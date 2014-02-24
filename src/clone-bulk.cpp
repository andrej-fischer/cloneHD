//clone-bulk.cpp

//own headers...
#include "emission.h"
#include "log-space.h"
#include "common-functions.h"
#include "clone.h"

using namespace std;



void Clone::allocate_bulk_mean(){//mean only...
  if (snvEmit->is_set==0) abort();
  //allocate prior...
  bulk_prior_mean = new double * [snvEmit->nSamples];
  for (int s=0; s<snvEmit->nSamples; s++){
    bulk_prior_mean[s] = new double [ snvEmit->nSites[s] ];
  }
  //allocate posterior...
  bulk_post_mean  = new double ** [nTimes];
  for (int t=0; t<nTimes; t++){
    bulk_post_mean[t] = new double * [snvEmit->nSamples];   
    for (int s=0; s<snvEmit->nSamples; s++){
      bulk_post_mean[t][s] = new double [snvEmit->nSites[s]];
    }
  }
  //set pointers to prior...
  bulk_mean = new double ** [nTimes];
  for (int t=0; t<nTimes; t++){
    bulk_mean[t] = new double * [snvEmit->nSamples];   
    for (int s=0; s<snvEmit->nSamples; s++){
      bulk_mean[t][s] = bulk_prior_mean[s];
    }
  }
}

void Clone::allocate_bulk_dist(){//distribution and mean...
  if (snvEmit->is_set==0) abort();
  //allocate prior...
  bulk_prior      = new gsl_matrix * [snvEmit->nSamples];
  bulk_prior_mean = new double * [snvEmit->nSamples];
  for (int s=0; s<snvEmit->nSamples; s++){
    bulk_prior[s] = gsl_matrix_calloc( snvEmit->nSites[s], bulkGrid + 1);
    bulk_prior_mean[s] = new double [ snvEmit->nSites[s] ];
  }
  //allocate posterior...
  bulk_post       = new gsl_matrix ** [nTimes];
  bulk_post_mean  = new double ** [nTimes];
  for (int t=0; t<nTimes; t++){
    bulk_post[t] = new gsl_matrix * [snvEmit->nSamples];   
    bulk_post_mean[t] = new double * [snvEmit->nSamples];   
    for (int s=0; s<snvEmit->nSamples; s++){
      bulk_post[t][s] = gsl_matrix_calloc( snvEmit->nSites[s], bulkGrid + 1);
      bulk_post_mean[t][s] = new double [snvEmit->nSites[s]];
    }
  }
  //set pointers to prior...
  bulk_dist = new gsl_matrix ** [nTimes];
  bulk_mean = new double ** [nTimes];
  for (int t=0; t<nTimes; t++){
    bulk_dist[t] = new gsl_matrix * [snvEmit->nSamples];   
    bulk_mean[t] = new double * [snvEmit->nSamples];   
    for (int s=0; s<snvEmit->nSamples; s++){
      bulk_dist[t][s] = bulk_prior[s];
      bulk_mean[t][s] = bulk_prior_mean[s];
    }
  }
}

void Clone::set_bulk_to_post(){
  for (int t=0; t<nTimes; t++){
    for (int s=0; s<snvEmit->nSamples; s++){
      bulk_mean[t][s] = bulk_post_mean[t][s];
      if (bulk_post != NULL)  bulk_dist[t][s] = bulk_post[t][s];
    }
  }
}

void Clone::set_bulk_to_prior(){
  for (int t=0; t<nTimes; t++){
    for (int s=0; s<snvEmit->nSamples; s++){
      bulk_mean[t][s] = bulk_prior_mean[s];
      if (bulk_prior != NULL)  bulk_dist[t][s] = bulk_prior[s];
    }
  }
}

// Bayesian update of the SNV bulk
void Clone::update_bulk(int sample){
  if (!snvEmit->is_set) abort();
  if (bulk_mean == NULL) abort();
  gsl_vector * bpost=NULL, *flat=NULL, *emit=NULL;
  gsl_vector_view bprior;
  bpost = gsl_vector_alloc(bulkGrid+1);
  flat  = gsl_vector_alloc(bulkGrid+1);
  gsl_vector_set_all(flat,1.0);
  //get SNV posterior...
  alpha_snv[sample]=NULL;
  gamma_snv[sample]=NULL;
  Clone::get_snv_posterior(sample);
  //update the bulk
  for (int time=0; time<nTimes; time++){   
    unsigned int n,N;
    for (int idx=0; idx<snvEmit->nSites[sample]; idx++){
      n = snvEmit->reads[time][sample][idx];
      N = snvEmit->depths[time][sample][idx];
      if (N>0){
	emit = snvEmit->EmitProb[N][n];
	if (bulk_prior == NULL){ 
	  Clone::get_bulk_post_dist( flat, bpost, emit, time, sample, idx);
	  bulk_post_mean[time][sample][idx] 
	    = 0.5*(bulk_prior_mean[sample][idx] + snvEmit->xgrid[gsl_vector_max_index(bpost)]);
	}
	else{
	  bprior = gsl_matrix_row( bulk_prior[sample], idx);
	  Clone::get_bulk_post_dist( &bprior.vector, bpost, emit, time, sample, idx);
	  gsl_matrix_set_row( bulk_post[time][sample], idx, bpost);
	  bulk_post_mean[time][sample][idx] = get_mean( bpost, 0.0, 1.0);
	}
      }
      else{//no observation
	bulk_post_mean[time][sample][idx] = bulk_prior_mean[sample][idx];
	if (bulk_prior != NULL){
	  bprior = gsl_matrix_row( bulk_prior[sample], idx);
	  gsl_matrix_set_row( bulk_post[time][sample], idx, &bprior.vector);
	}
      }
    }
  }
  //clean up...
  gsl_matrix_free(gamma_snv[sample]);
  gamma_snv[sample] = NULL;
  if (bpost!=NULL) gsl_vector_free(bpost);
  if (flat!=NULL)  gsl_vector_free(flat);
}

// emission probability
void Clone::get_bulk_post_dist( gsl_vector * bprior, gsl_vector * bpost, gsl_vector * emit, int time, int sample, int idx){
  double dx = snvEmit->dx;
  double prob=0, x=0, y=0, val=0 , f1=0, f2=0, nu=0;
  int old=-1, i1=0, i2=0;
  int chr = snvEmit->chr[sample];
  double ncn = (double) normal_copy[chr];
  int evt = snvEmit->event_of_idx[sample][idx];
  double total_cn = 0;
  if (snvEmit->mean_tcn != NULL){
    total_cn = snvEmit->mean_tcn[time][sample][evt];
  }
  else{
    total_cn = tcn[chr][time][level_of[chr]];
  }
  gsl_vector_set_zero(bpost);
  for (int i=0; i <= bulkGrid; i++){
    y = double(i)*dx * (1.0 - purity[time]) * ncn / total_cn;
    prob=0.0;
    old = -1;
    for (int j=0; j<nLevels; j++){
      x = y + clone_spectrum[time][j] / total_cn;
      if (x > 1.0) continue;
      i1 = int(x/dx);
      nu = x/dx - double(i1);
      if (i1 != old){
	f1 = emit->data[i1];
	i2 = (i1<bulkGrid) ? i1+1 : i1;
	f2 = emit->data[i2];
	old=i1;
      }
      val = (1.0-nu)*f1 + nu*f2;
      val *= gsl_matrix_get( gamma_snv[sample], evt, j);
      prob += val;
    }
    if (prob != prob) abort();
    gsl_vector_set( bpost, i, prob);
  }
  gsl_vector_mul( bpost, bprior);
  double norm = gsl_blas_dasum(bpost) - 0.5*(bpost->data[0] + bpost->data[bulkGrid]);
  norm *= 1.0/double(bulkGrid);
  if (norm<= 0.0) abort();
  gsl_vector_scale( bpost, 1.0/norm);
}

//get the minimum mean bulk freq for all segments
void Clone::get_bulk_min(){
  if (bulk_min==NULL){//allocate
    bulk_min = new double ** [nTimes];
    for (int t=0; t<nTimes; t++){
      bulk_min[t] = new double * [snvEmit->nSamples];
      for (int s=0; s<snvEmit->nSamples; s++){
	bulk_min[t][s] = new double [snvEmit->nEvents[s]];
      }
    }
  }
  for (int t=0; t<nTimes; t++){
    for (int s=0; s<snvEmit->nSamples; s++){
      for (int evt=0; evt < snvEmit->nEvents[s]; evt++){
	if (bulk_fix >= 0.0){
	  bulk_min[t][s][evt] = bulk_fix;
	}
	else if (bulk_mean != NULL){
	  bulk_min[t][s][evt] = 1.1;
	  int idxi = snvEmit->idx_of_event[s][evt];
	  int idxf = (evt<snvEmit->nEvents[s]-1) ? snvEmit->idx_of_event[s][evt+1]-1 : snvEmit->nSites[s]-1;
	  for (int idx=idxi; idx<=idxf; idx++){
	    bulk_min[t][s][evt] = min( bulk_min[t][s][evt], bulk_mean[t][s][idx]);
	  }
	  bulk_min[t][s][evt] = max( bulk_min[t][s][evt], 1.0e-3);
	}
	else{
	  abort();
	}
      }
    }
  }
}
