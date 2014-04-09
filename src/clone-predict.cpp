//clone-predict.cpp

//own headers...
#include "emission.h"
#include "log-space.h"
#include "clone.h"

using namespace std;


void Clone::set_TransMat_cna(){
  if (TransMat_cna==NULL){
    TransMat_cna = new gsl_matrix * [cnaEmit->nSamples];
    for (int s=0; s<cnaEmit->nSamples;s++) TransMat_cna[s] = NULL;
  }
  for (int s=0; s<cnaEmit->nSamples;s++){
    if (TransMat_cna[s]!=NULL) gsl_matrix_free(TransMat_cna[s]);
    TransMat_cna[s] = gsl_matrix_alloc(nLevels,nLevels);
    set_TransMat_cna( TransMat_cna[s], cnaEmit->chr[s]);
  }
}

// only one clone can change its state,
// or clones in the same state can change in parallel
void Clone::set_TransMat_cna( gsl_matrix * Trans, int chr){
  double norm,p;
  gsl_vector_view row;
  int jumps,cni,cnf;
  for (int i=0; i<nLevels; i++){
    for (int j=0; j<nLevels; j++){
      jumps=0;
      p=1.0;
      for(int k=0; k < nClones; k++){
	/*if (copynumber[j][k] > maxtcn_per_clone[chr][k]){
	  jumps = 2;
	  break;
	  }*/
	if( copynumber[i][k] != copynumber[j][k]){
	  if ( jumps==0 ){
	    cni = copynumber[i][k];
	    cnf = copynumber[j][k];
	    jumps++;
	  }
	  //else if (cni != copynumber[i][k] || cnf != copynumber[j][k]){
	  else if (cni - cnf != copynumber[i][k] - copynumber[j][k]){
	    jumps++;
	  }
	  //if (copynumber[j][k] == 0) p*= 0.01;
	}
      }
      if (jumps <= 1){
	gsl_matrix_set( Trans, i, j, p);
      }
      else{
	gsl_matrix_set( Trans, i, j, 0.0);
      }
    }
    row  = gsl_matrix_row(Trans,i);
    norm = gsl_blas_dasum(&row.vector);
    if (norm <= 0) abort();
    gsl_vector_scale(&row.vector,1.0/norm);
  }
}


void Clone::set_TransMat_snv(){
  if (TransMat_snv==NULL){
    TransMat_snv = new gsl_matrix * [snvEmit->nSamples];
    for (int s=0; s<snvEmit->nSamples;s++) TransMat_snv[s] = NULL;
  }
  for (int s=0; s<snvEmit->nSamples;s++){
    if (TransMat_snv[s]!=NULL) gsl_matrix_free(TransMat_snv[s]);
    TransMat_snv[s] = gsl_matrix_alloc(nLevels,nLevels);
    set_TransMat_snv( TransMat_snv[s], snvEmit->chr[s]);
  }
}


// only one clone can change its state,
void Clone::set_TransMat_snv(gsl_matrix * Trans, int chr){
  double norm;
  gsl_vector_view row;
  int jumps;
  for (int i=0; i<nLevels; i++){
    for (int j=0; j<nLevels; j++){
      jumps=0;
      for(int k=0; k < nClones; k++){
	if (copynumber[j][k] > maxtcn_per_clone[chr][k]){
	  jumps = 2;
	  break;
	}
	if( copynumber[i][k]    != copynumber[j][k] 
	    && copynumber[i][k] <= maxtcn_per_clone[chr][k]){
	  jumps++;
	}
      }
      if (jumps <= 1){
	gsl_matrix_set(Trans,i, j, 1.0);
      }
      else{
	gsl_matrix_set(Trans,i, j, 0.0);
      }
    }
    row  = gsl_matrix_row(Trans,i);
    norm = gsl_blas_dasum(&row.vector);
    if (norm <= 0){
      cout<<"ERROR\n";
      abort();
    }
    gsl_vector_scale(&row.vector,1.0/norm);
  }
}


// predict step with transition matrix...
void Clone::predict( gsl_vector * prior, gsl_vector * post, Emission * myEmit, double pj, gsl_matrix * T){
  if (pj == 0.0){
    gsl_vector_memcpy( prior, post);//no jump possible
  }
  else{
    if (myEmit->log_space) for (int l=0;l<nLevels;l++) post->data[l] = exp(post->data[l]);
    gsl_vector_memcpy( prior, post);
    gsl_blas_dgemv( CblasTrans, pj, T, post, 1.0-pj, prior);
    if (myEmit->log_space){
      for (int l=0;l<nLevels;l++){
	prior->data[l] = prior->data[l]>0.0 ? log(prior->data[l]) : logzero;
      }
    }
  }
}



// predict step with convex mixing...
void Clone::predict( gsl_vector * prior, gsl_vector * post, Emission * myEmit, double pj, gsl_vector * flat){
  if (pj==0.0){
    gsl_vector_memcpy(prior,post);
  }
  else{
    gsl_vector_memcpy( prior, flat);
    if(pj<1.0){//convex combination
      if(myEmit->log_space){
	gsl_vector_add_constant(prior,log(pj));
	gsl_vector_add_constant(post,log(1.0-pj));
	log_vector_add(prior,post);
      }
      else{
	gsl_vector_scale(prior,pj);
	gsl_vector_scale(post,1.0-pj);
	gsl_vector_add(prior,post);
      }
    }
  }
}

void Clone::apply_maxtcn_mask(gsl_vector * prior, int chr, int log_space){
  for (int l=0;l<nLevels;l++){
    for (int j=0; j<nClones; j++){
      if ( copynumber[l][j] > maxtcn_per_clone[chr][j] ){
	prior->data[l] = 0.0;
	break;
      }
    }
  }
  double norm = gsl_blas_dasum(prior);
  if (norm <= 0.0) abort();
  gsl_vector_scale(prior,1.0/norm);
  if(log_space){
    for (int l=0;l<nLevels;l++){
      prior->data[l] = prior->data[l]>0.0 ? log(prior->data[l]) : logzero;
    }
  }
}

