//clone-prior.cpp

//own headers...
#include "emission.h"
#include "log-space.h"
#include "clone.h"

using namespace std;


// CNA prior (only used for chr entry or w/o correlations)...
void Clone::set_cna_prior( gsl_vector * prior, int sample){
  // cna_pen_norm: penalty for being  different from the normal copy number
  // cna_pen_diff: penalty for having different copynumbers across clones
  // cna_pen_zero: penalty for zero total copies
  if (cnaEmit->is_set == 0) abort();
  if (nClones==0){
    gsl_vector_set_all(prior,1.0);
  }
  else{
    std::set<int> cns;
    int chr = cnaEmit->chr[sample];
    int ncn = normal_copy[chr];
    for (int i=0; i<nLevels; i++){
      cns.clear();
      double p = 1.0;
      for (int j=0; j<nClones; j++){
	p *= pow( cna_pen_norm, abs(copynumber[i][j] - ncn));
	if (copynumber[i][j] == 0) p *= cna_pen_zero;
	if (copynumber[i][j] > maxtcn_per_clone[chr][j]) p = 0.0;
	cns.insert( copynumber[i][j] );
      }
      p *= pow( cna_pen_diff, (int) cns.size());
      gsl_vector_set( prior, i, p);
    } 
    //normalize and logify...
    double norm = gsl_blas_dasum(prior);
    if (norm <= 0.0 || norm != norm) abort();
    gsl_vector_scale( prior, 1.0 / norm);
    if (cnaEmit->log_space){
      for (int l=0; l<nLevels;l++){ 
	prior->data[l] = prior->data[l] > 0.0 ? log(prior->data[l]) : logzero;
      }
    }
  }
}

//used in SNV-only mode, w/o correlation and w/o cn-info...
void Clone::initialize_snv_prior_param(){// SNV prior, conditional on max-tcn
  if (initial_snv_prior_param != NULL) gsl_matrix_free(initial_snv_prior_param);
  initial_snv_prior_param = NULL;
  if (nClones == 0) return;
  initial_snv_prior_param = gsl_matrix_calloc( maxtcn+1, maxtcn+1);
  gsl_matrix_set( initial_snv_prior_param, 0, 0, snv_fpr);
  double p = snv_pen_high;// penalty for higher genotypes
  for (int cn=1; cn <= maxtcn; cn++){
    if ( all_maxtcn.count(cn) == 0 ) continue;
    gsl_vector_view subrow = gsl_matrix_subrow( initial_snv_prior_param, cn, 0, cn+1);
    gsl_vector_set( &subrow.vector, 0, p);//P0=P1>P2>... or P00=P10=P10=P11>P20...
    for (int i=1; i<=cn; i++) gsl_vector_set( &subrow.vector, i, pow( p, i));
    gsl_vector_scale( &subrow.vector, 1.0 / gsl_blas_dasum(&subrow.vector) );
  }
  //set above fixed priors
  Clone::set_snv_prior(initial_snv_prior_param);
}


//SNV-only mode, w/o cn-info...
void Clone::set_snv_prior( gsl_matrix * snv_prior_param){
  double fpr = gsl_matrix_get( snv_prior_param, 0, 0);
  snv_prior.clear();
  std::map<int, vector<int> >::iterator it;
  gsl_vector * mem = gsl_vector_calloc(nLevels);    
  for (it=maxtcn_per_clone.begin(); it != maxtcn_per_clone.end(); it++){
    int chr = it->first;
    mem->data[0] = 0;
    for (int i=1; i<nLevels; i++){
      double p=1.0;
      for (int j=0; j<nClones; j++){
	int limit = maxtcn_per_clone[chr][j];
	if (copynumber[i][j] <= limit){
	  p *= limit==0 ? 1.0 : gsl_matrix_get( snv_prior_param, limit, copynumber[i][j]);
	}
	else{
	  p = 0.0;
	  break;
	}
      }
      //gsl_vector_set( snv_prior[chr], i, p);
      mem->data[i] = p;
    }
    //normalize and logify...
    double norm = gsl_blas_dasum(mem);
    if (norm<=0) abort();
    gsl_vector_scale( mem, (1.0-fpr) / norm);
    gsl_vector_set( mem, 0, fpr);
    if (snvEmit->log_space){
      for (int l=0; l<nLevels; l++){
	double val = gsl_vector_get( mem, l);
	gsl_vector_set( mem, l, val>0 ? log(val) : logzero);
      }
    }
    snv_prior[chr] = gsl_vector_calloc(nLevels);
    gsl_vector_memcpy(snv_prior[chr],mem);
  }
  gsl_vector_free(mem);
}

//CNA + BAF (+SNV) mode...
void Clone::set_baf_prior_map(){
  if ( baf_prior_map == NULL){
    baf_prior_map = gsl_matrix_alloc( maxtcn+1, maxtcn+1);
  }
  gsl_matrix_set_zero(baf_prior_map); 
  double f = 0;
  for (int cn=0; cn <= maxtcn; cn++){
    for (int bcn=0; bcn <= cn; bcn++){//penalty for complex chromosome status 
      f = pow( baf_pen_comp, int( fabs(double(bcn) - 0.5*double(cn))));
      gsl_matrix_set( baf_prior_map, bcn, cn, f);
    }
    //normalize...
    gsl_vector_view col = gsl_matrix_column( baf_prior_map,cn);
    double norm = gsl_blas_dasum(&col.vector);
    if (norm <= 0.0) abort();
    gsl_vector_scale( &col.vector, 1.0 / norm);
  }
}

//CNA (+BAF) + SNV mode...
void Clone::set_snv_prior_map(){//either via BAF or else via CNA
  if (nClones == 0) abort();
  if (cnaEmit->is_set == 0) abort();
  //allocate
  if (bafEmit->is_set){//via CNA + BAF posterior...
    if ( snv_prior_from_cna_baf_map == NULL){
      snv_prior_from_cna_baf_map = new gsl_matrix * [maxtcn+1];
      for (int cn=0; cn <= maxtcn; cn++){ 
	snv_prior_from_cna_baf_map[cn] = gsl_matrix_alloc( maxtcn+1, maxtcn+1);
      }
    }
    for (int cn=0; cn <= maxtcn; cn++){//cn=total cn
      gsl_matrix_set_zero( snv_prior_from_cna_baf_map[cn]); 
      for (int j=0; j<=cn; j++){//j=minor cn
	for (int i=0; i<= cn; i++){//penalty for multiple hit SNVs
	  double pen = pow( snv_pen_mult, max( 0, i - max(j,cn-j)) );
	  gsl_matrix_set( snv_prior_from_cna_baf_map[cn], i, j, pen);
	}
	//normalize...
	gsl_vector_view col = gsl_matrix_column( snv_prior_from_cna_baf_map[cn], j);
	double norm = gsl_blas_dasum(&col.vector);
	if (norm <=0.0) abort();
	gsl_vector_scale( &col.vector, 1.0 / norm);
      }
    } 
  }
  //via CNA posterior only...
  if ( snv_prior_from_cna_map == NULL){//allocate
    snv_prior_from_cna_map = gsl_matrix_alloc( maxtcn+1, maxtcn+1);
  }
  gsl_matrix_set_zero( snv_prior_from_cna_map);  
  double pen = snvEmit->connect ? 1.0 : snv_pen_high;// penalty for high genotypes 
  for (int cn=0; cn <= maxtcn; cn++){
    gsl_matrix_set( snv_prior_from_cna_map, 0, cn, pen);
    for (int i=1; i<=cn; i++){
      gsl_matrix_set( snv_prior_from_cna_map, i, cn, pow(pen,i));
    }
    //normalize...
    gsl_vector_view col = gsl_matrix_column( snv_prior_from_cna_map, cn);
    double norm = gsl_blas_dasum(&col.vector);
    if (norm <=0.0) abort();
    gsl_vector_scale( &col.vector, 1.0 / norm);
  }
}




// CNA + BAF (+SNV) mode
void Clone::get_baf_prior_from_cna_post(gsl_vector * prior, gsl_vector * post){
  gsl_vector * post_per_clone  = gsl_vector_calloc( nClones*(maxtcn+1) );
  gsl_matrix * prior_per_clone = gsl_matrix_calloc( nClones, maxtcn+1);
  gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, post, 0.0, post_per_clone);
  gsl_vector_view po_pc,pr_pc;
  for (int i=0; i<nClones; i++){
    po_pc = gsl_vector_subvector( post_per_clone, i*(maxtcn+1), maxtcn+1);
    pr_pc = gsl_matrix_row( prior_per_clone, i);
    gsl_blas_dgemv( CblasNoTrans, 1.0, baf_prior_map, &po_pc.vector, 0.0, &pr_pc.vector);
  }
  gsl_vector_set_all( prior, 1.0);
  for (int i=0; i<nLevels; i++){
    for (int j=0; j<nClones; j++){
      prior->data[i] *= gsl_matrix_get( prior_per_clone, j, copynumber[i][j]);
    }
  }
  if (bafEmit->log_space){//log-transform?
    for (int l=0; l<nLevels; l++){
      prior->data[l] = prior->data[l] > 0.0 ? log(prior->data[l]) : logzero;
    }
  }
  gsl_matrix_free(prior_per_clone);
  gsl_vector_free(post_per_clone);
}


// CNA + SNV mode...
void Clone::get_snv_prior_from_cna_post(gsl_vector * prior, gsl_vector * cnapost){
  gsl_vector * cnapostpc = gsl_vector_calloc( nClones*(maxtcn+1));
  gsl_matrix * snv_prpc  = gsl_matrix_calloc( nClones, maxtcn+1);
  gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, cnapost, 0.0, cnapostpc);
  gsl_vector_view cnappc,prpc;
  for (int i=0; i<nClones; i++){
    cnappc = gsl_vector_subvector( cnapostpc, i*(maxtcn+1), maxtcn+1);
    prpc = gsl_matrix_row( snv_prpc, i);
    gsl_blas_dgemv( CblasNoTrans, 1.0, snv_prior_from_cna_map, &cnappc.vector, 0.0, &prpc.vector);
  }
  Clone::apply_snv_prpc( prior, snv_prpc, cnapost->data[0]);
  gsl_matrix_free(snv_prpc);
  gsl_vector_free(cnapostpc);
}




//CNA + BAF + SNV mode...
void Clone::get_snv_prior_from_cna_baf_post(gsl_vector * prior, gsl_vector * cnapost, gsl_vector * bafpost){
  gsl_vector * cnapostpc  = gsl_vector_calloc( nClones*(maxtcn+1));
  gsl_vector * bafpostpc  = gsl_vector_calloc( nClones*(maxtcn+1));
  gsl_matrix * snv_prpc   = gsl_matrix_calloc( nClones, maxtcn+1);
  gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, cnapost, 0.0, cnapostpc);
  gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, bafpost, 0.0, bafpostpc);
  gsl_vector_view prpc;
  gsl_vector * bafppc = gsl_vector_alloc(maxtcn+1);
  for (int i=0; i<nClones; i++){
    prpc   = gsl_matrix_row( snv_prpc, i);
    for (int cn=0; cn<=maxtcn; cn++){
      gsl_vector_set_zero(bafppc); 
      double norm=0;
      for (int j=0; j<=cn;j++){
	bafppc->data[j] = gsl_vector_get( bafpostpc, i*(maxtcn+1)+j) + 1.0e-10;
	norm += bafppc->data[j];
      }
      if (norm <= 0.0) abort();
      gsl_vector_scale(bafppc,1.0/norm);
      double pcna = gsl_vector_get( cnapostpc, i*(maxtcn+1) + cn);
      gsl_blas_dgemv( CblasNoTrans, pcna, snv_prior_from_cna_baf_map[cn], bafppc, 1.0, &prpc.vector);
    }
  }
  Clone::apply_snv_prpc( prior, snv_prpc, cnapost->data[0]);
  gsl_matrix_free(snv_prpc);
  gsl_vector_free(cnapostpc);
  gsl_vector_free(bafpostpc);
  gsl_vector_free(bafppc);
}


//CNA (+BAF) + SNV mode...
void Clone::apply_snv_prpc( gsl_vector * prior, gsl_matrix * snv_prpc,  double pzero){
  gsl_vector_set_all( prior, 1.0);
  for (int level=0; level<nLevels; level++){
    for (int j=0; j<nClones; j++){
      prior->data[level] *= gsl_matrix_get( snv_prpc, j, copynumber[level][j]);
    }
  }
  //normalize...
  if ( !snvEmit->connect ){
    prior->data[0] = 0.0;
    double norm = gsl_blas_dasum(prior);
    if (norm > 0.0 ) gsl_vector_scale(prior, (1.0-snv_fpr)*(1.0-pzero) / norm);
    prior->data[0] = 1.0 - (1.0-snv_fpr)*(1.0-pzero);//SNV false positive rate
  }
  else{
    double norm = gsl_blas_dasum(prior);
    if (norm <=0.0 || norm != norm) abort();
    gsl_vector_scale(prior, 1.0 / norm);
  }
  //log-transform?
  if (snvEmit->log_space){
    for (int l=0; l<nLevels; l++){
      prior->data[l] = prior->data[l] > 0.0 ? log(prior->data[l]) : logzero;
    }
  }
}




//get mean total copy number...
void Clone::get_mean_tcn(int sample){//only ever used for cnaEmit
  int chr = cnaEmit->chr[sample];
  if (nClones == 0){   
    for (int t=0; t<nTimes; t++){
      for (int evt=0; evt < cnaEmit->nEvents[sample]; evt++){
       cnaEmit->mean_tcn[t][sample][evt] = tcn[chr][t][0];
      }
    }
  }
  else{//nClones > 0
    if (gamma_cna[sample] == NULL) abort();
    double mn;
    for (int t=0; t<nTimes; t++){
      for (int evt=0; evt < cnaEmit->nEvents[sample]; evt++){
       gsl_vector_view TCN  = gsl_vector_view_array( tcn[chr][t], nLevels);
       gsl_vector_view post = gsl_matrix_row( gamma_cna[sample], evt);
       gsl_blas_ddot( &TCN.vector, &post.vector, &mn);
       if (mn<=0.0){
	 cout<<"ERROR\n";
	 abort();
       }
       cnaEmit->mean_tcn[t][sample][evt] = mn;
      }
    }
  }
}

void Clone::map_mean_tcn( Emission * fromEmit, int fromSample, Emission * toEmit){
  int fromChr = fromEmit->chr[fromSample];
  if (toEmit->chrs.count(fromChr) == 0) abort();
  if (toEmit->mean_tcn==NULL) abort();
  int toSample = toEmit->idx_of[fromChr];
  for (int evt=0; evt < toEmit->nEvents[toSample]; evt++){
    int idx = toEmit->idx_of_event[toSample][evt];
    int fromEvt = toEmit->Event_of_idx[toSample][idx];
    for (int t=0; t<nTimes; t++){//mean total copynumber
      toEmit->mean_tcn[t][toSample][evt] = fromEmit->mean_tcn[t][fromSample][fromEvt];
    }
  }
}



void Clone::get_avail_cn(Emission * myEmit, int sample){
  int chr = myEmit->chr[sample];
  if (nClones==0){
    for (int t=0; t<nTimes; t++){
      for (int evt=0; evt < myEmit->nEvents[sample]; evt++){
       for (int cn=0; cn<=maxtcn; cn++){
         if (myEmit==cnaEmit){
           myEmit->av_cn[t][sample][evt][cn] = (cn <= normal_copy[chr]) ? 1.0 : 0.0;
         }
         else if (myEmit==bafEmit){
           myEmit->av_cn[t][sample][evt][cn] = (cn <= 1) ? 1.0 : 0.0;
         }
       }
      }
    }
  }
  else{//nClones > 0
    gsl_matrix * gamma = NULL;
    if (myEmit == cnaEmit) gamma = gamma_cna[sample];
    if (myEmit == bafEmit) gamma = gamma_baf[sample];
    if (gamma== NULL) abort();
    double val=0;
    gsl_vector * post_per_clone = gsl_vector_alloc((maxtcn+1)*nClones);
    int * cnest = new int [nClones];//conservative estimate
    for (int t=0; t<nTimes; t++){
      for (int evt=0; evt < myEmit->nEvents[sample]; evt++){
        gsl_vector_view post = gsl_matrix_row(gamma,evt);
        gsl_blas_dgemv(CblasNoTrans, 1.0, margin_map, &post.vector, 0.0, post_per_clone);
        for (int j=0; j<nClones; j++){
          val=0;
          for (int cn=0; cn<=maxtcn; cn++){
            val += gsl_vector_get( post_per_clone, j*(maxtcn+1) + cn);
            if (val > 0.99){
              cnest[j] = cn;
              break;
            }
          }
        }
        for (int cn=0; cn<=maxtcn; cn++){
          double av = 0.0;
          if (myEmit == cnaEmit && cn <= normal_copy[chr]) av += 1.0-purity[t];
          if (myEmit == bafEmit && cn <= 1) av += 1.0-purity[t];
          for (int j=0;j<nClones;j++){
            if (cn<=cnest[j]) av += gsl_matrix_get(freqs,t,j);
          }
          myEmit->av_cn[t][sample][evt][cn] = av;
        }
      }
    }
    gsl_vector_free(post_per_clone);
  }
}


void Clone::get_snv_prior_from_av_cn(gsl_vector * prior, int sample, int evt){
  if (snvEmit->av_cn==NULL) abort();
  int snvChr = snvEmit->chr[sample];
  int found=0;
  prior->data[0] = 0;
  for (int l=1; l<nLevels; l++){
    found=0;
    prior->data[l] = (snv_prior[snvChr])->data[l];
    for (int t=0;t<nTimes;t++){//genotype not available?
      for (int cn=0; cn<=maxtcn; cn++){
       if (cn_usage[t][cn][l] > snvEmit->av_cn[t][sample][evt][cn]){
         found = 1;
	 prior->data[l] *= snv_pen_mult;//multiple hit SNV
         break;
       }
       if (found) break;
      }
      if (found) break;
    }
  }
  double norm = gsl_blas_dasum(prior);
  if (norm <= 0.0) abort();
  gsl_vector_scale(prior, (1.0-snv_fpr)/norm);
  prior->data[0] = snv_fpr;
}
