//clone-llh.cpp

//own headers...
#include "emission.h"
#include "log-space.h"
#include "clone.h"

using namespace std;



double Clone::get_all_total_llh(){
  if (bafEmit->is_set || snvEmit->is_set){
    alpha_cna = new gsl_matrix * [cnaEmit->nSamples];
    gamma_cna = new gsl_matrix * [cnaEmit->nSamples];
    for ( int s=0; s < cnaEmit->nSamples; s++){
      alpha_cna[s] = gsl_matrix_calloc( cnaEmit->nEvents[s], nLevels);
      gamma_cna[s] = gsl_matrix_calloc( cnaEmit->nEvents[s], nLevels);
    }
    save_cna_alpha = 1;
    if (bafEmit->is_set && snvEmit->is_set){
      alpha_baf = new gsl_matrix * [bafEmit->nSamples];
      gamma_baf = new gsl_matrix * [bafEmit->nSamples];
      for ( int s=0; s < bafEmit->nSamples; s++){
	alpha_baf[s] = gsl_matrix_calloc( bafEmit->nEvents[s], nLevels);
	gamma_baf[s] = gsl_matrix_calloc( bafEmit->nEvents[s], nLevels);
      }
      save_baf_alpha = 1;
    }
    else{
      save_baf_alpha = 0;
    }
  }
  else{
    save_cna_alpha = 0;
  }
  get_gofs = 0;
  save_snv_alpha = 0;
  total_llh     = 0.0;
  cna_total_llh = 0.0;
  baf_total_llh = 0.0;
  snv_total_llh = 0.0;
  double * llhs = NULL;
  int sample;
#ifdef _OPENMP
  int nt = min( cnaEmit->nSamples, omp_get_max_threads());
#pragma omp parallel for schedule( dynamic, 1) default(shared) num_threads(nt)
#endif
  for ( sample=0; sample < cnaEmit->nSamples; sample++){
    double llh=0,ent=-1;   
    Clone::do_cna_Fwd( sample, llh, llhs);
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      cna_total_llh += llh;
    }
    if (save_cna_alpha==1){
      int cnaChr = cnaEmit->chr[sample];
      Clone::do_cna_Bwd( sample, ent);
      Clone::get_mean_tcn(sample);
      if ( bafEmit->is_set && bafEmit->chrs.count(cnaChr) == 1){
	Clone::map_mean_tcn( cnaEmit, sample, bafEmit);	
	if ( snvEmit->is_set && snvEmit->chrs.count(cnaChr) == 1){
	  int bafsample = bafEmit->idx_of[cnaEmit->chr[sample]];
	  Clone::map_mean_tcn( bafEmit, bafsample, snvEmit);
	}
      }
      else if ( snvEmit->is_set && snvEmit->chrs.count(cnaChr) == 1){
	Clone::map_mean_tcn( cnaEmit, sample, snvEmit);
      }
      gsl_matrix_free(alpha_cna[sample]);
      alpha_cna[sample] = NULL;
    }
  }//END PARALLEL FOR
  //
  // BAF
  if ( bafEmit->is_set ){
#ifdef _OPENMP
    int nt = min( bafEmit->nSamples, omp_get_max_threads());
#pragma omp parallel for schedule( dynamic, 1) default(shared) num_threads(nt)
#endif
    for ( sample=0; sample < bafEmit->nSamples; sample++){//START PARALLEL FOR
      double llh=0,ent=-1;
      Clone::do_baf_Fwd( sample, llh, llhs);
#ifdef _OPENMP
#pragma omp critical
#endif
      {
	baf_total_llh += llh;
      }
      if (save_baf_alpha==1){
	Clone::do_baf_Bwd( sample, ent);
	gsl_matrix_free(alpha_baf[sample]);
	alpha_baf[sample] = NULL;
      }
    }//END PARALLEL FOR
  }
  //
  // SNV
  if( snvEmit->is_set ){
#ifdef _OPENMP
    int nt = min( snvEmit->nSamples,  omp_get_max_threads());
#pragma omp parallel for schedule( dynamic, 1) default(shared) num_threads(nt)
#endif
    for ( sample=0; sample < snvEmit->nSamples; sample++){//START PARALLEL FOR
      double llh = 0.0;
      Clone::do_snv_Fwd(sample, llh, llhs);
#ifdef _OPENMP
#pragma omp critical
#endif
      {
	snv_total_llh += llh;
      }
    }//END PARALLEL FOR
  }
  //cleanup...
  if (cnaEmit->is_set && save_cna_alpha==1){
    for ( sample=0; sample < cnaEmit->nSamples; sample++){
      gsl_matrix_free(gamma_cna[sample]);
    }
    delete [] alpha_cna;
    delete [] gamma_cna;
    alpha_cna = NULL;
    gamma_cna = NULL;
  }
  if (bafEmit->is_set && save_baf_alpha==1){
    for ( sample=0; sample < bafEmit->nSamples; sample++){
      gsl_matrix_free(gamma_baf[sample]);
    }
    delete [] alpha_baf;
    delete [] gamma_baf;
    alpha_baf = NULL;
    gamma_baf = NULL;
  }
  total_llh = cna_total_llh + baf_total_llh + snv_total_llh;
  return(total_llh);
}




double Clone::get_cna_total_llh(){
  int sample;
  save_cna_alpha = 0;
  cna_total_llh  = 0.0;
  double * llhs = NULL;
#ifdef _OPENMP
  int nt = min( cnaEmit->nSamples,  omp_get_max_threads());
#pragma omp parallel for schedule( dynamic, 1) default(shared) num_threads(nt)
#endif
  for ( sample=0; sample < cnaEmit->nSamples; sample++){
    double llh;
    Clone::do_cna_Fwd( sample, llh, llhs);
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      cna_total_llh += llh;
    }
  }//END PARALLEL FOR
  return(cna_total_llh);
}


double Clone::get_baf_total_llh(){
  if ( nClones > 0 && cnaEmit->is_set && gamma_cna == NULL ) abort();
  save_baf_alpha = 0;
  baf_total_llh  = 0.0;
  double * llhs = NULL;
  int sample;
#ifdef _OPENMP
  int nt = min( bafEmit->nSamples,  omp_get_max_threads());
#pragma omp parallel for schedule( dynamic, 1) default(shared) num_threads(nt)
#endif
  for ( sample=0; sample< bafEmit->nSamples; sample++){
    double llh;
    Clone::do_baf_Fwd( sample, llh, llhs);
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      baf_total_llh += llh;
    }
  }
  return(baf_total_llh);
}




double Clone::get_snv_total_llh(){
  if ( nClones > 0 && cnaEmit->is_set && gamma_cna == NULL ) abort();
  int sample;
  save_snv_alpha = 0;
  snv_total_llh  = 0.0;
  double * llhs  = NULL;
#ifdef _OPENMP
  int nt = min( snvEmit->nSamples,  omp_get_max_threads());
#pragma omp parallel for schedule( dynamic, 1) default(shared) num_threads(nt)
#endif
  for ( sample=0; sample< snvEmit->nSamples; sample++){
    double llh=0;
    Clone::do_snv_Fwd( sample, llh, llhs);
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      snv_total_llh += llh;
    }
  }
  return(snv_total_llh);
}




double Clone::get_cna_posterior(int sample){
  double llh=0;
  double * llhs = NULL;
  save_cna_alpha = 1;
  //set fw-bw arrays    
  if (alpha_cna[sample] == NULL) 
    alpha_cna[sample] = gsl_matrix_calloc( cnaEmit->nEvents[sample], nLevels);
  if (gamma_cna[sample] == NULL) 
    gamma_cna[sample] = gsl_matrix_calloc( cnaEmit->nEvents[sample], nLevels);
  Clone::do_cna_Fwd( sample, llh, llhs);
  Clone::do_cna_Bwd( sample, cna_total_ent);
  gsl_matrix_free(alpha_cna[sample]);
  alpha_cna[sample] = NULL;
  return(llh);
}

double Clone::get_baf_posterior(int sample){
  if (cnaEmit->is_set){
    int cna_sample = cnaEmit->idx_of[bafEmit->chr[sample]];
    if ( nClones > 0 && (gamma_cna == NULL || gamma_cna[cna_sample] == NULL)){
      abort();
    }
  }
  double llh=0;
  double * llhs = NULL;
  save_baf_alpha = 1;
  //set fw-bw arrays    
  if (alpha_baf[sample] == NULL) 
    alpha_baf[sample] = gsl_matrix_calloc( bafEmit->nEvents[sample], nLevels);
  if (gamma_baf[sample] == NULL) 
    gamma_baf[sample] = gsl_matrix_calloc( bafEmit->nEvents[sample], nLevels);
  Clone::do_baf_Fwd( sample, llh, llhs);
  Clone::do_baf_Bwd( sample, baf_total_ent);
  gsl_matrix_free(alpha_baf[sample]);
  alpha_baf[sample] = NULL;
  return(llh);
}


double Clone::get_snv_posterior(int sample){
  if (cnaEmit->is_set){
    int cna_sample = cnaEmit->idx_of[snvEmit->chr[sample]];
    if ( nClones > 0 && (gamma_cna == NULL || gamma_cna[cna_sample] == NULL)){
      abort();
    }
  }
  double llh=0;
  double * llhs = NULL;
  save_snv_alpha = 1;
  //set fw-bw arrays    
  if (alpha_snv[sample] == NULL) 
    alpha_snv[sample] = gsl_matrix_calloc( snvEmit->nEvents[sample], nLevels);
  if (gamma_snv[sample] == NULL) 
    gamma_snv[sample] = gsl_matrix_calloc( snvEmit->nEvents[sample], nLevels);
  Clone::do_snv_Fwd( sample, llh, llhs);
  Clone::do_snv_Bwd( sample, snv_total_ent);
  gsl_matrix_free(alpha_snv[sample]);
  alpha_snv[sample] = NULL;
  return(llh);
}
