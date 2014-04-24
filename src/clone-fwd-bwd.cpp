//clone-fwd-bwd.cpp

//own headers...
#include "emission.h"
#include "log-space.h"
#include "clone.h"

using namespace std;




double Clone::entropy(gsl_vector * x){
  double H = 0.0;
  for (int i=0; i<(int) x->size; i++){
    if (x->data[i] > 0.0) H -= x->data[i] * log(x->data[i]);
  }
  return(H);
}

void Clone::scale_prior( gsl_vector*& prior, int n){
  gsl_vector_scale( prior, double(n));//log-space!
  double norm = log_vector_norm(prior);
  gsl_vector_add_constant( prior, -norm);
}

void Clone::combine_prior( gsl_vector*& prior, gsl_vector*& mem, int n){
  double norm;
  if ( gsl_vector_max(prior) <= 0.0){//log-space
    if (n > 1) Clone::scale_prior(mem,n);
    gsl_vector_add( prior, mem);
    norm = log_vector_norm(prior);
    gsl_vector_add_constant(prior,-norm);
  }
  else{
    gsl_vector_mul( prior, mem);
    norm = gsl_blas_dasum(prior);
    if (norm!=norm || norm < 0.0) abort();
    gsl_vector_scale( prior, 1.0/norm);
  }
}

//*** CNA FWD/BWD ***************************************************
void Clone::do_cna_Fwd( int sample, double& llh, double*& llhs){
  gsl_vector * entry = gsl_vector_alloc(nLevels);
  gsl_vector * mem   = gsl_vector_alloc(nLevels);
  gsl_vector * prior = gsl_vector_alloc(nLevels);
  gsl_vector * post  = gsl_vector_alloc(nLevels);
  gsl_matrix * Trans = NULL;
  if (nClones>0){//preparations...
    Clone::set_cna_prior( entry, sample);
    if (cnaEmit->connect){
      Trans = gsl_matrix_alloc( nLevels, nLevels);
      gsl_matrix_memcpy(Trans,TransMat_cna[sample]);
    }
  }
  double norm = 0.0, pj=1.0;
  int idx=0, nidx=0, last_evt = cnaEmit->nEvents[sample]-1;
  llh = 0.0;
  for (int evt=0; evt <= last_evt ; evt++){
    idx = cnaEmit->idx_of_event[sample][evt];
    nidx = (evt < last_evt) ? cnaEmit->idx_of_event[sample][evt+1] : cnaEmit->nSites[sample];
    //***PREDICT STEP***
    if (nClones > 0){
      if (cnaEmit->connect && evt > 0){//connect to the left...
	pj = cnaEmit->pjump[sample][idx];
	Clone::predict( prior, post, cnaEmit, pj, Trans);	
	gsl_vector_memcpy( mem, entry);
	Clone::combine_prior( prior, mem, nidx-idx);
      }
      else{
	gsl_vector_memcpy( prior, entry);
	if (nidx-idx > 1) Clone::scale_prior(prior,nidx-idx);
      }
    }
    else{//nClones == 0
      gsl_vector_set_all( prior, cnaEmit->log_space ? 0.0 : 1.0);
    }
    //***UPDATE***
    norm = Clone::update( prior, post, cnaEmit, sample, evt, llhs);
    llh += norm;
    if (save_cna_alpha == 1) gsl_matrix_set_row( alpha_cna[sample], evt, post);
  }
  // cleanup    
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(entry);
  gsl_vector_free(mem);
  if (Trans != NULL) gsl_matrix_free(Trans);
}

void Clone::do_cna_Bwd(int sample, double& ent){
  if (alpha_cna[sample] == NULL || gamma_cna[sample] == NULL) abort();
  gsl_vector * prior = gsl_vector_alloc(nLevels);
  gsl_vector * post  = gsl_vector_alloc(nLevels);
  gsl_vector * entry = gsl_vector_alloc(nLevels);
  gsl_vector * mem   = gsl_vector_alloc(nLevels);
  gsl_matrix * Trans = NULL;
  //int cnaChr = cnaEmit->chr[sample];
  if (nClones>0){
    Clone::set_cna_prior( entry, sample);
    if (cnaEmit->connect){
      Trans = gsl_matrix_alloc(nLevels,nLevels);
      gsl_matrix_memcpy(Trans,TransMat_cna[sample]);
    }
  }
  double * llhs = NULL;
  int idx=0, nidx=0, nloci=1;
  gsl_vector_view alph;
  double pj = 1.0, norm=0;
  int last_evt = cnaEmit->nEvents[sample]-1;
  int last_idx = cnaEmit->idx_of_event[sample][last_evt];
  for (int evt = last_evt; evt >= 0; evt--){
    idx   = cnaEmit->idx_of_event[sample][evt];
    nidx  = (evt < last_evt) ? cnaEmit->idx_of_event[sample][evt+1] : cnaEmit->nSites[sample];
    nloci = nidx-idx;
    //***PREDICTION STEP***
    if (nClones>0){
      if ( cnaEmit->connect && evt < last_evt){//connect to the right... 
	pj = cnaEmit->pjump[sample][last_idx];
	last_idx = idx;
	Clone::predict( prior, post, cnaEmit, pj, Trans);	
	gsl_vector_memcpy(mem,entry);
	Clone::combine_prior( prior, mem, nloci);
      }
      else{
	gsl_vector_memcpy(prior,entry);
	if (nloci > 1) Clone::scale_prior(prior,nloci);
      }
    }
    else{//nClones==0
      gsl_vector_set_all( prior, cnaEmit->log_space ? 0.0 : 1.0);
    }
    //***GET POSTERIOR***
    gsl_vector_memcpy( mem, prior);
    alph = gsl_matrix_row( alpha_cna[sample], evt);
    if (cnaEmit->log_space){//multiply with forward posterior...
      gsl_vector_add( mem, &alph.vector);
      norm = log_vector_norm(mem);
      if (norm != norm) abort();
      gsl_vector_add_constant( mem, -norm);
      for (int l=0; l<nLevels; l++) mem->data[l] = exp(mem->data[l]);
    }
    else{
      gsl_vector_mul( mem, &alph.vector);
      norm = gsl_blas_dasum(mem);
      if (norm <= 0.0 || norm != norm) abort();
      gsl_vector_scale(mem,1.0/norm);
    }//multiply done.
    gsl_matrix_set_row( gamma_cna[sample], evt, mem);
    if (get_gofs){
      Clone::get_cna_gof(mem,sample,evt);
      for (int t=0;t<nTimes;t++) cna_gofs[t] += cna_all_gofs[t][sample][evt];
    }
    if (ent >= 0.0) ent += double(nloci) * Clone::entropy(mem);
    //***UPDATE STEP*** (normalization term not needed here)
    Clone::update( prior, post, cnaEmit, sample, evt, llhs);
  }
  // cleanup    
  gsl_vector_free(entry);
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(mem);
  if (Trans!= NULL) gsl_matrix_free(Trans);
}




//***BAF FWD-BWD*********************************************************
void Clone::do_baf_Fwd( int sample, double& llh, double*& llhs){
  int cnaSample = 0;
  int bafChr = bafEmit->chr[sample];
  if (!cnaEmit->is_set) abort();
  cnaSample = cnaEmit->idx_of[bafChr];
  if (nClones > 0){
    if (gamma_cna == NULL) abort();
    if (gamma_cna[cnaSample] == NULL) abort();
  }
  gsl_vector * prior = gsl_vector_alloc(nLevels);
  gsl_vector * post  = gsl_vector_alloc(nLevels);
  gsl_vector * mem   = gsl_vector_alloc(nLevels);
  gsl_vector * Prior = gsl_vector_alloc(nLevels);
  gsl_vector * flat = gsl_vector_alloc(nLevels);
  gsl_vector_set_all( flat, 1.0/double(nLevels));
  if (nClones>0) Clone::apply_maxtcn_mask( flat, bafChr, bafEmit->log_space);
  gsl_vector_memcpy(prior,flat);
  int idx=0, cna_evt=0, last_cna_evt =-1, nidx=0, nLoci=1;
  gsl_vector_view cna_post;
  double pj=0.0, norm  = 0.0;
  int last_evt = bafEmit->nEvents[sample]-1;
  llh = 0.0;
  for (int evt=0; evt <= last_evt; evt++){
    idx   = bafEmit->idx_of_event[sample][evt];
    nidx  = (evt < last_evt) ? bafEmit->idx_of_event[sample][evt+1] : bafEmit->nSites[sample];
    nLoci = nidx-idx;
    //***PREDICT STEP***
    if (nClones > 0){    
      if (bafEmit->connect && evt > 0){//connect to the left
	pj = bafEmit->pjump[sample][idx];
	Clone::predict( prior, post, bafEmit, pj, flat);
      }
      //***CONSISTENCY WITH CNA***
      cna_evt = bafEmit->Event_of_idx[sample][idx];
      if (cna_evt != last_cna_evt){//new BAF prior from CNA post
	cna_post = gsl_matrix_row( gamma_cna[cnaSample], cna_evt);
	get_baf_prior_from_cna_post( Prior, &cna_post.vector);
	last_cna_evt = cna_evt;
      }     
      if (bafEmit->connect){
	gsl_vector_memcpy( mem, Prior);
	Clone::combine_prior( prior, mem, nLoci);
      }
      else{
	gsl_vector_memcpy( prior, Prior);
	if (nLoci > 1) Clone::scale_prior( prior, nLoci);
      }
    }
    else{//nClones == 0
      gsl_vector_set_all( prior, bafEmit->log_space ? 0.0 : 1.0);
    }
    //***UPDATE STEP***
    norm = Clone::update( prior, post, bafEmit, sample, evt, llhs);
    llh += norm;
    if (save_baf_alpha == 1){
      gsl_matrix_set_row( alpha_baf[sample], evt, post);
    }
  }
  // cleanup    
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(mem);
  gsl_vector_free(flat);
  gsl_vector_free(Prior);
}


void Clone::do_baf_Bwd( int sample, double& ent){
  if (alpha_baf[sample] == NULL || gamma_baf[sample] == NULL) abort();
  int cnaSample = 0;
  int bafChr = bafEmit->chr[sample];
  if (!cnaEmit->is_set) abort();
  cnaSample = cnaEmit->idx_of[bafChr];
  if ( nClones>0 ){
    if ( gamma_cna == NULL ) abort();
    if ( gamma_cna[cnaSample] == NULL) abort();
  }
  gsl_vector * prior = gsl_vector_alloc(nLevels);
  gsl_vector * Prior = gsl_vector_alloc(nLevels);
  gsl_vector * post  = gsl_vector_alloc(nLevels);
  gsl_vector * mem   = gsl_vector_alloc(nLevels);
  gsl_vector * flat = gsl_vector_alloc(nLevels);
  gsl_vector_set_all(flat,1.0/double(nLevels));
  if (nClones>0) Clone::apply_maxtcn_mask( flat, bafChr, bafEmit->log_space);
  gsl_vector_memcpy(prior,flat);
  gsl_vector_view cna_post;
  double * llhs = NULL;
  gsl_vector_view alph;
  double pj = 0.0, norm=0;
  int last_evt = bafEmit->nEvents[sample] - 1;
  int last_idx = bafEmit->idx_of_event[sample][last_evt];
  int idx = 0, nidx=0, nLoci=1, cna_evt=0, last_cna_evt=-1;
  for (int evt = last_evt; evt >= 0; evt--){
    idx   = bafEmit->idx_of_event[sample][evt];
    nidx  = (evt < last_evt) ? bafEmit->idx_of_event[sample][evt+1] : bafEmit->nSites[sample];
    nLoci = nidx-idx;
    //***PREDICTION STEP***
    if (nClones > 0 ){ 
      if (bafEmit->connect && evt < last_evt){//connect to the right...
	pj = bafEmit->pjump[sample][last_idx];
	Clone::predict( prior, post, bafEmit, pj, flat);
	last_idx = idx;
      }
      //***CONSISTENCY WITH CNA***
      cna_evt = bafEmit->Event_of_idx[sample][idx];
      if (cna_evt != last_cna_evt){
	cna_post = gsl_matrix_row( gamma_cna[cnaSample], cna_evt);
	get_baf_prior_from_cna_post( Prior, &cna_post.vector);
	last_cna_evt = cna_evt;
      }      
      if (bafEmit->connect){
	gsl_vector_memcpy( mem, Prior);
	Clone::combine_prior( prior, mem, nLoci);
      }
      else{
	gsl_vector_memcpy( prior, Prior);
	if (nLoci > 1) Clone::scale_prior( prior, nLoci);
      }
    }
    else{//nClones==0
      gsl_vector_set_all( prior, bafEmit->log_space ? 0.0 : 1.0);
    }
    //***GET POSTERIOR***
    gsl_vector_memcpy( mem, prior);
    alph = gsl_matrix_row( alpha_baf[sample], evt);
    if (bafEmit->log_space){//multiply with fwd-posterior
      gsl_vector_add( mem, &alph.vector);
      norm = log_vector_norm(mem);
      if (norm != norm) abort();
      gsl_vector_add_constant(mem,-norm);
      for (int l=0; l<nLevels; l++) mem->data[l] = exp(mem->data[l]);
    }
    else{
      gsl_vector_mul( mem, &alph.vector);
      norm = gsl_blas_dasum(mem);
      if (norm <= 0.0 || norm != norm) abort();
      gsl_vector_scale( mem, 1.0/norm);
    }//multiply done
    //if (nClones > 0 && symmetrize_baf) Clone::sym_baf( mem, &cna_post.vector);
    gsl_matrix_set_row( gamma_baf[sample], evt, mem);
    if (get_gofs){
      Clone::get_baf_gof(mem,sample,evt);
      for (int t=0;t<nTimes;t++) baf_gofs[t] += baf_all_gofs[t][sample][evt];
    }
    if (ent >= 0) ent += double(nLoci)*Clone::entropy(mem);
    //***UPDATE STEP*** (normalization term not needed here)
    Clone::update( prior, post, bafEmit, sample, evt, llhs);
  }
  // cleanup    
  gsl_vector_free(prior);
  gsl_vector_free(Prior);
  gsl_vector_free(post);
  gsl_vector_free(mem);
  gsl_vector_free(flat);
}




//***SNV FWD-BWD***************************************************
void Clone::do_snv_Fwd(int sample, double& llh, double*& llhs){
  int snvChr = snvEmit->chr[sample];
  int cnaSample=-1, bafSample=-1;
  if (cnaEmit->is_set){
    if (cnaEmit->chrs.count(snvChr) == 0) abort();
    cnaSample = cnaEmit->idx_of[snvChr];
    if (bafEmit->is_set && bafEmit->chrs.count(snvChr) == 1){
      bafSample = bafEmit->idx_of[snvChr];
    }
  }
  if ( nClones > 0 && cnaEmit->is_set){//need CNA post
    if( gamma_cna == NULL || gamma_cna[cnaSample] == NULL) abort();
    if (bafEmit->is_set && bafSample >=0){//also need BAF post
      if( gamma_baf == NULL || gamma_baf[bafSample] == NULL) abort();
    }
  }
  gsl_vector * Prior = gsl_vector_alloc(nLevels);
  gsl_vector * prior = gsl_vector_alloc(nLevels);
  gsl_vector * post  = gsl_vector_alloc(nLevels);
  gsl_vector * mem   = gsl_vector_alloc(nLevels);
  gsl_matrix * Trans = NULL;
  if ( nClones>0 && snvEmit->connect ){
    Trans = gsl_matrix_alloc( nLevels, nLevels);
    gsl_matrix_memcpy(Trans,TransMat_snv[sample]);
  }
  if ( nClones > 0 && !cnaEmit->is_set && !snvEmit->connect && snvEmit->av_cn==NULL){
    gsl_vector_memcpy( Prior, snv_prior[snvChr]);
  }
  gsl_vector_set_all(prior,1.0/double(nLevels));
  if (nClones>0) Clone::apply_maxtcn_mask( prior, snvChr, snvEmit->log_space);
  gsl_vector_view cna_post, baf_post;
  double norm = 0.0, pj=0.0;
  int cna_evt=-1, last_cna_evt=-1;
  int baf_evt=-1, last_baf_evt=-1, baf_idx=0;
  int idx=0, nidx=0, nLoci=1, last_evt=(snvEmit->nEvents[sample]-1);
  llh = 0.0;
  for (int evt=0; evt <= last_evt; evt++){
    idx  = snvEmit->idx_of_event[sample][evt];
    nidx = (evt < last_evt) ? snvEmit->idx_of_event[sample][evt+1] : snvEmit->nSites[sample];
    nLoci = nidx-idx;
    //***PREDICT STEP***
    if (nClones > 0){
      if (snvEmit->connect && evt > 0){//connect to the left...
	pj = snvEmit->pjump[sample][idx];
	Clone::predict( prior, post, snvEmit, pj, Trans);
      }
      //***CONSISTENCY WITH CNA AND BAF***
      if (cnaEmit->is_set){//connect to CNA+BAF
	if (bafEmit->is_set && bafSample >= 0){//use CNA+BAF post
	  baf_evt = snvEmit->Event_of_idx[sample][idx];
	  baf_idx = bafEmit->idx_of_event[bafSample][baf_evt];
	  cna_evt = bafEmit->Event_of_idx[bafSample][baf_idx];
	}
	else{//use CNA post only
	  cna_evt = snvEmit->Event_of_idx[sample][idx];
	}
	if (cna_evt != last_cna_evt || baf_evt != last_baf_evt){//new segment, new prior
	  cna_post = gsl_matrix_row( gamma_cna[cnaSample], cna_evt);
	  if (bafEmit->is_set && bafSample >= 0){
	    baf_post = gsl_matrix_row( gamma_baf[bafSample], baf_evt);
	    get_snv_prior_from_cna_baf_post( Prior, &cna_post.vector, &baf_post.vector);
	  }
	  else{
	    get_snv_prior_from_cna_post( Prior, &cna_post.vector);
	  }
	  last_cna_evt = cna_evt;
	  last_baf_evt = baf_evt;
	}	
	if (snvEmit->connect){
	  gsl_vector_memcpy( mem, Prior);
	  Clone::combine_prior( prior, mem, nLoci);
	}
	else{
	  gsl_vector_memcpy( prior, Prior);
	  if (nidx-idx > 1) Clone::scale_prior( prior, nLoci);
	}
      }
      else if ( !snvEmit->connect ){
	if (snvEmit->av_cn != NULL && (evt == 0 || snvEmit->mean_tcn[0][sample][evt-1] != snvEmit->mean_tcn[0][sample][evt])){
	  Clone::get_snv_prior_from_av_cn( Prior, sample, evt);
	}
	gsl_vector_memcpy( prior, Prior);
      }
    }
    else{//nClones == 0
      gsl_vector_set_all( prior, snvEmit->log_space ? 0.0 : 1.0);
    }
    //***UPDATE STEP***
    norm = Clone::update( prior, post, snvEmit, sample, evt, llhs);
    llh += norm;
    if (save_snv_alpha == 1) gsl_matrix_set_row( alpha_snv[sample], evt, post);
  }
  // cleanup    
  gsl_vector_free(mem);
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(Prior);
  if (Trans != NULL) gsl_matrix_free(Trans);
}


void Clone::do_snv_Bwd( int sample, double& ent){
  if (alpha_snv[sample] == NULL || gamma_snv[sample] == NULL) abort();
  int snvChr = snvEmit->chr[sample];
  int cnaSample = -1,bafSample = -1;
  if (cnaEmit->is_set){
    if (cnaEmit->chrs.count(snvChr) == 0) abort();
    cnaSample = cnaEmit->idx_of[snvChr];
    if (bafEmit->is_set && bafEmit->chrs.count(snvChr) == 1){
      bafSample = bafEmit->idx_of[snvChr];
    }
  }
  if ( nClones >0 && cnaEmit->is_set){
    if( gamma_cna == NULL || gamma_cna[cnaSample] == NULL) abort();
    if (bafEmit->is_set && bafSample >= 0){
      if( gamma_baf == NULL || gamma_baf[bafSample] == NULL) abort();
    }
  }
  gsl_vector * Prior = gsl_vector_alloc(nLevels);
  gsl_vector * prior    = gsl_vector_alloc(nLevels);
  gsl_vector * post     = gsl_vector_alloc(nLevels);
  gsl_vector * mem      = gsl_vector_alloc(nLevels);
  gsl_matrix * Trans = NULL;
  if ( nClones>0 && snvEmit->connect){
    Trans = gsl_matrix_alloc( nLevels, nLevels);
    gsl_matrix_memcpy(Trans,TransMat_snv[sample]);
  }
  if ( nClones > 0 && !cnaEmit->is_set && !snvEmit->connect && snvEmit->av_cn==NULL){
    gsl_vector_memcpy( Prior, snv_prior[snvChr]);
  }
  gsl_vector_set_all(prior,1.0/double(nLevels));
  if (nClones>0) Clone::apply_maxtcn_mask( prior, snvChr, snvEmit->log_space);
  ent = 0.0;
  gsl_vector_view alph;
  gsl_vector_view cna_post,baf_post;
  double pj = 1.0, norm=0;
  double * llhs = NULL;
  int cna_evt=-1, last_cna_evt=-1;
  int baf_evt=-1, last_baf_evt=-1, baf_idx=0;
  int idx=0, nidx=0, nLoci=1;
  int last_evt = snvEmit->nEvents[sample]-1;
  int last_idx = snvEmit->idx_of_event[sample][last_evt];
  for (int evt = last_evt; evt >= 0 ; evt--){
    idx   = snvEmit->idx_of_event[sample][evt];
    nidx  = (evt < last_evt) ? snvEmit->idx_of_event[sample][evt+1] : snvEmit->nSites[sample];
    nLoci = nidx-idx;
    //***PREDICTION STEP***
    if (nClones > 0){
      if ( snvEmit->connect && evt < last_evt){//connect to the right...
	pj = snvEmit->pjump[sample][last_idx];
	last_idx = idx;
	Clone::predict( prior, post, snvEmit, pj, Trans);
      }
      if (cnaEmit->is_set){//connect to CNA...
	if (bafEmit->is_set && bafSample >= 0){//use BAF post
	  baf_evt = snvEmit->Event_of_idx[sample][idx];
	  baf_idx = bafEmit->idx_of_event[bafSample][baf_evt];
	  cna_evt = bafEmit->Event_of_idx[bafSample][baf_idx];
	}
	else{//use CNA post only
	  cna_evt = snvEmit->Event_of_idx[sample][idx];
	}
	if (cna_evt != last_cna_evt || baf_evt != last_baf_evt){
	  cna_post = gsl_matrix_row( gamma_cna[cnaSample], cna_evt);
	  if (bafEmit->is_set && bafSample >= 0){
	    baf_post = gsl_matrix_row( gamma_baf[bafSample], baf_evt);
	    get_snv_prior_from_cna_baf_post( Prior, &cna_post.vector, &baf_post.vector);
	  }
	  else{
	    get_snv_prior_from_cna_post( Prior, &cna_post.vector);
	  }
	  last_cna_evt = cna_evt;
	  last_baf_evt = baf_evt;
	}	
	if (snvEmit->connect){
	  gsl_vector_memcpy( mem, Prior);
	  Clone::combine_prior( prior, mem, nLoci);
	}
	else{
	  gsl_vector_memcpy( prior, Prior);
	  if (nidx-idx > 1) Clone::scale_prior( prior, nLoci);
	}
      }
      else if ( !snvEmit->connect ){
	if (snvEmit->av_cn != NULL && (evt == last_evt || snvEmit->mean_tcn[0][sample][evt+1] != snvEmit->mean_tcn[0][sample][evt])){
	  Clone::get_snv_prior_from_av_cn( Prior, sample, evt);
	}
	gsl_vector_memcpy( prior, Prior);
      }
    }  
    else{// nClones == 0
      gsl_vector_set_all( prior, snvEmit->log_space ? 0.0 : 1.0);
    }
    //***GET POSTERIOR***
    gsl_vector_memcpy( mem, prior);
    alph = gsl_matrix_row( alpha_snv[sample], evt);
    if (snvEmit->log_space){
      gsl_vector_add( mem, &alph.vector);
      norm = log_vector_norm(mem);
      if (norm != norm) abort();
      gsl_vector_add_constant(mem,-norm);
      for (int l=0; l<nLevels; l++) mem->data[l] = exp(mem->data[l]);
    }
    else{
      gsl_vector_mul( mem, &alph.vector);
      norm = gsl_blas_dasum(mem);
      if (norm <= 0.0 || norm != norm) abort();
      gsl_vector_scale( mem, 1.0/norm);
    }
    gsl_matrix_set_row( gamma_snv[sample], evt, mem);
    if (get_gofs){
      Clone::get_snv_gof(mem,sample,evt);
      for (int t=0;t<nTimes;t++) snv_gofs[t] += snv_all_gofs[t][sample][evt];
    }
    if (ent==0) ent += double(nLoci)*Clone::entropy(mem);
    //***UPDATE STEP*** (normalization term not needed here)
    Clone::update( prior, post, snvEmit, sample, evt, llhs);
  }
  // cleanup    
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(mem);
  gsl_vector_free(Prior);
  if (Trans!= NULL) gsl_matrix_free(Trans);
}


//*** GODNESS OF FITS **********************************************************
void Clone::allocate_all_gofs(){
  if (cnaEmit->is_set) cna_all_gofs = new double ** [nTimes];
  if (bafEmit->is_set) baf_all_gofs = new double ** [nTimes];
  if (snvEmit->is_set) snv_all_gofs = new double ** [nTimes];
  for (int t=0;t<nTimes;t++){
    if (cnaEmit->is_set){
      cna_all_gofs[t] = new double * [cnaEmit->nSamples];
      for (int s=0;s<cnaEmit->nSamples; s++)
	cna_all_gofs[t][s] = new double [cnaEmit->nEvents[s]];
    }
    if (bafEmit->is_set){
      baf_all_gofs[t] = new double * [bafEmit->nSamples];
      for (int s=0;s<bafEmit->nSamples; s++)
	baf_all_gofs[t][s] = new double [bafEmit->nEvents[s]];
    }
    if (snvEmit->is_set){
      snv_all_gofs[t] = new double * [snvEmit->nSamples];
      for (int s=0;s<snvEmit->nSamples; s++)
	snv_all_gofs[t][s] = new double [snvEmit->nEvents[s]];
    }
  }
}



void Clone::get_cna_gof(gsl_vector * post, int s, int evt){
  int cnaChr = cnaEmit->chr[s];
  for (int t=0; t<nTimes; t++) cna_all_gofs[t][s][evt] = 0;
  int first = cnaEmit->idx_of_event[s][evt];     
  int last 
    = (evt < cnaEmit->nEvents[s]-1) 
    ? cnaEmit->idx_of_event[s][evt+1] - 1 
    : cnaEmit->nSites[s]-1; 
  unsigned int n,N;
  double xobs,x,b,g;
  for (int idx=first; idx<=last; idx++){
    b = (cnaEmit->bias == NULL) ? 1.0 : cnaEmit->bias[s][idx];
    for (int t=0; t<nTimes; t++){
      n = cnaEmit->reads[t][s][idx];
      N = cnaEmit->depths[t][s][idx];
      if (N==0) continue;
      xobs = double(n) / double(N);
      g = 0;
      for (int l=0; l<nLevels; l++){
	x = tcn[cnaChr][t][l] * b * mass->data[t];
	g += post->data[l] * fabs(xobs-x);
      }
      cna_all_gofs[t][s][evt] += g;
    }
  }
}



void Clone::get_baf_gof(gsl_vector * post, int s, int evt){
  int bafChr = bafEmit->chr[s];
  for (int t=0; t<nTimes; t++) baf_all_gofs[t][s][evt] = 0;
  int first = bafEmit->idx_of_event[s][evt];     
  int last 
    = (evt < bafEmit->nEvents[s]-1) 
    ? bafEmit->idx_of_event[s][evt+1] - 1 
    : bafEmit->nSites[s]-1; 
  unsigned int n,N;
  double xobs,x,mntcn,g;
  for (int idx=first; idx<=last; idx++){
    for (int t=0; t<nTimes; t++){
      n = bafEmit->reads[t][s][idx];
      N = bafEmit->depths[t][s][idx];
      if (N==0) continue;
      xobs = double(n)/double(N);
      xobs = min(xobs,1.0-xobs);
      g = 0;
      mntcn = bafEmit->mean_tcn[t][s][evt];
      for (int l=0; l<nLevels; l++){
	x = tcn[bafChr][t][l] / mntcn;
	x = min(x,1.0-x);
	g += post->data[l] * fabs(xobs-x);
      }
      baf_all_gofs[t][s][evt] += g;
    }
  }
}

void Clone::get_snv_gof(gsl_vector * post, int s, int evt){
  int snvChr = snvEmit->chr[s];
  for (int t=0; t<nTimes; t++) snv_all_gofs[t][s][evt] = 0;
  int first = snvEmit->idx_of_event[s][evt];     
  int last 
    = (evt < snvEmit->nEvents[s]-1) 
    ? snvEmit->idx_of_event[s][evt+1] - 1 
    : snvEmit->nSites[s]-1; 
  unsigned int n,N;
  double xobs=0,x=0,y=0,b=0,mntcn=1,g=0;
  for (int idx=first; idx<=last; idx++){
    for (int t=0; t<nTimes; t++){
      n = snvEmit->reads[t][s][idx];
      N = snvEmit->depths[t][s][idx];
      if (N==0) continue;
      xobs = double(n)/double(N);
      mntcn 
	= (snvEmit->mean_tcn == NULL) 
	? tcn[snvChr][t][level_of[snvChr]] 
	: snvEmit->mean_tcn[t][s][evt];
      if (bulk_mean != NULL){
	b = (bulk_fix >= 0.0) ? bulk_fix : bulk_mean[t][s][idx];
	y =  b * (1.0-purity[t]) * double(normal_copy[snvChr]) / mntcn;
      }
      g=0;
      for (int l=0; l<nLevels; l++){
	x = y + clone_spectrum[t][l] / mntcn; 
	g += post->data[l] * fabs(xobs-x);
      }
      snv_all_gofs[t][s][evt] += g;
    }
  }
}



/*
void Clone::sym_baf( gsl_vector * bafPost, gsl_vector * cnaPost){
  if (map1==NULL && map2==NULL){
    map1 = new gsl_matrix * [maxtcn+1];
    for (int c=0; c<=maxtcn; c++){
      map1[c] = gsl_matrix_calloc(maxtcn+1,maxtcn+1);
      for (int i=0; i<=maxtcn; i++){
	if (i>c) break;
	for (int j=0; j<=maxtcn; j++){
	  if (j>c) break;
	  double val 
	    = (i==j && 2*j == c) 
	    ? 1.0 
	    : ( (i==j || c-i==j) ? 0.5 : 0.0);
	  gsl_matrix_set( map1[c], i, j, val);
	}
      } 
    }
    map2 = new gsl_matrix * [nLevels];
    for (int l=0; l<nLevels; l++){
      map2[l] = gsl_matrix_calloc(nLevels,nLevels);
      for (int i=0; i<nLevels; i++){
	for (int j=0; j<nLevels; j++){
	  double val = 1.0;
	  for (int n=0; n<nClones; n++){
	    val *= gsl_matrix_get( map1[copynumber[l][n]], copynumber[i][n], copynumber[j][n]);
	  }
	  gsl_matrix_set(map2[l],i,j,val);
	} 
      }
    }
  }//allocated
  gsl_vector * mem1 = gsl_vector_calloc(nLevels);
  gsl_vector * mem2 = gsl_vector_calloc(nLevels);
  for (int l=0; l<nLevels; l++){
    gsl_blas_dgemv(CblasNoTrans,1.0,map2[l],bafPost,0.0,mem2);
    double norm = gsl_blas_dasum(mem2);
    if (norm <= 0 && cnaPost->data[l] > 1.0e-4){
      abort();    
    }
    else if (norm > 0){
      gsl_vector_scale(mem2, cnaPost->data[l] / norm);
      gsl_vector_add(mem1, mem2);
    }
  }
  double norm = gsl_blas_dasum(mem1);
  if (norm<=0) abort();
  gsl_vector_scale(mem1,1.0/norm);
  gsl_vector_memcpy(bafPost,mem1);
  gsl_vector_free(mem1);
  gsl_vector_free(mem2);
}
*/
