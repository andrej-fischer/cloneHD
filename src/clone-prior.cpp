//clone-prior.cpp

//own headers...
#include "emission.h"
#include "log-space.h"
#include "clone.h"

using namespace std;


// CNA prior (only used for chr entry or w/o correlations)...
void Clone::set_cna_prior( gsl_vector * prior, int sample){
  double pncn = 0.5; //penalty for being  different from the normal copy number
  double pdif = 0.01;//penalty for having different copynumbers in clones
  if (cnaEmit->is_set == 0) abort();
  if (nClones==0){
    gsl_matrix_set_all(prior,1.0);
  }
  else{
    std::set<int> cns;
    int ncn = normal_copy[ cnaEmit->chr[sample] ];
    for (int i=0; i<nLevels; i++){
      cns.clear();
      for (int j=0; j<nClones; j++) cns.insert( copynumber[i][j] );
      double p = pow( pdif, (int) cns.size());
      for (int j=0; j<nClones; j++) p *= pow( pncn, abs(copynumber[i][j] - ncn) );
      gsl_vector_set( prior, i, p);
    } 
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
void Clone::initialize_cn_prior_snv(){// SNV prior, conditional on max-cn
  if (init_cn_prior_snv != NULL) gsl_matrix_free(init_cn_prior_snv);
  init_cn_prior_snv = gsl_matrix_calloc( maxcn+1, maxcn+1);
  gsl_matrix_set( init_cn_prior_snv, 0, 0, snv_fpr);//default:1.0e-4
  double p   = 0.5;// initial penalty for higher genotypes
  for (int cn=1; cn <= maxcn; cn++){
    if ( !maxcns.exists(cn) ) continue;
    gsl_vector_view subrow = gsl_matrix_subrow( init_cn_prior_snv, cn, 0, cn+1);
    gsl_vector_set( &subrow.vector, 0, p);//P0=P1>P2>... or P00=P10=P10=P11>P20...
    for (int i=1; i<=cn; i++) gsl_vector_set( &subrow.vector, i, pow( p, i));
    gsl_vector_scale( &subrow.vector, 1.0 / gsl_blas_dasum(&subrow.vector) );
  }
  //set above fixed priors
  if (cn_prior_snv != NULL) gsl_matrix_free(cn_prior_snv);
  cn_prior_snv = gsl_matrix_allocate(maxcn+1,nLevels);
  Clone::set_cn_prior_snv(init_cn_prior_snv);
}

//SNV-only mode, w/o cn-info...
void Clone::set_cn_prior_snv( gsl_matrix * snv_prior_by_maxcn){
  if (snvEmit->is_set == 0) abort();
  if ((int) cn_prior_snv->size2 != nLevels) abort();
  //cn==0
  gsl_matrix_set( cn_prior_snv, 0, 0, 1.0);
  double val;
  if (snvEmit->log_space){//log-transform?
    for (int l=0; l<nLevels; l++){
      val = gsl_matrix_get(cn_prior_snv,0,l);
      gsl_matrix_set( cn_prior_snv, 0, l, val>0.0 ? log(val) : logzero);
    }
  }
  //cn>0
  for (int cn=1; cn<=maxcn; cn++){//if total c.n. = cn...
    gsl_vector_view row = gsl_matrix_row(prior_per_clone,cn);
    if ( gsl_vector_max(&row.vector) <= 0.0) continue;
    gsl_vector * pr = gsl_vector_calloc(nLevels);
    for (int i=0; i<nLevels; i++){
      double p=1.0;
      for (int j=0; j<nClones; j++){
	if (copynumber[i][j] <= cn){
	  p *= gsl_matrix_get( prior_per_clone, cn, copynumber[i][j]);
	}
	else{
	  p = 0.0;
	  break;
	}
      }
      gsl_vector_set( pr, i, p);
    }
    // set all-zero-combination to SNV false positive rate...
    if ( !snvEmit->connect ){//CHECK
      double p0 = gsl_matrix_get( prior_per_clone, 0, 0);
      gsl_vector_set( pr, 0, p0);
    }
    //normalize
    double norm = gsl_blas_dasum(pr);
    if (norm <= 0.0 ) abort(); 
    gsl_vector_scale( pr, 1.0/norm);
    if (snvEmit->log_space){//log-transform
      for (int l=0; l<nLevels; l++){
	pr->data[l] = pr->data[l] > 0.0 ? log(pr->data[l]) : logzero;
      }
    }
    cn_prior_snv.insert( pair<int,gsl_vector*>(cn,pr));
  }
}

//CNA + BAF (+SNV) mode...
void Clone::set_baf_prior_map(){
  if ( baf_prior_map == NULL){
    baf_prior_map = gsl_matrix_alloc(maxcn+1,maxcn+1);
  }
  gsl_matrix_set_zero( baf_prior_map ); 
  double p = baf_pen;//penalty for complex chromosome status (default:1.0)
  double f = 0;
  for (int cn=0; cn <= maxcn; cn++){
    for (int bcn=0; bcn <= cn; bcn++){
      f = pow( p, int( fabs(double(bcn) - 0.5*double(cn)) ));
      gsl_matrix_set( baf_prior_map, bcn, cn, f);
    }
    //normalize...
    gsl_vector_view col = gsl_matrix_column( baf_prior_map,cn);
    double norm = gsl_blas_dasum(&col.vector);
    if (norm <= 0.0 || norm != norm) abort();
    gsl_vector_scale( &col.vector, 1.0 / norm);
  }
}

//CNA (+BAF) + SNV mode...
void Clone::set_snv_prior_map(){//either via BAF or else via CNA
  if (nClones == 0) abort();
  if (cnaEmit->is_set == 0) abort();
  //allocate
  if ( snv_prior_from_cna_baf_map == NULL){
    snv_prior_from_cna_baf_map = new gsl_matrix * [maxcn+1];
    for (int cn=0; cn <= maxcn; cn++){ 
      snv_prior_from_cna_baf_map[cn] = gsl_matrix_alloc( maxcn+1, maxcn+1);
    }
  }
  //via CNA + BAF posterior...
  double p = snv_pen;//penalty for SNVs in cn higher than max BAF cn (multiple hits)
  for (int cn=0; cn <= maxcn; cn++){
    gsl_matrix_set_zero( snv_prior_from_cna_baf_map[cn]); 
    for (int j=0; j<=cn; j++){
      for (int i=0; i<= cn; i++){
	double pen = pow( p, max( 0, i - max(j,cn-j)) );
	gsl_matrix_set( snv_prior_from_cna_baf_map[cn], i, j, pen);
      }
      //normalize...
      gsl_vector_view col = gsl_matrix_column( snv_prior_from_cna_baf_map[cn], j);
      double norm = gsl_blas_dasum(&col.vector);
      if (norm <=0 || norm != norm) abort();
      gsl_vector_scale( &col.vector, 1.0 / norm);
    }
  } 
  //via CNA posterior only...
  if ( snv_prior_from_cna_map == NULL){//allocate
    snv_prior_from_cna_map = gsl_matrix_alloc( maxcn+1, maxcn+1);
  }
  gsl_matrix_set_zero( snv_prior_from_cna_map);  
  p = (snvEmit->connect) ? 1.0 : snv_pen;// penalty for high genotypes 
  for (int cn=0; cn <= maxcn; cn++){
    for (int i=0; i <= cn; i++){
      gsl_matrix_set( snv_prior_from_cna_map, i, cn, pow(p,i));
    }
    //normalize...
    gsl_vector_view col = gsl_matrix_column( snv_prior_from_cna_map, cn);
    double norm = gsl_blas_dasum(&col.vector);
    if (norm <=0 || norm != norm) abort();
    gsl_vector_scale( &col.vector, 1.0 / norm);
  }
}




// CNA + BAF (+SNV) mode
void Clone::get_baf_prior_from_cna_post(gsl_vector * prior, gsl_vector * post){
  gsl_vector * post_per_clone  = gsl_vector_calloc( nClones*(maxcn+1) );
  gsl_matrix * prior_per_clone = gsl_matrix_calloc( nClones, maxcn+1);
  gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, post, 0.0, post_per_clone);
  gsl_vector_view po_pc,pr_pc;
  for (int i=0; i<nClones; i++){
    po_pc = gsl_vector_subvector( post_per_clone, i*(maxcn+1), maxcn+1);
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
      prior->data[l] = prior->data[l] > 0.0 ? log(prior->data[l]) : - 1.0e10;
    }
  }
  gsl_matrix_free(prior_per_clone);
  gsl_vector_free(post_per_clone);
}

// CNA + SNV mode...
void Clone::get_snv_prior_from_cna_post(gsl_vector * prior, gsl_vector * cnapost){
  gsl_vector * cnapostpc = gsl_vector_calloc( nClones*(maxcn+1));
  gsl_matrix * snv_prpc  = gsl_matrix_calloc( nClones, maxcn+1);
  gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, cnapost, 0.0, cnapostpc);
  gsl_vector_view cnappc,prpc;
  for (int i=0; i<nClones; i++){
    cnappc = gsl_vector_subvector( cnapostpc, i*(maxcn+1), maxcn+1);
    prpc = gsl_matrix_row( snv_prpc, i);
    gsl_blas_dgemv( CblasNoTrans, 1.0, snv_prior_from_cna_map, &cnappc.vector, 0.0, &prpc.vector);
  }
  Clone::apply_snv_prpc( prior, snv_prpc, cnapost->data[0]);
  gsl_matrix_free(snv_prpc);
  gsl_vector_free(cnapostpc);
}




//CNA + BAF + SNV mode...
void Clone::get_snv_prior_from_cna_baf_post(gsl_vector * prior, gsl_vector * cnapost, gsl_vector * bafpost){
  gsl_vector * cnapostpc  = gsl_vector_calloc( nClones*(maxcn+1));
  gsl_vector * bafpostpc  = gsl_vector_calloc( nClones*(maxcn+1));
  gsl_matrix * snv_prpc   = gsl_matrix_calloc( nClones, maxcn+1);
  gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, cnapost, 0.0, cnapostpc);
  gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, bafpost, 0.0, bafpostpc);
  gsl_vector_view prpc;
  gsl_vector * bafppc = gsl_vector_alloc(maxcn+1);
  for (int i=0; i<nClones; i++){
    //bafppc = gsl_vector_subvector( bafpostpc, i*(maxcn+1), maxcn+1);    
    prpc   = gsl_matrix_row( snv_prpc, i);
    for (int cn=0; cn<=maxcn; cn++){
      gsl_vector_set_zero(bafppc); 
      double norm=0;
      for (int j=0; j<=cn;j++){
	bafppc->data[j] = gsl_vector_get( bafpostpc, i*(maxcn+1)+j) + 1.0e-10;
	norm += bafppc->data[j];
      }
      gsl_vector_scale(bafppc,1.0/norm);
      double pcna = gsl_vector_get( cnapostpc, i*(maxcn+1) + cn);
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
void Clone::apply_snv_prpc( gsl_vector * prior, gsl_matrix * snv_prpc,  double pc0){
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
    if (norm > 0.0 ) gsl_vector_scale(prior, (1.0-snv_fpr)*(1.0-pc0) / norm);
    prior->data[0] = 1.0 - (1.0-snv_fpr)*(1.0-pc0);//SNV false positive rate
  }
  else{
    double norm = gsl_blas_dasum(prior);
    if (norm <=0.0 || norm!= norm) abort();
    gsl_vector_scale(prior, 1.0 / norm);
  }
  //log-transform?
  if (snvEmit->log_space){
    for (int l=0; l<nLevels; l++){
      prior->data[l] = prior->data[l] > 0.0 ? log(prior->data[l]) : - 1.0e10;
    }
  }
}




//get mean total copy number...
void Clone::get_phi(int sample){//only ever used for cnaEmit
  if (nClones == 0){
    double ncn = normal_copy[cnaEmit->chr[sample]];
    for (int t=0; t<nTimes; t++){
      for (int evt=0; evt < cnaEmit->nEvents[sample]; evt++){
	cnaEmit->phi[t][sample][evt] = double(ncn);
      }
    }
    for (int evt=0; evt<cnaEmit->nEvents[sample]; evt++){
      cnaEmit->cnmax[sample][evt] = ncn;
    }
  }
  else{
    if (gamma_cna[sample] == NULL) abort();
    double val=0;
    gsl_vector * post1 = gsl_vector_alloc(nLevels);
    gsl_vector * post2 = gsl_vector_alloc(maxcn+1);
    int ncn = normal_copy[ cnaEmit->chr[sample] ];
    for (int t=0; t<nTimes; t++){
      double nrml = (1.0-purity[t]) * double(ncn);
      for (int evt=0; evt < cnaEmit->nEvents[sample]; evt++){
	cnaEmit->phi[t][sample][evt] = 0.0;
	if (t==0){
	  gsl_vector_set_zero(post2);
	}
	for (int l=0; l<nLevels; l++){
	  val = gsl_matrix_get( gamma_cna[sample], evt, l);
	  post1->data[l] = val;
	  cnaEmit->phi[t][sample][evt] += val * (nrml + clone_spectrum[t][l]);
	  if (t==0){
	    int mx = *std::max_element( copynumber[l], copynumber[l] + nClones);
	    post2->data[mx] += val;
	  }
	}
	if (t==0){//get the maximum copynumber across clones...
	  cnaEmit->cnmax[sample][evt] = 0;
	  double p=0.0;
	  for (int cn=0; cn<=maxcn; cn++){
	    p += post2->data[cn];
	    if (p > 0.99){//conservative estimate of the highest copynumber (with pError < 1%)
	      cnaEmit->cnmax[sample][evt] = cn;
	      break;
	    }
	  }
	}
      }
    }
    gsl_vector_free(post1);
    gsl_vector_free(post2);
  }
}

void Clone::map_phi( Emission * fromEmit, int from_sample, Emission * toEmit){
  int fromChr = fromEmit->chr[from_sample];
  if ( toEmit->chrs.count(fromChr) == 0 ) abort();
  int sample = toEmit->idx_of[fromChr];
  if (nClones == 0){
    double ncn = normal_copy[toEmit->chr[sample]];
    for (int t=0; t<nTimes; t++){
      for (int evt=0; evt<toEmit->nEvents[sample]; evt++){
	toEmit->phi[t][sample][evt] = double(ncn);
      }
    }
    for (int evt=0; evt<toEmit->nEvents[sample]; evt++){
      toEmit->cnmax[sample][evt] = ncn;
    }
  }
  else{
    for (int evt=0; evt < toEmit->nEvents[sample]; evt++){
      int idx = toEmit->idx_of_event[sample][evt];
      int from_evt = toEmit->Event_of_idx[sample][idx];
      for (int t=0; t<nTimes; t++){//mean total copynumber
	toEmit->phi[t][sample][evt] = fromEmit->phi[t][from_sample][from_evt];
      }
      //max total copynumber
      toEmit->cnmax[sample][evt] = fromEmit->cnmax[from_sample][from_evt];
    }
  }
}



