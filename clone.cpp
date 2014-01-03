//clone.cpp


#define PI 3.1415926

//own headers...
#include "emission.h"
#include "common-functions.h"
#include "minimization.h"
#include "log-space.h"
#include "clone.h"

using namespace std;

// Constructor
Clone::Clone(){
  is_set=0;
  allocated = 0;
  nTimes    = 0;
  nClones   = 0;
  nLevels   = 0;
  //
  freqs     = NULL;
  purity    = NULL;
  marginals  = NULL;
  margin_map = NULL;
  baf_prior_map  = NULL;
  snp_prior_map  = NULL;
  maxcn     = 4;
  logn_set  = 0;
  nFreq     = 0;
  max_nFreq = 4;
  TransMat_cnv  = NULL;
  TransMat_snp  = NULL;
  mass      = NULL;
  log_mass  = NULL;
  nmean     = NULL;
  got_gamma = 0;
  mass_candidates = NULL;
  normal_copy      = NULL;
  copynumber       = NULL;
  copynumber_post  = NULL;
  init_cn_prior_snp= NULL;
  clone_spectrum   = NULL;
  cn2_post  = NULL;
  bulk_fix  = -1.0;
  baf_pen = 1.0;
  snp_pen = 0.01;
  snp_fpr = 1.0e-4;
  tcn    = NULL;//total copy number
  log_tcn= NULL;//log thereof
  bafEmitLog = NULL;
  cnvEmitLog = NULL;
  snpEmitLog = NULL;
  //posteriors
  alpha_cnv=NULL;
  gamma_cnv=NULL;
  alpha_baf=NULL;
  gamma_baf=NULL; 
  alpha_snp=NULL;
  gamma_snp=NULL;
  //
  cnv_total_llh=0;
  baf_total_llh=0;
  snp_total_llh=0;
  //bulk variables
  bulk_prior      = NULL;
  bulk_post       = NULL;
  bulk_dist       = NULL;
  bulk_prior_mean = NULL;
  bulk_post_mean  = NULL;
  bulk_mean       = NULL; 
  bulk_min = NULL;
}

void Clone::allocate( Emission * cnv, Emission * baf, Emission * snp){
  cnvEmit = cnv;
  bafEmit = baf;
  snpEmit = snp;
  if (cnvEmit->is_set) nTimes  = cnvEmit->nTimes;
  if (bafEmit->is_set) nTimes  = bafEmit->nTimes;
  if (snpEmit->is_set) nTimes  = snpEmit->nTimes;
  total_loci = 0;
  if (cnvEmit->is_set) total_loci += cnvEmit->total_loci;
  if (bafEmit->is_set) total_loci += bafEmit->total_loci;
  if (snpEmit->is_set) total_loci += snpEmit->total_loci;
  //the normal copy number of human DNA
  normal_copy = new int [25];
  normal_copy[0] = 0;
  for (int c=1; c<=23; c++){
    normal_copy[c] = 2;
  }
  normal_copy[24] = 0;//female by default
  //check whether there is a Y chromosome...
  int male=0;
  if (cnvEmit->is_set && cnvEmit->maxchr >= 24 && cnvEmit->idx_of[24] >= 0) male=1;
  if (bafEmit->is_set && bafEmit->maxchr >= 24 && bafEmit->idx_of[24] >= 0) male=1;
  if (snpEmit->is_set && snpEmit->maxchr >= 24 && snpEmit->idx_of[24] >= 0) male=1;
  if (male){
    normal_copy[23] = 1;
    normal_copy[24] = 1;
  }
  if (cnvEmit->is_set){
    mass     = gsl_vector_alloc(nTimes);
    log_mass = gsl_vector_alloc(nTimes);
    nmean    = gsl_vector_alloc(nTimes);
    Clone::set_logn();
    Clone::get_nmean();
    if (bafEmit->is_set) bafEmit->map_idx_to_Event(cnvEmit);
    if (snpEmit->is_set){
      if (bafEmit->is_set){
	snpEmit->map_idx_to_Event(bafEmit);
      }
      else{
	snpEmit->map_idx_to_Event(cnvEmit);
      }
    }
  }
  //
  purity = new double [nTimes];
  clone_spectrum = new double * [nTimes];
  for (int t=0; t<nTimes; t++) clone_spectrum[t] = NULL;
  min_purity = gsl_vector_calloc(nTimes);
  if (cnvEmit->is_set){
    tcn     = new double ** [nTimes];
    log_tcn = new double ** [nTimes];
    for (int t=0; t<nTimes; t++){
      tcn[t]       = new double * [cnvEmit->nSamples];
      log_tcn[t]   = new double * [cnvEmit->nSamples];
      for (int s=0; s<cnvEmit->nSamples; s++){
	tcn[t][s]     = NULL;
	log_tcn[t][s] = NULL;
      }
    }
  }
  allocated = 1;//done
}



void Clone::allocate_bulk_mean(){//mean only...
  if (snpEmit->is_set==0) abort();
  //allocate prior...
  bulk_prior_mean = new double * [snpEmit->nSamples];
  for (int s=0; s<snpEmit->nSamples; s++){
    bulk_prior_mean[s] = new double [ snpEmit->nSites[s] ];
  }
  //allocate posterior...
  bulk_post_mean  = new double ** [nTimes];
  for (int t=0; t<nTimes; t++){
    bulk_post_mean[t] = new double * [snpEmit->nSamples];   
    for (int s=0; s<snpEmit->nSamples; s++){
      bulk_post_mean[t][s] = new double [snpEmit->nSites[s]];
    }
  }
  //set pointers to prior...
  bulk_mean = new double ** [nTimes];
  for (int t=0; t<nTimes; t++){
    bulk_mean[t] = new double * [snpEmit->nSamples];   
    for (int s=0; s<snpEmit->nSamples; s++){
      bulk_mean[t][s] = bulk_prior_mean[s];
    }
  }
}


void Clone::allocate_bulk_dist(){//distribution and mean...
  if (snpEmit->is_set==0) abort();
  //allocate prior...
  bulk_prior      = new gsl_matrix * [snpEmit->nSamples];
  bulk_prior_mean = new double * [snpEmit->nSamples];
  for (int s=0; s<snpEmit->nSamples; s++){
    bulk_prior[s] = gsl_matrix_calloc( snpEmit->nSites[s], snpEmit->gridSize+1);
    bulk_prior_mean[s] = new double [ snpEmit->nSites[s] ];
  }
  //allocate posterior...
  bulk_post       = new gsl_matrix ** [nTimes];
  bulk_post_mean  = new double ** [nTimes];
  for (int t=0; t<nTimes; t++){
    bulk_post[t] = new gsl_matrix * [snpEmit->nSamples];   
    bulk_post_mean[t] = new double * [snpEmit->nSamples];   
    for (int s=0; s<snpEmit->nSamples; s++){
      bulk_post[t][s] = gsl_matrix_calloc( snpEmit->nSites[s], snpEmit->gridSize+1);
      bulk_post_mean[t][s] = new double [snpEmit->nSites[s]];
    }
  }
  //set pointers to prior...
  bulk_dist = new gsl_matrix ** [nTimes];
  bulk_mean = new double ** [nTimes];
  for (int t=0; t<nTimes; t++){
    bulk_dist[t] = new gsl_matrix * [snpEmit->nSamples];   
    bulk_mean[t] = new double * [snpEmit->nSamples];   
    for (int s=0; s<snpEmit->nSamples; s++){
      bulk_dist[t][s] = bulk_prior[s];
      bulk_mean[t][s] = bulk_prior_mean[s];
    }
  }
}

void Clone::set_bulk_to_post(){
  for (int t=0; t<nTimes; t++){
    for (int s=0; s<snpEmit->nSamples; s++){
      bulk_mean[t][s] = bulk_post_mean[t][s];
      if (bulk_post != NULL)  bulk_dist[t][s] = bulk_post[t][s];
    }
  }
}

void Clone::set_bulk_to_prior(){
  for (int t=0; t<nTimes; t++){
    for (int s=0; s<snpEmit->nSamples; s++){
      bulk_mean[t][s] = bulk_prior_mean[s];
      if (bulk_prior != NULL)  bulk_dist[t][s] = bulk_prior[s];
    }
  }
}

//Destructor
Clone::~Clone(){
  if (allocated == 1){
    delete [] purity;
    for (int l=0; l<nLevels; l++) delete [] copynumber[l]; 
    delete [] copynumber; 
    for (int t=0; t<nTimes; t++) delete [] clone_spectrum[t];
    delete [] clone_spectrum;
    delete [] normal_copy;
    if (cnvEmit->is_set){
      for (int t=0; t<nTimes; t++){
	for (int s=0; s<cnvEmit->nSamples; s++){
	  if ( tcn[t][s] != NULL) delete [] tcn[t][s];
	  if ( log_tcn[t][s] != NULL) delete [] log_tcn[t][s];
	}
	delete [] tcn[t];
	delete [] log_tcn[t];
      }
      delete [] tcn;
      delete [] log_tcn;
    }
    //
    if (mass!=NULL) gsl_vector_free(mass);
    if (log_mass!=NULL) gsl_vector_free(log_mass);
    if (nmean!=NULL) gsl_vector_free(nmean);
    if (min_purity!=NULL) gsl_vector_free(min_purity);
  }
  if (freqs != NULL) gsl_matrix_free(freqs);
  if (margin_map != NULL) gsl_matrix_free(margin_map);
  if (baf_prior_map != NULL) gsl_matrix_free(baf_prior_map);
  if (snp_prior_map != NULL){
    for (int cn=0; cn<=maxcn; cn++) gsl_matrix_free(snp_prior_map[cn]);
    delete [] snp_prior_map;
  }
  if (TransMat_cnv != NULL) gsl_matrix_free(TransMat_cnv);
  if (TransMat_snp != NULL) gsl_matrix_free(TransMat_snp);
  if (logn_set == 1){
    logn.clear();
    loggma.clear();
  }
}



// set frequency-dependent variables
void Clone::set(const gsl_matrix * freq){
  if( freq == NULL){//no clone
    nFreq   = 0;
    nClones = 0;
    freqs = NULL;
    if( nLevels == 0 || nClones != 0){//first time
      Clone::set_copynumbers();
      Clone::set_margin_map();
    }
  }// new number of clones?
  else if ( freq != NULL && (int) freq->size2 != nFreq ){
    nFreq = (int) freq->size2;
    nClones = nFreq;
    if (freqs != NULL) gsl_matrix_free(freqs);
    freqs = gsl_matrix_alloc(nTimes,nFreq);
    // set all the copynumbers per level
    Clone::set_copynumbers();
    Clone::set_margin_map();
    if (cnvEmit->is_set){
      Clone::set_TransMat_cnv();
    }
    if (bafEmit->is_set) 
      Clone::set_baf_prior_map();
    if (snpEmit->is_set){
      if (cnvEmit->is_set) Clone::set_snp_prior_map();
      Clone::set_TransMat_snp();
    }
  }
  if (freq != NULL){
    if (freqs == NULL) abort();
    gsl_matrix_memcpy( freqs, freq);
  }
  // set the clonal spectrum of possible levels
  Clone::set_clone_spectrum(freq); 
  for (int t=0; t<nTimes; t++){
    purity[t] = 0.0;
    if (freq != NULL){
      for (int f=0; f<nFreq;f++){
	purity[t] += gsl_matrix_get( freq, t, f);
      }
    }
  }
  is_set = 1;
}


void Clone::set_copynumbers(){
  //clear old...
  if (copynumber != NULL && nLevels != 0){
    for (int l=0; l<nLevels; l++){
      delete [] copynumber[l];
    }
    delete [] copynumber;
  }
  if (nClones == 0){// no clone
    copynumber       = new int * [1];
    copynumber[0]    = new int [1]; 
    copynumber[0][0] = 0;
    nLevels = 1;
  }
  else{//one or more clones
    int * offset = new int [nFreq]; 
    for (int f=0; f< nFreq; f++){
      offset[f]   = 0;
    }
    nLevels = 1;// the total number of levels
    for (int f=0;f<nFreq;f++){
      nLevels *= maxcn + 1;
    }
    int nl = nLevels;
    for (int f=0;f<nFreq;f++){
      nl = (int) double(nl) / double(maxcn + 1);
      offset[f] = nl;
    }
    // get the copynumbers for each level
    copynumber = new int * [nLevels];
    for (int l=0; l<nLevels; l++){
      copynumber[l] = new int [nFreq]; 
    }
    for (int f=0; f<nFreq; f++){
      for (int l=0; l<nLevels; l++){
	copynumber[l][f] = int(l/offset[f]) % (maxcn+1);
      }
    }
    delete [] offset;
  }
}



// set the spectrum of possible total clonal allele frequencies for diploids
void Clone::set_clone_spectrum(const gsl_matrix * freq){
  if (freq != NULL && nTimes != (int) freq->size1) abort();
  if (freq == NULL){//no clones
    nLevels = 1;
    for (int t=0; t<nTimes; t++){
      clone_spectrum[t]    = new double [1];
      clone_spectrum[t][0] = 0.0;
    }
  }
  else{//one or more clones
    double cf;
    for (int t=0; t<nTimes; t++){
      if (clone_spectrum[t] != NULL) delete [] clone_spectrum[t];
      clone_spectrum[t] = new double [nLevels];
      for (int l=0; l<nLevels; l++){      
	cf=0.0;
	for (int f=0;f<nFreq;f++){
	  cf += copynumber[l][f] * gsl_matrix_get(freq,t,f);
	}
	clone_spectrum[t][l] = cf;
      }
    }
  }
}


void Clone::get_complexity(){
  if (cnvEmit->is_set){
    complexity = double(nTimes)*( double(nClones) + 1.0) + pow( double(maxcn+1), nClones);
    double size = double(cnvEmit->total_loci);
    if (bafEmit->is_set) size += double(bafEmit->total_loci);
    if (snpEmit->is_set) size += double(snpEmit->total_loci);
    complexity *= log(size);
  }
  else if (bafEmit->is_set){
    complexity = double(nTimes)*(double(nClones) + 1.0) + pow(double(maxcn+1), nClones);
    double size = double(bafEmit->total_loci);
    complexity *= log(size);
  }
  else if (snpEmit->is_set){
    double * cnmax_freq = new double [maxcn+1];
    for (int cn=0; cn<=maxcn;cn++) cnmax_freq[cn]=0.0;
    double norm=0.0;
    if ( snpEmit->cnmax==NULL ){
      for (int s=0; s<snpEmit->nSamples; s++){
	cnmax_freq[ normal_copy[ snpEmit->chr[s] ] ] += double(snpEmit->nSites[s]);
	norm += double(snpEmit->nSites[s]);
      }  
    }
    else{
      for (int s=0; s<snpEmit->nSamples; s++){
	for (int evt=0; evt<snpEmit->nEvents[s];evt++){
	  cnmax_freq[snpEmit->cnmax[s][evt]] += 1.0;
	}
	norm += double(snpEmit->nEvents[s]);
      }
    }
    for (int cn=0; cn<=maxcn;cn++){
      cnmax_freq[cn] /= norm;
    }
    complexity = 0.0;
    for (int cn=1; cn<=maxcn;cn++){
      complexity += cnmax_freq[cn] * pow( double(cn+1), nClones);
      if (!snpEmit->connect && snpEmit->cnmax_seen.count(cn) == 1){
	complexity += double(cn+1);//d.o.f. of cn-prior
      }
    }
    if (!snpEmit->connect) complexity += 1.0;
    complexity += double(nTimes)*(double(nClones) + 1.0);
    complexity *= log(double(snpEmit->total_loci));
    delete [] cnmax_freq;
  }
}



void Clone::set_cn_prior_cnv( gsl_vector * prior, int sample){
  double pncn = 0.5; //penalty for being  different from the normal copy number
  double pdif = 0.01;//penalty for having different copynumbers in clones
  if (cnvEmit->is_set == 0) abort();
  std::set<int> cns;
  int ncn = normal_copy[ cnvEmit->chr[sample] ];
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
  if (cnvEmit->log_space){
    for (int l=0; l<nLevels;l++){ 
      prior->data[l] = prior->data[l]> 0.0 ? log(prior->data[l]) : -1.0e6;
    }
  }
}





void Clone::initialize_cn_prior_snp(){//SNP prior, conditional on max cn
  if (init_cn_prior_snp != NULL) gsl_matrix_free(init_cn_prior_snp);
  init_cn_prior_snp = gsl_matrix_calloc(maxcn+1,maxcn+1);
  gsl_matrix_set( init_cn_prior_snp, 0, 0, snp_fpr);
  int nprior = 1;
  double p   = 0.5;// penalty for higher genotypes
  for (int cn=1; cn <= maxcn; cn++){
    gsl_vector_view subrow = gsl_matrix_subrow( init_cn_prior_snp, cn, 0, cn+1);
    if (!snpEmit->connect ){
      if( snpEmit->cnmax_seen.count(cn) != 0){//un-correlated SNP data
	gsl_vector_set( &subrow.vector, 0, p);
	for (int i=1; i<=cn; i++) gsl_vector_set( &subrow.vector, i, pow( p, i));
	gsl_vector_scale( &subrow.vector, 1.0 / gsl_blas_dasum(&subrow.vector) );
      }
    }
    else{//flat prior for correlated SNP data (Transition matrix will be used)
      gsl_vector_set_all( &subrow.vector, 1.0/double(cn+1));
    }
    nprior++;
  }
  printf("Using these SNP copynumber priors per clone (combinations have multiplicative priors)\n");
  for (int i=0; i<(int) init_cn_prior_snp->size1; i++){
    for (int j=0; j<(int) init_cn_prior_snp->size2; j++){
      printf("%.3f ", gsl_matrix_get( init_cn_prior_snp, i, j));
    }
    cout<<endl;
  }
  //set above fixed priors
  Clone::set_cn_prior_snp(init_cn_prior_snp);
}



//SNP prior if only max cn is known...
void Clone::set_cn_prior_snp( gsl_matrix * prior_per_clone){
  if (snpEmit->is_set == 0) abort();
  // delete old ones...
  map<int,gsl_vector*>::iterator it;
  if ( !cn_prior_snp.empty() ){
    for (it=cn_prior_snp.begin(); it != cn_prior_snp.end(); ++it){
      gsl_vector_free(cn_prior_snp[it->first]);
    }
    cn_prior_snp.clear();
  }
  //set new ones...
  //cn==0
  gsl_vector * pr = gsl_vector_calloc(nLevels); 
  gsl_vector_set_all( pr, 0.0);
  gsl_vector_set( pr, 0, 1.0);
  if (snpEmit->log_space){//log-transform?
    for (int l=0; l<nLevels; l++){
      pr->data[l] = (pr->data[l] > 0.0) ? log(pr->data[l]) : - 1.0e6;
    }
  }
  cn_prior_snp.insert(pair<int,gsl_vector*>(0,pr));
  //cn>0
  for (int cn=1; cn<=maxcn; cn++){
    gsl_vector_view row = gsl_matrix_row(prior_per_clone,cn);
    if ( gsl_vector_max(&row.vector) <= 0.0) continue;
    gsl_vector * pr = gsl_vector_calloc(nLevels);
    gsl_vector_set_all( pr, 0.0);
    for (int i=1; i<nLevels; i++){
      double p=1.0;
      for (int j=0; j<nClones; j++){
	if (copynumber[i][j] <= cn){
	  p *= gsl_matrix_get( prior_per_clone, cn, copynumber[i][j]);//soft/hard?
	}
	else{
	  p = 0.0;
	  break;
	}
      }
      gsl_vector_set( pr, i, p);
    }
    // scale all-zero combination to sth rather small...
    if ( !snpEmit->connect ){//CHECK
      double p0 = gsl_matrix_get( prior_per_clone, 0, 0);
      gsl_vector_set( pr, 0, p0);
    }
    //normalize
    double norm = gsl_blas_dasum(pr);
    if (norm <= 0.0 ) abort(); 
    gsl_vector_scale( pr, 1.0/norm);
    if (snpEmit->log_space){//log-transform
      for (int l=0; l<nLevels; l++){
	pr->data[l] = pr->data[l] > 0.0 ? log(pr->data[l]) : - 1.0e6;
      }
    }
    cn_prior_snp.insert( pair<int,gsl_vector*>(cn,pr));
  }
}



void Clone::set_margin_map(){
  if (margin_map != NULL) gsl_matrix_free(margin_map);
  if (nClones == 0){
    margin_map = gsl_matrix_calloc( 1, 1);
    gsl_matrix_set_all(margin_map,1.0);
  }
  else{
    int ct = 0;
    for (int f=0; f<nFreq; f++) ct += maxcn+1;
    margin_map = gsl_matrix_calloc( ct, nLevels);
    ct=0;
    for (int f=0; f<nFreq; f++){
      for (int i=0; i<maxcn+1; i++){
	for (int j=0; j<nLevels; j++){
	  if (copynumber[j][f] == i){
	    gsl_matrix_set( margin_map, ct, j, 1.0);
	  }
	}
	ct++;
      }
    }
  }
}





void Clone::set_baf_prior_map(){
  if ( baf_prior_map == NULL){
    baf_prior_map = gsl_matrix_alloc(maxcn+1,maxcn+1);
  }
  gsl_matrix_set_zero( baf_prior_map ); 
  double p = baf_pen;//penalty for complex chromosome status
  double f = 0;
  for (int cn=0; cn <= maxcn; cn++){
    for (int bafcn=0; bafcn<=cn; bafcn++){
      f = pow( p, int( fabs(double(bafcn) - 0.5*double(cn)) ));
      gsl_matrix_set( baf_prior_map, bafcn, cn, f);
    }
    //normalize...
    gsl_vector_view col = gsl_matrix_column( baf_prior_map,cn);
    double norm = gsl_blas_dasum(&col.vector);
    if (norm <= 0.0 || norm != norm) abort();
    gsl_vector_scale( &col.vector, 1.0 / norm);
  }
}




void Clone::set_snp_prior_map(){//either via BAF or else via CNV
  if ( snp_prior_map == NULL){
    snp_prior_map = new gsl_matrix * [maxcn+1];
    for (int cn=0; cn <= maxcn; cn++){ 
      snp_prior_map[cn] = gsl_matrix_alloc( maxcn+1, maxcn+1);
    }
  }
  if (nClones == 0) abort();
  int cnf=0;
  if ( cnvEmit->is_set && bafEmit->is_set ){//via CNV and BAF posterior
    double p = snp_pen;//penalty for SNPs in cn higher than max BAF cn (multiple hits)
    for (int cn=0; cn <= maxcn; cn++){
      gsl_matrix_set_zero( snp_prior_map[cn]); 
      for (int j=0; j<= maxcn; j++){
	for (int i=0; i<= cn; i++){
	  double pen = pow( p, max( 0, i - max(j,cn-j)) );
	  gsl_matrix_set( snp_prior_map[cn], i, j, pen);
	}
      }
    } 
    cnf=maxcn;
  }
  else if (cnvEmit->is_set){//via CNV posterior only
    gsl_matrix_set_zero( snp_prior_map[0]);  
    double p = (snpEmit->connect) ? 1.0 : snp_pen;// penalty for high genotypes 
    for (int j=0; j<= maxcn; j++){
      for (int i=0; i<= j; i++){
	gsl_matrix_set( snp_prior_map[0], i, j, pow(p,i));
      }
    }
  }
  //normalize...
  for (int cn=0; cn <= cnf; cn++){
    for (int j=0; j<= maxcn; j++){
      gsl_vector_view col = gsl_matrix_column( snp_prior_map[cn], j);
      double norm = gsl_blas_dasum(&col.vector);
      if (norm <=0 || norm != norm) abort();
      gsl_vector_scale( &col.vector, 1.0 / norm);
    }
  }
}



void Clone::get_baf_prior_from_cnv_post(gsl_vector * prior, gsl_vector * post){
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



void Clone::get_snp_prior_from_cnv_post(gsl_vector * prior, gsl_vector * cnvpost){
  gsl_vector * cnvpostpc = gsl_vector_calloc( nClones*(maxcn+1));
  gsl_matrix * snp_prpc  = gsl_matrix_calloc( nClones, maxcn+1);
  gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, cnvpost, 0.0, cnvpostpc);
  gsl_vector_view cnvppc,prpc;
  for (int i=0; i<nClones; i++){
    cnvppc = gsl_vector_subvector( cnvpostpc, i*(maxcn+1), maxcn+1);
    prpc = gsl_matrix_row( snp_prpc, i);
    gsl_blas_dgemv( CblasNoTrans, 1.0, snp_prior_map[0], &cnvppc.vector, 0.0, &prpc.vector);
  }
  Clone::apply_snp_prpc(prior,snp_prpc);
  gsl_matrix_free(snp_prpc);
  gsl_vector_free(cnvpostpc);
}




//if both CNV and BAF are available to inform SNP...
void Clone::get_snp_prior_from_cnv_baf_post(gsl_vector * prior, gsl_vector * cnvpost, gsl_vector * bafpost){
  gsl_vector * cnvpostpc  = gsl_vector_calloc( nClones*(maxcn+1));
  gsl_vector * bafpostpc  = gsl_vector_calloc( nClones*(maxcn+1));
  gsl_matrix * snp_prpc   = gsl_matrix_calloc( nClones, maxcn+1);
  gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, cnvpost, 0.0, cnvpostpc);
  gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, bafpost, 0.0, bafpostpc);
  gsl_vector_view bafppc, prpc;
  for (int i=0; i<nClones; i++){
    bafppc = gsl_vector_subvector( bafpostpc, i*(maxcn+1), maxcn+1);   
    prpc   = gsl_matrix_row( snp_prpc, i);
    for (int cn=0; cn<=maxcn; cn++){
      double pcnv = gsl_vector_get( cnvpostpc, i*(maxcn+1) + cn);
      gsl_blas_dgemv( CblasNoTrans, pcnv, snp_prior_map[cn], &bafppc.vector, 1.0, &prpc.vector);
    }
  }
  Clone::apply_snp_prpc( prior, snp_prpc);
  gsl_matrix_free(snp_prpc);
  gsl_vector_free(cnvpostpc);
  gsl_vector_free(bafpostpc);
}
  


void Clone::apply_snp_prpc( gsl_vector * prior, gsl_matrix * snp_prpc){
  gsl_vector_set_all( prior, 1.0);
  int level=0;
  if (snpEmit->connect == 0){
    prior->data[0] = 1.0e-4;// 000... is not part of the somatic SNP data
    level = 1;
  }
  while (level<nLevels){
    for (int j=0; j<nClones; j++){
      prior->data[level] *= gsl_matrix_get( snp_prpc, j, copynumber[level][j]);
    }
    level++;
  }
  double norm = gsl_blas_dasum(prior);
  if (norm <=0.0 || norm!= norm) abort();
  gsl_vector_scale(prior,1.0/norm);
  if (snpEmit->log_space){//log-transform?
    for (int l=0; l<nLevels; l++){
      prior->data[l] = prior->data[l] > 0.0 ? log(prior->data[l]) : - 1.0e10;
    }
  }
}



void Clone::get_cnv_marginals(){
  if (cnvEmit->is_set == 0){
    cout<<"ERROR-1 in Clone::get_cnv_marginals()\n";
    exit(1);
  }
  if (marginals != NULL) gsl_matrix_free(marginals);
  int ct = 0;
  for (int f=0; f<nFreq; f++) ct += maxcn+1;
  marginals = gsl_matrix_calloc( cnvEmit->nSamples, ct);
  gsl_vector * post = gsl_vector_calloc(nLevels);
  gsl_vector * mem  = gsl_vector_calloc(nLevels);
  if (copynumber_post != NULL) gsl_matrix_free(copynumber_post);
  if (cn2_post != NULL) gsl_vector_free(cn2_post);
  copynumber_post = gsl_matrix_calloc( cnvEmit->nSamples, nLevels);
  cn2_post        = gsl_vector_calloc(nLevels);
  alpha_cnv = new gsl_matrix * [cnvEmit->nSamples];
  gamma_cnv = new gsl_matrix * [cnvEmit->nSamples];
  for (int s=0; s< cnvEmit->nSamples; s++){
    gsl_vector_view cn_row = gsl_matrix_row(copynumber_post,s);
    gsl_vector_view mg_row = gsl_matrix_row(marginals,s);
    alpha_cnv[s] = NULL;
    gamma_cnv[s] = NULL;
    //get the posterior here
    Clone::get_cnv_posterior(s);
    int ncn = normal_copy[ cnvEmit->chr[s] ];
    int evt=0,last_evt=-1;
    for (int idx=0; idx < cnvEmit->nSites[s]; idx++){
      evt = cnvEmit->event_of_idx[s][idx];
      if (evt != last_evt){
	gsl_matrix_get_row( post, gamma_cnv[s], evt);
	last_evt = evt;
      }
      if ( ncn == 2) gsl_vector_add( cn2_post, post);
      gsl_vector_add( &cn_row.vector, post);
      gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, post, 1.0, &mg_row.vector);
    }
    gsl_matrix_free(gamma_cnv[s]);
    gamma_cnv[s] = NULL;
    gsl_vector_add( mem, &cn_row.vector);
  }
  delete [] alpha_cnv;
  delete [] gamma_cnv;
  alpha_cnv = NULL;
  gamma_cnv = NULL;
  gsl_vector_free(post);
  //normalize cn2 posterior
  double norm = gsl_blas_dasum( cn2_post);
  gsl_vector_scale( cn2_post, 1.0/norm);
  //normalize posterior
  norm = gsl_blas_dasum(mem);
  for (int s=0; s<cnvEmit->nSamples; s++){
    gsl_vector_view cn_row = gsl_matrix_row(copynumber_post,s);
    gsl_vector_scale(&cn_row.vector,1.0/norm);
  }
  gsl_vector_free(mem);
  //normalize marginals
  gsl_vector_view part;
  ct = 0;
  for (int f=0; f<nClones; f++){
    norm = 0.0;
    for (int s=0; s<cnvEmit->nSamples; s++){
      part = gsl_matrix_subrow( marginals, s, ct, maxcn+1);
      norm += gsl_blas_dasum(&part.vector);
    }
    gsl_matrix_view pt = gsl_matrix_submatrix( marginals, 0, ct, cnvEmit->nSamples, maxcn+1);
    if (norm > 0.0){
      gsl_matrix_scale(&pt.matrix, 1.0 / norm);
    }
    else{
      cout<<"ERROR-2 in Clone::get_marginals().\n";
      exit(1);
    }
    ct += maxcn+1;
  }
}


void Clone::get_mass_candidates(){//only in normal_copy = 2 samples!
  Clone::get_cnv_marginals();
  levels_sorted.clear();
  for (int i=0; i<nLevels; i++) levels_sorted.push_back(i);
  SortDesc sd;
  sd.arg = cn2_post->data;
  std::sort( levels_sorted.begin(), levels_sorted.end(), sd);
  if (mass_candidates != NULL) gsl_matrix_free(mass_candidates);
  mass_candidates = gsl_matrix_calloc(nLevels,nTimes);
  double z;
  for (int i=0; i<nLevels; i++){
    int level = levels_sorted[i];
    for (int t=0; t<nTimes; t++){
      z = 2.0 * (1.0 - purity[t]) + clone_spectrum[t][level];
      z *= mass->data[t] / 2.0;
      gsl_matrix_set( mass_candidates, i, t, z);
    }
  }
}




void Clone::set_logn(){//only needed for cnvEmit
  for (int t=0; t<nTimes; t++){
    for (int s=0; s<cnvEmit->nSamples; s++){
      for (int l=0; l<cnvEmit->nSites[s]; l++){
	unsigned int N = cnvEmit->depths[t][s][l];
	if (N==0) continue;//no observation
	unsigned int n = cnvEmit->reads[t][s][l];
	if( logn.count(n) == 0 ){
	  logn[n]   = log(double(n));
	  loggma[n+1] = gsl_sf_lngamma(double(n+1));
	}
	if( logn.count(N) == 0 ){
	  logn[N]   = log(double(N));
	  loggma[N+1] = gsl_sf_lngamma(double(N+1));
	}
	if (cnvEmit->mode == 4){//for negative binomial model
	  int n1 = int(cnvEmit->shape*double(N));
	  if ( loggma.count(n1) == 0){
	    loggma[n1] = gsl_sf_lngamma(double(n1));
	  }
	  int n2 = n + n1;
	  if ( loggma.count(n2) == 0){
	    loggma[n2] = gsl_sf_lngamma(double(n2));
	  }
	}
      }
    }
  }
  for (int s=0; s<cnvEmit->nSamples; s++){
    int ncn = normal_copy[cnvEmit->chr[s]];
    if( logn.count(ncn) == 0 ){
      logn[ncn] = log(double(ncn));
    }
  }
  logn_set = 1;
}



void Clone::get_nmean(){//only needed for cnvEmit
  cnvEmit->minRate = 1.0e6;
  cnvEmit->maxRate = 0;
  for (int t=0; t<nTimes; t++){
    int ct = 0;
    double mn = 0;
    for (int s=0; s<cnvEmit->nSamples; s++){
      for (int l=0; l< cnvEmit->nSites[s]; l++){
	unsigned int n = cnvEmit->reads[t][s][l];
	unsigned int N = cnvEmit->depths[t][s][l];
	if (N==0) continue;
	if (N>0 && n>0){
	  cnvEmit->minRate = min(cnvEmit->minRate, double(n)/double(N));
	  cnvEmit->maxRate = max(cnvEmit->maxRate, double(n)/double(N));
	}
	mn += double(n) / double(N);
	ct++;
      }
    }
    //set mass to mean value
    nmean->data[t] = mn / double(ct);
  }
}


// only one clone can change its state,
// or clones in the same state can change in parallel
void Clone::set_TransMat_cnv(){
  if (TransMat_cnv != NULL) gsl_matrix_free(TransMat_cnv);
  TransMat_cnv = gsl_matrix_calloc(nLevels,nLevels);
  double norm;
  gsl_vector_view row;
  int jumps,cni,cnf;
  for (int i=0; i<nLevels; i++){
    for (int j=0; j<nLevels; j++){
      jumps=0;
      for(int k=0; k < nFreq; k++){
	if( copynumber[i][k] != copynumber[j][k]){
	  if ( jumps==0 ){
	    cni = copynumber[i][k];
	    cnf = copynumber[j][k];
	    jumps++;
	  }
	  else if (cni != copynumber[i][k] || cnf != copynumber[j][k]){
	    //else if (cni - cnf != copynumber[i][k] - copynumber[j][k]){
	    jumps++;
	  }
	}
      }
      if (jumps <= 1){
	gsl_matrix_set(TransMat_cnv,i, j, 1.0);
      }
      else{
	gsl_matrix_set(TransMat_cnv,i, j, 0.0);
      }
    }
    row  = gsl_matrix_row(TransMat_cnv,i);
    norm = gsl_blas_dasum(&row.vector);
    gsl_vector_scale(&row.vector,1.0/norm);
  }
}




// only one clone can change its state,
void Clone::set_TransMat_snp(){
  if (TransMat_snp != NULL) gsl_matrix_free(TransMat_snp);
  TransMat_snp = gsl_matrix_calloc(nLevels,nLevels);
  double norm;
  gsl_vector_view row;
  int jumps;
  for (int i=0; i<nLevels; i++){
    for (int j=0; j<nLevels; j++){
      jumps=0;
      for(int k=0; k < nFreq; k++){
	if( copynumber[i][k] != copynumber[j][k]){
	  jumps++;
	}
      }
      if (jumps <= 1){
	gsl_matrix_set(TransMat_snp,i, j, 1.0);
      }
      else{
	gsl_matrix_set(TransMat_snp,i, j, 0.0);
      }
    }
    row  = gsl_matrix_row(TransMat_snp,i);
    norm = gsl_blas_dasum(&row.vector);
    gsl_vector_scale(&row.vector,1.0/norm);
  }
}



void Clone::set_mass(gsl_vector * m){
  if ((int) m->size != nTimes) abort();
  gsl_vector_memcpy(mass,m);
  for (int t=0;t<nTimes;t++) log_mass->data[t] = log(mass->data[t]);
}


double Clone::get_all_total_llh(){
  if (bafEmit->is_set || snpEmit->is_set){
    alpha_cnv = new gsl_matrix * [cnvEmit->nSamples];
    gamma_cnv = new gsl_matrix * [cnvEmit->nSamples];
    for ( int s=0; s < cnvEmit->nSamples; s++){
      alpha_cnv[s] = gsl_matrix_calloc( cnvEmit->nEvents[s], nLevels);
      gamma_cnv[s] = gsl_matrix_calloc( cnvEmit->nEvents[s], nLevels);
    }
    save_cnv_alpha = 1;
    if (bafEmit->is_set && snpEmit->is_set){
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
    save_cnv_alpha = 0;
  }
  save_snp_alpha = 0;
  total_llh     = 0.0;
  cnv_total_llh = 0.0;
  baf_total_llh = 0.0;
  snp_total_llh = 0.0;
  int sample;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
  for ( sample=0; sample < cnvEmit->nSamples; sample++){
    double llh,ent;
    Clone::do_cnv_Fwd( sample, llh);
#pragma omp critical
    {
      cnv_total_llh += llh;
    }
    if (save_cnv_alpha==1){
      Clone::do_cnv_Bwd( sample, ent);
      Clone::get_phi(sample);
      if ( bafEmit->is_set ){
	Clone::map_phi( cnvEmit, sample, bafEmit);	
	if ( snpEmit->is_set ){
	  int bafsample = bafEmit->idx_of[cnvEmit->chr[sample]];
	  Clone::map_phi( bafEmit, bafsample, snpEmit);
	}
      }
      else{
	if ( snpEmit->is_set ) Clone::map_phi( cnvEmit, sample, snpEmit);
      }
      gsl_matrix_free(alpha_cnv[sample]);
      alpha_cnv[sample] = NULL;
    }
  }//END PARALLEL FOR
  //
  //BAF
  if ( bafEmit->is_set ){
#pragma omp parallel for schedule( dynamic, 1) default(shared)
    for ( sample=0; sample < bafEmit->nSamples; sample++){//START PARALLEL FOR
      double llh=0,ent=0;	
      Clone::do_baf_Fwd( sample, llh);
#pragma omp critical
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
  //SNP
  if( snpEmit->is_set ){
#pragma omp parallel for schedule( dynamic, 1) default(shared)
    for ( sample=0; sample < snpEmit->nSamples; sample++){//START PARALLEL FOR
      double llh = 0.0;
      Clone::do_snp_Fwd(sample, llh);
#pragma omp critical
      {
	snp_total_llh += llh;
      }
    }//END PARALLEL FOR
  }
  //cleanup...
  if (save_cnv_alpha==1){
    for ( sample=0; sample < cnvEmit->nSamples; sample++){
      gsl_matrix_free(gamma_cnv[sample]);
    }
    delete [] alpha_cnv;
    delete [] gamma_cnv;
    alpha_cnv = NULL;
    gamma_cnv = NULL;
  }
  if (save_baf_alpha==1){
    for ( sample=0; sample < bafEmit->nSamples; sample++){
      gsl_matrix_free(gamma_baf[sample]);
    }
    delete [] alpha_baf;
    delete [] gamma_baf;
    alpha_baf = NULL;
    gamma_baf = NULL;
  }
  total_llh = cnv_total_llh + baf_total_llh + snp_total_llh;
  return(total_llh);
}

double Clone::get_cnv_total_llh(){
  int sample;
  save_cnv_alpha = 0;
  cnv_total_llh  = 0.0;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
  for ( sample=0; sample < cnvEmit->nSamples; sample++){
    double llh;
    Clone::do_cnv_Fwd( sample, llh);
#pragma omp critical
    {
      cnv_total_llh += llh;
    }
  }//END PARALLEL FOR
  return(cnv_total_llh);
}


double Clone::get_baf_total_llh(){
  if ( nClones > 0 && cnvEmit->is_set && gamma_cnv == NULL ) abort();
  save_baf_alpha = 0;
  baf_total_llh = 0.0;
  int sample;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
  for ( sample=0; sample< bafEmit->nSamples; sample++){
    double llh;
    Clone::do_baf_Fwd( sample, llh);
#pragma omp critical
    {
      baf_total_llh += llh;
    }
  }
  return(baf_total_llh);
}


double Clone::get_snp_total_llh(){
  if ( nClones > 0 && cnvEmit->is_set && gamma_cnv == NULL ) abort();
  int sample;
  save_snp_alpha = 0;
  snp_total_llh  = 0.0;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
  for ( sample=0; sample< snpEmit->nSamples; sample++){
    double llh;
    Clone::do_snp_Fwd( sample, llh);
#pragma omp critical
    {
      snp_total_llh += llh;
    }
  }
  return(snp_total_llh);
}




double Clone::get_cnv_posterior(int sample){
  double llh=0, entropy=0;
  save_cnv_alpha = 1;
  //set fw-bw arrays    
  if (alpha_cnv[sample] == NULL) 
    alpha_cnv[sample] = gsl_matrix_calloc( cnvEmit->nEvents[sample], nLevels);
  if (gamma_cnv[sample] == NULL) 
    gamma_cnv[sample] = gsl_matrix_calloc( cnvEmit->nEvents[sample], nLevels);
  Clone::do_cnv_Fwd( sample, llh);
  Clone::do_cnv_Bwd( sample, entropy);
  gsl_matrix_free(alpha_cnv[sample]);
  alpha_cnv[sample] = NULL;
  return(llh);
}

double Clone::get_baf_posterior(int sample){
  if (cnvEmit->is_set){
    int cnv_sample = cnvEmit->idx_of[bafEmit->chr[sample]];
    if ( nClones > 0 && (gamma_cnv == NULL || gamma_cnv[cnv_sample] == NULL)){
      abort();
    }
  }
  double llh=0, entropy=0;
  save_baf_alpha = 1;
  //set fw-bw arrays    
  if (alpha_baf[sample] == NULL) 
    alpha_baf[sample] = gsl_matrix_calloc( bafEmit->nEvents[sample], nLevels);
  if (gamma_baf[sample] == NULL) 
    gamma_baf[sample] = gsl_matrix_calloc( bafEmit->nEvents[sample], nLevels);
  Clone::do_baf_Fwd( sample, llh);
  Clone::do_baf_Bwd( sample, entropy);
  gsl_matrix_free(alpha_baf[sample]);
  alpha_baf[sample] = NULL;
  return(llh);
}


double Clone::get_snp_posterior(int sample){
  if (cnvEmit->is_set){
    int cnv_sample = cnvEmit->idx_of[snpEmit->chr[sample]];
    if ( nClones > 0 && (gamma_cnv == NULL || gamma_cnv[cnv_sample] == NULL)){
      abort();
    }
  }
  double llh=0, ent=0;
  save_snp_alpha = 1;
  //set fw-bw arrays    
  if (alpha_snp[sample] == NULL) 
    alpha_snp[sample] = gsl_matrix_calloc( snpEmit->nEvents[sample], nLevels);
  if (gamma_snp[sample] == NULL) 
    gamma_snp[sample] = gsl_matrix_calloc( snpEmit->nEvents[sample], nLevels);
  Clone::do_snp_Fwd( sample, llh);
  Clone::do_snp_Bwd( sample, ent);
  gsl_matrix_free(alpha_snp[sample]);
  alpha_snp[sample] = NULL;
  return(llh);
}


void Clone::do_cnv_Fwd( int sample, double& llh){
  Clone::set_tcn(sample);//pre-compute total copy number
  gsl_vector * prior    = gsl_vector_alloc(nLevels);
  gsl_vector * post     = gsl_vector_alloc(nLevels);
  gsl_matrix * Trans = NULL;
  double flat = (cnvEmit->log_space) ? -log(double(nLevels)) : 1.0/double(nLevels);
  gsl_vector_set_all(prior, flat);
  gsl_vector_set_all(post,  flat);  
  if (nClones>0){//preparations...
    Clone::set_cn_prior_cnv( prior, sample);
    if (cnvEmit->connect){
      Trans = gsl_matrix_alloc( nLevels, nLevels);
      gsl_matrix_memcpy( Trans, TransMat_cnv);
    }
  }
  else{
    gsl_vector_set_all(prior,flat);
  }
  gsl_vector_memcpy( post, prior);
  double norm = 0.0;
  double pj = 1.0;
  int idx=0;
  llh = 0.0;
  for (int evt=0; evt < cnvEmit->nEvents[sample]; evt++){
    idx = cnvEmit->idx_of_event[sample][evt];
    //***PREDICT STEP***
    if (nClones > 0 && cnvEmit->connect && evt > 0){//connect to the left...
      pj = cnvEmit->pjump[sample][idx];
      Clone::predict( prior, post, cnvEmit, pj, Trans);
    }//connect done.
    else if (nClones == 0){
      gsl_vector_set_all( prior, flat);
    }
    //***UPDATE***
    norm = Clone::update( prior, post, cnvEmit, sample, evt);
    llh += norm;
    if (save_cnv_alpha == 1) gsl_matrix_set_row( alpha_cnv[sample], evt, post);
  }
  // cleanup    
  gsl_vector_free(prior);
  gsl_vector_free(post);
  if (Trans != NULL) gsl_matrix_free(Trans);
}



void Clone::do_baf_Fwd( int sample, double& llh){
  int cnv_sample = 0;
  if (cnvEmit->is_set){
    cnv_sample = cnvEmit->idx_of[bafEmit->chr[sample]];
    if (nClones > 0){
      if (gamma_cnv == NULL) abort();
      if (gamma_cnv[cnv_sample] == NULL) abort();
    }
  }
  gsl_vector * prior = gsl_vector_alloc(nLevels);
  gsl_vector * mem   = gsl_vector_alloc(nLevels);
  gsl_vector * post  = gsl_vector_alloc(nLevels);
  gsl_vector * cn_prior  = gsl_vector_alloc(nLevels);
  double flat = (bafEmit->log_space) ? -log(double(nLevels)) : 1.0/double(nLevels);
  gsl_vector_set_all( cn_prior, flat);
  gsl_vector_set_all( prior,    flat);
  gsl_vector_set_all( post,     flat);
  double norm  = 0.0;
  int idx=0, cnv_evt=0, last_cnv_evt =-1, nidx=0;
  gsl_vector_view cnv_post;
  double pj=0.0;
  int last_evt = bafEmit->nEvents[sample]-1;
  llh = 0.0;
  for (int evt=0; evt <= last_evt; evt++){
    idx = bafEmit->idx_of_event[sample][evt];
    //***PREDICT STEP***
    if (nClones > 0){    
      if (bafEmit->connect && evt > 0){//connect to the left...
	pj = bafEmit->pjump[sample][idx];
	Clone::predict( prior, post, bafEmit, pj, flat);
      }//...connect done
      if (cnvEmit->is_set){//connect with CNV	
	cnv_evt = bafEmit->Event_of_idx[sample][idx];
	if (cnv_evt != last_cnv_evt){
	  cnv_post = gsl_matrix_row( gamma_cnv[cnv_sample], cnv_evt);
	  get_baf_prior_from_cnv_post( cn_prior, &cnv_post.vector);
	  last_cnv_evt = cnv_evt;
	}
	if (bafEmit->connect) gsl_vector_memcpy( mem, prior);
	gsl_vector_memcpy( prior, cn_prior);
	nidx = (evt < last_evt) ? bafEmit->idx_of_event[sample][evt+1] : bafEmit->nSites[sample];
	if ( nidx-idx > 1){//exponentiate prior for all observations in this block
	  gsl_vector_scale( prior, double(nidx-idx));//log-space!
	  norm = log_vector_norm(prior);
	  gsl_vector_add_constant(prior,-norm);
	}  
	//multiply two priors and rescale...
	if (bafEmit->connect){
	  if (bafEmit->log_space){
	    gsl_vector_add(prior,mem);
	    norm = log_vector_norm(prior);
	    gsl_vector_add_constant(prior,-norm);
	  }
	  else{
	    gsl_vector_mul(prior,mem);
	    norm = gsl_blas_dasum(prior);
	    if (norm!=norm || norm < 0.0) abort();
	    gsl_vector_scale( prior, 1.0/norm);
	  }
	}
      }
    }
    else{//just one element!
      gsl_vector_set_all( prior, flat);
    }
    //***UPDATE STEP***
    norm = Clone::update( prior, post, bafEmit, sample, evt);
    llh += norm;
    if (save_baf_alpha == 1) gsl_matrix_set_row( alpha_baf[sample], evt, post);
  }
  // cleanup    
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(mem);
  gsl_vector_free(cn_prior);
}


void Clone::do_snp_Fwd(int sample, double& llh){
  int baf_sample = (bafEmit->is_set) ? bafEmit->idx_of[snpEmit->chr[sample]] : 0;
  int cnv_sample = (cnvEmit->is_set) ? cnvEmit->idx_of[snpEmit->chr[sample]] : 0;
  if ( nClones >0 ){
    if (cnvEmit->is_set)
      if( gamma_cnv == NULL || gamma_cnv[cnv_sample] == NULL) abort();
    if (bafEmit->is_set)
      if( gamma_baf == NULL || gamma_baf[baf_sample] == NULL) abort();
  }
  gsl_vector * cn_prior = gsl_vector_alloc(nLevels);
  gsl_vector * prior    = gsl_vector_alloc(nLevels);
  gsl_vector * post     = gsl_vector_alloc(nLevels);
  gsl_vector * mem      = gsl_vector_alloc(nLevels);
  double flat = (snpEmit->log_space) ? -log(double(nLevels)) : 1.0/double(nLevels);
  gsl_vector_set_all( cn_prior, flat);
  gsl_vector_set_all( prior,    flat);
  gsl_vector_set_all( post,     flat);
  gsl_matrix * Trans = NULL;
  gsl_vector_view cnv_post, baf_post;
  gsl_vector_memcpy( post, prior);
  int cn=0, ocn=-1;
  if ( nClones>0 && snpEmit->connect == 1){
    Trans = gsl_matrix_alloc( nLevels, nLevels);
    gsl_matrix_memcpy( Trans, TransMat_snp);
  }
  double norm = 0.0;
  double pj   = 1.0;
  int ncn     = normal_copy[ snpEmit->chr[sample] ];
  int cnv_evt=-1, last_cnv_evt=-1;
  int baf_evt=-1, last_baf_evt=-1, baf_idx=0;
  int idx=0, nidx=0, last_evt=snpEmit->nEvents[sample]-1;
  llh = 0.0;
  for (int evt=0; evt <= last_evt; evt++){
    idx = snpEmit->idx_of_event[sample][evt];
    //***PREDICT STEP***
    if (nClones > 0){
      if (snpEmit->connect && evt > 0){//connect to the left...
	pj = snpEmit->pjump[sample][idx];
	Clone::predict( prior, post, snpEmit, pj, Trans);
      }//...connect to left done
      if (cnvEmit->is_set){//connect to CNV...
	if (bafEmit->is_set){
	  baf_evt = snpEmit->Event_of_idx[sample][idx];
	  baf_idx = bafEmit->idx_of_event[baf_sample][baf_evt];
	  cnv_evt = bafEmit->Event_of_idx[baf_sample][baf_idx];
	}
	else{
	  cnv_evt = snpEmit->Event_of_idx[sample][idx];
	}
	if (cnv_evt != last_cnv_evt || baf_evt != last_baf_evt){//new segment, new prior
	  cnv_post = gsl_matrix_row( gamma_cnv[cnv_sample], cnv_evt);
	  if (bafEmit->is_set){
	    baf_post = gsl_matrix_row( gamma_baf[baf_sample], baf_evt);
	    get_snp_prior_from_cnv_baf_post( cn_prior, &cnv_post.vector, &baf_post.vector);
	  }
	  else{
	    get_snp_prior_from_cnv_post( cn_prior, &cnv_post.vector);
	  }
	  last_cnv_evt = cnv_evt;
	  last_baf_evt = baf_evt;
	}
	if (snpEmit->connect) gsl_vector_memcpy( mem, prior);
	gsl_vector_memcpy( prior, cn_prior);
	nidx = (evt < last_evt) ? snpEmit->idx_of_event[sample][evt+1] : snpEmit->nSites[sample];
	if ( nidx-idx > 1){//exponentiate prior for all observations in this block
	  gsl_vector_scale( prior, double(nidx-idx));//log-space!
	  norm = log_vector_norm(prior);
	  gsl_vector_add_constant(prior,-norm);
	}  
	// multiply the two priors and rescale...
	if (snpEmit->connect){
	  if (snpEmit->log_space){
	    gsl_vector_add(prior,mem);
	    norm = log_vector_norm(prior);
	    gsl_vector_add_constant(prior,-norm);
	  }
	  else{
	    gsl_vector_mul(prior,mem);
	    norm = gsl_blas_dasum(prior);
	    if (norm!=norm || norm <= 0.0) abort();
	    gsl_vector_scale(prior,1.0/norm);
	  }
	}
      }//connect to CNV done.
      else if (!cnvEmit->is_set && !snpEmit->connect){//cnv not set, no correlations...
	cn = (snpEmit->cnmax == NULL) ? ncn : snpEmit->cnmax[sample][evt];
	if (cn != ocn){
	  gsl_vector_memcpy( cn_prior, cn_prior_snp[cn]);
	  ocn = cn;
	}//or flat over all levels?
	gsl_vector_memcpy( prior, cn_prior);
      }
    }
    else if (nClones == 0){//just one element!
      gsl_vector_set_all( prior, flat);
    }
    //***UPDATE STEP***
    norm = Clone::update( prior, post, snpEmit, sample, evt);
    llh += norm;
    if (save_snp_alpha == 1) gsl_matrix_set_row( alpha_snp[sample], evt, post);
  }
  // cleanup    
  gsl_vector_free(mem);
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(cn_prior);
  if (Trans != NULL) gsl_matrix_free(Trans);
}



void Clone::do_cnv_Bwd(int sample, double& ent){
  if (alpha_cnv[sample] == NULL || gamma_cnv[sample] == NULL) abort();
  gsl_vector * prior    = gsl_vector_alloc(nLevels);
  gsl_vector * post     = gsl_vector_alloc(nLevels);
  gsl_vector * mem      = gsl_vector_alloc(nLevels);
  double flat = (cnvEmit->log_space) ? -log(double(nLevels)) : 1.0/double(nLevels);
  gsl_vector_set_all( prior, flat);
  gsl_vector_set_all( post,  flat);
  gsl_matrix * Trans = NULL;
  ent = 0.0;
  if (nClones>0){
    Clone::set_cn_prior_cnv( prior, sample);
    if (cnvEmit->connect){
      Trans = gsl_matrix_alloc(nLevels,nLevels);
      gsl_matrix_memcpy( Trans, TransMat_cnv);
    }
  }
  gsl_vector_memcpy( post, prior);
  gsl_vector_view alph;
  double pj = 1.0;
  double norm;
  int idx = 0;
  int last_evt = cnvEmit->nEvents[sample]-1;
  int last_idx = cnvEmit->idx_of_event[sample][last_evt];
  for (int evt = last_evt; evt >= 0; evt--){
    idx = cnvEmit->idx_of_event[sample][evt];
    //***PREDICTION STEP***
    if ( cnvEmit->connect && nClones > 0 &&  evt < last_evt){//connect to the right... 
      pj = cnvEmit->pjump[sample][last_idx];
      last_idx = idx;
      Clone::predict( prior, post, cnvEmit, pj, Trans);
    }//connect done.
    else{
      gsl_vector_set_all( prior, flat);
    }
    //***GET POSTERIOR***
    gsl_vector_memcpy( mem, prior);
    alph = gsl_matrix_row( alpha_cnv[sample], evt);
    if (cnvEmit->log_space){//multiply with forward posterior...
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
    gsl_matrix_set_row( gamma_cnv[sample], evt, mem);
    //ent += Clone::entropy(mem);
    //***UPDATE STEP*** (normalization term not needed here)
    Clone::update( prior, post, cnvEmit, sample, evt);
  }
  // cleanup    
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(mem);
  if (Trans!= NULL) gsl_matrix_free(Trans);
}



void Clone::do_baf_Bwd( int sample, double& ent){
  if (alpha_baf[sample] == NULL || gamma_baf[sample] == NULL) abort();
  int cnv_sample = 0;
  if (cnvEmit->is_set){
    cnv_sample = cnvEmit->idx_of[bafEmit->chr[sample]];
    if ( nClones>0 ){
      if ( gamma_cnv == NULL ) abort();
      if ( gamma_cnv[cnv_sample] == NULL) abort();
    }
  }
  gsl_vector * prior    = gsl_vector_alloc(nLevels);
  gsl_vector * cn_prior = gsl_vector_alloc(nLevels);
  gsl_vector * post     = gsl_vector_alloc(nLevels);
  gsl_vector * mem      = gsl_vector_alloc(nLevels);
  double flat = (bafEmit->log_space) ? -log(double(nLevels)) : 1.0/double(nLevels);
  gsl_vector_set_all(prior,flat);
  gsl_vector_set_all(cn_prior,flat);
  gsl_vector_set_all(post,flat);
  gsl_vector_view cnv_post;
  double norm;
  ent = 0.0;
  gsl_vector_view alph;
  double pj = 0.0;
  int last_evt = bafEmit->nEvents[sample] - 1;
  int last_idx = bafEmit->idx_of_event[sample][last_evt];
  int idx = 0,nidx=0;
  int cnv_evt=0,last_cnv_evt=-1;
  for (int evt = last_evt; evt >= 0; evt--){
    idx = bafEmit->idx_of_event[sample][evt];
    //***PREDICTION STEP***
    if (nClones > 0 ){ 
      if (bafEmit->connect && evt < last_evt){//connect to the right...
	pj = bafEmit->pjump[sample][last_idx];
	last_idx = idx;
	Clone::predict( prior, post, bafEmit, pj, flat);
      }//...connect done.   
      if (cnvEmit->is_set){//connect with CNV...	
	cnv_evt = bafEmit->Event_of_idx[sample][idx];
	if (cnv_evt != last_cnv_evt){
	  cnv_post = gsl_matrix_row( gamma_cnv[cnv_sample], cnv_evt);
	  get_baf_prior_from_cnv_post( cn_prior, &cnv_post.vector);
	  last_cnv_evt = cnv_evt;
	}
	if (bafEmit->connect) gsl_vector_memcpy( mem, prior);
	gsl_vector_memcpy( prior, cn_prior);
	nidx = (evt < last_evt) ? bafEmit->idx_of_event[sample][evt+1] : bafEmit->nSites[sample];
	if ( nidx-idx > 1){//exponentiate prior for all observations in this block
	  gsl_vector_scale( prior, double(nidx-idx));//log-space!
	  norm = log_vector_norm(prior);
	  gsl_vector_add_constant(prior,-norm);
	}  
	//multiply two priors and rescale...
	if (bafEmit->connect){
	  if (bafEmit->log_space){
	    gsl_vector_add(prior,mem);
	    norm = log_vector_norm(prior);
	    gsl_vector_add_constant(prior,-norm);
	  }
	  else{
	    gsl_vector_mul(prior,mem);
	    norm = gsl_blas_dasum(prior);
	    if (norm!=norm || norm <0.0) abort();
	    gsl_vector_scale(prior,1.0/norm);
	  }
	}
      }//...connect with CNV done.
    }
    else if (nClones == 0){//just one element!
      gsl_vector_set_all( prior, flat);
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
    gsl_matrix_set_row( gamma_baf[sample], evt, mem);
    //ent += Clone::entropy(mem);
    //***UPDATE STEP*** (normalization term not needed here)
    Clone::update( prior, post, bafEmit, sample, evt);
  }
  // cleanup    
  gsl_vector_free(prior);
  gsl_vector_free(cn_prior);
  gsl_vector_free(post);
  gsl_vector_free(mem);
}




void Clone::do_snp_Bwd( int sample, double& ent){
  if (alpha_snp[sample] == NULL || gamma_snp[sample] == NULL) abort();
  int baf_sample = (bafEmit->is_set) ? bafEmit->idx_of[snpEmit->chr[sample]] : 0;
  int cnv_sample = (cnvEmit->is_set) ? cnvEmit->idx_of[snpEmit->chr[sample]] : 0;
  if ( nClones >0 ){
    if (cnvEmit->is_set)
      if( gamma_cnv == NULL || gamma_cnv[cnv_sample] == NULL) abort();
    if (bafEmit->is_set)
      if( gamma_baf == NULL || gamma_baf[baf_sample] == NULL) abort();
  }
  gsl_vector * cn_prior = gsl_vector_alloc(nLevels);
  gsl_vector * prior    = gsl_vector_alloc(nLevels);
  gsl_vector * post     = gsl_vector_alloc(nLevels);
  gsl_vector * mem      = gsl_vector_alloc(nLevels);
  double flat = (snpEmit->log_space) ? -log(double(nLevels)) : 1.0/double(nLevels);
  gsl_vector_set_all(prior, flat);
  gsl_vector_set_all(cn_prior, flat);
  gsl_vector_set_all(post,  flat);
  gsl_matrix * Trans = NULL;
  ent = 0.0;
  if (nClones>0 && snpEmit->connect){
    Trans = gsl_matrix_alloc( nLevels, nLevels);
    gsl_matrix_memcpy( Trans, TransMat_snp);
  }
  gsl_vector_view alph;
  gsl_vector_view cnv_post,baf_post;
  double pj = 1.0, norm=0;
  int ncn = normal_copy[snpEmit->chr[sample]];
  int cn=0, ocn=-1;
  int cnv_evt=-1, last_cnv_evt=-1;
  int baf_evt=-1, last_baf_evt=-1, baf_idx=0;
  int idx=0, nidx=0;
  int last_evt = snpEmit->nEvents[sample]-1;
  int last_idx = snpEmit->idx_of_event[sample][last_evt];
  for (int evt = last_evt; evt >= 0 ; evt--){
    idx = snpEmit->idx_of_event[sample][evt];
    //***PREDICTION STEP***
    if (nClones > 0){
      if ( snpEmit->connect && evt < last_evt){//connect to the right...
	pj = snpEmit->pjump[sample][last_idx];
	last_idx = idx;
	Clone::predict( prior, post, snpEmit, pj, Trans);
      }//...connect to right done
      if (cnvEmit->is_set){//connect to CNV...
	if (bafEmit->is_set){
	  baf_evt = snpEmit->Event_of_idx[sample][idx];
	  baf_idx = bafEmit->idx_of_event[baf_sample][baf_evt];
	  cnv_evt = bafEmit->Event_of_idx[baf_sample][baf_idx];
	}
	else{
	  cnv_evt = snpEmit->Event_of_idx[sample][idx];
	}
	if (cnv_evt != last_cnv_evt || baf_evt != last_baf_evt){
	  cnv_post = gsl_matrix_row( gamma_cnv[cnv_sample], cnv_evt);
	  if (bafEmit->is_set){
	    baf_post = gsl_matrix_row( gamma_baf[baf_sample], baf_evt);
	    get_snp_prior_from_cnv_baf_post( cn_prior, &cnv_post.vector, &baf_post.vector);
	  }
	  else{
	    get_snp_prior_from_cnv_post( cn_prior, &cnv_post.vector);
	  }
	  last_cnv_evt = cnv_evt;
	  last_baf_evt = baf_evt;
	}
	if (snpEmit->connect) gsl_vector_memcpy( mem, prior);
	gsl_vector_memcpy( prior, cn_prior);
	nidx = (evt < last_evt) ? snpEmit->idx_of_event[sample][evt+1] : snpEmit->nSites[sample];
	if ( nidx-idx > 1){//exponentiate prior for all observations in this block
	  gsl_vector_scale( prior, double(nidx-idx));//log-space!
	  norm = log_vector_norm(prior);
	  gsl_vector_add_constant(prior,-norm);
	}  
	// multiply two priors and rescale...
	if (snpEmit->connect){
	  if (snpEmit->log_space){
	    gsl_vector_add(prior,mem);
	    norm = log_vector_norm(prior);
	    gsl_vector_add_constant(prior,-norm);
	  }
	  else{
	    gsl_vector_mul(prior,mem);
	    norm = gsl_blas_dasum(prior);
	    if (norm!=norm || norm <0.0) abort();
	    gsl_vector_scale(prior,1.0/norm);
	  }
	}
      }//connect to CNV done.
      else if (!cnvEmit->is_set && !snpEmit->connect){//cnv not set, no correlations...
	cn = (snpEmit->cnmax == NULL) ? ncn : snpEmit->cnmax[sample][evt];
	if (cn != ocn){
	  gsl_vector_memcpy( cn_prior, cn_prior_snp[cn]);
	  ocn = cn;
	}//or flat over all levels?
	gsl_vector_memcpy( prior, cn_prior);
      }
    }  
    else if (nClones == 0){//just one element!
      gsl_vector_set_all( prior, flat);
    }
    //***GET POSTERIOR***
    gsl_vector_memcpy( mem, prior);
    alph = gsl_matrix_row( alpha_snp[sample], evt);
    if (snpEmit->log_space){
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
    gsl_matrix_set_row( gamma_snp[sample], evt, mem);
    //ent += Clone::entropy(mem);
    //***UPDATE STEP*** (normalization term not needed here)
    Clone::update( prior, post, snpEmit, sample, evt);
  }
  // cleanup    
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(mem);
  gsl_vector_free(cn_prior);
  if (Trans!= NULL) gsl_matrix_free(Trans);
}


double Clone::entropy(gsl_vector * x){
  double H = 0.0;
  for (int i=0; i<(int) x->size; i++){
    if (x->data[i] > 0.0) H -= x->data[i] * log(x->data[i]);
  }
  return(H);
}

// predict step...
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
	prior->data[l] = prior->data[l]>0.0 ? log(prior->data[l]) : -1.0e6;
      }
    }
  }
}

void Clone::predict( gsl_vector * prior, gsl_vector * post, Emission * myEmit, double pj, double flat){
  if (pj==0.0){
    gsl_vector_memcpy(prior,post);
  }
  else{
    gsl_vector_set_all(prior,flat);
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



// general update step...
double Clone::update( gsl_vector * prior, gsl_vector * post, Emission * myEmit, int sample, int evt){ 
  if (myEmit == cnvEmit){
    Clone::update_cnv( post, sample, evt);
  }
  else if (myEmit == bafEmit){
    Clone::update_baf( post, sample, evt);
  }
  else if (myEmit == snpEmit){
    Clone::update_snp( prior, post, sample, evt);
  }
  // update step: multiply prior with the likelihood and scale to get the posterior...
  if (myEmit->log_space){
    gsl_vector_add( post, prior);
    double norm = log_vector_norm( post );
    if ( norm != norm ) abort();
    gsl_vector_add_constant( post, -norm);
    return(norm);
  }
  else{
    if (gsl_vector_isnull(post)){
      gsl_vector_memcpy( post, prior);
      return(-1.0e10);
    }
    else{
      gsl_vector_mul( post, prior);
      double norm = gsl_blas_dasum(post);
      if ( norm != norm || norm <= 0.0) abort();
      gsl_vector_scale( post, 1.0 / norm);
      return( log(norm) );
    }
  }
}


//to be precomputed before each CNV-Fwd run!
void Clone::set_tcn(int s){//only used for cnvEmit...
  if (cnvEmit->is_set==0) abort();
  for (int t=0; t<nTimes; t++){
    int ncn = normal_copy[ cnvEmit->chr[s] ];
    if ( tcn[t][s] != NULL)     delete [] tcn[t][s];
    if ( log_tcn[t][s] != NULL) delete [] log_tcn[t][s];
    tcn[t][s]     = new double [nLevels];
    log_tcn[t][s] = new double [nLevels];
    for (int l=0; l<nLevels; l++){
      tcn[t][s][l] = double(ncn)*(1.0-purity[t]) + clone_spectrum[t][l];
      log_tcn[t][s][l] = log(tcn[t][s][l] + 1.0e-10);
    }
  }
}





// update step for CNV...
void Clone::update_cnv( gsl_vector * post, int sample, int evt){ 
  if ( cnvEmit->coarse_grained ){
    Clone::update_cnv_event( post, sample, evt);
  }
  else{
    if (nClones==0){
      Clone::update_cnv_site_noclone( post, sample, evt);
    }
    else{
      Clone::update_cnv_site_wclone( post, sample, evt);
    }
  }
}


void Clone::update_cnv_event( gsl_vector * post, int sample, int evt){
  gsl_vector_set_all( post, 0.0);//log-space!
  gsl_vector_view emit;
  double x,dx,val;
  for (int t=0; t<nTimes; t++){
    emit = gsl_matrix_row( cnvEmitLog[t][sample], evt);
    dx   = (cnv_xmax[t] - cnv_xmin[t]) / double(cnvEmit->gridSize);
    for (int level=0; level<nLevels; level++){
      x = tcn[t][sample][level] * mass->data[t];
      if ( x < cnv_xmin[t] || x > cnv_xmax[t]){//outside range...
	val = -1.0e10;
      }
      else{//inside range...
	val = Clone::get_interpolation( x, cnv_xmin[t], cnv_xmax[t], dx, &emit.vector);
      }
      post->data[level] += val;
    }
  }
}



void Clone::update_cnv_site_noclone( gsl_vector * post, int sample, int site){//no-clone special case...
  gsl_vector_set_all( post, 1.0);
  double val=0,rate=0;
  int ncn = normal_copy[cnvEmit->chr[sample]];
  double b  =  (cnvEmit->bias != NULL) ? cnvEmit->bias[sample][site] : -1.0;
  double lb =  (cnvEmit->log_bias != NULL) ? cnvEmit->log_bias[sample][site] : 0.0;
  for (int time=0; time<nTimes; time++){
    int N = cnvEmit->depths[time][sample][site];
    int n = cnvEmit->reads[time][sample][site];   
    double rnd = cnvEmit->rnd_emit / (double(N)*(cnvEmit->maxRate - cnvEmit->minRate) + 1.0);
    double nrnd = (1.0 - cnvEmit->rnd_emit);
    if (cnvEmit->mode == 3){//Poisson
      val  = - mass->data[time] * double(N*ncn);
      if (b>0.0) val *= b;
      val -= loggma[n+1];
      rate = log_mass->data[time] + logn[N] + logn[ncn];
      if (b>0.0) rate += lb;
      val += double(n) * rate;    
    }
    else if (cnvEmit->mode == 4){//Negative Binomial
      double n1   = double(N)*cnvEmit->shape;//rounding (error?)
      double nrml = double(ncn);
      val  = loggma[int(n1)+n] - loggma[int(n1)] - loggma[n+1] + n1*cnvEmit->log_shape;
      rate = mass->data[time] * nrml;
      if (b>0.0) rate *= b;
      val += double(n) * log(rate);
      val -= (double(n) + n1) * log(rate + cnvEmit->shape);
      if (val != val) abort();
    }
    val = exp(val);
    if (rnd > 0.0) val = val*nrnd + rnd;
    post->data[0] *= val; 
  }
}

//normal case with one or more clones
void Clone::update_cnv_site_wclone( gsl_vector * post, int sample, int site){
  gsl_vector_set_all( post, 1.0);
  double b  = (cnvEmit->bias != NULL) ? cnvEmit->bias[sample][site] : -1.0;
  double lb = (cnvEmit->log_bias != NULL) ? cnvEmit->log_bias[sample][site] : 0.0;
  for (int time=0; time<nTimes; time++){
    int N = cnvEmit->depths[time][sample][site];
    int n = cnvEmit->reads[time][sample][site];
    double rnd  = cnvEmit->rnd_emit / (double(N)*(cnvEmit->maxRate - cnvEmit->minRate) + 1.0);
    double nrnd = 1.0 - cnvEmit->rnd_emit;
    double val=0;
    if (cnvEmit->mode==3){//Poisson
      double pre1 = double(n) * (log_mass->data[time] + logn[N]) - loggma[n+1];
      double pre2 = mass->data[time] * double(N);
      if (b > 0.0){//modulation
	pre1 += double(n) * lb;
	pre2 *= b;
      }
      for (int l=0; l<nLevels; l++){
	val = pre1 - pre2*tcn[time][sample][l] + double(n)*log_tcn[time][sample][l];
	val = exp(val);
	if (rnd > 0.0) val = val*nrnd + rnd;
	post->data[l] *= val;
      }
    }
    else if (cnvEmit->mode == 4){//Negative Binomial
      double n1 = double(N)*cnvEmit->shape;//rounding (error?)
      double pre = loggma[int(n1)+n] - loggma[int(n1)] - loggma[n+1] + n1*cnvEmit->log_shape;
      double rate;
      for (int l=0; l<nLevels; l++){
	rate = mass->data[time] * tcn[time][sample][l];
	val  = pre +  double(n) * (log_mass->data[time] + log_tcn[time][sample][l]);
	if (b > 0.0){
	  rate *= b;
	  val += double(n)*lb;
	}
	val  = val - (double(n) + n1) * log(rate + cnvEmit->shape);
	val = exp(val);
	if (rnd > 0.0) val = val*nrnd + rnd;
	post->data[l] *= val;
      }
    }
  }
}


void Clone::update_baf( gsl_vector * post, int sample, int evt){
  if (bafEmit->coarse_grained){
    Clone::update_baf_event( post, sample, evt);
  }
  else{
    Clone::update_baf_site( post, sample, evt);
  }
}

void Clone::update_baf_event( gsl_vector * post, int sample, int evt){
  gsl_vector_set_all( post, 0.0);//log-space!
  double x,y,val,total_cn;
  double dx = bafEmit->dx;
  gsl_vector_view emit;
  for (int t=0; t<nTimes; t++){
    total_cn = (bafEmit->phi == NULL) ? 2.0 : bafEmit->phi[t][sample][evt];//only diploid chromosomes!
    y = (1.0 - purity[t]) / total_cn;
    emit = gsl_matrix_row( bafEmitLog[t][sample], evt);
    for (int level=0; level<nLevels; level++){
      x = y + clone_spectrum[t][level] / total_cn;
      if ( x < 0.0 || x > 1.0){//outside range...
	val = -1.0e6;
      }
      else{//inside range...
	val = Clone::get_interpolation( x, 0.0, 1.0, dx, &emit.vector);
      }
      post->data[level] += val;
    }
  }
}


void Clone::update_baf_site( gsl_vector * post, int sample, int site){
  gsl_vector_set_all( post, 1.0);
  double x, y, total_cn, val;
  unsigned int N, n1, n2;
  double dx = bafEmit->dx;
  int evt   = bafEmit->event_of_idx[sample][site];
  gsl_vector * emit1=NULL, * emit2=NULL;
  for (int t=0; t<nTimes; t++){
    N  = bafEmit->depths[t][sample][site];
    n1 = bafEmit->reads[t][sample][site];
    emit1 = bafEmit->EmitLog[N][n1];
    n2 = N-n1;
    if (n2 != n1) emit2 = bafEmit->EmitLog[N][n2];
    double rnd = bafEmit->rnd_emit / double(N+1);
    total_cn = (bafEmit->phi == NULL) ? 2.0 : bafEmit->phi[t][sample][evt];
    y = (1.0-purity[t]) / total_cn;
    for (int level=0; level<nLevels; level++){
      x = y + clone_spectrum[t][level] / total_cn;
      if (x > 1.0){//outside range...
	val = 0.0;
      }
      else if (x<1.0 && x>0.0){//inside range...
	val = exp( Clone::get_interpolation( x, 0.0, 1.0, dx, emit1) );
	if (n1 != n2) val += exp( Clone::get_interpolation( x, 0.0, 1.0, dx, emit2) );
      }
      else if ( x==1.0 ){//edges...
	if (n1==N || n2==N){
	  val = 1.0;
	}
	else{
	  val = 0.0;
	}
      }
      else if ( x==0.0 ){
	if (n1==0 || n2==0){
	  val = 1.0;
	}
	else{
	  val = 0.0;
	}
      }
      else{
	abort();
      }
      if (bafEmit->rnd_emit > 0.0) val = val*(1.0-bafEmit->rnd_emit) + rnd;
      if (val!=val || val < 0.0) abort();
      post->data[level] *= val;
      //printf("x=%e n=%i N=%i p=%e\n", x,n1,N,val);
      //exit(0);
    }
  }
}


void Clone::update_snp( gsl_vector * prior, gsl_vector * post, int sample, int evt){
  if (snpEmit->coarse_grained){
    Clone::update_snp_event( post, sample, evt);
  }
  else if ( bulk_fix >= 0.0 || (bulk_mean != NULL && bulk_dist == NULL)){
    Clone::update_snp_fixed( prior, post, sample, evt);
  }
  else{
    Clone::update_snp_nfixed( prior, post, sample, evt);
  }
}

void Clone::update_snp_event( gsl_vector * post, int sample, int evt){
  gsl_vector_set_all( post, 0.0);//log-space!
  double val,total_cn;
  double ncn = double(normal_copy[snpEmit->chr[sample]]);
  for (int t=0; t<nTimes; t++){
    total_cn = (snpEmit->phi == NULL) ? ncn : snpEmit->phi[t][sample][evt];
    double b = (1.0 - purity[t]) * ncn / total_cn;
    for (int level=0; level<nLevels; level++){
      double a = clone_spectrum[t][level] / total_cn;
      if (a > 1.0){//outside range...
	int idx  = snpEmit->idx_of_event[sample][evt];
	int nidx = (evt < snpEmit->nEvents[sample]-1) ? snpEmit->idx_of_event[sample][evt+1] : snpEmit->nSites[sample];
	val = -1.0e6 * double(nidx-idx);
      }
      else{//inside range...
	val = Clone::get_interpolation( a, 0.0, 1.0, b, 0.0, 1.0/bulk_min[t][sample][evt], snpEmitLog[t][sample][evt]);
      }
      post->data[level] += val;
    }
  }
  /*
  for (int level=0; level<nLevels; level++) printf("%.2e ", post->data[level]);
  cout<<endl;
  */
}

//Dirac-Delta bulk prior distribution...
void Clone::update_snp_fixed( gsl_vector * prior, 
			      gsl_vector * post, 
			      int sample, 
			      int site){
  gsl_vector_set_all( post, 1.0);
  int ncn = normal_copy[snpEmit->chr[sample]];
  double x,y,val=0;
  unsigned int N,n;
  int evt = snpEmit->event_of_idx[sample][site];
  double dx = snpEmit->dx;
  gsl_vector * emit = NULL;
  for (int t=0; t<nTimes; t++){
    N = snpEmit->depths[t][sample][site];
    n = snpEmit->reads[t][sample][site];
    if (N==0){
      continue;
    }
    emit = snpEmit->EmitLog[N][n];
    double rnd  = snpEmit->rnd_emit / double(N+1);
    double nrnd = 1.0 - snpEmit->rnd_emit;
    double total_cn = (snpEmit->phi == NULL) ? double(ncn) : snpEmit->phi[t][sample][evt];
    double bfix = (bulk_fix >= 0.0) ? bulk_fix : bulk_mean[t][sample][site];
    y =  bfix * (1.0-purity[t]) * double(ncn) / total_cn;
    for (int level=0; level<nLevels; level++){
      if (prior->data[level] <= 0.0){
	post->data[level] = 0.0;
	continue;
      }
      x = y + clone_spectrum[t][level] / total_cn; 
      if (snp_err>0.0){
	x = snp_err + (1.0-snp_err)*x;
      }
      if (x > 1.0){//outside range...
	val = 0.0;
      }
      else if (x<1.0 && x>0.0){//inside range...
	val = Clone::get_interpolation( x, 0.0, 1.0, dx, emit);
	val = exp(val);
      }
      else if ( x==1.0 ){//edges...
	val = (n==N) ? 1.0 : 0.0;
      }
      else if ( x==0.0 ){
	val = (n==0) ? 1.0 : 0.0;
      }
      else{
	abort();
      }
      if (val != val || val < 0.0) abort();
      if (rnd > 0.0) val = val*nrnd + rnd;
      post->data[level] *= val;
    }
  }
  /*
  for (int level=0; level<nLevels; level++) 
    printf("%.2e ", post->data[level] > 0.0 ? log(post->data[level]) : -1.0e6);
  cout<<endl;
  */
}




 



//non-trivial case with one or more clones
void Clone::update_snp_nfixed( gsl_vector * prior, 
			       gsl_vector * post, 
			       int sample, 
			       int site){
  gsl_vector_set_all(post, 1.0);
  int ncn = normal_copy[snpEmit->chr[sample]];
  double val=0;
  int evt = snpEmit->event_of_idx[sample][site];
  for (int time=0; time<nTimes; time++){//product over time points
    unsigned int N = snpEmit->depths[time][sample][site];
    unsigned int n = snpEmit->reads[time][sample][site];
    gsl_vector * emit = snpEmit->EmitProb[N][n];
    double rnd  = snpEmit->rnd_emit / double(N+1);
    double nrnd = 1.0 - snpEmit->rnd_emit;
    //evaluate likelihood function at a finite grid of clonal frequencies
    int nPts = (nLevels <= 100 ) ? nLevels-1 : 100;
    //grid for linear interpolation//***TEST***
    gsl_vector * mem = gsl_vector_calloc(nPts+1);
    double total_cn = (snpEmit->phi == NULL) ? double(ncn) : snpEmit->phi[time][sample][evt];
    double maxq = 0.0;
    for (int l=0; l<nLevels; l++){
      if (clone_spectrum[time][l] / total_cn <= 1.0){
	maxq = max( maxq, clone_spectrum[time][l] / total_cn);
      }
    }
    double dq = maxq / double(nPts);
    gsl_vector_view blk = gsl_matrix_row( bulk_dist[time][sample], site);
    double b = (1.0-purity[time])*double(ncn) / total_cn;
    for (int pt=0; pt<=nPts; pt++){
      double a = (nPts == nLevels-1) ? clone_spectrum[time][pt] / total_cn : double(pt)*dq;
      if( N == 0 ){
	val =  1.0;// no observation -> all states equally likely!
      }
      else if (a > 1.0){
	val = 0.0;
      }
      else{
	if ( nPts != nLevels-1 ){
	  val = Clone::trapezoidal( &blk.vector, a, b, emit, 1);//log
	}
	else if (prior->data[pt] > 0.0){
	  val = Clone::trapezoidal( &blk.vector, a, b, emit, 0);//linear
	}
      }	
      gsl_vector_set( mem, pt, val);
    }
    // mutiply into posterior...
    if (nPts == nLevels-1){// we have computed lh at all the correct points
      if (rnd > 0.0){
	gsl_vector_scale(mem, nrnd);
	gsl_vector_add_constant( mem, rnd);
      }
      gsl_vector_mul( post, mem);
    }
    else{//need to approximate lh via linear interpolation...
      double a,nu,f1,f2;
      int low;	
      // remaining levels...
      for (int level=0; level<nLevels; level++){
	if (prior->data[level] <= 0.0) continue;
	a = clone_spectrum[time][level] / total_cn;	 
	if (a>1.0){
	  val = 0.0;
	}
	else{
	  low = int(a/dq);
	  nu = a/dq - double(low);
	  f1 = gsl_vector_get( mem, low);
	  f2 = (low<nPts) ? gsl_vector_get( mem, low + 1) : f1;
	  val = (1.0-nu)*f1 + nu*f2;//linear interpolation for log-value
	  val = exp(val);
	  if (snpEmit->rnd_emit > 0.0) val = val*(1.0-snpEmit->rnd_emit) + rnd;
	}
	post->data[level] *= val;
      }
    }
    gsl_vector_free(mem);
  }
}

//1-D interpolation...
double Clone::get_interpolation(double x, double xmin, double xmax, double dx, gsl_vector * emit){
  int idx   = int((x-xmin)/dx);
  double nu = (x-xmin)/dx - double(idx);
  double f1,f2;
  int gs = (int) emit->size - 1;
  if (idx > 0 && idx < gs){//use log-emission probability
    f1 = emit->data[idx];
    f2 = emit->data[idx+1];
  }
  else if (idx==0){//else use linear extra-polation of log
    if (x>xmin){
      f1 = emit->data[1];
      f2 = emit->data[2];
      nu -= 1.0;
    }
    else{
      f1=-1.0e6;
      f2=-1.0e6;
    }
  }
  else{
    if (x<xmax){
      f1 = emit->data[gs-2];
      f2 = emit->data[gs-1];
      nu += 1.0;
    }
    else{
      f1=-1.0e6;
      f2=-1.0e6;
    }
  }
  double val = (1.0-nu)*f1 + nu*f2;//linear interpolation (in log-space)
  return(val);
}

//2-D interpolation...
double Clone::get_interpolation(double x, double xmin, double xmax,
				double y, double ymin, double ymax,
				gsl_matrix * emit){
  double dx = (xmax-xmin) / (double(emit->size1 - 1));
  double dy = (ymax-ymin) / (double(emit->size2 - 1));
  int ix1 = int((x-xmin)/dx);
  int iy1 = int((y-ymin)/dy);
  double f,f1,f2,f3,f4;
  double t = (x-(xmin+ix1*dx)) / dx;
  double u = (y-(ymin+iy1*dy)) / dy;  
  f1 = gsl_matrix_get(emit, ix1,   iy1);
  f2 = gsl_matrix_get(emit, ix1+1, iy1);
  f3 = gsl_matrix_get(emit, ix1+1, iy1+1);
  f4 = gsl_matrix_get(emit, ix1,   iy1+1);  
  f  = (1.0-t)*(1.0-u)*f1 + t*(1.0-u)*f2 + t*u*f3 + (1.0-t)*u*f4;
  return(f);
}



//only ever used on snpEmit!
double Clone::trapezoidal( gsl_vector * blk, double a, double b, gsl_vector * emit, int get_log){
  int olow = -1;
  double dx = snpEmit->dx;
  double dy = dx * b;
  double y = a - dy;
  double f1=0,f2=0,nu=0,pre=0;
  int low,high;
  int gs = snpEmit->gridSize;
  if ( gsl_vector_isnull(blk) ) abort();
  double val = 0.0, dval=0;
  for (int j=0; j <= gs; j++){
    y += dy;
    if (y > 1.0) break;
    low  = int( y / dx);
    if ( low != olow){
      f1 = emit->data[low];
      if (low < gs){
	high = low+1;
	f2   = emit->data[high];
      }
      else{
	high = low;
	f2   = f1;
      }
      olow = low;
    }
    // linear interpolation of 1-D distribution between two grid-points...
    pre = (j==0 || j==gs) ? 0.5 : 1.0;
    nu = y / dx - double(low);
    dval = pre * blk->data[j] * ( (1.0-nu)*f1 + nu*f2 );
    val += dval;
    if (dval < 1.0e-6*val) break;
  }
  val *= dx;
  if (val != val) abort();
  if (get_log==1) val = (val > 0.0) ? log(val) : -1.0e6;
  return(val);
}




//get mean total copy number...
void Clone::get_phi(int sample){//only ever used for cnvEmit
  if (nClones == 0){
    double ncn = normal_copy[cnvEmit->chr[sample]];
    for (int t=0; t<nTimes; t++){
      for (int evt=0; evt < cnvEmit->nEvents[sample]; evt++){
	cnvEmit->phi[t][sample][evt] = double(ncn);
      }
    }
    for (int evt=0; evt<cnvEmit->nEvents[sample]; evt++){
      cnvEmit->cnmax[sample][evt] = ncn;
    }
  }
  else{
    if (gamma_cnv[sample] == NULL) abort();
    double val=0;
    gsl_vector * post1 = gsl_vector_alloc(nLevels);
    gsl_vector * post2 = gsl_vector_alloc(maxcn+1);
    int ncn = normal_copy[ cnvEmit->chr[sample] ];
    for (int t=0; t<nTimes; t++){
      double nrml = (1.0-purity[t]) * double(ncn);
      for (int evt=0; evt < cnvEmit->nEvents[sample]; evt++){
	cnvEmit->phi[t][sample][evt] = 0.0;
	if (t==0){
	  gsl_vector_set_zero(post2);
	}
	for (int l=0; l<nLevels; l++){
	  val = gsl_matrix_get( gamma_cnv[sample], evt, l);
	  post1->data[l] = val;
	  cnvEmit->phi[t][sample][evt] += val * (nrml + clone_spectrum[t][l]);
	  if (t==0){
	    int mx = *std::max_element( copynumber[l], copynumber[l] + nClones);
	    post2->data[mx] += val;
	  }
	}
	if (t==0){//get the maximum copynumber across clones...
	  cnvEmit->cnmax[sample][evt] = 0;
	  double p=0.0;
	  for (int cn=0; cn<=maxcn; cn++){
	    p += post2->data[cn];
	    if (p > 0.99){//conservative estimate of the highest copynumber (with pError < 1%)
	      cnvEmit->cnmax[sample][evt] = cn;
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

void Clone::map_phi( Emission * fromEmit, int from_sample, Emission * toEmit){//only ever done for snpEmit and bafEmit
  //fill
  int sample = toEmit->idx_of[fromEmit->chr[from_sample]];
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
      for (int t=0; t<nTimes; t++){//total copynumber
	toEmit->phi[t][sample][evt] = fromEmit->phi[t][from_sample][from_evt];
      }
      //max copynumber
      toEmit->cnmax[sample][evt] = fromEmit->cnmax[from_sample][from_evt];
    }
  }
}



void Clone::get_cnvEmitLog(){//only ever used for BAF data
  if (cnvEmitLog != NULL) abort();
  cnv_xmin = new double [nTimes];
  cnv_xmax = new double [nTimes];
  cnvEmitLog = new gsl_matrix ** [nTimes];
  double lnrnd = (cnvEmit->rnd_emit>0.0) ? log(1.0-cnvEmit->rnd_emit) : 0.0;
  for (int t=0; t<nTimes; t++){
    cnvEmit->get_log = 1;
    cnvEmit->init_range(t);
    cnv_xmin[t] = cnvEmit->xmin;
    cnv_xmax[t] = cnvEmit->xmax;
    cnvEmit->set_EmitProb(t);
    cnvEmitLog[t] = new gsl_matrix * [cnvEmit->nSamples];
    for (int s=0; s<cnvEmit->nSamples; s++){
      cnvEmitLog[t][s] = gsl_matrix_calloc( cnvEmit->nEvents[s], cnvEmit->gridSize+1);
      int evt;
#pragma omp parallel for schedule( dynamic, 1) default(shared)      
      for ( evt=0; evt<cnvEmit->nEvents[s]; evt++){
	gsl_vector * mem = gsl_vector_calloc(cnvEmit->gridSize+1);
	gsl_vector * rnd = gsl_vector_calloc(cnvEmit->gridSize+1);
	unsigned int n,N;
	double r=0;
	gsl_vector_view emit = gsl_matrix_row( cnvEmitLog[t][s], evt);
	int idxi = cnvEmit->idx_of_event[s][evt];
	int idxf = (evt<cnvEmit->nEvents[s]-1) ? cnvEmit->idx_of_event[s][evt+1] - 1: cnvEmit->nSites[s]-1;
	for(int idx=idxi; idx<=idxf; idx++){
	  n = cnvEmit->reads[t][s][idx];
	  N = cnvEmit->depths[t][s][idx];
	  if (N>0){//observation?
	    if (cnvEmit->bias == NULL){
	      gsl_vector_memcpy( mem, cnvEmit->EmitLog[N][n]);
	    }
	    else{
	      cnvEmit->get_eprob_wBias( mem, cnvEmit->EmitLog[N][n], cnvEmit->bias[s][idx], n, N, 1);
	    }	  
	    if (cnvEmit->rnd_emit > 0.0){//random?
	      gsl_vector_add_constant( mem, lnrnd);
	      r = double(N)*(cnvEmit->maxRate - cnvEmit->minRate);
	      if (cnvEmit->bias != NULL) r *= cnvEmit->bias[s][idx];
	      r =  cnvEmit->rnd_emit / (r+1.0);
	      gsl_vector_set_all( rnd, log(r) );
	      log_vector_add( mem, rnd);
	    }
	    gsl_vector_add( &emit.vector, mem);
	  }
	}
	gsl_vector_free(mem);
	gsl_vector_free(rnd);
      }
    }  
  }
}


void Clone::get_bafEmitLog(){//only ever used for BAF data
  if (bafEmitLog != NULL) abort();
  bafEmitLog = new gsl_matrix ** [nTimes];
  double lrnd  = (bafEmit->rnd_emit>0.0) ? log(bafEmit->rnd_emit) : 0.0;
  double lnrnd = (bafEmit->rnd_emit>0.0) ? log(1.0-bafEmit->rnd_emit) : 0.0;
  for (int t=0; t<nTimes; t++){
    bafEmitLog[t] = new gsl_matrix * [bafEmit->nSamples];
    for (int s=0; s<bafEmit->nSamples; s++){
      bafEmitLog[t][s] = gsl_matrix_calloc( bafEmit->nEvents[s], bafEmit->gridSize+1);      
      int evt;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
      for (evt=0; evt<bafEmit->nEvents[s]; evt++){
	unsigned int n,N;
	gsl_vector * mem = gsl_vector_calloc(bafEmit->gridSize+1);
	gsl_vector * rnd = gsl_vector_calloc(bafEmit->gridSize+1);
	gsl_vector_view emit = gsl_matrix_row( bafEmitLog[t][s], evt);
	int idxi = bafEmit->idx_of_event[s][evt];
	int idxf = (evt<bafEmit->nEvents[s]-1) ? bafEmit->idx_of_event[s][evt+1] - 1: bafEmit->nSites[s]-1;
	for(int idx=idxi; idx<=idxf; idx++){
	  n = bafEmit->reads[t][s][idx];
	  N = bafEmit->depths[t][s][idx];
	  if (N>0){//observation?
	    gsl_vector_memcpy( mem, bafEmit->EmitLog[N][n]);
	    if (n != N-n) log_vector_add( mem, bafEmit->EmitLog[N][N-n]);	  
	    if (bafEmit->rnd_emit > 0.0){//random?
	      gsl_vector_add_constant( mem, lnrnd);
	      gsl_vector_set_all( rnd, lrnd - log(double(N+1)) );
	      log_vector_add( mem, rnd);
	    }
	    gsl_vector_add( &emit.vector, mem);
	  }
	}
	gsl_vector_free(mem);
	gsl_vector_free(rnd);
      }      
    }  
  }  
  //cout<<"done."<<endl;
}





void Clone::get_snpEmitLog(){//only ever used for SNP data
  if (!snpEmit->is_set) abort();
  if (snpEmitLog==NULL){//first allocation
    snpEmitLog = new gsl_matrix *** [nTimes];
    for (int t=0;t<nTimes; t++){
      snpEmitLog[t] = new gsl_matrix ** [snpEmit->nSamples];
      for (int s=0; s<snpEmit->nSamples; s++){
	snpEmitLog[t][s] = new gsl_matrix * [snpEmit->nEvents[s]];
	for (int evt=0; evt< snpEmit->nEvents[s]; evt++){
	  snpEmitLog[t][s][evt] = gsl_matrix_calloc(snpEmit->gridSize+1,snpEmit->gridSize+1);
	}
      }
    }
  }
  Clone::get_bulk_min();
  double dx = snpEmit->dx;
  double lnrnd = log(1.0 - snpEmit->rnd_emit);
  for (int t=0; t<nTimes; t++){   
    for (int s=0; s<snpEmit->nSamples; s++){      
      int evt;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
      for (evt=0; evt<snpEmit->nEvents[s]; evt++){
	gsl_matrix * mem     = gsl_matrix_calloc(snpEmit->gridSize+1,snpEmit->gridSize+1);
	gsl_matrix * rnd_vec = gsl_matrix_calloc(snpEmit->gridSize+1,snpEmit->gridSize+1);
	gsl_matrix_set_zero(snpEmitLog[t][s][evt]);
	int idxi = snpEmit->idx_of_event[s][evt];
	int idxf = (evt<snpEmit->nEvents[s]-1) ? snpEmit->idx_of_event[s][evt+1] - 1: snpEmit->nSites[s]-1;
	for(int idx=idxi; idx<=idxf; idx++){
	  double val=0;
	  unsigned int N = snpEmit->depths[t][s][idx];
	  if (N==0) continue;//no observation
	  unsigned int n = snpEmit->reads[t][s][idx];
	  gsl_vector * emit = (bulk_prior==NULL) ? snpEmit->EmitLog[N][n] : snpEmit->EmitProb[N][n];
	  double lrnd  = (snpEmit->rnd_emit>0.0) ? log(snpEmit->rnd_emit/double(N+1)) : 0.0;	  
	  //evaluate likelihood function at a finite grid of clonal frequencies
	  double bfix = (bulk_fix >= 0.0) ? bulk_fix : bulk_mean[t][s][idx];
	  gsl_vector_view blk;
	  if (bulk_prior != NULL) blk = gsl_matrix_row( bulk_dist[t][s], idx);
	  for (int i=0; i<=snpEmit->gridSize; i++){
	    double a = double(i) / double(snpEmit->gridSize);
	    for (int j=0; j<=snpEmit->gridSize; j++){
	      double b = double(j) / (double(snpEmit->gridSize) * bulk_min[t][s][evt]);
	      if( b > (1.0-a) / bulk_min[t][s][evt]){
		val = -1.0e6;
	      }
	      else if (bulk_prior==NULL){//bulk point estimate...	      
		double x =  a + bfix * b;
		if (x > 1.0){//outside range...
		  val = -1.0e6;
		}
		else if (x>0.0 && x<1.0){//inside range...
		  val = Clone::get_interpolation( x, 0.0, 1.0, dx, emit);
		}
		else if ( x==1.0 ){//edges...
		  val = (n==N) ? 0.0 : -1.0e6;
		}
		else if ( x==0.0 ){
		  val = (n==0) ? 0.0 : -1.0e6;
		}
	      }
	      else{//bulk full distribution...		
		val = Clone::trapezoidal( &blk.vector, a, b, emit, 1);//log-space
	      }	
	      gsl_matrix_set( mem, i, j, val);
	    }
	  }
	  // mutiply into posterior...
	  if (snpEmit->rnd_emit > 0.0){
	    gsl_matrix_add_constant( mem, lnrnd);
	    gsl_matrix_set_all( rnd_vec, lrnd);
	    log_matrix_add( mem, rnd_vec);
	  }
	  gsl_matrix_add( snpEmitLog[t][s][evt], mem);
	}
	gsl_matrix_free(mem);
	gsl_matrix_free(rnd_vec);
      }
    }  
  }
}




// Bayesian update of the SNP bulk
void Clone::update_bulk(int sample){
  if (!snpEmit->is_set) abort();
  if (bulk_mean == NULL) abort();
  gsl_vector * bpost=NULL, *flat=NULL, *emit=NULL;
  gsl_vector_view bprior;
  bpost = gsl_vector_alloc(snpEmit->gridSize+1);
  flat  = gsl_vector_alloc(snpEmit->gridSize+1);
  gsl_vector_set_all(flat,1.0);
  //get SNP posterior...
  alpha_snp[sample]=NULL;
  gamma_snp[sample]=NULL;
  Clone::get_snp_posterior(sample);
  //update the bulk
  for (int time=0; time<nTimes; time++){   
    unsigned int n,N;
    for (int idx=0; idx<snpEmit->nSites[sample]; idx++){
      n = snpEmit->reads[time][sample][idx];
      N = snpEmit->depths[time][sample][idx];
      if (N>0){
	emit = snpEmit->EmitProb[N][n];
	if (bulk_prior == NULL){ 
	  Clone::get_bulk_post_dist( flat, bpost, emit, time, sample, idx);
	  bulk_post_mean[time][sample][idx] 
	    = 0.5*(bulk_prior_mean[sample][idx] + snpEmit->xgrid[gsl_vector_max_index(bpost)]);
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
  gsl_matrix_free(gamma_snp[sample]);
  gamma_snp[sample] = NULL;
  if (bpost!=NULL) gsl_vector_free(bpost);
  if (flat!=NULL)  gsl_vector_free(flat);
}

// emission probability
void Clone::get_bulk_post_dist( gsl_vector * bprior, gsl_vector * bpost, gsl_vector * emit, int time, int sample, int idx){
  double dx = snpEmit->dx;
  double prob=0, x=0, y=0, val=0 , f1=0, f2=0, nu=0;
  int old=-1, i1=0, i2=0;
  int gS = snpEmit->gridSize;
  double ncn = (double) normal_copy[snpEmit->chr[sample]];
  int evt = snpEmit->event_of_idx[sample][idx];
  double total_cn = (snpEmit->phi==NULL) ? ncn : snpEmit->phi[time][sample][evt];
  gsl_vector_set_zero(bpost);
  for (int i=0; i <= gS; i++){
    y = double(i)*dx * (1.0 - purity[time]) * ncn;
    prob=0.0;
    old = -1;
    for (int j=0; j<nLevels; j++){
      x = (y + clone_spectrum[time][j]) / total_cn;
      if (x > 1.0) continue;
      i1 = int(x/dx);
      nu = x/dx - double(i1);
      if (i1 != old){
	f1 = emit->data[i1];
	i2 = (i1<gS) ? i1+1 : i1;
	f2 = emit->data[i2];
	old=i1;
      }
      val = (1.0-nu)*f1 + nu*f2;
      val *= gsl_matrix_get( gamma_snp[sample], evt, j);
      prob += val;
    }
    if (prob != prob) abort();
    gsl_vector_set( bpost, i, prob);
  }
  gsl_vector_mul(bpost,bprior);
  double norm = gsl_blas_dasum(bpost) - 0.5*(bpost->data[0] + bpost->data[gS]);
  norm *= dx;
  gsl_vector_scale(bpost,1.0/norm);
}


void Clone::get_bulk_min(){
  if (bulk_min==NULL){
    bulk_min = new double ** [nTimes];
    for (int t=0; t<nTimes; t++){
      bulk_min[t] = new double * [snpEmit->nSamples];
      for (int s=0; s<snpEmit->nSamples; s++){
	bulk_min[t][s] = new double [snpEmit->nEvents[s]];
      }
    }
  }
  for (int t=0; t<nTimes; t++){
    for (int s=0; s<snpEmit->nSamples; s++){
      for (int evt=0; evt < snpEmit->nEvents[s]; evt++){
	bulk_min[t][s][evt] = 1.1;
	int idxi = snpEmit->idx_of_event[s][evt];
	int idxf = (evt<snpEmit->nEvents[s]-1) ? snpEmit->idx_of_event[s][evt+1]-1 : snpEmit->nSites[s]-1;
	for (int idx=idxi; idx<=idxf; idx++){
	  bulk_min[t][s][evt] = min( bulk_min[t][s][evt], bulk_mean[t][s][idx]);
	}
      }
    }
  }
}
