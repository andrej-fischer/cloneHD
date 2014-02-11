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
  snv_prior_from_cna_baf_map  = NULL;
  snv_prior_from_cna_map      = NULL;
  maxcn     = 4;
  logn_set  = 0;
  nFreq     = 0;
  max_nFreq = 4;
  TransMat_cna  = NULL;
  TransMat_snv  = NULL;
  mass      = NULL;
  log_mass  = NULL;
  nmean     = NULL;
  got_gamma = 0;
  mass_candidates = NULL;
  normal_copy      = NULL;
  copynumber       = NULL;
  copynumber_post  = NULL;
  init_cn_prior_snv= NULL;
  clone_spectrum   = NULL;
  majcn_post  = NULL;
  bulk_fix  = -1.0;
  baf_pen = 1.0;
  snv_pen = 0.01;
  snv_fpr = 1.0e-4;
  tcn    = NULL;//total copy number
  log_tcn= NULL;//log thereof
  bafEmitLog = NULL;
  cnaEmitLog = NULL;
  snvEmitLog = NULL;
  //posteriors
  alpha_cna=NULL;
  gamma_cna=NULL;
  alpha_baf=NULL;
  gamma_baf=NULL; 
  alpha_snv=NULL;
  gamma_snv=NULL;
  //
  cna_total_llh=0;
  baf_total_llh=0;
  snv_total_llh=0;
  //bulk variables
  bulk_prior      = NULL;
  bulk_post       = NULL;
  bulk_dist       = NULL;
  bulk_prior_mean = NULL;
  bulk_post_mean  = NULL;
  bulk_mean       = NULL; 
  bulk_min = NULL;
  logzero = -1.0e6;
}

void Clone::allocate( Emission * cna, Emission * baf, Emission * snv, const char * chr_fn){
  cnaEmit = cna;
  bafEmit = baf;
  snvEmit = snv;
  if (cnaEmit->is_set) nTimes  = cnaEmit->nTimes;
  if (bafEmit->is_set) nTimes  = bafEmit->nTimes;
  if (snvEmit->is_set) nTimes  = snvEmit->nTimes;
  total_loci = 0;
  if (cnaEmit->is_set) total_loci += cnaEmit->total_loci;
  if (bafEmit->is_set) total_loci += bafEmit->total_loci;
  if (snvEmit->is_set) total_loci += snvEmit->total_loci;
  // normal copy number
  Clone::set_normal_copy(chr_fn);
  if (cnaEmit->is_set){
    mass     = gsl_vector_alloc(nTimes);
    log_mass = gsl_vector_alloc(nTimes);
    nmean    = gsl_vector_alloc(nTimes);
    Clone::set_logn();
    Clone::get_nmean();
    if (bafEmit->is_set){//CNA + BAF (+SNV)
      for (int s=0; s<bafEmit->nSamples; s++){//all BAF are mapped to CNA
	bafEmit->map_idx_to_Event(cnaEmit,s);
      }
      if (snvEmit->is_set){//if BAF exists, use it for SNVs in autosomes only! 
	for (int s=0; s<snvEmit->nSamples; s++){
	  int snvChr = snvEmit->chr[s];
	  if ( bafEmit->chrs.count(snvChr) == 1){//BAF exists -> autosome 
	    snvEmit->map_idx_to_Event( bafEmit, s);
	  }
	  else{//BAF dne -> hapl. sex chr -> use CNA
	    snvEmit->map_idx_to_Event( cnaEmit, s);
	  }
	}
      }
    }
    else if (snvEmit->is_set){// CNA + SNV
      for (int s=0; s<snvEmit->nSamples; s++){
	snvEmit->map_idx_to_Event(cnaEmit,s);
      }
    }
  }
  //
  purity = new double [nTimes];
  clone_spectrum = new double * [nTimes];
  for (int t=0; t<nTimes; t++) clone_spectrum[t] = NULL;
  min_purity = gsl_vector_calloc(nTimes);
  if (cnaEmit->is_set){
    tcn     = new double ** [nTimes];
    log_tcn = new double ** [nTimes];
    for (int t=0; t<nTimes; t++){
      tcn[t]       = new double * [cnaEmit->nSamples];
      log_tcn[t]   = new double * [cnaEmit->nSamples];
      for (int s=0; s<cnaEmit->nSamples; s++){
	tcn[t][s]     = NULL;
	log_tcn[t][s] = NULL;
      }
    }
  }
  allocated = 1;//done
}

void Clone::set_normal_copy(const char * chr_fn){
  normal_copy = new int [100];
  normal_copy[0] = 0;
  for (int c=0; c<100; c++) normal_copy[c] = 0;
  if (chr_fn == NULL){//if --chr [file] was not given, this is used
    for (int c=1; c<=23; c++) normal_copy[c] = 2;
    normal_copy[24] = 0;//female by default
    //check whether there is a Y chromosome...
    int male=0;
    if (cnaEmit->is_set && cnaEmit->maxchr >= 24 && cnaEmit->idx_of[24] >= 0) male=1;
    if (bafEmit->is_set && bafEmit->maxchr >= 24 && bafEmit->idx_of[24] >= 0) male=1;
    if (snvEmit->is_set && snvEmit->maxchr >= 24 && snvEmit->idx_of[24] >= 0) male=1;
    if (male){
      normal_copy[23] = 1;
      normal_copy[24] = 1;
    }
    maj_ncn = 2;//majority normal copy number
  }
  else{
    Clone::get_normal_copy(chr_fn);
  }
  //test normal copy numbers
  int found_err=0;
  if (cnaEmit->is_set){
    for (int s=0; s<cnaEmit->nSamples; s++){
      if (normal_copy[cnaEmit->chr[s]] == 0) found_err=1;
    }
  }
  if (bafEmit->is_set){
    for (int s=0; s<bafEmit->nSamples; s++){
      if (normal_copy[bafEmit->chr[s]] == 0) found_err=1;
    }
  }
  if (snvEmit->is_set){
    for (int s=0; s<snvEmit->nSamples; s++){
      if (normal_copy[snvEmit->chr[s]] == 0) found_err=1;
    }
  }
  if (found_err){
    cout<<"Normal copy number for some chromosomes not given."<<endl;
    cout<<"Use --chr [file] to set normal copy numbers."<<endl;
    exit(1);
  }
}

void Clone::get_normal_copy(const char * chr_fn){
  if (chr_fn==NULL) abort();
  ifstream ifs;
  string line;
  stringstream line_ss;
  ifs.open( chr_fn, ios::in);
  if (ifs.fail()){
    printf("ERROR: file %s cannot be opened.\n", chr_fn);
    exit(1);
  }
  int chr = 0, ncn=0;
  map<int,int> tally;
  while( ifs.good() ){
    line.clear();
    getline( ifs, line);
    if (line.empty()) break;
    if (line[0] == '#') continue;
    line_ss.clear();
    line_ss.str(line);
    line_ss >> chr >> ncn;
    if (chr<0 || chr >= 100) abort();
    normal_copy[chr] = ncn;
    if (tally.count(ncn) == 0){
      tally.insert(std::pair<int,int>(ncn,0));
    }
    if (cnaEmit->is_set){
      //if (chr > cnaEmit->maxchr || cnaEmit->idx_of[chr] < 0) continue;
      if ( cnaEmit->chrs.count(chr) == 0 ) continue;
      tally[ncn] += cnaEmit->nSites[cnaEmit->idx_of[chr]];
    }
    else if (snvEmit->is_set){
      //if (chr > snvEmit->maxchr || snvEmit->idx_of[chr] < 0) continue;
      if ( snvEmit->chrs.count(chr) == 0 ) continue;
      tally[ncn] += snvEmit->nSites[snvEmit->idx_of[chr]];
    }
  }
  ifs.close();
  int mx=0;
  maj_ncn=2;
  map<int,int>::iterator it;
  for (it=tally.begin(); it != tally.end(); it++){
    if (it->second > mx){
      maj_ncn = it->first;
      mx = it->second;
    }
  }
}




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
    bulk_prior[s] = gsl_matrix_calloc( snvEmit->nSites[s], snvEmit->gridSize+1);
    bulk_prior_mean[s] = new double [ snvEmit->nSites[s] ];
  }
  //allocate posterior...
  bulk_post       = new gsl_matrix ** [nTimes];
  bulk_post_mean  = new double ** [nTimes];
  for (int t=0; t<nTimes; t++){
    bulk_post[t] = new gsl_matrix * [snvEmit->nSamples];   
    bulk_post_mean[t] = new double * [snvEmit->nSamples];   
    for (int s=0; s<snvEmit->nSamples; s++){
      bulk_post[t][s] = gsl_matrix_calloc( snvEmit->nSites[s], snvEmit->gridSize+1);
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

//Destructor
Clone::~Clone(){
  if (allocated == 1){
    delete [] purity;
    for (int l=0; l<nLevels; l++) delete [] copynumber[l]; 
    delete [] copynumber; 
    for (int t=0; t<nTimes; t++) delete [] clone_spectrum[t];
    delete [] clone_spectrum;
    delete [] normal_copy;
    if (cnaEmit->is_set){
      for (int t=0; t<nTimes; t++){
	for (int s=0; s<cnaEmit->nSamples; s++){
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
  if (snv_prior_from_cna_baf_map != NULL){
    for (int cn=0; cn<=maxcn; cn++) gsl_matrix_free(snv_prior_from_cna_baf_map[cn]);
    delete [] snv_prior_from_cna_baf_map;
  }
  if (snv_prior_from_cna_map != NULL){
    gsl_matrix_free(snv_prior_from_cna_map);
  }
  if (TransMat_cna != NULL) gsl_matrix_free(TransMat_cna);
  if (TransMat_snv != NULL) gsl_matrix_free(TransMat_snv);
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
    if (cnaEmit->is_set){
      Clone::set_TransMat_cna();
    }
    if (bafEmit->is_set) 
      Clone::set_baf_prior_map();
    if (snvEmit->is_set){
      if (cnaEmit->is_set) Clone::set_snv_prior_map();
      Clone::set_TransMat_snv();
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
  if (cnaEmit->is_set){
    complexity = double(nTimes)*( double(nClones) + 1.0) + pow( double(maxcn+1), nClones);
    double size = double(cnaEmit->total_loci);
    if (bafEmit->is_set) size += double(bafEmit->total_loci);
    if (snvEmit->is_set) size += double(snvEmit->total_loci);
    complexity *= log(size);
  }
  //TODO: if maxcn_mask is used???
  //
  /* *** BAF only not supported ***
    else if (bafEmit->is_set){
    complexity = double(nTimes)*(double(nClones) + 1.0) + pow(double(maxcn+1), nClones);
    double size = double(bafEmit->total_loci);
    complexity *= log(size);
    }
  */
  else if (snvEmit->is_set){
    double * cnmax_freq = new double [maxcn+1];
    for (int cn=0; cn<=maxcn;cn++) cnmax_freq[cn]=0.0;
    double norm=0.0;
    if ( snvEmit->cnmax == NULL ){
      for (int s=0; s<snvEmit->nSamples; s++){
	cnmax_freq[ normal_copy[ snvEmit->chr[s] ] ] += double(snvEmit->nSites[s]);
	norm += double(snvEmit->nSites[s]);
      }  
    }
    else{
      for (int s=0; s<snvEmit->nSamples; s++){
	for (int evt=0; evt<snvEmit->nEvents[s];evt++){
	  cnmax_freq[snvEmit->cnmax[s][evt]] += 1.0;
	}
	norm += double(snvEmit->nEvents[s]);
      }
    }
    for (int cn=0; cn<=maxcn;cn++){
      cnmax_freq[cn] /= norm;
    }
    complexity = 0.0;
    for (int cn=1; cn<=maxcn;cn++){
      complexity += cnmax_freq[cn] * pow( double(cn+1), nClones);
      if (!snvEmit->connect && snvEmit->cnmax_seen.count(cn) == 1){
	complexity += double(cn+1);//d.o.f. of cn-prior
      }
    }
    if (!snvEmit->connect) complexity += 1.0;
    complexity += double(nTimes)*(double(nClones) + 1.0);
    complexity *= log(double(snvEmit->total_loci));
    delete [] cnmax_freq;
  }
}



void Clone::set_cn_prior_cna( gsl_vector * prior, int sample){
  double pncn = 0.5; //penalty for being  different from the normal copy number
  double pdif = 0.01;//penalty for having different copynumbers in clones
  if (cnaEmit->is_set == 0) abort();
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
      prior->data[l] = prior->data[l]> 0.0 ? log(prior->data[l]) : logzero;
    }
  }
}


//SNV-mode only...
void Clone::initialize_cn_prior_snv(){// SNV prior, conditional on max-cn
  if (init_cn_prior_snv != NULL) gsl_matrix_free(init_cn_prior_snv);
  init_cn_prior_snv = gsl_matrix_calloc(maxcn+1,maxcn+1);
  gsl_matrix_set( init_cn_prior_snv, 0, 0, snv_fpr);//snv_fpr: hopefully small (default:1.0e-4)
  int nprior = 1;
  double p   = 0.5;// initial penalty for higher genotypes
  for (int cn=1; cn <= maxcn; cn++){
    gsl_vector_view subrow = gsl_matrix_subrow( init_cn_prior_snv, cn, 0, cn+1);
    if (!snvEmit->connect ){
      if( snvEmit->cnmax_seen.count(cn) != 0){//un-correlated SNV data
	gsl_vector_set( &subrow.vector, 0, p);
	for (int i=1; i<=cn; i++) gsl_vector_set( &subrow.vector, i, pow( p, i));
	gsl_vector_scale( &subrow.vector, 1.0 / gsl_blas_dasum(&subrow.vector) );
      }
    }
    else{//flat prior for correlated SNV data (Transition matrix will be used)
      gsl_vector_set_all( &subrow.vector, 1.0/double(cn+1));
    }
    nprior++;
  }
  printf("Using these SNV copynumber priors per clone (combinations have multiplicative priors)\n");
  for (int i=0; i<(int) init_cn_prior_snv->size1; i++){
    for (int j=0; j<(int) init_cn_prior_snv->size2; j++){
      printf("%.3f ", gsl_matrix_get( init_cn_prior_snv, i, j));
    }
    cout<<endl;
  }
  //set above fixed priors
  Clone::set_cn_prior_snv(init_cn_prior_snv);
}



//SNV-mode only...
void Clone::set_cn_prior_snv( gsl_matrix * prior_per_clone){
  if (snvEmit->is_set == 0) abort();
  // delete old ones...
  map<int,gsl_vector*>::iterator it;
  if ( !cn_prior_snv.empty() ){
    for (it=cn_prior_snv.begin(); it != cn_prior_snv.end(); ++it){
      gsl_vector_free(cn_prior_snv[it->first]);
    }
    cn_prior_snv.clear();
  }
  //set new ones...
  //cn==0
  gsl_vector * pr = gsl_vector_calloc(nLevels); 
  gsl_vector_set_all( pr, 0.0);
  gsl_vector_set( pr, 0, 1.0);
  if (snvEmit->log_space){//log-transform?
    for (int l=0; l<nLevels; l++){
      pr->data[l] = (pr->data[l] > 0.0) ? log(pr->data[l]) : logzero;
    }
  }
  cn_prior_snv.insert(pair<int,gsl_vector*>(0,pr));
  //cn>0
  for (int cn=1; cn<=maxcn; cn++){//if total c.n. = cn...
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
  if ( snv_prior_from_cna_baf_map == NULL){//allocate
    snv_prior_from_cna_baf_map = new gsl_matrix * [maxcn+1];
    for (int cn=0; cn <= maxcn; cn++){ 
      snv_prior_from_cna_baf_map[cn] = gsl_matrix_alloc( maxcn+1, maxcn+1);
    }
  }
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
  //via CNA posterior only
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
  Clone::apply_snv_prpc(prior,snv_prpc);
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
  Clone::apply_snv_prpc( prior, snv_prpc);
  gsl_matrix_free(snv_prpc);
  gsl_vector_free(cnapostpc);
  gsl_vector_free(bafpostpc);
  gsl_vector_free(bafppc);
}
  
//CNA (+BAF) + SNV mode...
void Clone::apply_snv_prpc( gsl_vector * prior, gsl_matrix * snv_prpc){
  gsl_vector_set_all( prior, 1.0);
  int level=0;
  if (snvEmit->connect == 0){
    prior->data[0] = snv_fpr;//SNV false positive rate
    level = 1;
  }
  while (level<nLevels){
    for (int j=0; j<nClones; j++){
      prior->data[level] *= gsl_matrix_get( snv_prpc, j, copynumber[level][j]);
    }
    level++;
  }
  //normalize...
  double norm = gsl_blas_dasum(prior);
  if (norm <=0.0 || norm!= norm) abort();
  gsl_vector_scale(prior,1.0/norm);
  if (snvEmit->log_space){//log-transform?
    for (int l=0; l<nLevels; l++){
      prior->data[l] = prior->data[l] > 0.0 ? log(prior->data[l]) : - 1.0e10;
    }
  }
}


//probability weight in chromosomes with majority normal copy number
void Clone::get_cna_marginals(){
  if (cnaEmit->is_set == 0){
    cout<<"ERROR-1 in Clone::get_cna_marginals()\n";
    exit(1);
  }
  if (marginals != NULL) gsl_matrix_free(marginals);
  int ct = 0;
  for (int f=0; f<nFreq; f++) ct += maxcn+1;
  marginals = gsl_matrix_calloc( cnaEmit->nSamples, ct);
  gsl_vector * post = gsl_vector_calloc(nLevels);
  gsl_vector * mem  = gsl_vector_calloc(nLevels);
  if (copynumber_post != NULL) gsl_matrix_free(copynumber_post);
  if (majcn_post != NULL) gsl_vector_free(majcn_post);
  copynumber_post = gsl_matrix_calloc( cnaEmit->nSamples, nLevels);
  majcn_post      = gsl_vector_calloc(nLevels);
  alpha_cna = new gsl_matrix * [cnaEmit->nSamples];
  gamma_cna = new gsl_matrix * [cnaEmit->nSamples];
  for (int s=0; s< cnaEmit->nSamples; s++){
    gsl_vector_view cn_row = gsl_matrix_row(copynumber_post,s);
    gsl_vector_view mg_row = gsl_matrix_row(marginals,s);
    alpha_cna[s] = NULL;
    gamma_cna[s] = NULL;
    //get the posterior here
    Clone::get_cna_posterior(s);
    int ncn = normal_copy[ cnaEmit->chr[s] ];
    int evt=0,last_evt=-1;
    for (int idx=0; idx < cnaEmit->nSites[s]; idx++){
      evt = cnaEmit->event_of_idx[s][idx];
      if (evt != last_evt){
	gsl_matrix_get_row( post, gamma_cna[s], evt);
	last_evt = evt;
      }
      if ( ncn == maj_ncn) gsl_vector_add( majcn_post, post);
      gsl_vector_add( &cn_row.vector, post);
      gsl_blas_dgemv( CblasNoTrans, 1.0, margin_map, post, 1.0, &mg_row.vector);
    }
    gsl_matrix_free(gamma_cna[s]);
    gamma_cna[s] = NULL;
    gsl_vector_add( mem, &cn_row.vector);
  }
  delete [] alpha_cna;
  delete [] gamma_cna;
  alpha_cna = NULL;
  gamma_cna = NULL;
  gsl_vector_free(post);
  //normalize majority ncn posterior
  double norm = gsl_blas_dasum( majcn_post);
  if (norm <= 0.0) abort();
  gsl_vector_scale( majcn_post, 1.0/norm);
  //normalize posterior
  norm = gsl_blas_dasum(mem);
  if (norm <= 0.0) abort();
  for (int s=0; s<cnaEmit->nSamples; s++){
    gsl_vector_view cn_row = gsl_matrix_row(copynumber_post,s);
    gsl_vector_scale(&cn_row.vector,1.0/norm);
  }
  gsl_vector_free(mem);
  //normalize marginals
  gsl_vector_view part;
  ct = 0;
  for (int f=0; f<nClones; f++){
    norm = 0.0;
    for (int s=0; s<cnaEmit->nSamples; s++){
      part = gsl_matrix_subrow( marginals, s, ct, maxcn+1);
      norm += gsl_blas_dasum(&part.vector);
    }
    gsl_matrix_view pt = gsl_matrix_submatrix( marginals, 0, ct, cnaEmit->nSamples, maxcn+1);
    if (norm <= 0.0) abort();
    gsl_matrix_scale( &pt.matrix, 1.0 / norm);
    ct += maxcn+1;
  }
}

//only in samples where normal copy number (ncn) == majority ncn
void Clone::get_mass_candidates(){
  Clone::get_cna_marginals();
  levels_sorted.clear();
  for (int i=0; i<nLevels; i++) levels_sorted.push_back(i);
  SortDesc sd;
  sd.arg = majcn_post->data;
  std::sort( levels_sorted.begin(), levels_sorted.end(), sd);
  if (mass_candidates != NULL) gsl_matrix_free(mass_candidates);
  mass_candidates = gsl_matrix_calloc(nLevels,nTimes);
  double z;
  for (int i=0; i<nLevels; i++){
    int level = levels_sorted[i];
    for (int t=0; t<nTimes; t++){
      z = double(maj_ncn) * (1.0 - purity[t]) + clone_spectrum[t][level];
      z *= mass->data[t] / double(maj_ncn);
      gsl_matrix_set( mass_candidates, i, t, z);
    }
  }
}




void Clone::set_logn(){//only needed for cnaEmit
  for (int t=0; t<nTimes; t++){
    for (int s=0; s<cnaEmit->nSamples; s++){
      for (int l=0; l<cnaEmit->nSites[s]; l++){
	unsigned int N = cnaEmit->depths[t][s][l];
	if (N==0) continue;//no observation
	unsigned int n = cnaEmit->reads[t][s][l];
	if( logn.count(n) == 0 ){
	  logn[n]   = log(double(n));
	  loggma[n+1] = gsl_sf_lngamma(double(n+1));
	}
	if( logn.count(N) == 0 ){
	  logn[N]   = log(double(N));
	  loggma[N+1] = gsl_sf_lngamma(double(N+1));
	}
	if (cnaEmit->mode == 4){//for negative binomial model
	  int n1 = int(cnaEmit->shape*double(N));
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
  for (int s=0; s<cnaEmit->nSamples; s++){
    int ncn = normal_copy[cnaEmit->chr[s]];
    if( logn.count(ncn) == 0 ){
      logn[ncn] = log(double(ncn));
    }
  }
  logn_set = 1;
}



void Clone::get_nmean(){//only needed for cnaEmit
  cnaEmit->minRate = 1.0e6;
  cnaEmit->maxRate = 0;
  for (int t=0; t<nTimes; t++){
    int ct = 0;
    double mn = 0;
    for (int s=0; s<cnaEmit->nSamples; s++){
      for (int l=0; l< cnaEmit->nSites[s]; l++){
	unsigned int n = cnaEmit->reads[t][s][l];
	unsigned int N = cnaEmit->depths[t][s][l];
	if (N==0) continue;
	if (N>0 && n>0){
	  cnaEmit->minRate = min(cnaEmit->minRate, double(n)/double(N));
	  cnaEmit->maxRate = max(cnaEmit->maxRate, double(n)/double(N));
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
void Clone::set_TransMat_cna(){
  if (TransMat_cna != NULL) gsl_matrix_free(TransMat_cna);
  TransMat_cna = gsl_matrix_calloc(nLevels,nLevels);
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
	gsl_matrix_set(TransMat_cna,i, j, 1.0);
      }
      else{
	gsl_matrix_set(TransMat_cna,i, j, 0.0);
      }
    }
    row  = gsl_matrix_row(TransMat_cna,i);
    norm = gsl_blas_dasum(&row.vector);
    gsl_vector_scale(&row.vector,1.0/norm);
  }
}




// only one clone can change its state,
void Clone::set_TransMat_snv(){
  if (TransMat_snv != NULL) gsl_matrix_free(TransMat_snv);
  TransMat_snv = gsl_matrix_calloc(nLevels,nLevels);
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
	gsl_matrix_set(TransMat_snv,i, j, 1.0);
      }
      else{
	gsl_matrix_set(TransMat_snv,i, j, 0.0);
      }
    }
    row  = gsl_matrix_row(TransMat_snv,i);
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
  save_snv_alpha = 0;
  total_llh     = 0.0;
  cna_total_llh = 0.0;
  baf_total_llh = 0.0;
  snv_total_llh = 0.0;
  int sample;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
  for ( sample=0; sample < cnaEmit->nSamples; sample++){
    double llh,ent;   
    Clone::do_cna_Fwd( sample, llh);
#pragma omp critical
    {
      cna_total_llh += llh;
    }
    if (save_cna_alpha==1){
      int cnaChr = cnaEmit->chr[sample];
      Clone::do_cna_Bwd( sample, ent);
      Clone::get_phi(sample);
      if ( bafEmit->is_set && bafEmit->chrs.count(cnaChr) == 1){
	Clone::map_phi( cnaEmit, sample, bafEmit);	
	if ( snvEmit->is_set && snvEmit->chrs.count(cnaChr) == 1){
	  int bafsample = bafEmit->idx_of[cnaEmit->chr[sample]];
	  Clone::map_phi( bafEmit, bafsample, snvEmit);
	}
      }
      else if ( snvEmit->is_set && snvEmit->chrs.count(cnaChr) == 1){
	Clone::map_phi( cnaEmit, sample, snvEmit);
      }
      gsl_matrix_free(alpha_cna[sample]);
      alpha_cna[sample] = NULL;
    }
  }//END PARALLEL FOR
  //
  // BAF
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
  // SNV
  if( snvEmit->is_set ){
#pragma omp parallel for schedule( dynamic, 1) default(shared)
    for ( sample=0; sample < snvEmit->nSamples; sample++){//START PARALLEL FOR
      double llh = 0.0;
      Clone::do_snv_Fwd(sample, llh);
#pragma omp critical
      {
	snv_total_llh += llh;
      }
    }//END PARALLEL FOR
  }
  //cleanup...
  if (save_cna_alpha==1){
    for ( sample=0; sample < cnaEmit->nSamples; sample++){
      gsl_matrix_free(gamma_cna[sample]);
    }
    delete [] alpha_cna;
    delete [] gamma_cna;
    alpha_cna = NULL;
    gamma_cna = NULL;
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
  total_llh = cna_total_llh + baf_total_llh + snv_total_llh;
  return(total_llh);
}

double Clone::get_cna_total_llh(){
  int sample;
  save_cna_alpha = 0;
  cna_total_llh  = 0.0;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
  for ( sample=0; sample < cnaEmit->nSamples; sample++){
    double llh;
    Clone::do_cna_Fwd( sample, llh);
#pragma omp critical
    {
      cna_total_llh += llh;
    }
  }//END PARALLEL FOR
  return(cna_total_llh);
}


double Clone::get_baf_total_llh(){
  if ( nClones > 0 && cnaEmit->is_set && gamma_cna == NULL ) abort();
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


double Clone::get_snv_total_llh(){
  if ( nClones > 0 && cnaEmit->is_set && gamma_cna == NULL ) abort();
  int sample;
  save_snv_alpha = 0;
  snv_total_llh  = 0.0;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
  for ( sample=0; sample< snvEmit->nSamples; sample++){
    double llh;
    Clone::do_snv_Fwd( sample, llh);
#pragma omp critical
    {
      snv_total_llh += llh;
    }
  }
  return(snv_total_llh);
}




double Clone::get_cna_posterior(int sample){
  double llh=0, entropy=0;
  save_cna_alpha = 1;
  //set fw-bw arrays    
  if (alpha_cna[sample] == NULL) 
    alpha_cna[sample] = gsl_matrix_calloc( cnaEmit->nEvents[sample], nLevels);
  if (gamma_cna[sample] == NULL) 
    gamma_cna[sample] = gsl_matrix_calloc( cnaEmit->nEvents[sample], nLevels);
  Clone::do_cna_Fwd( sample, llh);
  Clone::do_cna_Bwd( sample, entropy);
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


double Clone::get_snv_posterior(int sample){
  if (cnaEmit->is_set){
    int cna_sample = cnaEmit->idx_of[snvEmit->chr[sample]];
    if ( nClones > 0 && (gamma_cna == NULL || gamma_cna[cna_sample] == NULL)){
      abort();
    }
  }
  double llh=0, ent=0;
  save_snv_alpha = 1;
  //set fw-bw arrays    
  if (alpha_snv[sample] == NULL) 
    alpha_snv[sample] = gsl_matrix_calloc( snvEmit->nEvents[sample], nLevels);
  if (gamma_snv[sample] == NULL) 
    gamma_snv[sample] = gsl_matrix_calloc( snvEmit->nEvents[sample], nLevels);
  Clone::do_snv_Fwd( sample, llh);
  Clone::do_snv_Bwd( sample, ent);
  gsl_matrix_free(alpha_snv[sample]);
  alpha_snv[sample] = NULL;
  return(llh);
}


void Clone::do_cna_Fwd( int sample, double& llh){
  Clone::set_tcn(sample);//pre-compute total copy number
  gsl_vector * prior    = gsl_vector_alloc(nLevels);
  gsl_vector * post     = gsl_vector_alloc(nLevels);
  gsl_matrix * Trans = NULL;
  double flat = (cnaEmit->log_space) ? -log(double(nLevels)) : 1.0/double(nLevels);
  gsl_vector_set_all(prior, flat);
  gsl_vector_set_all(post,  flat);  
  if (nClones>0){//preparations...
    Clone::set_cn_prior_cna( prior, sample);
    if (cnaEmit->connect){
      Trans = gsl_matrix_alloc( nLevels, nLevels);
      gsl_matrix_memcpy( Trans, TransMat_cna);
    }
  }
  else{
    gsl_vector_set_all(prior,flat);
  }
  gsl_vector_memcpy( post, prior);
  double norm = 0.0;
  double pj = 1.0;
  int idx=0;
  int cnaChr = cnaEmit->chr[sample];
  int mx = maxcn_mask.count(cnaChr) ? maxcn_mask[cnaChr]: maxcn;
  llh = 0.0;
  for (int evt=0; evt < cnaEmit->nEvents[sample]; evt++){
    idx = cnaEmit->idx_of_event[sample][evt];
    //***PREDICT STEP***
    if (nClones > 0 && cnaEmit->connect && evt > 0){//connect to the left...
      pj = cnaEmit->pjump[sample][idx];
      if (mx < maxcn){//some levels shut down
	Clone::predict( prior, post, cnaEmit, pj, Trans, mx);
      }
      else{//all levels open
	Clone::predict( prior, post, cnaEmit, pj, Trans);
      }
    }//connect done.
    else if (nClones == 0){
      gsl_vector_set_all( prior, flat);
    }
    //***UPDATE***
    norm = Clone::update( prior, post, cnaEmit, sample, evt);
    llh += norm;
    if (save_cna_alpha == 1) gsl_matrix_set_row( alpha_cna[sample], evt, post);
  }
  // cleanup    
  gsl_vector_free(prior);
  gsl_vector_free(post);
  if (Trans != NULL) gsl_matrix_free(Trans);
}



void Clone::do_baf_Fwd( int sample, double& llh){
  int cna_sample = 0;
  if (cnaEmit->is_set){
    cna_sample = cnaEmit->idx_of[bafEmit->chr[sample]];
    if (nClones > 0){
      if (gamma_cna == NULL) abort();
      if (gamma_cna[cna_sample] == NULL) abort();
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
  int idx=0, cna_evt=0, last_cna_evt =-1, nidx=0;
  gsl_vector_view cna_post;
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
      if (cnaEmit->is_set){//connect with CNA	
	cna_evt = bafEmit->Event_of_idx[sample][idx];
	if (cna_evt != last_cna_evt){
	  cna_post = gsl_matrix_row( gamma_cna[cna_sample], cna_evt);
	  get_baf_prior_from_cna_post( cn_prior, &cna_post.vector);
	  last_cna_evt = cna_evt;
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


void Clone::do_snv_Fwd(int sample, double& llh){
  int snvChr = snvEmit->chr[sample];
  int cnaSample = -1;
  if (cnaEmit->is_set){
    if (cnaEmit->chrs.count(snvChr) == 0) abort();
    cnaSample = cnaEmit->idx_of[snvEmit->chr[sample]];
  }
  int bafSample = -1;  
  if (bafEmit->is_set && bafEmit->chrs.count(snvChr) == 1){
    bafSample = bafEmit->idx_of[snvEmit->chr[sample]];
  }
  if ( nClones >0 ){
    if (cnaEmit->is_set)//need CNA post
      if( gamma_cna == NULL || gamma_cna[cnaSample] == NULL) abort();
    if (bafEmit->is_set && bafSample >=0 )//also need BAF post
      if( gamma_baf == NULL || gamma_baf[bafSample] == NULL) abort();
  }
  gsl_vector * cn_prior = gsl_vector_alloc(nLevels);
  gsl_vector * prior    = gsl_vector_alloc(nLevels);
  gsl_vector * post     = gsl_vector_alloc(nLevels);
  gsl_vector * mem      = gsl_vector_alloc(nLevels);
  double flat = (snvEmit->log_space) ? -log(double(nLevels)) : 1.0/double(nLevels);
  gsl_vector_set_all( cn_prior, flat);
  gsl_vector_set_all( prior,    flat);
  gsl_vector_set_all( post,     flat);
  gsl_matrix * Trans = NULL;
  gsl_vector_view cna_post, baf_post;
  gsl_vector_memcpy( post, prior);
  int cn=0, ocn=-1;
  if ( nClones>0 && snvEmit->connect == 1){
    Trans = gsl_matrix_alloc( nLevels, nLevels);
    gsl_matrix_memcpy( Trans, TransMat_snv);
  }
  int mx = maxcn_mask.count(snvChr) ? maxcn_mask[snvChr]: maxcn;
  double norm = 0.0;
  double pj   = 1.0;
  int ncn     = normal_copy[ snvEmit->chr[sample] ];
  int cna_evt=-1, last_cna_evt=-1;
  int baf_evt=-1, last_baf_evt=-1, baf_idx=0;
  int idx=0, nidx=0, last_evt=snvEmit->nEvents[sample]-1;
  llh = 0.0;
  for (int evt=0; evt <= last_evt; evt++){
    idx = snvEmit->idx_of_event[sample][evt];
    //***PREDICT STEP***
    if (nClones > 0){
      if (snvEmit->connect && evt > 0){//connect to the left...
	pj = snvEmit->pjump[sample][idx];
	if (mx<maxcn){
	  Clone::predict( prior, post, snvEmit, pj, Trans, mx);
	}
	else{
	  Clone::predict( prior, post, snvEmit, pj, Trans);
	}
      }//...connect to left done
      if (cnaEmit->is_set){//connect to CNA...
	if (bafEmit->is_set && bafSample >= 0){//use BAF post
	  baf_evt = snvEmit->Event_of_idx[sample][idx];
	  baf_idx = bafEmit->idx_of_event[bafSample][baf_evt];
	  cna_evt = bafEmit->Event_of_idx[bafSample][baf_idx];
	}
	else{
	  cna_evt = snvEmit->Event_of_idx[sample][idx];
	}
	if (cna_evt != last_cna_evt || baf_evt != last_baf_evt){//new segment, new prior
	  cna_post = gsl_matrix_row( gamma_cna[cnaSample], cna_evt);
	  if (bafEmit->is_set && bafSample >= 0){
	    baf_post = gsl_matrix_row( gamma_baf[bafSample], baf_evt);
	    get_snv_prior_from_cna_baf_post( cn_prior, &cna_post.vector, &baf_post.vector);
	  }
	  else{
	    get_snv_prior_from_cna_post( cn_prior, &cna_post.vector);
	  }
	  last_cna_evt = cna_evt;
	  last_baf_evt = baf_evt;
	}
	if (snvEmit->connect) gsl_vector_memcpy( mem, prior);
	gsl_vector_memcpy( prior, cn_prior);
	nidx = (evt < last_evt) ? snvEmit->idx_of_event[sample][evt+1] : snvEmit->nSites[sample];
	if ( nidx-idx > 1){//exponentiate prior for all observations in this block
	  gsl_vector_scale( prior, double(nidx-idx));//log-space!
	  norm = log_vector_norm(prior);
	  gsl_vector_add_constant(prior,-norm);
	}  
	// multiply the two priors and rescale...
	if (snvEmit->connect){
	  if (snvEmit->log_space){
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
      }//connect to CNA done.
      else if (!cnaEmit->is_set && !snvEmit->connect){//cna not set, no correlations...
	cn = (snvEmit->cnmax == NULL) ? ncn : snvEmit->cnmax[sample][evt];
	if (cn != ocn){
	  gsl_vector_memcpy( cn_prior, cn_prior_snv[cn]);
	  ocn = cn;
	}//or flat over all levels?
	gsl_vector_memcpy( prior, cn_prior);
      }
    }
    else if (nClones == 0){//just one element!
      gsl_vector_set_all( prior, flat);
    }
    //***UPDATE STEP***
    norm = Clone::update( prior, post, snvEmit, sample, evt);
    llh += norm;
    if (save_snv_alpha == 1) gsl_matrix_set_row( alpha_snv[sample], evt, post);
  }
  // cleanup    
  gsl_vector_free(mem);
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(cn_prior);
  if (Trans != NULL) gsl_matrix_free(Trans);
}



void Clone::do_cna_Bwd(int sample, double& ent){
  if (alpha_cna[sample] == NULL || gamma_cna[sample] == NULL) abort();
  gsl_vector * prior    = gsl_vector_alloc(nLevels);
  gsl_vector * post     = gsl_vector_alloc(nLevels);
  gsl_vector * mem      = gsl_vector_alloc(nLevels);
  double flat = (cnaEmit->log_space) ? -log(double(nLevels)) : 1.0/double(nLevels);
  gsl_vector_set_all( prior, flat);
  gsl_vector_set_all( post,  flat);
  gsl_matrix * Trans = NULL;
  ent = 0.0;
  if (nClones>0){
    Clone::set_cn_prior_cna( prior, sample);
    if (cnaEmit->connect){
      Trans = gsl_matrix_alloc(nLevels,nLevels);
      gsl_matrix_memcpy( Trans, TransMat_cna);
    }
  }
  gsl_vector_memcpy( post, prior);
  int cnaChr = cnaEmit->chr[sample];
  int mx = maxcn_mask.count(cnaChr) ? maxcn_mask[cnaChr]: maxcn;
  gsl_vector_view alph;
  double pj = 1.0;
  double norm;
  int idx = 0;
  int last_evt = cnaEmit->nEvents[sample]-1;
  int last_idx = cnaEmit->idx_of_event[sample][last_evt];
  for (int evt = last_evt; evt >= 0; evt--){
    idx = cnaEmit->idx_of_event[sample][evt];
    //***PREDICTION STEP***
    if ( cnaEmit->connect && nClones > 0 &&  evt < last_evt){//connect to the right... 
      pj = cnaEmit->pjump[sample][last_idx];
      last_idx = idx;
      if (mx < maxcn){//some levels shut down
	Clone::predict( prior, post, cnaEmit, pj, Trans, mx);
      }
      else{//all levels open
	Clone::predict( prior, post, cnaEmit, pj, Trans);
      }
    }//connect done.
    else{
      gsl_vector_set_all( prior, flat);
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
    //ent += Clone::entropy(mem);
    //***UPDATE STEP*** (normalization term not needed here)
    Clone::update( prior, post, cnaEmit, sample, evt);
  }
  // cleanup    
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(mem);
  if (Trans!= NULL) gsl_matrix_free(Trans);
}



void Clone::do_baf_Bwd( int sample, double& ent){
  if (alpha_baf[sample] == NULL || gamma_baf[sample] == NULL) abort();
  int cna_sample = 0;
  if (cnaEmit->is_set){
    cna_sample = cnaEmit->idx_of[bafEmit->chr[sample]];
    if ( nClones>0 ){
      if ( gamma_cna == NULL ) abort();
      if ( gamma_cna[cna_sample] == NULL) abort();
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
  gsl_vector_view cna_post;
  double norm;
  ent = 0.0;
  gsl_vector_view alph;
  double pj = 0.0;
  int last_evt = bafEmit->nEvents[sample] - 1;
  int last_idx = bafEmit->idx_of_event[sample][last_evt];
  int idx = 0,nidx=0;
  int cna_evt=0,last_cna_evt=-1;
  for (int evt = last_evt; evt >= 0; evt--){
    idx = bafEmit->idx_of_event[sample][evt];
    //***PREDICTION STEP***
    if (nClones > 0 ){ 
      if (bafEmit->connect && evt < last_evt){//connect to the right...
	pj = bafEmit->pjump[sample][last_idx];
	last_idx = idx;
	Clone::predict( prior, post, bafEmit, pj, flat);
      }//...connect done.   
      if (cnaEmit->is_set){//connect with CNA...	
	cna_evt = bafEmit->Event_of_idx[sample][idx];
	if (cna_evt != last_cna_evt){
	  cna_post = gsl_matrix_row( gamma_cna[cna_sample], cna_evt);
	  get_baf_prior_from_cna_post( cn_prior, &cna_post.vector);
	  last_cna_evt = cna_evt;
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
      }//...connect with CNA done.
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




void Clone::do_snv_Bwd( int sample, double& ent){
  if (alpha_snv[sample] == NULL || gamma_snv[sample] == NULL) abort();
  int snvChr = snvEmit->chr[sample];
  int cnaSample = -1;
  if (cnaEmit->is_set){
    if (cnaEmit->chrs.count(snvChr) == 0) abort();
    cnaSample = cnaEmit->idx_of[snvChr];
  }
  int bafSample = -1;
  if (bafEmit->is_set && bafEmit->chrs.count(snvChr) == 1){
    bafSample = bafEmit->idx_of[snvChr];
  }
  if ( nClones >0 ){
    if (cnaEmit->is_set)
      if( gamma_cna == NULL || gamma_cna[cnaSample] == NULL) abort();
    if (bafEmit->is_set && bafSample >= 0)
      if( gamma_baf == NULL || gamma_baf[bafSample] == NULL) abort();
  }
  gsl_vector * cn_prior = gsl_vector_alloc(nLevels);
  gsl_vector * prior    = gsl_vector_alloc(nLevels);
  gsl_vector * post     = gsl_vector_alloc(nLevels);
  gsl_vector * mem      = gsl_vector_alloc(nLevels);
  double flat = (snvEmit->log_space) ? -log(double(nLevels)) : 1.0/double(nLevels);
  gsl_vector_set_all(prior, flat);
  gsl_vector_set_all(cn_prior, flat);
  gsl_vector_set_all(post,  flat);
  gsl_matrix * Trans = NULL;
  ent = 0.0;
  if (nClones>0 && snvEmit->connect){
    Trans = gsl_matrix_alloc( nLevels, nLevels);
    gsl_matrix_memcpy( Trans, TransMat_snv);
  }
  int mx = maxcn_mask.count(snvChr) ? maxcn_mask[snvChr]: maxcn;
  gsl_vector_view alph;
  gsl_vector_view cna_post,baf_post;
  double pj = 1.0, norm=0;
  int ncn = normal_copy[snvEmit->chr[sample]];
  int cn=0, ocn=-1;
  int cna_evt=-1, last_cna_evt=-1;
  int baf_evt=-1, last_baf_evt=-1, baf_idx=0;
  int idx=0, nidx=0;
  int last_evt = snvEmit->nEvents[sample]-1;
  int last_idx = snvEmit->idx_of_event[sample][last_evt];
  for (int evt = last_evt; evt >= 0 ; evt--){
    idx = snvEmit->idx_of_event[sample][evt];
    //***PREDICTION STEP***
    if (nClones > 0){
      if ( snvEmit->connect && evt < last_evt){//connect to the right...
	pj = snvEmit->pjump[sample][last_idx];
	last_idx = idx;
	if (mx<maxcn){
	  Clone::predict( prior, post, snvEmit, pj, Trans, mx);
	}
	else {
	  Clone::predict( prior, post, snvEmit, pj, Trans);
	}
      }//...connect to right done
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
	    get_snv_prior_from_cna_baf_post( cn_prior, &cna_post.vector, &baf_post.vector);
	  }
	  else{
	    get_snv_prior_from_cna_post( cn_prior, &cna_post.vector);
	  }
	  last_cna_evt = cna_evt;
	  last_baf_evt = baf_evt;
	}
	if (snvEmit->connect) gsl_vector_memcpy( mem, prior);
	gsl_vector_memcpy( prior, cn_prior);
	nidx = (evt < last_evt) ? snvEmit->idx_of_event[sample][evt+1] : snvEmit->nSites[sample];
	if ( nidx-idx > 1){//exponentiate prior for all observations in this block
	  gsl_vector_scale( prior, double(nidx-idx));//log-space!
	  norm = log_vector_norm(prior);
	  gsl_vector_add_constant(prior,-norm);
	}  
	// multiply two priors and rescale...
	if (snvEmit->connect){
	  if (snvEmit->log_space){
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
      }//connect to CNA done.
      else if (!cnaEmit->is_set && !snvEmit->connect){//CNA not set, no correlations...
	cn = (snvEmit->cnmax == NULL) ? ncn : snvEmit->cnmax[sample][evt];
	if (cn != ocn){
	  gsl_vector_memcpy( cn_prior, cn_prior_snv[cn]);
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
    //ent += Clone::entropy(mem);
    //***UPDATE STEP*** (normalization term not needed here)
    Clone::update( prior, post, snvEmit, sample, evt);
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

//predict step with transition matrix and maximum cn constraint...
void Clone::predict( gsl_vector * prior, gsl_vector * post, Emission * myEmit, double pj, gsl_matrix * T, int mxcn){
  if (pj == 0.0){
    gsl_vector_memcpy( prior, post);//no jump possible
  }
  else{
    if (myEmit->log_space) for (int l=0;l<nLevels;l++) post->data[l] = exp(post->data[l]);
    gsl_vector_memcpy( prior, post);
    gsl_blas_dgemv( CblasTrans, pj, T, post, 1.0-pj, prior);
    for (int l=0;l<nLevels;l++){
      for (int j=0; j<nClones; j++){
	if (copynumber[l][j] > mxcn){
	  prior->data[l] = 0.0;
	  break;
	}
      }
    }
    double norm = gsl_blas_dasum(prior);
    if (norm <= 0.0) abort();
    gsl_vector_scale(prior,1.0/norm);
    if (myEmit->log_space){
      for (int l=0;l<nLevels;l++){
	prior->data[l] = prior->data[l]>0.0 ? log(prior->data[l]) : logzero;
      }
    }
  }
}

// predict step with convex mixing...
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
  if (myEmit == cnaEmit){
    Clone::update_cna( prior, post, sample, evt);
  }
  else if (myEmit == bafEmit){
    Clone::update_baf( prior, post, sample, evt);
  }
  else if (myEmit == snvEmit){
    Clone::update_snv( prior, post, sample, evt);
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


//to be precomputed before each CNA-Fwd run!
void Clone::set_tcn(int s){//only used for cnaEmit...
  if (cnaEmit->is_set==0) abort();
  for (int t=0; t<nTimes; t++){
    int ncn = normal_copy[ cnaEmit->chr[s] ];
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





// update step for CNA...
void Clone::update_cna( gsl_vector * prior, gsl_vector * post, int sample, int evt){ 
  if ( cnaEmit->coarse_grained ){
    Clone::update_cna_event( prior, post, sample, evt);
  }
  else{
    if (nClones==0){
      Clone::update_cna_site_noclone( post, sample, evt);
    }
    else{
      Clone::update_cna_site_wclone( prior, post, sample, evt);
    }
  }
}


void Clone::update_cna_event( gsl_vector * prior, gsl_vector * post, int sample, int evt){
  gsl_vector_set_all( post, 0.0);//log-space!
  gsl_vector_view emit;
  double x,dx,val;
  for (int t=0; t<nTimes; t++){
    emit = gsl_matrix_row( cnaEmitLog[t][sample], evt);
    dx   = (cna_xmax[t] - cna_xmin[t]) / double(cnaEmit->gridSize);
    for (int level=0; level<nLevels; level++){
      if (prior->data[level] <= logzero) continue;
      x = tcn[t][sample][level] * mass->data[t];
      if ( x < cna_xmin[t] || x > cna_xmax[t]){//outside range...
	val = -1.0e10;
      }
      else{//inside range...
	val = Clone::get_interpolation( x, cna_xmin[t], cna_xmax[t], dx, &emit.vector);
      }
      post->data[level] += val;
    }
  }
}



void Clone::update_cna_site_noclone( gsl_vector * post, int sample, int site){//no-clone special case...
  gsl_vector_set_all( post, 1.0);
  double val=0,rate=0;
  int ncn = normal_copy[cnaEmit->chr[sample]];
  double b  =  (cnaEmit->bias != NULL) ? cnaEmit->bias[sample][site] : -1.0;
  double lb =  (cnaEmit->log_bias != NULL) ? cnaEmit->log_bias[sample][site] : 0.0;
  for (int time=0; time<nTimes; time++){
    int N = cnaEmit->depths[time][sample][site];
    int n = cnaEmit->reads[time][sample][site];   
    double rnd = cnaEmit->rnd_emit / (double(N)*(cnaEmit->maxRate - cnaEmit->minRate) + 1.0);
    double nrnd = (1.0 - cnaEmit->rnd_emit);
    if (cnaEmit->mode == 3){//Poisson
      val  = - mass->data[time] * double(N*ncn);
      if (b>0.0) val *= b;
      val -= loggma[n+1];
      rate = log_mass->data[time] + logn[N] + logn[ncn];
      if (b>0.0) rate += lb;
      val += double(n) * rate;    
    }
    else if (cnaEmit->mode == 4){//Negative Binomial
      double n1   = double(N)*cnaEmit->shape;//rounding (error?)
      double nrml = double(ncn);
      val  = loggma[int(n1)+n] - loggma[int(n1)] - loggma[n+1] + n1*cnaEmit->log_shape;
      rate = mass->data[time] * nrml;
      if (b>0.0) rate *= b;
      val += double(n) * log(rate);
      val -= (double(n) + n1) * log(rate + cnaEmit->shape);
      if (val != val) abort();
    }
    val = exp(val);
    if (rnd > 0.0) val = val*nrnd + rnd;
    post->data[0] *= val; 
  }
}

//normal case with one or more clones
void Clone::update_cna_site_wclone( gsl_vector * prior, gsl_vector * post, int sample, int site){
  gsl_vector_set_all( post, 1.0);//NOT logspace!
  double b  = (cnaEmit->bias != NULL) ? cnaEmit->bias[sample][site] : -1.0;
  double lb = (cnaEmit->log_bias != NULL) ? cnaEmit->log_bias[sample][site] : 0.0;
  for (int time=0; time<nTimes; time++){
    int N = cnaEmit->depths[time][sample][site];
    int n = cnaEmit->reads[time][sample][site];
    double rnd  = cnaEmit->rnd_emit / (double(N)*(cnaEmit->maxRate - cnaEmit->minRate) + 1.0);
    double nrnd = 1.0 - cnaEmit->rnd_emit;
    double val=0;
    if (cnaEmit->mode==3){//Poisson
      double pre1 = double(n) * (log_mass->data[time] + logn[N]) - loggma[n+1];
      double pre2 = mass->data[time] * double(N);
      if (b > 0.0){//modulation
	pre1 += double(n) * lb;
	pre2 *= b;
      }
      for (int l=0; l<nLevels; l++){
	if (prior->data[l] <= 0.0) continue;
	val = pre1 - pre2*tcn[time][sample][l] + double(n)*log_tcn[time][sample][l];
	val = exp(val);
	if (rnd > 0.0) val = val*nrnd + rnd;
	post->data[l] *= val;
      }
    }
    else if (cnaEmit->mode == 4){//Negative Binomial
      double n1 = double(N)*cnaEmit->shape;//rounding (error?)
      double pre = loggma[int(n1)+n] - loggma[int(n1)] - loggma[n+1] + n1*cnaEmit->log_shape;
      double rate;
      for (int l=0; l<nLevels; l++){
	if (prior->data[l] <= 0.0) continue;
	rate = mass->data[time] * tcn[time][sample][l];
	val  = pre +  double(n) * (log_mass->data[time] + log_tcn[time][sample][l]);
	if (b > 0.0){
	  rate *= b;
	  val += double(n)*lb;
	}
	val  = val - (double(n) + n1) * log(rate + cnaEmit->shape);
	val = exp(val);
	if (rnd > 0.0) val = val*nrnd + rnd;
	post->data[l] *= val;
      }
    }
  }
}


void Clone::update_baf( gsl_vector * prior, gsl_vector * post, int sample, int evt){
  if (bafEmit->coarse_grained){
    Clone::update_baf_event( prior, post, sample, evt);
  }
  else{
    Clone::update_baf_site( prior, post, sample, evt);
  }
}

void Clone::update_baf_event( gsl_vector * prior, gsl_vector * post, int sample, int evt){
  gsl_vector_set_all( post, 0.0);//log-space!
  double x,y,val,total_cn;
  double dx = bafEmit->dx;
  gsl_vector_view emit;
  for (int t=0; t<nTimes; t++){
    total_cn = (bafEmit->phi == NULL) ? 2.0 : bafEmit->phi[t][sample][evt];//only diploid chromosomes!
    y = (1.0 - purity[t]) / total_cn;
    emit = gsl_matrix_row( bafEmitLog[t][sample], evt);
    for (int level=0; level<nLevels; level++){
      if (prior->data[level] <= logzero) continue;
      x = y + clone_spectrum[t][level] / total_cn;
      if ( x < 0.0 || x > 1.0){//outside range...
	val = logzero;
      }
      else{//inside range...
	val = Clone::get_interpolation( x, 0.0, 1.0, dx, &emit.vector);
      }
      post->data[level] += val;
    }
  }
}


void Clone::update_baf_site( gsl_vector * prior, gsl_vector * post, int sample, int site){
  gsl_vector_set_all( post, 1.0);//not log-space!
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
      if (prior->data[level] <= 0.0) continue;
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


void Clone::update_snv( gsl_vector * prior, gsl_vector * post, int sample, int evt){
  if (snvEmit->coarse_grained){
    Clone::update_snv_event( prior, post, sample, evt);
  }
  else if ( bulk_fix >= 0.0 || (bulk_mean != NULL && bulk_dist == NULL)){
    Clone::update_snv_fixed( prior, post, sample, evt);
  }
  else{
    Clone::update_snv_nfixed( prior, post, sample, evt);
  }
}

void Clone::update_snv_event( gsl_vector * prior, gsl_vector * post, int sample, int evt){
  gsl_vector_set_all( post, 0.0);//log-space!
  double val,total_cn;
  double ncn = double(normal_copy[snvEmit->chr[sample]]);
  for (int t=0; t<nTimes; t++){
    total_cn = (snvEmit->phi == NULL) ? ncn : snvEmit->phi[t][sample][evt];
    double b = (1.0 - purity[t]) * ncn / total_cn;
    for (int level=0; level<nLevels; level++){
      if (prior->data[level] <= logzero) continue;
      double a = clone_spectrum[t][level] / total_cn;
      if (a > 1.0){//outside range...
	int idx  = snvEmit->idx_of_event[sample][evt];
	int nidx = (evt < snvEmit->nEvents[sample]-1) ? snvEmit->idx_of_event[sample][evt+1] : snvEmit->nSites[sample];
	val = logzero * double(nidx-idx);
      }
      else{//inside range...
	val = Clone::get_interpolation( a, 0.0, 1.0, b, 0.0, 1.0/bulk_min[t][sample][evt], snvEmitLog[t][sample][evt]);
      }
      post->data[level] += val;
    }
  }
}

//Dirac-Delta bulk prior distribution (maybe fixed at zero)...
void Clone::update_snv_fixed( gsl_vector * prior, gsl_vector * post, int sample, int site){
  gsl_vector_set_all( post, 1.0);
  int ncn = normal_copy[snvEmit->chr[sample]];
  double x,y,val=0;
  unsigned int N,n;
  int evt = snvEmit->event_of_idx[sample][site];
  double dx = snvEmit->dx;
  gsl_vector * emit = NULL;
  for (int t=0; t<nTimes; t++){
    N = snvEmit->depths[t][sample][site];
    n = snvEmit->reads[t][sample][site];
    if (N==0){
      continue;
    }
    emit = snvEmit->EmitLog[N][n];
    double rnd  = snvEmit->rnd_emit / double(N+1);
    double nrnd = 1.0 - snvEmit->rnd_emit;
    double total_cn = (snvEmit->phi == NULL) ? double(ncn) : snvEmit->phi[t][sample][evt];
    double bfix = (bulk_fix >= 0.0) ? bulk_fix : bulk_mean[t][sample][site];
    y =  bfix * (1.0-purity[t]) * double(ncn) / total_cn;//for cancer app, this is usually 0.0!
    for (int level=0; level<nLevels; level++){
      if (prior->data[level] <= 0.0) continue;
      x = y + clone_spectrum[t][level] / total_cn; 
      if (snv_err>0.0){//frequency of erroneous reads
	x = snv_err + (1.0-snv_err)*x;
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
}




//non-trivial case with one or more clones
void Clone::update_snv_nfixed( gsl_vector * prior, gsl_vector * post, int sample, int site){
  gsl_vector_set_all(post, 1.0);
  int ncn = normal_copy[snvEmit->chr[sample]];
  double val=0;
  int evt = snvEmit->event_of_idx[sample][site];
  for (int time=0; time<nTimes; time++){//product over time points
    unsigned int N = snvEmit->depths[time][sample][site];
    unsigned int n = snvEmit->reads[time][sample][site];
    gsl_vector * emit = snvEmit->EmitProb[N][n];
    double rnd  = snvEmit->rnd_emit / double(N+1);
    double nrnd = 1.0 - snvEmit->rnd_emit;
    //evaluate likelihood function at a finite grid of clonal frequencies
    int nPts = (nLevels <= 100 ) ? nLevels-1 : 100;
    //grid for linear interpolation//***TEST***
    gsl_vector * mem = gsl_vector_calloc(nPts+1);
    double total_cn = (snvEmit->phi == NULL) ? double(ncn) : snvEmit->phi[time][sample][evt];
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
	  if (snvEmit->rnd_emit > 0.0) val = val*(1.0-snvEmit->rnd_emit) + rnd;
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
      f1=logzero;
      f2=logzero;
    }
  }
  else{
    if (x<xmax){
      f1 = emit->data[gs-2];
      f2 = emit->data[gs-1];
      nu += 1.0;
    }
    else{
      f1=logzero;
      f2=logzero;
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



//only ever used on snvEmit!
double Clone::trapezoidal( gsl_vector * blk, double a, double b, gsl_vector * emit, int get_log){
  int olow = -1;
  double dx = snvEmit->dx;
  double dy = dx * b;
  double y = a - dy;
  double f1=0,f2=0,nu=0,pre=0;
  int low,high;
  int gs = snvEmit->gridSize;
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
  if (get_log==1) val = (val > 0.0) ? log(val) : logzero;
  return(val);
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



void Clone::get_cnaEmitLog(){//only ever used for BAF data
  if (cnaEmitLog != NULL) abort();
  cna_xmin = new double [nTimes];
  cna_xmax = new double [nTimes];
  cnaEmitLog = new gsl_matrix ** [nTimes];
  double lnrnd = (cnaEmit->rnd_emit>0.0) ? log(1.0-cnaEmit->rnd_emit) : 0.0;
  for (int t=0; t<nTimes; t++){
    cnaEmit->get_log = 1;
    cnaEmit->init_range(t);
    cna_xmin[t] = cnaEmit->xmin;
    cna_xmax[t] = cnaEmit->xmax;
    cnaEmit->set_EmitProb(t);
    cnaEmitLog[t] = new gsl_matrix * [cnaEmit->nSamples];
    for (int s=0; s<cnaEmit->nSamples; s++){
      cnaEmitLog[t][s] = gsl_matrix_calloc( cnaEmit->nEvents[s], cnaEmit->gridSize+1);
      int evt;
#pragma omp parallel for schedule( dynamic, 1) default(shared)      
      for ( evt=0; evt<cnaEmit->nEvents[s]; evt++){
	gsl_vector * mem = gsl_vector_calloc(cnaEmit->gridSize+1);
	gsl_vector * rnd = gsl_vector_calloc(cnaEmit->gridSize+1);
	unsigned int n,N;
	double r=0;
	gsl_vector_view emit = gsl_matrix_row( cnaEmitLog[t][s], evt);
	int idxi = cnaEmit->idx_of_event[s][evt];
	int idxf = (evt<cnaEmit->nEvents[s]-1) ? cnaEmit->idx_of_event[s][evt+1] - 1: cnaEmit->nSites[s]-1;
	for(int idx=idxi; idx<=idxf; idx++){
	  n = cnaEmit->reads[t][s][idx];
	  N = cnaEmit->depths[t][s][idx];
	  if (N>0){//observation?
	    if (cnaEmit->bias == NULL){
	      gsl_vector_memcpy( mem, cnaEmit->EmitLog[N][n]);
	    }
	    else{
	      cnaEmit->get_eprob_wBias( mem, cnaEmit->EmitLog[N][n], cnaEmit->bias[s][idx], n, N, 1);
	    }	  
	    if (cnaEmit->rnd_emit > 0.0){//random?
	      gsl_vector_add_constant( mem, lnrnd);
	      r = double(N)*(cnaEmit->maxRate - cnaEmit->minRate);
	      if (cnaEmit->bias != NULL) r *= cnaEmit->bias[s][idx];
	      r =  cnaEmit->rnd_emit / (r+1.0);
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





void Clone::get_snvEmitLog(){//only ever used for SNV data
  if (!snvEmit->is_set) abort();
  if (snvEmitLog==NULL){//first allocation
    snvEmitLog = new gsl_matrix *** [nTimes];
    for (int t=0;t<nTimes; t++){
      snvEmitLog[t] = new gsl_matrix ** [snvEmit->nSamples];
      for (int s=0; s<snvEmit->nSamples; s++){
	snvEmitLog[t][s] = new gsl_matrix * [snvEmit->nEvents[s]];
	for (int evt=0; evt< snvEmit->nEvents[s]; evt++){
	  snvEmitLog[t][s][evt] = gsl_matrix_calloc(snvEmit->gridSize+1,snvEmit->gridSize+1);
	}
      }
    }
  }
  Clone::get_bulk_min();
  double dx = snvEmit->dx;
  double lnrnd = log(1.0 - snvEmit->rnd_emit);
  for (int t=0; t<nTimes; t++){   
    for (int s=0; s<snvEmit->nSamples; s++){      
      int evt;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
      for (evt=0; evt<snvEmit->nEvents[s]; evt++){
	gsl_matrix * mem     = gsl_matrix_calloc(snvEmit->gridSize+1,snvEmit->gridSize+1);
	gsl_matrix * rnd_vec = gsl_matrix_calloc(snvEmit->gridSize+1,snvEmit->gridSize+1);
	gsl_matrix_set_zero(snvEmitLog[t][s][evt]);
	int idxi = snvEmit->idx_of_event[s][evt];
	int idxf = (evt<snvEmit->nEvents[s]-1) ? snvEmit->idx_of_event[s][evt+1] - 1: snvEmit->nSites[s]-1;
	for(int idx=idxi; idx<=idxf; idx++){
	  double val=0;
	  unsigned int N = snvEmit->depths[t][s][idx];
	  if (N==0) continue;//no observation
	  unsigned int n = snvEmit->reads[t][s][idx];
	  gsl_vector * emit = (bulk_prior==NULL) ? snvEmit->EmitLog[N][n] : snvEmit->EmitProb[N][n];
	  double lrnd  = (snvEmit->rnd_emit>0.0) ? log(snvEmit->rnd_emit/double(N+1)) : 0.0;	  
	  //evaluate likelihood function at a finite grid of clonal frequencies
	  double bfix = (bulk_fix >= 0.0) ? bulk_fix : bulk_mean[t][s][idx];
	  gsl_vector_view blk;
	  if (bulk_prior != NULL) blk = gsl_matrix_row( bulk_dist[t][s], idx);
	  for (int i=0; i<=snvEmit->gridSize; i++){
	    double a = double(i) / double(snvEmit->gridSize);
	    for (int j=0; j<=snvEmit->gridSize; j++){
	      double b = double(j) / (double(snvEmit->gridSize) * bulk_min[t][s][evt]);
	      if( b > (1.0-a) / bulk_min[t][s][evt]){
		val = logzero;
	      }
	      else if (bulk_prior==NULL){//bulk point estimate...	      
		double x =  a + bfix * b;
		if (x > 1.0){//outside range...
		  val = logzero;
		}
		else if (x>0.0 && x<1.0){//inside range...
		  val = Clone::get_interpolation( x, 0.0, 1.0, dx, emit);
		}
		else if ( x==1.0 ){//edges...
		  val = (n==N) ? 0.0 : logzero;
		}
		else if ( x==0.0 ){
		  val = (n==0) ? 0.0 : logzero;
		}
	      }
	      else{//bulk full distribution...		
		val = Clone::trapezoidal( &blk.vector, a, b, emit, 1);//log-space
	      }	
	      gsl_matrix_set( mem, i, j, val);
	    }
	  }
	  // mutiply into posterior...
	  if (snvEmit->rnd_emit > 0.0){
	    gsl_matrix_add_constant( mem, lnrnd);
	    gsl_matrix_set_all( rnd_vec, lrnd);
	    log_matrix_add( mem, rnd_vec);
	  }
	  gsl_matrix_add( snvEmitLog[t][s][evt], mem);
	}
	gsl_matrix_free(mem);
	gsl_matrix_free(rnd_vec);
      }
    }  
  }
}




// Bayesian update of the SNV bulk
void Clone::update_bulk(int sample){
  if (!snvEmit->is_set) abort();
  if (bulk_mean == NULL) abort();
  gsl_vector * bpost=NULL, *flat=NULL, *emit=NULL;
  gsl_vector_view bprior;
  bpost = gsl_vector_alloc(snvEmit->gridSize+1);
  flat  = gsl_vector_alloc(snvEmit->gridSize+1);
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
  int gS = snvEmit->gridSize;
  double ncn = (double) normal_copy[snvEmit->chr[sample]];
  int evt = snvEmit->event_of_idx[sample][idx];
  double total_cn = (snvEmit->phi==NULL) ? ncn : snvEmit->phi[time][sample][evt];
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
      val *= gsl_matrix_get( gamma_snv[sample], evt, j);
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
      bulk_min[t] = new double * [snvEmit->nSamples];
      for (int s=0; s<snvEmit->nSamples; s++){
	bulk_min[t][s] = new double [snvEmit->nEvents[s]];
      }
    }
  }
  for (int t=0; t<nTimes; t++){
    for (int s=0; s<snvEmit->nSamples; s++){
      for (int evt=0; evt < snvEmit->nEvents[s]; evt++){
	bulk_min[t][s][evt] = 1.1;
	int idxi = snvEmit->idx_of_event[s][evt];
	int idxf = (evt<snvEmit->nEvents[s]-1) ? snvEmit->idx_of_event[s][evt+1]-1 : snvEmit->nSites[s]-1;
	for (int idx=idxi; idx<=idxf; idx++){
	  bulk_min[t][s][evt] = min( bulk_min[t][s][evt], bulk_mean[t][s][idx]);
	}
      }
    }
  }
}
