//clone.cpp

//own headers...
#include "emission.h"
#include "common-functions.h"
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
  freqs      = NULL;
  purity     = NULL;
  min_purity = NULL;
  marginals  = NULL;
  margin_map = NULL;
  baf_prior_map  = NULL;
  snv_prior_from_cna_baf_map  = NULL;
  snv_prior_from_cna_map      = NULL;
  maxtcn    = 4;
  logn_set  = 0;
  nFreq     = 0;
  mass      = NULL;
  log_mass  = NULL;
  nmean     = NULL;
  got_gamma = 0;
  mass_candidates = NULL;
  copynumber       = NULL;
  copynumber_post  = NULL;
  initial_snv_prior_param= NULL;
  learn_priors = 0;
  cn_usage = NULL;
  clone_spectrum   = NULL;
  majcn_post  = NULL;
  bulk_fix  = -1.0;
  cna_pen_norm = 1.0;
  cna_pen_zero = 1.0;
  cna_pen_diff = 1.0;
  baf_pen_comp = 1.0;
  snv_pen_mult = 0.01;
  snv_pen_high = 0.5;
  snv_fpr = 1.0e-4;
  // pre-computed
  bafEmitLog = NULL;
  cnaEmitLog = NULL;
  snvEmitLog = NULL;
  TransMat_cna=NULL;
  TransMat_snv=NULL;
  // (fwd) posteriors
  alpha_cna=NULL;
  gamma_cna=NULL;
  alpha_baf=NULL;
  gamma_baf=NULL; 
  alpha_snv=NULL;
  gamma_snv=NULL;
  // llh's
  cna_total_llh=0;
  baf_total_llh=0;
  snv_total_llh=0;
  // bulk variables
  bulk_prior      = NULL;
  bulk_post       = NULL;
  bulk_dist       = NULL;
  bulk_prior_mean = NULL;
  bulk_post_mean  = NULL;
  bulk_mean       = NULL; 
  bulk_min = NULL;
  //grids
  cnaGrid  = 300;
  bafGrid  = 100;
  snvGrid  = 100;
  bulkGrid = 100;
  logzero = -1.0e6;
  cna_llhs = NULL;
  baf_llhs = NULL;
  snv_llhs = NULL;
  cna_gofs = NULL;
  baf_gofs = NULL;
  snv_gofs = NULL;
  get_gofs = 0;
  bafSymMap = NULL;
}



//Destructor
Clone::~Clone(){
  if (allocated == 1){
    delete [] purity;
    for (int l=0; l<nLevels; l++) delete [] copynumber[l]; 
    delete [] copynumber; 
    for (int t=0; t<nTimes; t++) delete [] clone_spectrum[t];
    delete [] clone_spectrum;
    normal_copy.clear();
    //
    if (mass!=NULL) gsl_vector_free(mass);
    if (log_mass!=NULL) gsl_vector_free(log_mass);
    if (nmean!=NULL) gsl_vector_free(nmean);
    if (min_purity!=NULL) gsl_vector_free(min_purity);
    delete [] cna_llhs;
    delete [] baf_llhs;
    delete [] snv_llhs;
    delete [] cna_gofs;
    delete [] baf_gofs;
    delete [] snv_gofs;
  }
  if (freqs != NULL) gsl_matrix_free(freqs);
  if (margin_map != NULL) gsl_matrix_free(margin_map);
  if (baf_prior_map != NULL) gsl_matrix_free(baf_prior_map);
  if (snv_prior_from_cna_baf_map != NULL){
    for (int cn=0; cn<=maxtcn; cn++) gsl_matrix_free(snv_prior_from_cna_baf_map[cn]);
    delete [] snv_prior_from_cna_baf_map;
  }
  if (snv_prior_from_cna_map != NULL){
    gsl_matrix_free(snv_prior_from_cna_map);
  }
  logn.clear();
  loggma.clear();
  if (initial_snv_prior_param != NULL) gsl_matrix_free(initial_snv_prior_param);
  snv_prior.clear();
}


//allocator...
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
  //all chr present
  chrs.clear();
  if (cnaEmit->is_set){
    for (int s=0;s<cnaEmit->nSamples;s++) 
      chrs.insert(cnaEmit->chr[s]);
  }
  if (bafEmit->is_set){
    for (int s=0;s<bafEmit->nSamples;s++)
      chrs.insert(bafEmit->chr[s]);
  }
  if (snvEmit->is_set){
    for (int s=0;s<snvEmit->nSamples;s++) 
      chrs.insert(snvEmit->chr[s]);
  }
  if (cnaEmit->is_set){
    mass     = gsl_vector_alloc(nTimes);
    log_mass = gsl_vector_alloc(nTimes);
    nmean    = gsl_vector_alloc(nTimes);
    Clone::set_logn();
    Clone::get_nmean();
    if (bafEmit->is_set){//CNA + BAF (+SNV)
      for (int s=0; s<bafEmit->nSamples; s++){//all BAF are mapped to CNA
	bafEmit->map_idx_to_Event( cnaEmit, s);
      }
      if (snvEmit->is_set){//if BAF exists, use it for SNVs in autosomes only! 
	for (int s=0; s<snvEmit->nSamples; s++){
	  int snvChr = snvEmit->chr[s];
	  if ( bafEmit->chrs.count(snvChr) == 1){//BAF exists -> autosome 
	    snvEmit->map_idx_to_Event( bafEmit, s);
	  }
	  else{//BAF does not exist -> hapl. sex chr -> use CNA
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
  //allocate tcn...  
  maxChr = *(chrs.rbegin());
  tcn     = new double ** [maxChr+1];
  log_tcn = new double ** [maxChr+1];
  for (int chr=0; chr<=maxChr; chr++){
    if (chrs.count(chr)){
      tcn[chr] = new double * [nTimes];
      log_tcn[chr] = new double * [nTimes];
      for (int t=0; t<nTimes; t++){    
	tcn[chr][t]  = NULL;
	log_tcn[chr][t] = NULL;
      }
    }
    else{
      tcn[chr] = NULL;
      log_tcn[chr] = NULL;
    }
  }
  cna_llhs = new double [nTimes];
  baf_llhs = new double [nTimes];
  snv_llhs = new double [nTimes];
  cna_gofs = new double [nTimes];
  baf_gofs = new double [nTimes];
  snv_gofs = new double [nTimes];
  allocated = 1;//done
}

void Clone::set_normal_copy(const char * chr_fn){
  normal_copy.clear();
  if (chr_fn == NULL){//if --chr [file] was not given
    for (int chr=1; chr<=23; chr++) normal_copy.insert(pair<int,int>(chr,2));
    normal_copy.insert(pair<int,int>(24,0));//female by default
    int male=0;//check whether there is a Y chromosome...
    if (cnaEmit->is_set && cnaEmit->chrs.count(24) == 1) male=1;
    if (snvEmit->is_set && snvEmit->chrs.count(24) == 1) male=1;
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
      if (normal_copy.count(cnaEmit->chr[s]) == 0) found_err=1;
    }
  }
  if (bafEmit->is_set){
    for (int s=0; s<bafEmit->nSamples; s++){
      if (normal_copy.count(bafEmit->chr[s]) == 0) found_err=1;
      if (normal_copy[bafEmit->chr[s]] != 2)       found_err=2;
    }
  }
  if (snvEmit->is_set){
    for (int s=0; s<snvEmit->nSamples; s++){
      if (normal_copy.count(snvEmit->chr[s]) == 0) found_err=1;
    }
  }
  if (found_err==1){
    cout<<"Normal copy number for some chromosomes not given."<<endl;
    cout<<"Use --chr [file] to set normal copy numbers."<<endl;
    exit(1);
  }
  else if (found_err==2){
    cout<<"Normal copy number for some BAF data chromosomes is not 2."<<endl;
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
    if (normal_copy.count(chr) == 1) abort();
    normal_copy.insert(pair<int,int>(chr, ncn));
    if (tally.count(ncn) == 0){
      tally.insert(std::pair<int,int>(ncn,0));
    }
    if (cnaEmit->is_set){
      if ( cnaEmit->chrs.count(chr) == 0 ) continue;
      tally[ncn] += cnaEmit->nSites[cnaEmit->idx_of[chr]];
    }
    else if (snvEmit->is_set){
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

// set frequency-dependent variables
void Clone::set(const gsl_matrix * freq){
  int isnew=0;
  if( freq == NULL){// no clone
    nFreq   = 0;
    nClones = 0;
    freqs = NULL;
    if( nLevels == 0 || nClones != 0) isnew=1;
  }// new number of clones?
  else if ( freq != NULL && (int) freq->size2 != nFreq ){
    nFreq = (int) freq->size2;
    nClones = nFreq;
    if (freqs != NULL) gsl_matrix_free(freqs);
    freqs = gsl_matrix_alloc(nTimes,nFreq);
    isnew=1;
  }
  if (isnew){
    // set all the copynumbers per level
    Clone::set_copynumbers();
    Clone::set_margin_map();
    Clone::set_maxtcn_per_clone();
    Clone::set_all_levels();   
    if (nClones>0 && cnaEmit->is_set){
      if (bafEmit->is_set)  Clone::set_baf_prior_map();
      if (snvEmit->is_set)  Clone::set_snv_prior_map();
    }  
    if (nClones>0){
      if (cnaEmit->is_set) Clone::set_TransMat_cna();
      if (snvEmit->is_set && snvEmit->connect) Clone::set_TransMat_snv();
    }
    Clone::allocate_tcn();
    if (nClones>0 && snvEmit->is_set && snvEmit->av_cn != NULL && !snvEmit->connect){
      Clone::allocate_cn_usage();
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
  Clone::set_tcn();
  if (nClones>0 && snvEmit->is_set && snvEmit->av_cn != NULL && !snvEmit->connect){
    Clone::set_cn_usage();
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
      nLevels *= maxtcn + 1;
    }
    int nl = nLevels;
    for (int f=0;f<nFreq;f++){
      nl = (int) double(nl) / double(maxtcn + 1);
      offset[f] = nl;
    }
    // get the copynumbers for each level
    copynumber = new int * [nLevels];
    for (int l=0; l<nLevels; l++){
      copynumber[l] = new int [nFreq]; 
    }
    for (int f=0; f<nFreq; f++){
      for (int l=0; l<nLevels; l++){
	copynumber[l][f] = int(l/offset[f]) % (maxtcn+1);
      }
    }
    delete [] offset;
  }
}


void Clone::set_all_levels(){
  level_of.clear();
  std::map<int, vector<int> >::iterator it;
  for ( it=maxtcn_per_clone.begin(); it != maxtcn_per_clone.end(); it++){
    int chr = it->first;
    if (nClones==0){
      level_of[chr] = 0;
    }
    else{
      int lev=0;
      for (int j=0; j<nClones;j++){
	lev += maxtcn_per_clone[chr][j] * (j<nClones-1 ? pow(maxtcn+1,nClones-1-j) : 1);
      }
      level_of[chr] = lev;
    }
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

void Clone::set_maxtcn_per_clone(){
  maxtcn_per_clone.clear();
  std::map<int, vector<int> >::iterator it;
  for (it=maxtcn_input.begin(); it != maxtcn_input.end(); it++){
    int chr = it->first;
    vector<int> mxpc;
    if (nClones==0){
      mxpc.push_back(0);
      maxtcn_per_clone.insert(std::pair< int, vector<int> >(chr,mxpc));
    }
    else{
      for (int j=0;j<nClones;j++){
	if ( j < (int) maxtcn_input[chr].size()){
	  mxpc.push_back( maxtcn_input[chr][j] );
	}
	else{
	  mxpc.push_back( maxtcn_input[chr].back() );
	}
      }
      maxtcn_per_clone.insert(std::pair< int, vector<int> >(chr,mxpc));
    }
  }
}

void Clone::allocate_tcn(){
  for (int chr=0;chr<=maxChr; chr++){
    if (tcn[chr] == NULL) continue;
    for (int t=0;t<nTimes;t++){
      if (tcn[chr][t] != NULL) delete [] tcn[chr][t];
      if (log_tcn[chr][t] != NULL) delete [] log_tcn[chr][t];
      tcn[chr][t] = new double [nLevels];
      log_tcn[chr][t] = new double [nLevels];
    }
  }
}

void Clone::set_tcn(){
  for (int chr=0; chr<=maxChr; chr++){
    if( tcn[chr] == NULL ) continue;
    double ncn = double(normal_copy[chr]);
    for (int t=0; t<nTimes; t++){    
      for (int l=0; l<nLevels; l++){
	tcn[chr][t][l]  = ncn*(1.0-purity[t]) + clone_spectrum[t][l];
	log_tcn[chr][t][l] = log(tcn[chr][t][l] + 1.0e-10);
      }
    }
  }
}

void Clone::allocate_cn_usage(){
  if (cn_usage == NULL){
    cn_usage = new double ** [nTimes];
    for (int t=0;t<nTimes;t++){
      cn_usage[t] = new double * [maxtcn+1];
      for (int cn=0; cn<=maxtcn; cn++) cn_usage[t][cn] = NULL;
    }
  }
  for (int t=0;t<nTimes;t++){
    for (int cn=0; cn<=maxtcn; cn++){
      if (cn_usage[t][cn] != NULL) delete [] cn_usage[t][cn];
      cn_usage[t][cn] = new double [nLevels];
    }
  }
}

void Clone::set_cn_usage(){
  for (int t=0; t<nTimes; t++){
    for (int cn=0;cn<=maxtcn;cn++){
      for (int l=0;l<nLevels;l++) cn_usage[t][cn][l] = 0;
    }
    for (int l=0;l<nLevels;l++){
      for (int j=0;j<nClones;j++){
	cn_usage[t][copynumber[l][j]][l] += gsl_matrix_get(freqs,t,j);
      }
    }
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
    for (int f=0; f<nFreq; f++) ct += maxtcn+1;
    margin_map = gsl_matrix_calloc( ct, nLevels);
    ct=0;
    for (int f=0; f<nFreq; f++){
      for (int i=0; i<=maxtcn; i++){
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



void Clone::get_complexity(){
  double cnaC=0, bafC=0, snvC=0;
  double cnaN=0, bafN=0, snvN=0;
  std::map<int, vector<int> >::iterator it;
  for (it=maxtcn_per_clone.begin(); it != maxtcn_per_clone.end(); it++){
    int chr = it->first;
    double val=1.0;
    for (int j=0; j<nClones; j++){
      val *= (double) maxtcn_per_clone[chr][j] + 1;
    }
    //val -= 1.0;
    if (cnaEmit->is_set && cnaEmit->chrs.count(chr)==1){
      cnaC += val*double(cnaEmit->nSites[cnaEmit->idx_of[chr]]);
      cnaN += double(cnaEmit->nSites[cnaEmit->idx_of[chr]]);      
    }
    if (bafEmit->is_set && bafEmit->chrs.count(chr)==1){
      bafC += val*double(bafEmit->nSites[bafEmit->idx_of[chr]]);
      bafN += double(bafEmit->nSites[bafEmit->idx_of[chr]]);      
    }
    if (snvEmit->is_set && snvEmit->chrs.count(chr)==1){
      snvC += val*double(snvEmit->nSites[snvEmit->idx_of[chr]]);
      snvN += double(snvEmit->nSites[snvEmit->idx_of[chr]]);      
    }
  }
  complexity = (double) nTimes*( nClones + cnaEmit->is_set);
  double size=0;
  if (cnaEmit->is_set){
    complexity += cnaC / cnaN;
    size += cnaN*double(nTimes);
    if (bafEmit->is_set) size += bafN*double(nTimes);
    if (snvEmit->is_set) size += snvN*double(nTimes);
  }
  else if (snvEmit->is_set){
    complexity += snvC / snvN;
    size += snvN*double(nTimes);;
    if (learn_priors){
      std::set<int>::iterator iter;
      for (iter=all_maxtcn.begin(); iter!=all_maxtcn.end(); iter++){
	complexity += *iter;
      }
      complexity += 1.0;//for fpr
    }
  }
  complexity *= log(size); 
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

void Clone::set_mass(gsl_vector * m){
  if ((int) m->size != nTimes) abort();
  gsl_vector_memcpy(mass,m);
  for (int t=0;t<nTimes;t++) log_mass->data[t] = log(mass->data[t]);
}



//probability weight in chromosomes with majority normal copy number
void Clone::get_cna_marginals(){
  if (cnaEmit->is_set == 0){
    cout<<"ERROR-1 in Clone::get_cna_marginals()\n";
    exit(1);
  }
  if (marginals != NULL) gsl_matrix_free(marginals);
  int ct = 0;
  for (int f=0; f<nFreq; f++) ct += maxtcn+1;
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
      part = gsl_matrix_subrow( marginals, s, ct, maxtcn+1);
      norm += gsl_blas_dasum(&part.vector);
    }
    gsl_matrix_view pt = gsl_matrix_submatrix( marginals, 0, ct, cnaEmit->nSamples, maxtcn+1);
    if (norm <= 0.0) abort();
    gsl_matrix_scale( &pt.matrix, 1.0 / norm);
    ct += maxtcn+1;
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


void Clone::set_bafSymMap(){
  if (bafSymMap != NULL) abort();
  bafSymMap = new gsl_matrix * [maxtcn+1];
  for (int cn=0; cn<=maxtcn; cn++){
    bafSymMap[cn] = gsl_matrix_calloc(maxtcn+1,maxtcn+1);    
    for (int i=0; i<=cn; i++){
      for (int j=0; j<=cn; j++){
	double val 
	  = ( j==i && 2*j == cn ) 
	  ? 1.0 
	  : ( (j==i || j==cn-i) ? 0.5 : 0.0);
	gsl_matrix_set( bafSymMap[cn], i, j, val);
      }
    } 
  }
}
