//cloneHD-functions.cpp

#include "common-functions.h"
#include "cloneHD-functions.h"

#define PI 3.1415926
#define LOG2 0.693147

using namespace std;

//general function to read in mean, variance and distribution...
void get_track(const char * track_fn, 
	       gsl_matrix **& distribution,  
	       double **& mean,
	       double **& var,
	       Emission * myEmit
	       ){
  ifstream ifs;
  string line;
  stringstream line_ss;
  ifs.open( track_fn, ios::in);
  if (ifs.fail()){
    printf("ERROR in get_track(): file %s cannot be opened.\n", track_fn);
    exit(1);
  }
  if (mean == NULL){
    printf("ERROR-1 in get_track().\n");
    exit(1);
  }
  int chr = 0, old_chr = -1, l=0, locus;
  double mn, sd, jp;
  while( ifs.good() ){
    line.clear();
    getline( ifs, line);
    if (line.empty()) break;
    if (line[0] == '#') continue;
    line_ss.clear();
    line_ss.str(line);
    line_ss >> chr >> locus; 
    l = (chr == old_chr) ? l+1 : 0;
    if (chr != old_chr && (chr > myEmit->maxchr || myEmit->idx_of[chr] < 0) ){
      printf("ERROR 1a in get_track(): unexpected chromosome.\n");
      cout<<line<<endl;
      exit(1);
    }
    // exit(0);
    old_chr = chr;
    if( (int) myEmit->loci[ myEmit->idx_of[chr] ][l] != locus){
      printf("ERROR 2 in get_track()\n");
      cout<<line<<endl;
      printf( "%i, %i, %i vs %i\n", myEmit->idx_of[chr], l, myEmit->loci[ myEmit->idx_of[chr] ][l], locus);
      exit(1);
    }
    line_ss >> mn;
    mean[ myEmit->idx_of[chr] ][l] = mn;
    if (var != NULL) {
      line_ss >> sd;
      var[ myEmit->idx_of[chr] ][l]  = pow(sd,2);
    }
    if (distribution == NULL) continue;
    double p;
    line_ss >> jp;
    for (int i=0; i < (int) (distribution[ myEmit->idx_of[chr] ])->size2; i++){
      if (line_ss.good() != true){
	printf("ERROR 3 in get_track()\n");
	exit(1);
      }
      line_ss >> p;
      gsl_matrix_set( distribution[ myEmit->idx_of[chr] ], l, i, p);
    }
  }
  ifs.close();
}


// get the maximum total c.n. mask per chromosome
void get_maxcn_mask(const char * mask_fn, Clone * myClone, int maxcn_gw){
  myClone->maxcn_mask.clear();
  myClone->all_maxcn.clear();
  myClone->maxcn = maxcn_gw;
  //use explicit upper limit per chromosome, and mxcn_gw elsewhere
  if (mask_fn != NULL){
    ifstream ifs;
    string line;
    stringstream line_ss;
    ifs.open( mask_fn, ios::in);
    if (ifs.fail()){
      printf("ERROR: file %s cannot be opened.\n", mask_fn);
      exit(1);
    }
    int chr = 0, chrMxcn=0;   
    while( ifs.good() ){
      line.clear();
      getline( ifs, line);
      if (line.empty()) break;
      if (line[0] == '#') continue;
      line_ss.clear();
      line_ss.str(line);
      line_ss >> chr;
      if (chr<0 || chr >= 100) abort();
      if ( myClone->maxcn_mask.count(chr) == 0){
	vector<int> mx;
	while( line_ss.good() ){
	  line_ss >> chrMxcn;
	  mx.push_back(chrMxcn);
	  myClone->all_maxcn.insert(chrMxcn);
	  myClone->maxcn = max( myClone->maxcn, chrMxcn);
	}
	myClone->maxcn_mask.insert(std::pair<int, vector<int> >(chr,mx) );
      }
      else{
	printf("ERROR: in file %s: chr %i appears twice.\n", mask_fn, chr);
	exit(1);
      }
    }
    ifs.close();
  }
  //insert for all chromosomes not fixed the limit maxcn_gw
  if ( myClone->cnaEmit->is_set){
    for (int s=0; s< myClone->cnaEmit->nSamples; s++){
      int cnaChr = myClone->cnaEmit->chr[s];
      if ( myClone->maxcn_mask.count(cnaChr) == 0){
	vector<int> mx;
	mx.push_back(maxcn_gw);
	myClone->maxcn_mask.insert(std::pair<int, vector<int> >(cnaChr,mx));
	myClone->all_maxcn.insert(maxcn_gw);
      }
    }
  }
  else if ( myClone->snvEmit->is_set){
    for (int s=0; s< myClone->snvEmit->nSamples; s++){
      int snvChr = myClone->snvEmit->chr[s];
      if ( myClone->maxcn_mask.count(snvChr) == 0){
	vector<int> mx;
	mx.push_back(maxcn_gw);
	myClone->maxcn_mask.insert(std::pair<int, vector<int> >(snvChr,mx));
	myClone->all_maxcn.insert(maxcn_gw);
      }
    }
  }
  else{
    abort();
  }
}



//get total copynumber tracks from file...
void get_mean_tcn( const char * cn_fn, Clone * myClone, Emission * myEmit){
  ifstream ifs;
  string line;
  stringstream line_ss;
  ifs.open( mntcn_fn, ios::in);
  if (ifs.fail()){
    printf("ERROR file %s cannot be opened.\n", mntcn_fn);
    exit(1);
  }
  int chr=0, old_chr = -1;
  int cn_locusi=0, cn_locusf=0, nloci=0, locus=-1;
  int evt=0;
  double mcn, x;
  double * buff = new double [myEmit->nTimes];
  int sample=-1, perevt=0, persite=0;
  int wait=0;
  while( ifs.good() ){
    line.clear();
    getline( ifs, line);
    if (line.empty()) break;
    if (line[0] == '#') continue;
    line_ss.clear();
    line_ss.str(line);
    if (old_chr == -1 && locus == -1){//check format
      int cols=0;
      while(line_ss >> x) cols++;
      if (cols == 1 + 3 + myEmit->nTimes){
	perevt = 1;
      }
      else if (cols == 1 + 1 + myEmit->nTimes){
	persite = 1;
      }
      else{
	printf("ERROR: check format in %s.\n", mntcn_fn);
	exit(1);
      }
      line_ss.clear();//reset
      line_ss.str(line);
    }
    line_ss >> chr >> cn_locusi;
    if (perevt){
      line_ss >> nloci >> cn_locusf;
    }
    else{
      cn_locusf = cn_locusi;
      nloci = 1;
    }
    if (chr != old_chr){//new chromosome
      if(myEmit->chrs.count(chr) == 0){
	wait=1;
      }
      else{
	wait=0;
      }
      if (old_chr != -1 ){//not the first new chr
	sample = myEmit->idx_of[old_chr];
	if ( evt < myEmit->nEvents[sample]){//remaining entries from last chr
	  int oevt = evt-1;
	  for (int e=evt; e<myEmit->nEvents[sample]; e++){
	    for (int t=0; t<myEmit->nTimes; t++){
	      myEmit->mean_tcn[t][sample][e] = myEmit->mean_tcn[t][sample][oevt];
	    }
	  }
	}
	evt = 0;
      }
    }
    old_chr = chr;
    if (wait) continue;
    sample  = myEmit->idx_of[chr];
    if ( evt >= myEmit->nEvents[sample]) continue;//chromosome is complete!
    locus = (int) myEmit->loci[sample][myEmit->idx_of_event[sample][evt]];//current target locus
    if( cn_locusf <  locus) continue;
    // now we are above or at next event...
    for (int t=0; t<myEmit->nTimes; t++){
      if (line_ss.good() != true) abort();
      line_ss >> mcn;
      buff[t] = mcn;
    }
    //fill
    while(locus <= cn_locusf){
      for (int t=0; t<myEmit->nTimes; t++) myEmit->mean_tcn[t][sample][evt] = buff[t];
      evt++;
      if ( evt >= myEmit->nEvents[sample]) break;
      locus = (int) myEmit->loci[sample][myEmit->idx_of_event[sample][evt]];
    }
  }
  //last bit
  sample = myEmit->idx_of[old_chr];
  if (evt < myEmit->nEvents[sample]){
    int oevt = evt-1;
    for (int e=evt; e<myEmit->nEvents[sample]; e++){
      for (int t=0; t<myEmit->nTimes; t++){
	myEmit->mean_tcn[t][sample][e] = myEmit->mean_tcn[t][sample][oevt];
      }
    }
  }
  //done
  ifs.close();
}




//get copy number availability
void get_avail_cn( const char * avcn_fn, Clone * myClone, Emission * myEmit){
  ifstream ifs;
  string line;
  stringstream line_ss;
  ifs.open( avcn_fn, ios::in);
  if (ifs.fail()){
    printf("ERROR file %s cannot be opened.\n", avcn_fn);
    exit(1);
  }
  int chr=0, old_chr = -1;
  int cn_locusi=0, cn_locusf=0, nloci=0, locus=-1;
  int evt=0;
  double frac, x;
  int nEl = (myClone->maxcn+1)*myClone->nTimes;
  double ** buff = new double [myClone->nTimes];
  for (int t=0; t<myClone->nTimes; t++) buff[t] = new double [myClone->maxcn+1];
  int sample=-1, perevt=0, persite=0;
  int wait=0;
  while( ifs.good() ){
    line.clear();
    getline( ifs, line);
    if (line.empty()) break;
    if (line[0] == '#') continue;
    line_ss.clear();
    line_ss.str(line);
    if (old_chr == -1 && locus == -1){//check format
      int cols=0;
      while(line_ss >> x) cols++;
      if (cols == 1 + 3 + nEl){
	perevt = 1;
      }
      else if (cols == 1 + 1 + nEl){
	persite = 1;
      }
      else{
	printf("ERROR: check format in %s.\n", avcn_fn);
	exit(1);
      }
      line_ss.clear();//reset
      line_ss.str(line);
    }
    line_ss >> chr >> cn_locusi;
    if (perevt){
      line_ss >> nloci >> cn_locusf;
    }
    else{
      cn_locusf = cn_locusi;
      nloci = 1;
    }
    if (chr != old_chr){//new chromosome
      if(myEmit->chrs.count[chr] == 0){
	wait = 1;
      }
      else{
	wait = 0;
      }
      if (old_chr != -1 ){//not the first new chr
	sample = myEmit->idx_of[old_chr];
	if ( evt < myEmit->nEvents[sample]){//remaining from last chr
	  int oevt = evt-1;
	  for (int e=evt; e<myEmit->nEvents[sample]; e++){
	    for (int t=0; t<myEmit->nTimes; t++){
	      for (int cn=0; cn<=myClone->maxcn; cn++){
		myEmit->av_cn[t][sample][e][cn] = myEmit->av_cn[t][sample][oevt][cn];
	      }
	    }
	  }
	}
	evt = 0;
      }
    }
    old_chr = chr;
    if (wait) continue;
    sample  = myEmit->idx_of[chr];
    if ( evt >= myEmit->nEvents[sample]) continue;//chromosome is complete!
    locus = (int) myEmit->loci[sample][myEmit->idx_of_event[sample][evt]];//current target locus
    if( cn_locusf <  locus) continue;
    // now we are above or at next event...
    for (int t=0; t<myEmit->nTimes; t++){
      for (int cn=0; cn<=myClone->maxcn; cn++){
	if (line_ss.good() != true) abort();
	line_ss >> frac;
	buff[t][cn] = frac;
      }
    }
    // now fill
    while(locus <= cn_locusf){
      for (int t=0; t<myEmit->nTimes; t++){
	for (int cn=0; cn<=myClone->maxcn; cn++){
	  myEmit->av_cn[t][sample][evt][cn] = buff[t][cn];
	}
      }
      evt++;
      if ( evt >= myEmit->nEvents[sample]) break;
      locus = (int) myEmit->loci[sample][myEmit->idx_of_event[sample][evt]];
    }
  }
  //last bit
  sample = myEmit->idx_of[old_chr];
  if (evt < myEmit->nEvents[sample]){
    int oevt = evt-1;
    for (int e=evt; e<myEmit->nEvents[sample]; e++){
      for (int t=0; t<myEmit->nTimes; t++){
	for (int cn=0; cn<=myClone->maxcn; cn++){
	  myEmit->av_cn[t][sample][e][cn] = myEmit->av_cn[t][sample][oevt][cn];
	}
      }
    }
  }
  //done
  ifs.close();
}





//read in frequencies to act as lower limit in minimisations...
void get_purity( const char * purity_fn, gsl_vector *& purity){
  printf("Using data in %s for sample purities (lower bounds for lambdas)...\n", purity_fn);
  ifstream ifs;
  string line;
  stringstream line_ss;
  ifs.open( purity_fn, ios::in);
  if (ifs.fail()){
    printf("ERROR-1 in get_purity(): file %s cannot be opened.\n", purity_fn);
    exit(1);
  }
  double x;
  int ct=0;
  while( ifs.good() ){
    line.clear();
    getline( ifs, line);
    if (line.empty()) break;
    line_ss.clear();
    line_ss.str(line);
    line_ss >> x;
    if (ct >= (int) purity->size){
      cout<<"ERROR-2 in get_purity()\n"<<endl;
      exit(1);
    }
    else{
      purity->data[ct] = x;
    }
    ct++;
  }
  ifs.close();
  if (ct != (int) purity->size){
    cout<<"ERROR-3 in get_purity()\n"<<endl;
    exit(1);
  }
}


//read in clone frequencies and mass parameters from file...
void get_fixed_clones(gsl_matrix *& clones, gsl_vector *& mass, const char * clones_fn, int nT){
  printf("Using data in %s as input...\n", clones_fn);
  ifstream ifs;
  string line;
  stringstream line_ss;
  ifs.open( clones_fn, ios::in);
  if (ifs.fail()){
    printf("ERROR: file %s cannot be opened.\n", clones_fn);
    exit(1);
  }
  //
  string first;
  double x;
  int rows=0,cols=0,ocols=0;
  int with_mass = 0;
  while( ifs.good() ){//get dimensions...
    line.clear();
    getline( ifs, line);
    if (line.empty()) break;
    if (line[0] == '#') continue;
    line_ss.clear();
    line_ss.str(line);
    cols = 0;
    while(line_ss >> x){
      if (cols==0 && x>1.0) with_mass = 1;
      cols++;
    }
    if (ocols > 0 && cols != ocols){
      cout<<"ERROR in get_fixed_clones(): column sizes not consistent\n";
      exit(1);
    }
    else{
      ocols = cols;
    }
    rows++;
  }
  ifs.close();
  //allocate...
  if (rows < nT){
    printf("ERROR: not enough lines in file %s for the given data.\n", clones_fn);
    exit(1);
  }
  else if ( 0 != rows % nT){
    printf("ERROR: no. lines in file %s incompatible with the given data.\n", clones_fn);
    exit(1);
  }
  if ( cols - with_mass > 0 ) clones = gsl_matrix_alloc( rows, cols - with_mass);
  if (with_mass==1) mass = gsl_vector_alloc(rows);
  if (mass==NULL && clones==NULL){
    printf("ERROR: file %s seems corrupted (expect matrix).\n", clones_fn);
    exit(1);
  }
  ifs.open( clones_fn, ios::in);
  rows=0;
  while( ifs.good() ){
    line.clear();
    getline( ifs, line);
    if (line.empty()) break;
    if (line[0] == '#') continue;
    line_ss.clear();
    line_ss.str(line);
    cols = 0;
    while(line_ss >> x){
      if (with_mass == 1){
	if (cols == 0){
	  gsl_vector_set( mass, rows, x);
	}
	else{
	  gsl_matrix_set( clones, rows, cols-1, x);
	}
      }
      else{
	gsl_matrix_set( clones, rows, cols, x);
      }
      cols++;
    }
    rows++;
  }
  ifs.close();
  rows = (mass != NULL) ? (int) mass->size : (int) clones->size1;
  if (rows == nT){//print for single input...
    for (int i=0; i<rows; i++){
      if (mass != NULL) printf("%e ", mass->data[i]);
      if (clones!=NULL){
	for (int j=0; j<(int)clones->size2; j++){
	  printf("%e ", gsl_matrix_get(clones,i,j));
	}
      }
      cout<<endl;
    }
    cout<<endl;
  }
}



//***posterior jump probability track***
void get_jump_probability(  Clone * myClone, cmdl_opts& opts){
  double ** vardummy = NULL;
  gsl_matrix ** distdummy = NULL;
  if (myClone->cnaEmit->is_set){//***CNA JUMPS***
    if(opts.cna_jumps_fn != NULL){//1. either use external jump probability track
      get_track( opts.cna_jumps_fn, distdummy, myClone->cnaEmit->pjump, vardummy, myClone->cnaEmit);
      for (int s=0; s< myClone->cnaEmit->nSamples; s++){
	myClone->cnaEmit->coarse_grain_jumps( s, opts.min_jump, 5);
      }
      myClone->cnaEmit->get_events_via_jumps();
    }
    else if (opts.cna_jump >= 0.0){//2. or constant jump probability per base
      myClone->cnaEmit->set_pjump(opts.cna_jump);
    }
    else{
      abort();
    }
    myClone->cnaEmit->allocate_mean_tcn();
    if (myClone->bafEmit->is_set){//CNA+BAF: map CNA evts to BAF
      for (int s=0; s<myClone->bafEmit->nSamples; s++){
	myClone->bafEmit->map_idx_to_Event( myClone->cnaEmit, s);
      }
      if (myClone->snvEmit->is_set){//CNA+BAF+SNV: map CNA evts to SNV outside of autosomes
	for (int s=0; s<myClone->snvEmit->nSamples; s++){
	  int snvChr = myClone->snvEmit->chr[s];
	  if (myClone->bafEmit->chrs.count(snvChr) == 0){
	    myClone->snvEmit->map_idx_to_Event( myClone->cnaEmit, s);
	  }
	}
      }
    }
    else if (myClone->snvEmit->is_set){//CNA+SNV (no BAF): map CNA events to SNV
      for (int s=0; s<myClone->snvEmit->nSamples; s++){
	myClone->snvEmit->map_idx_to_Event( myClone->cnaEmit, s);
      }
    }
  }
  if ( myClone->bafEmit->is_set ){//***BAF JUMPS***
    if (opts.baf_jumps_fn != NULL){//1. either external jump probability track
      get_track( opts.baf_jumps_fn, distdummy, myClone->bafEmit->pjump, vardummy, myClone->bafEmit);
      for (int s=0; s< myClone->bafEmit->nSamples; s++){// ignore improbable jump events
	myClone->bafEmit->coarse_grain_jumps( s, opts.min_jump, 5);
      }
      myClone->bafEmit->get_events_via_jumps();
      if ( myClone->cnaEmit->is_set && opts.cna_jumps_fn != NULL ){//allow transit at CNA jumps
	myClone->bafEmit->add_break_points_via_jumps( myClone->cnaEmit, opts.min_jump);
      }
      myClone->bafEmit->get_events_via_jumps();
    }
    else if ( myClone->cnaEmit->is_set && opts.cna_jumps_fn != NULL ){//2. map the CNA jumps to BAF
      myClone->bafEmit->map_jumps(myClone->cnaEmit);
      myClone->bafEmit->get_events_via_jumps();
    }
    else if (opts.baf_jump >= 0.0){//or 3. constant jump probability per base
      myClone->bafEmit->set_pjump(opts.baf_jump);
    }
    else{
      abort();
    }
    myClone->bafEmit->allocate_mean_tcn();
    if (myClone->snvEmit->is_set){//CNA + BAF + SNV: map BAF to SNV for autosomes only
      for (int s=0; s<myClone->snvEmit->nSamples; s++){
	int snvChr = myClone->snvEmit->chr[s];
	if (myClone->bafEmit->chrs.count(snvChr) == 1){
	  myClone->snvEmit->map_idx_to_Event( myClone->bafEmit, s);
	}
      }
    }
  }
  if ( myClone->snvEmit->is_set ){//***SNV JUMPS***
    if ( opts.snv_jumps_fn != NULL ){
      get_track( opts.snv_jumps_fn, distdummy, myClone->snvEmit->pjump, vardummy, myClone->snvEmit);
      for (int s=0; s< myClone->snvEmit->nSamples; s++){//ignore improbable jump events
	myClone->snvEmit->coarse_grain_jumps( s, opts.min_jump, 5);
      }
      myClone->snvEmit->get_events_via_jumps();
      if ( opts.cna_jumps_fn != NULL ){//allow SNV to jump at CNA jump sites
	myClone->snvEmit->add_break_points_via_jumps( myClone->cnaEmit, opts.min_jump);
      }
      myClone->snvEmit->get_events_via_jumps();
    }
    else if (opts.snv_jump >= 0.0){
      myClone->snvEmit->set_pjump(opts.snv_jump);
    }    
    if (myClone->cnaEmit->is_set){
      myClone->snvEmit->allocate_mean_tcn();
    }
    else{
      if (opts.mntcn_fn != NULL) myClone->snvEmit->allocate_mean_tcn();
      if (opts.avcn_fn  != NULL) myClone->snvEmit->allocate_av_cn(myClone->maxcn);
    }
  }
}





//***BIAS FIELD***
void get_bias_field( Clone * myClone, cmdl_opts& opts){
  myClone->cnaEmit->allocate_bias();
  get_bias( opts.bias_fn, myClone->cnaEmit);
  //normalize the bias field...
  double bias_mean=0.0,norm=0.0;
  for (int s=0; s< myClone->cnaEmit->nSamples; s++){
    int cnaChr = myClone->cnaEmit->chr[s];
    int ncn    = myClone->normal_copy[cnaChr];
    for (int i=0; i< myClone->cnaEmit->nSites[s]; i++){
      bias_mean += myClone->cnaEmit->bias[s][i] / double(ncn);
    }
    norm += double(myClone->cnaEmit->nSites[s]);
  }
  bias_mean /= norm;
  for (int s=0; s< myClone->cnaEmit->nSamples; s++){
    int cnaChr = myClone->cnaEmit->chr[s];
    int ncn = myClone->normal_copy[cnaChr];
    for (int i=0; i< myClone->cnaEmit->nSites[s]; i++){
      myClone->cnaEmit->bias[s][i] /= bias_mean * double(ncn);
      myClone->cnaEmit->log_bias[s][i] = log(myClone->cnaEmit->bias[s][i]);
      for (int t=1; t<myClone->cnaEmit->nTimes; t++){
	myClone->cnaEmit->bias[s][i]     = myClone->cnaEmit->bias[s][i];
	myClone->cnaEmit->log_bias[s][i] = myClone->cnaEmit->log_bias[s][i];
      }
    }
  }
}



void print_llh_for_set(gsl_matrix * clones, gsl_vector * mass, Clone * myClone, cmdl_opts& opts){
  if (clones == NULL){
    cout<<"ERROR: all parameters have to be specified.\n";
    exit(1);
  }
  char clonal_out[1024];
  sprintf( clonal_out, "%s.llh-values.txt", opts.pre);
  FILE * clonal_fp = fopen(clonal_out,"w");
  fprintf(clonal_fp, "# cna-llh baf-llh snv-llh total-llh\n");
  int rows = (int) clones->size1;
  myClone->nClones = (int) clones->size2;
  int nC = myClone->nClones;
  int nT = myClone->nTimes;
  int nPts = rows / nT; 
  gsl_matrix_view nclones;
  gsl_vector_view nmass;
  printf("Printing log-likelihood values for %i parameter sets to %s\n", nPts, clonal_out);
  if (myClone->cnaEmit->is_set && mass != NULL){        
    for (int i=0; i<nPts; i++){
      printf("\r%i", i+1);
      cout<<flush;
      nclones = gsl_matrix_submatrix( clones, i*nT, 0, nT, nC);
      nmass   = gsl_vector_subvector( mass, i*nT, nT); 
      myClone->set( &nclones.matrix );		
      myClone->set_mass( &nmass.vector );
      myClone->get_all_total_llh();
      fprintf( clonal_fp, "%.8e %.8e %.8e %.8e\n", 
	       myClone->cna_total_llh , myClone->baf_total_llh , myClone->snv_total_llh,
	       myClone->cna_total_llh + myClone->baf_total_llh + myClone->snv_total_llh);
    }
  }
  else if ( myClone->snvEmit->is_set){   
    fprintf(clonal_fp, "\n");  
    for (int i=0; i<nPts; i++){
      printf("\r%i", i+1);
      cout<<flush;
      nclones = gsl_matrix_submatrix( clones, i*nT, 0, nT, nC);
      myClone->set(&nclones.matrix);		
      myClone->get_snv_total_llh();
      fprintf( clonal_fp, "%.8e %.8e %.8e %.8e\n", 
	       0. , 0. , myClone->snv_total_llh,
	       myClone->snv_total_llh);
    }
  }
  cout<<"...done"<<endl;
  fclose(clonal_fp);
}
