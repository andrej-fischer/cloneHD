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
    for (int i=0; i <= myEmit->gridSize; i++){
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




//get total copynumber tracks from file...
int get_phi( const char * phi_fn, Emission * myEmit){
  myEmit->cnmax_seen.clear();
  ifstream ifs;
  string line;
  stringstream line_ss;
  ifs.open( phi_fn, ios::in);
  if (ifs.fail()){
    printf("ERROR file %s cannot be opened.\n", phi_fn);
    exit(1);
  }
  int chr = 0, old_chr = -1;
  int cn_locusi=0, cn_locusf=0, nloci=0, locus=-1;
  int evt=0;
  double tcn,x;
  int mxcn;
  double * tcn_buff = new double [myEmit->nTimes];
  int global_mx=0;
  int sample=-1,pevt=0;
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
      if (cols == 1 + 3 + myEmit->nTimes + 1) pevt=1;
      line_ss.clear();//reset
      line_ss.str(line);
    }
    line_ss >> chr >> cn_locusi;
    if (pevt){
      line_ss >> nloci >> cn_locusf;
    }
    else{
      cn_locusf = cn_locusi;
    }
    if (chr != old_chr && (chr > myEmit->maxchr || myEmit->idx_of[chr] < 0)){
      printf("ERROR in file %s chromosome %i not present in data.\n", phi_fn, chr);
      cout<<line<<endl;
      exit(1);
    }
    if (old_chr != -1 && chr != old_chr){//new chromosome
      sample = myEmit->idx_of[old_chr];
      if ( evt < myEmit->nEvents[sample]){//remaining from last chr
	int oevt = evt-1;
	for (int e=evt; e<myEmit->nEvents[sample]; e++){
	  for (int t=0; t<myEmit->nTimes; t++){
	    myEmit->phi[t][sample][e] = myEmit->phi[t][sample][oevt];
	  }
	  myEmit->cnmax[sample][e] = myEmit->cnmax[sample][oevt];
	}
      }
      evt = 0;
    }
    old_chr = chr;
    sample  = myEmit->idx_of[chr];
    if ( evt >= myEmit->nEvents[sample]) continue;//chromosome is complete!
    locus = (int) myEmit->loci[sample][myEmit->idx_of_event[sample][evt]];//current target locus
    if( cn_locusf <  locus) continue;
    // now we are above or at next event...
    for (int t=0; t<myEmit->nTimes; t++){
      if (line_ss.good() != true) abort();
      line_ss >> tcn;
      tcn_buff[t] = tcn;
    }
    line_ss >> mxcn;
    global_mx = max(global_mx,mxcn);
    myEmit->cnmax_seen.insert(mxcn);
    //fill
    while(locus <= cn_locusf){
      for (int t=0; t<myEmit->nTimes; t++) myEmit->phi[t][sample][evt] = tcn_buff[t];
      myEmit->cnmax[sample][evt] = mxcn;
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
	myEmit->phi[t][sample][e] = myEmit->phi[t][sample][oevt];
      }
      myEmit->cnmax[sample][e] = myEmit->cnmax[sample][oevt];
    }
  }
  //done
  ifs.close();
  return(global_mx);
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
  if (myClone->cnvEmit->is_set){//***CNV JUMPS***
    if(opts.cnv_jumps_fn != NULL){//1. either use external jump probability track
      get_track( opts.cnv_jumps_fn, distdummy, myClone->cnvEmit->pjump, vardummy, myClone->cnvEmit);
      for (int s=0; s< myClone->cnvEmit->nSamples; s++){
	myClone->cnvEmit->coarse_grain_jumps( s, opts.min_jump, 5);
      }
      myClone->cnvEmit->get_events_via_jumps();
    }
    else if (opts.cnv_jump >= 0.0){//2. or constant jump probability per base
      myClone->cnvEmit->set_pjump(opts.cnv_jump);
    }
    else{
      abort();
    }
    myClone->cnvEmit->allocate_phi();
    myClone->cnvEmit->allocate_cnmax();
    if (myClone->bafEmit->is_set){
      myClone->bafEmit->map_idx_to_Event(myClone->cnvEmit);
    }
    else{
      if (myClone->snpEmit->is_set) myClone->snpEmit->map_idx_to_Event(myClone->cnvEmit);
    }
  }
  if ( myClone->bafEmit->is_set ){//***BAF JUMPS***
    if (opts.baf_jumps_fn != NULL){//1. either external jump probability track
      get_track( opts.baf_jumps_fn, distdummy, myClone->bafEmit->pjump, vardummy, myClone->bafEmit);
      for (int s=0; s< myClone->bafEmit->nSamples; s++){// ignore improbable jump events
	myClone->bafEmit->coarse_grain_jumps( s, opts.min_jump, 5);
      }
      myClone->bafEmit->get_events_via_jumps();
      if ( myClone->cnvEmit->is_set && opts.cnv_jumps_fn != NULL ){//allow transit at CNV jumps
	myClone->bafEmit->add_break_points_via_jumps( myClone->cnvEmit, opts.min_jump);
      }
      myClone->bafEmit->get_events_via_jumps();
    }
    else if ( myClone->cnvEmit->is_set && opts.cnv_jumps_fn != NULL ){//2. map the CNV jumps to BAF
      myClone->bafEmit->map_jumps(myClone->cnvEmit);
      myClone->bafEmit->get_events_via_jumps();
    }
    else if (opts.baf_jump >= 0.0){//or 3. constant jump probability per base
      myClone->bafEmit->set_pjump(opts.baf_jump);
    }
    else{
      abort();
    }
    if (myClone->cnvEmit->is_set || opts.cn_fn != NULL){
      myClone->bafEmit->allocate_phi();
      myClone->bafEmit->allocate_cnmax();
    }
    if (myClone->snpEmit->is_set) myClone->snpEmit->map_idx_to_Event(myClone->bafEmit);
  }
  if ( myClone->snpEmit->is_set ){//***SNP JUMPS***
    if ( opts.snp_jumps_fn != NULL ){
      get_track( opts.snp_jumps_fn, distdummy, myClone->snpEmit->pjump, vardummy, myClone->snpEmit);
      for (int s=0; s< myClone->snpEmit->nSamples; s++){//ignore improbable jump events
	myClone->snpEmit->coarse_grain_jumps( s, opts.min_jump, 5);
      }
      myClone->snpEmit->get_events_via_jumps();
      if ( opts.cnv_jumps_fn != NULL ){
	myClone->snpEmit->add_break_points_via_jumps( myClone->cnvEmit, opts.min_jump);
      }
      myClone->snpEmit->get_events_via_jumps();
    }
    else if (opts.snp_jump >= 0.0){
      myClone->snpEmit->set_pjump(opts.snp_jump);
    }
    if (myClone->cnvEmit->is_set || opts.cn_fn != NULL){
      myClone->snpEmit->allocate_phi();
      myClone->snpEmit->allocate_cnmax();
    }
  }
}





//***BIAS FIELD***
void get_bias_field( Clone * myClone, cmdl_opts& opts){
  myClone->cnvEmit->allocate_bias();
  get_bias( opts.bias_fn, myClone->cnvEmit);
  //normalize the bias field...
  double bias_mean=0.0,norm=0.0;
  for (int s=0; s< myClone->cnvEmit->nSamples; s++){
    int ncn = myClone->normal_copy[ myClone->cnvEmit->chr[s]];
    for (int i=0; i< myClone->cnvEmit->nSites[s]; i++){
      bias_mean += myClone->cnvEmit->bias[s][i] / double(ncn);
    }
    norm += double(myClone->cnvEmit->nSites[s]);
  }
  bias_mean /= norm;
  for (int s=0; s< myClone->cnvEmit->nSamples; s++){
    int ncn = myClone->normal_copy[ myClone->cnvEmit->chr[s]];
    for (int i=0; i< myClone->cnvEmit->nSites[s]; i++){
      myClone->cnvEmit->bias[s][i] /= bias_mean * double(ncn);
      myClone->cnvEmit->log_bias[s][i] = log(myClone->cnvEmit->bias[s][i]);
      for (int t=1; t<myClone->cnvEmit->nTimes; t++){
	myClone->cnvEmit->bias[s][i]     = myClone->cnvEmit->bias[s][i];
	myClone->cnvEmit->log_bias[s][i] = myClone->cnvEmit->log_bias[s][i];
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
  if (myClone->cnvEmit->is_set && mass != NULL){        
    for (int i=0; i<nPts; i++){
      printf("\r%i", i+1);
      cout<<flush;
      nclones = gsl_matrix_submatrix( clones, i*nT, 0, nT, nC);
      nmass   = gsl_vector_subvector( mass,   i*nT, nT); 
      /*for (int t=0;t<nT;t++){
	printf("%.4f ", gsl_vector_get(&nmass.vector,t));
	for (int n=0; n<nC;n++){
	  printf("%.4f ", gsl_matrix_get(&nclones.matrix,t,n));
	}
	cout<<endl;
      }
      */
      myClone->set( &nclones.matrix );		
      myClone->set_mass( &nmass.vector );
      myClone->get_all_total_llh();
      fprintf( clonal_fp, "%.8e %.8e %.8e %.8e\n", 
	       myClone->cnv_total_llh , myClone->baf_total_llh , myClone->snp_total_llh,
	       myClone->cnv_total_llh + myClone->baf_total_llh + myClone->snp_total_llh);
    }
  }
  else if ( myClone->snpEmit->is_set){   
    fprintf(clonal_fp, "\n");  
    for (int i=0; i<nPts; i++){
      printf("\r%i", i+1);
      cout<<flush;
      nclones = gsl_matrix_submatrix( clones, i*nT, 0, nT, nC);
      myClone->set(&nclones.matrix);		
      myClone->get_snp_total_llh();
      fprintf( clonal_fp, "%.8e %.8e %.8e %.8e\n", 
	       0. , 0. , myClone->snp_total_llh,
	       myClone->snp_total_llh);
    }
  }
  cout<<"...done"<<endl;
  fclose(clonal_fp);
}




//top-level clone inference function (checks dimensions and chooses sub-programmes)...
int infer_clones( gsl_matrix * Clones, gsl_vector * Mass, Clone * myClone, cmdl_opts& opts){
  //***check input***
  if (Mass != NULL && myClone->cnvEmit->is_set){
    if (Clones != NULL && Clones->size1 != Mass->size){
      cout<<"ERROR-1a in infer_clones()\n";
      exit(1);
    }
    if (myClone->nTimes != (int) Mass->size){
      cout<<"ERROR-1b in infer_clones()\n";
      exit(1);
    }
  }
  if (Clones != NULL && (int) Clones->size1 != myClone->nTimes){
    cout<<"ERROR-1c in infer_clones()\n";
    exit(1);
  }
  if (Clones != NULL){
    opts.force = (int) Clones->size2;
    opts.nmax  = opts.force+1;
  }
  //output file pointer 
  char clonal_out[1024];
  sprintf( clonal_out, "%s.clonal.txt", opts.pre);
  FILE * clonal_fp = fopen(clonal_out,"w");
  fprintf(clonal_fp, "# n cnv-llh baf-llh snv-llh total-llh total-bic\n");
  //containers for the best estimates, given n...
  gsl_vector ** best_mass   = new gsl_vector * [opts.nmax+1];
  gsl_matrix ** best_clones = new gsl_matrix * [opts.nmax+1];
  gsl_matrix ** best_priors = new gsl_matrix * [opts.nmax+1];
  for (int n=0; n<=opts.nmax; n++){
    best_mass[n]   = NULL;
    best_clones[n] = NULL;
    best_priors[n] = NULL;
  }
  //exit(0);
  double bic=0, llh=0, max_llh=0,cnv_llh=0,baf_llh=0,snp_llh=0, max_bic=0; 
  int steps=0,btrial=0,bestn=0;
  int nT = myClone->nTimes;
  // *** NO CLONE SCENARIO (n==0) ***
  if (opts.force <= 0){
    myClone->nClones = 0;
    myClone->set(NULL);
    llh = 0.0;
    if (myClone->cnvEmit->is_set){// ***CNV***
      //with CNV data, the mass must be found even for n=0...
      best_mass[0] = gsl_vector_alloc(nT);
      if ( Mass == NULL){//need to get masses...
	cnv_llh = cnv_only_mass_noclones(  best_mass[0], myClone, 0, steps);
	myClone->set_mass(best_mass[0]);	
      }
      else if ( Mass != NULL){//masses are fixed...
	gsl_vector_memcpy( best_mass[0],  Mass);
	myClone->set_mass( best_mass[0]);
	cnv_llh = myClone->get_cnv_total_llh();
      }   
      //set trivial (normal) total copynumber...  
      for (int s=0; s<myClone->cnvEmit->nSamples; s++){
	myClone->get_phi(s);
	if (myClone->bafEmit->is_set) myClone->map_phi( myClone->cnvEmit, s, myClone->bafEmit);
	if (myClone->snpEmit->is_set) myClone->map_phi( myClone->cnvEmit, s, myClone->snpEmit);
      }
    }
    if (myClone->bafEmit->is_set){//*** BAF - NO CLONE ***
      baf_llh = myClone->get_baf_total_llh();
    }
    if (myClone->snpEmit->is_set){//*** SNP - NO CLONE ***
      if (myClone->bulk_mean != NULL){
	snp_bulk_update(myClone);
	myClone->set_bulk_to_post();
      }
      snp_llh = myClone->get_snp_total_llh();
      if (myClone->bulk_prior != NULL) myClone->set_bulk_to_prior();
    }
    llh = cnv_llh + baf_llh + snp_llh;   
    max_bic = 2.0*llh;
    if (myClone->cnvEmit->is_set){//no. parameters
      double complexity = double(myClone->nTimes)*log( double(myClone->total_loci) );
      max_bic -= complexity;
      printf("complexity = %f\n", complexity);
    }
    report_results( cnv_llh, baf_llh, snp_llh, steps, best_mass[0], NULL);
    printf("\nno clone best total llh = %+.6e, bic = %+.6e\n", llh, max_bic);
    cout<<endl;
    fprintf(clonal_fp, "0 %.10e %.10e %.10e %.10e %.10e\n", 
	    cnv_llh, baf_llh, snp_llh, llh, max_bic);
  }
  //exit(0);
  // *** ONE OR MORE CLONES (n>0) ***
  for (int n=1; n <= opts.nmax; n++){
    if (opts.force >= 0 && opts.force != n) continue;//if n is forced, skip!
    myClone->nClones = n;
    best_clones[n] = gsl_matrix_alloc( nT, n);//clones always have to be found
    gsl_matrix * clones = gsl_matrix_alloc( nT, n);
    gsl_vector * mass   = NULL;
    gsl_matrix * priors = NULL;
    if (myClone->cnvEmit->is_set){//mass only for CNV data...
      mass         = gsl_vector_calloc(nT);
      best_mass[n] = gsl_vector_calloc(nT);
    }
    if (myClone->snpEmit->is_set){//SNP cn-priors only for SNP data...
      priors         = gsl_matrix_calloc(myClone->maxcn+1,myClone->maxcn+1);
      best_priors[n] = gsl_matrix_calloc(myClone->maxcn+1,myClone->maxcn+1);
    }
    for (int trial=0; trial<opts.trials; trial++){//independent trials...
      printf("\nTrial %2i of %2i:...\n", trial+1, opts.trials);
      //get clones/mass and llh value...
      double cl=0.0,bl=0.0,sl=0.0;
      get_clones( clones, Clones, mass, Mass, priors, myClone, opts, cl, bl, sl);
      llh = cl + bl + sl;
      //test llh value...
      if (trial==0 || llh > max_llh){
	max_llh = llh;
	cnv_llh = cl;
	baf_llh = bl;
	snp_llh = sl;
	btrial  = trial;//trial with the best solution
	// keep results...	
	gsl_matrix_memcpy( best_clones[n], clones);
	if (mass != NULL)   gsl_vector_memcpy( best_mass[n],   mass);
	if (priors != NULL) gsl_matrix_memcpy( best_priors[n], priors);
      }
      cout<<endl;
      if (Clones != NULL && Mass != NULL) break;
    }
    // get BIC value...
    myClone->get_complexity();
    printf("complexity = %f\n", myClone->complexity);
    bic = 2.0*max_llh - myClone->complexity;
    printf("%i clone best total llh (trial %i) = %+.6e, bic = %+.6e\n\n", n, btrial+1, max_llh, bic);
    fprintf(clonal_fp, "%i %.10e %.10e %.10e %.10e %.10e\n", 
	    n, cnv_llh, baf_llh, snp_llh, max_llh, bic);
    // test BIC value...
    if (bic > max_bic || (opts.force > 0 && opts.force == n) ){
      bestn   = n;
      max_bic = bic;     
    }
    //local cleanup...
    gsl_matrix_free(clones);
    if (mass != NULL)   gsl_vector_free(mass);
    if (priors != NULL) gsl_matrix_free(priors);
  }
  //print all results to file
  for (int n=0; n<=opts.nmax; n++){
    if (best_clones[n] == NULL) continue;
    fprintf(clonal_fp, "# %i clones\n", n);
    for (int t=0; t<myClone->nTimes; t++){
      if (best_mass[n] != NULL){
	fprintf(clonal_fp, "%.3f ", gsl_vector_get( best_mass[n], t));
      }
      if (best_clones[n] != NULL){
	for (int f=0; f<n; f++){
	  fprintf(clonal_fp, "%.6f ", gsl_matrix_get( best_clones[n], t, f));
	}
      }
      fprintf(clonal_fp, "\n");
    }
    /*if (best_priors[n] != NULL){//print priors???
      for (int i=0; i< (int) best_priors[n]->size1; i++){
      for (int j=0; j< (int) best_priors[n]->size2; j++){
      fprintf(clonal_fp, "%.3e ", gsl_matrix_get( best_priors[n], i, j));
      }
      fprintf(clonal_fp, "\n");
      }
      }
    */
  }
  //insert the best solutions...
  myClone->nClones = bestn;
  myClone->set(best_clones[bestn]);
  if (best_mass[bestn] != NULL){
    myClone->set_mass(best_mass[bestn]);
  }
  if (best_priors[bestn] != NULL && bestn > 0){
    myClone->set_cn_prior_snp( best_priors[bestn] );
  }
  /*if (myClone->cnvEmit->is_set==0 && myClone->bafEmit->is_set){
    myClone->set_cn_prior_baf();
  }
  */
  //cleanup
  for (int n=0; n<=opts.nmax; n++){
    if (best_mass[n] != NULL)   gsl_vector_free(best_mass[n]);
    if (best_clones[n] != NULL) gsl_matrix_free(best_clones[n]);
    if (best_priors[n] != NULL) gsl_matrix_free(best_priors[n]);
  }
  delete [] best_mass;
  delete [] best_clones;
  delete [] best_priors;
  fclose(clonal_fp);
  //done
  return(bestn);
}

// get max-llh clones/mass for fixed n>0...
double get_clones( gsl_matrix *& clones, 
		   gsl_matrix *& Clones, 
		   gsl_vector *& mass, 
		   gsl_vector *& Mass, 
		   gsl_matrix *& priors, 
		   Clone * myClone,
		   cmdl_opts& opts,
		   double& cnv_llh,
		   double& baf_llh,
		   double& snp_llh
		   ){
  int nC = myClone->nClones;
  //default set
  gsl_matrix_set_all( clones, 1.0/double(nC));
  myClone->set( clones );
  double llh=0;
  if (myClone->cnvEmit->is_set){
    llh = get_clones_cnv( clones, Clones, mass, Mass, myClone, opts, cnv_llh, baf_llh, snp_llh);
  }
  else if (myClone->bafEmit->is_set){
    llh = get_clones_baf( clones, Clones, myClone, opts);
    baf_llh = llh;
  }
  else if (myClone->snpEmit->is_set){
    if (myClone->snpEmit->connect){
      llh = get_clones_snp_wcorr( clones, Clones, myClone, opts);
    }
    else{
      llh = get_clones_snp_ncorr( clones, Clones, priors, myClone, opts); 
    }
    snp_llh = llh;
  }
  return(llh);
}



double get_clones_cnv( gsl_matrix *& clones, 
		       gsl_matrix *& Clones, 
		       gsl_vector *& mass, 
		       gsl_vector *& Mass, 
		       Clone * myClone,
		       cmdl_opts& opts,
		       double& cnv_llh,
		       double& baf_llh,
		       double& snp_llh
		       ){
  int nT = myClone->nTimes;
  int nC = myClone->nClones;
  int nL = myClone->nLevels;
  double cl=0,bl=0,sl=0,llh=0;
  int steps;
  gsl_matrix * nclones = gsl_matrix_alloc(nT,nC);
  gsl_vector * mem = gsl_vector_alloc(nC);
  if ( myClone->snpEmit->is_set && myClone->bulk_mean != NULL) myClone->set_bulk_to_prior();
  if (Mass==NULL && Clones == NULL){//nothing is fixed
    //*** STEP 1: FIND FIRST ESTIMATES OF MASS AND CLONES ***    
    //random initial values
    for (int t=0; t<nT; t++){
      set_random_start_freq( mem, (myClone->min_purity)->data[t]);
      gsl_matrix_set_row( clones, t, mem);
      double p = (1.5 - double(rand()) / double(RAND_MAX));
      mass->data[t] = 0.5*double(myClone->nmean->data[t]) * p;
    }
    llh = cnv_clones_mass( clones, mass, myClone, 0, steps, cnv_llh, baf_llh, snp_llh);
    report_results( cnv_llh, baf_llh, snp_llh, steps, mass, clones);
    gsl_matrix * candidate_masses = gsl_matrix_calloc(nL,nT);
    gsl_vector * levels = gsl_vector_calloc(nL);
    gsl_vector * nmass  = gsl_vector_calloc(nT);
    int repeat = 1;
    int iter   = 0;
    double nllh;
    while(repeat==1){
      repeat=0;
      iter++;
      //*** STEP 2: FIND/REFINE CANDIDATE MASSES (GIVEN FIXED CLONES/MASS) ***	
      get_candidate_masses( clones, mass, myClone, candidate_masses, levels, opts.min_occ);
      //*** STEP 3: FOR EACH CANDIDATE MASS, FIND NEW CLONES ***
      for (int i=0; i<nL; i++){
	//get new estimates for clones with fixed mass
	gsl_vector_view cmass = gsl_matrix_row(candidate_masses,i);
	if ( gsl_vector_isnull(&cmass.vector)) continue;
	//report...
	int level = levels->data[i];
	printf("\rGauging masses via state ");
	for (int j=0; j<nC; j++) printf("%i", myClone->copynumber[level][j]);
	printf(" (%6.3f %%) to: ", myClone->cn2_post->data[level] * 100.0);
	for (int t=0; t<nT; t++){
	  printf("%.3e ", gsl_matrix_get(candidate_masses, i, t));
	}
	cout<<endl;
	//fix candidate mass...
	myClone->set_mass( &cmass.vector );
	//find new clones, given this mass...
	for (int t=0; t<nT; t++){//random initial clones
	  set_random_start_freq( mem, myClone->min_purity->data[t]);
	  gsl_matrix_set_row( nclones, t, mem);
	}
	steps=0;
	nllh = cnv_clones_fixed_mass( nclones, myClone, opts.restarts, steps, cl, bl, sl);
	if (nllh > llh + 1.0 ){
	  llh = nllh;
	  cnv_llh = cl;
	  baf_llh = bl;
	  snp_llh = sl;
	  gsl_vector_memcpy( mass, &cmass.vector);
	  gsl_matrix_memcpy( clones, nclones);
	  printf("\r");
	  report_results( cl,bl,sl, steps, mass, clones);
	  repeat=1; 
	}
	//
	if (repeat==1){//***STEP 4: REFINE MASSES, then GOTO 2 ***
	  gsl_matrix_memcpy( nclones, clones);
	  gsl_vector_memcpy( nmass, mass);
	  steps=0;
	  nllh = cnv_clones_mass( nclones, nmass, myClone, 0, steps, cl, bl, sl);
	  if ( nllh > llh ){
	    llh = nllh;
	    cnv_llh = cl;
	    baf_llh = bl;
	    snp_llh = sl;
	    gsl_vector_memcpy( mass, nmass);
	    gsl_matrix_memcpy( clones, nclones);
	    printf("\r");
	    report_results( cl, bl, sl, steps, mass, clones);
	  }
	  break;
	}
      }
    }
    gsl_matrix_free(candidate_masses);
    gsl_vector_free(levels);
    gsl_vector_free(nmass);
    //SNP bulk update, if required...
    if (opts.bulk_updates>0){
      for (int it=0; it<opts.bulk_updates; it++){
	snp_bulk_update(myClone);
	myClone->set_bulk_to_post();
	if (myClone->bulk_prior==NULL) break;
	if (myClone->snpEmit->coarse_grained) myClone->get_snpEmitLog();
	gsl_matrix_memcpy( nclones, clones);
	gsl_vector_memcpy( nmass, mass);
	nllh = cnv_clones_mass( nclones, nmass, myClone, 0, steps, cl, bl, sl);
	if (nllh > llh + 1.0){
	  llh = nllh;
	  cnv_llh = cl;
	  baf_llh = bl;
	  snp_llh = sl;
	  gsl_vector_memcpy( mass, nmass);
	  gsl_matrix_memcpy( clones, nclones);
	  report_results( cl,bl,sl, steps, mass, clones);
	}
	else{
	  break;
	}
      }
    }
  }
  else if (Mass == NULL && Clones != NULL){// clones are fixed, mass is not
    gsl_matrix_memcpy( clones, Clones);
    myClone->set(Clones);
    gsl_vector * nmass = gsl_vector_alloc(mass->size);
    for (int run=0; run<opts.restarts; run++){      
      for (int t=0; t<nT; t++){
	double p = (1.5 - double(rand()) / double(RAND_MAX));
	nmass->data[t] = 0.5*double(myClone->nmean->data[t]) * p;
      }
      steps=0;
      double cl,bl,sl;
      double nllh = cnv_mass_fixed_clones( nmass, myClone, 0, steps, cl,bl,sl);
      if (run==0 || nllh > llh + 1.0){	
	cnv_llh=cl;
	baf_llh=bl;
	snp_llh=sl;
	llh = nllh;
	gsl_vector_memcpy(mass,nmass);
	printf("%2i of %3i\n",run+1,opts.restarts);
	report_results( cnv_llh, baf_llh, snp_llh, steps, mass, clones);
      }
    }
    gsl_vector_free(nmass);
    myClone->set_mass(mass);
    //SNP bulk update, if required...
    if (opts.bulk_updates>0){
      snp_bulk_update(myClone);
      myClone->set_bulk_to_post();
      if (myClone->snpEmit->coarse_grained) myClone->get_snpEmitLog();
      //get CNV posterior and SNP llh...
      myClone->alpha_cnv = new gsl_matrix * [myClone->cnvEmit->nSamples];
      myClone->gamma_cnv = new gsl_matrix * [myClone->cnvEmit->nSamples];
      int s;
      snp_llh=0; 
#pragma omp parallel for schedule( dynamic, 1) default(shared)
      for ( s=0; s<myClone->snpEmit->nSamples; s++){  
	int cnv_sample = myClone->cnvEmit->idx_of[myClone->snpEmit->chr[s]];
	myClone->alpha_cnv[cnv_sample]=NULL;
	myClone->gamma_cnv[cnv_sample]=NULL;
	myClone->get_cnv_posterior(cnv_sample);
	myClone->map_phi( myClone->cnvEmit, cnv_sample, myClone->snpEmit);
	double l;
	myClone->do_snp_Fwd( s, l);
	#pragma omp critical
	{
	  snp_llh += l;
	}
	gsl_matrix_free(myClone->gamma_cnv[cnv_sample]);  
      }
      delete [] myClone->gamma_cnv;
      delete [] myClone->alpha_cnv;
      report_results( cnv_llh, baf_llh, snp_llh, steps, mass, clones);	
    }
  }
  else if (Mass != NULL && Clones == NULL){// mass is fixed, clones are not
    gsl_vector_memcpy( mass, Mass);
    myClone->set_mass(Mass);
    for (int t=0; t<nT; t++){	
      set_random_start_freq( mem, (myClone->min_purity)->data[t]);
      gsl_matrix_set_row( clones, t, mem);
    }
    steps=0;
    llh = cnv_clones_fixed_mass( clones, myClone, opts.restarts, steps, cnv_llh, baf_llh, snp_llh);
    report_results( cnv_llh, baf_llh, snp_llh, steps, mass, clones);
    //SNP bulk update, if required...
    if (opts.bulk_updates>0){
      for (int it=0; it<opts.bulk_updates; it++){
	snp_bulk_update(myClone);
	myClone->set_bulk_to_post();
	if (myClone->bulk_prior==NULL) break;
	if (myClone->snpEmit->coarse_grained) myClone->get_snpEmitLog();
	gsl_matrix_memcpy( nclones, clones);
	double nllh = cnv_clones_fixed_mass( nclones, myClone, 0, steps, cl, bl, sl);
	if (nllh > llh ){
	  llh = nllh;
	  cnv_llh = cl;
	  baf_llh = bl;
	  snp_llh = sl;
	  gsl_matrix_memcpy( clones, nclones);
	  report_results( cl,bl,sl, steps, mass, clones);
	}
	else{
	  break;
	}
      }
    }
  }
  else if (Mass != NULL && Clones != NULL){//all is fixed
    gsl_matrix_memcpy( clones, Clones);
    myClone->set(Clones);
    gsl_vector_memcpy( mass, Mass);
    myClone->set_mass(Mass);
    //llh = cnv_llh_all_fixed( myClone);
    llh = myClone->get_all_total_llh();
    cnv_llh = myClone->cnv_total_llh;
    baf_llh = myClone->baf_total_llh;
    snp_llh = myClone->snp_total_llh;
    report_results(  cnv_llh, baf_llh, snp_llh, 0, mass, clones);
    //SNP bulk update, if required...
    if (opts.bulk_updates>0){
      snp_bulk_update(myClone);
      myClone->set_bulk_to_post();
      if (myClone->snpEmit->coarse_grained) myClone->get_snpEmitLog();
      //get CNV posterior and SNP llh...
      myClone->alpha_cnv = new gsl_matrix * [myClone->cnvEmit->nSamples];
      myClone->gamma_cnv = new gsl_matrix * [myClone->cnvEmit->nSamples];
      int s;
      snp_llh = 0.0;
      myClone->save_snp_alpha = 0; 
#pragma omp parallel for schedule( dynamic, 1) default(shared)
      for ( s=0; s<myClone->snpEmit->nSamples; s++){  
	int cnv_sample = myClone->cnvEmit->idx_of[myClone->snpEmit->chr[s]];
	myClone->alpha_cnv[cnv_sample]=NULL;
	myClone->gamma_cnv[cnv_sample]=NULL;
	myClone->get_cnv_posterior(cnv_sample);
	myClone->map_phi( myClone->cnvEmit, cnv_sample, myClone->snpEmit);
	double l=0;
	myClone->do_snp_Fwd(s,l);
#pragma omp critical
	{
	  snp_llh += l;
	}
	gsl_matrix_free(myClone->gamma_cnv[cnv_sample]);  
      }
      delete [] myClone->gamma_cnv;
      delete [] myClone->alpha_cnv;
      report_results( cnv_llh, baf_llh, snp_llh, 0, mass, clones);	
    }
  }
  gsl_matrix_free(nclones);
  gsl_vector_free(mem);
  return(llh);
}

double get_clones_baf( gsl_matrix *& clones, 
		       gsl_matrix *& Clones, 
		       Clone * myClone,
		       cmdl_opts& opts
		       ){
  int nT = myClone->nTimes;
  int nC = myClone->nClones;	
  double llh;
  //myClone->set_cn_prior_baf();
  if (Clones == NULL){
    gsl_vector * mem = gsl_vector_alloc(nC);
    gsl_matrix * nclones = gsl_matrix_alloc(nT,nC);
    for (int t=0; t<nT; t++){
      set_random_start_freq( mem, myClone->min_purity->data[t]);
      gsl_matrix_set_row( nclones, t, mem);
    }
    int steps=0;
    llh = baf_clones( nclones, myClone, opts.restarts, steps);
    report_results( 0, llh, 0, steps, NULL, nclones);
    myClone->set(nclones);
    gsl_matrix_memcpy(clones,nclones);
    gsl_vector_free(mem);
    gsl_matrix_free(nclones);
  }
  else{
    gsl_matrix_memcpy( clones, Clones);
    myClone->set(Clones);
    llh = myClone->get_baf_total_llh();
  }
  return(llh);
}


//***SNP ONLY CLONE INFERENCE***
double get_clones_snp_ncorr( gsl_matrix *& clones, 
			     gsl_matrix *& Clones, 
			     gsl_matrix *& priors, 
			     Clone * myClone,
			     cmdl_opts& opts
			     ){
  int nT = myClone->nTimes;
  int nC = myClone->nClones;
  int steps;
  double llh=0;
  gsl_vector * mem = gsl_vector_alloc(nC);
  //STEP 1: learn clones with fixed priors...
  myClone->initialize_cn_prior_snp();// initialize SNP copynumber prior
  gsl_matrix_memcpy( priors, myClone->init_cn_prior_snp);
  if ( Clones == NULL ){
    for (int t=0; t<nT; t++){
      set_random_start_freq( mem, (myClone->min_purity)->data[t]);
      gsl_matrix_set_row( clones, t, mem);
    }
    steps=0;
    llh = snp_clones_fixed_priors( clones, myClone, opts.restarts, steps);
    report_results( 0, 0, llh, steps, NULL, clones);
    myClone->set(clones);
  }
  else{
    gsl_matrix_memcpy( clones, Clones);
    myClone->set(Clones);
  }     
  //STEP 2: for uncorrelated SNP data, learn cn-priors as well...
  if (opts.learn_priors){
    gsl_matrix * npriors = gsl_matrix_alloc(priors->size1,priors->size2);
    gsl_matrix_memcpy(npriors,priors);
    if (Clones == NULL){//STEP 3a: learn clones and priors jointly...
      gsl_matrix * nclones = gsl_matrix_alloc(nT,nC);
      gsl_matrix_memcpy(nclones,clones);
      steps=0;
      double nllh = snp_clones_priors( nclones, npriors, myClone, opts.restarts, steps);
      if (nllh>llh){
	gsl_matrix_memcpy(clones,nclones);
	gsl_matrix_memcpy(priors,npriors);
	llh=nllh;
	report_results( 0, 0, llh, steps, NULL, clones);
	printf("Found these copynumber priors per clone\n");
	for (int i=0; i< (int) priors->size1; i++){
	  for (int j=0; j< (int) priors->size2; j++){
	    printf("%.3f ", gsl_matrix_get( priors, i, j));
	  }
	  cout<<endl;
	}
      }
      myClone->set_cn_prior_snp(priors);
      myClone->set(clones);
      gsl_matrix_free(nclones);
    }
    else{//STEP 4b: learn priors with fixed clones...
      myClone->set(Clones);
      gsl_matrix_memcpy( clones, Clones);
      llh = snp_priors_fixed_clones( priors, myClone, 0, steps);
      report_results( 0, 0, llh, steps, NULL, clones);
      printf("Found these copynumber priors per clone\n");
      for (int i=0; i< (int) priors->size1; i++){
	for (int j=0; j< (int) priors->size2; j++){
	  printf("%.3f ", gsl_matrix_get( priors, i, j));
	}
	cout<<endl;
      }
      myClone->set_cn_prior_snp(priors);
    }
    gsl_matrix_free(npriors);
  }
  gsl_vector_free(mem);
  return(llh);
}


double get_clones_snp_wcorr( gsl_matrix *& clones, 
			     gsl_matrix *& Clones, 
			     Clone * myClone,
			     cmdl_opts& opts
			     ){
  int nT = myClone->nTimes;
  int nC = myClone->nClones;
  int steps;
  double llh=0;
  gsl_vector * mem = gsl_vector_alloc(nC);
  gsl_matrix * nclones = gsl_matrix_alloc(nT,nC);
  //STEP 1: re-set bulk to original prior...
  if ( myClone->bulk_mean != NULL) myClone->set_bulk_to_prior();
  //STEP 2: learn clones with prior bulk...
  if ( Clones == NULL ){
    for (int t=0; t<nT; t++){
      set_random_start_freq( mem, (myClone->min_purity)->data[t]);
      gsl_matrix_set_row( clones, t, mem);
    }
    steps=0;
    llh = snp_clones_fixed_priors( clones, myClone, 0, steps);
    report_results( 0., 0., llh, steps, NULL, clones);
    myClone->set(clones);
  }
  else{
    gsl_matrix_memcpy( clones, Clones);
    myClone->set(Clones);
    llh = myClone->get_snp_total_llh();
    report_results( 0., 0., llh, 0, NULL, clones);
  }     
  //STEP 3: iterative bulk update
  if( myClone->bulk_mean != NULL && opts.bulk_updates > 0){
    gsl_matrix_memcpy(nclones,clones);
    for (int it=0; it < opts.bulk_updates; it++){
      snp_bulk_update(myClone);
      myClone->set_bulk_to_post();
      if (myClone->snpEmit->coarse_grained) myClone->get_snpEmitLog();
      if (Clones != NULL){
	llh = myClone->get_snp_total_llh();
	report_results( 0, 0, llh, 0, NULL, clones);
	break;
      }
      else if (myClone->bulk_prior!=NULL){
	double nllh = snp_clones_fixed_priors( nclones, myClone, 0, steps);
	if (nllh > llh+1.0){
	  gsl_matrix_memcpy(clones,nclones);
	  llh = nllh;
	  report_results( 0, 0, llh, steps, NULL, clones);
	}
	else{
	  break;
	}
      }
    }
  }
  return(llh);
}

double cnv_only_mass_noclones( gsl_vector *& best_mass, Clone * myClone, int restarts, int& steps){	
  gsl_vector * mass  = gsl_vector_calloc(myClone->nTimes);
  gsl_vector * range = gsl_vector_calloc(myClone->nTimes);
  double llh, max_llh;
  gsl_vector ** simplex = NULL;
  gsl_vector * lower    = NULL;
  Q_par Qpar;
  Qpar.myClone  = myClone;
  Qpar.nSimplex = 0;
  Qpar.simplexD.clear();
  Qpar.clones_fixed = 1;//to NULL
  Qpar.mass_fixed   = 0;
  Qpar.prior_fixed  = 1;//does not apply
  Qpar.cnv=1;
  Qpar.baf=0;
  Qpar.snp=0;
  void * param = static_cast<void*>(&Qpar);
  steps = 0;
  for (int trial=0; trial<10; trial++){//try several times to get mass estimates
    for (int t=0; t<myClone->nTimes;t++){
      double p = (1.5 - double(rand()) / double(RAND_MAX));
      mass->data[t] = 0.5*double(myClone->nmean->data[t]) * p;
    }
    int s = 0;  
    llh = - find_optimum_wrestarts( 0, simplex, lower, mass, range, param, &Q, 1.0e-3, restarts, s, 0);
    steps += s;
    if ( trial==0 || llh > max_llh ){
      max_llh = llh;
      gsl_vector_memcpy( best_mass, mass);
    }
  }
  myClone->set_mass(best_mass);
  gsl_vector_free(mass);
  gsl_vector_free(range);
  return(max_llh);
}





//find maxLLH mass and clones for CNV data only...
double cnv_clones_fixed_mass( gsl_matrix*& clones, Clone * myClone, int restarts, int& steps,
			      double& cnv_llh, double& baf_llh, double& snp_llh){
  int nT = myClone->nTimes;
  int nC = myClone->nClones;
  Q_par Qpar;
  Qpar.myClone = myClone; 
  gsl_vector ** simplex = new gsl_vector * [nT];
  gsl_vector * lower    = gsl_vector_alloc(nT);
  gsl_vector * other    = NULL;
  gsl_vector * range    = NULL;
  Qpar.nSimplex = nT;
  Qpar.simplexD.clear();
  for (int t=0; t<nT; t++){
    simplex[t] = gsl_vector_alloc(nC);
    gsl_matrix_get_row( simplex[t], clones, t);
    lower->data[t] = (myClone->min_purity)->data[t];
    Qpar.simplexD.push_back(nC);
  }
  Qpar.cnv = 1;
  Qpar.baf = (myClone->bafEmit->is_set) ? 1 : 0;  
  Qpar.snp = (myClone->snpEmit->is_set) ? 1 : 0;
  Qpar.clones_fixed = 0;
  Qpar.mass_fixed   = 1;
  Qpar.prior_fixed  = 1;
  void * param = static_cast<void*>(&Qpar);
  // get mass and clones...
  steps=0;
  double llh = -find_optimum_wrestarts( nT, simplex, lower, other, range, param, &Q, 1.0e-3, restarts, steps, 0);
  for (int t=0; t<nT; t++) gsl_matrix_set_row( clones, t, simplex[t]);
  //cleanup...
  for (int t=0; t<nT; t++) gsl_vector_free(simplex[t]);
  delete [] simplex;
  gsl_vector_free(lower);
  myClone->set(clones);
  llh = myClone->get_all_total_llh();
  //myClone->get_cnv_total_llh();
  cnv_llh = myClone->cnv_total_llh;
  baf_llh = myClone->baf_total_llh;
  snp_llh = myClone->snp_total_llh;
  return(llh);
}


double cnv_mass_fixed_clones( gsl_vector*& mass, Clone * myClone, int restarts, int& steps, 
			      double& cnv_llh, double& baf_llh, double& snp_llh){
  int nT = myClone->nTimes;
  Q_par Qpar;
  Qpar.myClone = myClone; 
  gsl_vector ** simplex = NULL;
  gsl_vector * lower    = NULL;
  gsl_vector * other    = gsl_vector_alloc(nT);
  gsl_vector * range    = gsl_vector_calloc(nT);;
  Qpar.nSimplex = 0;
  Qpar.simplexD.clear();
  for (int t=0; t<nT; t++) other->data[t] = mass->data[t];
  Qpar.cnv = 1;
  Qpar.baf = (myClone->bafEmit->is_set) ? 1 : 0;  
  Qpar.snp = (myClone->snpEmit->is_set) ? 1 : 0;
  Qpar.clones_fixed = 1;
  Qpar.mass_fixed   = 0;
  Qpar.prior_fixed  = 1;//does not apply here
  void * param = static_cast<void*>(&Qpar);
  // get mass and clones...
  steps=0;
  double llh = -find_optimum_wrestarts( 0, simplex, lower, other, range, param, &Q, 1.0e-3, restarts, steps, 0);
  gsl_vector_memcpy( mass, other);
  //cleanup...
  gsl_vector_free(other);
  gsl_vector_free(range);
  myClone->set_mass(mass);
  llh = myClone->get_all_total_llh();
  cnv_llh = myClone->cnv_total_llh;
  baf_llh = myClone->baf_total_llh;
  snp_llh = myClone->snp_total_llh;
  return(llh);
}




//learn both using all data...
double cnv_clones_mass( gsl_matrix*& clones, gsl_vector*& mass, Clone * myClone, int restarts, int& steps,
			double& cnv_llh, double& baf_llh, double& snp_llh){
  int nT = myClone->nTimes;
  int nC = myClone->nClones;
  Q_par Qpar;
  Qpar.myClone = myClone; 
  gsl_vector ** simplex = new gsl_vector * [nT];
  gsl_vector * lower    = gsl_vector_alloc(nT);
  gsl_vector * other    = gsl_vector_alloc(nT);
  gsl_vector * range    = gsl_vector_calloc(nT);;
  Qpar.nSimplex = nT;
  Qpar.simplexD.clear();
  for (int t=0; t<nT; t++){
    simplex[t] = gsl_vector_alloc(nC);
    gsl_matrix_get_row(simplex[t], clones, t);
    lower->data[t] = (myClone->min_purity)->data[t];
    Qpar.simplexD.push_back(nC);
  }
  for (int t=0; t<nT; t++) other->data[t] = mass->data[t];
  Qpar.cnv = 1;
  Qpar.baf = (myClone->bafEmit->is_set) ? 1 : 0;  
  Qpar.snp = (myClone->snpEmit->is_set) ? 1 : 0;
  Qpar.clones_fixed = 0;
  Qpar.mass_fixed   = 0;
  Qpar.prior_fixed  = 1;//does not apply here
  void * param = static_cast<void*>(&Qpar);
  // get mass and clones...
  steps=0;
  double llh = -find_optimum_wrestarts( nT, simplex, lower, other, range, param, &Q, 1.0e-3, restarts, steps, 0);
  //copy results...
  for (int t=0; t<nT; t++) gsl_matrix_set_row( clones, t, simplex[t]);
  gsl_vector_memcpy( mass, other);
  //cleanup...
  for (int t=0; t<nT; t++) gsl_vector_free(simplex[t]);
  delete [] simplex;
  gsl_vector_free(other);
  gsl_vector_free(lower);
  gsl_vector_free(range);
  myClone->set(clones);
  myClone->set_mass(mass);
  llh = myClone->get_all_total_llh();
  cnv_llh = myClone->cnv_total_llh;
  baf_llh = myClone->baf_total_llh;
  snp_llh = myClone->snp_total_llh;
  return(llh);
}



// get candidate masses from the heights in the data, one of which is all-normal (cn=2)
void get_candidate_masses( gsl_matrix * clones, 
			   gsl_vector * mass, 
			   Clone * myClone, 
			   gsl_matrix*& candidate_masses,
			   gsl_vector*& levels,
			   double min_occ){
  int nT = myClone->nTimes;
  int nL = myClone->nLevels;
  gsl_matrix_set_zero(candidate_masses);
  gsl_vector_set_zero(levels);
  myClone->set(clones);
  myClone->set_mass(mass);
  myClone->get_mass_candidates();//this gets all the candidate masses (sorted by occupancy) 
  int ct=0;
  for (int i=0; i< nL; i++){//loop through the candidate masses...
    int level = myClone->levels_sorted[i];
    if (myClone->cn2_post->data[level] < min_occ || i>=15) break;//***THRESHOLD OCCUPANCY***
    //test whether new masses are sufficiently different from everything sofar...
    int found=0;
    for (int j=0; j<i; j++){
      double maxdiff=0.0,md=0,m1,m2;
      for (int t=0; t<nT; t++){
	m1 = gsl_matrix_get(myClone->mass_candidates, i, t);
	m2 = gsl_matrix_get(myClone->mass_candidates, j, t);
	md = fabs( m1 - m2) / m1;
	maxdiff = max(md,maxdiff);
      }
      if ( maxdiff < 0.01){//less than 1% relative difference
	found=1;
	break;
      };
    }
    if (found==1) continue;
    gsl_vector_view row = gsl_matrix_row( myClone->mass_candidates, i);
    gsl_matrix_set_row( candidate_masses, ct, &row.vector);
    levels->data[ct] = level;
    ct++;
  }
}

/*
//evaluate LLH if all is fixed for CNV data *WITH* BAF/SNP...
double cnv_llh_all_fixed( Clone * myClone ){
  //allocate alpha/gamma for cnv
  myClone->alpha_cnv = new gsl_matrix * [myClone->cnvEmit->nSamples];
  myClone->gamma_cnv = new gsl_matrix * [myClone->cnvEmit->nSamples];
  for ( int s=0; s< myClone->cnvEmit->nSamples; s++){
    myClone->alpha_cnv[s] = gsl_matrix_alloc( myClone->cnvEmit->nEvents[s], myClone->nLevels);
    myClone->gamma_cnv[s] = gsl_matrix_alloc( myClone->cnvEmit->nEvents[s], myClone->nLevels);;
  }
  int sample;
  //CNV contribution....
  myClone->save_cnv_alpha = 1;
  myClone->cnv_total_llh  = 0.0;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
  for ( sample=0; sample< myClone->cnvEmit->nSamples; sample++){//START PARALLEL FOR
    //1. get posterior for cnv
    double llh=0,ent=0;
    myClone->do_cnv_Fwd( sample, llh);
#pragma omp critical
    {
      myClone->cnv_total_llh += llh;    
    }
    myClone->do_cnv_Bwd( sample, ent);
    gsl_matrix_free(myClone->alpha_cnv[sample]);
    myClone->alpha_cnv[sample] = NULL;
    //2. get phi for cnv (total copynumber)
    myClone->get_phi(sample);
  }//END PARALLEL FOR
  //
  myClone->total_llh = myClone->cnv_total_llh;
  //
  //local cleanup...
  delete [] myClone->alpha_cnv;
  myClone->alpha_cnv = NULL;
  //BAF contributions...
  myClone->save_baf_alpha = 0;
  myClone->baf_total_llh  = 0.0;
  if ( myClone->bafEmit->is_set){
#pragma omp parallel for schedule( dynamic, 1) default(shared)
    for ( sample=0; sample< myClone->cnvEmit->nSamples; sample++){//START PARALLEL FOR
      double llh_baf = 0.0;
      int cnv_chr    = myClone->cnvEmit->chr[sample];
      int baf_sample = myClone->bafEmit->idx_of[cnv_chr];
      if (baf_sample >= 0){  
	myClone->map_phi( myClone->bafEmit, sample);	
	myClone->do_baf_Fwd( baf_sample, llh_baf);
#pragma omp critical
	{
	  myClone->baf_total_llh += llh_baf;
	}
      }
    }//END PARALLEL FOR
    myClone->total_llh += myClone->baf_total_llh;
  }
  //SNP contributions...
  myClone->snp_total_llh = 0.0;
  myClone->save_snp_alpha = 0;
  if( myClone->snpEmit->is_set){
#pragma omp parallel for schedule( dynamic, 1) default(shared)
    for ( sample=0; sample< myClone->cnvEmit->nSamples; sample++){//START PARALLEL FOR
      double llh_snp = 0.0;
      int cnv_chr = myClone->cnvEmit->chr[sample];
      int snp_sample = myClone->snpEmit->idx_of[cnv_chr];
      if( snp_sample >= 0){
	myClone->map_phi( myClone->snpEmit, sample);	
	myClone->do_snp_Fwd( snp_sample, llh_snp);
#pragma omp critical
	{
	  myClone->snp_total_llh += llh_snp;
	}
      }
    }//END PARALLEL FOR
    myClone->total_llh += myClone->snp_total_llh;
  }
  //cleanup CNV posterior...
  for ( sample=0; sample< myClone->cnvEmit->nSamples; sample++){
    gsl_matrix_free( myClone->gamma_cnv[sample]);
  }
  delete [] myClone->gamma_cnv;
  myClone->gamma_cnv = NULL;
  //
  return(myClone->total_llh);
}
*/

double baf_clones( gsl_matrix*& clones, Clone * myClone, int restarts, int& steps){
  int nT = myClone->nTimes;
  int nC = myClone->nClones;
  Q_par Qpar;
  Qpar.myClone = myClone; 
  gsl_vector ** simplex = new gsl_vector * [nT];
  gsl_vector * lower    = gsl_vector_alloc(nT);
  gsl_vector * other    = NULL;
  gsl_vector * range    = NULL;
  Qpar.nSimplex = nT;
  Qpar.simplexD.clear();
  for (int t=0; t<nT; t++){
    simplex[t] = gsl_vector_alloc(nC);
    gsl_matrix_get_row( simplex[t], clones, t);
    lower->data[t] = myClone->min_purity->data[t];
    Qpar.simplexD.push_back(nC);
  }
  Qpar.cnv = 0;
  Qpar.baf = 1;
  Qpar.snp = 0;
  Qpar.clones_fixed = 0;
  Qpar.mass_fixed   = 1;//does not apply here
  Qpar.prior_fixed  = 1;//does not apply here
  void * param = static_cast<void*>(&Qpar);
  // get mass and clones...
  steps=0;
  double llh = -find_optimum_wrestarts( nT, simplex, lower, other, range, param, &Q, 1.0e-3, restarts, steps,0);
  //copy results...
  for (int t=0; t<nT; t++) gsl_matrix_set_row( clones, t, simplex[t]);
  myClone->set(clones);
  llh = myClone->get_baf_total_llh();
  //cleanup...
  for (int t=0; t<nT; t++) gsl_vector_free(simplex[t]);
  delete [] simplex;
  return(llh);
}


double snp_clones_fixed_priors( gsl_matrix*& clones, Clone * myClone, int restarts, int& steps){
  int nT = myClone->nTimes;
  int nC = myClone->nClones;
  Q_par Qpar;
  Qpar.myClone = myClone; 
  gsl_vector ** simplex = new gsl_vector * [nT];
  gsl_vector * lower    = gsl_vector_alloc(nT);
  gsl_vector * other    = NULL;
  gsl_vector * range    = NULL;
  Qpar.nSimplex = nT;
  Qpar.simplexD.clear();
  for (int t=0; t<nT; t++){
    simplex[t] = gsl_vector_alloc(nC);
    gsl_matrix_get_row(simplex[t],clones,t);
    lower->data[t] = (myClone->min_purity)->data[t];
    Qpar.simplexD.push_back(nC);
  }
  Qpar.cnv = 0;
  Qpar.baf = 0;
  Qpar.snp = 1;
  Qpar.clones_fixed = 0;
  Qpar.mass_fixed   = 1;// does not apply here
  Qpar.prior_fixed  = 1;// priors are fixed!
  void * param = static_cast<void*>(&Qpar);
  // get mass and clones...
  steps=0;
  double llh = - find_optimum_wrestarts( nT, simplex, lower, other, range, 
					 param, &Q, 1.0e-3, restarts, steps, 0);
  for (int t=0; t<nT; t++) gsl_matrix_set_row( clones, t, simplex[t]);
  myClone->set(clones);
  llh = myClone->get_snp_total_llh();
  //cleanup...
  for (int t=0; t<nT; t++) gsl_vector_free(simplex[t]);
  delete [] simplex;
  gsl_vector_free(lower);
  return(llh);
}




double snp_priors_fixed_clones( gsl_matrix*& priors, Clone * myClone, int restarts, int& steps){
  Q_par Qpar;
  Qpar.myClone = myClone;
  int nP=0;
  for (int cn=0; cn<(int) priors->size1; cn++){
    gsl_vector_view row = gsl_matrix_row(priors,cn);
    if ( gsl_vector_max(&row.vector) <= 0.0) continue;
    nP++;
  }
  gsl_vector ** simplex = new gsl_vector * [nP];
  gsl_vector * lower    = gsl_vector_alloc(nP);
  gsl_vector * other    = NULL;
  gsl_vector * range    = NULL;
  Qpar.nSimplex = nP;
  Qpar.simplexD.clear();
  int ct=0;
  for (int cn=0; cn<(int) priors->size1; cn++){
    gsl_vector_view row = gsl_matrix_row(priors,cn);
    if ( gsl_vector_max(&row.vector) <= 0.0) continue;
    simplex[ct] = gsl_vector_alloc(cn+1);
    for (int i=0; i<=cn; i++) gsl_vector_set( simplex[ct], i, gsl_matrix_get(priors,cn,i));
    lower->data[ct] = (cn==0) ? 0.0 : 1.0;
    Qpar.simplexD.push_back(cn+1);
    ct++;
  }
  Qpar.cnv = 0;
  Qpar.baf = 0;
  Qpar.snp = 1;
  Qpar.clones_fixed = 1;
  Qpar.mass_fixed   = 1;// does not apply here
  Qpar.prior_fixed  = 0;
  void * param = static_cast<void*>(&Qpar);
  // get mass and clones...
  steps=0;
  double llh = - find_optimum_wrestarts( nP, simplex, lower, other, range,
					 param, &Q, 1.0e-3, restarts, steps, 0);
  for (int p=0; p<nP; p++){
    int cn = Qpar.simplexD[p] - 1;
    for (int j=0; j<=cn; j++) gsl_matrix_set( priors, cn, j, (simplex[p])->data[j]);
  }
  myClone->set_cn_prior_snp(priors);
  llh = myClone->get_snp_total_llh();
  //cleanup...
  for (int p=0; p<nP; p++) gsl_vector_free(simplex[p]);
  delete [] simplex;
  gsl_vector_free(lower);
  return(llh);
}



double snp_clones_priors( gsl_matrix*& clones, gsl_matrix*& priors, Clone * myClone, int restarts, int& steps){
  int nT = myClone->nTimes;
  int nC = myClone->nClones;
  Q_par Qpar;
  Qpar.myClone = myClone;
  int nP=0;
  for (int cn=0; cn<(int) priors->size1; cn++){
    gsl_vector_view row = gsl_matrix_row(priors,cn);
    if ( gsl_vector_max(&row.vector) <= 0.0) continue;
    nP++;
  }
  gsl_vector ** simplex = new gsl_vector * [nT+nP];
  gsl_vector * lower    = gsl_vector_alloc(nT+nP);
  gsl_vector * other    = NULL;
  gsl_vector * range    = NULL;
  Qpar.nSimplex = nT+nP;
  Qpar.simplexD.clear();
  int ct=0;
  for (int t=0; t<nT; t++){
    simplex[ct] = gsl_vector_alloc(nC);
    gsl_matrix_get_row(simplex[t],clones,t);
    lower->data[ct] = (myClone->min_purity)->data[t];
    Qpar.simplexD.push_back(nC);
    ct++;
  }
  for (int cn=0; cn<(int) priors->size1; cn++){
    gsl_vector_view row = gsl_matrix_row(priors,cn);
    if ( gsl_vector_max(&row.vector) <= 0.0) continue;
    simplex[ct] = gsl_vector_alloc(cn+1);
    for (int i=0; i<=cn; i++) gsl_vector_set( simplex[ct], i, gsl_matrix_get(priors,cn,i));
    lower->data[ct] = (cn==0) ? 0.0 : 1.0;
    Qpar.simplexD.push_back(cn+1);
    ct++;
  }
  Qpar.cnv = 0;
  Qpar.baf = 0;
  Qpar.snp = 1;
  Qpar.clones_fixed = 0;
  Qpar.mass_fixed   = 1;// does not apply here
  Qpar.prior_fixed  = 0;
  void * param = static_cast<void*>(&Qpar);
  // get mass and clones...
  steps=0;
  double llh = - find_optimum_wrestarts( nT+nP, simplex, lower, other, range, 
					 param, &Q, 1.0e-3, restarts, steps, 0);
  for (int t=0; t<nT; t++){
    gsl_matrix_set_row(clones,t,simplex[t]);
    ct++;
  }
  for (int p=0; p<nP; p++){
    int cn = Qpar.simplexD[nT+p] - 1;
    for (int j=0; j<=cn; j++) gsl_matrix_set( priors, cn, j, (simplex[nT+p])->data[j]);
  }
  myClone->set(clones);
  myClone->set_cn_prior_snp(priors);
  llh = myClone->get_snp_total_llh();
  //cleanup...
  for (int i=0; i<nT+nP; i++) gsl_vector_free(simplex[i]);
  delete [] simplex;
  gsl_vector_free(lower);
  return(llh);
}




double Q( const gsl_vector * x, void * p){
  Q_par * Qpar = static_cast<Q_par*> (p); 
  int nT = Qpar->myClone->nTimes;
  int nC = Qpar->myClone->nClones;
  gsl_matrix * clones   = NULL;
  gsl_matrix * prior    = NULL;
  gsl_vector ** simplex = NULL;
  gsl_vector * lower    = NULL;
  gsl_vector * other    = NULL;
  gsl_vector * range    = NULL; 
  //allocate
  if ( Qpar->nSimplex > 0){
    if ((int) Qpar->simplexD.size() != Qpar->nSimplex) abort();
    simplex = new gsl_vector * [Qpar->nSimplex];
    for (int i=0; i<Qpar->nSimplex; i++){
      simplex[i] = gsl_vector_alloc( Qpar->simplexD[i] );
    }
    lower = gsl_vector_alloc(Qpar->nSimplex);
  }
  if ( Qpar->mass_fixed == 0){
    other = gsl_vector_calloc( nT);
    range = gsl_vector_calloc( nT);
  }
  //set lower
  int ct=0;
  if (Qpar->clones_fixed == 0){
    for (int t=0; t<nT; t++){
      lower->data[t] = gsl_vector_get(Qpar->myClone->min_purity,t);
      ct++;
    }
  }
  if (Qpar->prior_fixed == 0){
    lower->data[ct] = 0.0;
    ct++;
    while (ct < Qpar->nSimplex){
      lower->data[ct] = 1.0;
      ct++;
    }
  }
  //unmap 
  int err = arg_unmap( x, Qpar->nSimplex, simplex, lower, other, range);
  double LLH = 0.0;
  if (err==0){//inside range?
    ct=0;
    if ( Qpar->clones_fixed == 0 ){//set clones
      clones = gsl_matrix_alloc(nT,nC);
      for (int t=0; t<nT; t++){
	gsl_matrix_set_row( clones, t, simplex[t]);
	ct++;
      }
      Qpar->myClone->set(clones);
    }
    if ( Qpar->prior_fixed == 0 ){//set prior
      prior = gsl_matrix_calloc( Qpar->myClone->maxcn + 1, Qpar->myClone->maxcn + 1);
      for (int i=ct; i<Qpar->nSimplex; i++){
	int cn = (int) (simplex[i])->size - 1;
	for (int j=0; j<=cn; j++) gsl_matrix_set( prior, cn, j, gsl_vector_get(simplex[i],j));
      }
      Qpar->myClone->set_cn_prior_snp(prior);
    }
    if (Qpar->mass_fixed == 0){//set mass
      Qpar->myClone->set_mass(other);
    }
    //get total LLH
    if ( Qpar->cnv == 1 && Qpar->baf == 0 && Qpar->snp == 0 ){
      //get only the total llh for CNV *without* BAF/SNP
      LLH = Qpar->myClone->get_cnv_total_llh();
    }
    else if (Qpar->cnv == 1){//get llh for CNV *with* BAF/SNP
      //LLH = cnv_llh_all_fixed(Qpar->myClone);
      LLH = Qpar->myClone->get_all_total_llh();
    }
    else if (Qpar->baf == 1){
      // get the total llh for BAF *only*
      LLH = Qpar->myClone->get_baf_total_llh();
    }
    else if (Qpar->snp == 1){
      // get the total llh for SNP *only*
      LLH = Qpar->myClone->get_snp_total_llh();
    }
    else{
      abort();
    }
  }
  // clean up
  if (clones != NULL) gsl_matrix_free(clones);
  if (prior != NULL)  gsl_matrix_free(prior);
  if (simplex != NULL){
    for (int t=0; t<Qpar->nSimplex; t++) gsl_vector_free(simplex[t]);
    delete [] simplex;
  }
  if (other != NULL) gsl_vector_free(other);
  if (range != NULL) gsl_vector_free(range);
  if (lower != NULL) gsl_vector_free(lower);
  return( err==1 ? 1.0e20 : - LLH );
}




void snp_bulk_update(Clone * myClone){
  myClone->alpha_snp = new gsl_matrix * [myClone->snpEmit->nSamples];
  myClone->gamma_snp = new gsl_matrix * [myClone->snpEmit->nSamples];
  if (myClone->cnvEmit->is_set){
    myClone->alpha_cnv = new gsl_matrix * [myClone->cnvEmit->nSamples];
    myClone->gamma_cnv = new gsl_matrix * [myClone->cnvEmit->nSamples];
  }
  int cnv_sample=0;
  int s;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
  for ( s=0; s<myClone->snpEmit->nSamples; s++){  
    if (myClone->cnvEmit->is_set){
      cnv_sample = myClone->cnvEmit->idx_of[myClone->snpEmit->chr[s]];
      myClone->alpha_cnv[cnv_sample]=NULL;
      myClone->gamma_cnv[cnv_sample]=NULL;
      myClone->get_cnv_posterior(cnv_sample);
      myClone->map_phi( myClone->cnvEmit, cnv_sample, myClone->snpEmit);
    }
    myClone->update_bulk(s);
    if (myClone->cnvEmit->is_set) gsl_matrix_free(myClone->gamma_cnv[cnv_sample]);
  }  
  myClone->set_bulk_to_post();
  //clean up...
  if (myClone->cnvEmit->is_set){
    delete [] myClone->gamma_cnv;
    delete [] myClone->alpha_cnv;
    myClone->gamma_cnv = NULL;
    myClone->alpha_cnv = NULL;
  }
  delete [] myClone->gamma_snp;
  delete [] myClone->alpha_snp;
  myClone->gamma_snp = NULL;
  myClone->alpha_snp = NULL;	
}

void set_random_start_freq(gsl_vector *& freq, double lower){
  if (freq->size == 1){
    double p = double(rand()) / double(RAND_MAX);
    p += lower * (1.0 - p);
    gsl_vector_set(freq,0,p);
  }
  else{
    list<double> parts;
    list<double> draw;
    for (int j=0; j<(int) freq->size-1; j++){
      parts.push_back(double(rand()) / double(RAND_MAX));
    }
    parts.push_back(1.0);
    parts.sort();
    double p = double(rand()) / double(RAND_MAX);
    p += pow(lower,2) * (1.0 - p);
    p  = sqrt(p);
    double last=0.0;
    for (int j=0; j<(int) freq->size; j++){
      //gsl_vector_set( freq, j, p * (parts.front()-last));
      draw.push_back(p * (parts.front()-last));
      last = parts.front();
      parts.pop_front();
    }
    parts.clear();
    draw.sort();
    draw.reverse();
    for (int j=0; j<(int) freq->size; j++){
      gsl_vector_set( freq, j, draw.front());
      //printf("%.1e ", draw.front());
      draw.pop_front();
    }
    // cout<<endl;
    draw.clear();
  }
}

void report_results(double cllh, double bllh, double sllh, int steps, gsl_vector * mass, gsl_matrix * freq){
  if (cllh != 0.0) printf("cnv-llh = %e ", cllh);
  if (bllh != 0.0) printf("baf-llh = %e ", bllh);
  if (sllh != 0.0) printf("snv-llh = %e ", sllh);
  printf("tot-llh = %e ", cllh+bllh+sllh);
  printf("(%i)\n", steps);
  if (mass != NULL || freq != NULL){
    int nT = (mass != NULL) ? (int) mass->size : (int) freq->size1;
    for (int t=0; t<nT; t++){
      if (mass != NULL) printf("%.3f ", gsl_vector_get(mass,t));
      if (freq != NULL){
	for (int i=0;i<(int)freq->size2;i++){
	  printf("%.6f ", gsl_matrix_get(freq,t,i));
	}
      }
      cout<<endl;
    }
  }
}
