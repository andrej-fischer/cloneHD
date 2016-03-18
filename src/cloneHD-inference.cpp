//cloneHD-inference.cpp

#include "common-functions.h"
#include "cloneHD-functions.h"
#include "cloneHD-inference.h"
#include "clone.h"
#include "emission.h"
#include "minimization.h"

#define PI 3.1415926
#define LOG2 0.693147

using namespace std;


void report_results(double cllh, double bllh, double sllh, int steps, gsl_vector * mass, gsl_matrix * freq){
  if (cllh != 0.0) printf("cna-llh = %.8e ", cllh);
  if (bllh != 0.0) printf("baf-llh = %.8e ", bllh);
  if (sllh != 0.0) printf("snv-llh = %.8e ", sllh);
  printf("tot-llh = %.8e ", cllh+bllh+sllh);
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


//random points within simplex
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
      draw.push_back(p * (parts.front()-last));
      last = parts.front();
      parts.pop_front();
    }
    parts.clear();
    draw.sort();
    draw.reverse();
    for (int j=0; j<(int) freq->size; j++){
      gsl_vector_set( freq, j, draw.front());
      draw.pop_front();
    }
    draw.clear();
  }
}


//top-level clone inference function (checks dimensions and chooses sub-programmes)...
int infer_clones( gsl_matrix * Clones, gsl_vector * Mass, Clone * myClone, cmdl_opts& opts){
  Emission * cnaEmit = myClone->cnaEmit;
  Emission * bafEmit = myClone->bafEmit;
  Emission * snvEmit = myClone->snvEmit; 
  //***check input***
  if (Mass != NULL && cnaEmit->is_set){
    if (Clones != NULL && Clones->size1 != Mass->size) abort();
    if (myClone->nTimes != (int) Mass->size) abort();
  }
  if (Clones != NULL && (int) Clones->size1 != myClone->nTimes) abort();
  if (Clones != NULL){
    opts.force = (int) Clones->size2;
    opts.nmax  = opts.force+1;
  }
  //output file pointer 
  char clonal_out[1024];
  sprintf( clonal_out, "%s.summary.txt", opts.pre);
  FILE * clonal_fp = fopen(clonal_out,"w");
  if (clonal_fp == NULL){
    printf("ERROR: file %s could not be opened for writing.\n", clonal_out);
    exit(1);
  }
  fprintf( clonal_fp, "# n cna-llh baf-llh snv-llh total-llh total-bic\n");
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
  double bic=0, llh=0, max_llh=0, cna_llh=0, baf_llh=0, snv_llh=0, max_bic=0; 
  int steps=0,btrial=0,bestn=0;
  int nT = myClone->nTimes;
  // *** NO CLONE SCENARIO (n==0) ***
  if (opts.force <= 0){
    myClone->nClones = 0;
    myClone->set(NULL);
    llh = 0.0;
    if (cnaEmit->is_set){// ***CNA***
      //with CNA data, the mass must be found even for n=0...
      best_mass[0]   = gsl_vector_alloc(nT);
      best_clones[0] = gsl_matrix_calloc(nT,1);
      if ( Mass == NULL ){//need to get masses...
       cna_llh = cna_only_mass_noclones(  best_mass[0], myClone, 0, steps);
       myClone->set_mass(best_mass[0]);	
      }
      else if ( Mass != NULL ){//masses are fixed...
       gsl_vector_memcpy( best_mass[0],  Mass);
       myClone->set_mass( best_mass[0]);
       cna_llh = myClone->get_cna_total_llh();
      }   
      //set trivial (normal) mean total copynumber...  
      for (int cnaSample=0; cnaSample<cnaEmit->nSamples; cnaSample++){
       myClone->get_mean_tcn(cnaSample);
       int cnaChr = cnaEmit->chr[cnaSample];
       if (bafEmit->is_set && bafEmit->chrs.count(cnaChr) == 1){
         myClone->map_mean_tcn( cnaEmit, cnaSample, bafEmit);
         if (snvEmit->is_set && snvEmit->chrs.count(cnaChr) == 1){
           int bafSample = bafEmit->idx_of[cnaChr];
           myClone->map_mean_tcn( bafEmit, bafSample, snvEmit);
         }
       }
       else if (snvEmit->is_set && snvEmit->chrs.count(cnaChr) == 1){
         myClone->map_mean_tcn( cnaEmit, cnaSample, snvEmit);
       }
      }
    }
    if (bafEmit->is_set){//*** BAF - NO CLONE ***
      baf_llh = myClone->get_baf_total_llh();
    }
    if (snvEmit->is_set){//*** SNV - NO CLONE ***
      if (myClone->bulk_mean != NULL){
       snv_bulk_update(myClone);
       myClone->set_bulk_to_post();
      }
      snv_llh = myClone->get_snv_total_llh();
      if (myClone->bulk_mean != NULL)
      myClone->set_bulk_to_prior();
    }
    llh = cna_llh + baf_llh + snv_llh;//n==0 total LLH
    myClone->get_complexity();
    max_bic = 2.0*llh -  myClone->complexity;//n==0 BIC...
    printf("complexity = %f\n", myClone->complexity);
    /*if (cnaEmit->is_set){//no. parameters
      double complexity = double(myClone->nTimes)*log( double(myClone->total_loci) );
      max_bic -= complexity;   
      }*/
    report_results( cna_llh, baf_llh, snv_llh, steps, best_mass[0], NULL);
    printf("\nno clone best total llh = %+.6e, bic = %+.6e\n", llh, max_bic);
    cout<<endl;
    fprintf(clonal_fp, "0 %.10e %.10e %.10e %.10e %.10e\n", 
	    cna_llh, baf_llh, snv_llh, llh, max_bic);
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
    if (cnaEmit->is_set){//mass only for CNA data...
      mass         = gsl_vector_calloc(nT);
      best_mass[n] = gsl_vector_calloc(nT);
    }
    else if (snvEmit->is_set && myClone->learn_priors){//SNV cn-priors only for SNV data...
      priors         = gsl_matrix_calloc(myClone->maxtcn+1,myClone->maxtcn+1);
      best_priors[n] = gsl_matrix_calloc(myClone->maxtcn+1,myClone->maxtcn+1);
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
        cna_llh = cl;
        baf_llh = bl;
        snv_llh = sl;
	btrial  = trial;//trial with the best solution
	// keep results, sorted by sizw in t==1...	
        std::vector<int> idx;
        for (int i = 0; i < n; ++i) idx.push_back(i);
        SortDesc sd;
        gsl_vector_view row = gsl_matrix_row(clones,0);
        sd.arg = (&row.vector)->data;
        std::sort( idx.begin(), idx.end(), sd);
        for (int i = 0; i < nT; ++i){
          for (int j = 0; j < n; ++j){
            gsl_matrix_set(best_clones[n],i,j,gsl_matrix_get(clones,i,idx[j]));
          }
        }
        // gsl_matrix_memcpy( best_clones[n], clones);
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
    printf("%i clone best total llh (trial %i) = %+.8e, bic = %+.8e\n\n", n, btrial+1, max_llh, bic);
    fprintf(clonal_fp, "%i %.10e %.10e %.10e %.10e %.10e\n", n, cna_llh, baf_llh, snv_llh, max_llh, bic);
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
	  fprintf(clonal_fp, "%.4f ", gsl_matrix_get( best_clones[n], t, f));
	}
      }
      fprintf(clonal_fp, "\n");
    }
  }
  //insert the best solutions...
  myClone->nClones = bestn;
  myClone->set(best_clones[bestn]);
  if (best_mass[bestn] != NULL){
    myClone->set_mass(best_mass[bestn]);
  }
  if (!cnaEmit->is_set && snvEmit->is_set && !snvEmit->connect && bestn > 0){
    myClone->initialize_snv_prior_param();
  }
  if (best_priors[bestn] != NULL && bestn > 0){
    myClone->set_snv_prior( best_priors[bestn] );
  }
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
		   double& cna_llh,
		   double& baf_llh,
		   double& snv_llh
		   ){
  int nC = myClone->nClones;
  //default set
  gsl_matrix_set_all( clones, 1.0/double(nC));
  myClone->set( clones );
  double llh=0;
  if (myClone->cnaEmit->is_set){//CNA + (BAF + SNV)
    llh = get_clones_cna( clones, Clones, mass, Mass, myClone, opts, cna_llh, baf_llh, snv_llh);
  }
  else if (myClone->snvEmit->is_set){//SNV-only
    if (myClone->snvEmit->connect){
      llh = get_clones_snv_wcorr( clones, Clones, myClone, opts);
    }
    else{
      llh = get_clones_snv_ncorr( clones, Clones, priors, myClone, opts); 
    }
    snv_llh = llh;
  }
  return(llh);
}


double get_clones_cna( gsl_matrix *& clones, 
		       gsl_matrix *& Clones, 
		       gsl_vector *& mass, 
		       gsl_vector *& Mass, 
		       Clone * myClone,
		       cmdl_opts& opts,
		       double& cna_llh,
		       double& baf_llh,
		       double& snv_llh
		       ){
  Emission * cnaEmit = myClone->cnaEmit;
  //  Emission * bafEmit = myClone->bafEmit;
  Emission * snvEmit = myClone->snvEmit; 
  int nT = myClone->nTimes;
  int nC = myClone->nClones;
  int nL = myClone->nLevels;
  double cl=0,bl=0,sl=0,llh=0;
  int steps;
  gsl_matrix * nclones = gsl_matrix_alloc(nT,nC);
  gsl_vector * mem = gsl_vector_alloc(nC);
  if ( snvEmit->is_set && myClone->bulk_mean != NULL) myClone->set_bulk_to_prior();
  if (Mass==NULL && Clones == NULL){//nothing is fixed
    //*** STEP 1: FIND FIRST ESTIMATES OF MASS AND CLONES ***    
    //random initial values
    for (int t=0; t<nT; t++){
      set_random_start_freq( mem, (myClone->min_purity)->data[t]);
      gsl_matrix_set_row( clones, t, mem);
      double p = (1.5 - double(rand()) / double(RAND_MAX));
      mass->data[t] = 0.5*double(myClone->nmean->data[t]) * p;
    }
    int restarts = opts.mass_gauging ? 0 : opts.restarts;
    llh = cna_clones_mass( clones, mass, myClone, restarts, steps, cna_llh, baf_llh, snv_llh);
    report_results( cna_llh, baf_llh, snv_llh, steps, mass, clones);
    gsl_matrix * candidate_masses = gsl_matrix_calloc(nL,nT);
    gsl_vector * levels = gsl_vector_calloc(nL);
    gsl_vector * nmass  = gsl_vector_calloc(nT);
    int repeat = opts.mass_gauging;
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
	printf(" (%6.3f %%) to: ", myClone->majcn_post->data[level] * 100.0);
	for (int t=0; t<nT; t++){
	  printf("%.3e ", gsl_matrix_get(candidate_masses, i, t));
	}
	cout<<endl;
	//fix candidate masses...
	myClone->set_mass( &cmass.vector );
	//find new clones, given this mass...
	for (int t=0; t<nT; t++){//random initial clones
	  set_random_start_freq( mem, myClone->min_purity->data[t]);
	  gsl_matrix_set_row( nclones, t, mem);
	}
	steps=0;
	nllh = cna_clones_fixed_mass( nclones, myClone, opts.restarts, steps, cl, bl, sl);
	if (nllh > llh + 1.0 ){
	  llh = nllh;
	  cna_llh = cl;
	  baf_llh = bl;
	  snv_llh = sl;
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
	  nllh = cna_clones_mass( nclones, nmass, myClone, 0, steps, cl, bl, sl);
	  if ( nllh > llh ){
	    llh = nllh;
	    cna_llh = cl;
	    baf_llh = bl;
	    snv_llh = sl;
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
    //SNV bulk update, if required...
    if (opts.bulk_updates>0){
      for (int it=0; it<opts.bulk_updates; it++){
	snv_bulk_update(myClone);
	myClone->set_bulk_to_post();
	if (myClone->bulk_prior==NULL) break;
	if (snvEmit->coarse_grained) myClone->get_snvEmitLog();
	gsl_matrix_memcpy( nclones, clones);
	gsl_vector_memcpy( nmass, mass);
	nllh = cna_clones_mass( nclones, nmass, myClone, 0, steps, cl, bl, sl);
	if (nllh > llh + 1.0){
	  llh = nllh;
	  cna_llh = cl;
	  baf_llh = bl;
	  snv_llh = sl;
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
      double nllh = cna_mass_fixed_clones( nmass, myClone, 0, steps, cl,bl,sl);
      if (run==0 || nllh > llh + 1.0){	
	cna_llh=cl;
	baf_llh=bl;
	snv_llh=sl;
	llh = nllh;
	gsl_vector_memcpy(mass,nmass);
	printf("%2i of %3i\n",run+1,opts.restarts);
	report_results( cna_llh, baf_llh, snv_llh, steps, mass, clones);
      }
    }
    gsl_vector_free(nmass);
    myClone->set_mass(mass);
    //SNV bulk update (just one!), if required...
    if (opts.bulk_updates>0){
      snv_bulk_update(myClone);
      myClone->set_bulk_to_post();
      if (snvEmit->coarse_grained) myClone->get_snvEmitLog();
      //get CNA posterior and SNV llh...
      myClone->alpha_cna = new gsl_matrix * [cnaEmit->nSamples];
      myClone->gamma_cna = new gsl_matrix * [cnaEmit->nSamples];
      int s;
      snv_llh=0; 
      double*llhs;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic, 1) default(shared)
#endif
      for ( s=0; s<snvEmit->nSamples; s++){  
	int cna_sample = cnaEmit->idx_of[snvEmit->chr[s]];
	myClone->alpha_cna[cna_sample]=NULL;
	myClone->gamma_cna[cna_sample]=NULL;
	myClone->get_cna_posterior(cna_sample);
	myClone->get_mean_tcn( cna_sample);
	myClone->map_mean_tcn( cnaEmit, cna_sample, snvEmit);
	double l;
	myClone->do_snv_Fwd( s, l, llhs);
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  snv_llh += l;
	}
	gsl_matrix_free(myClone->gamma_cna[cna_sample]);  
      }
      delete [] myClone->gamma_cna;
      delete [] myClone->alpha_cna;
      report_results( cna_llh, baf_llh, snv_llh, steps, mass, clones);	
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
    llh = cna_clones_fixed_mass( clones, myClone, opts.restarts, steps, cna_llh, baf_llh, snv_llh);
    report_results( cna_llh, baf_llh, snv_llh, steps, mass, clones);
    //SNV bulk update, if required...
    if (opts.bulk_updates>0){
      for (int it=0; it<opts.bulk_updates; it++){
	snv_bulk_update(myClone);
	myClone->set_bulk_to_post();
	if (myClone->bulk_prior==NULL) break;
	if (myClone->snvEmit->coarse_grained) myClone->get_snvEmitLog();
	gsl_matrix_memcpy( nclones, clones);
	double nllh = cna_clones_fixed_mass( nclones, myClone, 0, steps, cl, bl, sl);
	if (nllh > llh ){
	  llh = nllh;
	  cna_llh = cl;
	  baf_llh = bl;
	  snv_llh = sl;
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
    //llh = cna_llh_all_fixed( myClone);
    llh = myClone->get_all_total_llh();
    cna_llh = myClone->cna_total_llh;
    baf_llh = myClone->baf_total_llh;
    snv_llh = myClone->snv_total_llh;
    report_results(  cna_llh, baf_llh, snv_llh, 0, mass, clones);
    //SNV bulk update, if required...
    if (opts.bulk_updates>0){
      snv_bulk_update(myClone);
      myClone->set_bulk_to_post();
      if (snvEmit->coarse_grained) myClone->get_snvEmitLog();
      //get CNA posterior and SNV llh...
      myClone->alpha_cna = new gsl_matrix * [cnaEmit->nSamples];
      myClone->gamma_cna = new gsl_matrix * [cnaEmit->nSamples];
      int s;
      snv_llh = 0.0;
      double * llhs = NULL;
      myClone->save_snv_alpha = 0; 
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic, 1) default(shared)
#endif
      for ( s=0; s<snvEmit->nSamples; s++){  
	int cna_sample = cnaEmit->idx_of[snvEmit->chr[s]];
	myClone->alpha_cna[cna_sample]=NULL;
	myClone->gamma_cna[cna_sample]=NULL;
	myClone->get_cna_posterior(cna_sample);
	myClone->get_mean_tcn(cna_sample);
	myClone->map_mean_tcn( cnaEmit, cna_sample, snvEmit);
	double l=0;
	myClone->do_snv_Fwd(s,l,llhs);
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  snv_llh += l;
	}
	gsl_matrix_free(myClone->gamma_cna[cna_sample]);  
      }
      delete [] myClone->gamma_cna;
      delete [] myClone->alpha_cna;
      report_results( cna_llh, baf_llh, snv_llh, 0, mass, clones);	
    }
  }
  gsl_matrix_free(nclones);
  gsl_vector_free(mem);
  return(llh);
}


//***SNV ONLY INFERENCE***
double get_clones_snv_ncorr( gsl_matrix *& clones, 
			     gsl_matrix *& Clones, 
			     gsl_matrix *& priors, 
			     Clone * myClone,
			     cmdl_opts& opts
			     ){//no correlations
  Emission * snvEmit = myClone->snvEmit;
  int nT = myClone->nTimes;
  int nC = myClone->nClones;
  int steps;
  double llh=0;
  gsl_vector * mem = gsl_vector_alloc(nC);
  myClone->initialize_snv_prior_param();
  //STEP 1: learn clones with fixed priors...
  if (myClone->learn_priors && snvEmit->av_cn == NULL){
    //myClone->initialize_snv_prior_param();
    printf("Using these SNV copynumber prior parameters\n");
    for (int i=0; i<=myClone->maxtcn; i++){
      for (int j=0; j<=myClone->maxtcn; j++){
	printf("%.3f ", gsl_matrix_get( myClone->initial_snv_prior_param, i, j));
      }
      cout<<endl;
    }
    if (priors==NULL) abort();
    gsl_matrix_memcpy( priors, myClone->initial_snv_prior_param);
  }
  if ( Clones == NULL ){
    for (int t=0; t<nT; t++){
      set_random_start_freq( mem, (myClone->min_purity)->data[t]);
      gsl_matrix_set_row( clones, t, mem);
    }
    steps=0;
    llh = snv_clones_fixed_priors( clones, myClone, opts.restarts, steps);
    report_results( 0, 0, llh, steps, NULL, clones);
    myClone->set(clones);
  }
  else{
    gsl_matrix_memcpy( clones, Clones);
    myClone->set(Clones);
    llh = myClone->get_snv_total_llh();
    report_results( 0, 0, llh, 0, NULL, Clones);
  }     
  //STEP 2: for uncorrelated SNV data, without availability constraint,
  // learn clones and SNV prior parameter...
  if (myClone->learn_priors){
    if (snvEmit->av_cn != NULL) abort();
    gsl_matrix * npriors = gsl_matrix_alloc( priors->size1, priors->size2);
    //irrelevant columns of 'priors' should already be set to zero
    gsl_matrix_memcpy( npriors, priors);
    if (Clones == NULL){//STEP 3a: learn clones and priors jointly...
      gsl_matrix * nclones = gsl_matrix_alloc(nT,nC);
      gsl_matrix_memcpy( nclones, clones);
      steps=0;
      double nllh = snv_clones_priors( nclones, npriors, myClone, opts.restarts, steps);
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
      myClone->set_snv_prior(priors);
      myClone->set(clones);
      gsl_matrix_free(nclones);
    }
    else{//STEP 4b: learn priors with fixed clones...
      myClone->set(Clones);
      gsl_matrix_memcpy( clones, Clones);
      llh = snv_priors_fixed_clones( priors, myClone, 0, steps);
      report_results( 0, 0, llh, steps, NULL, clones);
      printf("Found these copynumber priors per clone\n");
      for (int i=0; i< (int) priors->size1; i++){
	for (int j=0; j< (int) priors->size2; j++){
	  printf("%.3f ", gsl_matrix_get( priors, i, j));
	}
	cout<<endl;
      }
      myClone->set_snv_prior(priors);
    }
    gsl_matrix_free(npriors);
  }
  gsl_vector_free(mem);
  return(llh);
}



//SNV-only, w/ correlation
double get_clones_snv_wcorr( gsl_matrix *& clones, 
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
    llh = snv_clones_fixed_priors( clones, myClone, opts.restarts, steps);
    report_results( 0., 0., llh, steps, NULL, clones);
    myClone->set(clones);
  }
  else{
    gsl_matrix_memcpy( clones, Clones);
    myClone->set(Clones);
    llh = myClone->get_snv_total_llh();
    report_results( 0., 0., llh, 0, NULL, clones);
  }     
  //STEP 3: iterative bulk update
  if( myClone->bulk_mean != NULL && opts.bulk_updates > 0){
    gsl_matrix_memcpy(nclones,clones);
    for (int it=0; it < opts.bulk_updates; it++){
      snv_bulk_update(myClone);
      myClone->set_bulk_to_post();
      if (myClone->snvEmit->coarse_grained) myClone->get_snvEmitLog();
      if (Clones != NULL){
	llh = myClone->get_snv_total_llh();
	report_results( 0, 0, llh, 0, NULL, clones);
	break;
      }
      else if (myClone->bulk_prior!=NULL){
	double nllh = snv_clones_fixed_priors( nclones, myClone, 0, steps);
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

double cna_only_mass_noclones( gsl_vector *& best_mass, Clone * myClone, int restarts, int& steps){	
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
  Qpar.cna=1;
  Qpar.baf=0;
  Qpar.snv=0;
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





//find maxLLH mass and clones for CNA data only...
double cna_clones_fixed_mass( gsl_matrix*& clones, Clone * myClone, int restarts, int& steps,
			      double& cna_llh, double& baf_llh, double& snv_llh){
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
  Qpar.cna = 1;
  Qpar.baf = (myClone->bafEmit->is_set) ? 1 : 0;  
  Qpar.snv = (myClone->snvEmit->is_set) ? 1 : 0;
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
  //myClone->get_cna_total_llh();
  cna_llh = myClone->cna_total_llh;
  baf_llh = myClone->baf_total_llh;
  snv_llh = myClone->snv_total_llh;
  return(llh);
}


double cna_mass_fixed_clones( gsl_vector*& mass, Clone * myClone, int restarts, int& steps, 
			      double& cna_llh, double& baf_llh, double& snv_llh){
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
  Qpar.cna = 1;
  Qpar.baf = (myClone->bafEmit->is_set) ? 1 : 0;  
  Qpar.snv = (myClone->snvEmit->is_set) ? 1 : 0;
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
  cna_llh = myClone->cna_total_llh;
  baf_llh = myClone->baf_total_llh;
  snv_llh = myClone->snv_total_llh;
  return(llh);
}




//learn both using all data...
double cna_clones_mass( gsl_matrix*& clones, gsl_vector*& mass, Clone * myClone, int restarts, int& steps,
			double& cna_llh, double& baf_llh, double& snv_llh){
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
  Qpar.cna = 1;
  Qpar.baf = (myClone->bafEmit->is_set) ? 1 : 0;  
  Qpar.snv = (myClone->snvEmit->is_set) ? 1 : 0;
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
  cna_llh = myClone->cna_total_llh;
  baf_llh = myClone->baf_total_llh;
  snv_llh = myClone->snv_total_llh;
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
    //***THRESHOLD OCCUPANCY***
    if (myClone->majcn_post->data[level] < min_occ || i>=15 ) break;
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



double snv_clones_fixed_priors( gsl_matrix*& clones, Clone * myClone, int restarts, int& steps){
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
  Qpar.cna = 0;
  Qpar.baf = 0;
  Qpar.snv = 1;
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
  llh = myClone->get_snv_total_llh();
  //cleanup...
  for (int t=0; t<nT; t++) gsl_vector_free(simplex[t]);
  delete [] simplex;
  gsl_vector_free(lower);
  return(llh);
}


double snv_priors_fixed_clones( gsl_matrix*& priors, Clone * myClone, int restarts, int& steps){
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
  Qpar.cna = 0;
  Qpar.baf = 0;
  Qpar.snv = 1;
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
  myClone->set_snv_prior(priors);
  llh = myClone->get_snv_total_llh();
  //cleanup...
  for (int p=0; p<nP; p++) gsl_vector_free(simplex[p]);
  delete [] simplex;
  gsl_vector_free(lower);
  return(llh);
}


double snv_clones_priors( gsl_matrix*& clones, gsl_matrix*& priors, Clone * myClone, int restarts, int& steps){
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
  Qpar.cna = 0;
  Qpar.baf = 0;
  Qpar.snv = 1;
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
  myClone->set_snv_prior(priors);
  llh = myClone->get_snv_total_llh();
  //cleanup...
  for (int i=0; i<nT+nP; i++) gsl_vector_free(simplex[i]);
  delete [] simplex;
  gsl_vector_free(lower);
  return(llh);
}


double Q( const gsl_vector * x, void * p){
  Q_par * Qpar = static_cast<Q_par*> (p); 
  Clone * myClone = Qpar->myClone;
  int nT = myClone->nTimes;
  int nC = myClone->nClones;
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
      lower->data[t] = gsl_vector_get(myClone->min_purity,t);
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
      myClone->set(clones);
    }
    if ( Qpar->prior_fixed == 0 ){//set prior
      prior = gsl_matrix_calloc(myClone->maxtcn+1, myClone->maxtcn+1);
      for (int i=ct; i<Qpar->nSimplex; i++){
	int cn = (int) (simplex[i])->size - 1;
	for (int j=0; j<=cn; j++) gsl_matrix_set( prior, cn, j, gsl_vector_get(simplex[i],j));
      }
      myClone->set_snv_prior(prior);
    }
    if (Qpar->mass_fixed == 0){//set mass
      myClone->set_mass(other);
    }
    //get total LLH
    if ( Qpar->cna == 1 && Qpar->baf == 0 && Qpar->snv == 0 ){
      //get only the total llh for CNA *without* BAF/SNV
      LLH = myClone->get_cna_total_llh();
    }
    else if (Qpar->cna == 1){//get llh for CNA *with* BAF/SNV
      //LLH = cna_llh_all_fixed(Qpar->myClone);
      LLH = myClone->get_all_total_llh();
    }
    else if (Qpar->snv == 1){
      // get the total llh for SNV *only*
      LLH = myClone->get_snv_total_llh();
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




void snv_bulk_update(Clone * myClone){
  Emission * cnaEmit = myClone->cnaEmit;
  Emission * snvEmit = myClone->snvEmit;  
  myClone->alpha_snv = new gsl_matrix * [snvEmit->nSamples];
  myClone->gamma_snv = new gsl_matrix * [snvEmit->nSamples];
  if (cnaEmit->is_set){
    myClone->alpha_cna = new gsl_matrix * [cnaEmit->nSamples];
    myClone->gamma_cna = new gsl_matrix * [cnaEmit->nSamples];
  }
  int cna_sample=0;
  int s;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic, 1) default(shared)
#endif
  for ( s=0; s<snvEmit->nSamples; s++){  
    int snvChr = snvEmit->chr[s];
    if (cnaEmit->is_set && cnaEmit->chrs.count(snvChr) == 1){
      cna_sample = cnaEmit->idx_of[snvChr];
      myClone->alpha_cna[cna_sample] = NULL;
      myClone->gamma_cna[cna_sample] = NULL;
      myClone->get_cna_posterior(cna_sample);
      myClone->get_mean_tcn( cna_sample);
      myClone->map_mean_tcn( cnaEmit, cna_sample, snvEmit);
    }
    myClone->update_bulk(s);
    if (myClone->gamma_cna != NULL && myClone->gamma_cna[cna_sample] != NULL){
      gsl_matrix_free(myClone->gamma_cna[cna_sample]);
    }
  }  
  myClone->set_bulk_to_post();
  //clean up...
  if (cnaEmit->is_set){
    delete [] myClone->gamma_cna;
    delete [] myClone->alpha_cna;
    myClone->gamma_cna = NULL;
    myClone->alpha_cna = NULL;
  }
  delete [] myClone->gamma_snv;
  delete [] myClone->alpha_snv;
  myClone->gamma_snv = NULL;
  myClone->alpha_snv = NULL;	
}

