/*
  ******************************************************************************

Copyright (c) 11/12/13  Genome Research Ltd.

Author: Andrej Fischer (af7[at]sanger.ac.uk)

This file is part of cloneHD.

cloneHD is free software: you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation; either version 3 
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  
If not, see <http://www.gnu.org/licenses/>.


******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <vector>
#include <list>


// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"

//own headers...
#include "cloneHD-functions.h"
#include "common-functions.h"


using namespace std;


// *** MAIN START***
int main (int argc, const char * argv[]){
  cmdl_opts opts;
  get_opts( argc, argv, opts);
  vector<int> chrs;
  vector<int> idx_of;
  vector<int> nSites;
  int nTimes=0, nT=0;
  //*** EMITTED DATA OBJECTS ***
  Emission cnaEmit, bafEmit, snvEmit;
  if (opts.cna_fn != NULL){//CNA DATA
    get_dims( opts.cna_fn, nTimes, chrs, nSites);
    cnaEmit.mode = (opts.cna_shape > 0.0) ? 4 : 3;
    cnaEmit.shape = opts.cna_shape;
    cnaEmit.log_shape = (opts.cna_shape > 0.0) ? log(opts.cna_shape) : 0.0;
    cnaEmit.rnd_emit  = opts.cna_rnd;
    if (opts.cna_jump >= 0.0)      cnaEmit.connect = 1;
    if (opts.cna_jumps_fn != NULL) cnaEmit.connect = 1;
    cnaEmit.set( nTimes, chrs, nSites, opts.grid);
    get_data( opts.cna_fn, &cnaEmit);
  }
  if (opts.baf_fn != NULL){//BAF DATA
    get_dims( opts.baf_fn, nT, chrs, nSites);
    if (nTimes > 0 && nT != nTimes){
      cout<<"ERROR-1a in cloneHD main(): incompatible sample sizes in CNA and BAF data.\n";
      exit(1);
    }
    nTimes = nT;
    bafEmit.mode = (opts.baf_shape > 0.0) ? 2 : 1;//shape?
    bafEmit.shape = opts.baf_shape;
    bafEmit.log_shape = (opts.baf_shape > 0.0) ? log(opts.baf_shape) : 0.0;
    bafEmit.rnd_emit  = opts.baf_rnd;
    bafEmit.reflect = 1;//reflect at 0.5 (n == N-n)
    bafEmit.get_log = 1;
    if (opts.baf_jump >= 0.0)    bafEmit.connect = 1;
    if (opts.baf_jumps_fn!=NULL) bafEmit.connect = 1;
    if (opts.cna_jumps_fn!=NULL) bafEmit.connect = 1;
    bafEmit.set( nTimes, chrs, nSites, opts.grid);
    get_data( opts.baf_fn, &bafEmit);
    bafEmit.set_EmitProb(-1);
  }
  if (opts.snv_fn != NULL){//SNV DATA
    get_dims( opts.snv_fn, nT, chrs, nSites);
    if (nTimes>0 && nT != nTimes){
      cout<<"ERROR-1b in cloneHD:main(): incompatible sample sizes CNA and SNV data.\n";
      exit(1);
    }
    nTimes = nT;
    snvEmit.mode = (opts.snv_shape > 0.0) ? 2 : 1;
    snvEmit.shape = opts.snv_shape;
    snvEmit.log_shape = (opts.snv_shape > 0.0) ? log(opts.snv_shape) : 0.0;
    snvEmit.rnd_emit = opts.snv_rnd;
    snvEmit.get_log = 1;
    snvEmit.reflect = 0;
    if (opts.snv_jump >= 0.0)    snvEmit.connect = 1;
    if (opts.snv_jumps_fn!=NULL) snvEmit.connect = 1;
    snvEmit.set( nTimes, chrs, nSites, opts.grid);
    get_data( opts.snv_fn, &snvEmit);
    snvEmit.set_EmitProb(-1);
  }
  //*** ANNOUNCE ***
  printf("\ncloneHD: probabilistic inference of sub-clonality using...\n\n");
  if (cnaEmit.is_set){
    printf("CNA data in %s: %i sites in %i chr across %i samples\n", 
	   opts.cna_fn, cnaEmit.total_loci, cnaEmit.nSamples, nTimes);
  }
  if (bafEmit.is_set){
    printf("BAF data in %s: %i sites in %i chr across %i samples\n", 
	   opts.baf_fn, bafEmit.total_loci, bafEmit.nSamples, nTimes);
  }
  if (snvEmit.is_set){
    printf("SNV data in %s: %i sites in %i chr across %i samples\n", 
	   opts.snv_fn, snvEmit.total_loci, snvEmit.nSamples, nTimes);
  }
  cout<<endl;
  // *** ALLOCATE CLONE ***
  Clone myClone;
  myClone.allocate( &cnaEmit, &bafEmit, &snvEmit, opts.chr_fn);
  myClone.baf_pen  = (opts.baf_pen > 0.0)  ? opts.baf_pen : 1.0;// BAF penalty for complex chr status
  myClone.snv_fpr  = (opts.snv_fpr > 0.0)  ? opts.snv_fpr : 1.0e-4;//SNV false-positive rate
  myClone.snv_err  = (opts.snv_err >= 0.0) ? opts.snv_err : 0.0;   //SNV frequency of false positives
  myClone.snv_pen  = (opts.snv_pen > 0.0)  ? opts.snv_pen : (bafEmit.is_set ? 0.01 : 0.5);//penalty for high SNV genotypes
  myClone.bulk_fix = opts.bulk_fix;
  // *** GET MAXCN MASK ***
  get_maxcn_mask( opts.maxcn_mask_fn, &myClone, opts.maxcn);
  if (myClone.maxcn > opts.maxcn)
    printf("Found new higher maximum copy number %i in %s.\n", myClone.maxcn, opts.maxcn_mask_fn);
  // *** GET SNV BULK PRIOR ***
  if ( snvEmit.is_set && opts.bulk_fn != NULL ){
    if (opts.bulk_mean)  myClone.allocate_bulk_mean();
    if (opts.bulk_prior) myClone.allocate_bulk_dist();
    printf("Using data in %s as SNV bulk prior...\n", opts.bulk_fn);
    double ** vardummy  = NULL;
    gsl_matrix ** distdummy = NULL;
    if (opts.bulk_mean){
      get_track( opts.bulk_fn, distdummy, myClone.bulk_prior_mean, vardummy, &snvEmit);
    }
    else if (opts.bulk_prior){
      get_track( opts.bulk_fn, myClone.bulk_prior, myClone.bulk_prior_mean, vardummy, &snvEmit);
    }
    else{
      abort();
    }
  }
  //*** GET JUMP PROBABILITY TRACKS and COLLAPSE TO EVENTS***
  get_jump_probability( &myClone, opts);
  if ( snvEmit.is_set){
    // *** GET TOTAL MEAN COPYNUMBER TRACKS ***  
    if( opts.mntcn_fn != NULL ){
      get_mean_tcn( opts.mntcn_fn, &myClone, &snvEmit);
    } 
    // *** GET AVAILABLE COPYNUMBER TRACKS ***  
    if ( snvEmit.is_set && opts.avcn_fn != NULL ){
      get_avail_cn( opts.avcn_fn, &myClone, &snvEmit);
    }
    //check whether this is still needed
    /*if ( cnaEmit.is_set==0){
      snvEmit.cnmax_seen.clear();
      for (int s=0; s<snvEmit.nSamples; s++){//all chromosomes are normal???
      snvEmit.cnmax_seen.insert(myClone.normal_copy[snvEmit.chr[s]]);
      }
      }
    */
  }
  //*** GET READ DEPTH BIAS FIELD ***
  if (cnaEmit.is_set && opts.bias_fn != NULL){
    get_bias_field( &myClone, opts);
  }
  //*** PREPARE COARSE-GRAINED DATA ***
  if (cnaEmit.is_set && opts.cna_jumps_fn != NULL){
    cnaEmit.log_space      = 1;
    cnaEmit.coarse_grained = 1;
    printf("Collapsed CNA data to %5i sites based on potential jump events.\n", cnaEmit.total_events);
    cout<<"Precomputing for CNA..."<<flush;
    myClone.get_cnaEmitLog();
    cout<<"done."<<endl;
  }
  if (bafEmit.is_set && (opts.cna_jumps_fn != NULL || opts.baf_jumps_fn != NULL)){
    bafEmit.log_space      = 1;
    bafEmit.coarse_grained = 1;
    printf("Collapsed BAF data to %5i sites based on potential jump events.\n", bafEmit.total_events);
    cout<<"Precomputing for BAF..."<<flush;
    myClone.get_bafEmitLog();
    cout<<"done."<<endl;
  }
  if (snvEmit.is_set && opts.snv_jumps_fn != NULL){
    snvEmit.log_space      = 1;
    snvEmit.coarse_grained = 1;
    printf("Collapsed SNV data to %5i sites based on potential jump events.\n", snvEmit.total_events);
    cout<<"Precomputing for SNV..."<<flush;
    myClone.get_snvEmitLog();//TBD
    cout<<"done."<<endl;
  }
  cout<<endl;
  // get purities...
  if (opts.purity_fn != NULL){
    get_purity( opts.purity_fn, myClone.min_purity);
  }
  // get user pre-defined clones
  gsl_matrix * clones = NULL;
  gsl_vector * mass   = NULL;
  if (opts.clones_fn != NULL) get_fixed_clones( clones, mass, opts.clones_fn, nTimes);
  int bestn=0, rows=0;
  if (mass != NULL   && (int) mass->size > nTimes)    rows = (int) mass->size;
  if (clones != NULL && (int) clones->size1 > nTimes) rows = (int) clones->size1;
  if (rows > nTimes){//print LLH's for predefined parameter values...
    print_llh_for_set( clones, mass, &myClone, opts);
    return(0);//conditional end of program
  }
  else{
    // ****** INFERENCE STARTS HERE ******
    bestn = infer_clones( clones, mass, &myClone, opts);
    // ****** INFERENCE COMPLETED  ******
  }
  printf("cloneHD found support for %i sub-clone(s) in the data.\n", bestn);
  bafEmit.connect = 0;
  // NOTE: at the end, best solution is inserted
  //
  // *** PRINT CLONAL RESULTS ***
  // margin-map to get posterior per clone...
  char buff[1024]; 
  sprintf( buff, "%s.clonal.margin_map.txt", opts.pre);
  FILE * margin_fp = fopen( buff, "w");
  for (int i=0; i<(int) myClone.margin_map->size1; i++){
    for (int j=0; j<(int) myClone.margin_map->size2; j++){
      fprintf(margin_fp, "%.1f ", gsl_matrix_get(myClone.margin_map,i,j));
    }
    fprintf(margin_fp, "\n");
  }
  fclose(margin_fp);
  // used total copynumber...
  FILE * baf_utcn_fp=NULL, * snv_utcn_fp=NULL; 
  if ( bafEmit.is_set ){
    sprintf( buff, "%s.baf.used-tcn.txt", opts.pre);
    baf_utcn_fp = fopen( buff, "w");
    fprintf( baf_utcn_fp, "#sample site used-total-copynumber\n");
  }
  if (snvEmit.is_set){
    sprintf( buff, "%s.snv.used-tcn.txt", opts.pre);
    snv_utcn_fp = fopen( buff, "w");
    fprintf( snv_utcn_fp, "#sample site used-total-copynumber\n");
  }
  // print CNA posterior distributions..
  if (cnaEmit.is_set){
    sprintf( buff, "%s.cna.posterior.txt", opts.pre);
    FILE * cna_fp = fopen( buff, "w");
    print_clonal_header( cna_fp, &myClone, &cnaEmit, opts);
    sprintf( buff,"%s.copynumber.txt", opts.pre); 
    FILE * phi_fp =  fopen(buff,"w");
    fprintf( phi_fp, "#sample site total-copynumber max-copynumber\n");
    //allocate space for the posterior
    myClone.alpha_cna = new gsl_matrix * [cnaEmit.nSamples];
    myClone.gamma_cna = new gsl_matrix * [cnaEmit.nSamples];
    for (int s=0; s < cnaEmit.nSamples; s++){//print each chromosome...
      myClone.alpha_cna[s] = NULL;
      myClone.gamma_cna[s] = NULL;
      myClone.get_cna_posterior(s);
      print_posterior( cna_fp, &myClone, &cnaEmit, s, opts);
      //get, map and print total copynumber...
      myClone.get_phi(s);
      print_phi( phi_fp, &myClone, &cnaEmit, s, opts);
      int cnaChr = cnaEmit.chr[s];
      int bafSample = -1;
      int snvSample = -1;
      if (bafEmit.is_set && bafEmit.chrs.count(cnaChr) == 1){
	myClone.map_phi( &cnaEmit, s, &bafEmit);
	bafSample = bafEmit.idx_of[cnaChr];
	if (snvEmit.is_set && snvEmit.chrs.count(cnaChr) == 1){
	  myClone.map_phi( &bafEmit, bafSample, &snvEmit);
	}
      }
      else if (snvEmit.is_set && snvEmit.chrs.count(cnaChr) == 1){
	myClone.map_phi( &cnaEmit, s, &snvEmit);
	snvSample = snvEmit.idx_of[cnaChr];
      }     
      //print used total-copynumber tracks...
      if ( bafSample >= 0 && bafEmit.phi != NULL){
	print_phi( baf_utcn_fp,  &myClone, &bafEmit, bafSample, opts);
      }
      if (snvSample >= 0 && snvEmit.phi != NULL){
	print_phi( snv_utcn_fp,  &myClone, &snvEmit, snvSample, opts);
      }
    }
    fclose(cna_fp);
    fclose(phi_fp);
    if (bafEmit.is_set){//print BAF posterior...
      sprintf( buff, "%s.baf.posterior.txt", opts.pre);
      FILE * baf_fp = fopen( buff, "w");
      print_clonal_header( baf_fp, &myClone, &bafEmit, opts);
      myClone.alpha_baf = new gsl_matrix * [myClone.bafEmit->nSamples];
      myClone.gamma_baf = new gsl_matrix * [myClone.bafEmit->nSamples];    
      for (int s=0; s < bafEmit.nSamples; s++){
	myClone.alpha_baf[s] = NULL;
	myClone.gamma_baf[s] = NULL;
	myClone.get_baf_posterior(s);
	print_posterior( baf_fp, &myClone, &bafEmit, s, opts);
      }
      fclose(baf_fp);
    }
    if (snvEmit.is_set){//print SNV posterior...
      sprintf( buff, "%s.snv.posterior.txt", opts.pre);
      FILE * snv_fp = fopen( buff, "w");
      print_clonal_header( snv_fp, &myClone, &snvEmit, opts);
      myClone.alpha_snv = new gsl_matrix * [myClone.snvEmit->nSamples];
      myClone.gamma_snv = new gsl_matrix * [myClone.snvEmit->nSamples];   
      for (int s=0; s < snvEmit.nSamples; s++){
	myClone.alpha_snv[s] = NULL;
	myClone.gamma_snv[s] = NULL;
	myClone.get_snv_posterior(s);
	print_posterior( snv_fp, &myClone, &snvEmit, s, opts);
	//cleanup
	gsl_matrix_free(myClone.gamma_snv[s]);
	myClone.gamma_snv[s] = NULL;
      }
      fclose(snv_fp);
      delete [] myClone.gamma_snv;
      delete [] myClone.alpha_snv;
    }
    // clean up posterior...
    for (int s=0; s < cnaEmit.nSamples; s++) gsl_matrix_free(myClone.gamma_cna[s]);
    delete []  myClone.gamma_cna;
    myClone.gamma_cna = NULL;
    if (bafEmit.is_set){
      for (int s=0; s < bafEmit.nSamples; s++) gsl_matrix_free(myClone.gamma_baf[s]);
      delete [] myClone.gamma_baf;
      delete [] myClone.alpha_baf;
    }
  } 
  else if (snvEmit.is_set){//SNV only...
    sprintf( buff, "%s.snv.posterior.txt", opts.pre);
    FILE * snv_fp = fopen( buff, "w");
    print_clonal_header( snv_fp, &myClone, &snvEmit, opts);
    myClone.alpha_snv = new gsl_matrix * [myClone.snvEmit->nSamples];
    myClone.gamma_snv = new gsl_matrix * [myClone.snvEmit->nSamples];    
    for (int s=0; s < snvEmit.nSamples; s++){
      myClone.alpha_snv[s] = NULL;
      myClone.gamma_snv[s] = NULL;
      myClone.get_snv_posterior(s);
      print_posterior( snv_fp, &myClone, &snvEmit, s, opts);
      //cleanup
      gsl_matrix_free( myClone.gamma_snv[s]);
      myClone.gamma_snv[s] = NULL;
      //total copynumber...
      print_phi( snv_utcn_fp,  &myClone, &snvEmit, s, opts);
    }
    fclose(snv_fp);
    delete [] myClone.gamma_snv;
    delete [] myClone.alpha_snv;
  }
  if (snv_utcn_fp != NULL) fclose(snv_utcn_fp);
  if (baf_utcn_fp != NULL) fclose(baf_utcn_fp);
  //SNV BULK
  if (snvEmit.is_set && myClone.bulk_mean != NULL){
    if (opts.bulk_updates > 0) myClone.set_bulk_to_post();
    //print bulk posterior mean/distribution//TBD
    for (int t=0; t<nTimes; t++){
      sprintf( buff, "%s.bulk.posterior.%i.txt", opts.pre, t+1);
      FILE * bulk_fp = fopen( buff, "w");
      fprintf( bulk_fp, "#sample site mean");
      if (myClone.bulk_prior!=NULL){
	fprintf( bulk_fp, " std-dev dist");
      }
      fprintf( bulk_fp, "\n");
      for (int s=0; s < snvEmit.nSamples; s++){
	for (int idx=0; idx<snvEmit.nSites[s]; idx++){
	  double bmean = myClone.bulk_mean[t][s][idx];
	  fprintf( bulk_fp, "%i %i %.3f", snvEmit.chr[s], snvEmit.loci[s][idx], bmean);
	  if (myClone.bulk_prior!=NULL){
	    gsl_vector_view bdist = gsl_matrix_row( myClone.bulk_dist[t][s], idx);
	    double var = get_var( &bdist.vector, 0.0, 1.0, bmean);
	    fprintf( bulk_fp, " %.3f", sqrt(var));
	    for (int i=0; i<=snvEmit.gridSize; i++){
	      fprintf( bulk_fp, " %.3f", (&bdist.vector)->data[i]);
	    }
	  }
	  fprintf( bulk_fp, "\n");
	}	
      }
      fclose(bulk_fp);
    }
  }
  // all done...
  return (0);
}
// *** MAIN END ***




// get command line arguments...
void get_opts( int argc, const char ** argv, cmdl_opts& opts){
  default_opts(opts);
  int opt_idx = 1;
  string opt_switch;
  while ( opt_idx < argc && (argv[opt_idx][0] == '-')){
    opt_switch = argv[opt_idx];
    if ( opt_switch.compare("--print-options") == 0){
      print_opts();
      exit(0);
    }
    opt_idx++;
    if (opt_idx==argc) break;
    if ( argv[opt_idx][0] == '-' && argv[opt_idx][1] == '-') continue;
    if ( opt_switch.compare("--cna") == 0){
      opts.cna_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--baf") == 0){
      opts.baf_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--snv") == 0){
      opts.snv_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--clones") == 0){
      opts.clones_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--bulk-mean") == 0 ){
      opts.bulk_fn = argv[opt_idx];
      opts.bulk_mean=1;
    }
    else if ( opt_switch.compare("--bulk-prior") == 0 ){
      opts.bulk_fn = argv[opt_idx];
      opts.bulk_prior=1;
    }
    else if ( opt_switch.compare("--bulk-updates") == 0 ){
      opts.bulk_updates = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--mean-tcn") == 0){
      opts.mntcn_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--avail-cn") == 0){
      opts.avcn_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--bias") == 0){
      opts.bias_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--pre") == 0){
      opts.pre = argv[opt_idx];
    }
    else if ( opt_switch.compare("--chr") == 0){
      opts.chr_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--grid") == 0){
      opts.grid = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--seed") == 0){
      opts.seed = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--maxcn") == 0){
      opts.maxcn = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--maxcn-mask") == 0){
      opts.maxcn_mask_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--nmax") == 0){
      opts.nmax = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--trials") == 0){
      opts.trials = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--restarts") == 0){
      opts.restarts = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--cna-rnd") == 0){
      opts.cna_rnd = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--baf-rnd") == 0){
      opts.baf_rnd = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--snv-rnd") == 0){
      opts.snv_rnd = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--snv-err") == 0){
      opts.snv_err = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--snv-fpr") == 0){
      opts.snv_fpr = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--cna-jump") == 0){
      opts.cna_jump = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--baf-jump") == 0){
      opts.baf_jump = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--snv-jump") == 0){
      opts.snv_jump = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--bulk-sigma") == 0){
      opts.bulk_sigma = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--cna-shape") == 0){
      opts.cna_shape = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--baf-shape") == 0){
      opts.baf_shape = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--snv-shape") == 0){
      opts.snv_shape = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--baf-pen") == 0){
      opts.baf_pen = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--snv-pen") == 0){
      opts.snv_pen = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--purity") == 0){
      opts.purity_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--cna-jumps") == 0){
      opts.cna_jumps_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--baf-jumps") == 0){
      opts.baf_jumps_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--snv-jumps") == 0){
      opts.snv_jumps_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--force") == 0){
      opts.force = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--print-all") == 0){
      opts.print_all = atoi(argv[opt_idx]);
      if (opts.print_all>0) opts.print_all = 1;
    }
    else if ( opt_switch.compare("--mass-gauging") == 0){
      opts.mass_gauging = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--bulk-fix") == 0){
      opts.bulk_fix = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--min-occ") == 0){
      opts.min_occ = atof(argv[opt_idx]);
    }  
    else if ( opt_switch.compare("--min-jump") == 0){
      opts.min_jump = atof(argv[opt_idx]);
    } 
    else if ( opt_switch.compare("--learn-priors") == 0){
      opts.learn_priors = atoi(argv[opt_idx]);
      if (opts.learn_priors>0) opts.learn_priors = 1;
    }   
    else{
      cout<<"Usage:"<<endl;
      cout<<"cloneHD --print-options"<<endl;
      exit(0);
    }
    opt_idx++;
    opt_switch.clear();
  }
  test_opts(opts);
  srand(opts.seed);
}

void default_opts(cmdl_opts& opts){
  opts.cna_fn    = NULL;
  opts.baf_fn    = NULL;
  opts.snv_fn    = NULL;
  opts.bulk_fn   = NULL;
  opts.clones_fn = NULL;
  opts.bias_fn   = NULL;
  opts.mntcn_fn  = NULL;
  opts.avcn_fn   = NULL;
  opts.chr_fn    = NULL;
  opts.purity_fn = NULL;
  opts.cna_jumps_fn   = NULL;
  opts.baf_jumps_fn   = NULL;
  opts.snv_jumps_fn   = NULL;
  opts.maxcn_mask_fn   = NULL;
  opts.pre      = "./out";
  opts.grid     = 300;
  opts.cna_jump = -1.0;//negative values mean NO correlations
  opts.baf_jump = -1.0;
  opts.snv_jump = -1.0;
  opts.cna_rnd = 0.0;//random error rate
  opts.baf_rnd = 0.0;
  opts.snv_rnd = 0.0;
  opts.snv_err = -1.0;
  opts.snv_fpr = -1.0;
  opts.cna_shape = -1.0;//shape parameters
  opts.baf_shape = -1.0;
  opts.snv_shape = -1.0;
  opts.baf_pen   = -1.0;
  opts.snv_pen   = -1.0;
  opts.bulk_fix  = -1.0;//constant bulk for SNV
  opts.bulk_sigma = -1.0;
  opts.force    = -1;
  opts.trials   = 1;
  opts.nmax     = 3;
  opts.seed     = 123456 * (int(time(NULL)) % 10) + (int(time(NULL)) % 1000);
  opts.maxcn    = 4;
  opts.min_occ  = 0.01;
  opts.min_jump = 0.01;
  opts.print_all = 0;
  opts.bulk_mean    = 0;
  opts.bulk_prior   = 0;
  opts.bulk_updates = 0;
  opts.restarts = 10;
  opts.learn_priors = 1;
  opts.mass_gauging = 1;
}

void test_opts(cmdl_opts& opts){
  // *** CHECK COMMAND LINE ARGUMENTS ***
  if ( opts.cna_fn == NULL && opts.baf_fn == NULL && opts.snv_fn == NULL){
    cout<<"One of --cna [file], --baf [file] and --snv [file] must be given.\n";
    exit(1);
  }
  if (opts.snv_fn != NULL){
    if ( opts.snv_jumps_fn != NULL || opts.snv_jump >= 0.0){//with SNV persistence
      if ( opts.bulk_fn == NULL && opts.bulk_fix < 0.0){
	cout<<"With --snv [file] with correlations, one of --bulk-(prior/mean) [file] or --bulk-fix [double] must be given.\n";
	exit(1);
      }
      opts.snv_err = 0.0;
    }
    else{//no SNV persistence
      opts.bulk_fix = 0.0;
    }
    if (opts.bulk_fn != NULL && opts.bulk_fix >= 0.0){
      cout<<"Only one of --bulk-(mean/prior) [file] and --bulk-fix [double] can be used.\n";
      exit(1);
    }
    if (opts.bulk_fn != NULL && opts.bulk_mean == 1 && opts.bulk_prior == 1){
      cout<<"Only one of --bulk-mean or --bulk-prior [file] can be used.\n";
      exit(1);
    }
  }
  if (opts.cna_fn != NULL && (opts.mntcn_fn != NULL || opts.avcn_fn != NULL )){
    cout<<"--mean-tcn [file] and --avail-cn [file] cannot be used with --cna [file].\n";
    exit(1);
  }
  if ( opts.cna_fn != NULL && opts.cna_jump < 0.0 && opts.cna_jumps_fn == NULL ){
    cout<<"With --cna [file], --cna-jump [double] or --cna-jumps [file] must be given.\n";
    exit(1);
  }
  if (opts.bulk_fix == 0.0 && opts.snv_rnd == 0.0){
    opts.snv_rnd = 1.0e-9;
  }
  //some settings...
  if (opts.clones_fn != NULL){
    opts.nmax = 100;//THIS CAUSES A BUG! FIX IT!
  }
  if (opts.force > 0){
    opts.nmax = opts.force;
  }
}

void print_opts(){
  cout<<endl<<"./build/cloneHD --cna [file] --snv [file] --baf [file] --pre [string:./out] --clones [file] --purity $purity --chr [file] --bias [file] --mean-tcn [file] --avail-cn [file] --grid [int:300] --seed [int:time(0)] --trials [int:1] --restarts [int:10] --nmax [int:3] --force [int] --maxcn [int:4] --maxcn-mask [file] --cna-jump [double] --baf-jump [double] --snv-jump [double] --cna-jumps [files] --baf-jumps [file] --snv-jumps [file] --cna-rnd [double:0] --baf-rnd [double:0] --snv-rnd [double:0] --snv-err [double] --snv-fpr [double] --cna-shape [double:inf] --baf-shape [double:inf] --snv-shape [double:inf] --baf-pen [double:1.0] --snv-pen [double:0.01] --min-occ [double:0.01] --min-jump [double:0.01] --learn-priors [0/1:0] --mass-gauging [0/1:1] ";
  //cout<<"--bulk-prior [file] --bulk-mean [file] --bulk-fix [double:0] --bulk-sigma [double] --bulk-updates [int:0]";
  cout<<endl;
  exit(0);
}

void print_clonal_header(FILE * fp, Clone * myClone, Emission * myEmit, cmdl_opts& opts){
  fprintf( fp,"#copynumbers:\n");
  for (int k=0; k < myClone->nLevels; k++){
    for (int f=0; f<myClone->nClones; f++){
      fprintf( fp, "%i", myClone->copynumber[k][f]);
    }
    fprintf( fp," ");
  }
  fprintf( fp,"\n#clonal frequencies:\n");
  for (int t=0; t<myClone->nTimes; t++){
    for (int k=0; k<myClone->nLevels; k++){
      fprintf( fp, "%.3e ", myClone->clone_spectrum[t][k]);
    }
    fprintf( fp, "\n");
  }
  if (opts.print_all || myEmit->coarse_grained == 0 ){//HERE!!!
    fprintf( fp,"#sample locus PostDist\n");
  }
  else{
    fprintf( fp,"#sample first-locus nloci last-locus PostDist\n");
  }
}


void print_posterior( FILE * post_fp, Clone * myClone, Emission * myEmit, int s, cmdl_opts& opts){
  for (int evt=0; evt < myEmit->nEvents[s]; evt++){   
    int first = myEmit->idx_of_event[s][evt];     
    int last  = (evt < myEmit->nEvents[s]-1) ?  myEmit->idx_of_event[s][evt+1] - 1 : myEmit->nSites[s]-1; 
    gsl_matrix * post=NULL;
    if ( myEmit == myClone->cnaEmit) post =  myClone->gamma_cna[s];
    if ( myEmit == myClone->bafEmit) post =  myClone->gamma_baf[s];
    if ( myEmit == myClone->snvEmit) post =  myClone->gamma_snv[s];
    if( opts.print_all || myEmit->coarse_grained == 0 ){
      for (int idx=first; idx<=last; idx++){
      	fprintf( post_fp, "%i %6i", myEmit->chr[s], myEmit->loci[s][idx]);
        for (int j=0; j< myClone->nLevels; j++){
          fprintf( post_fp, " %.3f", gsl_matrix_get( post, evt, j));
      	}
      	fprintf( post_fp,"\n");
      }
    }
    else{
      fprintf( post_fp, "%i %6i %6i %6i", myEmit->chr[s], 
	       myEmit->loci[s][first], last-first+1, myEmit->loci[s][last]
	       );
      for (int j=0; j< myClone->nLevels; j++){
	fprintf( post_fp, " %.3f", gsl_matrix_get( post, evt, j));
      }
      fprintf( post_fp,"\n");
    }
  }
}

void print_mean_tcn( FILE * mntcn_fp, Clone * myClone, Emission * myEmit, int s, cmdl_opts& opts){
  double ncn = double(myClone->normal_copy[myEmit->chr[s]]);
  for (int evt=0; evt < myEmit->nEvents[s]; evt++){   
    int first = myEmit->idx_of_event[s][evt];     
    int last  = (evt < myEmit->nEvents[s]-1) ?  myEmit->idx_of_event[s][evt+1] - 1 : myEmit->nSites[s]-1; 
    if( opts.print_all || myEmit->coarse_grained == 0 ){
      for (int idx=first; idx<=last; idx++){
      	fprintf( mntcn_fp, "%i %6i", myEmit->chr[s], myEmit->loci[s][idx]);
	for (int t=0; t<myEmit->nTimes; t++){
	  double tcn = (myEmit->mntcn==NULL) ? ncn : myEmit->mntcn[t][s][evt];
	  fprintf( mntcn_fp, " %.3f", tcn);
	}
	//int mcn = (myEmit->cnmax==NULL) ? int(ncn) : myEmit->cnmax[s][evt];//???
	//fprintf( phi_fp, " %i\n", mcn);
	fprintf( phi_fp, "\n");
      }
    }
    else{
      fprintf( phi_fp, "%i %6i %6i %6i", myEmit->chr[s], 
	       myEmit->loci[s][first], last-first+1, myEmit->loci[s][last]
	       );
      for (int t=0; t<myEmit->nTimes; t++){
	double tcn = (myEmit->mntcn==NULL) ? ncn : myEmit->mntcn[t][s][evt];
	fprintf( mntcn_fp, " %.3f", tcn);
      }
      //int mcn = (myEmit->cnmax==NULL) ? int(ncn) : myEmit->cnmax[s][evt];//???
      //fprintf( phi_fp, " %i\n", mcn);
      fprintf( phi_fp, "\n");
    }
  }
}

//TODO
//void print_avail_cn( FILE * avcn_fp, Clone * myClone, int s, cmdl_opts& opts){}
