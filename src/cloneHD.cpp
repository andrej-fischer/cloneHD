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
#include "cloneHD-inference.h"
#include "common-functions.h"
#include "clone.h"
#include "emission.h"


using namespace std;



// *** MAIN START***
int main (int argc, const char * argv[]){
  cmdl_opts opts;
  get_opts( argc, argv, opts);
  int nTimes=0, nT=0;
  //*** EMITTED DATA OBJECTS ***
  Emission cnaEmit, bafEmit, snvEmit;
  if (opts.cna_fn != NULL) get_cna_data( &cnaEmit, opts, nTimes);
  if (opts.baf_fn != NULL) get_baf_data( &bafEmit, opts, nTimes, nT);
  if (opts.snv_fn != NULL) get_snv_data( &snvEmit, opts, nTimes, nT);
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
  myClone.baf_pen  = opts.baf_pen;//BAF penalty for complex chr status
  myClone.snv_fpr  = opts.snv_fpr;//SNV false-positive rate
  myClone.snv_fpf  = opts.snv_fpf;//SNV frequency of false positives
  myClone.snv_pen  = (bafEmit.is_set ? opts.snv_pen : 0.5);//penalty for high SNV genotypes
  myClone.bulk_fix = opts.bulk_fix;
  myClone.cnaGrid  = opts.cnaGrid;
  myClone.bafGrid  = opts.bafGrid;
  myClone.snvGrid  = opts.snvGrid;
  myClone.bulkGrid = opts.bulkGrid;
  // *** GET MAX-TCN INFO ***
  get_maxtcn_input( opts.maxtcn_fn, opts.maxtcn, &myClone);
  /*std::map<int, vector<int> >::iterator it;
  for (it=myClone.maxtcn_input.begin(); it != myClone.maxtcn_input.end(); it++){
    printf("%i ",it->first);
    for (int j=0; j<it->second.size(); j++) printf("%i ", it->second[j]);
    cout<<endl;
  }
  exit(1);
  */
  // *** GET SNV BULK PRIOR ***
  if ( snvEmit.is_set && opts.bulk_fn != NULL ){
    printf("Using data in %s as SNV bulk prior...\n", opts.bulk_fn);
    get_snv_bulk_prior( &myClone, opts);
  }
  //*** GET JUMP PROBABILITY TRACKS and COLLAPSE TO EVENTS***
  get_jump_probability( &myClone, opts);
  //...now all segments are fixed and mean_tcn/av_cn allocated.
  if ( snvEmit.is_set && !cnaEmit.is_set ){//for SNV only
    // *** GET TOTAL MEAN COPYNUMBER TRACK ***  
    if( opts.mntcn_fn != NULL ){
      get_mean_tcn( opts.mntcn_fn, &myClone, &snvEmit);
    } 
    // *** GET AVAILABLE COPYNUMBER TRACK ***  
    if ( opts.avcn_fn != NULL ){
      get_avail_cn( opts.avcn_fn, &myClone, &snvEmit);
    }
  }
  //*** GET READ DEPTH BIAS FIELD ***
  if (cnaEmit.is_set && opts.bias_fn != NULL){
    get_bias_field( &myClone, opts);
  }
  //*** PREPARE COARSE-GRAINED DATA ***
  if (cnaEmit.is_set && opts.cna_jumps_fn != NULL){
    cnaEmit.log_space      = 1;
    cnaEmit.coarse_grained = 1;
    printf("Collapsed CNA data to %5i segments based on potential jump events.\n", cnaEmit.total_events);
    cout<<"Precomputing for CNA..."<<flush;
    myClone.get_cnaEmitLog();
    cout<<"done."<<endl;
  }
  if (bafEmit.is_set && (opts.cna_jumps_fn != NULL || opts.baf_jumps_fn != NULL)){
    bafEmit.log_space      = 1;
    bafEmit.coarse_grained = 1;
    printf("Collapsed BAF data to %5i segments based on potential jump events.\n", bafEmit.total_events);
    cout<<"Precomputing for BAF..."<<flush;
    myClone.get_bafEmitLog();
    cout<<"done."<<endl;
  }
  if (snvEmit.is_set && opts.snv_jumps_fn != NULL){
    snvEmit.log_space      = 1;
    snvEmit.coarse_grained = 1;
    printf("Collapsed SNV data to %5i segments based on potential jump events.\n", snvEmit.total_events);
    cout<<"Precomputing for SNV..."<<flush;
    myClone.get_snvEmitLog();
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
    return(0);
  }
  else{
    // ****** INFERENCE STARTS HERE ******
    bestn = infer_clones( clones, mass, &myClone, opts);
    printf("cloneHD found support for %i sub-clone(s) in the data.\n", bestn);
    // ****** INFERENCE COMPLETED ********
  }
  print_all_results( &myClone, opts);
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
    else if ( opt_switch.compare("--max-tcn") == 0){
      if ( isdigit(argv[opt_idx][0]) ){
	opts.maxtcn = atoi(argv[opt_idx]);
      }
      else{
	opts.maxtcn_fn = argv[opt_idx];
      }
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
    else if ( opt_switch.compare("--cna-grid") == 0){
      opts.cnaGrid = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--baf-grid") == 0){
      opts.bafGrid = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--snv-grid") == 0){
      opts.snvGrid = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--bulk-grid") == 0){
      opts.bulkGrid = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--seed") == 0){
      opts.seed = atoi(argv[opt_idx]);
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
    else if ( opt_switch.compare("--snv-fpfreq") == 0){
      opts.snv_fpf = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--snv-fprate") == 0){
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
      if (opts.print_all > 0) opts.print_all = 1;
    }
    else if ( opt_switch.compare("--mass-gauging") == 0){
      opts.mass_gauging = atoi(argv[opt_idx]);
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
      cout<<"Usage see:"<<endl;
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
  opts.maxtcn_fn  = NULL;
  opts.avcn_fn   = NULL;
  opts.chr_fn    = NULL;
  opts.purity_fn = NULL;
  opts.cna_jumps_fn   = NULL;
  opts.baf_jumps_fn   = NULL;
  opts.snv_jumps_fn   = NULL;
  opts.pre      = "./out";
  opts.cnaGrid     = 300;
  opts.bafGrid     = 100;
  opts.snvGrid     = 100;
  opts.bulkGrid    = 100;
  opts.cna_jump = -1.0;//negative values mean NO correlations
  opts.baf_jump = -1.0;
  opts.snv_jump = -1.0;
  opts.cna_rnd = 0.0;//random error rate
  opts.baf_rnd = 0.0;
  opts.snv_rnd = 0.0;
  opts.snv_fpf = 0.01;
  opts.snv_fpr = 0.001;
  opts.cna_shape = -1.0;//shape parameters
  opts.baf_shape = -1.0;
  opts.snv_shape = -1.0;
  opts.baf_pen   = 1.0;
  opts.snv_pen   = 0.01;
  opts.bulk_fix   = -1.0;//constant bulk for SNV
  opts.bulk_sigma = -1.0;
  opts.force    = -1;
  opts.trials   = 1;
  opts.nmax     = 3;
  opts.seed     = 123456 * (int(time(NULL)) % 10) + (int(time(NULL)) % 1000);
  opts.maxtcn    = 4;
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
    cout<<"ERROR: One of --cna [file], --baf [file] and --snv [file] must be given.\n";
    exit(1);
  }
  if (opts.snv_fn != NULL){
    if ( opts.snv_jumps_fn != NULL || opts.snv_jump >= 0.0){//with SNV persistence
      if ( opts.bulk_fn == NULL && opts.bulk_fix < 0.0){
	cout<<"ERROR: With --snv [file] with correlations, one of --bulk-(prior/mean) [file] or --bulk-fix [double] must be given.\n";
	exit(1);
      }
      opts.snv_fpf = 0.0;//there are no false positives
    }
    else{//no SNV persistence
      opts.bulk_fix = 0.0;
    }
    if (opts.bulk_fn != NULL && opts.bulk_fix >= 0.0){
      cout<<"ERROR: Only one of --bulk-(mean/prior) [file] and --bulk-fix [double] can be used.\n";
      exit(1);
    }
    if (opts.bulk_fn != NULL && opts.bulk_mean == 1 && opts.bulk_prior == 1){
      cout<<"ERROR: Only one of --bulk-mean or --bulk-prior [file] can be used.\n";
      exit(1);
    }
  }
  if (opts.cna_fn != NULL && (opts.mntcn_fn != NULL || opts.avcn_fn != NULL )){
    cout<<"ERROR: --mean-tcn [file] and --avail-cn [file] cannot be used with --cna [file].\n";
    exit(1);
  }
  if ( opts.cna_fn != NULL && opts.cna_jump < 0.0 && opts.cna_jumps_fn == NULL ){
    cout<<"ERROR: With --cna [file], --cna-jump [double] or --cna-jumps [file] must be given.\n";
    exit(1);
  }
  if (opts.bulk_fix == 0.0 && opts.snv_rnd == 0.0){
    opts.snv_rnd = 1.0e-9;
  }
  if (opts.force > 0){
    opts.nmax = opts.force;
  }
}

void print_opts(){
  cout<<endl<<"./build/cloneHD --cna [file] --snv [file] --baf [file] --pre [string:./out] --clones [file] --purity [file] --chr [file] --bias [file] --max-tcn [file/int:4] --mean-tcn [file] --avail-cn [file] --seed [int:time(0)] --trials [int:1] --restarts [int:10] --nmax [int:3] --force [int] --cna-jump [double] --baf-jump [double] --snv-jump [double] --cna-jumps [files] --baf-jumps [file] --snv-jumps [file] --cna-rnd [double:0] --baf-rnd [double:0] --snv-rnd [double:0] --snv-fpfreq [double:0.01] --snv-fprate [double:0.001] --cna-shape [double:inf] --baf-shape [double:inf] --snv-shape [double:inf] --baf-pen [double:1.0] --snv-pen [double:0.01] --min-occ [double:0.01] --min-jump [double:0.01] --learn-priors [0/1:0] --mass-gauging [0/1:1] --bulk-prior [file] --bulk-mean [file] --bulk-fix [double:0] --bulk-sigma [double] --bulk-updates [int:0] --cna-grid [int:300] --snv-grid [int:100] --baf-grid [int:100]";
  cout<<endl;
  exit(0);
}
