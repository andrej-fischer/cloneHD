/*
  filterHD.cpp
 
  Input: a file with four columns:
  chr locus reads depth

  Each chromosome is treated as an independent sample from the same process.
  The hidden emission rate state sequence is modeled as a jump-difusion trajectory.
  Its value x at any point in time is the mean rate of a integer-valued emission process that 
  generates the observed data.

  Output:
  The posterior mean and standard deviation of the emission rate.
  The posterior jump probabiliy for each transition.
  On request (--dist), the posterior distribution at each point (may be large files!).

  Required Arguments:
  --data        The input observed data file
  --mode        Emission model:
                1 Binomial
		2 Beta-Binomial
		3 Poisson
		4 Negative Binomial

  Optional Arguments:
  --pre         The prefix to put before all output files
  --grid [int]  The grid size for the distributions (partitions [0,1] into [grid] bins).
  --dist        Prints all the posterior distributions as well.
  --jump [double]  To fix the jump probability per base
  --sigma [double] To fix the diffusion constant
  --rnd [double]   To fix the random emission rate
  --shape [double] To fix the shape parameter in mode 2 and 4
  --bias [file]    A bias field modulating Poisson emissions. Use filterHD output file style 
                   (posterior distribution not needed). Bias field is assumed not to have jumps.
                   With bias field, filterHD will set --sigma 0.0.

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
#include "emission.h"
#include "jump-diffusion.h"
#include "minimization.h"
#include "common-functions.h"

using namespace std;


struct cmdl_opts{
  const char * data_fn;
  const char * bias_fn;
  const char * pre;
  int grid, dist, nojump, mode, seed, filter_pVal,filter_shortSeg, jumps, reflect;
  double sigma, jump, rnd_emit, shape,xmin, xmax, min_jump;
  double sigma_i, jump_i, rnd_emit_i, shape_i;
};


//*** OWN FUNCTIONS ***
void get_opts( int argc, const char ** argv, cmdl_opts& opts);
void print_opts();
double Q( const gsl_vector * x, void * p);
double find_JD_parameters(JumpDiffusion * myJD, cmdl_opts& opts);
struct fpar{
  JumpDiffusion * myJD;
  vector<int> to_opt;
};
void init_parameter(gsl_vector*& var, gsl_vector*& range, vector<int>& to_opt, JumpDiffusion * myJD, cmdl_opts& opts);

// *** MAIN START***
int main (int argc, const char * argv[]){
  cmdl_opts opts;
  get_opts( argc, argv, opts);
  srand(opts.seed);
  vector<int> chrs;
  vector<int> nSites;
  int nTimes;
  get_dims( opts.data_fn, nTimes, chrs, nSites);
  int nSamples = (int) chrs.size();
  int total_nLoci=0;
  for (int s=0; s<nSamples; s++){
    total_nLoci += nSites[s];
  }
  //announce:
  printf("\nfilterHD: Fitting jump-diffusion model to data at %i loci ", total_nLoci);
  printf("in %i segment(s) in %i sample(s) with a ", nSamples, nTimes);
  if (opts.mode == 1){
    printf("binomial ");
  }
  else if (opts.mode == 2){
    printf("beta-binomial ");
  }
  else if (opts.mode == 3){
    printf("poisson ");
  }
  else if (opts.mode == 4){
    printf("negative-binomial ");
  }
  printf("emission model:\n");
  //the data and emission model object
  Emission myEmit;
  myEmit.mode    = opts.mode;
  myEmit.shape   = opts.shape;
  myEmit.get_log = 1;
  myEmit.reflect = opts.reflect;
  myEmit.set( nTimes, chrs, nSites, opts.grid);
  get_data( opts.data_fn, &myEmit);
  // *** BIAS FIELD *** emission bias field
  if (opts.bias_fn != NULL){
    myEmit.allocate_bias();
    get_bias( opts.bias_fn, &myEmit);
  }
  // *** POSTERIOR JUMP TRACK ACROSS ALL TIME POINTS ***
  double ** jumps = NULL;
  if (opts.jumps==1){
    jumps = new double * [nSamples];
    for (int s=0; s<nSamples; s++){
      jumps[s] = new double [nSites[s]];
      for (int l=0; l<nSites[s]; l++){
	jumps[s][l] = 1.0;
      }
    }
  }
  int ** mask = NULL;
  if (opts.filter_pVal || opts.filter_shortSeg > 0){
    mask = new int * [nSamples];
    for (int s=0; s<nSamples; s++){
      mask[s] = new int [nSites[s]];
      for (int l=0; l<nSites[s]; l++){
	mask[s][l] = 1;
      }
    }
  }
  //***Jump-Diffusion filtering of each sample***
  for (int t=0; t<nTimes; t++){
    printf("\nFiltering sample %i of %i:\n", t+1, nTimes);
    //the jump diffusion propagation object
    JumpDiffusion myJD( &myEmit, t);
    //find maximum-likelihood estimates of all parameters
    double llh = find_JD_parameters( &myJD, opts);
    printf("Filtered sample %i of %i: llh = %.5e, --jump %.3e --sigma %.3e --rnd %.3e",
	   t+1, nTimes, llh, myJD.jump, myJD.sigma, myJD.rnd_emit);
    if (opts.mode == 2 || opts.mode==4) printf(" --shape %.3e", myJD.myEmit->shape);
    cout<<endl;
    if ( (opts.mode == 2 || opts.mode==4) && myJD.myEmit->shape > 1.0e3){
      printf("With --shape %.3e, you might consider choosing mode %i.\n", 
	     myJD.myEmit->shape, opts.mode==2 ? 1 : 3);
    } 
    //filter out data points which are not compatible with the emission model
    if ( opts.reflect == 0 && (opts.filter_pVal|| opts.filter_shortSeg > 0) ){
      double crit=0.0;
      for (int s=0; s < myJD.nSamples; s++) crit += double(myJD.nSites[s]);
      crit = 10.0/crit; 
      for (int s=0; s < myJD.nSamples; s++){
	myJD.get_posterior(s);
	double * mean = new double [myJD.nSites[s]];
	double msd    = 0.0;
	for (int l=0; l < myJD.nSites[s]; l++){
	  gsl_vector_view post =  gsl_matrix_row(myJD.gamma[s],l);
	  mean[l] = get_mean( &post.vector, myJD.myEmit->xmin, myJD.myEmit->xmax);
	  msd += sqrt( get_var( &post.vector, myJD.myEmit->xmin, myJD.myEmit->xmax, mean[l]) );
	}
	msd /= double(myJD.nSites[s]);
	if (opts.filter_pVal){
	  for (int l=0; l < myJD.nSites[s]; l++){
	    //gsl_vector_view post =  gsl_matrix_row(myJD.gamma[s],l);
	    //double mean = get_mean( &post.vector, myJD.myEmit->xmin, myJD.myEmit->xmax);
	    double pval = myEmit.get_pval( t, s, l, mean[l]);
	    if (pval < crit) mask[s][l] = 0;
	  }
	}
	gsl_matrix_free(myJD.gamma[s]);
	myJD.gamma[s] = NULL;
	if (opts.filter_shortSeg > 0){
	  int last=0;
	  for (int l=0; l < myJD.nSites[s]; l++){
	    if ( l>0 && fabs(mean[l]-mean[l-1]) > 4.0*msd){
	      if ( l < last + opts.filter_shortSeg) 
		for (int i=last;i<l;i++) mask[s][i] = 0;
	      last = l;
	    }
	  }
	}
	delete [] mean;
      }
    }
    // print posterior information to file
    char buff[1024];
    sprintf(buff,"%s.posterior.%i.txt", opts.pre, t+1);
    FILE * total_fp = fopen(buff,"w");
    fprintf(total_fp, "#sample site mean std-dev jump-prob");
    fprintf(total_fp, " posterior %.5e %.5e\n", myJD.myEmit->xmin,  myJD.myEmit->xmax);
    int uidx  = opts.reflect ? int(0.5*double( myJD.gridSize)) :  myJD.gridSize;
    double mx = opts.reflect ? 0.5 : myJD.myEmit->xmax;
    double mn = myJD.myEmit->xmin;
    gsl_vector * post = gsl_vector_alloc(uidx+1);
    for (int s=0; s < myJD.nSamples; s++){
      //get posterior distribution with the ML parameters
      myJD.get_posterior(s);
      for (int l=0; l < myJD.nSites[s]; l++){
	if (opts.reflect){//distribution in lower half
	  gsl_vector_view lower = gsl_matrix_subrow( myJD.gamma[s], l, 0, uidx+1);
	  gsl_vector_memcpy( post, &lower.vector);
	  double norm = gsl_blas_dasum(post);
	  norm = (norm - 0.5*(post->data[0] + post->data[uidx])) * myJD.myEmit->dx;
	  if (norm <= 0.0) abort();
	  gsl_vector_scale(post,1.0/norm);
	}
	else{
	  gsl_matrix_get_row( post, myJD.gamma[s], l);
	}
	double mean = get_mean( post, mn, mx);
	double var  = get_var(  post, mn, mx, mean);
	if (opts.jumps==1){
	  jumps[s][l] *= exp(myJD.pnojump[s][l]);
	}
	fprintf(total_fp, "%i %6i %.2e %.2e %.2e", 
		chrs[s], myJD.loci[s][l], mean, sqrt(var), exp(myJD.pjump[s][l]));
	if(opts.dist==1){// full posterior distribution? LARGE!
	  for (int i=0; i <= myJD.gridSize; i++){
	    double p = gsl_matrix_get( myJD.gamma[s], l, i);
	    fprintf(total_fp, " %.2e", p); 
	  }
	}
	fprintf(total_fp,"\n");
      }
      gsl_matrix_free(myJD.gamma[s]);
      myJD.gamma[s] = NULL;
    }
    fclose(total_fp);
    gsl_vector_free(post);
    // *** RESET ***
    myEmit.delete_old_Emit();
    myEmit.range_set = 0;
    myEmit.reset_mask();
  }
  // *** PRINT JUMPS ***
  if (opts.jumps == 1){
    for (int s=0; s < nSamples; s++){
      for (int l=0; l < nSites[s]; l++){
	myEmit.pjump[s][l] = 1.0 - jumps[s][l];
      }
      if (opts.min_jump > 0.0){
	myEmit.coarse_grain_jumps( s, opts.min_jump, 5);
      }
      delete [] jumps[s];
    }
    delete [] jumps;
    char buff[1024];
    sprintf(buff,"%s.jumps.txt", opts.pre);
    FILE * jumps_fp = fopen(buff,"w");
    fprintf(jumps_fp, "#sample site jump-prob\n");
    for (int s=0; s < nSamples; s++){
      for (int l=0; l < nSites[s]; l++){
	fprintf(jumps_fp, "%i %6i %.2e\n", chrs[s], myEmit.loci[s][l], myEmit.pjump[s][l]);
      }
    }
    fclose(jumps_fp);
  }
  //FILTERED?
  if (opts.filter_pVal || opts.filter_shortSeg > 0){
    //filter short segments for reflected data
    if ( opts.reflect && opts.filter_shortSeg > 0){
      for (int s=0; s < myEmit.nSamples; s++){
	int last=0;
	for (int l=0; l < myEmit.nSites[s]; l++){
	  if ( l>0 && myEmit.pjump[s][l] > 0.01 ){
	    if ( l < last + opts.filter_shortSeg) 
	      for (int i=last;i<l;i++) mask[s][i] = 0;
	    last = l;
	  }
	}
      }
    }
    //print
    char buff[1024];  
    sprintf(buff,"%s.filtered.txt", opts.pre);
    FILE * filtered_fp = fopen(buff,"w");
    for (int s=0; s < myEmit.nSamples; s++){
      for (int l=0; l < myEmit.nSites[s]; l++){
	if( mask[s][l] == 1 ){
	  fprintf( filtered_fp, "%i %6i", chrs[s], myEmit.loci[s][l]);
	  for (int t=0;t<nTimes;t++) 
	    fprintf( filtered_fp, " %3i %3i", myEmit.reads[t][s][l], myEmit.depths[t][s][l]);
	  fprintf( filtered_fp, "\n");
	}
      }
      delete [] mask[s];
    }
    delete [] mask;
    fclose(filtered_fp);
  }
  //done
  return (0);
}
// *** MAIN END ***


// get command line arguments...
void get_opts( int argc, const char ** argv, cmdl_opts& opts){
  int opt_idx = 1;
  string opt_switch;
  opts.data_fn  = NULL;
  opts.bias_fn  = NULL;
  opts.pre      = "./out";
  opts.grid     = 100;
  opts.dist     = 0;
  opts.sigma    = -1.0;
  opts.jump     = -1.0;
  opts.rnd_emit = -1.0;
  opts.shape    = -1.0;
  opts.sigma_i    = -1.0;
  opts.jump_i     = -1.0;
  opts.rnd_emit_i = -1.0;
  opts.shape_i    = -1.0;
  opts.mode     = 0;
  opts.xmin = -1.0;
  opts.xmax = -1.0;
  opts.seed = (int) time(NULL);
  opts.filter_pVal     = 0;
  opts.filter_shortSeg = 0;
  opts.jumps  = 0;
  opts.reflect = 0;
  opts.min_jump = 0.0;
  while ( opt_idx < argc && (argv[opt_idx][0] == '-')){
    opt_switch = argv[opt_idx];
    if ( opt_switch.compare("--print-options") == 0){
      print_opts();
      exit(0);
    }
    opt_idx++;
    if (opt_idx==argc) break;
    if ( argv[opt_idx][0] == '-') continue;
    if ( opt_switch.compare("--data") == 0){//the input data
      opts.data_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--bias") == 0){//the input data
      opts.bias_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--pre") == 0){//the prefix for all output files
      opts.pre = argv[opt_idx];
    }
    else if ( opt_switch.compare("--grid") == 0){//the size of the grid for the continuous distributions
      opts.grid = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--seed") == 0){//random seed
      opts.seed = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--mode") == 0){//emission model: see above
      opts.mode = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--sigma") == 0){//diffusion constant
      opts.sigma = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--rnd") == 0){//random emission rate
      opts.rnd_emit = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--jump") == 0){//jump probability per base
      opts.jump = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--shape") == 0){//shape parameter for mode 2/4
      opts.shape = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--sigmai") == 0){//diffusion constant
      opts.sigma_i = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--rndi") == 0){//random emission rate
      opts.rnd_emit_i = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--jumpi") == 0){//jump probability per base
      opts.jump_i = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--shapei") == 0){//shape parameter for mode 2/4
      opts.shape_i = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--xmin") == 0){//shape parameter for mode 2/4
      opts.xmin = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--xmax") == 0){//shape parameter for mode 2/4
      opts.xmax = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--dist") == 0){//whether to print posterior
      opts.dist = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--filter-pVal") == 0){//whether to filter out some data points
      opts.filter_pVal = 1;
    }
    else if ( opt_switch.compare("--filter-shortSeg") == 0){//whether to filter out some data points
      opts.filter_shortSeg = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--jumps") == 0){//whether to filter out some data points
      opts.jumps = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--min-jump") == 0){//whether to filter out some data points
      opts.min_jump = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--reflect") == 0){//whether to filter out some data points
      opts.reflect = atoi(argv[opt_idx]);
    }
    else {
      cout << "Usage: filterHD --print-options"<<endl;
      exit(1);
    }
    opt_switch.clear();
    opt_idx++;
  }
  if (opts.mode==0){
    cout<<"ERROR: choose emission mode with --mode [1,2,3,4]\n";
    exit(1);
  }
  if (opts.bias_fn != NULL){
    if( opts.filter_pVal || opts.filter_shortSeg){
      cout<<"ERROR: --bias [file] and --filter-pVal or --filter-shortSeg [int] cannot be used together\n";
      exit(1);
    }
  }
  if(opts.reflect==1 && (opts.mode==3||opts.mode==4)){
     cout<<"ERROR: --reflect [0/1] can only be used in mode 1 and 2.\n";
     exit(1);
  }
  if (opts.filter_shortSeg > 0 && opts.min_jump == 0.0 && opts.reflect == 1) opts.min_jump = 1.0e-5;
}

void print_opts(){
  cout<<"# All filterHD command line options to be used in shell scripts"<<endl;
  cout<<"# the input data"<<endl;
  cout<<"data=[file]"<<endl;
  cout<<"# the emission mode (1=Binomia, 2=Beta-Binomial, 3=Poisson, 4=Negative-Binomial)"<<endl;
  cout<<"mode=[1/2/3/4]"<<endl;
  cout<<"# the output prefix (default: ./out)"<<endl;
  cout<<"pre=[string]"<<endl;
  cout<<"# the file to use for a bias field (filterHD posterior output)"<<endl;
  cout<<"bias=[file]"<<endl;
  //
  cout<<"# fix the jump probability (per base)"<<endl;
  cout<<"jump=[double]"<<endl;
  cout<<"# fix the diffusion constant (per base)"<<endl;
  cout<<"sigma=[double]"<<endl;
  cout<<"# fix the random error rate (optional)"<<endl;
  cout<<"rnd=[double]"<<endl;
  cout<<"# fix the shape parameter (optional)"<<endl;
  cout<<"shape=[double]"<<endl;
  //
  cout<<"# initialize the jump probability (per base)"<<endl;
  cout<<"jumpi=[double]"<<endl;
  cout<<"# initialize the diffusion constant (per base)"<<endl;
  cout<<"sigmai=[double]"<<endl;
  cout<<"# initialize the random error rate"<<endl;
  cout<<"rndi=[double]"<<endl;
  cout<<"# initialize the shape parameter"<<endl;
  cout<<"shapei=[double]"<<endl;
  //
  cout<<"# the grid size"<<endl;
  cout<<"grid=[int:100]"<<endl;
  cout<<"# whether to print posterior distributions"<<endl;
  cout<<"dist=[0/1:0]"<<endl;
  cout<<"# whether to print posterior jump probabilities"<<endl;
  cout<<"jumps=[0/1:0]"<<endl;
  cout<<"# whether to filter sites that are random "<<endl;
  cout<<"filter-pVal=[0/1:0]"<<endl;
  cout<<"# whether to filter sites in short segments "<<endl;
  cout<<"filter-shortSeg=[int:0]"<<endl;
  //
  cout<<"# fix lower end of range for hidden data"<<endl;
  cout<<"xmin=[double]"<<endl;
  cout<<"# fix upper end of range for hidden data"<<endl;
  cout<<"xmax=[double]"<<endl;
  //
  cout<<"cmd=\"filterHD --data $data --mode $mode --bias $bias --pre $pre --grid $grid --jump $jump --sigma --sigma --rnd $rnd --shape $shape --jumps $jumps --dist $dist --filter-pVal $filterpVal --filter-shortSeg $filtershortSeg --jumpi $jumpi --sigmai $sigmai --rndi $rndi --shapei $shapei --xmin $xmin --xmax $xmax"<<endl;
  cout<<"echo $cmd"<<endl;
  cout<<"$cmd"<<endl;
  exit(0);
}



// *** LEARN PARAMETERS OF THE JUMP-DIFFUSION MODEL ***
double find_JD_parameters(JumpDiffusion * myJD, cmdl_opts& opts){
  vector<int> to_opt;
  int nvar=0;
  double llh=0;
  if (opts.xmin >= 0.0 && opts.xmax >= 0.0){
    myJD->myEmit->xmin = opts.xmin;
    myJD->myEmit->xmax = opts.xmax;
    myJD->myEmit->ymin = opts.xmin;
    myJD->myEmit->ymax = opts.xmax;
    myJD->myEmit->set_grid();
  }
  else{
    myJD->myEmit->init_range(myJD->time);//get the initial range of rates
  }
  printf("Initial range is %.3e < x < %.3e", myJD->myEmit->xmin, myJD->myEmit->xmax);
  if (opts.bias_fn != NULL){
    printf(", %.3e < y < %.3e", myJD->myEmit->ymin, myJD->myEmit->ymax);
  }
  cout<<endl;
  if ( opts.jump < 0.0){// jump probability
    to_opt.push_back(0);
    nvar++;
  }
  else{
    myJD->jump = opts.jump;
    myJD->set_pstay();
  }
  if (opts.sigma < 0.0){// diffusion constant
    to_opt.push_back(1);
    nvar++;
  }
  else{
    myJD->sigma = opts.sigma;
    myJD->get_DiffProp();
  }
  if (opts.rnd_emit < 0.0){//random error rate
    to_opt.push_back(2);
    nvar++;
  }
  else{
    myJD->rnd_emit = opts.rnd_emit;
  }
  if ( opts.shape < 0.0 && (opts.mode == 2 || opts.mode == 4)){//shape parameter
    to_opt.push_back(3);
    nvar++;
  }
  else{
    myJD->myEmit->shape = opts.shape;
  }
  if(nvar>0){
    gsl_vector * var   = gsl_vector_calloc(nvar);
    gsl_vector * range = gsl_vector_calloc(nvar);
    //set initial values
    init_parameter(var, range,  to_opt, myJD, opts);
    fpar myfpar;
    myfpar.myJD    = myJD;
    myfpar.to_opt  = to_opt;
    void * param = static_cast<void*>(&myfpar);
    // get the ML estimates and ML value
    int steps = 0;
    gsl_vector ** simplex = NULL;
    gsl_vector * lower    = NULL;
    //header...
    printf("%-4s ", "eval");
    for (int i=0; i<nvar; i++){
      if (to_opt[i] == 0){//jump probability
	printf("%-11s ", "jump");
      }
      else if(to_opt[i] == 1){//diffusion constant
	printf("%-11s ", "sigma");
      }
      else if(to_opt[i] == 2){//random rate
	printf("%-11s ", "rnd");
      }
      else if(to_opt[i] == 3){//shape parameter
	printf("%-11s ", "shape");
      }
    }
    printf("-llh\n");
    llh = - find_local_optimum( 0, simplex, lower, var, range,  param, &Q, 1.0e-3, steps, 1);
    //adapt the range if needed
    if ((opts.mode==3 || opts.mode==4) && (opts.xmin < 0.0 && opts.xmax < 0.0)  && opts.bias_fn == NULL){
      int redo = myJD->adapt_range();
      if (redo==1){
	printf("Adapted range to %.3e < x < %.3e\n", myJD->myEmit->xmin, myJD->myEmit->xmax);
	//init_parameter( var, range,  to_opt, myJD, opts);
	llh = - find_local_optimum( 0, simplex, lower, var, range, param, &Q, 1.0e-3, steps, 1);
      }
    }
    //set the ML values into the objects
    for (int i=0; i<nvar; i++){
      if (to_opt[i] == 0){//jump probability
	myJD->jump = var->data[i];
	myJD->set_pstay();
      }
      else if(to_opt[i] == 1){//diffusion constant
	myJD->sigma = var->data[i];
	myJD->get_DiffProp();
      }
      else if(to_opt[i] == 2){//random rate
	myJD->rnd_emit = var->data[i];
      }
      else if(to_opt[i] == 3){//shape parameter
	myJD->myEmit->shape = var->data[i];
	myJD->myEmit->set_EmitProb(myJD->time);
      }
    }
    gsl_vector_free(var);
    gsl_vector_free(range);
  }
  else{
    llh = myJD->get_total_llh();
    //adapt the range if needed
    if ( (opts.mode==3 || opts.mode==4) && (opts.xmin < 0.0 && opts.xmax < 0.0) && opts.bias_fn == NULL){
      myJD->adapt_range();
      printf("Adapted range to %.3e < x < %.3e\n", myJD->myEmit->xmin, myJD->myEmit->xmax);
      llh = myJD->get_total_llh();
    }
  }
  return(llh);
}



void init_parameter(gsl_vector*& var, gsl_vector*& range, vector<int>& to_opt, JumpDiffusion * myJD, cmdl_opts& opts){
  int nvar = (int) var->size;
  for (int i=0; i<nvar; i++){
    if (to_opt[i] == 0){//jump probability
      var->data[i]   = (opts.jump_i > 0.0) ? opts.jump_i : 1.0e-5;
      range->data[i] = 1.0;
    }
    else if(to_opt[i] == 1){//diffusion constant
      if (opts.sigma_i > 0.0){
	var->data[i]   = opts.sigma_i;
      }
      else{
	var->data[i]   = 0.1*myJD->myEmit->dx / sqrt(myJD->myEmit->median_dist);
      }
	range->data[i] = 0.0;
    }
    else if(to_opt[i] == 2){//random rate
      var->data[i]   = (opts.rnd_emit_i > 0.0) ? opts.rnd_emit_i : 1.0e-5;
      range->data[i] = 1.0;
    }
    else if(to_opt[i] == 3){//shape parameter
      var->data[i]   = (opts.shape_i > 0.0) ? opts.shape_i : 100.0;
      range->data[i] = 0.0;
    }
  }
}



double Q( const gsl_vector * x, void * p){
  //JumpDiffusion * myJD = static_cast<JumpDiffusion*> (p);
  fpar * myfpar = static_cast<fpar*> (p);
  int nvar = (int) (myfpar->to_opt).size();
  gsl_vector * var   = gsl_vector_alloc(nvar);
  gsl_vector * range = gsl_vector_alloc(nvar);
  for (int i=0; i<nvar; i++){
    if ((myfpar->to_opt)[i] == 0){//jump probability in [0,1]
      range->data[i] = 1.0;
    }
    else if((myfpar->to_opt)[i] == 1){//diffusion constant in [0,\infty]
      range->data[i] = 0.0;
    }
    else if((myfpar->to_opt)[i] == 2){//random rate in [0,1]
      range->data[i] = 1.0;
    }
    else if((myfpar->to_opt)[i] == 3){//shape parameter in [0,\infty]
      range->data[i] = 0.0;
    }
  }
  gsl_vector ** simplex = NULL;
  gsl_vector * lower    = NULL;
  int err = arg_unmap( x, 0, simplex, lower, var, range);
  if (err==1){
    gsl_vector_free(var);
    gsl_vector_free(range);
    return(1.0e20);
  }
  //set the ML values into the objects
  for (int i=0; i<nvar; i++){
    if ((myfpar->to_opt)[i] == 0){//jump probability
      myfpar->myJD->jump = var->data[i];
      myfpar->myJD->set_pstay();
    }
    else if((myfpar->to_opt)[i] == 1){//diffusion constant
      myfpar->myJD->sigma = var->data[i];
      myfpar->myJD->get_DiffProp();
    }
    else if((myfpar->to_opt)[i] == 2){//random rate
      myfpar->myJD->rnd_emit = var->data[i];
    }
    else if((myfpar->to_opt)[i] == 3){//shape parameter
      myfpar->myJD->myEmit->shape = var->data[i];
      myfpar->myJD->myEmit->set_EmitProb(myfpar->myJD->time);
    }
  }
  // DO FWD TO GET LLH
  double llh = (myfpar->myJD)->get_total_llh();
  gsl_vector_free(var);
  gsl_vector_free(range);
  return(-llh);
}

