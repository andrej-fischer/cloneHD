//cloneHD-functions.h

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

//own headers
#include "clone.h"
#include "minimization.h"
#include "emission.h"

using namespace std;


struct cmdl_opts{
  const char * cna_fn;
  const char * baf_fn;
  const char * snv_fn;
  const char * pre;
  const char * bias_fn;
  const char * mntcn_fn;
  const char * maxtcn_fn;
  const char * avcn_fn; 
  const char * chr_fn;
  const char * bulk_fn;
  const char * clones_fn;
  const char * purity_fn;
  const char * cna_jumps_fn;
  const char * baf_jumps_fn;
  const char * snv_jumps_fn;
  //
  int force, trials, restarts, nmax, seed, maxtcn, print_all, learn_priors;
  int mass_gauging;
  double cna_jump, baf_jump, snv_jump;
  double cna_shape, baf_shape, snv_shape;
  double cna_rnd, baf_rnd, snv_rnd;
  double baf_pen, snv_pen;
  double snv_fpr, snv_fpf;
  double bulk_fix, bulk_sigma, bulk_rnd;
  double min_occ,min_jump;
  int bulk_mean, bulk_prior, bulk_updates;
  int cnaGrid, bafGrid, snvGrid, bulkGrid;
};


void get_opts( int argc, const char ** argv, cmdl_opts& opts);
void read_opts( const char * opts_fn, cmdl_opts& opts);
void default_opts(cmdl_opts& opts);
void test_opts(cmdl_opts& opts);
void print_opts();
void print_clonal_header( FILE * fp,  Clone * myClone, Emission * myEmit, cmdl_opts& opts);
void print_posterior( FILE * cna_fp,  Clone * myClone, Emission * myEmit,  int s, cmdl_opts& opts);
void print_mean_tcn( FILE * mntcn_fp, Clone * myClone, Emission * cnaEmit, int s, cmdl_opts& opts);
void print_avail_cn( FILE * avcn_fp,  Clone * myClone, Emission * cnaEmit, int s, cmdl_opts& opts);

void get_dims( const char * data_fn,
	       int& nTimes,
	       vector<int>& chrs,
	       vector<int>& nSites
	       );
void get_data( const char * data_fn, 
	       Emission * myEmit
	       );
void get_bulk_prior( gsl_matrix **& bulk, 
		     double **& bulk_prior_mean,
		     Emission * myEmit,
		     cmdl_opts& opts
		     );
void get_track(const char * track_fn,
	       gsl_matrix **& distribution, 
	       double **& mean, 
	       double **& var, 
	       Emission * myEmit
	       );

void get_maxcn_mask(const char * maxcn_mask_fn, Clone * myClone, int maxcn_gw);
void get_mean_tcn( const char * mtcn_fn, Clone * myClone, Emission * myEmit);
void get_avail_cn( const char * avcn_fn, Clone * myClone, Emission * myEmit);
void get_purity( const char * purity_fn, gsl_vector *& purity);
void get_fixed_clones( gsl_matrix *& clones, gsl_vector *& mass, const char * clones_fn, int nTimes);
void get_jump_probability(  Clone * myClone, cmdl_opts& opts);
void get_bias_field( Clone * myClone, cmdl_opts& opts);
void print_llh_for_set(gsl_matrix * clones, gsl_vector * mass, Clone * myClone, cmdl_opts& opts);
