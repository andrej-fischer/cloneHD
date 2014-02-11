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
  const char * avcn_fn; 
  const char * chr_fn;
  const char * bulk_fn;
  const char * clones_fn;
  const char * purity_fn;
  const char * cna_jumps_fn;
  const char * baf_jumps_fn;
  const char * snv_jumps_fn;
  const char * maxcn_mask_fn;
  //
  int grid, force, trials, restarts, nmax, seed, maxcn, print_all, learn_priors;
  int mass_gauging;
  double cna_jump, baf_jump, snv_jump;
  double cna_shape, baf_shape, snv_shape;
  double cna_rnd, baf_rnd, snv_rnd, snv_err;
  double baf_pen, snv_pen, snv_fpr;
  double bulk_fix, bulk_sigma, bulk_rnd;
  double min_occ,min_jump;
  int bulk_mean, bulk_prior, bulk_updates;
};


void get_opts( int argc, const char ** argv, cmdl_opts& opts);
void read_opts( const char * opts_fn, cmdl_opts& opts);
void default_opts(cmdl_opts& opts);
void test_opts(cmdl_opts& opts);
void print_opts();
void print_clonal_header( FILE * fp, Clone * myClone,  Emission * myEmit, cmdl_opts& opts);
void print_posterior( FILE * cna_fp, Clone * myClone, Emission * myEmit, int s, cmdl_opts& opts);
void print_phi( FILE * phi_fp, Clone * myClone, Emission * cnaEmit, int s, cmdl_opts& opts);

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
void get_avail_tcn( const char * avcn_fn, Clone * myClone, Emission * myEmit);
void get_purity( const char * purity_fn, gsl_vector *& purity);
void get_fixed_clones( gsl_matrix *& clones, gsl_vector *& mass, const char * clones_fn, int nTimes);
void get_jump_probability(  Clone * myClone, cmdl_opts& opts);
void get_bias_field( Clone * myClone, cmdl_opts& opts);

void print_llh_for_set(gsl_matrix * clones, gsl_vector * mass, Clone * myClone, cmdl_opts& opts);

int infer_clones( gsl_matrix * Clones, gsl_vector * Mass, Clone * myClone, cmdl_opts& opts);

double get_clones( gsl_matrix *& clones, 
		   gsl_matrix *& Clones, 
		   gsl_vector *& mass, 
		   gsl_vector *& Mass, 
		   gsl_matrix *& priors, 
		   Clone * myClone,
		   cmdl_opts& opts,
		   double& cl,
		   double& bl,
		   double& sl);

double get_clones_cna( gsl_matrix *& clones, 
		       gsl_matrix *& Clones, 
		       gsl_vector *& mass, 
		       gsl_vector *& Mass, 
		       Clone * myClone,
		       cmdl_opts& opts,
		       double& cl,
		       double& bl,
		       double& sl);

double get_clones_baf( gsl_matrix *& clones, 
		       gsl_matrix *& Clones, 
		       Clone * myClone,
		       cmdl_opts& opts
		       );

double get_clones_snv_ncorr( gsl_matrix *& clones, 
			     gsl_matrix *& Clones, 
			     gsl_matrix *& priors, 
			     Clone * myClone,
			     cmdl_opts& opts
			     );

double get_clones_snv_wcorr( gsl_matrix *& clones, 
			     gsl_matrix *& Clones, 
			     Clone * myClone,
			     cmdl_opts& opts
			     );

double cna_only_mass_noclones( gsl_vector *& mass, Clone * myClone, int restarts, int& steps);
//double cna_only_clones_mass( gsl_matrix*& clones, gsl_vector*& mass, Clone * myClone, int& steps);
double cna_clones_fixed_mass( gsl_matrix*& clones, Clone * myClone, int restarts, 
			      int& steps, double& cl, double& bl, double& sl);
double cna_mass_fixed_clones(gsl_vector*& mass, Clone * myClone, int restarts, 
			     int& steps, double& cl, double& bl, double& sl);

double cna_clones_mass( gsl_matrix*& clones, gsl_vector*& mass, Clone * myClone, int restarts, 
			int& steps, double& cl, double& bl, double& sl);


void get_candidate_masses( gsl_matrix * clones, 
			   gsl_vector * mass, 
			   Clone * myClone, 
			   gsl_matrix*& candidate_masses,
			   gsl_vector*& levels,
			   double min_occ);

double cna_llh_all_fixed(Clone * myClone);

double baf_clones( gsl_matrix*& clones, Clone* myClone, int restarts, int& steps);

double snv_clones_fixed_priors( gsl_matrix*& clones, Clone * myClone, int restarts, int& steps);
void snv_iterative_bulk_update( double& llh, gsl_matrix*& clones, Clone * myClone, int iter);
double snv_priors_fixed_clones( gsl_matrix*& priors, Clone * myClone, int restarts, int& steps);
double snv_clones_priors( gsl_matrix*& clones, gsl_matrix*& priors, Clone * myClone, int restarts, int& steps);
void snv_bulk_update(Clone * myClone);

void set_random_start_freq(gsl_vector *& freq, double lower);
void report_results( double cl, double bl, double sl, int steps, gsl_vector * mass, gsl_matrix * freq);

double Q( const gsl_vector * x, void * p);
struct Q_par{
  Clone * myClone;
  int nSimplex;
  vector<int> simplexD;
  int clones_fixed;
  int mass_fixed;
  int prior_fixed;
  int cna,baf,snv;
};


