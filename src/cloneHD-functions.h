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
class Clone;
class Emission;

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
  double cna_pen_zero, cna_pen_norm, cna_pen_diff;
  double baf_pen_comp;
  double snv_pen_high, snv_pen_mult;
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
void print_usage();

void print_all_results( Clone * myClone, cmdl_opts& opts);
void print_posterior_header( FILE * fp,  Clone * myClone, Emission * myEmit, cmdl_opts& opts);
void print_posterior( FILE * fp,  Clone * myClone, Emission * myEmit,  int s, cmdl_opts& opts);
void print_perclone_header( FILE * fp,  Clone * myClone, Emission * myEmit, cmdl_opts& opts);
void print_perclone_posterior( FILE ** fp,  Clone * myClone, Emission * myEmit,  int s, cmdl_opts& opts);
void print_mean_tcn( FILE * mntcn_fp, Clone * myClone, Emission * cnaEmit, int s, cmdl_opts& opts);
void print_avail_cn( FILE * avcn_fp,  Clone * myClone, Emission * cnaEmit, int s, cmdl_opts& opts);
void print_gof( Clone * myClone, Emission * myEmit, cmdl_opts& opts);

void get_cna_data(Emission * cnaEmit, cmdl_opts& opts, int& nTimes);
void get_baf_data(Emission * bafEmit, cmdl_opts& opts, int& nTimes, int& nT);
void get_snv_data(Emission * snvEmit, cmdl_opts& opts, int& nTimes, int& nT);
void get_snv_bulk_prior( Clone * myClone, cmdl_opts& opts);
void get_track(const char * fn, gsl_matrix **& dist, double **& mn, double **& var, Emission * myEmit);
void match_jumps(const char * jumps_fn, Emission * myEmit);
void get_maxtcn_input(const char * maxtcn_fn, int maxtcn_gw, Clone * myClone);
void get_mean_tcn( const char * mtcn_fn, Clone * myClone, Emission * myEmit);
void get_avail_cn( const char * avcn_fn, Clone * myClone, Emission * myEmit);
void get_purity( const char * purity_fn, gsl_vector *& purity);
void get_fixed_clones( gsl_matrix *& clones, gsl_vector *& mass, const char * clones_fn, int nTimes);
void get_jump_probability(  Clone * myClone, cmdl_opts& opts);
void get_bias_field( Clone * myClone, cmdl_opts& opts);
void print_llh_for_set(gsl_matrix * clones, gsl_vector * mass, Clone * myClone, cmdl_opts& opts);
