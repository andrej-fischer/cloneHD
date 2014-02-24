//cloneHD-inference.h

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

class Clone;
class Emission;

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
