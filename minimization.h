//minimization.h

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
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_linalg.h"

using namespace std;


// The  numerical minimization routine to get the clone frequencies...
double find_local_optimum( int nSimplex,
			   gsl_vector**& simplex, 
			   gsl_vector * lower,
			   gsl_vector * other,
			   gsl_vector * range,
			   void * params,
			   double (*obj_fn)( const gsl_vector * x, void * p),
			   double prec,
			   int& steps,
			   int verbose
			   );

double find_optimum_wrestarts(int nSimplex,
			      gsl_vector**& simplex,
			      gsl_vector * lower,
			      gsl_vector * other,
			      gsl_vector * range,
			      void * params,
			      double (*obj_fn)( const gsl_vector * x, void * p),
			      double prec,
			      int restarts,
			      int& steps,
			      int verbose
			      );

void spherical_random_step_uniform( double ri, double& rf, double lower,
				    const gsl_vector * anglei, gsl_vector*& anglef,
				    double eps);
void simplex_random_step_uniform(const gsl_vector*simplexi, gsl_vector*& simplexf,
				 double lower, double eps);

void arg_map( int nSimplex, gsl_vector**& simplex, gsl_vector * lower, const gsl_vector * other, const gsl_vector * range, gsl_vector ** x);
int arg_unmap( const gsl_vector * x, int nSimplex, gsl_vector**& simplex, gsl_vector * lower, gsl_vector * other, const gsl_vector * range);

double spherical_to_simplex( double radial, const gsl_vector * angle, gsl_vector *& simplex, int getLJD);

void simplex_to_spherical( const gsl_vector * simplex, double& radial, gsl_vector*& angle);

double logify( double x,  double R);
double delogify(double y, double R);

/*
double simulated_annealing(
			   gsl_matrix * freqs, // set of points in or on simplex
			   gsl_vector * other, // other arguments
			   gsl_vector * range,// range of other arguments
			   void * params,
			   double (*obj_fn)( const gsl_vector * x, void * p),
			   int& steps
			   );
int accepted(double dE, double T);
*/
