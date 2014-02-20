//common-functions.h

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
class Emission;

using namespace std;

void get_dims( const char * data_fn, int& nTimes, vector<int>& chrs, vector<int>& nSites, int keep);
void get_data( const char * data_fn, Emission * myEmit);
void get_bias( const char * bias_fn, Emission * myEmit);
double get_mean( gsl_vector * dist, double xmin, double xmax);
double get_var(  gsl_vector * dist, double xmin, double xmax, double mean);
