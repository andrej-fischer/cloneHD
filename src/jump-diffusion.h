//jump-diffusion.h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_multimin.h"

using namespace std;

class JumpDiffusion{
public:
  JumpDiffusion(Emission * emit, int time);
  ~JumpDiffusion();
  Emission * myEmit; 
  int nSamples;
  int time;
  int mode;
  double sigma, jump, rnd_emit;
  int Fwd_done, Bwd_done, wTotal, save_alpha;
  int gridSize;
  int * nSites;
  unsigned int ** dist;
  unsigned int ** loci;
  unsigned int ** mask;  
  double ** pstay;
  double ** pjump;
  double ** pnojump;
  double ** bias;
  void get_EmitProb(int read, int depth, double * xgrid, gsl_vector * eprob);// emission probability
  gsl_vector * proposal;
  gsl_matrix ** alpha;
  gsl_matrix ** gamma;
  gsl_matrix ** total;
  void set_pstay();
  int pstay_set;
  double do_Fwd(int sample);
  void do_Bwd(int sample);
  //
  gsl_matrix ** DiffProp;
  int set_DiffProp(gsl_matrix * propagator, double variance);
  void get_DiffProp();
  int DiffProp_set;
  void reset_DiffProp();
  vector<int> is_identity;
  map<unsigned int,int> position;
  //
  //void set_DiffProp_Log(gsl_matrix * propagator, double variance);
  int predict(gsl_vector * prior, gsl_vector * post, gsl_matrix*& DiffProp, gsl_matrix**& DP_pt, int sampe, int site);
  double total_llh;
  double get_total_llh();
  void get_posterior(int sample);
  int adapt_range();
};

