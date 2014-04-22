//emission.h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <set>
//#include <unordered_map>
#include <vector>
#include <list>
#include <algorithm>

// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_psi.h"
#include "gsl/gsl_statistics_double.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_cdf.h"



using namespace std;

class Emission{
public:
  Emission();
  void set(int ntimes, vector<int>& chrs, vector<int>& nsites, int grid);
  ~Emission();
  void clear();
  void delete_old_Emit();
  int is_set;
  void set_dist();
  int dist_set;
  int connect;
  double median_dist;
  map<unsigned int, int> dist_count;
  map<unsigned int, int> frequent_dist;
  int get_log, get_der, get_mv;
  unsigned int nmax, Nmax;
  //unordered_map< unsigned int, unordered_map< unsigned int, gsl_vector*> > EmitProb;
  //unordered_map< unsigned int, unordered_map< unsigned int, gsl_vector*> > EmitLog;
  map< unsigned int, map< unsigned int, gsl_vector*> > EmitProb;
  map< unsigned int, map< unsigned int, gsl_vector*> > EmitLog;
  double shape, log_shape, rnd_emit;
  double minRate, maxRate;
  //
  int mode, reflect, log_space;
  void set_EmitProb(int time);
  void binomial(int N, int n);
  void beta_binomial(int N, int n);
  void poisson(int N, int n);
  void negative_binomial(int N, int n);
  double get_single_EmitLog(double x, unsigned int n, unsigned int N);
  void get_eprob_wBias( gsl_vector * eprob, gsl_vector * emit, double b, unsigned int n, unsigned int N, int get_log);
  int EmitProb_set;
  //
  int nTimes, nSamples;
  int gridSize;
  double dx,xmin,xmax;
  double dy,ymin,ymax;
  double * xgrid;
  double * ygrid;
  unsigned int *** reads;
  unsigned int *** depths;
  unsigned int ** loci;
  unsigned int ** mask;
  unsigned int ** dist;
  unsigned int *** nObs;
  void get_nObs();
  int * nSites;
  int * chr;
  std::set<int> chrs;
  int * idx_of;
  int maxchr;
  double ** bias;
  double ** log_bias;
  void allocate_bias();
  void allocate_mean_tcn();
  void allocate_av_cn(int maxcn);
  int total_loci, total_events;
  unsigned int total_dist;
  void set_grid();
  void init_range(int time);
  int range_set;
  void reset_mask();
  double get_pval(int time, int sample, int site, double mean);
  void coarse_grain_jumps( int sample, double plow, int range);
  double ** pjump;
  void set_pjump(double jump);
  double *** mean_tcn;//mean total copy number
  double **** av_cn;//copy number availability
  void init_events();
  int * nEvents;
  unsigned int ** Event_of_idx;// map from idx to cnv-event 
  void map_idx_to_Event(Emission * Emit, int sample);
  void map_jumps(Emission * Emit);
  void add_break_points_via_jumps(Emission * Emit, double pmin);
  void get_events_via_jumps();
  unsigned int ** idx_of_event;// map from event to idx
  unsigned int ** event_of_idx;// map from idx to event
  int idx_to_Event_mapped;
  int coarse_grained;
};


bool value_comparer(std::map<int,double>::value_type &i1, std::map<int,double>::value_type &i2);
