//clone.h

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
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

//#include <unordered_map>

// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"

//sort by values in array
struct SortAsc{
  double * arg; 
  bool operator() (int i, int j) {return(arg[i]<arg[j]);}
};
struct SortDesc{
  double * arg; 
  bool operator() (int i, int j) {return(arg[i]>arg[j]);}
};


using namespace std;

//forward class declaration
class Emission;

class Clone{
 public:
  Clone();
  ~Clone();
  Emission * cnaEmit, * bafEmit, * snvEmit;
  // functions
  // *** in clone.cpp ********************************************************************************
  void allocate(Emission * cnaEmit, Emission * bafEmit, Emission * snvEmit, const char * chr_fn);
  void clean();
  int nTimes, nClones;
  int allocated, is_set;
  int total_loci;
  std::set<int> chrs;
  int maxChr;
  // copy number combination states and frequencies
  int ** copynumber;
  void set_copynumbers();
  gsl_matrix * freqs;
  double * purity;
  gsl_vector * min_purity;
  double ** clone_spectrum;
  int nLevels, nFreq;
  void set(const gsl_matrix * freq);
  void set_clone_spectrum(const gsl_matrix * freq);
  double logzero;
  //normal copy number per chr
  std::map<int,int> normal_copy;
  void get_normal_copy(const char * chr_fn);
  void set_normal_copy(const char * chr_fn);
  int maj_ncn;
  // maximum total copy number per chr
  int maxtcn;
  std::map<int, vector<int> > maxtcn_input;
  std::map<int, vector<int> > maxtcn_per_clone;
  std::set<int> all_maxtcn;
  void set_maxtcn_per_clone();
  void set_all_levels();
  std::map<int,int> level_of;
  double *** tcn;
  double *** log_tcn;
  void allocate_tcn();
  void set_tcn();
  // penalties and parameters
  double cna_pen_zero, cna_pen_norm, cna_pen_diff;
  double baf_pen_comp;
  double snv_pen_high, snv_pen_mult;
  double snv_fpr,snv_fpf;
  // pre computed variables (consider unordered_map<>)
  map< unsigned int, double> logn;
  map< unsigned int, double> loggma;
  int logn_set;
  void set_logn();
  // masses
  void set_mass(gsl_vector * mass);
  gsl_vector * mass, * log_mass;
  gsl_vector * nmean;
  void get_nmean();
  // mass-gauging
  gsl_matrix * margin_map;
  void get_cna_marginals();
  gsl_matrix * marginals;
  gsl_matrix * mass_candidates;
  void get_mass_candidates();
  vector<int> levels_sorted;
  gsl_matrix * copynumber_post;
  gsl_vector * majcn_post;
  gsl_matrix ** bafSymMap;
  void set_bafSymMap();
  // *** in clone-prior.cpp **************************************************************************
  int learn_priors;
  gsl_matrix * baf_prior_map;
  gsl_matrix ** snv_prior_from_cna_baf_map;
  gsl_matrix *  snv_prior_from_cna_map;
  void set_margin_map();
  void set_baf_prior_map();
  void set_snv_prior_map();
  void get_baf_prior_from_cna_post(gsl_vector * prior, gsl_vector * post);
  void get_snv_prior_from_cna_post(gsl_vector * prior, gsl_vector * post);
  void get_snv_prior_from_cna_baf_post(gsl_vector * prior, gsl_vector * cnapost, gsl_vector * bafpost);
  void apply_snv_prpc( gsl_vector * prior, gsl_matrix * snv_prpc, double pc0);
  std::map<int,gsl_vector*> snv_prior;
  void set_cna_prior( gsl_vector * prior, int sample);
  void set_snv_prior( gsl_matrix * prior_param );
  gsl_matrix * initial_snv_prior_param;
  void initialize_snv_prior_param();
  // mean total c.n. and available c.n.
  void get_mean_tcn(int sample);//for cna only
  void map_mean_tcn( Emission * fromEmit, int from_sample,  Emission * toEmit);//from cna/baf to baf/snv
  void get_avail_cn( Emission * myEmit, int sample);//for cna/baf
  void get_snv_prior_from_av_cn( gsl_vector * prior, int sample, int evt);
  double *** cn_usage;
  void allocate_cn_usage();
  void set_cn_usage();
  // *** in clone-predict.cpp ************************************************************************  
  gsl_matrix ** TransMat_cna;
  gsl_matrix ** TransMat_snv;
  void set_TransMat_cna();
  void set_TransMat_snv();
  void set_TransMat_cna(gsl_matrix * Trans, int chr);
  void set_TransMat_snv(gsl_matrix * Trans, int chr);
  void predict( gsl_vector * prior, gsl_vector * post, Emission * myEmit, double pj, gsl_matrix * T);
  void predict( gsl_vector * prior, gsl_vector * post, Emission * myEmit, double pj, gsl_vector * flat);
  void apply_maxtcn_mask( gsl_vector * prior, int chr, int log_space);
  // *** in clone-fwd-bwd.cpp ************************************************************************
  void combine_prior(gsl_vector*& prior, gsl_vector*& mem, int n);
  void scale_prior(gsl_vector*& prior, int n);
  void do_cna_Fwd( int sample, double& llh, double*& llhs);
  void do_cna_Bwd( int sample, double& ent);
  void do_baf_Fwd( int sample, double& llh, double*& llhs);
  void do_baf_Bwd( int sample, double& ent);
  void do_snv_Fwd( int sample, double& llh, double*& llhs);
  void do_snv_Bwd( int sample, double& ent);
  int got_gamma, save_cna_alpha, save_baf_alpha, save_snv_alpha;
  gsl_matrix ** alpha_cna, ** alpha_baf, ** alpha_snv;
  gsl_matrix ** gamma_cna, ** gamma_baf, ** gamma_snv;
  double entropy(gsl_vector * x);
  int get_gofs;
  void allocate_all_gofs();
  void get_cna_gof(gsl_vector * post, int sample, int evt);
  void get_baf_gof(gsl_vector * post, int sample, int evt);
  void get_snv_gof(gsl_vector * post, int sample, int evt);
  double * cna_gofs, * baf_gofs, * snv_gofs;
  double *** cna_all_gofs, *** baf_all_gofs, *** snv_all_gofs;
  //void sym_baf( gsl_vector * bafPost, gsl_vector * cnvPost);
  //gsl_matrix ** map1, ** map2;
  //int symmetrize_baf;
  // *** in clone-llh.cpp ****************************************************************************
  double get_cna_posterior(int sample);
  double get_baf_posterior(int sample);
  double get_snv_posterior(int sample);
  // log-likelihoods
  double total_llh, total_entropy; 
  double cna_total_llh, baf_total_llh, snv_total_llh;
  double * cna_llhs, * baf_llhs, * snv_llhs;
  double cna_total_ent, baf_total_ent, snv_total_ent;
  double get_all_total_llh();
  double get_cna_total_llh();
  double get_baf_total_llh();
  double get_snv_total_llh();
  // *** in clone-update.cpp *************************************************************************
  // update step CNA
  double update( gsl_vector * prior, gsl_vector * post, Emission * myEmit, int sample, int site, double*& llhs);
  void update_cna( gsl_vector * prior, gsl_vector * post, int sample, int site, gsl_matrix * Post);
  void update_cna_event( gsl_vector * prior, gsl_vector * post, int sample, int evt, gsl_matrix * Post);
  void update_cna_site_noclone( gsl_vector * post, int sample, int site, gsl_matrix * Post);
  void update_cna_site_wclone( gsl_vector * prior, gsl_vector * post, int sample, int site, gsl_matrix * Post);
  // update step BAF
  void update_baf( gsl_vector * prior, gsl_vector * post, int sample, int evt, gsl_matrix * Post);
  void update_baf_event( gsl_vector * prior, gsl_vector * post, int sample, int evt, gsl_matrix * Post);
  void update_baf_site( gsl_vector * prior, gsl_vector * post, int sample, int site, gsl_matrix * Post);
  // update step SNV
  void update_snv( gsl_vector * prior, gsl_vector * post, int sample, int evt, gsl_matrix * Post);
  void update_snv_event( gsl_vector * prior, gsl_vector * post, int sample, int evt, gsl_matrix * Post);
  void update_snv_site_ncorr( gsl_vector * prior, gsl_vector * post, int sample, int site, gsl_matrix * Post);
  void update_snv_site_fixed( gsl_vector * prior, gsl_vector * post, int sample, int site, gsl_matrix * Post);
  void update_snv_site_nfixed( gsl_vector * prior, gsl_vector * post, int sample, int site, gsl_matrix * Post);	
  // 1D and 2D interpolation
  double get_interpolation(double x, double xmin, double xmax, double dx, gsl_vector * emit);
  double get_interpolation(double x, double xmin, double xmax,
			   double y, double ymin, double ymax,
			   gsl_matrix * emit);
  double trapezoidal( gsl_vector * blk, double a, double b, gsl_vector * emit, int get_log);
  //precomputed log-emission probabilities for BAF update
  gsl_matrix *** bafEmitLog;
  gsl_matrix *** cnaEmitLog;
  gsl_matrix **** snvEmitLog;
  void get_bafEmitLog();
  void get_cnaEmitLog();
  void get_snvEmitLog();
  double * cna_xmin;
  double * cna_xmax;
  int snvGrid, bafGrid, cnaGrid, bulkGrid;
  // BIC complecity penalty
  double complexity;
  void get_complexity();
  // *** in clone-bulk.cpp ***************************************************************************
  double bulk_fix;
  void allocate_bulk_dist();
  void allocate_bulk_mean();
  void update_bulk( int time, int sample);
  gsl_matrix **  bulk_prior;
  gsl_matrix *** bulk_post;
  gsl_matrix *** bulk_dist;
  double **  bulk_prior_mean;
  double *** bulk_post_mean;
  double *** bulk_mean;
  double ***  bulk_min;
  void set_bulk_to_prior();
  void set_bulk_to_post();
  void get_bulk_min();
  //SNV bulk updates
  void update_bulk(int sample);
  void get_bulk_post_dist( gsl_vector * bprior, gsl_vector * bpost, gsl_vector * emit, int time, int sample, int idx);
};
