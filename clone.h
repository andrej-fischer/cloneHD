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
class Bulk;

class Clone{
 public:
  Clone();
  ~Clone();
  Emission * cnaEmit, * bafEmit, * snvEmit;
  Bulk * myBulk;
  void allocate(Emission * cnaEmit, Emission * bafEmit, Emission * snvEmit, const char * chr_fn);
  void get_normal_copy(const char * chr_fn);
  void set_normal_copy(const char * chr_fn);
  void clean();
  int nTimes, nClones, maxcn;
  int allocated, is_set;
  int maj_ncn;
  double bulk_fix,snv_err,baf_pen,snv_pen, snv_fpr;
  int total_loci;
  //SNV bulk contribution
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
  //
  gsl_matrix * margin_map;
  gsl_matrix * baf_prior_map;
  gsl_matrix ** snv_prior_map;
  void set_margin_map();
  void set_baf_prior_map();
  void set_snv_prior_map();
  void get_cna_marginals();
  gsl_matrix * marginals;
  gsl_matrix * mass_candidates;
  void get_mass_candidates();
  vector<int> levels_sorted;
  void get_baf_prior_from_cna_post(gsl_vector * prior, gsl_vector * post);
  void get_snv_prior_from_cna_post(gsl_vector * prior, gsl_vector * post);
  void get_snv_prior_from_cna_baf_post(gsl_vector * prior, gsl_vector * cnapost, gsl_vector * bafpost);
  void apply_snv_prpc( gsl_vector * prior, gsl_matrix * snv_prpc);
  //
  int ** copynumber;
  void set_copynumbers();
  int * normal_copy; //copy number of normal human DNA
  std::map<int,gsl_vector*> cn_prior_snv;
  std::map<int,gsl_vector*> cn_prior_baf;
  void set_cn_prior_cna( gsl_vector * cn_prior, int sample);
  void set_cn_prior_snv( gsl_matrix * prior_per_clone);
  //void set_cn_prior_baf();
  gsl_matrix * init_cn_prior_snv;
  void initialize_cn_prior_snv();
  gsl_matrix * copynumber_post;
  gsl_vector * majcn_post;
  //
  gsl_matrix * TransMat_snv, * TransMat_cna;
  void set_TransMat_cna();
  void set_TransMat_snv();
  //
  void set_mass(gsl_vector * mass);
  gsl_vector * mass, * log_mass;
  gsl_vector * nmean;
  //
  //unordered_map< unsigned int, double> logn;
  //unordered_map< unsigned int, double> loggma;
  map< unsigned int, double> logn;
  map< unsigned int, double> loggma;
  int logn_set;
  void set_logn();
  unsigned int max_nFreq;
  void get_nmean();
  //
  gsl_matrix * freqs;
  double * purity;
  gsl_vector * min_purity;
  double ** clone_spectrum;
  int nLevels, nFreq;
  void set(const gsl_matrix * freq);
  void set_clone_spectrum(const gsl_matrix * freq);
  //
  void get_event_map(Emission * myEmit);
  void get_phi(int sample);//for cna
  void map_phi( Emission * fromEmit, int from_sample,  Emission * toEmit);//from cna/baf to baf/snv
  int got_gamma, save_cna_alpha,save_baf_alpha,save_snv_alpha;
  //
  void do_cna_Fwd( int sample, double& llh);
  void do_cna_Bwd( int sample, double& ent);
  void do_baf_Fwd( int sample, double& llh);
  void do_baf_Bwd( int sample, double& ent);
  void do_snv_Fwd( int sample, double& llh);
  void do_snv_Bwd( int sample, double& ent);
  //
  gsl_matrix ** alpha_cna, ** alpha_baf, ** alpha_snv;
  gsl_matrix ** gamma_cna, ** gamma_baf, ** gamma_snv;
  double total_llh, total_entropy; 
  double cna_total_llh, baf_total_llh, snv_total_llh;
  double get_all_total_llh();
  double get_cna_total_llh();
  double get_baf_total_llh();
  double get_snv_total_llh();
  double get_cna_posterior(  int sample);
  double get_baf_posterior(  int sample);
  double get_snv_posterior(  int sample);
  double entropy(gsl_vector * x);
  //
  void predict( gsl_vector * prior, gsl_vector * post, Emission * myEmit, double pj, gsl_matrix * T);
  void predict( gsl_vector * prior, gsl_vector * post, Emission * myEmit, double pj, double flat);
  double update( gsl_vector * prior, gsl_vector * post, Emission * myEmit, int sample, int site);
  //
  void update_cna( gsl_vector * post, int sample, int site);
  void update_cna_event( gsl_vector * post, int sample, int evt);
  void update_cna_site_noclone( gsl_vector * post, int sample, int site);
  void update_cna_site_wclone( gsl_vector * post, int sample, int site);
  //
  void update_baf( gsl_vector * post, int sample, int evt);
  void update_baf_event( gsl_vector * post, int sample, int evt);
  void update_baf_site( gsl_vector * post, int sample, int site);
  //
  void update_snv( gsl_vector * prior, gsl_vector * post, int sample, int site);
  void update_snv_event( gsl_vector * post,int sample,int evt);
  void update_snv_fixed(  gsl_vector * prior, gsl_vector * post, int sample, int site);
  void update_snv_nfixed( gsl_vector * prior, gsl_vector * post, int sample, int site);	  
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
  //precomputed variables for cna update 
  double *** tcn;
  double *** log_tcn;
  void set_tcn(int sample);
  double complexity;
  void get_complexity();
  //SNV bulk updates
  void update_bulk(int sample);
  void get_bulk_post_dist( gsl_vector * bprior, gsl_vector * bpost, gsl_vector * emit, int time, int sample, int idx);
  
};



