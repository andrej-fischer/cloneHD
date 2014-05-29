//clone-update.cpp

//own headers...
#include "emission.h"
#include "log-space.h"
#include "clone.h"

using namespace std;

// general update step...
double Clone::update( gsl_vector * prior, gsl_vector * post, Emission * myEmit, int sample, int evt, double*& llhs){ 
  gsl_matrix * Post = NULL;
  gsl_vector * mem  = NULL;
  if (llhs != NULL){
    Post = gsl_matrix_alloc( nTimes, nLevels);
    mem  = gsl_vector_alloc( nLevels);
  }
  if (myEmit == cnaEmit){
    Clone::update_cna( prior, post, sample, evt, Post);
  }
  else if (myEmit == bafEmit){
    Clone::update_baf( prior, post, sample, evt, Post);
  }
  else if (myEmit == snvEmit){
    Clone::update_snv( prior, post, sample, evt, Post);
  }
  double norm=0;
  // update step: multiply prior with the likelihood and scale to get the posterior...
  if (myEmit->log_space){
    if (Post!=NULL){
      for (int t=0; t<nTimes; t++){
	gsl_matrix_get_row(mem,Post,t);
	gsl_vector_add( mem, prior);
	norm = log_vector_norm( mem );
	if ( norm != norm ) abort();
	llhs[t] += norm;
      }
      gsl_matrix_free(Post);
      gsl_vector_free(mem);
    }
    gsl_vector_add( post, prior);
    norm = log_vector_norm( post );
    if ( norm != norm ) abort();
    gsl_vector_add_constant( post, -norm);
    return(norm);
  }
  else{
    if (gsl_vector_isnull(post)){      
      if (Post != NULL){
	for (int t=0; t<nTimes; t++) llhs[t] += logzero;
	gsl_matrix_free(Post);
	gsl_vector_free(mem);
      }
      gsl_vector_memcpy( post, prior);
      return(logzero);
    }
    else{
      if (Post!=NULL){
	for (int t=0; t<nTimes; t++){
	  gsl_matrix_get_row(mem,Post,t);
	  gsl_vector_mul( mem, prior);
	  norm = gsl_blas_dasum(mem);
	  if ( norm != norm || norm <= 0.0) abort();
	  llhs[t] += log(norm);
	}
	gsl_matrix_free(Post);
	gsl_vector_free(mem);
      }
      gsl_vector_mul( post, prior);
      norm = gsl_blas_dasum(post);
      if ( norm != norm || norm <= 0.0){
	cout<<"ERROR\n";
	abort();
      }
      gsl_vector_scale( post, 1.0 / norm);
      return( log(norm) );
    }
  }
}


//*** CNA UPDATE ************************************************************************
void Clone::update_cna( gsl_vector * prior, gsl_vector * post, int sample, int evt, gsl_matrix * Post){ 
  if ( cnaEmit->coarse_grained ){
    Clone::update_cna_event( prior, post, sample, evt, Post);
  }
  else{
    if (nClones==0){
      Clone::update_cna_site_noclone( post, sample, evt, Post);
    }
    else{
      Clone::update_cna_site_wclone( prior, post, sample, evt, Post);
    }
  }
}

void Clone::update_cna_event( gsl_vector * prior, gsl_vector * post, int sample, int evt, gsl_matrix * Post){
  gsl_vector_set_all( post, 0.0);//log-space!
  if (Post != NULL) gsl_matrix_set_all(Post,0.0);
  gsl_vector_view emit;
  double x,dx,val;
  int chr = cnaEmit->chr[sample];
  for (int time=0; time<nTimes; time++){
    emit = gsl_matrix_row( cnaEmitLog[time][sample], evt);
    dx   = (cna_xmax[time] - cna_xmin[time]) / double(cnaGrid);
    for (int level=0; level<nLevels; level++){
      if (prior->data[level] <= logzero) continue;
      x = tcn[chr][time][level] * mass->data[time];
      if ( x < cna_xmin[time] || x > cna_xmax[time]){//outside range...
	val = logzero;
      }
      else{//inside range...
	val = Clone::get_interpolation( x, cna_xmin[time], cna_xmax[time], dx, &emit.vector);
      }
      post->data[level] += val;
      if (Post!=NULL) gsl_matrix_set(Post,time,level,val);
    }
  }
}

void Clone::update_cna_site_noclone( gsl_vector * post, int sample, int site, gsl_matrix * Post){//no-clone special case...
  gsl_vector_set_all( post, 1.0);
  if (Post != NULL) gsl_matrix_set_all(Post, 1.0);
  double val=0,rate=0;
  int ncn = normal_copy[cnaEmit->chr[sample]];
  double b  =  (cnaEmit->bias != NULL) ? cnaEmit->bias[sample][site] : -1.0;
  double lb =  (cnaEmit->log_bias != NULL) ? cnaEmit->log_bias[sample][site] : 0.0;
  for (int time=0; time<nTimes; time++){
    int N = cnaEmit->depths[time][sample][site];
    int n = cnaEmit->reads[time][sample][site];   
    double rnd = cnaEmit->rnd_emit / (double(N)*(cnaEmit->maxRate - cnaEmit->minRate) + 1.0);
    double nrnd = (1.0 - cnaEmit->rnd_emit);
    if (cnaEmit->mode == 3){//Poisson
      val  = - mass->data[time] * double(N*ncn);
      if (b>0.0) val *= b;
      val -= loggma[n+1];
      rate = log_mass->data[time] + logn[N] + logn[ncn];
      if (b>0.0) rate += lb;
      val += double(n) * rate;    
    }
    else if (cnaEmit->mode == 4){//Negative Binomial
      double n1   = double(N)*cnaEmit->shape;//rounding (error?)
      double nrml = double(ncn);
      val  = loggma[int(n1)+n] - loggma[int(n1)] - loggma[n+1] + n1*cnaEmit->log_shape;
      rate = mass->data[time] * nrml;
      if (b>0.0) rate *= b;
      val += double(n) * log(rate);
      val -= (double(n) + n1) * log(rate + cnaEmit->shape);
      if (val != val) abort();
    }
    val = exp(val);
    if (rnd > 0.0) val = val*nrnd + rnd;
    post->data[0] *= val; 
    if (Post!=NULL) gsl_matrix_set(Post,time,0,val);
  }
}

void Clone::update_cna_site_wclone( gsl_vector * prior, gsl_vector * post, int sample, int site, gsl_matrix * Post){
  gsl_vector_set_all( post, 1.0);//NOT logspace!
  if (Post != NULL) gsl_matrix_set_all(Post, 1.0);
  double b  = (cnaEmit->bias != NULL) ? cnaEmit->bias[sample][site] : -1.0;
  double lb = (cnaEmit->log_bias != NULL) ? cnaEmit->log_bias[sample][site] : 0.0;
  int chr = cnaEmit->chr[sample];
  for (int time=0; time<nTimes; time++){
    int N = cnaEmit->depths[time][sample][site];
    int n = cnaEmit->reads[time][sample][site];
    double rnd  = cnaEmit->rnd_emit / (double(N)*(cnaEmit->maxRate - cnaEmit->minRate) + 1.0);
    double nrnd = 1.0 - cnaEmit->rnd_emit;
    double val=0;
    if (cnaEmit->mode==3){//Poisson
      double pre1 = double(n) * (log_mass->data[time] + logn[N]) - loggma[n+1];
      double pre2 = mass->data[time] * double(N);
      if (b > 0.0){//modulation
	pre1 += double(n) * lb;
	pre2 *= b;
      }
      for (int level=0; level<nLevels; level++){
	if (prior->data[level] <= 0.0) continue;
	val = pre1 - pre2*tcn[chr][time][level] + double(n)*log_tcn[chr][time][level];
	val = exp(val);
	if (rnd > 0.0) val = val*nrnd + rnd;
	post->data[level] *= val;
	if (Post!=NULL) gsl_matrix_set(Post,time,level,val);
      }
    }
    else if (cnaEmit->mode == 4){//Negative Binomial
      double n1 = double(N)*cnaEmit->shape;//rounding (error?)
      double pre = loggma[int(n1)+n] - loggma[int(n1)] - loggma[n+1] + n1*cnaEmit->log_shape;
      double rate;
      for (int level=0; level<nLevels; level++){
	if (prior->data[level] <= 0.0) continue;
	rate = mass->data[time] * tcn[chr][time][level];
	val  = pre +  double(n) * (log_mass->data[time] + log_tcn[chr][time][level]);
	if (b > 0.0){
	  rate *= b;
	  val += double(n)*lb;
	}
	val  = val - (double(n) + n1) * log(rate + cnaEmit->shape);
	val = exp(val);
	if (rnd > 0.0) val = val*nrnd + rnd;
	post->data[level] *= val;
	if (Post!=NULL) gsl_matrix_set(Post,time,level,val);
      }
    }
  }
}

//*** BAF UPDATE ************************************************************************
void Clone::update_baf( gsl_vector * prior, gsl_vector * post, int sample, int evt, gsl_matrix * Post){
  if (bafEmit->coarse_grained){
    Clone::update_baf_event( prior, post, sample, evt, Post);
  }
  else{
    Clone::update_baf_site( prior, post, sample, evt, Post);
  }
}

void Clone::update_baf_event( gsl_vector * prior, gsl_vector * post, int sample, int evt, gsl_matrix * Post){
  if (bafEmit->mean_tcn==NULL) abort();
  gsl_vector_set_all( post, 0.0);//log-space!
  if (Post != NULL) gsl_matrix_set_all(Post, 0.0);
  double x,y,val;
  double dx = 1.0 / double(bafGrid);
  gsl_vector_view emit;
  for (int time=0; time<nTimes; time++){
    double mntcn = bafEmit->mean_tcn[time][sample][evt];
    y = (1.0 - purity[time]) / mntcn;
    emit = gsl_matrix_row( bafEmitLog[time][sample], evt);
    for (int level=0; level<nLevels; level++){
      if (prior->data[level] <= logzero) continue;
      x = y + clone_spectrum[time][level] / mntcn;
      if ( x < 0.0 || x > 1.0){//outside range...
	val = logzero;
      }
      else{//inside range...
	val = Clone::get_interpolation( x, 0.0, 1.0, dx, &emit.vector);
      }
      post->data[level] += val;
      if (Post!=NULL) gsl_matrix_set(Post,time,level,val);
    }
  }
}


void Clone::update_baf_site( gsl_vector * prior, gsl_vector * post, int sample, int site, gsl_matrix * Post){
  gsl_vector_set_all( post, 1.0);//not log-space!
  if (Post != NULL) gsl_matrix_set_all(Post, 1.0);
  double x, y, mntcn, val;
  unsigned int N, n1, n2;
  double dx = bafEmit->dx;
  int evt   = bafEmit->event_of_idx[sample][site];
  gsl_vector * emit1=NULL, * emit2=NULL;
  for (int time=0; time<nTimes; time++){
    N  = bafEmit->depths[time][sample][site];
    n1 = bafEmit->reads[time][sample][site];
    emit1 = bafEmit->EmitLog[N][n1];
    n2 = N-n1;
    if (n2 != n1) emit2 = bafEmit->EmitLog[N][n2];
    double rnd = bafEmit->rnd_emit / double(N+1);
    mntcn =  bafEmit->mean_tcn[time][sample][evt];
    y = (1.0-purity[time]) / mntcn;
    for (int level=0; level<nLevels; level++){
      if (prior->data[level] <= 0.0) continue;
      x = y + clone_spectrum[time][level] / mntcn;
      if (x > 1.0){//outside range...
	val = 0.0;
      }
      else if (x<1.0 && x>0.0){//inside range...
	val = exp( Clone::get_interpolation( x, 0.0, 1.0, dx, emit1) );
	if (n1 != n2) val += exp( Clone::get_interpolation( x, 0.0, 1.0, dx, emit2) );
      }
      else if ( x==1.0 ){//edges...
	if (n1==N || n2==N){
	  val = 1.0;
	}
	else{
	  val = 0.0;
	}
      }
      else if ( x==0.0 ){
	if (n1==0 || n2==0){
	  val = 1.0;
	}
	else{
	  val = 0.0;
	}
      }
      else{
	abort();
      }
      if (bafEmit->rnd_emit > 0.0) val = val*(1.0-bafEmit->rnd_emit) + rnd;
      if (val!=val || val < 0.0) abort();
      post->data[level] *= val;
      if (Post!=NULL) gsl_matrix_set(Post,time,level,val);
    }
  }
}


//*** SNV UPDATE ************************************************************************
void Clone::update_snv( gsl_vector * prior, gsl_vector * post, int sample, int evt, gsl_matrix * Post){
  if (snvEmit->coarse_grained){
    Clone::update_snv_event( prior, post, sample, evt, Post);
  }
  else if ( !snvEmit->connect ){
    Clone::update_snv_site_ncorr( prior, post, sample, evt, Post);
  }
  else if ( bulk_fix >= 0.0 || (bulk_mean != NULL && bulk_dist == NULL)){
    Clone::update_snv_site_fixed( prior, post, sample, evt, Post);
  }
  else{
    Clone::update_snv_site_nfixed( prior, post, sample, evt, Post);
  }
}

void Clone::update_snv_event( gsl_vector * prior, gsl_vector * post, int sample, int evt, gsl_matrix * Post){
  gsl_vector_set_all( post, 0.0);//log-space!
  if (Post != NULL) gsl_matrix_set_all(Post, 0.0);
  double val=0;
  int chr = snvEmit->chr[sample];
  for (int time=0; time<nTimes; time++){
    double mntcn 
      = (snvEmit->mean_tcn == NULL) 
      ? tcn[chr][time][level_of[chr]] 
      : snvEmit->mean_tcn[time][sample][evt];
    double b = (1.0 - purity[time]) * double(normal_copy[chr]) / mntcn;
    for (int level=0; level<nLevels; level++){
      if (prior->data[level] <= logzero) continue;
      double a = clone_spectrum[time][level] / mntcn;
      if (a<0.0){
	cout<<"here\n";
	abort();
      }
      if (a > 1.0){//outside range...
	int idx  = snvEmit->idx_of_event[sample][evt];
	int nidx 
	  = (evt < snvEmit->nEvents[sample]-1) 
	  ? snvEmit->idx_of_event[sample][evt+1] 
	  : snvEmit->nSites[sample];
	val = logzero * double(nidx-idx);
      }
      else if (bulk_fix == 0.0){
	gsl_vector_view col = gsl_matrix_column( snvEmitLog[time][sample][evt], 0);
	val = Clone::get_interpolation( a, 0.0, 1.0, 1.0/double((&col.vector)->size-1), &col.vector);
      }
      else{//inside range...
	val = Clone::get_interpolation( a, 0.0, 1.0, b, 0.0, 1.0/bulk_min[time][sample][evt], snvEmitLog[time][sample][evt]);
      }
      post->data[level] += val;
      if (Post!=NULL) gsl_matrix_set(Post,time,level,val);
    }
  }
}

void Clone::update_snv_site_ncorr( gsl_vector * prior, gsl_vector * post, int sample, int site, gsl_matrix * Post){
  gsl_vector_set_all( post, 1.0);
  if (Post != NULL) gsl_matrix_set_all(Post, 1.0);
  int chr = snvEmit->chr[sample];
  double x,val=0;
  unsigned int N,n;
  int evt = snvEmit->event_of_idx[sample][site];
  double dx = snvEmit->dx;
  gsl_vector * emit = NULL;
  for (int time=0; time<nTimes; time++){
    N = snvEmit->depths[time][sample][site];
    n = snvEmit->reads[time][sample][site];
    if (N==0){
      continue;
    }
    emit = snvEmit->EmitLog[N][n];
    double rnd  = snvEmit->rnd_emit / double(N+1);
    double nrnd = 1.0 - snvEmit->rnd_emit;
    double mntcn 
      = (snvEmit->mean_tcn != NULL) 
      ? snvEmit->mean_tcn[time][sample][evt]
      : tcn[chr][time][level_of[chr]];
    for (int level=0; level<nLevels; level++){
      if (prior->data[level] <= 0.0) continue;
      if (level==0){//allele frequency of false positives
	x = snv_fpf; 
      }
      else{
	x = clone_spectrum[time][level] / mntcn; 
      }
      if (x > 1.0){//outside range...
	val = 0.0;
      }
      else if (x<1.0 && x>0.0){//inside range...
	val = Clone::get_interpolation( x, 0.0, 1.0, dx, emit);
	val = exp(val);
      }
      else if ( x==1.0 ){//edges...
	val = (n==N) ? 1.0 : 0.0;
      }
      else if ( x==0.0 ){
	val = (n==0) ? 1.0 : 0.0;
      }
      else{
	abort();
      }
      if (val != val || val < 0.0) abort();
      if (rnd > 0.0) val = val*nrnd + rnd;
      post->data[level] *= val;
      if (Post!=NULL) gsl_matrix_set(Post,time,level,val);
    }
  }
}

//Dirac-Delta bulk prior distribution (maybe fixed at zero)...
void Clone::update_snv_site_fixed( gsl_vector * prior, gsl_vector * post, int sample, int site, gsl_matrix * Post){
  gsl_vector_set_all( post, 1.0);
  if (Post != NULL) gsl_matrix_set_all(Post, 1.0);
  int chr = snvEmit->chr[sample];
  double x,y,val=0;
  unsigned int N,n;
  int evt = snvEmit->event_of_idx[sample][site];
  double dx = snvEmit->dx;
  gsl_vector * emit = NULL;
  for (int time=0; time<nTimes; time++){
    N = snvEmit->depths[time][sample][site];
    n = snvEmit->reads[time][sample][site];
    if (N==0){
      continue;
    }
    emit = snvEmit->EmitLog[N][n];
    double rnd  = snvEmit->rnd_emit / double(N+1);
    double nrnd = 1.0 - snvEmit->rnd_emit;
    double mntcn 
      = (snvEmit->mean_tcn == NULL) 
      ? tcn[chr][time][level_of[chr]] 
      : snvEmit->mean_tcn[time][sample][evt];
    double bfix = (bulk_fix >= 0.0) ? bulk_fix : bulk_mean[time][sample][site];
    y =  bfix * (1.0-purity[time]) * double(normal_copy[chr]) / mntcn;
    for (int level=0; level<nLevels; level++){
      if (prior->data[level] <= 0.0) continue;
      x = y + clone_spectrum[time][level] / mntcn; 
      if (x > 1.0){//outside range...
	val = 0.0;
      }
      else if (x<1.0 && x>0.0){//inside range...
	val = Clone::get_interpolation( x, 0.0, 1.0, dx, emit);
	val = exp(val);
      }
      else if ( x==1.0 ){//edges...
	val = (n==N) ? 1.0 : 0.0;
      }
      else if ( x==0.0 ){
	val = (n==0) ? 1.0 : 0.0;
      }
      else{
	abort();
      }
      if (val != val || val < 0.0) abort();
      if (rnd > 0.0) val = val*nrnd + rnd;
      post->data[level] *= val;
      if (Post!=NULL) gsl_matrix_set(Post,time,level,val);
    }
  }
}

//HERE

//non-trivial case with one or more clones
void Clone::update_snv_site_nfixed( gsl_vector * prior, gsl_vector * post, int sample, int site, gsl_matrix * Post){
  gsl_vector_set_all(post, 1.0);
  if (Post != NULL) gsl_matrix_set_all(Post, 1.0);
  int chr = snvEmit->chr[sample];
  double val=0;
  int evt = snvEmit->event_of_idx[sample][site];
  for (int time=0; time<nTimes; time++){//product over time points
    unsigned int N = snvEmit->depths[time][sample][site];
    unsigned int n = snvEmit->reads[time][sample][site];
    gsl_vector * emit = snvEmit->EmitProb[N][n];
    double rnd  = snvEmit->rnd_emit / double(N+1);
    double nrnd = 1.0 - snvEmit->rnd_emit;
    //evaluate likelihood function at a finite grid of clonal frequencies
    int nPts = (nLevels <= 100 ) ? nLevels-1 : 100;
    //grid for linear interpolation//***TEST***
    gsl_vector * mem = gsl_vector_calloc(nPts+1);
    double mntcn 
      = (snvEmit->mean_tcn == NULL) 
      ? tcn[chr][time][level_of[chr]] 
      : snvEmit->mean_tcn[time][sample][evt];
    double maxq = 0.0;
    for (int l=0; l<nLevels; l++){
      if (clone_spectrum[time][l] / mntcn <= 1.0){
	maxq = max( maxq, clone_spectrum[time][l] / mntcn);
      }
    }
    double dq = maxq / double(nPts);
    gsl_vector_view blk = gsl_matrix_row( bulk_dist[time][sample], site);
    double b = (1.0-purity[time])*double(normal_copy[chr]) / mntcn;
    for (int pt=0; pt<=nPts; pt++){
      double a = (nPts == nLevels-1) ? clone_spectrum[time][pt] / mntcn : double(pt)*dq;
      if( N == 0 ){
	val =  1.0;// no observation -> all states equally likely!
      }
      else if (a > 1.0){
	val = 0.0;
      }
      else{
	if ( nPts != nLevels-1 ){
	  val = Clone::trapezoidal( &blk.vector, a, b, emit, 1);//log
	}
	else if (prior->data[pt] > 0.0){
	  val = Clone::trapezoidal( &blk.vector, a, b, emit, 0);//linear
	}
      }	
      gsl_vector_set( mem, pt, val);
    }
    // mutiply into posterior...
    if (nPts == nLevels-1){// we have computed lh at all the correct points
      if (rnd > 0.0){
	gsl_vector_scale(mem, nrnd);
	gsl_vector_add_constant( mem, rnd);
      }
      gsl_vector_mul( post, mem);
      if (Post != NULL){
	gsl_vector_view row = gsl_matrix_row(Post,time);
	gsl_vector_memcpy(&row.vector,mem);
      }
    }
    else{//need to approximate lh via linear interpolation...
      double a,nu,f1,f2;
      int low;	
      // remaining levels...
      for (int level=0; level<nLevels; level++){
	if (prior->data[level] <= 0.0) continue;
	a = clone_spectrum[time][level] / mntcn;	 
	if (a>1.0){
	  val = 0.0;
	}
	else{
	  low = int(a/dq);
	  nu = a/dq - double(low);
	  f1 = gsl_vector_get( mem, low);
	  f2 = (low<nPts) ? gsl_vector_get( mem, low + 1) : f1;
	  val = (1.0-nu)*f1 + nu*f2;//linear interpolation for log-value
	  val = exp(val);
	  if (snvEmit->rnd_emit > 0.0) val = val*(1.0-snvEmit->rnd_emit) + rnd;
	}
	post->data[level] *= val;
	if (Post!=NULL) gsl_matrix_set(Post,time,level,val);
      }
    }
    gsl_vector_free(mem);
  }
}

//only ever used on snvEmit!
double Clone::trapezoidal( gsl_vector * blk, double a, double b, gsl_vector * emit, int get_log){
  int olow = -1;
  int gsb = (int) blk->size  - 1;
  int gse = (int) emit->size - 1;
  double dxb = 1.0 / double(gsb);
  double dxe = 1.0 / double(gse);
  double y;
  double f1=0,f2=0,nu=0,pre=0;
  int low,high;
  //both blk and emit are defined on [0,1].
  if ( gsl_vector_isnull(blk) ) abort();
  double val = 0.0, dval=0;
  for (int j=0; j <= gsb; j++){
    y = a + double(j)*dxb*b;
    if (y > 1.0) break;
    low  = int( y / dxe);
    if ( low != olow){
      f1 = emit->data[low];
      if (low < gse){
	high = low+1;
	f2   = emit->data[high];
      }
      else{
	high = low;
	f2   = f1;
      }
      olow = low;
    }
    // linear interpolation of 1-D distribution between two grid-points...
    pre = (j==0 || j==gsb) ? 0.5 : 1.0;
    nu = y / dxe - double(low);
    dval = pre * blk->data[j] * ( (1.0-nu)*f1 + nu*f2 );
    val += dval;
    //if (dval < 1.0e-6*val) break;
  }
  val *= dxb;
  if (val != val) abort();
  if (get_log==1) val = (val > 0.0) ? log(val) : logzero;
  return(val);
}



//*** INTERPOLATION ****************************************************************************
double Clone::get_interpolation(double x, double xmin, double xmax, double dx, gsl_vector * emit){
  int idx   = int((x-xmin)/dx);
  double nu = (x-xmin)/dx - double(idx);
  double f1,f2;
  int gs = (int) emit->size - 1;
  if (idx >= 0 && idx <= gs){//use log-emission probability
    f1 = emit->data[idx];
    f2 = emit->data[min(gs,idx+1)];
  }
  else{
    abort();
  }
  /*else if (idx==0){//else use linear extra-polation of log
    if (x>xmin){
      f1 = emit->data[1];
      f2 = emit->data[2];
      nu -= 1.0;
    }
    else{
      f1=logzero;
      f2=logzero;
    }
  }
  else{
    if (x<xmax){
      f1 = emit->data[gs-2];
      f2 = emit->data[gs-1];
      nu += 1.0;
    }
    else{
      f1=logzero;
      f2=logzero;
    }
    }*/
  double val = (1.0-nu)*f1 + nu*f2;//linear interpolation (in log-space)
  return(val);
}

//2-D interpolation...
double Clone::get_interpolation(double x, double xmin, double xmax,
				double y, double ymin, double ymax,
				gsl_matrix * emit){
  double dx = (xmax-xmin) / (double(emit->size1 - 1));
  double dy = (ymax-ymin) / (double(emit->size2 - 1));
  int ix1 = int((x-xmin)/dx);
  int iy1 = int((y-ymin)/dy);
  double f,f1,f2,f3,f4;
  double t = (x-(xmin+ix1*dx)) / dx;
  double u = (y-(ymin+iy1*dy)) / dy;  
  f1 = gsl_matrix_get(emit, ix1,   iy1);
  f2 = gsl_matrix_get(emit, ix1+1, iy1);
  f3 = gsl_matrix_get(emit, ix1+1, iy1+1);
  f4 = gsl_matrix_get(emit, ix1,   iy1+1);  
  f  = (1.0-t)*(1.0-u)*f1 + t*(1.0-u)*f2 + t*u*f3 + (1.0-t)*u*f4;
  return(f);
}


//*** EmitLog's for precomputed emissions ******************************************************
void Clone::get_cnaEmitLog(){//only ever used for BAF data
  if (cnaEmitLog != NULL) abort();
  cna_xmin = new double [nTimes];
  cna_xmax = new double [nTimes];
  cnaEmitLog = new gsl_matrix ** [nTimes];
  double lnrnd = (cnaEmit->rnd_emit>0.0) ? log(1.0-cnaEmit->rnd_emit) : 0.0;
  for (int t=0; t<nTimes; t++){
    cnaEmit->get_log = 1;
    cnaEmit->init_range(t);
    cna_xmin[t] = cnaEmit->xmin;
    cna_xmax[t] = cnaEmit->xmax;
    cnaEmit->set_EmitProb(t);
    cnaEmitLog[t] = new gsl_matrix * [cnaEmit->nSamples];
    for (int s=0; s<cnaEmit->nSamples; s++){
      cnaEmitLog[t][s] = gsl_matrix_calloc( cnaEmit->nEvents[s], cnaEmit->gridSize+1);
      int evt;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic, 1) default(shared)  
#endif    
      for ( evt=0; evt<cnaEmit->nEvents[s]; evt++){
	gsl_vector * mem = gsl_vector_calloc(cnaEmit->gridSize+1);
	gsl_vector * rnd = gsl_vector_calloc(cnaEmit->gridSize+1);
	unsigned int n,N;
	double r=0;
	gsl_vector_view emit = gsl_matrix_row( cnaEmitLog[t][s], evt);
	int idxi = cnaEmit->idx_of_event[s][evt];
	int idxf = (evt<cnaEmit->nEvents[s]-1) ? cnaEmit->idx_of_event[s][evt+1] - 1: cnaEmit->nSites[s]-1;
	for(int idx=idxi; idx<=idxf; idx++){
	  n = cnaEmit->reads[t][s][idx];
	  N = cnaEmit->depths[t][s][idx];
	  if (N>0){//observation?
	    if (cnaEmit->bias == NULL){
	      gsl_vector_memcpy( mem, cnaEmit->EmitLog[N][n]);
	    }
	    else{
	      cnaEmit->get_eprob_wBias( mem, cnaEmit->EmitLog[N][n], cnaEmit->bias[s][idx], n, N, 1);
	    }	  
	    if (cnaEmit->rnd_emit > 0.0){//random?
	      gsl_vector_add_constant( mem, lnrnd);
	      r = double(N)*(cnaEmit->maxRate - cnaEmit->minRate);
	      if (cnaEmit->bias != NULL) r *= cnaEmit->bias[s][idx];
	      r =  cnaEmit->rnd_emit / (r+1.0);
	      gsl_vector_set_all( rnd, log(r) );
	      log_vector_add( mem, rnd);
	    }
	    gsl_vector_add( &emit.vector, mem);
	  }
	}
	gsl_vector_free(mem);
	gsl_vector_free(rnd);
      }
    }  
  }
}

void Clone::get_bafEmitLog(){//only ever used for BAF data
  if (bafEmitLog != NULL) abort();
  bafEmitLog = new gsl_matrix ** [nTimes];
  double lrnd  = (bafEmit->rnd_emit>0.0) ? log(bafEmit->rnd_emit) : 0.0;
  double lnrnd = (bafEmit->rnd_emit>0.0) ? log(1.0-bafEmit->rnd_emit) : 0.0;
  for (int t=0; t<nTimes; t++){
    bafEmitLog[t] = new gsl_matrix * [bafEmit->nSamples];
    for (int s=0; s<bafEmit->nSamples; s++){
      bafEmitLog[t][s] = gsl_matrix_calloc( bafEmit->nEvents[s], bafEmit->gridSize+1);      
      int evt;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic, 1) default(shared)
#endif
      for (evt=0; evt<bafEmit->nEvents[s]; evt++){
	unsigned int n,N;
	gsl_vector * mem = gsl_vector_calloc(bafEmit->gridSize+1);
	gsl_vector * rnd = gsl_vector_calloc(bafEmit->gridSize+1);
	gsl_vector_view emit = gsl_matrix_row( bafEmitLog[t][s], evt);
	int idxi = bafEmit->idx_of_event[s][evt];
	int idxf = (evt<bafEmit->nEvents[s]-1) ? bafEmit->idx_of_event[s][evt+1] - 1: bafEmit->nSites[s]-1;
	for(int idx=idxi; idx<=idxf; idx++){
	  n = bafEmit->reads[t][s][idx];
	  N = bafEmit->depths[t][s][idx];
	  if (N>0){//observation?
	    gsl_vector_memcpy( mem, bafEmit->EmitLog[N][n]);
	    if (n != N-n) log_vector_add( mem, bafEmit->EmitLog[N][N-n]);	  
	    if (bafEmit->rnd_emit > 0.0){//random?
	      gsl_vector_add_constant( mem, lnrnd);
	      gsl_vector_set_all( rnd, lrnd - log(double(N+1)) );
	      log_vector_add( mem, rnd);
	    }
	    gsl_vector_add( &emit.vector, mem);
	  }
	}
	gsl_vector_free(mem);
	gsl_vector_free(rnd);
      }      
    }  
  }  
  //cout<<"done."<<endl;
}

void Clone::get_snvEmitLog(){//only ever used for SNV data
  if (!snvEmit->is_set) abort();
  int s1 = snvEmit->gridSize+1;
  int s2 = (bulk_fix == 0.0) ? 1 : (snvEmit->gridSize+1);
  if (snvEmitLog==NULL){//first allocation
    snvEmitLog = new gsl_matrix *** [nTimes];
    for (int t=0;t<nTimes; t++){
      snvEmitLog[t] = new gsl_matrix ** [snvEmit->nSamples];
      for (int s=0; s<snvEmit->nSamples; s++){
	snvEmitLog[t][s] = new gsl_matrix * [snvEmit->nEvents[s]];
	for (int evt=0; evt< snvEmit->nEvents[s]; evt++){
	  snvEmitLog[t][s][evt] = gsl_matrix_calloc( s1, s2);
	}
      }
    }
  }
  Clone::get_bulk_min();
  //double dx = snvEmit->dx;
  double lnrnd = log(1.0 - snvEmit->rnd_emit);
  cout<<endl;
  for (int t=0; t<nTimes; t++){   
    for (int s=0; s<snvEmit->nSamples; s++){   
      printf("\r%2i/%2i %2i/%2i...",t+1,nTimes,s+1,snvEmit->nSamples);
      cout<<flush;
      int evt;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic, 1) default(shared)
#endif
      for (evt=0; evt<snvEmit->nEvents[s]; evt++){
	gsl_matrix * mem     = gsl_matrix_calloc(s1,s2);
	gsl_matrix * rnd_vec = gsl_matrix_calloc(s1,s2);
	gsl_matrix_set_zero(snvEmitLog[t][s][evt]);
	int idxi = snvEmit->idx_of_event[s][evt];
	int idxf = (evt<snvEmit->nEvents[s]-1) ? snvEmit->idx_of_event[s][evt+1] - 1: snvEmit->nSites[s]-1;
	for(int idx=idxi; idx<=idxf; idx++){
	  double val=0;
	  unsigned int N = snvEmit->depths[t][s][idx];
	  if (N==0) continue;//no observation
	  unsigned int n = snvEmit->reads[t][s][idx];
	  gsl_vector * emit = (bulk_prior==NULL) ? snvEmit->EmitLog[N][n] : snvEmit->EmitProb[N][n];
	  double dx = 1.0/double(emit->size-1);
	  double lrnd  = (snvEmit->rnd_emit>0.0) ? log(snvEmit->rnd_emit/double(N+1)) : 0.0;	  
	  //evaluate likelihood function at a finite grid of clonal frequencies
	  double bfix = (bulk_fix >= 0.0) ? bulk_fix : bulk_mean[t][s][idx];
	  gsl_vector_view blk;
	  if (bulk_prior != NULL) blk = gsl_matrix_row( bulk_dist[t][s], idx);
	  //printf("\nbmin = %e\n", bulk_min[t][s][evt]);
	  for (int i=0; i<=snvEmit->gridSize; i++){
	    double a = double(i) / double(snvEmit->gridSize);
	    if (bulk_fix == 0.0){//bulk fixed at zero	      
	      if (a>0.0 && a<1.0){//inside range...
		val = emit->data[i];
	      }
	      else if ( a==1.0 ){//edges...
		val = (n==N) ? 0.0 : logzero;
	      }
	      else if ( a==0.0 ){
		val = (n==0) ? 0.0 : logzero;
	      }
	      gsl_matrix_set( mem, i, 0, val);
	    }
	    else{
	      for (int j=0; j<=snvEmit->gridSize; j++){
		double b = double(j) / (double(snvEmit->gridSize) * bulk_min[t][s][evt]);
		if (bulk_prior==NULL){//bulk point estimate...	      
		  double x =  a + bfix * b;
		  if (x > 1.0){//outside range...
		    val = logzero;
		  }
		  else if (x>0.0 && x<1.0){//inside range...
		    val = Clone::get_interpolation( x, 0.0, 1.0, dx, emit);
		  }
		  else if ( x==1.0 ){//edges...
		    val = (n==N) ? 0.0 : logzero;
		  }
		  else if ( x==0.0 ){
		    val = (n==0) ? 0.0 : logzero;
		  }
		}
		else{//bulk full distribution...		
		  val = Clone::trapezoidal( &blk.vector, a, b, emit, 1);//log-space
		}	
		gsl_matrix_set( mem, i, j, val);
		//printf("%e ", val);
	      }
	      //cout<<endl;
	    }
	  }
	  // mutiply into posterior...
	  if (snvEmit->rnd_emit > 0.0){
	    gsl_matrix_add_constant( mem, lnrnd);
	    gsl_matrix_set_all( rnd_vec, lrnd);
	    log_matrix_add( mem, rnd_vec);
	  }
	  gsl_matrix_add( snvEmitLog[t][s][evt], mem);
	}
	gsl_matrix_free(mem);
	gsl_matrix_free(rnd_vec);
      }
    }  
  }
  cout<<endl;
}

