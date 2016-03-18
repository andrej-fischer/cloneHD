//emission.cpp

//own headers...
#include "emission.h"

#define PI 3.1415926

using namespace std;

// Constructor
Emission::Emission(){
  is_set = 0;
  dx=0;
  EmitProb_set = 0;
  dist_set = 0;
  range_set = 0;
  shape    = 1.0;
  reflect  = 0;
  mode     = 0;// 1: binomial, 2: beta-binomial, 3: Poisson, 4: negative binomial
  xmin     = 0;
  xmax     = 0;
  ymin     = 0;
  ymax     = 0;
  xgrid = NULL;
  ygrid = NULL;
  maxRate  = 0; 
  get_log  = 0;//get log of emission probability?
  get_der  = 0;//get derivatives?
  get_mv   = 0;//get mean and variance?
  bias     = NULL;
  log_bias = NULL;
  mask     = NULL;
  pjump          = NULL;
  mean_tcn       = NULL;
  av_cn          = NULL;
  idx_of_event   = NULL;
  event_of_idx   = NULL;
  nEvents        = NULL;
  Event_of_idx   = NULL;
  log_space      = 0;
  coarse_grained = 0;
  idx_to_Event_mapped=0;
  connect=0;
  nObs = NULL;
}


// real constructor
void Emission::set(int ntimes, vector<int>& Chrs, vector<int>& nsites, int grid){
  if ((int) Chrs.size() != (int) nsites.size()){
    cout<<"FATAL ERROR (use gdb to locate).\n";
    abort();
  }
  nSamples  = (int) nsites.size();
  nTimes    = ntimes;
  nSites  = new int [nSamples]; 
  chr     = new int [nSamples];
  chrs.clear();
  total_loci   = 0;
  maxchr = 0;
  for (int s=0; s<nSamples; s++){
    chr[s]     = Chrs[s];
    chrs.insert(chr[s]);
    maxchr      = max( maxchr, chr[s]);
    nSites[s]   = nsites[s];
    total_loci += nSites[s];
  }
  idx_of = new int [maxchr+1];
  for (int i=0; i<=maxchr; i++)  idx_of[i] = -1;
  for (int s=0; s<nSamples; s++) idx_of[chr[s]] = s;
  // set xgrid...
  gridSize = grid;
  xgrid  = new double [gridSize+1];
  ygrid  = new double [gridSize+1];
  reads  = new unsigned int ** [nTimes];
  depths = new unsigned int ** [nTimes];
  for (int t=0; t<nTimes; t++){
    reads[t]  = new unsigned int * [nSamples];
    depths[t] = new unsigned int * [nSamples];
    for (int s=0; s<nSamples; s++){
      reads[t][s]    = new unsigned int [nSites[s]];
      depths[t][s]   = new unsigned int [nSites[s]];
    }
  }
  loci = new unsigned int * [nSamples];
  mask = new unsigned int * [nSamples];
  pjump = new double * [nSamples];
  for (int s=0; s<nSamples; s++){
    loci[s] = new unsigned int [nSites[s]];
    mask[s] = new unsigned int [nSites[s]];
    pjump[s] = new double [nSites[s]];
    for (int l=0; l<nSites[s];l++) pjump[s][l] = 0.0;
  }
  Emission::reset_mask();
  Emission::init_events();
  is_set=1;
}

void Emission::init_events(){
  nEvents = new int [nSamples];//coarse-grained events
  Event_of_idx   = new unsigned int * [nSamples];// map from idx to cnv-event 
  idx_of_event   = new unsigned int * [nSamples];// maps within
  event_of_idx   = new unsigned int * [nSamples];// maps within
  total_events   = 0;
  for (int s=0; s<nSamples; s++){
    nEvents[s]    = nSites[s];//first, everything is an event
    total_events += nSites[s];
  }
  for (int s=0; s<nSamples; s++){//first, everything is an event
    event_of_idx[s] = new unsigned int [nSites[s]];
    Event_of_idx[s] = new unsigned int [nSites[s]];
    idx_of_event[s] = new unsigned int [nSites[s]];
    for (int idx=0; idx<nSites[s]; idx++){
      Event_of_idx[s][idx] = 0;
      event_of_idx[s][idx] = idx;
      idx_of_event[s][idx] = idx;
    }
  }
}


//map each observations to an event in another data track...
void Emission::map_idx_to_Event(Emission * Emit, int sample){
  if (Emit->is_set == 0) abort(); 
  if (Emit->chrs.count(chr[sample]) == 0) abort();
  int Sample = Emit->idx_of[chr[sample]];
  int Event = 0;
  int Idx   = Emit->idx_of_event[Sample][Event];
  int Locus = Emit->loci[Sample][Idx];
  int idx   = 0;
  int locus = loci[sample][idx];
  while(locus < Locus){//left over-hang
    Event_of_idx[sample][idx] = 0;
    idx++;
    if (idx == nSites[sample]) break;  
    locus = loci[sample][idx];
  }
  for ( Event = 1; Event < Emit->nEvents[Sample]; Event++){
    if (idx == nSites[sample]) break; 
    Idx   = Emit->idx_of_event[Sample][Event];
    Locus = Emit->loci[Sample][Idx];
    while(locus < Locus){
      Event_of_idx[sample][idx]= Event - 1;
      idx++;
      if (idx == nSites[sample]) break;  
      locus = loci[sample][idx];
    }     
  }
  while( idx < nSites[sample]){//right over-hang
    Event_of_idx[sample][idx] = Emit->nEvents[Sample] - 1;
    idx++;
  }
}

void Emission::map_jumps(Emission * Emit){
  for (int s=0; s<nSamples; s++){
    int Sample = Emit->idx_of[chr[s]];
    for (int idx=0; idx<nSites[s]; idx++) pjump[s][idx]=0;
    int old = Event_of_idx[s][0];
    for (int idx=0; idx<nSites[s]; idx++){
      int Evt = Event_of_idx[s][idx];
      if (Evt != old){
	old = Evt;
	int Idx = Emit->idx_of_event[Sample][Evt];
	pjump[s][idx] = Emit->pjump[Sample][Idx];
      }
    }
  }
}

void Emission::add_break_points_via_jumps(Emission * Emit, double pmin){
  for (int s=0; s<nSamples; s++){
    int Sample = Emit->idx_of[chr[s]];
    int old    = Event_of_idx[s][0];
    for (int idx=0; idx<nSites[s]; idx++){
      int Evt = Event_of_idx[s][idx];
      if (Evt != old){
	old = Evt;
	int Idx = Emit->idx_of_event[Sample][Evt];
	double pj = pjump[s][idx];
	double PJ = Emit->pjump[Sample][Idx];
	if (pj>0.0){
	  pj = 1.0 - (1.0-pj)*(1.0-PJ);
	}
	else{
	  pj = pmin;
	}
	pjump[s][idx] = pj;
      }
    }
  }
}



void Emission::get_events_via_jumps(){
  total_events = 0;
  for (int s=0; s<nSamples; s++){
    nEvents[s] = 1;
    for (int idx=1; idx<nSites[s]; idx++){
      if ( pjump[s][idx] > 0.0) nEvents[s]++;
    }
    total_events += nEvents[s];
    if (idx_of_event[s] != NULL)   delete [] idx_of_event[s];
    idx_of_event[s] = new unsigned int [nEvents[s]];
    int evt=0;
    event_of_idx[s][0] = 0;
    idx_of_event[s][0] = 0;
    for (int idx=1; idx<nSites[s]; idx++){
      if ( pjump[s][idx] > 0.0){
	evt++;
	idx_of_event[s][evt] = idx;
      }
      event_of_idx[s][idx] = evt;
    }
  }
}

void Emission::allocate_bias(){//only once...
  if (bias != NULL || log_bias != NULL) abort(); 
  bias     = new double * [nSamples];
  log_bias = new double * [nSamples];
  for (int s=0; s<nSamples; s++){
    bias[s]     = new double [nSites[s]];
    log_bias[s] = new double [nSites[s]];
  }
}


void Emission::allocate_mean_tcn(){
  if (mean_tcn != NULL) abort();
  mean_tcn = new double ** [nTimes];
  for (int t=0; t<nTimes; t++){
    mean_tcn[t] = new double * [nSamples];
    for (int s=0; s<nSamples; s++){
      if (nEvents[s] == 0){ 
	mean_tcn[t][s] = NULL;
      }
      else{
	mean_tcn[t][s] = new double [nEvents[s]];
      }
    }
  }
}


void Emission::allocate_av_cn(int maxtcn){//only once!
  if (av_cn != NULL) abort();
  av_cn = new double *** [nTimes];
  for (int t=0; t<nTimes; t++){
    av_cn[t] = new double ** [nSamples];
    for (int s=0; s<nSamples; s++){
      if (nEvents[s] == 0){ 
	av_cn[t][s] = NULL;
      }
      else{
	av_cn[t][s] = new double * [nEvents[s]];
	for (int evt=0; evt<nEvents[s]; evt++){
	  av_cn[t][s][evt] = new double [maxtcn+1];
	}
      }
    }
  }
}

void Emission::get_nObs(){
  if (nObs!=NULL) abort();
  nObs = new unsigned int ** [nTimes];
  for (int t=0; t<nTimes; t++){
    nObs[t] = new unsigned int * [nSamples];
    for (int s=0; s<nSamples; s++){
      if (nEvents[s] == 0){
	nObs[t][s] = NULL; 
	continue;
      }
      nObs[t][s] = new unsigned int [nEvents[s]];
    }
  }
  for (int s=0; s<nSamples; s++){
    for (int evt=0; evt<nEvents[s]; evt++){
      int first = idx_of_event[s][evt];     
      int last 
	= (evt < nEvents[s]-1) 
	? idx_of_event[s][evt+1] - 1 
	: nSites[s] - 1; 
      unsigned int n,N;
      for (int t=0; t<nTimes; t++){
	nObs[t][s][evt] = 0;
	for (int idx=first; idx<=last; idx++){
	  n = reads[t][s][idx];
	  N = depths[t][s][idx];
	  if (N==0){
	    if (n>0) abort();
	    continue;
	  }
	  nObs[t][s][evt]++;
	}
      }
    }
  }
}


//Destructor
Emission::~Emission(){
  if (pjump != NULL){
    for (int s=0; s<nSamples; s++) delete [] pjump[s];
    delete [] pjump;
  }
  if (mean_tcn != NULL){
    for (int t=0; t<nTimes; t++){
      for (int s=0; s<nSamples; s++){
	if (mean_tcn[t][s] != NULL) delete [] mean_tcn[t][s];
      }
      delete [] mean_tcn[t];
    }
    delete [] mean_tcn;
  }
  if (av_cn != NULL){
    for (int t=0; t<nTimes; t++){
      for (int s=0; s<nSamples; s++){
	if (av_cn[t][s] != NULL){
	  for (int evt=0; evt<nEvents[s]; evt++){
	    delete [] av_cn[t][s][evt];
	  }
	  delete [] av_cn[t][s];
	}
      }
      delete [] av_cn[t];
    }
    delete [] av_cn;
  }
  if (bias != NULL){
    for (int s=0; s<nSamples; s++) delete [] bias[s];
    delete [] bias;
  }
  if (log_bias != NULL){
    for (int s=0; s<nSamples; s++) delete [] log_bias[s];
    delete [] log_bias;
  }
  if (is_set==1){
    Emission::clear();
  }
}

void Emission::clear(){
  if (is_set==0) abort();
  for (int t=0; t<nTimes; t++){
    for (int s=0; s<nSamples; s++){
      delete [] reads[t][s];
      delete [] depths[t][s];
    }
    delete reads[t];
    delete depths[t];
  }
  for (int s=0; s<nSamples; s++){
    delete [] loci[s];
    delete [] mask[s];
    delete [] dist[s];
  }
  delete [] reads;
  delete [] depths;
  delete [] loci;
  delete [] mask;
  delete [] dist;
  delete [] xgrid;
  delete [] ygrid;
  delete [] nSites;
  delete [] nEvents;
  for (int s=0; s<nSamples; s++){
    delete [] event_of_idx[s];
    if (idx_of_event[s] != NULL)  delete [] idx_of_event[s];
  }
  delete [] event_of_idx;
  delete [] idx_of_event;
  if (EmitProb_set==1) Emission::delete_old_Emit();
  is_set   = 0;
  dist_set = 0;
}





void Emission::delete_old_Emit(){
  //unordered_map< unsigned int, unordered_map< unsigned int, gsl_vector*> >::iterator it1;
  map< unsigned int, map< unsigned int, gsl_vector*> >::iterator it1;
  for (it1 = EmitProb.begin(); it1 != EmitProb.end(); ++it1){
    //unordered_map< unsigned int, gsl_vector* >::iterator it2;
    map< unsigned int, gsl_vector* >::iterator it2;
    for (it2 = (EmitProb[it1->first]).begin(); it2 != (EmitProb[it1->first]).end(); ++it2){
      gsl_vector_free( (EmitProb[it1->first])[it2->first]);
    }
    (EmitProb[it1->first]).clear();
  }
  EmitProb.clear();
  //
  if (get_log==1){
    for (it1 = EmitLog.begin(); it1 != EmitLog.end(); ++it1){
      //unordered_map< unsigned int, gsl_vector* >::iterator it2;
      map< unsigned int, gsl_vector* >::iterator it2;
      for (it2 = (EmitLog[it1->first]).begin(); it2 != (EmitLog[it1->first]).end(); ++it2){
	gsl_vector_free( (EmitLog[it1->first])[it2->first]);
      }
      (EmitLog[it1->first]).clear();
    }
    EmitLog.clear();
  }
  EmitProb_set = 0;
}


void Emission::reset_mask(){
  for (int s=0; s<nSamples; s++){
     for (int l=0; l<nSites[s]; l++){
       mask[s][l] = 1;
     }
  }
}



void Emission::set_dist(){
  if (loci == NULL){
    cout<<"ERROR-1 in Emission::set_dist()\n";
    exit(1);
  }
  dist = new unsigned int * [nSamples];
  total_dist=0.0;
  vector<double> distances;
  for (int s=0; s<nSamples; s++){
    dist[s] = new unsigned int [nSites[s]];
    dist[s][0] = 0;
    for (int l=1; l <nSites[s]; l++){
      if (mask[s][l] == 1){
	int k=l;
	unsigned distance = 0;
	while (k>0){
	  distance += loci[s][k] - loci[s][k-1];
	  if (mask[s][k-1]==1) break;
	  k--;
	}
	if (k==0 && mask[s][k] == 0){
	  dist[s][l] = 0;
	}
	else{
	  dist[s][l] = distance;
	  distances.push_back(double(distance));
	  dist_count[distance] += 1;
	  total_dist += distance;
	}
      }
      else{
	dist[s][l] = 0;
      }
    }
  }
  gsl_sort(distances.data(),1,distances.size());
  median_dist = gsl_stats_quantile_from_sorted_data ( distances.data(), 1, distances.size(), 0.5);
  distances.clear();
  // get the distribution of distances...
  map<unsigned int,int>::iterator it;
  for (it = dist_count.begin(); it != dist_count.end(); ++it){
    if (it->second > 1){
      frequent_dist.insert(pair<unsigned int,int>(it->first,it->second));
    }
  }
  dist_set = 1;
}


void Emission::set_grid(){
  dx = (xmax-xmin) / double(gridSize);
  dy = (ymax-ymin) / double(gridSize);
  for (int i=0; i<= gridSize; i++){
    xgrid[i] = xmin + double(i) * dx;
    ygrid[i] = ymin + double(i) * dy;
  }
}


//get initial values for xmin, xmax etc.
void Emission::init_range(int time){
  if (mode==1 || mode==2){
    xmin = 0.0;
    ymin = 0.0;
    xmax = 1.0;
    ymax = 1.0;
  }
  else if (mode==3 || mode==4){
    vector<double> allx;
    vector<double> ally;
    double x,y;
    maxRate=0;
    minRate=1.0e6;
    for (int t=0; t<nTimes; t++){
      if (time >= 0 && time != t) continue;
      for (int s=0; s<nSamples; s++){
	for (int i=0; i < nSites[s]; i++){
	  if (mask[s][i] == 0) continue;
	  if ( depths[t][s][i] > 0 ){
	    y = double(reads[t][s][i]) / double(depths[t][s][i]);
	    x = (bias != NULL) ? y / bias[s][i] : y;
	    ally.push_back( y );
	    if (bias != NULL)  allx.push_back( x );
	    maxRate = max( maxRate, x);
	    minRate = min( minRate, x);
	  }
	}
      }
    }
    //
    gsl_sort( ally.data(), 1, ally.size());
    double eps = 1.0e-3;
    ymin = 0.01 + gsl_stats_quantile_from_sorted_data ( ally.data(), 1, ally.size(), eps);
    ymax = gsl_stats_quantile_from_sorted_data ( ally.data(), 1, ally.size(), 1.0-eps);
    if (bias != NULL){
      eps = 1.0e-3;
      gsl_sort( allx.data(), 1, allx.size());
      xmin = 0.01 + gsl_stats_quantile_from_sorted_data ( allx.data(), 1, allx.size(), eps);
      xmax = gsl_stats_quantile_from_sorted_data ( allx.data(), 1, allx.size(), 1.0-eps);
    }
    else{
      xmin = ymin;
      xmax = ymax;
    }
  }
  Emission::set_grid();
  range_set=1;
}




//emission probability as a function of total freq x
void Emission::set_EmitProb(int time){//if time < 0, do for all!
  if (mode == 0){
    printf("ERROR-1 in Emission::set_EmitProb(): mode not set.\n");
    exit(1);
  }
  if (mode==3||mode==4) reflect=0;
  if (range_set == 0) Emission::init_range(time);
  //delete old stuff
  if (EmitProb_set == 1) Emission::delete_old_Emit();
  //
  for (int t=0; t<nTimes; t++){
    if(time >= 0 && t != time) continue;
    for (int s=0; s<nSamples; s++){
      unsigned int n1,n2,N;
      for (int i=0; i < nSites[s]; i++){
	//if (mask[s][i] == 0) continue;
	N  = depths[t][s][i];
	n1 = reads[t][s][i];
	n2 = (reflect == 1) ? N-n1 : n1;	
	//test
	if ( N==0 ){//no observation!
	  if (n1==0) continue;
	  abort();
	}
	else if (EmitProb.count(N) == 1 && EmitProb[N].count(n1) == 1){
	  continue;
	}
	else{//allocate
	  EmitProb[N][n1] = gsl_vector_alloc(gridSize+1);
	  if (n1!=n2) EmitProb[N][n2] = gsl_vector_alloc(gridSize+1);
	  if (get_log == 1){//logarithm of likelihood?
	    EmitLog[N][n1] = gsl_vector_alloc(gridSize+1);
	    if (n1!=n2) EmitLog[N][n2] = gsl_vector_alloc(gridSize+1);
	  }
	}//now get the actual values...
	if (mode == 1){
	  Emission::binomial(N,n1);
	  if (n1!=n2) Emission::binomial(N,n2);
	}
	else if (mode == 2){
	  Emission::beta_binomial(N,n1);
	  if (n1!=n2) Emission::beta_binomial(N,n2);
	}
	else if (mode==3){
	  Emission::poisson(N,n1);
	}
	else if (mode==4){
	  Emission::negative_binomial(N,n1);
	}
	else{
	  printf("ERROR-5 in Emission::set_EmitProb(): mode (%i) not set.\n", mode);
	  exit(1);
	}
      }
    }
  }
  EmitProb_set = 1;
}

//*** BINOMIAL ***
void Emission::binomial(int N, int n){
  double x,p;
  double pre = gsl_sf_lngamma(double(N+1));
  pre -= gsl_sf_lngamma(double(n+1));
  pre -= gsl_sf_lngamma(double(N-n+1));
  for (int j=0; j<=gridSize; j++){
    x = double(j)*dx;
    if (x==0.0){
      p = (n==0) ? 1.0 : 0.0;
    }
    else if (x==1.0){
      p = (n==N) ? 1.0 : 0.0;
    }
    else if (N==0){
      if (n==0){
	p = 1.0;
      }
      else{
	printf("ERROR-2 in Emission::set_EmitProb()\n");
	exit(1);
      }
    }
    else{
      p = gsl_ran_binomial_pdf(n, x, N);
    }
    if (p<0.0 || p!= p){
      printf("ERROR-3 in Emission::set_EmitProb(): %e\n", p);
    }
    gsl_vector_set( EmitProb[N][n], j, p);
    //logarithm of emission probability...
    if (get_log==1){
      if (x>0.0 && x< 1.0){
	gsl_vector_set(EmitLog[N][n], j, pre + double(n)*log(x) + double(N-n)*log(1.0-x));
      }
      else if ( (x==0.0 && n==0) || (x==1.0 && n==N) ){
	gsl_vector_set(EmitLog[N][n], j, 0.0);
      }
      else{
	gsl_vector_set(EmitLog[N][n], j, -1.0e6);
      }
    }	  
  }
}


//*** BETA-BINOMIAL ***
void Emission::beta_binomial(int N, int n){
  double x,p;
  double p0 = gsl_sf_lngamma(double(N+1)) - gsl_sf_lngamma(double(n+1)) - gsl_sf_lngamma(double(N-n+1)); 
  for (int j=0; j<=gridSize; j++){
    x = double(j)*dx;
    if (x==0.0){
      gsl_vector_set(EmitProb[N][n], j, n==0 ? 1.0 : 0.0);
      if (get_log == 1) gsl_vector_set(EmitLog[N][n], j, n==0 ? 0.0 : -1.0e6);
    }
    else if (x==1.0){
      gsl_vector_set(EmitProb[N][n], j, n==N ? 1.0 : 0.0);
      if (get_log == 1) gsl_vector_set(EmitLog[N][n], j, n==N ? 0.0 : -1.0e6);
    }
    else{
      p  = p0 + gsl_sf_lnbeta(double(n) + shape*x, double(N-n) + shape*(1.0-x));
      p -= gsl_sf_lnbeta( shape*x,  shape*(1.0-x));
      if (p > 0.0) abort();
      if (get_log==1){
	gsl_vector_set(EmitLog[N][n], j, p);
      }
      gsl_vector_set( EmitProb[N][n], j, exp(p));
    } 
  }	
}

//POISSON
void Emission::poisson(int N, int n){
  double y,p;
  double pre = gsl_sf_lngamma( double(n) + 1.0);
  if (pre != pre) abort();
  for (int j=0; j<=gridSize; j++){
    y = ygrid[j];
    if (y>0.0){
      p  = -y*double(N) + n*(log(y)+log(double(N))) - pre;
    }
    else{
      p = (n==0) ? 0.0 : -1.0e6;
    }
    if (get_log==1){
      gsl_vector_set(EmitLog[N][n], j, p);
    }
    gsl_vector_set(EmitProb[N][n], j, exp(p));
  }
}



//negative binomial emission model
void Emission::negative_binomial(int N, int n){
  //TBD: N==0, x==0
  double y,p;
  double r = double(N)*shape;
  double pre = gsl_sf_lngamma( double(n) + r) - gsl_sf_lngamma(double(n+1)) - gsl_sf_lngamma(r);
  pre += r * log(r);
  //printf("n=%i N=%i shape=%e\n", n, N, shape);
  for (int j=0; j<=gridSize; j++){
    y = ygrid[j];
    if (y>0.0){
      p  = pre - (r + double(n))*log(r+y*double(N)) + double(n)*log(y*double(N));
    }
    else{
      p = (n==0) ? 0.0 : -1.0e6;
    }
    if (get_log==1){
      gsl_vector_set(EmitLog[N][n], j, p);
    }
    gsl_vector_set(EmitProb[N][n], j, exp(p));
    //printf("%.4f %.3e\n", y, p);
  }
  //exit(0);
}



double Emission::get_single_EmitLog(double x, unsigned int n, unsigned int N){
  double p=0,pre=0;
  if (mode==1){        
    if ( (x==0.0 && n==0) || (x==1.0 && n==N) ){
      return(0.0);
    }
    else if (x>0.0 && x< 1.0){
      pre  = gsl_sf_lngamma(double(N+1));
      pre -= gsl_sf_lngamma(double(n+1));
      pre -= gsl_sf_lngamma(double(N-n+1));
      return( pre + double(n)*log(x) + double(N-n)*log(1.0-x) );
    }
    else{
      return(-1.0e6);
    }
  }
  else if (mode==2){    
    if ( (x==0.0 && n==0) || (x==1.0 && n==N) ){
      return(0.0);
    }
    else if (x>0.0 && x< 1.0){
      pre = gsl_sf_lngamma(double(N+1)) - gsl_sf_lngamma(double(n+1)) - gsl_sf_lngamma(double(N-n+1)); 
      p  = pre + gsl_sf_lnbeta(double(n) + shape*x, double(N-n) + shape*(1.0-x));
      p -= gsl_sf_lnbeta( shape*x,  shape*(1.0-x));
      if (p>0.0){
	abort();
      }
      return(p);
    }
    else{
      return(-1.0e6);
    }
  }
  else if (mode == 3){
    if (x>0.0){
      pre = gsl_sf_lngamma( double(n) + 1.0);
      p  = -x*double(N) + n*(log(x)+log(double(N))) - pre;
      return(p);
    }
    else if (n==0 && (x==0.0 || N==0)){
      return(0.0);
    }
    else{
      return(-1.0e6);
    }
  }
  else if (mode == 4){
    double r = double(N)*shape;
    pre = gsl_sf_lngamma( double(n) + r) - gsl_sf_lngamma(double(n+1)) - gsl_sf_lngamma(r);
    pre += r * log(r);
    if (x>0.0){
      p  = pre - (r + double(n))*log(r+x*double(N)) + double(n)*log(x*double(N));
      return(p);
    }
    else if (n==0 && (x==0.0 || N==0)){
      return(0.0);
    }
    else{
      return(-1.0e6);
    }
  }
  else{
    return(0);
  }
}


void Emission::get_eprob_wBias( gsl_vector * eprob, 
				gsl_vector * emit, 
				double b, 
				unsigned int n, 
				unsigned int N, 
				int get_log){
  double y,nu,f,f1,f2;
  int idx;
  for (int i=0; i <= gridSize; i++){
    y = b*xgrid[i];
    if (y <= ymin){
      nu  = (y - ymin) / dy;
      f1  = emit->data[0];
      f2  = emit->data[1];
      if (f1<f2){
	f = (1.0-nu)*f1 + nu*f2;
      }
      else{
	f = get_single_EmitLog( y, n, N);
      }
    }
    else if (y >= ymax){
      nu  = (y - ymin) / dy - double(gridSize);
      f1  = emit->data[gridSize-1]; 
      f2  = emit->data[gridSize];
      if (f2<f1){
	f = (1.0-nu)*f1 + nu*f2;
      }
      else{
	f = get_single_EmitLog( y, n, N);
      }
    }
    else{
      idx = int( (y - ymin) / dy);
      nu  = (y - ymin) / dy - double(idx);
      f1  = emit->data[idx];
      f2  = emit->data[idx+1];
      f = (1.0-nu)*f1 + nu*f2;
    }   
    if (get_log==1){
      eprob->data[i] = f; 
    }
    else{
      eprob->data[i] = exp(f);
    }
  }
}




double Emission::get_pval(int time, int sample, int site, double mean){
  unsigned int n = reads[time][sample][site];
  unsigned int N = depths[time][sample][site];
  if (reflect) n = min(n,N-n);
  double pval=0;
  if (mode==1){//binomial
    if (mean==0.0){
      pval = (n==0) ? 1.0 : 0.0;
    }
    else if (mean < 1.0){
      double p = gsl_cdf_binomial_P(n,mean,N);
      double q = gsl_cdf_binomial_Q(n,mean,N);
      q += gsl_ran_binomial_pdf(n,mean,N);
      pval = min(p,q);
    }
    else{
      pval = (n==N) ? 1.0 : 0.0;
    }
  }
  else if (mode==2){//beta-binomial
    double pre = gsl_sf_lngamma(double(N+1)) - gsl_sf_lnbeta( shape*mean,  shape*(1.0-mean));
    double pn = pre - gsl_sf_lngamma(double(n+1)) - gsl_sf_lngamma(double(N-n+1)); 
    pn += gsl_sf_lnbeta(double(n) + shape*mean, double(N-n) + shape*(1.0-mean));
    double q=0.0,p=0.0;
    if (n >= int(0.5*double(N))){
      for (unsigned int k=n+1;k<=N;k++){
	double y = pre - gsl_sf_lngamma(double(k+1)) - gsl_sf_lngamma(double(N-k+1));
	y += gsl_sf_lnbeta(double(k) + shape*mean, double(N-k) + shape*(1.0-mean));
	q += exp(y);
      }
      p = 1.0 - q - pn;
    }
    else{
      for (unsigned int k=n-1;k>=0;k--){
	double y  = pre - gsl_sf_lngamma(double(k+1)) - gsl_sf_lngamma(double(N-k+1));
	y += gsl_sf_lnbeta(double(k) + shape*mean, double(N-k) + shape*(1.0-mean));
	p += exp(y);
      }
      q = 1.0 - p - pn;
    }
    pval = min(p,q);
  }
  else if (mode==3){//poisson
    double p = gsl_cdf_poisson_P(n,mean*double(N));
    double q = gsl_cdf_poisson_Q(n,mean*double(N));
    q += gsl_ran_poisson_pdf(n,mean*double(N));
    pval = min(p,q);
  }
  else if (mode==4){//negative binomial
    double p = gsl_cdf_negative_binomial_P( n, shape / (mean+shape), double(N)*shape);
    double q = gsl_cdf_negative_binomial_Q( n, shape / (mean+shape), double(N)*shape);
    q += gsl_ran_negative_binomial_pdf( n, shape / (mean+shape), double(N)*shape);
    pval = min(p,q);
  }
  return(pval);
}



void Emission::set_pjump(double jp){
  if ( dist == NULL || dist_set == 0) abort();
  if (pjump == NULL) abort();
  double y = (jp < 1.0) ? log(1.0-jp) : 0.0;
  for (int s=0; s<nSamples; s++){
    pjump[s][0] = 0.0;
    for (int idx=1; idx<nSites[s]; idx++){
      if (jp == 0.0){
	pjump[s][idx] = 0.0;
      }
      else if (jp < 1.0){
	pjump[s][idx] = 1.0 - exp( double(dist[s][idx]) * y );
      }
      else if (jp == 1.0){
	pjump[s][idx] = 1.0;
      }
      else{
	abort();
      }
    }
  }
}


void Emission::coarse_grain_jumps( int sample, double min_jump, int range){
  std::map<int,double>::iterator it, max_it, first_it, last_it;
  int L = nSites[sample];
  double * cg_pjump = new double [L];
  for (int i=0; i<L; i++) cg_pjump[i] = 0.0;
  std::map<int,double> remaining;
  std::map<int,int> idxOf;
  pjump[sample][0] = 0;
  int ct=0;
  for (int i=0; i<L; i++){
    if (pjump[sample][i] <= 0) continue;
    remaining.insert(std::pair<int,double>( ct, pjump[sample][i]));
    idxOf.insert(std::pair<int,int>(ct,i));
    ct++;
  }
  //get global minimum...
  first_it = remaining.begin();
  last_it  = remaining.end();
  if (first_it==last_it){//no jumps whatsoever!
    for (int i=0; i<L; i++) pjump[sample][i] = cg_pjump[i];
    delete [] cg_pjump;
    remaining.clear();
    idxOf.clear();
    return;
  }
  it = std::min_element( first_it, last_it, value_comparer);
  double gmin = it->second;
  double p=0,p0=0,pmin=0,ps=0;
  while (remaining.size() > 0){
    max_it = std::max_element( remaining.begin(), remaining.end(), value_comparer);
    p0  = max_it->second;
    if (p0 < min_jump) break;//pjump at peak is below min_jump
    ps = 1.0 - p0;
    first_it = max_it;
    last_it  = max_it;
    pmin  = p0;
    //above...
    it = max_it;
    int incl=0,done=0;
    while( it != remaining.end() ){
      it++;
      p = it->second;
      if ( p>10.0*pmin || it->first != last_it->first+1 || p<10.0*gmin){
	done = 1;
      }
      else if (incl < range){
	ps *= 1.0 - p;
	incl++;
      }
      pmin = min(pmin,p);
      last_it = it;
      if (done) break;
    }
    //below...
    it   = max_it;
    pmin = p0;
    incl = 0;
    done = 0;
    while(it != remaining.begin()){
      it--;
      p = it->second;
      if( p>10.0*pmin || it->first != first_it->first-1 || p<10.*gmin) done=1;
      else if (incl<range){
	ps *= 1.0 - p;
	incl++;
      }
      first_it = it;
      pmin  = min(pmin,p);
      if(done) break;
    }
    //set effective pjump
    cg_pjump[idxOf[max_it->first]] = 1.0-ps;
    //remove...
    remaining.erase( first_it, last_it);
  }
  //fill
  for (int i=0; i<L; i++) pjump[sample][i] = cg_pjump[i];
  delete [] cg_pjump;
  remaining.clear();
  idxOf.clear();
}


bool value_comparer(std::map<int,double>::value_type &i1, std::map<int,double>::value_type &i2){
  return i1.second < i2.second;
}

