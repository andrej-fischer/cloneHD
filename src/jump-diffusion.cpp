//jump-diffusion.cpp

//own headers...
#include "emission.h"
#include "log-space.h"
#include "jump-diffusion.h"
#include "common-functions.h"

#define PI 3.1415926

// Constructor
JumpDiffusion::JumpDiffusion( Emission * emit, int t){
  myEmit   = emit;
  time     = t;
  nSamples = myEmit->nSamples;
  nSites   = myEmit->nSites;
  dist     = myEmit->dist;
  loci     = myEmit->loci;
  mask     = myEmit->mask;
  mode     = myEmit->mode;
  //bias     = myEmit->bias;
  gridSize = myEmit->gridSize;
  jump      = -1.0;
  sigma     = -1.0;
  rnd_emit  = -1.0;
  wTotal    = 0;
  pstay     = new double * [nSamples];
  pjump     = new double * [nSamples];
  pnojump   = new double * [nSamples];
  // matrices...
  alpha  = new gsl_matrix * [nSamples];
  gamma  = new gsl_matrix * [nSamples];  
  proposal = gsl_vector_alloc(gridSize+1);
  gsl_vector_set_all( proposal, 1.0);// uniform proposal distribution on [0,1]
  for (int s=0; s<nSamples; s++){
    pstay[s]   = new double [nSites[s]];
    pjump[s]   = new double [nSites[s]];
    pnojump[s] = new double [nSites[s]];
    alpha[s] = NULL;
    gamma[s] = NULL;
  }
  pstay_set  = 0;
  save_alpha = 0;
  DiffProp = NULL;
  DiffProp_set = 0;
}

JumpDiffusion::~JumpDiffusion(){
  for (int s=0; s<nSamples; s++){
    delete [] pstay[s];
    delete [] pjump[s];
    delete [] pnojump[s];
    if (alpha[s] != NULL) gsl_matrix_free(alpha[s]);
    if (gamma[s] != NULL) gsl_matrix_free(gamma[s]);
  }
  delete [] alpha;
  delete [] pstay;
  delete [] pjump;
  delete [] pnojump;
  delete [] gamma;
  JumpDiffusion::reset_DiffProp();
}

double JumpDiffusion::get_total_llh(){
  save_alpha  = 0;
  if (myEmit->EmitProb_set == 0){
    myEmit->set_EmitProb(time);
  }
  gsl_vector_set_all(proposal, 1.0/(myEmit->xmax - myEmit->xmin));
  //compute the staying probabilities..
  if (pstay_set == 0)    JumpDiffusion::set_pstay();
  if (DiffProp_set == 0) JumpDiffusion::get_DiffProp();
  int sample;
  total_llh   = 0.0;
  // SAMPLES:
#ifdef _OPENMP
  int nt = min( nSamples, omp_get_max_threads());
#pragma omp parallel for schedule( dynamic, 1) default(shared) num_threads(nt)
#endif
  for (sample=0; sample<nSamples; sample++){
    double llh = JumpDiffusion::do_Fwd(sample);
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      total_llh += llh;
    }
  }
  return(total_llh);
}

void JumpDiffusion::get_posterior(int sample){
  if (alpha[sample] != NULL) gsl_matrix_free(alpha[sample]);
  if (gamma[sample] != NULL) gsl_matrix_free(gamma[sample]);
  alpha[sample] = gsl_matrix_alloc(nSites[sample],gridSize+1);
  gamma[sample] = gsl_matrix_alloc(nSites[sample],gridSize+1);
  save_alpha  = 1;
  if (myEmit->EmitProb_set == 0){
    myEmit->set_EmitProb(time);
  }
  gsl_vector_set_all(proposal, 1.0/(myEmit->xmax - myEmit->xmin));
  if (pstay_set == 0)    JumpDiffusion::set_pstay();
  if (DiffProp_set == 0) JumpDiffusion::get_DiffProp();
  JumpDiffusion::do_Fwd(sample);
  JumpDiffusion::do_Bwd(sample);
  gsl_matrix_free(alpha[sample]);
  alpha[sample] = NULL;
}

void JumpDiffusion::reset_DiffProp(){
  if (DiffProp != NULL){
    for (int i=0; i< (int) myEmit->frequent_dist.size(); i++){
      gsl_matrix_free(DiffProp[i]);
    }
    delete [] DiffProp;
    DiffProp = NULL;
  }
}

void JumpDiffusion::get_DiffProp(){
  if (myEmit->dist_set == 0){
    cout<<"ERROR-1 in JumpDiffusion::set_DiffProp()\n";
    exit(1);
  }
  // allocate...
  map<unsigned int,int>::iterator it;
  if (DiffProp == NULL){
    DiffProp = new gsl_matrix * [myEmit->frequent_dist.size()];
    int i=0;
    for (it = myEmit->frequent_dist.begin(); it != myEmit->frequent_dist.end(); ++it){
      DiffProp[i] = gsl_matrix_alloc(gridSize+1,gridSize+1);
      position.insert(pair<unsigned int,int>(it->first,i));
      i++;
    }
  }
  // now set the matrices...
  is_identity.clear();
  int i=0;
  for (it=myEmit->frequent_dist.begin(); it != myEmit->frequent_dist.end(); it++){
    int is_id = JumpDiffusion::set_DiffProp( DiffProp[i],  sigma*sqrt(double(it->first)));
    is_identity.push_back(is_id);
    i++;
  }
  DiffProp_set = 1;
}


int JumpDiffusion::set_DiffProp( gsl_matrix * propagator, double sd){
  double dx = myEmit->dx;
  if (3.0*sd <= dx) return(1);// with the current resolution, propagator is Dirac-Delta!
  int range = 3 * ceil(sd / dx);// this means, only fluctuations up to 3 sigma are possible
  if (2*range > gridSize+1){
    sd = dx * double(gridSize) / 6.0;
    range = 3 * ceil(sd / dx);
  }
  gsl_vector * gauss = gsl_vector_alloc(2*range+1);
  for (int i=0; i<2*range+1; i++){
    gsl_vector_set(gauss,i, gsl_ran_gaussian_pdf( double(i-range)*dx, sd));
  }
  gsl_matrix_set_zero(propagator);
  double val = 0, norm=0;
  gsl_vector_view row;
  for (int i=0; i<= gridSize; i++){
    norm = 0.0;
    for (int j=i-range; j<=i+range; j++){
      if ( j<0 || j>gridSize) continue;
      val = gsl_vector_get( gauss, j-i+range);
      if( i <= range && j < range - i + 1){
	val += gsl_vector_get( gauss, range - i - j);
      }
      if( i >= gridSize - range && j >= 2*gridSize - i - range){
	val += gsl_vector_get( gauss, 2*gridSize + range - i - j);
      }
      gsl_matrix_set(propagator,i,j,val);
      norm += (j==0 || j==gridSize) ? 0.5*val : val;
    }
    row = gsl_matrix_row(propagator,i);
    norm *= dx;
    if (norm <=0.0 ||norm!=norm){
      cout<<"ERROR in JumpDiffusion::set_DiffProp()\n";
      printf("sd=%e dx=%e\n",sd,dx);
      exit(1);
    }
    gsl_vector_scale(&row.vector,1.0/norm);
  }
  gsl_vector_free(gauss);
  return(0);
}


//compute the staying probabilities..
void JumpDiffusion::set_pstay(){
  int s;
  // SAMPLES:
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic, 1) default(shared)
#endif
  for (s=0; s<nSamples; s++){
    pstay[s][0] = 1.0;
    for (int l=1; l< nSites[s]; l++){
      if (dist[s][l] > 0){
	if (jump == 1.0){
	  pstay[s][l] = 0.0;
	}
	else if (jump > 0.0){
	  pstay[s][l] = pow( 1.0-jump, dist[s][l]);
	}
	else if (jump == 0.0){
	  pstay[s][l] = 1.0;
	}
      }
      else{
	pstay[s][l] = 1.0;
      }
    }
  }
  pstay_set = 1;
}

int JumpDiffusion::predict(gsl_vector * prior, gsl_vector * post, 
			   gsl_matrix*& DP, gsl_matrix**& DP_pt, int s, int l){
  int is_id=0;
  //check whether DiffProp is pre-computed
  if (  (myEmit->frequent_dist).count(dist[s][l]) == 0 ){//not found
    is_id = JumpDiffusion::set_DiffProp( DP, sigma*sqrt(double(dist[s][l])));
    DP_pt = &DP;
  }
  else{//exists
    if ( is_identity[ position[dist[s][l]] ] == 1 ){
      is_id = 1;
      DP_pt = &DP;
    }
    else{
      is_id = 0;
      DP_pt = &( DiffProp[ position[dist[s][l]] ] );
    }
  }
  //apply diffusion propagator
  if (is_id == 0){
    gsl_vector * mem = gsl_vector_alloc(gridSize+1);
    post->data[0]        *= 0.5;
    post->data[gridSize] *= 0.5;
    gsl_blas_dgemv(CblasTrans, myEmit->dx, *DP_pt, post, 0.0, mem);
    gsl_vector_memcpy( post, mem);
    gsl_vector_free(mem);  
  }
  else{//increase variance by hand (exact for gaussian distribution)
    double mn  = get_mean( post, myEmit->xmin, myEmit->xmax);
    double var = get_var(  post, myEmit->xmin, myEmit->xmax, mn);
    double e = var / ( var + 100*pow(sigma,2)*double(dist[s][l]));
    for (int i=0; i<=gridSize; i++){
      double p = post->data[i];
      if (p>0.0) post->data[i] = pow(p,e);
    }
    double norm = gsl_blas_dasum(post) - 0.5*(post->data[0] + post->data[gridSize]);
    norm *= myEmit->dx;
    gsl_vector_scale(post,1.0/norm);
  }
  //apply jump propagator
  gsl_vector_scale( prior, 1.0 - pstay[s][l]);
  gsl_vector_scale( post,  pstay[s][l]);
  gsl_vector_add( prior, post);
  return(is_id);
}



double JumpDiffusion::do_Fwd(int s){
  // prepare fwd...
  gsl_vector * eprob = gsl_vector_alloc(gridSize+1);
  gsl_vector * prior = gsl_vector_alloc(gridSize+1);
  gsl_vector * post  = gsl_vector_alloc(gridSize+1);
  //gsl_vector * mem   = gsl_vector_alloc(gridSize+1);
  gsl_matrix * DP    = gsl_matrix_alloc(gridSize+1,gridSize+1);
  gsl_matrix ** DP_pt = NULL;
  double norm;
  double llh=0.0;
  int get_log=0;
  // Forward Pass
  for (int l=0; l<nSites[s]; l++){
    if (mask[s][l] == 0) continue;
    int N = myEmit->depths[time][s][l];
    int n = myEmit->reads[time][s][l];
    gsl_vector_memcpy( prior, proposal);
    if (dist[s][l] > 0) JumpDiffusion::predict( prior, post, DP, DP_pt, s, l);
    //emission probability
    if (N>0){//if observation
      if (myEmit->bias != NULL){
	myEmit->get_eprob_wBias( eprob, myEmit->EmitLog[N][n], myEmit->bias[s][l], n, N, get_log);
      }
      else{
	gsl_vector_memcpy( eprob, myEmit->EmitProb[N][n] );
	if (myEmit->reflect && n != N-n){
	  gsl_vector_add( eprob, myEmit->EmitProb[N][N-n]);
	  //gsl_vector_scale(eprob,0.5);
	}
      }
      //random emission channel
      gsl_vector_scale( eprob, 1.0-rnd_emit);
      if (mode==1 || mode==2){
	gsl_vector_add_constant( eprob, rnd_emit / double(N+1));
      }
      else if (mode==3 || mode==4){
	double rnd = double(N)*(myEmit->maxRate - myEmit->minRate);
	if (myEmit->bias != NULL) rnd *= myEmit->bias[s][l];
	rnd =  rnd_emit / (rnd+1.0);
	gsl_vector_add_constant( eprob, rnd);
      }
      gsl_vector_mul( prior, eprob);// at this time it is the posterior!
      norm  = gsl_blas_dasum(prior);
      norm -= 0.5*(prior->data[0] + prior->data[gridSize]);
      norm *= myEmit->dx;
      if (norm <=0 || norm != norm){
	cout<<"ERROR\n";
	abort();
      }
      gsl_vector_scale( prior, 1.0 / norm);
      llh += log(norm);// get part of the total log-likelihood
    }
    gsl_vector_memcpy(post,prior);
    if (save_alpha == 1){
      gsl_matrix_set_row(alpha[s], l, post);// save forward variable    
    }
  }
  // clean up...
  gsl_vector_free(eprob);
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_matrix_free(DP);
  return(llh);
}


// Backward Pass...
void JumpDiffusion::do_Bwd(int s){
  // prepare bwd...
  gsl_vector * eprob = gsl_vector_alloc(gridSize+1);
  gsl_vector * prior = gsl_vector_alloc(gridSize+1);
  gsl_vector * post  = gsl_vector_alloc(gridSize+1);
  gsl_vector * beta  = gsl_vector_alloc(gridSize+1);
  gsl_vector * last_beta  = gsl_vector_alloc(gridSize+1);
  gsl_vector * last_eprob = gsl_vector_alloc(gridSize+1);
  gsl_vector * mem   = gsl_vector_alloc(gridSize+1);
  gsl_vector * mem2  = gsl_vector_alloc(gridSize+1);
  gsl_matrix * DP    = gsl_matrix_calloc(gridSize+1,gridSize+1);
  gsl_matrix ** DP_pt = NULL;
  double x,y,norm;
  int get_log=0;
  int is_id=1;
  int last=-1;
  for (int l = nSites[s]-1; l>=0; l--){
    if (mask[s][l]==0) continue;
    gsl_vector_memcpy(prior,proposal);
    if (last>0) is_id = JumpDiffusion::predict( prior, post, DP, DP_pt, s, last);
    gsl_vector_memcpy( beta, prior);
    // get gamma, i.e. the total posterior probability vector
    gsl_vector_view alph = gsl_matrix_row(alpha[s],l);
    gsl_vector_memcpy( post, beta);
    gsl_vector_mul( post, &alph.vector);
    norm  = gsl_blas_dasum(post);
    norm -= 0.5*(post->data[0] + post->data[gridSize]);
    norm *= myEmit->dx;
    gsl_vector_scale( post, 1.0 / norm);
    // posterior on-site sojourn probability.
    gsl_matrix_set_row( gamma[s], l, post);
    // emission probability
    int N = myEmit->depths[time][s][l];
    int n = myEmit->reads[time][s][l];
    if (N>0){
      if (myEmit->bias != NULL){
	myEmit->get_eprob_wBias( eprob, myEmit->EmitLog[N][n], myEmit->bias[s][l], n, N, get_log);
      }
      else{
	gsl_vector_memcpy( eprob, myEmit->EmitProb[N][n] );
	if (myEmit->reflect && n != N-n){
	  gsl_vector_add( eprob, myEmit->EmitProb[N][N-n]);
	}
      }
      // random emission channel
      gsl_vector_scale( eprob, 1.0-rnd_emit);
      if (mode==1 || mode==2){
	gsl_vector_add_constant( eprob, rnd_emit / double(N+1));
      }
      else if (mode==3 || mode==4){
	double rnd = double(N)*(myEmit->maxRate-myEmit->minRate);
	if (myEmit->bias != NULL) rnd *= myEmit->bias[s][l];
	rnd =  rnd_emit / (rnd+1.0);
	gsl_vector_add_constant( eprob, rnd);
      }
      // get posterior update for the next step...
      gsl_vector_mul( prior, eprob);// now it is the posterior!
      norm = gsl_blas_dasum(prior);
      norm -= 0.5*(prior->data[0] + prior->data[gridSize]);
      norm *= myEmit->dx;
      gsl_vector_scale(prior, 1.0 / norm);
    }
    gsl_vector_memcpy( post, prior);
    // posterior jump-probability...
    if (last>0){
      if (jump > 0.0 && jump < 1.0){
	gsl_vector_memcpy(mem,last_beta);
	gsl_vector_mul(mem,last_eprob);
	mem->data[0]        *= 0.5;
	mem->data[gridSize] *= 0.5;
	x = gsl_blas_dasum(mem);
	x *= myEmit->dx * (1.0-pstay[s][last]) / (myEmit->xmax - myEmit->xmin);
	if (is_id==0){
	  gsl_blas_dgemv(CblasNoTrans, myEmit->dx, *DP_pt, mem, 0.0, mem2);
	}
	else{
	  gsl_vector_memcpy(mem2,mem);
	}
	gsl_vector_mul(mem2,&alph.vector);
	mem2->data[0]        *= 0.5;
	mem2->data[gridSize] *= 0.5;
	y = gsl_blas_dasum(mem2);
	y *= myEmit->dx * pstay[s][last];
	// log(pjump) and log(1.0-pjump) for the transition l->l+1
	pjump[s][last]   = log(x) - log(x+y);
	pnojump[s][last] = log(y) - log(x+y);
      }
      else if (jump == 0.0){
	pjump[s][last]   = -1.0e3;
	pnojump[s][last] = 0.0;
      }
      else if (jump == 1.0){
	pjump[s][last]   = 0.0;
	pnojump[s][last] = -1.0e3;
      }
    }
    gsl_vector_memcpy(last_beta,beta);
    gsl_vector_memcpy(last_eprob,eprob);
    last = l;
  }
  pjump[s][0]   = -1.0e3;
  pnojump[s][0] = 0.0;
  //clean up...
  gsl_vector_free(eprob);
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(beta);
  gsl_vector_free(mem);
  gsl_vector_free(mem2);
  gsl_vector_free(last_beta);
  gsl_vector_free(last_eprob);
  gsl_matrix_free(DP);
}



int JumpDiffusion::adapt_range(){
  gsl_vector * cum = gsl_vector_calloc(gridSize+1);
  //gsl_vector_view gma;
  double mn = myEmit->xmax;
  double mx = myEmit->xmin;
  int low=0;
  int redo=0;
  double eps = 1.0e-4;
  for (int s=0; s<nSamples; s++){
    JumpDiffusion::get_posterior(s);
    for (int l=0; l<nSites[s]; l++){
      if (mask[s][l] == 0) continue;
      low=0;
      cum->data[0] = 0.0;
      for (int i=1; i<=gridSize;i++){
	cum->data[i]  = cum->data[i-1];
	cum->data[i] += 0.5*myEmit->dx*(gsl_matrix_get(gamma[s],l,i-1) + gsl_matrix_get(gamma[s],l,i));
	if (low==0 && cum->data[i] >= eps){
	  mn = min( mn, myEmit->xgrid[i-1]);
	  low=1;
	}
	if (cum->data[i] >= 1.0-eps){
	  mx = max( mx, myEmit->xgrid[i]);
	  break;	  
	}		
      }
    }
    gsl_matrix_free(gamma[s]);
    gamma[s] = NULL;
  }
  gsl_vector_free(cum);
  if (mx-mn < 10.0){
    double xc = 0.5*(mn+mx);
    mn = max(0.01,xc - 5.0);
    mx = xc + 5.0;
  }
  if (mn > myEmit->xmin || mx < myEmit->xmax){
    redo = 1;
    myEmit->xmin = mn;
    myEmit->xmax = mx;
    if (myEmit->bias == NULL){
      myEmit->ymin = mn;
      myEmit->ymax = mx;
    }
    myEmit->set_grid();
    myEmit->set_EmitProb(time);
    Fwd_done=0;
    Bwd_done=0;
  }
  return(redo);
}
