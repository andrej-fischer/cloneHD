//common-functions.cpp

#include "common-functions.h"
#include "emission.h"

#define PI 3.1415926
#define LOG2 0.693147

using namespace std;


//get general dimensions of a data set for cloneHD...
void get_dims( const char * data_fn, 
	       int& nTimes,
	       vector<int>& chrs,
	       vector<int>& nSites,
	       int keep
	       ){
  ifstream data_ifs;
  string line;
  stringstream line_ss;
  data_ifs.open( data_fn, ios::in);
  if (data_ifs.fail()){
    printf("ERROR: file %s cannot be opened.\n", data_fn);
    exit(1);
  }
  nSites.clear();
  chrs.clear();
  int ct=0,l,r,d;
  int chr=0,old=-1,nT=0;
  while( data_ifs.good()){
    line.clear();
    getline( data_ifs, line);
    if (line.empty()) break;
    if (line[0] == '#') continue;
    line_ss.clear();
    line_ss.str(line);
    //check first entry for nTimes
    if (old == -1 && ct == 0){
      line_ss >> chr >> l; 
      while(line_ss >> r >> d){
	nT++;
      }
      line_ss.clear();
      line_ss.str(line);      
    }
    line_ss >> chr >> l;
    if (chr != old ){//new chromosome encounter     
      if (ct>0){
	nSites.push_back(ct);
	chrs.push_back(old);
      }
      ct=0;
    }
    old=chr;
    r = 0;
    for( int t=0;t<nT; t++){
      line_ss >> r >> d;
      if (r>0) break;
    }
    if (keep || r>0) ct++;
  }
  if (ct>0){
    nSites.push_back(ct);
    chrs.push_back(old);
  }
  nTimes = nT;
  data_ifs.close();
}


// read in data: expects columns to be "chr location (depth reads)^x"
void get_data( const char * data_fn, Emission * myEmit){
  ifstream data_ifs;
  string line;
  stringstream line_ss;
  data_ifs.open( data_fn, ios::in);
  if (data_ifs.fail()){
    printf("ERROR: file %s cannot be opened.\n", data_fn);
    exit(1);
  }
  int ct=0,l;
  int chr=0,old=-1, sample=0;
  int d,r, keep=0, wait=0;
  //now collect all data...
  while( data_ifs.good()){
    line.clear();
    getline( data_ifs, line);
    if (line.empty()) break;
    if (line[0] == '#') continue;
    line_ss.clear();
    line_ss.str(line);
    line_ss >> chr >> l;//chromosome and locus
    if (chr != old){
      if (myEmit->chrs.count(chr) == 0){
	printf("WARNING: chr %2i in file %s will be ignored.\n", chr, data_fn);
	wait = 1;
      }
      else{
	sample = myEmit->idx_of[chr];
	ct  = 0;
	wait = 0;
      }
      old = chr;
    }
    if (wait) continue;
    if (ct >= myEmit->nSites[sample]) continue;
    keep = 0;
    for (int t=0; t<myEmit->nTimes; t++){
      myEmit->loci[sample][ct] = l;//set locus
      line_ss >> r >> d;//get read and depth
      if (d == 0 && r > 0){
	printf("ERROR: depth = 0 in chr %i locus %i\n", chr, l);
	cout<<line<<endl;
	exit(1);
      }
      //set read and depth
      myEmit->reads[t][  myEmit->idx_of[chr] ][ct] = r;
      myEmit->depths[t][ myEmit->idx_of[chr] ][ct] = d;
      if (r>0) keep=1;
    }
    if (keep || myEmit->connect) ct++;
  }  
  data_ifs.close();
  // set the distances between loci
  for (int t=0; t<myEmit->nTimes; t++){
    myEmit->set_dist();
  }
}


void get_bias(const char * bias_fn, Emission * myEmit){
  ifstream ifs;
  string line;
  stringstream line_ss;
  ifs.open( bias_fn, ios::in);
  if (ifs.fail()){
    printf("ERROR: file %s cannot be opened.\n", bias_fn);
    exit(1);
  }
  int  chr = 0, old_chr = -1, idx=0, blocus=0, next_locus=0,l1=0;
  double b=0,b1=0,nu=0;
  while( ifs.good() ){
    line.clear();
    getline( ifs, line);
    if (line.empty()) break;
    if (line[0] == '#') continue;
    line_ss.clear();
    line_ss.str(line);
    line_ss >> chr >> blocus; 
    if (chr != old_chr && (chr > (int) myEmit->maxchr || myEmit->idx_of[chr] < 0)){
      printf("ERROR 1 in get_bias()\n");
      cout<<line<<endl;
      exit(1);
    }
    if (chr != old_chr){ 
      if (old_chr >=0){
	while ( idx < myEmit->nSites[myEmit->idx_of[old_chr]] ){
	  myEmit->bias[myEmit->idx_of[old_chr]][idx] = b;
	  idx++;
	}
      }
      idx=0;
      next_locus = (int) myEmit->loci[ myEmit->idx_of[chr] ][idx];
      old_chr    = chr;
    }
    if (idx >= myEmit->nSites[ myEmit->idx_of[chr] ]) continue;    
    line_ss >> b;
    if (idx==0){
      while ( next_locus < blocus ){//left overhang
	myEmit->bias[myEmit->idx_of[chr]][idx] = b;
	idx++;
	if (idx < myEmit->nSites[myEmit->idx_of[chr]]){
	  next_locus = (int) myEmit->loci[ myEmit->idx_of[chr] ][idx];
	}
	if (idx >= myEmit->nSites[myEmit->idx_of[chr]]) break;    
      }
    }
    if ( blocus <= next_locus ){
      b1 = b;
      l1 = blocus;
    }
    if ( blocus < next_locus ){
      continue;
    }
    else if (blocus==next_locus){
      myEmit->bias[myEmit->idx_of[chr]][idx] = b;
      idx++;
      if (idx < myEmit->nSites[myEmit->idx_of[chr]]){
	next_locus = (int) myEmit->loci[ myEmit->idx_of[chr] ][idx];
      }
    }
    else if (blocus > next_locus){
      while ( next_locus <= blocus ){
	nu = double(next_locus-l1)/double(blocus-l1);
	myEmit->bias[myEmit->idx_of[chr]][idx] = b1*(1.0-nu) + b*nu;
	idx++;
	if (idx < myEmit->nSites[myEmit->idx_of[chr]]){
	  next_locus = (int) myEmit->loci[ myEmit->idx_of[chr] ][idx];
	}
	if (idx >= myEmit->nSites[myEmit->idx_of[chr]]) break;    
      } 
    }
  }
  ifs.close();
}



//***Mean and Variance function***
double get_mean(gsl_vector * dist, double xmin, double xmax){
  double mean=0.0,P1,P2;
  int n = (int) dist->size;
  double dx = (xmax - xmin) / double(n-1);
  for (int i=0; i < n-1; i++){
    P1 = gsl_vector_get(dist,i);
    P2 = gsl_vector_get(dist,i+1);
    mean += 3.0*(P1+P2)*(xmin+double(i)*dx) + (P1+2.0*P2)*dx;
  }
  mean = mean * dx / 6.0;
  return(mean);
}

double get_var(gsl_vector * dist, double xmin, double xmax, double mean){
  double var=0.0, P1, P2,dev;
  int n = (int) dist->size;
  double dx = (xmax - xmin) / double(n-1);
  for (int i=0; i<n-1; i++){
    P1 = gsl_vector_get(dist,i);
    P2 = gsl_vector_get(dist,i+1);
    dev = xmin + double(i)*dx - mean;
    var += (P1+3.0*P2)*dx*dx + 4.0*(P1+2.0*P2)*dev*dx + 6.0*(P1+P2)*pow(dev,2);
  }
  var = var*dx/12.0;
  return(var);
}
