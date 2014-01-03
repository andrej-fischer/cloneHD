//log_space.cpp

#define PI 3.1415926

//own headers...
#include "log-space.h"

// adding two vectors in log space...
void log_vector_add(gsl_vector * one, gsl_vector * two){
  for (int i=0; i<(int) one->size; i++){
    gsl_vector_set( one, i, log_add( one->data[i], two->data[i]));
  }
}

// adding two matrices in log space...
void log_matrix_add(gsl_matrix * one, gsl_matrix * two){
  for (int i=0; i<(int) one->size1; i++){
    gsl_vector_view r1=gsl_matrix_row(one,i);
    gsl_vector_view r2=gsl_matrix_row(two,i);
    log_vector_add(&r1.vector,&r2.vector);
  }
}

// adding two scalars in log space...
double log_add(double one, double two){
  if (one > two){
    if (one-two > 10.0){
      return( one + exp(two - one));
    }
    else{
      return( one + log(1.0 + exp(two - one)));
    }
  }
  else if (one < two){
    if (two-one > 10.0){
      return( two + exp(one - two)); 
    }
    else{
      return( two + log(1.0 + exp(one - two))); 
    }
  }
  else{
    return(one + log(2.0));
  }
}

double log_sub(double one, double two){
  if (two > one){
    printf("ERROR in log_sub(%e,%e)\n",one,two);
  }
  if (one-two > 10.0){
    return( one - exp(two - one));
  }
  else{
    return( one + log(1.0 - exp(two - one)));
  }
}


void log_vector_invert(gsl_vector * vec){
  double val;
  for (int i=0; i<(int) vec->size; i++){
    val = vec->data[i];
    vec->data[i] = log(1.0 - exp(val));
  }
}


// 1-norm of a vector in log-space
double log_vector_norm(const gsl_vector * x){
  double norm = 0;
  double max  = gsl_vector_max(x);
  double crit = max - 10.0;
  for (int i=0; i<(int) x->size; i++){
    if ( x->data[i] >  crit) 
      norm += exp(x->data[i] - max);
  }
  if (norm>1.1){
    norm = log(norm) + max;
  }
  else{//Taylor-expansion of log(1+x) = x - 0.5*x*x
    norm = norm - 1.0;
    norm = norm*(1.0 - 0.5*norm) + max;
  }
  return(norm);
}

void log_vector_normalize( gsl_vector * x){
  double norm = log_vector_norm(x);
  gsl_vector_add_constant(x,-norm);
}

// 1-norm of a vector in log-space
double log_matrix_norm(const gsl_matrix * M){
  double norm=0;
  double max = gsl_matrix_max(M);
  for (int i=0; i<(int) M->size1; i++){
    for (int j=0; j<(int) M->size2; j++){
      norm += exp(gsl_matrix_get(M,i,j) - max);
    }
  }
  norm = log(norm) + max;
  return(norm);
}
