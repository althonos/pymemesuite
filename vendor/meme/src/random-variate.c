/********************************************************************
 * FILE: random-variate.c
 * AUTHOR: Timothy Bailey
 * CREATE DATE: 26/07/2011
 * PROJECT: MEME suite
 * COPYRIGHT: 2019, UNR
 ********************************************************************/
#include "random-variate.h"
#include "utils.h"

#define CHUNK 1000

/*************************************************************************
 * Initialization for generating samples from Poisson variate.
 *************************************************************************/
RV_STATE_T *init_poisson(
  double mu 		// mean of Poisson
) {
  int x, min_x, max_x;
  double p;
  RV_STATE_T *state = (RV_STATE_T *)mm_malloc(sizeof(RV_STATE_T));

  double *pdf = NULL;
  double *cdf = NULL;
  Resize(pdf, CHUNK, double);
  Resize(cdf, CHUNK, double);

  // Get CDF of Poisson distribution.
  x = 0;
  pdf[0] = cdf[0] = p = exp(-mu);
  while (p > 0 && cdf[x] < 1) {
    x += 1;
    if ((x % CHUNK) == 0) {
      Resize(pdf, x+CHUNK, double);
      Resize(cdf, x+CHUNK, double);
    }
    p *= mu/x;
    pdf[x] = p;
    cdf[x] = cdf[x-1] + p;
    if (cdf[x] > 1) cdf[x] = 1;
  }
  myfree(pdf);

  state->mu = mu;
  state->sd = mu;
  state->cdf = cdf;
  state->min_x = 0;
  state->max_x = x;

  return(state);
} // init_poisson

/******************************************************************
 * f_inverse
 *
 * Find inverse of f() by binary search.  f() is stored in an
 * array indexed from x_min to x_max.
 * f() is an increasing function of x.
 *
******************************************************************/
static int f_inverse(
  double y,			// get x where y=f(x)
  int x_min,
  int x_max,
  double *f			// array of f(x) for x_min <= x <= x_max
) {
  int x_lo = x_min;
  int x_hi = x_max;
  int x_mid = (x_hi+x_lo)/2;

  double y_lo = f[x_lo]; 
  double y_hi = f[x_hi]; 
  double y_mid = f[x_mid]; 

  while (x_hi-x_lo > 1) {
    if (y < y_lo) {
      x_hi = x_lo;
    } else if (y_lo <= y && y <= y_mid) {
      x_hi = x_mid;
    } else if (y_mid <= y && y <= y_hi) {
      x_lo = x_mid;
    } else if (y_hi < y) {
      x_lo = x_hi;
    }
    x_mid = (x_hi + x_lo)/2;
    y_lo = f[x_lo];
    y_hi = f[x_hi];
    y_mid  = f[x_mid];
  }

  return (y <= y_mid ? x_mid : x_hi);
} // f_inverse

/*************************************************************************
 * Generate a sample from a Poisson variate.
 *************************************************************************/
int rand_poisson(
  RV_STATE_T *state			// created by init_poisson()
) {
  double y = drand_mt();		// random uniform in [0.0, 1.0)
  return(f_inverse(y, state->min_x, state->max_x, state->cdf));
} // rand_poisson

/*************************************************************************
 * Generate a sample from a Gaussian variate using the
 * Marsaglia polar method.
 *************************************************************************/
double rand_gaussian(
  double mu,
  double sigma
) {
  double v1, v2, rsq, x;
  rsq = 1;
  while (rsq >= 1 || rsq == 0) {
    // Generate point uniformly within (-1,+1) square and
    // test to see if it is in the unit circle, not at origin.
    v1 = 2 * drand_mt() - 1;
    v2 = 2 * drand_mt() - 1;
    rsq = (v1*v1) + (v2*v2);
  }
  x = v1 * sqrt((-2 * log(rsq)) / rsq);
  return(x*sigma + mu);
} // rand_gaussian
