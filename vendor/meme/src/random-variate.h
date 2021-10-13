/********************************************************************
 * FILE: random-variate.c
 * AUTHOR: Timothy Bailey
 * CREATE DATE: 26/07/2011
 * PROJECT: MEME suite
 * COPYRIGHT: 2019, UNR
 ********************************************************************/
#ifndef RANDOM_VARIATE_H
#define RANDOM_VARIATE_H

typedef struct rv_state {
  double mu;
  double sd;
  double *cdf;
  int min_x;
  int max_x;
} RV_STATE_T;

RV_STATE_T *init_poisson(
  double mu             // mean of Poisson
);

/*************************************************************************
 * Generate a sample from a Poisson variate.
 *************************************************************************/
int rand_poisson(
  RV_STATE_T *state                     // created by init_poisson()
);

/*************************************************************************
 * Generate a sample from a Gaussian variate using the
 * Marsaglia polar method.
 *************************************************************************/
double rand_gaussian(
  double mu,
  double sigma
);

#endif
