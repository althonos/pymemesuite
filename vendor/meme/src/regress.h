/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1995-2014, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/
 
#ifndef REGRESS_H 
#define REGRESS_H

extern double regress(
  int n,                        /* number of points */
  double *x,                    /* x values */
  double *y,                    /* y values */
  double *m,                    /* slope */
  double *b                     /* y intercept */
);

double w_regress(
  int n,                        /* number of points */
  double *x,                    /* x values */
  double *y,                    /* y values */
  double *w,                    /* weights */
  double *m,                    /* slope */
  double *b                     /* y intercept */
);

/******************************************************************************
*       pearson_correlation
*
*       Returns the sample Pearson correlation coefficient of two sets of points
*       and computes an estimate of its significance (that it is large, one-tailed)
*       using the Fisher transform.  Also computes the regression line
*       and its mean-squared error.
*
******************************************************************************/
double pearson_correlation(
  int n,                        /* number of points */
  double *x,                    /* x values */
  double *y,                    /* y values */
  double *m,                    /* slope */
  double *b,                    /* intercept */
  double *mse,                  /* mean-squared error */
  double *log_pv,               /* log of Fisher transform p-value */
  bool two_tailed,              // p-value is for two-tailed test; 1-tailed otherwise
  bool *valid_list              // true if point i is valid; ignore point i otherwise; or NULL
);

#endif
