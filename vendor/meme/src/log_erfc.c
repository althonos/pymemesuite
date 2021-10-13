/************************************************************************
*       Copyright                                                       *
*       (1999-2016) The Regents of the University of California.        *
*       All Rights Reserved.                                            *
*       Author: Timothy L. Bailey
************************************************************************/
#include <math.h>

//
// from gsl-2.4/specfunc/erfc.c
//
static double erfc8_sum(double x)
{
  /* estimates erfc(x) valid for 8 < x < 100 */
  /* This is based on index 5725 in Hart et al */

  static double P[] = {
      2.97886562639399288862,
      7.409740605964741794425,
      6.1602098531096305440906,
      5.019049726784267463450058,
      1.275366644729965952479585264,
      0.5641895835477550741253201704
  };
  static double Q[] = {
      3.3690752069827527677,
      9.608965327192787870698,
      17.08144074746600431571095,
      12.0489519278551290360340491,
      9.396034016235054150430579648,
      2.260528520767326969591866945,
      1.0
  };
  double num=0.0, den=0.0;
  int i;

  num = P[5];
  for (i=4; i>=0; --i) {
      num = x*num + P[i];
  }
  den = Q[6];
  for (i=5; i>=0; --i) {
      den = x*den + Q[i];
  }

  return num/den;
}

inline static double log_erfc8(double x)
{
  double e;
  e = erfc8_sum(x);
  e = log(e) - x*x;
  return e;
}

//
// Get the log of the complementary error function.
//
double log_erfc(double x)
{
  return (x > 8.0) ? log_erfc8(x) : log(erfc(x));
}
