#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include "log_erfc.h"
/***********************************************************************/
/*
        get_dtc_log_pv

	Get p-value of the mean distance between the site center
	and sequence center using the Bates distribution.
*/
/***********************************************************************/
double get_dtc_log_pv(
  double dtc,                           // (weighted) average distance to center
  double wN,                            // (weighted) number of sites
  int w,                                // motif width
  int slen                              // sequence length
) {
  bool has_zero = (w+slen) % 2 ? false : true;       // 0 dist possible if motif & sequence both odd or even length
  int a = has_zero ? 0 : 1;     // minimum distance between site and sequence centers; rounded up if not has_zero
  int b = (slen - w)/2.0 + 0.5; // maximum distance between site and sequence centers
  double mu = (a + b)/2.0;                                      // mean of Bates distribution
  double sd = wN ? sqrt( ((b-a)*(b-a)) / (12*wN) ) : 100;       // standard deviation of Bates distribution
  double z = (dtc - mu) / sd;           // Z-score of dtc
  double log_pv = log(0.5) + log_erfc(-z/sqrt(2.0));    // Pr(avg_dist <= dtc)
  return(log_pv);
} // get_dtc_log_pv
