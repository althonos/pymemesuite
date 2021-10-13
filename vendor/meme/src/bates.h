#ifndef BATES_H
#define BATES_H
double get_dtc_log_pv(
  double dtc,                           // (weighted) average distance to center
  double wN,                            // (weighted) number of sites
  int w,                                // motif width
  int slen                              // sequence length
);
#endif
