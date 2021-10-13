#ifndef DTC_H
#define DTC_H

#include "meme.h"

void init_dtc_pv_tables(
  int minw,             // minimum motif width
  int maxw,             // maximum motif width
  int slen,             // sequence length
  int maxn              // maximum number of seqs
);

double get_dtc_pv(
  double dtc,                   // average, weighted distance to center
  double n,                     // (weighted) number of sequences
  int maxd,                     // maximum possible distance to center
  bool has_zero              // true if zero distance is possible
);

#endif
