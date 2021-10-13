/*
 * $Id:$
 * 
 */
#ifndef PSP_H
#define PSP_H
#include "meme.h"

extern void read_psp_file(
  char *psp_filename,		/* file containing positional priors */
  DATASET *dataset,		// the dataset
  bool psp_revcomp,          // PSP file contains both strands
  bool meme_revcomp,         // MEME using both strands
  MOTYPE mtype			/* OOPS, ZOOPS or TCM? */
);
extern void psp_renormalize(
  DATASET *dataset, 		/* the dataset */
  int new_w,			/* new motif width */
  bool invcomp,		/* reverse complement? */
  MOTYPE mtype			/* OOPS, ZOOPS or TCM? */
);
extern void add_psp_to_log_not_o(
  DATASET *dataset,             /* the dataset */
  int w,                        /* motif width */
  bool invcomp,              /* reverse complement? */
  MOTYPE mtype                   /* model type */
);

#endif
