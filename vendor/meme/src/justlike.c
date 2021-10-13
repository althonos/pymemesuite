/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/
/* 7-13-99 tlb; remove dz */

/* justlike.c */
/*	
	Set the z-matrix to 1 for each motif start; to 0 otherwise.
	Motifs have been read in from a .motif file.
*/

#include "meme.h"

/**********************************************************************/
/*
	like_e_step

*/
/**********************************************************************/
double like_e_step(
  MODEL *model,					/* the model */
  DATASET *dataset 				/* the dataset */
)
{
  int i, j;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;
  int nsites = 0;				/* number of sites read in */
  int imotif = model->imotif;			/* motif number */
  MOTIF motif = dataset->motifs[imotif-1];	/* known motif */
  int w = motif.width;				/* width of known motif */

  /* set all z's to 0 except motif starts */
  for (i=0; i < n_samples; i++) {		/* sequence */
    SAMPLE *s = samples[i];
    Z_T *zi = s->z;                          // zi[j], j in [-lseq...+lseq] 
    char *sample_name = s->sample_name;
    int lseq = s->length;
    for (j=0; j<=lseq-w; j++) {			/* offset */
      if (hash_lookup(sample_name, j+1, motif.ht) != NULL) {
        Zi(j) = 1;				/* motif start */
        nsites++;
      } else {
        Zi(j) = 0;				/* not motif start */
      }
    }
  }

  /* calculate lambda for motif */
  model->lambda = (double) nsites/wps(dataset, w, 1); 

  return 0; 
} /* like_e_step */

