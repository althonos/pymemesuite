#ifndef SPAMO_OUTPUT_H
#define SPAMO_OUTPUT_H

#include <stdio.h>

#include "motif.h"
#include "spamo-matches.h"

/**************************************************************************
 * Outputs TSV for the histogram of the spacings
 **************************************************************************/
void output_secondary_motif_tsv(FILE *tsv_output, MOTIF_DB_T *primary_db,
    MOTIF_T *primary_motif, SECONDARY_MOTIF_T *parent, SECONDARY_MOTIF_T *smotif,
    int n_secondary_motifs, LINKLST_T *rmotifs);

/**************************************************************************
 * Create an output file and dump the sequence matches to file.
 **************************************************************************/
void output_sequence_matches(char *dir, int margin, int bin, double sigthresh,
    bool sig_only, RBTREE_T *sequences, 
    //MOTIF_T *primary_motif, SECONDARY_MOTIF_T *secondary_motif, int *matches);
    MOTIF_T *primary_motif, SECONDARY_MOTIF_T *secondary_motif, ARRAY_T **matches);

#endif
