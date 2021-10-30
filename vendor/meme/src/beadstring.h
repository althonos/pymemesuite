#ifndef BEADSTRING_H
#define BEADSTRING_H

#include "utils.h"          // Generic utilities. 
#include "alphabet.h"
#include "red-black-tree.h"

// Structure for tracking beadstring command line parameters.
typedef struct options {

  char* command_line;     /* Full command line */
  char* output_dirname;   /* Name of the output directory */
  char* meme_filename;    /* Name of file containg motifs. */
  char* model_filename;   /* Name of external file for HMM. */
  char* seq_filename;     /* Name of file containg sequences. */
  char* bg_filename;      /* Name of file file containg background freq. */
  char* score_filename;
  char* order_string;     /* Motif order and spacing. */

  bool allow_weak_motifs;
  bool allow_clobber;
  bool fim;              /* Represent spacers as free insertion modules? */
  bool local_scoring;
  bool motif_scoring;
  bool perform_search;
  bool synth;
  bool zselo;

  int       max_seq;
  int       pam_distance;
  int       paths;
  int       spacer_states;

  PROB_T  motif_e_threshold;   /* e-value threshold for motif inclusion. */
  PROB_T  motif_p_threshold;   /* p-value threshold for motif occurences. */
  PROB_T  e_threshold;   /* e-value threshold for reporting cluster. */
  PROB_T  p_threshold;   /* Turns on p-value scoring and sets zer point. */
  double  gap_extend;    /* penalty for extending a gap */
  double  gap_open;      /* penalty for opening a gap */
  double  motif_pseudo;  /* pseudocount  for motif frequencies */
  double  spacer_pseudo; /* Spacer (self-loop) pseudocount. */
  double  trans_pseudo;  /* Transition pseudocount. */
  double  progress;      /* interval between progress reports */

  ALPH_T *alphabet;    // Alphabet specified by MEME file.
  RBTREE_T *request_nums; // Indices of requested motifs.
  RBTREE_T *request_ids; // IDs of requested motifs
} OPTIONS_T;


#endif
