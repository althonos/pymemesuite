#ifndef PSSM_H
#define PSSM_H

// referred to by mhmm-state.h
typedef struct pssm PSSM_T;

#include "mhmm-state.h"
#include "prior-dist.h"

//
// Range of integral score values for a PSSM column.
//
// Increased from 100 for more accuracy in descaled log-odds scores. (TLB 31 May 2013)
#define PSSM_RANGE 1000

//
// Macros to convert between scaled score and raw score.
//
#define scaled_to_raw(x,w,scale,offset) (((x)/(scale)) + ((w)*(offset)))
#define raw_to_scaled(x,w,scale,offset) (nint(((x) - ((w)*offset)) * (scale)))

//
// Macros to encapsulate PSSM_T structure
#define pssm_scaled_to_raw(x, pssm) (scaled_to_raw((x), (pssm)->w, (pssm)->scale, (pssm)->offset))
#define get_pssm_w(pssm) ((pssm)->w)
#define get_pssm_alph(pssm) ((pssm)->alph)
#define get_pssm_alphsize(pssm) ((pssm)->alphsize)
#define get_pssm_scale(pssm) ((pssm)->scale)
#define get_pssm_offset(pssm) ((pssm)->offset)
#define get_pssm_pv_length(pssm) (get_array_length((pssm)->pv))
#define get_pssm_pv(score, pssm) (get_array_item((score), (pssm)->pv))
#define get_pssm_score(row, col, pssm) (get_matrix_cell((row), (col), (pssm)->matrix))
#define get_pssm_min_pv(pssm) (get_array_item((pssm)->max_score, (pssm)->pv))
#define get_pssm_max_raw_score(pssm) (pssm_scaled_to_raw((pssm)->max_score, (pssm)))

//
// PSSM
//
// This object was created for AMA because the "scale" and "offset"
// parameters need to be stored with each PSSM, but were 
// globals in pssm.c.  This object should be used in all the programs
// that use PSSMs, since their scale and offset can differ, and their
// cdfs and pdfs should be kept with them.
//
struct pssm {
  MATRIX_T *matrix;	// The PSSM score matrix.
  ALPH_T *alph;		// The alphabet of the pssm
  int alphsize;         // The size of the alphabet after hashing
  int w;		// Width of PSSM.
  bool matrix_is_log;	// True if matrix is log likelihood ratio.
  bool matrix_is_scaled;	// True if matrix is scaled.
  double scale;		// Scale factor for scores.
  double offset;	// Offset for scores.
  int range;		// Scaled scores in range [0..range].
  ARRAY_T *pv;		// P-value table for scores.
  int num_gc_bins;	// Number of entries in gc_pv list.  If > 1, then ->pv is NULL.
  ARRAY_T **gc_pv;	// P-value tables for different GC contents: [gc_bin, score].
  int min_score;	// Smallest index with non-zero pdf.
  int max_score;	// Largest index with non-zero pdf.
  double nolog_max_score;
  MOTIF_T *motif;	// may be NULL but can be useful e.g. for id
};

//
// PSSM_PAIR
//
// PSSMs for the negative and positive DNA motifs.
//
typedef struct pssm_pair {
  PSSM_T* pos_pssm;		// positive strand PSSM
  PSSM_T* neg_pssm;		// negative strand PSSM 
  // Stuff below here is for AMA:
  // The pv lookup table for the average of n scores will be
  // in row log_2(n), for n=1, 2, 4, ...
  int num_gc_bins;		// this is the number of n_pv_lookup tables
  MATRIX_T** gc_n_pv_lookup;	// pv[gcbin, log_2(n), score] lookup table
  ARRAY_T* scaled_to_ama;	// for speed
  MOTIF_T* motif;               // use with care in case motif deallocated
} PSSM_PAIR_T;

void set_up_pssms_and_pvalues (
  bool motif_scoring, // Motif scoring?
  double p_threshold, // Scale/offset PSSM and create table if > 0
  bool use_both_strands, // Compute PSSM for negative strand, too?
  bool allow_weak_motifs, // Allow motifs with min p-value < p_threshold?
  MHMM_T *the_hmm, // The Hidden Markov Model
  PRIOR_DIST_T *prior_dist, // Distribution of priors
  double alpha // PSP alpha parameter
);

void compute_motif_score_matrix(
  bool use_pvalues, // Returns scores as p-values, not log-odds.
  double p_threshold,	// Divide p-values by this.
  int* int_sequence,	// Sequence as integer indices into alphabet.
  size_t seq_length,	// Length of sequence.
  double *priors, // Priors for sequence
  size_t num_priors, // Number of priors
  double alpha, // Alpha parameter for calculating prior log odds
  MHMM_T*    the_hmm,		// The hmm.
  MATRIX_T** motif_score_matrix
);

void scale_score_matrix(
  MATRIX_T *matrix, // The matrix (IN/OUT)
  int rows, int cols, // Rows and columns to access in matrix (IN)
  PRIOR_DIST_T *prior_dist, // Distribution of priors (OPTIONAL IN)
  double alpha, // Fraction of all TFBS that are in the TFBS of interest (IN)
  int range, // The desired range (IN)
  double *offset_p, // Set to calculated range (OPTIONAL OUT)
  double *scale_p // Set to calculated scale (OPTIONAL OUT)
);

void scale_pssm(
  PSSM_T *pssm, // The PSSM. (IN/OUT)
  PRIOR_DIST_T *prior_dist, // Distribution of priors (IN)
  double alpha, // Fraction of all TFBS that are the TFBS of interest (IN)
  int range // The desired range. (IN) 
);

ARRAY_T *scale_prior_dist(
  ARRAY_T *priors, // Distribution of priors (IN/OUT)
  int range,	   // The desired range. (IN) 
  double scale,    // The desired scale. (IN)
  double offset    // The desired offset. (IN)
);

void get_pv_lookup_pos_dep(
  PSSM_T* pssm,			// The PSSM.
  MATRIX_T* background_matrix,  // The background model PSSM matrix.
  ARRAY_T* scaled_prior_dist    // Scaled distribution of priors.
);

void get_pv_lookup(
  PSSM_T* pssm,			// The PSSM.
  ARRAY_T* background, 		// The background model.
  ARRAY_T* scaled_prior_dist	// Scaled distribution of priors.
);

double get_unscaled_pssm_score(
  double score,
  PSSM_T* pssm
);

double get_scaled_pssm_score(
  double score,
  PSSM_T* pssm
);

PSSM_T* build_motif_pssm(
  MOTIF_T* motif,	           // motif frequencies p_ia (IN)
  ARRAY_T* bg_freqs,	       // background frequencies b_a for pssm (IN)
  ARRAY_T* pv_bg_freqs,	     // background frequencies b_a for p-values (IN)
  PRIOR_DIST_T* prior_dist,  // Distribution of priors. May be NULL (IN)
  double alpha,              // Scale factor for non-specific priors. 
                             // Unused if prior_dist is NULL.
  int range,		             // range of scaled scores is [0..w*range]
  int num_gc_bins,	         // create pv tables for this number of GC bins
			                       // instead of using the pv_bg_freqs
  bool no_log  	       // make likelihood ratio pssm
);

PSSM_T* build_matrix_pssm(
  ALPH_T* alph, // alphabet (IN)
  MATRIX_T* matrix, // pssm matrix (IN)
  ARRAY_T* bg_freqs, // background frequencies b_a (IN)
  PRIOR_DIST_T *prior_dist, // Prior distribution
  double alpha, // PSP alpha parameter
  int range // range of scaled scores is [0..w*range] (IN)
);

double get_ama_pv(
  double ama_score,                     // average likelihood ratio score
  int seqlen,                           // length of sequence scanned
  double seq_gc,                        // total GC content of sequence
  PSSM_PAIR_T* pssm_pair                // pssms for pos and neg motifs
);

PSSM_PAIR_T* create_pssm_pair(
  PSSM_T* pos_pssm,		// positive strand pssm
  PSSM_T* neg_pssm 		// negative strand pssm
);

void free_pssm_pair(
  PSSM_PAIR_T *pssm_pair
);

PSSM_T* allocate_pssm(
  ALPH_T* alph,
  int w, 
  int alphsize, // may be different from alph_size(alph, ALPH_SIZE) if PSSM is hashed
  int num_gc_bins
);

void free_pssm(
  PSSM_T* pssm
);

/**************************************************************************
 * sum up the largest values for each position in the pssm to generate 
 * the best possible score for a match. Output is unscaled. 
 **************************************************************************/
double pssm_best_match_score(
  PSSM_T* pssm
);

// Print a PSSM to standard error.
void print_pssm(
  PSSM_T* pssm
);

// Find the unscaled (raw) score corresponding most closely to a given p-value.
double pssm_pvalue_to_score(
  PSSM_T *pssm,                 // a pssm
  double pvalue                 // a p-value
);

#endif
