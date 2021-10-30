/**************************************************************************
 * FILE: pssm.c
 * AUTHOR: Timothy L. Bailey
 * CREATE DATE: 6/6/02
 * PROJECT: MHMM
 * COPYRIGHT: 1998-2002, WNG, 2001-2002, TLB
 * DESCRIPTION: Core sequence routines for PSSM scoring and statistics.
 **************************************************************************/
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>
#include "macros.h"
#include "motif.h"       // Needed for complement_dna_freqs.
#include "utils.h"       // Generic utilities.
#include "mhmm-state.h"  // HMM.
#include "matrix.h"      // Routines for floating point matrices.
#include "array.h"       // Routines for floating point arrays.
#include "alphabet.h"    // The alphabet.
#include "fasta-io.h"
#include "binary-search.h"
#include "pssm.h"        // PSSM routines.

// 2 gives the best speed.
static int hash_n = 1;  // The number of characters to hash into 1 for efficiency.

#define get_pair_bin_value(n, i) (((double) (i)) / ((n) - 1))
#define get_pair_bin(n, gc) (((n)-1) * (gc))
#define get_pair_bin_lower(n, gc) floor(((n)-1) * (gc))
#define get_pair_bin_delta(n) (1.0/((n) - 1))

/***********************************************************************
 * allocate memory for a pssm
 ***********************************************************************/
PSSM_T* allocate_pssm(ALPH_T* alph, int w, int alphsize, int num_gc_bins){
  PSSM_T* pssm = mm_malloc(sizeof(PSSM_T));

  pssm->matrix = allocate_matrix(w, alphsize);
  pssm->alph = alph_hold(alph);
  pssm->w = w;
  pssm->alphsize = alphsize;
  pssm->matrix_is_log = false;
  pssm->matrix_is_scaled = false;
  pssm->scale = 0;
  pssm->offset = 0;
  pssm->range = -1;
  pssm->pv = NULL;
  pssm->gc_pv = num_gc_bins <= 1 ? NULL : mm_calloc(num_gc_bins, sizeof(ARRAY_T*));
  pssm->num_gc_bins = num_gc_bins <= 1 ? 0 : num_gc_bins;
  pssm->min_score = 0;
  pssm->max_score = 0;
  return pssm;
}

/**************************************************************************
*
*       get_unscaled_pssm_score
*
**************************************************************************/
double get_unscaled_pssm_score(
  double score,
  PSSM_T* pssm
)
{
  return scaled_to_raw(
           score, 
           get_pssm_w(pssm), 
           get_pssm_scale(pssm), 
           get_pssm_offset(pssm)
    );
}

/**************************************************************************
*
*       get_scaled_pssm_score
*
**************************************************************************/
double get_scaled_pssm_score(
  double score,
  PSSM_T* pssm
)
{
  const int w = get_pssm_w(pssm);
  const double scale = get_pssm_scale(pssm);
  const double offset = get_pssm_offset(pssm);
  return raw_to_scaled(score, w, scale, offset);
}

/**************************************************************************
*       get_min_pvalue
*
*       Return the minimum p-value for a given pssm.
*
**************************************************************************/
static double get_min_pvalue(
  PSSM_T *pssm, // The PSSM.
  PRIOR_DIST_T *prior_dist, // Distribution of priors
  double alpha
)
{
  int i, j;
  int max_score;
  int r = pssm->w;
  int c = pssm->alphsize;
  double min_p_value;

  // Get the largest score in each row and sum them.
  assert(pssm->matrix_is_scaled == true);
  max_score = 0;
  for (i=0; i<r; i++) {
    double large = -BIG;
    for (j=0; j<c; j++) {
      double x = get_matrix_cell(i, j, pssm->matrix);
      large = MAX(large, x);
    }
    max_score += large;
  }

  if (prior_dist) {
    // Get raw score before adding log-odds prior
    // We add 1 to the width of the motif to account for the prior
    max_score = scaled_to_raw(max_score, pssm->w, pssm->scale, pssm->offset);
    double max_lo_prior = get_max_lo_prior(prior_dist, alpha);
    max_score += max_lo_prior;
    // Back to scaled score
    max_score = raw_to_scaled(max_score, pssm->w + 1, pssm->scale, pssm->offset);
  }

  min_p_value = get_array_item(max_score, pssm->pv);

  return(min_p_value);
} /* get_min_pvalue */

/**************************************************************************
*
        hash_pssm_matrix_pos

        Recursively create a single position of a hashed PSSM.

*
**************************************************************************/
static void hash_pssm_matrix_pos(
  MATRIX_T *pssm,               // pssm to hash
  MATRIX_T *hashed_pssm,        // hashed pssm
  int  pos,                     // position in pssm
  int  hashed_pos,              // position in hashed pssm
  int  n,                       // number of columns to hash together
  double score,                 // cumulative score; call with 0
  int index                     // cumulative index; call with 0
)
{
  int i;
  int alen = get_num_cols(pssm);        // alphabet length
  int w = get_num_rows(pssm);           // pssm width

  if (n==0) {                           // done, set hashed_pssm entry
    set_matrix_cell(hashed_pos, index, score, hashed_pssm);
  } else {                              // combine next column of pssm
    for (i=0; i<=alen; i++) {           // letters + blank
      // not past right edge of motif and not blank?
      double s = (pos<w && i!=alen) ? get_matrix_cell(pos, i, pssm) : 0;
      hash_pssm_matrix_pos(pssm,
                    hashed_pssm,
                    pos+1,              // position in old pssm
                    hashed_pos,         // position working on
                    n-1,                // positions remaining to hash
                    score+s,            // score so far
                    index*(alen+1)+i);  // hashed alphabet index so far
    } // leter
  }
} // hash_pssm_matrix_pos

/**************************************************************************
*
        hash_pssm

        Hash a PSSM (in place) where each n columns are combined into
        a single column on the alphabet hashed to the nth power.
*
**************************************************************************/
void hash_pssm(
  PSSM_T* pssm,                         // PSSM to hash
  int n                                 // hash n columns to 1
)
{
  int w = get_pssm_w(pssm);                     // width of pssm
  int alen = get_pssm_alphsize(pssm);           // length of alphabet
  int hashed_w = (w+n-1)/n;                     // width of hashed pssm
  int hashed_alen = pow(alen+1, n) + 1;         // length of hashed alphabet

  // Allocate the hashed PSSM.
  MATRIX_T* pssm_matrix = pssm->matrix;
  MATRIX_T* hashed_pssm_matrix = allocate_matrix(hashed_w, hashed_alen);

  int pos, hashed_pos;
  for (pos=hashed_pos=0; pos<w; pos+=n, hashed_pos++) { // position in pssm
    hash_pssm_matrix_pos(pssm_matrix, hashed_pssm_matrix, pos, hashed_pos, n, 0, 0);
  } // position in pssm

  free_matrix(pssm_matrix);
  pssm->matrix = hashed_pssm_matrix;
  pssm->w = hashed_w;
  pssm->alphsize = hashed_alen;

} // hash_pssm

/**************************************************************************
*
*       hash_sequence
*
*       Hash a sequence, compressing hash_n letters into 1.
*
*       Return the newly allocated sequence.
*
 *************************************************************************/
static int* hash_sequence(
  ALPH_T *alph,
  int *int_sequence,                            // Sequence in integer format.
  int seq_length,                               // Length of sequence.
  int hash_n                                    // Number of letters to compress to 1.
)
{
  int i, j;
  int base = alph_size_full(alph) + 1;     // Base to hash to.
  int* hashed_sequence = NULL;

  // Allocate the hashed sequence.
  mm_resize(hashed_sequence, seq_length, int);

  for(i=0; i<seq_length; i++) {
    int c = int_sequence[i];                    // Character in hashed alphabet.
    int* old_cp;                                // Pointer to unhashed character in int_sequence.
    if ((seq_length - i - hash_n) < 0) {        // Hash window is within sequence.
      for(j=1, old_cp=&(int_sequence[i+1]); j<hash_n; j++, old_cp++) {
        c = (base * c) + *old_cp;
      }
    } else {                                    // Hash window runs off sequence end.
      for(j=1, old_cp=&(int_sequence[i+1]); j<hash_n; j++, old_cp++) {
        c = (base * c);
        if (old_cp - int_sequence < seq_length) c += *old_cp;
      }
    }
    hashed_sequence[i] = c;                     // Record the hashed character.
  }

  return(hashed_sequence);
} // hash_sequence

/**************************************************************************
*
*       reverse_complement_pssm_matrix
*
*       Turn a pssm matrix into its own reverse complement.
*
 *************************************************************************/
static void reverse_complement_pssm (
  ALPH_T* alph,
  MATRIX_T* pssm_matrix
)
{
  int i;
  ARRAY_T* left_scores;
  ARRAY_T* right_scores;
  int length = get_num_rows(pssm_matrix);

  // Consider each row (position) in the motif.
  for (i = 0; i < (int)((length+1) / 2); i++) {
    left_scores = get_matrix_row(i, pssm_matrix);
    right_scores = get_matrix_row(length - (i + 1), pssm_matrix);

    // Compute reverse complements in both directions.
    complement_swap_freqs(alph, left_scores, right_scores);
  }
} // reverse_complement_pssm_matrix

/**************************************************************************
*       scale_score_matrix
*
*       Scale and round the scores in a matrix so that the score of a word
*       is in the range [0..w*range].
*
**************************************************************************/
void scale_score_matrix(MATRIX_T *matrix, int rows, int cols,
    PRIOR_DIST_T *prior_dist, double alpha, int range,
    double *offset_p, double *scale_p) {
  int i, j;
  double small, large, scale, offset;

  // Get the largest and smallest scores in the PSSM.
  small = BIG;
  large = -BIG;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      double x;
      x = get_matrix_cell(i, j, matrix);
      small = MIN(small, x);
      large = MAX(large, x);
    }
  }

  // Get the smallest and largest prior log-odds from the prior distribution
  // and use them to adjust small and large.
  if (prior_dist != NULL) {
    double min_lo_prior, max_lo_prior;
    min_lo_prior = get_min_lo_prior(prior_dist, alpha);
    max_lo_prior = get_max_lo_prior(prior_dist, alpha);
    small = MIN(small, min_lo_prior);
    large = MAX(large, max_lo_prior);
  }
  
  // Find offset and scale factors so that PSSM scores for words is in the 
  // range: [0..w*range]
  // To make LO=0 map back to 0, need offset*scale to be an integer.
  // So we make offset and scale integers. (TLB 31 May 2013)
  if (large == small) { small = large - 1; }    // In case all motif entries are the same.
  offset = small = floor(small);                // Make offset an integer.
  scale = floor(range/(large-small));           // Ensure scaled scores are <= range.

  // Scale and round the PSSM entries.
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      double raw_score, scaled_score;
      raw_score = get_matrix_cell(i, j, matrix);
      scaled_score = raw_to_scaled(raw_score, 1, scale, offset);
      set_matrix_cell(i, j, scaled_score, matrix);
    }
  }

  // Record scale and offset of scores
  if (scale_p) *scale_p = scale;
  if (offset_p) *offset_p = offset;
}

/**************************************************************************
*       scale_pssm
*
*       Scale and round the scores in a PSSM so that the score of a word
*       is in the range [0..w*range].
*
*       Returns the scaled PSSM.
*
**************************************************************************/
void scale_pssm(
  PSSM_T *pssm,                   // The PSSM. (IN/OUT)
  PRIOR_DIST_T *prior_dist, // Distribution of priors (IN)
  double alpha,             // Fraction of all TFBS that are the TFBS of interest
  int range                               // The desired range. (IN) 
)
{
  scale_score_matrix(pssm->matrix, pssm->w, pssm->alphsize,
      prior_dist, alpha, range, &(pssm->offset), &(pssm->scale));
  pssm->range = range;
  pssm->matrix_is_scaled = true;

} // scale_pssm

/**************************************************************************
*       get_scaled_lo_prior_dist
*
*       Takes a scaled distribution of priors and creates a scaled distribution of
*       log odds priors. The parameters for the scaling of the input priors are
*       in the PRIOR_DIST_T data structure. The output distribution of log odss
*       priors are scaled to be in the same range as the PSSM   log odds using
*       the input parameters pssm_range, pssm_scale, and pssm_offset.
*
*       Special handling is required for a uniform distribution of priors.
*       In that case the max_prior == min_prior, and the distribution only 
*       contains one bin.
*
* Returns a new array containing the scaled log odds priors
**************************************************************************/
ARRAY_T  *get_scaled_lo_prior_dist(
  PRIOR_DIST_T *prior_dist,
  PSSM_T *pssm,
  double alpha,                    
  int range
) {

  assert(prior_dist != NULL);
  // Alocate enought space for elements in [0 ... range]
  ARRAY_T *scaled_lo_prior_dist = allocate_array(range + 1);

  if (prior_dist != NULL) {

    ARRAY_T *dist_array = get_prior_dist_array(prior_dist);
    int len_prior_dist = get_array_length(dist_array);
    double max_prior = get_prior_dist_maximum(prior_dist);
    double min_prior = get_prior_dist_minimum(prior_dist);
    double prior_dist_scale = get_prior_dist_scale(prior_dist);
    double prior_dist_offset = get_prior_dist_offset(prior_dist);
    init_array(0.0L, scaled_lo_prior_dist);
    if (max_prior == min_prior) {
        // Special case for uniform priors
        double value = 1.0;
        double lo_prior = my_log2(alpha * max_prior / (1.0L - (alpha * max_prior)));
        // Convert lo_prior to PSSM scale
        int scaled_index = raw_to_scaled(lo_prior, 1, pssm->scale, pssm->offset);
        set_array_item(scaled_index, value, scaled_lo_prior_dist);
    }
    else {
      int prior_index = 0;
      for (prior_index = 0; prior_index < len_prior_dist; ++prior_index) {
        double value = get_array_item(prior_index, dist_array);
        // Convert index giving scaled prior to raw prior.
        double scaled_prior = ((double) prior_index) + 0.5L;
        double prior \
          = scaled_to_raw(scaled_prior, 1, prior_dist_scale, prior_dist_offset);
        double lo_prior = my_log2(alpha * prior / (1.0L - (alpha * prior)));
        // Scale raw lo_prior using parameters from PSSM.
        int scaled_index = raw_to_scaled(lo_prior, 1, pssm->scale, pssm->offset);
        if (scaled_index < range) {
          double old_value = get_array_item(scaled_index, scaled_lo_prior_dist);
          set_array_item(scaled_index, value + old_value, scaled_lo_prior_dist);
        }
      }
    }
  }

  return scaled_lo_prior_dist;

}

/**************************************************************************
*
*       get_pdf_table
*
*       Compute the pdf of a pssm.
*
*       Returns an array of pdf values:
*               pdf[x] = Pr(score == x)
*       for 0<=x<=range*w.
*
*       Assumes:
*               1) motif scores are non-negative, integral
*               2) background model is position-dependent, 0-order Markov
*		3) if a distribution of log-odds priors is provided the range
*		   is extended to range * (w + 1)
*
**************************************************************************/
static ARRAY_T* get_pdf_table(
  PSSM_T* pssm,                                   // The PSSM.
  MATRIX_T* background_matrix,  // Background model PSSM matrix.
  ARRAY_T* scaled_lo_prior_dist // Scaled distribution of log odds priors.
)
{
  int i, j, k;
  MATRIX_T* matrix = pssm->matrix;// The PSSM matrix.
  int w = pssm->w;// PSSM width.
  int alen = pssm->alphsize;
  if (alen == alph_size_full(pssm->alph)) { 
    // CONSIDER We need to review how ambiguity characters are used 
    // in probability calculation throughout the code.
    // We don't want to include the background probabilities for ambiguity 
    // characters in this calculation, only the primary charcters of the alphabet.
    // However, Motiph the motiph 'alphabet' actually includes all possible
    // columns of the mutiple alignment, so we skip this clause if 
    // pssm->alphsize > alph_size_full(alph)
    alen = alph_size_core(pssm->alph);
  }
  int range = pssm->range;// Maximum score in PSSM.
  int size = w * range + 1;
  if (scaled_lo_prior_dist != NULL) {
    // Having priors expands the range of possible scores
    size = ((w + 1) * range) + 1;
  }
  ARRAY_T* pdf_old = allocate_array(size);
  ARRAY_T* pdf_new = allocate_array(size);

  init_array(0, pdf_new);
  set_array_item(0, 1, pdf_new);                // Prob(0)

  if (scaled_lo_prior_dist != NULL) {
    // Use distribution of log odds priors to 
    // initialize starting probabilities.
    for (k = 0; k <= range; k++) {
      double prob = get_array_item(k, scaled_lo_prior_dist);
      set_array_item(k, prob, pdf_new);
    }
  }

  // Compute the pdf recursively.
  for (i=0; i<w; i++) {

    int max = i * range;
    if (scaled_lo_prior_dist != NULL) {
      // Having priors expands the range of possible scores
      max += range;
    }

    // get position dependent background model
    ARRAY_T* background = get_matrix_row(i, background_matrix);
    SWAP(ARRAY_T*, pdf_new, pdf_old)
    for (k=0; k<=max+range; k++) {
      set_array_item(k, 0, pdf_new);
    }
    for (j=0; j<alen; j++) {
      int s = (int) get_matrix_cell(i, j, matrix);
      for(k=0; k<=max; k++) {
        double old = get_array_item(k, pdf_old);
        if (old != 0) {
          double new = get_array_item(k+s, pdf_new) +
            (old * get_array_item(j, background));
          set_array_item(k+s, new, pdf_new);
        } // old
      } // k
    } // j
  } // i

  // Free space.
  free_array(pdf_old);

  // Return the pdf.
  double total = array_total(pdf_new);

  return(pdf_new);
} // get_pdf_table

/**************************************************************************
*
*       get_pv_lookup
*
*       Create a lookup table for motif scores.
*
*       Sets an array of p-values:
*               pssm->pv[x] = Pr(score >= x)
*       for 0<=x<=range*w.
*
*       Assumes:
*               1) motif scores are non-negative, integral
*               2) background model is 0-order Markov
*
**************************************************************************/
void get_pv_lookup(
  PSSM_T* pssm,                             // The PSSM.
  ARRAY_T* background,          // The background model.
  ARRAY_T* scaled_lo_prior_dist // Scaled distribution of log odds priors.
)
{
  int i;
  int w = pssm->w;              // PSSM width.
  int alen = pssm->alphsize;    // PSSM alphabet size.

  assert(pssm->range > 0);

  // Create a temporary position-dependent background where all positions
  // are the same.
  MATRIX_T* background_matrix = allocate_matrix(0, get_array_length(background));
  for (i=0; i<w; i++) grow_matrix(background, background_matrix);

  get_pv_lookup_pos_dep(pssm, background_matrix, scaled_lo_prior_dist);

  free_matrix(background_matrix);
} // get_pv_lookup

/**************************************************************************
*
*       get_pv_lookup_pos_dep
*
*       Create a lookup table for motif scores where the background
*       model is position-dependent and has the same width as the
*       motif.
*
*       Sets an array of p-values:
*               pssm->pv[x] = Pr(score >= x)
*       for 0<=x<=range*w.
*
*       Assumes:
*               1) motif scores are non-negative, integral
*               2) background model is 0-order Markov
*
**************************************************************************/
void get_pv_lookup_pos_dep(
  PSSM_T* pssm,                             // The PSSM.
  MATRIX_T* background_matrix,  // The background model PSSM matrix.
  ARRAY_T* scaled_lo_prior_dist // Scaled distribution of log odds priors.
)
{
  int i;
  int w = pssm->w;              // PSSM width.
  int range = pssm->range;      // Maximum score in PSSM.
  int size = (scaled_lo_prior_dist == NULL 
    ?  w * range + 1 
    :  (w + 1) * range + 1);

  // Free any old pv array.
  if (pssm->pv != NULL) { free_array(pssm->pv); }

  // Get the pdf.
  ARRAY_T* pv = pssm->pv 
    = get_pdf_table(pssm, background_matrix, scaled_lo_prior_dist);

  // Get 1-cdf from the pdf.
  for (i=size-2; i>=0; i--) {
    double p_iplus1 = get_array_item(i+1, pv);
    double p_i = get_array_item(i, pv);
    double p = p_i + p_iplus1;
    set_array_item(i, MIN(1.0, p), pv);
    if (pssm->max_score==0 && p_iplus1 > 0) pssm->max_score = i+1; 
    if (p_i > 0) pssm->min_score = i;
  }
} // get_pv_lookup_pos_dep

/************************************************************************
 * See .h file for set_up_pssms_and_pvalues.
 *
 *************************************************************************/
static void set_up_pssm_and_pvalue(
  double p_threshold,           // Scale/offset PSSM and create table if > 0
  bool use_both_strands,   // Compute PSSM for negative strand, too?
  bool allow_weak_motifs,  // Allow motifs with min p-value < p_threshold?
  int      i_state,             // Index of starting state of motif.
  MHMM_T *the_hmm,              // The HMM.
  PRIOR_DIST_T *prior_dist, // Distribution of priors
  double alpha // PSP alpha parameter
)
{
  int i, len;
  int i_motif = the_hmm->states[i_state].i_motif;       // Name of motif.
  MHMM_STATE_T *start_state = the_hmm->states + i_state;// Starting state of motif.
  ALPH_T *alph = the_hmm->alph;
  //
  // Create the PSSM for the motif, storing it in the first state.
  //
  MATRIX_T* pssm_matrix = NULL;
  for (i=0; i<the_hmm->num_states; i++) {       // state
    MHMM_STATE_T* state = the_hmm->states + i;  // Current state.
    if (state->type == START_MOTIF_STATE && state->i_motif == i_motif) {
      pssm_matrix = allocate_matrix(state->w_motif, alph_size_full(alph));
      set_matrix_row(state->i_position, state->emit_odds, pssm_matrix);
    } else if (state->type == MID_MOTIF_STATE && state->i_motif == i_motif) {
      set_matrix_row(state->i_position, state->emit_odds, pssm_matrix);
    } else if (state->type == END_MOTIF_STATE && state->i_motif == i_motif) {
      set_matrix_row(state->i_position, state->emit_odds, pssm_matrix);
    }
  } // state
  // Created unscaled pssm with no p-valuve table.
  start_state->pssm = build_matrix_pssm(alph, pssm_matrix, NULL, prior_dist, alpha, 0);

  //
  // Scale PSSM and get p-value lookup table.
  //
  if (p_threshold > 0) {

    // Convert background to non-log form.
    ARRAY_T *background = allocate_array(get_array_length(the_hmm->background));
    convert_to_from_log_array(false, the_hmm->background, background);

    // Scale pssm and get p-value table.
    scale_pssm(
      start_state->pssm, 
      prior_dist,
      alpha,
      PSSM_RANGE
    );
    // Scale prior dist
    ARRAY_T *scaled_prior_dist = NULL;
    if (prior_dist != NULL) {
      scaled_prior_dist = get_scaled_lo_prior_dist(
          prior_dist, 
          start_state->pssm,
          alpha, 
          PSSM_RANGE
        );
    }
    // Get p-value table.
    get_pv_lookup(
      start_state->pssm,
      background,
      scaled_prior_dist
    );

    // Check that p_threshold is not too small.
    start_state->min_pvalue = get_min_pvalue(
      start_state->pssm,
      prior_dist,
      alpha
    );
    if (start_state->min_pvalue > p_threshold) {
      if (allow_weak_motifs) {
        fprintf(
          stderr,
          "Warning: Weak motif %s: minimum p-value (%.2g) > p-value threshold (%.2g)\n",
          the_hmm->states[i_state].motif_id,
          start_state->min_pvalue,
          p_threshold
        );
      } else {
        die(
          "Weak motif %s: minimum p-value (%.2g) > p-value threshold (%.2g)",
          the_hmm->states[i_state].motif_id,
          start_state->min_pvalue,
          p_threshold
        );
      }
    }
    // Find minimum motif score whose p-value is < p_threshold.
    len = get_array_length(start_state->pssm->pv);
    for (i=0; i < len && get_array_item((int) i, start_state->pssm->pv) > p_threshold; i++);
    if (i == len) {
      fprintf(
        stderr,
        "Warning: Motif %s has no scores with p-value < p-value threshold (%.2g)\n",
        the_hmm->states[i_state].motif_id,
        p_threshold
      );
    }
    start_state->min_sig_score = i-1;
    free_array(background);
  }

  //
  // Create PSSM for negative strand.
  //
  if (use_both_strands) {
    reverse_complement_pssm(alph, pssm_matrix);
    start_state->npssm = build_matrix_pssm(alph, pssm_matrix, NULL, prior_dist, alpha, 0);
  }
  free_matrix(pssm_matrix);

  //
  // Create a hashed PSSM for faster scanning.
  //
  if (hash_n > 1) {
    hash_pssm(start_state->pssm, hash_n);
  }

} // set_up_pssm_and_pvalue

/************************************************************************
 * See .h file.
 *************************************************************************/
void set_up_pssms_and_pvalues (
  bool motif_scoring,      // Doing motif scoring?
  double p_threshold,           // Scale/offset PSSM and create table if > 0
  bool use_both_strands,   // Compute negative PSSM?
  bool allow_weak_motifs,  // Allow motifs with min p-value < p_threshold?
  MHMM_T *the_hmm,
  PRIOR_DIST_T *prior_dist, // Distribution of priors
  double alpha // PSP alpha parameter
)
{
  int i_state;
  int num_states = the_hmm->num_states;
  MHMM_STATE_T* this_state;

  // Look through the model for start motif states.
  if (motif_scoring) {
    for (i_state = 0; i_state < num_states; i_state++) {
      this_state = &(the_hmm->states[i_state]);

      if (this_state->type == START_MOTIF_STATE) {
        set_up_pssm_and_pvalue(
          p_threshold,
          use_both_strands,
          allow_weak_motifs,
          i_state,
          the_hmm,
          prior_dist,
          alpha
        );
      }
    }
  }

  //
  // Set up the list of "hot" states for faster DP.
  // This is all states if not in motif-scoring mode.
  // In motif-scoring mode, it does not include the following states:
  //    START_STATE
  //    START_MOTIF_STATE
  //    MID_MOTIF_STATE
  //
  the_hmm->hot_states = NULL;
  mm_resize(the_hmm->hot_states, num_states, int);
  for (i_state=the_hmm->num_hot_states=0; i_state<num_states; i_state++) { // state
    MHMM_STATE_T *this_state = &the_hmm->states[i_state];
    // Skip cold states.  Always skip the START_STATE.
    if (this_state->type == START_STATE ||
         (motif_scoring &&
           (
             this_state->type == START_MOTIF_STATE ||
             this_state->type == MID_MOTIF_STATE
           )
          )
        ) {
      continue;
    }
    the_hmm->hot_states[the_hmm->num_hot_states++] = i_state;
  } // state
  mm_resize(the_hmm->hot_states, the_hmm->num_hot_states, int);

} // set_up_pssms_and_pvalues

/**************************************************************************
 *      compute_motif_scores
 *
 *      Score each position in the given sequence with a given motif.
 *      Return an array of scores:
 *
 *              s[i] = score of motif starting at position i in sequence
 *
 *      The score is either the raw score determined by the motif PSSM,
 *      or minus the logarithm (base 2) of the p-value of the raw score.
 *
 *      Negative scores indicate that the best score is on the
 *      negative strand.  The negative strand is only scored if the
 *      npssm is non-null in the start state of the motif.
 *
 *      Positions i where motif won't fit (at right edge of sequence)
 *      have s[i] = 0.
 *
 **************************************************************************/
static void compute_motif_scores (
  bool use_pvalues,   // Returns scores as p-values, not log-odds.
  double    p_threshold,   // Divide p-values by this.
  int       i_state,       // Index of the start state of the motif.
  MHMM_T*   the_hmm,       // The HMM, in log form.
  int*      int_sequence,  // The sequence, in integer format.
  size_t    seq_length,    // The length of the sequence.
  double *priors, // Array of priors
  size_t num_priors, // Number of priors
  double alpha, // Alpha parameter for calculating prior log odds
  ARRAY_T*  score_array    // The scores (pre-allocated array).
)
{
  int j, k;
  size_t i;
  size_t length = seq_length;                              // Length of sequence.
  PSSM_T *pssm = the_hmm->states[i_state].pssm;
  MATRIX_T* pssm_matrix = pssm->matrix;        // PSSM for motif.
  int w = the_hmm->states[i_state].w_motif;             // The motif width.
  ARRAY_T* pv = pssm->pv;      // Lookup table for p-values.
  double min_sig_score = the_hmm->states[i_state].min_sig_score; // Smallest significant score.
  double off = use_pvalues? my_log2(p_threshold) : 0;   // Convert p_threshold to bits.

  //
  // Score motif at each position it fits.
  // Motif may not hang off edge of sequences.
  // First position in sequence is not part of sequence.
  // Last position in sequence is not part of sequence.
  //
  set_array_item(0, 0, score_array);            // Position 0 not part of sequence.
  for (i=1; length>w && i<length-w; i++) {  	// motif start (size_t is unsignd!)
    int* seq = int_sequence+i;                  // Letter at current position.
    double score = 0;                           // Score for match.
    // i : position in original sequence.
    // j : position in hashed sequence.
    // k : position in hashed pssm.
    // seq[j] : letter in hashed sequence.
    for (j=k=0; j<w; j+=hash_n, k++, seq+=hash_n) {
      score += get_matrix_cell(k, *seq, pssm_matrix);
    }
    if (priors != NULL) {
      // The sequence has 1 char of padding at begining and end
      // the priors don't
      double prior = alpha * priors[i - 1];
      if (!isnan(prior)) {
        double log_prior_odds_score = my_log2(prior / (1.0L - prior));
        double scaled_log_prior_odds_score = raw_to_scaled(log_prior_odds_score, 1, pssm->scale, pssm->offset);
        score += scaled_log_prior_odds_score;
      }
    }
    double motif_score = 0.0L;
    if (use_pvalues) {
      double pvalue = get_array_item((int) score, pv);
      double log2_pvalue = my_log2(pvalue);
      motif_score = (score>min_sig_score) ? off-log2_pvalue : off;
    }
    set_array_item(i, motif_score, score_array);
  } // motif start

  // Set score where motif would run off edge to 0.
  for (i = (length > w ? length - w : 0); i < length; i++) set_array_item(i, 0, score_array);

} // compute_motif_scores

/**************************************************************************
 * Compute a matrix, for a given sequence, of motif scores w.r.t. each
 * sequence position.  The matrix may be pre-allocated, and it will
 * grow if necessary.
 **************************************************************************/
void compute_motif_score_matrix(
  bool  use_pvalues,       // Returns scores as p-values, not log-odds.
  double     p_threshold,       // Divide p-values by this.
  int*       int_sequence,      // Sequence as integer indices into alphabet.
  size_t     seq_length,        // Length of sequence.
  double *priors, // Priors for sequence
  size_t num_priors, // Number of priors
  double alpha, // Alpha parameter for calculating prior log odds
  MHMM_T*    the_hmm,           // The hmm.
  MATRIX_T** motif_score_matrix
)
{
  MHMM_STATE_T* this_state;
  int num_states;
  int i_state;
  int* hashed_sequence = NULL;  // Hashed version of sequence for efficiency.

  // Allocate or re-allocate space, as necessary.
  if ((*motif_score_matrix != NULL) &&
      (get_num_cols(*motif_score_matrix) < seq_length)) {
    free_matrix(*motif_score_matrix);
    *motif_score_matrix = NULL;
  }
  if (*motif_score_matrix == NULL) {
    *motif_score_matrix = allocate_matrix(the_hmm->num_motifs, seq_length);
  }

  // Hash the sequence to compressed format.
  hashed_sequence = hash_n>1 ?
    hash_sequence(the_hmm->alph, int_sequence, seq_length, hash_n) :
    int_sequence;

  // Look through the model for start motif states.
  num_states = the_hmm->num_states;
  for (i_state = 0; i_state < num_states; i_state++) {
    this_state = &(the_hmm->states[i_state]);

    if (this_state->type == START_MOTIF_STATE) {

      // Compute the scores for this motif at all sequence positions.
      compute_motif_scores(use_pvalues,
                           p_threshold,
                           i_state,
                           the_hmm,
                           hashed_sequence,
                           seq_length,
         priors,
         num_priors,
         alpha,
                           get_matrix_row(this_state->i_motif,
                                          *motif_score_matrix));
    }
  }

  // Free the hashed sequence.
  if (hash_n>1) myfree(hashed_sequence);
} // compute_motif_score_matrix

/*************************************************************************
 *  Build a Position Specific Scoring Matrix (PSSM) for a motif.
 *  The matrix will contain one row for each position in the motif that 
 *  has not been trimmed and one column for each possible nucleotide.
 *  The matrix will be log likelihood ratio (or likelihood ratio)
 *  and be scaled to the scores will be in the range [0..w*range].
 *  The PSSM object will contain the score and offset as well as
 *  the p-value lookup table. It will contain a pointer to the original
 *  motif data.
 *
 *  Caller is responsible for freeing the PSSM.
 *************************************************************************/
PSSM_T* build_motif_pssm(
  MOTIF_T* motif, // motif frequencies p_ia (IN)
  ARRAY_T* bg_freqs, // background frequencies b_a for pssm (IN)
  ARRAY_T* pv_bg_freqs, // background frequencies b_a for p-values (IN)
  PRIOR_DIST_T* prior_dist, // Distribution of priors. May be NULL (IN)
  double alpha, // Scale factor for non-specific priors. Unused if prior_dist is NULL.
  int range, // range of scaled scores is [0..w*range]
  // create pv tables for this number of GC bins
  // instead of using the pv_bg_freqs
  int num_gc_bins, 
  bool no_log // make likelihood ratio 
)
{
  assert(motif != NULL);
  assert(bg_freqs != NULL);
  assert(pv_bg_freqs != NULL);
  assert (range > 0);

  ALPH_T *alph = get_motif_alph(motif);
  const int a_size = alph_size_core(alph);
  const int trim_left = get_motif_trim_left(motif);
  const int trim_right = get_motif_trim_right(motif);
  const int w = get_motif_length(motif) - trim_left - trim_right;

  PSSM_T* pssm = allocate_pssm(alph, w, a_size, num_gc_bins);
  pssm->matrix_is_log = no_log ? false : true;
  pssm->matrix_is_scaled = true;

  // Process in column major order to avoid re-calculating bg_likelihood
  MATRIX_T* saved_pssm_matrix = NULL;
  if (no_log) {
    saved_pssm_matrix 
      = allocate_matrix(get_num_rows(pssm->matrix), get_num_cols(pssm->matrix));
  }
  double total_bg_likelihood = 0.0;
  double total_fg_likelihood = 0.0;
  int alph_index;
  for (alph_index = 0; alph_index < a_size; alph_index++) {
    double bg_likelihood = get_array_item(alph_index, bg_freqs);
    total_bg_likelihood += bg_likelihood;
    int motif_position_index, pssm_position_index;
    for (motif_position_index = trim_left, pssm_position_index = 0;
      pssm_position_index < w;
      motif_position_index++, pssm_position_index++) {
      double fg_likelihood =
        get_matrix_cell(motif_position_index, alph_index, get_motif_freqs(motif));
      total_fg_likelihood += fg_likelihood;
      double odds = fg_likelihood / bg_likelihood;
      // Caution: we're avoiding the overhead of a function call here
      // by operating directly on the matrix data structure.
      if (no_log) saved_pssm_matrix->rows[motif_position_index]->items[alph_index] = odds;
      pssm->matrix->rows[pssm_position_index]->items[alph_index] = my_log2(odds);
    }
  }

  // Check the input.
  const double epsilon = 0.001;
  assert((total_bg_likelihood - 1.0) < epsilon);
  assert((total_fg_likelihood - (double) w) < epsilon);
  
  // Scale the pssm and set the scale and offset in the PSSM object.
  scale_pssm(pssm, prior_dist, alpha, range);
  ARRAY_T *scaled_lo_prior_dist = NULL;
  if (prior_dist != NULL) {
    scaled_lo_prior_dist 
      = get_scaled_lo_prior_dist(prior_dist, pssm, alpha, range);
  }

  // Get the p-value lookup table for the scaled, log version of the PSSM. 
  if (num_gc_bins <= 1) {                       // not using GC binning
    get_pv_lookup(pssm, pv_bg_freqs, scaled_lo_prior_dist);
  } else { // using complementary pair binning (aka "GC binning")
    // create pv tables for different complementary pair backgrounds
    assert(alph_size_pairs(alph) == 2);
    int i, i_gc_bin, x1, x2, y1, y2;
    // determine complementary pairs
    x1 = 0;
    x2 = alph_complement(alph, x1);
    y1 = (x2 == 1 ? 2 : 1);
    y2 = alph_complement(alph, y1);
    // Allocate the background array big enough to include ambigs.
    ARRAY_T *bg = allocate_array(alph_size_full(alph));
    for (i_gc_bin=0; i_gc_bin < num_gc_bins; i_gc_bin++) { // GC bin
      // Create background model with given GC content
      double pair_prob = get_pair_bin_value(num_gc_bins, i_gc_bin);
      set_array_item(x1, (1-pair_prob)/2.0, bg);
      set_array_item(x2, (1-pair_prob)/2.0, bg);
      set_array_item(y1, pair_prob/2.0, bg);
      set_array_item(y2, pair_prob/2.0, bg);
      // Extend the distribution to account for ambiguous characters.
      calc_ambigs(alph, false, bg);
      int min_score = i_gc_bin==0 ? 0 : pssm->min_score;
      get_pv_lookup(pssm, bg, scaled_lo_prior_dist);
      pssm->gc_pv[i_gc_bin] = pssm->pv;
      pssm->min_score = min_score;              // don't let min_score change
      pssm->pv = NULL;
    }
    free_array(bg);
  }

  free_array(scaled_lo_prior_dist);

  // Use unscaled matrix for scoring?
  if (no_log) { 
    copy_matrix(saved_pssm_matrix, pssm->matrix);
    free_matrix(saved_pssm_matrix);
    pssm->matrix_is_scaled = false;
    pssm->nolog_max_score = 1;
    // Get the largest score in each row and multiply them.
    int i, j;
    for (i=0; i<pssm->w; i++) {
      double large = -BIG;
      for (j=0; j<pssm->alphsize; j++) {
	double x = get_matrix_cell(i, j, pssm->matrix);
	large = MAX(large, x);
      }
      pssm->nolog_max_score *= large;
    }
  } // no_log

  return pssm;
}

/***********************************************************************
 * create pssm from a scoring matrix
 *
 *      Builds a PSSM from a matrix.
 *      Scales it if range > 0;
 *      Computes the distribution if bg_freqs is not NULL.
 *
 ***********************************************************************/
PSSM_T* build_matrix_pssm(
  ALPH_T *alphabet,          // alphabet (IN)
  MATRIX_T* matrix,     // pssm matrix (IN)
  ARRAY_T* bg_freqs,    // background frequencies b_a (IN)
  PRIOR_DIST_T *prior_dist, // Prior distribution
  double alpha, // PSP alpha parameter
  int range             // range of scaled scores is [0..w*range] (IN)
)
{
  assert(matrix != NULL);

  int w = get_num_rows(matrix);
  int asize = get_num_cols(matrix);// could be hashed so no guarentee that it's the same size as alph_size_core(alph)

  // create the pssm
  PSSM_T* pssm = allocate_pssm(alphabet, w, asize, 0);
  pssm->matrix_is_log = true;
  //FIXME what if the matrix is too large because of ambiguous characters?!
  copy_matrix(matrix, pssm->matrix); 

  // Scale the pssm and set the scale and offset in the PSSM object.
  if (range > 0) {
    scale_pssm(
      pssm, 
      prior_dist,
      alpha,
      range
    );  
  }

  // Get the pv lookup table if background given.
  if (bg_freqs != NULL) {
    get_pv_lookup(pssm, bg_freqs, NULL /* no prior distribution */ );
  }

  return(pssm);
} // build_matrix_pssm

/**********************************************************************
 average_two_pvs() 

 Compute the pv distribution for the average of two distributions.
 Distributions are in two rows of the pssm_pair's gc_n_pv_lookup table.
 **********************************************************************/
ARRAY_T* average_two_pvs(
  PSSM_PAIR_T* pssm_pair,               // pssms for pos and neg motifs
  const int r1,                         // row with first pv distr.
  const int r2,                         // row with second pv distr.
  const int gcbin                       // which gc_n_pv_lookup table to use
)
{
  MATRIX_T* n_pv_lookup = pssm_pair->gc_n_pv_lookup[gcbin];// pv(log2(n), score) table
  ARRAY_T* scaled_to_ama = pssm_pair->scaled_to_ama;    // precomputed map
  PSSM_T* pssm = pssm_pair->pos_pssm;
  const int w = pssm->w;                        // width of pssms
  const double scale = pssm->scale;             // scale of pssms
  const double offset = pssm->offset;           // offset of pssms
  int min_score = pssm->min_score;

  // create space for new distribution
  int ncols = get_num_cols(n_pv_lookup);// number of possible scores
  ARRAY_T* new_pv = allocate_array(ncols);

  int s1, s2;           // scaled log scores
  for (s1 = min_score; s1 < ncols - 1; s1++) {
    // get probability (pdf) of scaled score in row r1
    double p1 = get_matrix_cell(r1, s1, n_pv_lookup) - get_matrix_cell(r1, s1+1, n_pv_lookup); 
    if (p1 == 0) continue;
    double ama1 = get_array_item(s1, scaled_to_ama);
    if (ama1 == 0) ama1 = SMALL_POS;
    //for (s2=min_score; s2<ncols-1; s2++) {
    // Only do above the diagonal to save time if convolving a row with itself.
    for (s2 = (r1==r2 ? s1 : min_score); s2 < ncols-1; s2++) {
      // get probability (pdf) of scaled score in row r2
      double p2 = get_matrix_cell(r2, s2, n_pv_lookup) - get_matrix_cell(r2, s2+1, n_pv_lookup); 
      if (p2 == 0) continue;
      // convert scaled log scores to ama scores
      double ama2 = get_array_item(s2, scaled_to_ama);
      if (ama2 == 0) ama2 = SMALL_POS;
      // average the two ama scores
      double ama = (ama1 + ama2)/2;
      // convert new ama score to scaled log score
      double log_ama = my_log2(ama);
      int s = raw_to_scaled(log_ama, w, scale, offset);
      assert(s < ncols);
      // Add the joint probability to the new PDF.
      // Double it for all elements not on the diagonal if convolving
      // a row with itself.
      double joint = (r1 == r2 && s1 != s2) ? 2*p1*p2 : p1*p2;
      double sum = get_array_item(s, new_pv) + joint;
      set_array_item(s, sum, new_pv);
    }
  }

  // Convert the PDF to 1-CDF (p-values)
  //printf("pdf: irow %d\n", irow);
  //print_array(new_pv, 8, 6, true, stdout);
  int s;
  for (s=ncols-2; s>=min_score; s--) {
    double p_splus1 = get_array_item(s+1, new_pv);
    double p_s = get_array_item(s, new_pv);
    double p = p_s + p_splus1;
    set_array_item(s, MIN(1.0, p), new_pv);
  }

  return(new_pv);
}

/**********************************************************************
  get_ama_pv()

  Get the pv of an AMA score.
  Builds the pv tables for AMA scores "on the fly" from the
  pv tables for scaled, log likelihood ratio scores stored in the pssms.
 **********************************************************************/
double get_ama_pv(
  double ama_score,                     // average likelihood ratio score
  int seqlen,                           // length of sequence scanned
  double seq_gc,                        // total GC content of sequence
  PSSM_PAIR_T* pssm_pair                // pssms for pos and neg motifs
)
{
  assert(ama_score >= 0);
  assert(seqlen > 0);
  assert(seq_gc == -1 || (seq_gc >= 0 && seq_gc <= 1));
  assert(pssm_pair != NULL);

  PSSM_T* pos_pssm = pssm_pair->pos_pssm;
  PSSM_T* neg_pssm = pssm_pair->neg_pssm;
  assert(pos_pssm != NULL);
  int w = pos_pssm->w;                          // width of pssm
  double scale = pos_pssm->scale;
  double offset = pos_pssm->offset;
  int num_gc_bins = pssm_pair->num_gc_bins;
  // number of possible scores
  int ncols = num_gc_bins <= 1 ? get_array_length(pos_pssm->pv) : get_array_length(pos_pssm->gc_pv[0]);
  // Get number of pssm scores in average. 
  int n = seqlen - w + 1;               // number of scores in average
  int i_gc_bin;

  // Make sure pssms are compatible.
  if (neg_pssm != NULL) assert(scale == neg_pssm->scale);
  if (neg_pssm != NULL) assert(offset == neg_pssm->offset);

  // Convert score to scaled average likelihood ratio score.
  double log_odds = my_log2(ama_score);
  int scaled_score = get_scaled_pssm_score(log_odds, pos_pssm);
  assert(scaled_score < ncols);

  // Precalculate scaled-to-ama-score conversion table for speed.
  if (! pssm_pair->scaled_to_ama) {
    pssm_pair->scaled_to_ama = allocate_array(ncols);
    int s;
    for (s = 0; s < ncols; s++) {
      double raw_score = scaled_to_raw(s, w, scale, offset);
      set_array_item(s, EXP2(raw_score), pssm_pair->scaled_to_ama);
    }
  }

  //
  // Set up x_1 and x_2, y_1 and y_2 for bilinear interpolation of p-values
  //
  int n_1, n_2; 
  int gc_1, gc_2;
  // Sequence length bins:
  // FIXME: test if linear n-axis better
  double seq_n_bin = my_log2(n);        // make n-axis logarithmic
  n_1 = floor(seq_n_bin);
  n_2 = n_1 == seq_n_bin ? n_1 : n_1 + 1;       // avoid making extra table if on boundary
  // GC bins:
  double seq_gc_bin;
  if (num_gc_bins <= 1) {               // not using GC binning
    seq_gc_bin = 0;
    gc_1 = gc_2 = 0;
  } else {                              // using GC binning
    seq_gc_bin = get_pair_bin(num_gc_bins, seq_gc);
    gc_1 = get_pair_bin_lower(num_gc_bins, seq_gc);
    gc_2 = gc_1 == seq_gc_bin ? gc_1 : gc_1 + 1;        // avoid bin for gc > 1.0
  }

  //
  // First time a GC content is seen, create tables for n=1.
  // Create tables for GC contents corresponding to bins gc_1 and/or gc_2.
  //
  // The pv lookup table for the average of n scores will be
  // in row log_2(N), for N=1, 2, 4, ...
  if (!pssm_pair->gc_n_pv_lookup[gc_1] || !pssm_pair->gc_n_pv_lookup[gc_2]) {
    for (i_gc_bin=gc_1; i_gc_bin<=gc_2; i_gc_bin++) {
      if (verbosity > NORMAL_VERBOSE) {
        fprintf(stderr, "Creating pv table for gc_bin %d n= %d (%d)\n", i_gc_bin, 1, 0);
      }
      MATRIX_T* n_pv_lookup =
        pssm_pair->gc_n_pv_lookup[i_gc_bin] = allocate_matrix(2, ncols);
      if (num_gc_bins == 1) {           // Use pssm->pv table (based on background freqs)
        set_matrix_row(0, pos_pssm->pv, n_pv_lookup);
        if (neg_pssm) set_matrix_row(1, neg_pssm->pv, n_pv_lookup);
      } else {                          // Use pssm->gc_pv (based on GC binning)
        set_matrix_row(0, pos_pssm->gc_pv[i_gc_bin], n_pv_lookup);
        if (neg_pssm) set_matrix_row(1, neg_pssm->gc_pv[i_gc_bin], n_pv_lookup);
      }
      if (neg_pssm) {                   // Average pos and neg distributions
        ARRAY_T* new_pv = average_two_pvs(pssm_pair, 0, 1, i_gc_bin);
        set_matrix_row(0, new_pv, n_pv_lookup);
        free_array(new_pv);             // free scratch array
      }
      remove_matrix_row(1, n_pv_lookup);// scratch row
    } // i_gc_bin
  } // First time one of gc_1 or gc_2 was seen.

  // Do matrices contain the lookup table for the current n for both gc_1 and gc_2?
  int nrows_1 = get_num_rows(pssm_pair->gc_n_pv_lookup[gc_1]);
  int nrows_2 = get_num_rows(pssm_pair->gc_n_pv_lookup[gc_2]);
  if (nrows_1 <= n_2 || nrows_2 <= n_2) {
    for (i_gc_bin=gc_1; i_gc_bin<=gc_2; i_gc_bin++) {
      int nrows = get_num_rows(pssm_pair->gc_n_pv_lookup[i_gc_bin]);
      if (nrows > n_2) continue;                        // Table rows exist already for this GC.
      MATRIX_T* n_pv_lookup = pssm_pair->gc_n_pv_lookup[i_gc_bin];
      // Create pv lookup tables for n up to 2**n_2.
      int irow;
      for (irow=nrows; irow<=n_2; irow++) {
        // Create the next pv_table for n = 2**irow by convolving the pdf for irow-1 with itself.
        if (verbosity > NORMAL_VERBOSE) {
          fprintf(stderr, "Creating pv table for gc_bin %d n= %d (%d)\n", i_gc_bin, (int) pow(2,irow), irow);
        }
        ARRAY_T* new_pv = average_two_pvs(pssm_pair, irow-1, irow-1, i_gc_bin);
        // Add new row to n_pv_lookup matrix.
        grow_matrix(new_pv, n_pv_lookup);
      } // irow
    } // i_gc_bin
  } // using GC bins

  //
  // Compute p-value using (bi)linear interpolation
  //
  double pv1 = get_matrix_cell(n_1, scaled_score, pssm_pair->gc_n_pv_lookup[gc_1]);
  double pv2 = get_matrix_cell(n_2, scaled_score, pssm_pair->gc_n_pv_lookup[gc_1]);
  double r1 = n_1==n_2 ? pv1 : ((n_2 - seq_n_bin) * pv1) + ((seq_n_bin - n_1) * pv2);
  double pv;
  if (num_gc_bins == 1) {               // linear interpolation on N only
    pv = r1;
  } else {                              // bilinear interpolation on N and GC
    pv1 = get_matrix_cell(n_1, scaled_score, pssm_pair->gc_n_pv_lookup[gc_2]);
    pv2 = get_matrix_cell(n_2, scaled_score, pssm_pair->gc_n_pv_lookup[gc_2]);
    //double r2 = (((n_1+1) - seq_n_bin) * pv1) + ((seq_n_bin - n_1) * pv2);
    double r2 = n_1==n_2 ? pv1 : ((n_2 - seq_n_bin) * pv1) + ((seq_n_bin - n_1) * pv2);
    pv = gc_1==gc_2 ? r1 : ((gc_2 - seq_gc_bin) * r1) + ((seq_gc_bin - gc_1) * r2);
  }
 
  return(pv); 
}

/**********************************************************************
  create_pssm_pair()

  Note: Keeps pointers to PSSMs in the object, so they get freed
  when it does.
 **********************************************************************/
PSSM_PAIR_T* create_pssm_pair(
  PSSM_T* pos_pssm,             // positive strand pssm
  PSSM_T* neg_pssm              // negative strand pssm
)
{
  assert(pos_pssm);

  // Set the number of GC bins to 1 if not using them since
  // we use the first entry in the gc_n_pv_lookup array in that case.
  int num_gc_bins = MAX(1, pos_pssm->num_gc_bins);

  PSSM_PAIR_T* pssm_pair = mm_malloc(sizeof(PSSM_PAIR_T));

  pssm_pair->pos_pssm = pos_pssm; 
  pssm_pair->neg_pssm = neg_pssm; 
  pssm_pair->num_gc_bins = num_gc_bins;
  pssm_pair->gc_n_pv_lookup = mm_calloc(num_gc_bins, sizeof(MATRIX_T *));
  pssm_pair->scaled_to_ama = NULL;

  return(pssm_pair);
}

/**********************************************************************
  free_pssm_pair()

  Note: Frees the PSSMs.
 **********************************************************************/
void free_pssm_pair(
  PSSM_PAIR_T* pssm_pair
)
{
  free_pssm(pssm_pair->pos_pssm);
  free_pssm(pssm_pair->neg_pssm);
  int i;
  for (i=0; i<pssm_pair->num_gc_bins; i++) free_matrix(pssm_pair->gc_n_pv_lookup[i]);
  myfree(pssm_pair->gc_n_pv_lookup);
  free_array(pssm_pair->scaled_to_ama);
  free(pssm_pair);
}

/**************************************************************************
 * Free a given PSSM.
 **************************************************************************/
void free_pssm(
  PSSM_T* pssm
)
{
  if (pssm == NULL) {
    return;
  }

  free_matrix(pssm->matrix);
  free_array(pssm->pv);
  int i;
  for (i=0; i<pssm->num_gc_bins; i++) free_array(pssm->gc_pv[i]);
  myfree(pssm->gc_pv);
  alph_release(pssm->alph);

  myfree(pssm);
} // free_pssm


/**************************************************************************
 * Sum up the largest values for each position in the pssm to generate 
 * the best possible score for a match. Output is unscaled. 
 **************************************************************************/
double pssm_best_match_score(
  PSSM_T* pssm
)
{
  double score = 0.0L, position_score, best_score;
  int motif_offset, motif_length, alphabet_index, alphabet_size;
  // get the motif length
  motif_length = get_pssm_w(pssm);
  //get the alphabet length
  alphabet_size = get_pssm_alphsize(pssm); 
  // get the largest score from every position of the motif
  for (motif_offset = 0; motif_offset < motif_length; ++motif_offset) {
    best_score = get_matrix_cell(motif_offset, 0, pssm->matrix);
    for (alphabet_index = 1; alphabet_index < alphabet_size; ++alphabet_index) {
      position_score = get_matrix_cell(motif_offset, alphabet_index, pssm->matrix);
      if (position_score > best_score) best_score = position_score;
    }
    score += best_score;
  }
  
  return get_unscaled_pssm_score(score, pssm);
}

// Print a PSSM to standard error.
void print_pssm(
  PSSM_T* pssm
)
{
  int motif_offset, motif_length, alphabet_index, alphabet_size;
  // get the motif length
  motif_length = get_pssm_w(pssm);
  //get the alphabet length
  alphabet_size = get_pssm_alphsize(pssm); 
  fprintf(stderr, "PSSM:\n");
  for (motif_offset = 0; motif_offset < motif_length; ++motif_offset) {
    for (alphabet_index = 0; alphabet_index < alphabet_size; ++alphabet_index) {
      double score = get_matrix_cell(motif_offset, alphabet_index, pssm->matrix);
      fprintf(stderr, "%.3f ", score);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
} // print_pssm

// comparator for binary search
int double_bsearch_compar(
  const void *key,
  const void *entry
) {
  double sought_pv = *((double*)key);
  double found_pv = *((double*)entry);
  if (sought_pv < found_pv) {
    return +1;
  } else if (sought_pv > found_pv) {
    return -1;
  } else {
    return 0;
  }
}

// Get the unscaled ("raw") score corresponding to a given p-value using binary search.
double pssm_pvalue_to_score(
  PSSM_T *pssm,			// a pssm
  double pvalue			// a p-value
) 
{
  int len = get_pssm_pv_length(pssm);
  double *pv = mm_malloc(len * sizeof(double));

  // copy the p-values to an array
  int i;
  for (i=0; i<len; i++) {
    pv[i] = get_pssm_pv(i, pssm);
  }
  
  // perform binary search
  int scaled_score = binary_search(&pvalue, pv, len, sizeof(double), double_bsearch_compar);
  double unscaled_score = (scaled_score >= 0) ? 
    pssm_scaled_to_raw(scaled_score, pssm) :
    pssm_scaled_to_raw(-scaled_score-1, pssm);
    
  // free the array
  myfree(pv);

  // return the unscaled score
  return(unscaled_score);
} // pssm_pvalue_to_score
