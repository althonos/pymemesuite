/**************************************************************************
 * FILE: mhmmscan.c
 * AUTHOR: William Stafford Noble, Timothy L. Bailey, Charles Grant
 * CREATE DATE: 5/21/02
 * PROJECT: MEME Suite
 * COPYRIGHT: 1998-2002, WNG, 2001, TLB
 * DESCRIPTION: Search a database of sequences using a motif-based
 * HMM.  Similar to mhmms, but allows arbitrarily long sequences and
 * multiple matches per sequence.
 **************************************************************************/

#ifdef MAIN
#define DEFINE_GLOBALS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "matrix.h"      // Routines for floating point matrices. 
#include "alphabet.h"
#include "array.h"       // Routines for floating point arrays. 
#include "cisml.h"       // CisML output format
#include "dir.h"         // Paths to MEME directories
#include "dp.h"          // Dynamic programming routines.
#include "fasta-io.h"    // Read FASTA format sequences. 
#include "fitevd.h"      // Extreme value distribution routines. 
#include "gendb.h"
#include "io.h"
#include "log-hmm.h"     // HMM log/log-odds conversion. 
#include "mhmmscan.h"
#include "metameme.h"    // Global metameme functions. 
#include "mhmm-state.h"  // HMM data structure. 
#include "mhmms.h"       // Framework for database search.
#include "pssm.h"        // Position-specific scoring matrix.
#include "rdb-matrix.h"  // For reading background files 
#include "read-mhmm.h"   // HMM input/output. 
#include "utils.h"       // Generic utilities. 

// Number of matches required before we generate synthetic scores if
// -synth given.
// If fewer matches than this are found, we assume that EM (or fitevd)
// will be unable to accurately estimate the distribution.
// The current value of this was guessed at by assuming that:
//   1) estimation works well if 10:1 ratio of negative:positive matches
//   2) rarely will there be more than 1000 true positives
#define MIN_MATCHES 10000

// Global variables
int       dp_rows = 0;           // Size of the DP matrix.
int       dp_cols = 0; 
MATRIX_T* dp_matrix = NULL;      // Dynamic programming matrix.
MATRIX_T* trace_matrix = NULL;   // Traceback for Viterbi.
MATCH_T*  complete_match = NULL; // All repeated matches.
MATCH_T*  partial_match = NULL;  // One partial match.

void init_mhmmscan_options(MHMMSCAN_OPTIONS_T *options) {
  options->allow_clobber = true;
  options->allow_weak_motifs = false;
  options->beta_set = false;
  options->both_strands = false;
  options->egcost_set = false;
  options->gap_extend_set = false;
  options->gap_open_set = false;
  options->use_obs_gc = false;
  options->max_gap_set = false;
  options->min_score_set = false;
  options->motif_scoring = false;
  options->use_synth = false;
  options->pam_distance_set = false;
  options->print_fancy = false;
  options->print_header = true;
  options->print_params = true;
  options->print_time = true;
  options->sort_output = true;
  options->text_only = false;
  options->use_pvalues = false;
  options->zero_spacer_emit_lo = false;

  options->bg_filename = NULL;
  options->command_line = NULL;
  options->output_dirname = NULL;
  options->gff_filename = NULL;
  options->hmm_filename = NULL;
  options->motif_filename = NULL;
  options->program = NULL;
  options->sc_filename = NULL;
  options->seq_filename = NULL;

  options->max_chars = -1;
  options->max_gap = -1;
  options->max_total_width = -1;
  options->max_seqs = NO_MAX_SEQS;
  options->output_width = DEFAULT_OUTPUT_WIDTH;
  options->pam_distance = -1;  // Use default PAM distance.
  options->progress_every = DEFAULT_PROGRESS_EVERY;

  options->beta = -1.0;
  options->egcost = 0.0; 
  options->gap_open = -1.0;    // No gap open penalty.
  options->gap_extend = -1.0;  // No gap extension penalty.
  options->motif_pthresh = -1.0; // Don't do p-value scoring.
  options->p_thresh = 1.0; // Default to allowing all matches
  options->q_thresh = 1.0; // Default to allowing all matches

  options->dp_thresh = DEFAULT_DP_THRESHOLD;
  options->e_thresh = DEFAULT_E_THRESHOLD;

  // Set up path values for output files.
  options->HTML_FILENAME = "mhmmscan.html";
  options->TEXT_FILENAME = "mhmmscan.txt";
  options->mhmmscan_filename = "mhmmscan.xml";
  options->CISML_FILENAME = "cisml.xml";
}

/***********************************************************************
  Free memory allocated in options processing
 ***********************************************************************/
void cleanup_mhmmscan_options(MHMMSCAN_OPTIONS_T *options) {
  alph_release(options->alphabet);
  myfree(options->command_line);
  myfree(options->mhmmscan_path);
  myfree(options->cisml_path);
  myfree(options->html_path);
  myfree(options->text_path);
  myfree(options->gff_path);
}

/***********************************************************************
 * Print MHMMSCAN specific information to an XML file
 ***********************************************************************/
void print_mhmmscan_xml_file(
  MHMMSCAN_OPTIONS_T *options,
  SCORE_SET *score_set,
  EVD_SET evd_set,
  ALPH_T* alph,
  int num_motifs,
  MOTIF_T *motifs,
  ARRAY_T *bgfreq,
  int num_seqs,
  long num_residues
) {
  int i;

  assert(options != NULL);

  FILE *out = fopen(options->mhmmscan_path, "w");
  if (!out) {
    die("Couldn't open file %s for output.\n", options->mhmmscan_path);
  }
  fputs("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n", out);
  fputs("<!-- Begin document body -->\n", out);
  fprintf(out, "<%s version=\"" VERSION "\" release=\"" ARCHIVE_DATE "\">\n",
      options->program);
  /*
  // TODO This namespace stuff has been in the XML since 4.8.0 but it
  // breaks XML parsing so I'm going to assume it was unintentional
  fputs("  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"", out);
  fputs("\n", out);
  fputs("  xsi:schemaLocation=", out);
  fputs("  xmlns:mhmmscan=\"http://noble.gs.washington.edu/schema/fimo\"\n>\n", out);
  */
  fprintf(out, "<command-line>%s</command-line>\n", options->command_line);
  fprintf(
    out, 
    "<sequence-data num-sequences=\"%d\" num-residues=\"%ld\" />\n", 
    num_seqs,
    num_residues
  );
  //TODO FIXME make work with other alphabets
  fprintf(
    out,
    "<alphabet>%s</alphabet>\n",
    alph_is_builtin_dna(options->alphabet) ? "nucleotide" : "protein"
  );
  for (i = 0; i < num_motifs; i++) {
    MOTIF_T *motif = motif_at(motifs, i);
    MOTIF_T *rc_motif = NULL;
    char *motif_id = get_motif_id(motif);
    char *rc_motif_id = NULL;
    if (i < (num_motifs - 1)) {
      rc_motif = motif_at(motifs, i + 1);
      rc_motif_id = get_motif_id(rc_motif);
    }
    char *best_possible_match = get_best_possible_match(motif);
    if (rc_motif_id && strcmp(motif_id, rc_motif_id) == 0) {
      ++i; // Pair of identiical motif ids indicate forward/reverse pair.
      char *best_possible_rc_match = get_best_possible_match(rc_motif);
      fprintf(
        out,
        "<motif name=\"\%s\" width=\"%d\" best-possible-match=\"%s\""
        " best-possible-rc-match=\"%s\"/>\n",
        motif_id,
        get_motif_length(motif),
        best_possible_match,
        best_possible_rc_match
      );
      myfree(best_possible_rc_match);
    }
    else {
      fprintf(
        out,
        "<motif name=\"\%s\" width=\"%d\" best-possible-match=\"%s\"/>\n",
        motif_id,
        get_motif_length(motif),
        best_possible_match
      );
    }
    myfree(best_possible_match);
  }
  fprintf(
    out,
    "<background source=\"%s\">\n",
    options->bg_filename ? options->bg_filename : "non-redundant database"
  );
  int alphabet_size = alph_size_core(alph);
  for (i = 0; i < alphabet_size; i++) {
    fprintf(
      out,
      "<value letter=\"%c\">%1.3f</value>\n",
      alph_char(alph, i),
      get_array_item(i, bgfreq)
    );
  }
  fputs("</background>\n", out);
  fputs("<cisml-file>cisml.xml</cisml-file>\n", out);
  fputs("<mcast-parameters>\n", out);
  fprintf(out, "<min-match-score>%1.3f</min-match-score>\n", options->dp_thresh);
  fprintf(out, "<e-thresh>%1.3f</e-thresh>\n", options->e_thresh);
  fprintf(out, "<p-thresh>%1.3f</p-thresh>\n", options->p_thresh);
  fprintf(out, "<q-thresh>%1.3f</q-thresh>\n", options->q_thresh);
  fprintf(out, "<max-gap>%d</max-gap>\n", options->max_gap);
  fprintf(out, "<max-total-width>%d</max-total-width>\n", options->max_total_width);
  fprintf(out, "<cost-factor>%.3g</cost-factor>\n", options->egcost);
  fprintf(out, "<gap-open-cost>%.3g</gap-open-cost>\n", options->gap_open);
  fprintf(out, "<gap-extend-cost>%.3g</gap-extend-cost>\n", options->gap_extend);
  fputs("</mcast-parameters>\n", out);
  fputs("<e-value-parameters>\n", out);
  fprintf(out, "<distribution>exponential</distribution>\n");
  if (evd_set.evds != NULL) {
    fprintf(out, "<status>%s</status>\n", evd_set.msg);
    fprintf(out, "<mu1>%.3g</mu1>\n", evd_set.evds[0].mu1);
    fprintf(out, "<mu2>%.3g</mu2>\n", evd_set.evds[0].mu2);
    fprintf(out, "<sigma2>%.3g</sigma2>\n", evd_set.evds[0].sigma2);
    fprintf(out, "<c>%.3g</c>\n", evd_set.evds[0].c);
    fprintf(out, "<num-scores-used>%d</num-scores-used>\n", evd_set.evds[0].n);
  }
  else {
    // Not enough scores to produce fit
    fputs("<status>To few scores to estimate fit</status>\n", out);
    fputs("<mu1>NaN</mu1>\n", out);
    fputs("<mu2>NaN</mu2>\n", out);
    fputs("<sigma2>NaN</sigma2>\n", out);
    fputs("<c>NaN</c>\n", out);
    fputs("<num-scores-used>0</num-scores-used>\n", out);
  }
  fprintf(out, "<num-matches>%d</num-matches>\n", score_set->n);
  fprintf(out, "<db-length>%g</db-length>\n", score_set->total_length);
  fputs("</e-value-parameters>\n", out);
  fprintf(out, "</%s>\n", options->program);

  fclose(out);
}

/**************************************************************************
 * Strip the filename off a path+filename string.
 * This function modifies the given string.
 **************************************************************************/
static char* strip_filename(char* filename) {
  int i;

  // Don't bother if it's empty.
  if (strlen(filename) == 0) {
    return(filename);
  }

  // Look for the last slash.
  for (i = strlen(filename) - 1; i >= 0; i--) {
    if (filename[i] == '/') {
      break;
    }
  }

  // Replace the slash with a '\0'.
  filename[i] = '\0';
  return(filename);
}

/**************************************************************************
 * Tell the user how far we are.
 **************************************************************************/
void user_update(
   int changed,
   int num_seqs,
   int num_segments,
   int num_matches,
   int num_stored,
   int progress_every
) {
  if (verbosity >= NORMAL_VERBOSE) {
    if ((changed != 0) && (changed % progress_every == 0)) {
      fprintf(stderr, "Sequences: %d Segments: %d Matches: %d Stored: %d\n",
              num_seqs, num_segments, num_matches, num_stored);
    }
  }
} // user_update

/**************************************************************************
 * read_and_score
 *
 * Read and score sequences from the input stream.
 * Set the number of sequences read.
 * Return the score_set.
 *
 **************************************************************************/
static SCORE_SET *read_and_score(
  ALPH_T*   alph,           // Alphabet
  SCORE_SET *score_set,     // Set of scores.
  FILE*     seq_file,       // Open stream to read from.
  SEQ_T     *sequence,      // Used if seq_file is NULL.
  bool negatives_only, // Stream contains only negative examples.
  int       max_chars,      // Number of chars to read at once.
  int       max_gap,        // Maximum gap length to allow in matches.
  MHMM_T*   log_hmm,        // The HMM, with probs converted to logs.
  bool motif_scoring,  // Perform motif-scoring.
  bool use_pvalues,    // Use p-value scoring?
  PROB_T    motif_pthresh, // P-value threshold for motif hits.
  PROB_T    dp_thresh,      // Score threshold for DP.
  PROB_T    e_thresh,    // E-value threshold for scores.
  bool got_evd,        // Have distribution?
  EVD_SET   evd_set,        // EVD data structure.
  bool print_fancy,    // Print alignments?
  bool store_gff,      // Store GFF entries?
  int       output_width,   // Width of output, in chars.
  int       progress_every, // Show progress after every n iterations.
  int       *num_seqs       // Number of sequences read.
) {
  int num_segments = 0; // Number of sequence segments processed.
  int num_matches = 0;  // Number of matches found.
  int start_pos = 0;    // Search from here.

  // Allocate memory.
  complete_match = allocate_match();
  partial_match = allocate_match();

  // Set up for computing score distribution.
  if (score_set == NULL) {
    score_set = set_up_score_set(
      motif_pthresh, 
      dp_thresh, 
      max_gap, 
      negatives_only,
      log_hmm);
  }

  *num_seqs = 0;

  while (
    (seq_file && read_one_fasta_segment(alph, seq_file, max_chars, &sequence)) 
    || (sequence != NULL)
  ) {
    // Number of motifs x sequence length.
    MATRIX_T* motif_score_matrix = NULL;

    // Keep track of total database size for E-value calculation.
    score_set->total_length += get_seq_length(sequence);

    // Let the user know what's going on.
    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Scoring %s (length=%d) at position %d.\n",
        get_seq_name(sequence), get_seq_length(sequence),
        start_pos);
    }

    // Convert the sequence to alphabet-specific indices. 
    prepare_sequence(sequence, alph, false /* no hard masking */);

    /* Allocate the dynamic programming matrix. Rows correspond to
       states in the model, columns to positions in the sequence. */
    if ((dp_rows < log_hmm->num_states) 
      || (dp_cols < get_seq_length(sequence))
    ) {
      free_matrix(dp_matrix);
      free_matrix(trace_matrix);
      if (dp_rows < log_hmm->num_states)
        dp_rows = log_hmm->num_states;
      if (dp_cols < get_seq_length(sequence))
        dp_cols = get_seq_length(sequence);
      // (Add one column for repeated match algorithm.)
      dp_matrix = allocate_matrix(dp_rows, dp_cols + 1);
      trace_matrix = allocate_matrix(dp_rows, dp_cols + 1);
    }

    // Compute the motif scoring matrix.
    if (motif_scoring) {
      motif_score_matrix = allocate_matrix(
        log_hmm->num_motifs,
        get_seq_length(sequence)
      );
      compute_motif_score_matrix(
        use_pvalues,
        motif_pthresh,
        get_int_sequence(sequence),
        get_seq_length(sequence),
        NULL, /* Array of priors (no priors) */
        0, /* Count of priors (no priors) */
        1.0, /* Alpha parameter for calculating log-ratio priors (no priors) */
        log_hmm,
        &motif_score_matrix
      );
    }

    // Fill in the DP matrix.
    repeated_match_algorithm(
      dp_thresh,
      get_int_sequence(sequence),
      get_seq_length(sequence),
      log_hmm,
      motif_score_matrix,
      dp_matrix,
      trace_matrix,
      complete_match
    );

    // Find all matches in the matrix.
    while (find_next_match(
      is_complete(sequence),
      start_pos,
      dp_matrix, 
      log_hmm,
      complete_match,
      partial_match)
    ) {

      // If this match starts in the overlap region, put if off until
      // the next segment.
      if (!is_complete(sequence) 
        && (get_start_match(partial_match) 
          > (get_seq_length(sequence) - OVERLAP_SIZE))
      ) {
        break;
      }

      // If this match starts at the beginning, then it's part of a
      // longer match from the previous matrix.  Skip it.  
      if (get_start_match(partial_match) == start_pos) {
        // fprintf(stderr, "Skipping match at position %d.\n", start_pos);

        // Change the starting position so that next call to find_next_match
        // will start to right of last match.
        start_pos = get_end_match(partial_match) + 1;
        continue;
      }
      start_pos = get_end_match(partial_match) + 1;

      // Store the score, ID, length and comment for later printing. 
      store_sequence(
        alph,
        motif_scoring, 
        true, // mhmmscan format 
        output_width,
        sequence,
        partial_match, 
        e_thresh,
        dp_thresh,
        motif_pthresh,
        score_set,
        got_evd,
        evd_set,
        print_fancy,
        log_hmm,
        motif_score_matrix,
        store_gff,
        NULL
      );

      /* Calculate the initial E-value distribution if the required
       * number sequences has been saved.  This will allow the
       * descriptions of low-scoring sequences not to be stored.  The
       * distribution will be recomputed using all scores when all
       * sequences have been read.  */
      if (score_set->n == EVD_NUM_SEQS && got_evd == false) {
        evd_set = calc_distr(
          *score_set,                   // Set of scores.  
          use_pvalues ? D_EXP : D_EVD,  // Use exponential distribution?
          true                          // Use match E-values.
        );
        if (evd_set.n > 0) {
          got_evd = true;
        }
      }

      num_matches++;
      user_update(
        num_matches, 
        *num_seqs, 
        num_segments, 
        num_matches, 
        get_num_stored(),
        progress_every
      );
    }

    // Does the input file contain more of this sequence?
    if (!is_complete(sequence)) { 
      int size_to_remove;

      // Compute the width of the buffer.
      size_to_remove = get_seq_length(sequence) - OVERLAP_SIZE;

      // Adjust the sequence accordingly.
      remove_flanking_xs(sequence);
      shift_sequence(size_to_remove, sequence);
      if (SHIFT_DEBUG) {
        fprintf(stderr, "Retained %d bases.\n", get_seq_length(sequence));
      }

      // Remove the left-over portion from the total database size.
      score_set->total_length -= get_seq_length(sequence);

      // Next time, start looking for a match after the buffer.
      start_pos -= size_to_remove;
      if (start_pos < 0) {
        start_pos = 0;
      }

    } else {

      // Free the memory used by this sequence. 
      free_seq(sequence);
      sequence = NULL;
      start_pos = 0;

      *num_seqs += 1;
      user_update(
        *num_seqs, 
        *num_seqs, 
        num_segments, 
        num_matches, 
        get_num_stored(),
        progress_every
      );
    }

    // Free the motif score matrix.
    free_matrix(motif_score_matrix);
    num_segments++;
    user_update(
      num_segments, 
      *num_seqs, 
      num_segments, 
      num_matches,
      get_num_stored(),
      progress_every
    );

    if (seq_file == NULL) {
      break;
    }
  }

  //
  // Free structures.
  //
  if (0) { // FIXME
    free_matrix(dp_matrix);
    free_matrix(trace_matrix);
    free_match(complete_match);
    free_match(partial_match);
  }

  // Return the set of scores.
  return(score_set);
} // read_and_score

/**************************************************************************
* Generate synthetic sequences and score them if fewer than
* MIN_MATCHES matches was found.  Don't do this if no matches
* were found (of course).
*
* Each sequence is 100000bp long.
* Continues to generate sequences until enough matches are found
* so that the standard error of mu1 will be 5% (std err = 100 * 1/sqrt(n))
* or the maximum database size is reached.
***************************************************************************/
SCORE_SET *generate_synth_seq(
  SCORE_SET *score_set,
  MHMM_T *log_hmm,
  bool got_evd,
  EVD_SET evd_set,
  char *bg_filename,
  int max_chars,
  int max_gap,
  bool motif_scoring,
  bool use_pvalues,
  bool use_obs_gc,
  PROB_T motif_pthresh,
  PROB_T dp_thresh,
  int output_width,
  int progress_every
) {

  SCORE_SET *synth_score_set = NULL;
  int i;
  int want = 1e3;  // Desired number of matches (for mu1 std err= 3%).
  double need = MIN_SCORES;  // Minimum number of matches (for mu1 std err= 5%).
  int found = 0;  // Matches found.
  double maxbp = 1e7;  // Maximum number of bp to generate.
  int db_size = 0;  // Number of bp generated. 
  int len = 1e5;  // Desired sequence length.
  int nseqs = 0;  // Number of sequences generated.
  double min_gc, max_gc;   // Min and Max observed GC-contents.
  int x1, x2, y1, y2; // the indexes of the complementary pairs
  ALPH_T *alph = log_hmm->alph; // Alphabet

  // Get the minimum and maximum GC-contents of match regions.
  min_gc = 1;
  max_gc = 0;
  if (use_obs_gc) {
    assert(alph_size_pairs(alph) == 2);
    for (i=0; i<score_set->n; i++) {
      double gc = score_set->scores[i].gc;
      if (gc < min_gc) min_gc = gc;
      if (gc > max_gc) max_gc = gc;
    }
    // determine alphabet indexes for the complementary pairs
    x1 = 0;
    x2 = alph_complement(alph, x1);
    y1 = (x2 == 1 ? 2 : 1);
    y2 = alph_complement(alph, y1);
  } else {
    x1 = x2 = y1 = y2 = 0; // stop compiler complaining
  }

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Generating and scoring synthetic sequences...\n");
  }

  // Local variables. 
  double    start_time;  // Time at start of sequence processing. 
  double    end_time;    // Time at end of sequence processing. 
  int       num_seqs;    // Number of sequences processed. 
  int       dummy;       // Throw away variable.
  bool use_synth;   // Distribution based on synthetic scores.
  FILE*     out_stream;  // Output stream (possibly running mhmm2html).
  ARRAY_T *bg;
  mt_state prng;

  // set PRNG to constant seed
  mts_seed32new(&prng, 1234);

  // get a background
  if (bg_filename != NULL) {
    // load bg from file
    bg = load_markov_model(alph, NULL, bg_filename);
  } else {
    // get NRDB or uniform freqs
    bg = allocate_array(alph_size_core(alph));
  }

  // Loop until enough matches found or we get tired.
  while (found < want && db_size < maxbp) {  
    FILE* synth_seq_file = tmpfile();    // Use file not seq variable.
    if (use_obs_gc) {
      double gc, a_or_t, g_or_c;
      // Select a GC content randomly and set bg
      gc = min_gc + (mts_ldrand(&prng) * (max_gc - min_gc));
      a_or_t = (1 - gc)/2;
      g_or_c = gc/2;
      set_array_item(x1, a_or_t, bg); 
      set_array_item(x2, a_or_t, bg); 
      set_array_item(y1, g_or_c, bg); 
      set_array_item(y2, g_or_c, bg); 
      fprintf(stderr, "gc = %f\r", gc);
    }
    gendb(
      synth_seq_file, // output
      &prng, // random number generator
      alph, // alphabet
      bg, // background
      0.0, // no ambiguous symbols
      1, // # of sequences
      len, // min length
      len // max length
    );

    nseqs++; // seed
    db_size += len; // Size of db so far.
    if (synth_seq_file) {
      rewind(synth_seq_file);
    }
    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "seqs: %d size: %d\r", nseqs, db_size);
    }

    //
    // Read and score synthetic sequences.
    //
    int dummy; // Throw away variable.
    synth_score_set = read_and_score(
      alph,
      synth_score_set,
      synth_seq_file,
      NULL,
      true, // Negatives only.
      max_chars,
      max_gap,
      log_hmm,
      motif_scoring,
      use_pvalues,
      motif_pthresh,
      dp_thresh,
      -1, // e_thresh: save nothing
      got_evd,
      evd_set,
      false,
      false,
      output_width,
      progress_every, 
      &dummy
    );
    found = synth_score_set->n;
    fprintf(stderr, "seqs: %d size: %d matches: %d\r", nseqs, db_size, found);

    // Close file to delete sequences.
    if (synth_seq_file) {
      fclose(synth_seq_file);
    }
    // free_seq(seq);      // free'd by read_and_score

    // Quit if we're not getting there.
    if (found/need < db_size/maxbp) {
      fprintf(
        stderr, 
        "\nGiving up generating synthetic sequences: match probability too low!\n"
      );
      break;
    }
  } // Generate synthetic scores.

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(
      stderr, 
      "Generated %d synthetic sequences (%d characters) with %d matches.\n",
      nseqs, 
      db_size, 
      found
    );
  }

  return synth_score_set;
}

/******************************************************************************
 * Run the mhmmscan algorithm
******************************************************************************/
void mhmmscan(
  int argc, 
  char* argv[], 
  MHMMSCAN_OPTIONS_T *options, 
  MHMM_T *hmm, 
  FILE* seq_file
) {

  MHMM_T*   log_hmm;             // The HMM, with probs converted to logs. 
  EVD_SET   evd_set;                 // EVD data structure 
  bool got_evd = false;         // no EVD found yet 
  SCORE_SET *score_set = NULL;       // Set of scores for computing distribution.
  SCORE_SET *synth_score_set = NULL; // Synthetic sequence scores.

  ALPH_T*   alph = hmm->alph;        // Alphabet
  double    start_time;  // Time at start of sequence processing. 
  double    end_time;    // Time at end of sequence processing. 
  int       num_seqs;    // Number of sequences processed. 
  FILE*     out_stream;  // Output stream (possibly running mhmm2html).

  // Record CPU time. 
  myclock();

  FILE*     gff_file;
  if (options->gff_filename != NULL) {
    if (open_file(options->gff_filename, "w", false, "gff", "gff", &gff_file) == 0) {
      exit(1);
    }
  } else {
    gff_file = NULL;
  }

  alph = hmm->alph;
  if (!(alph_is_builtin_dna(alph) || alph_is_builtin_protein(alph))) 
    die("%s alphabet is unsupported.\n", alph_name(alph));

  //
  // Update options with info from the HMM.
  //
  if (options->max_gap_set) {

    options->dp_thresh = 1e-6;  // Very small number.

    // 
    // If using egcost, we must get model statistics
    // in order to set dp_thresh.
    //
    if (options->egcost > 0.0) {    // Gap cost a fraction of expected hit score.
      score_set = set_up_score_set(
        options->motif_pthresh, 
        options->dp_thresh, 
        options->max_gap, 
        false,
        hmm
      );
      options->dp_thresh = 
        (options->egcost * options->max_gap * score_set->e_hit_score) / score_set->egap;
    } 

    //
    // Set gap costs.
    //
    options->gap_open = options->gap_extend = options->dp_thresh/options->max_gap;
    options->zero_spacer_emit_lo = true;
  }

  //
  // Prepare the model for recognition.
  //
  if (options->pam_distance == -1) {
    options->pam_distance = alph_is_builtin_protein(alph) ? DEFAULT_PROTEIN_PAM : DEFAULT_DNA_PAM;
  }
  if (!options->beta_set) {
    options->beta = alph_is_builtin_dna(alph) ? DEFAULT_PROTEIN_BETA : DEFAULT_DNA_BETA;
  }
  free_array(hmm->background);
  hmm->background = get_background(alph, options->bg_filename);
  log_hmm = allocate_mhmm(hmm->alph, hmm->num_states);
  convert_to_from_log_hmm(
    true, // Convert to logs.
    options->zero_spacer_emit_lo,
    options->gap_open,
    options->gap_extend,
    hmm->background,
    options->sc_filename,
    options->pam_distance,
    options->beta,
    hmm, 
    log_hmm
  );

  // Set up PSSM matrices if doing motif_scoring
  // and pre-compute motif p-values if using p-values.
  // Set up the hot_states list.
  set_up_pssms_and_pvalues(
    options->motif_scoring,
    options->motif_pthresh,
    options->both_strands,
    options->allow_weak_motifs,
    log_hmm,
    NULL, // No priors
    1.0 // Alpha not used
  );

  // Decide how many characters to read.
  if (options->max_chars == -1) {
    options->max_chars = (int)(MAX_MATRIX_SIZE / log_hmm->num_states);
  }
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(
      stderr, 
      "Reading sequences in blocks of %d characters.\n", 
      options->max_chars
    );
  }

  start_time = myclock();

  //
  // Read and score the real sequences.
  //
  evd_set.n = 0;			// prevent compiler warning
  score_set = read_and_score(
    alph,
    NULL,      // No score set yet.
    seq_file,
    NULL,      // No sequence.
    false,      // Positives and negatives.
    options->max_chars,
    options->max_gap,
    log_hmm,
    options->motif_scoring,
    options->use_pvalues,
    options->motif_pthresh,
    options->dp_thresh,
    options->e_thresh,
    got_evd,
    evd_set,
    options->print_fancy,
    gff_file != NULL,
    options->output_width,
    options->progress_every, 
    &num_seqs
  );

  if (options->use_synth == true && score_set->n>0 && score_set->n<MIN_MATCHES) {  
    synth_score_set = generate_synth_seq(
        score_set,
        log_hmm,
        got_evd,
        evd_set,
        options->bg_filename,
        options->max_chars,
        options->max_gap,
        options->motif_scoring,
        options->use_obs_gc,
        options->use_pvalues,
        options->motif_pthresh,
        options->dp_thresh,
        options->output_width,
        options->progress_every
    );
  } else {
    synth_score_set = score_set;    // Use actual scores for estimation.
  } // Generate and score synthetic sequences.

  end_time = myclock();

  /***********************************************
   * Calculate the E-values and store them as the keys.
   ***********************************************/
  // Recalculate the score distribution using synthetic scores. 
  // If successful, calculate the E-values and store them as keys.
  evd_set = calc_distr(
    *synth_score_set,                     // Set of scores.
    options->use_pvalues ? D_EXP : D_EVD, // Use exponential distribution?
    true                                  // Use match E-values.
  );
  // Distribution based on synthetic scores?
  evd_set.negatives_only = options->use_synth; 
  if (evd_set.n >= 1) {      // Found a valid score distribution.
    int q, t, N;
    q = 1; // Ignore query length.
    // Get p-value multiplier.
    if (options->use_synth) {
      N = evd_set.non_outliers = get_n(*score_set, evd_set);
    } else {
      N = evd_set.non_outliers; // Use number of non-outliers.
    }
    // p-value multiplier for E-value.
    evd_set.N = N;      
    // Record number of real scores in evd_set.
    evd_set.nscores = score_set->n;  

    // Get sequence length.
    if (options->use_pvalues) { // Exponential.
      t = 0; // Ignore sequence length.
    } else { // EVD
      t = score_set->total_length/score_set->n;  // average length
    }

    calc_evalues(&evd_set, N, q, t);
    got_evd = true;
  }

  /***********************************************
   * Start the mhmm2html process, if requested.
   ***********************************************/
  if (options->text_only) {
    out_stream = stdout;
  } else {
    out_stream = open_command_pipe(
      "mhmm2html", // Program
      strip_filename(argv[0]), // Search directory
      "-test -", // Test arguments.
      "mhmm2html", // Expected reply.
      "-", // Real arguments.
      true, // Return stdout on failure.
      "Warning: mhmm2html not found.  Producing text output."
    );
  }

  /***********************************************
   * Print header information.
   ***********************************************/
  if (options->print_header) {
    write_header(
      "mhmmscan",
      "Database search results",
      hmm->description,
      hmm->motif_file,
      options->hmm_filename,
      options->seq_filename,
      out_stream
    );
    if (gff_file != NULL) {
      write_header(
        "mhmmscan",
        "Database search results",
        hmm->description,
        hmm->motif_file,
        options->hmm_filename,
        options->seq_filename,
        gff_file
      );
    }
  }

  /***********************************************
   * Sort and print the results.
   ***********************************************/
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "\nSorting the scores.\n");
  }
  sort_and_print_scores(
    options->print_fancy, 
    options->print_header, 
    got_evd,
    options->motif_scoring,
    true, // Print in mhmmscan format.
    options->max_seqs, // Maximum number of sequences to print
    options->output_width,
    options->e_thresh,
    options->sort_output,
    gff_file,
    out_stream
  );

  if (options->print_params) {
    print_parameters(
      argv, 
      argc,
      "mhmmscan",
      alph,
      options->hmm_filename,
      options->seq_filename,
      true, // Viterbi search.
      options->dp_thresh,
      options->motif_scoring,
      options->use_pvalues,
      options->motif_pthresh,
      options->e_thresh,
      options->both_strands,
      options->bg_filename,
      options->sc_filename,
      options->pam_distance,
      options->beta,
      options->zero_spacer_emit_lo,
      options->max_gap,
      options->egcost,
      options->gap_open,
      options->gap_extend,
      options->print_fancy,
      options->print_time,
      start_time,
      end_time,
      num_seqs,
      evd_set,
      *score_set,
      out_stream
    );
  }

  // Tie up loose ends.
  free_mhmm(log_hmm);
  if (options->gff_filename != NULL) {
    fclose(gff_file);
  }    
  pclose(out_stream);

}

#ifdef MAIN

VERBOSE_T verbosity = NORMAL_VERBOSE;

#include "simple-getopt.h"  // Cmd line processing

/***********************************************************************
  Process command line options
 ***********************************************************************/
void process_command_line(
  int argc,
  char* argv[],
  MHMMSCAN_OPTIONS_T *options
) {

    // Define command line options.
    cmdoption const cmd_options[] = {
      {"allow-weak-motifs", NO_VALUE},
      {"bg-file", REQUIRED_VALUE},
      {"blocksize", REQUIRED_VALUE},
      {"e-thresh", REQUIRED_VALUE},
      {"eg-cost", REQUIRED_VALUE},
      {"fancy", NO_VALUE},
      {"gap-open", REQUIRED_VALUE},
      {"gap-extend", REQUIRED_VALUE},
      {"gff", REQUIRED_VALUE},
      {"maxseqs", REQUIRED_VALUE},
      {"max-gap", REQUIRED_VALUE},
      {"min-score", REQUIRED_VALUE},
      {"motif-scoring", NO_VALUE},
      {"noheader", NO_VALUE},
      {"noparams", NO_VALUE},
      {"nosort", NO_VALUE},
      {"notime", NO_VALUE},
      {"p-thresh", REQUIRED_VALUE},
      {"pam", REQUIRED_VALUE},
      {"progress", REQUIRED_VALUE},
      {"pseudo-weight", REQUIRED_VALUE},
      {"quiet", NO_VALUE},
      {"score-file", REQUIRED_VALUE},
      {"use-obs-gc-synth", NO_VALUE},
      {"synth", NO_VALUE},
      {"text", NO_VALUE},
      {"verbosity", REQUIRED_VALUE},
      {"width", REQUIRED_VALUE},
      {"zselo", NO_VALUE}
    };
    const int num_options = sizeof(cmd_options) / sizeof(cmdoption);
    int option_index = 0;

    // Define the usage message.
    char *usage = 
      "Usage: mhmmscan [options] <HMM> <database>\n"
      "\n"
      "   Options:\n"
      "     --allow-weak-motifs\n"
      "     --bg-file <file>\n"
      "     --blocksize <int>\n"
      "     --e-thresh <E-value> (default=10.0)\n"
      "     --fancy\n"
      "     --gff <file>\n"
      "     --max-gap <int>\n"
      "     --maxseqs <int>\n"
      "     --noheader\n"
      "     --noparams\n"
      "     --nosort\n"
      "     --notime\n"
      "     --p-thresh <p-value>\n"
      "     --progress <int>\n"
      "     --quiet\n"
      "     --use-obs-gc-synth\n"
      "     --synth\n"
      "     --text\n"
      "     --verbosity 1|2|3|4|5 (default=2)\n"
      "     --width <int> (default=79)\n"
      "\n"
      "   Advanced options:\n"
      "     --eg-cost <fraction>\n"
      "     --gap-extend <cost>\n"
      "     --gap-open <cost>\n"
      "     --min-score <score>\n"
      "     --motif-scoring\n"
      "     --pam <distance> (default=250 [protein] 1 [DNA])\n"
      "     --pseudo-weight <weight> (default=10)\n"
      "     --score-file <file>\n"
      "     --zselo\n"
      "\n";

    simple_setopt(argc, argv, num_options, cmd_options);
       
    // Parse the command line.
    while (1) { 
      int c = 0;
      char* option_name = NULL;
      char* option_value = NULL;
      const char* message = NULL;

      // Read the next option, and break if we're done.
      c = simple_getopt(&option_name, &option_value, &option_index);
      if (c == 0) {
        break;
      } else if (c < 0) {
        simple_getopterror(&message);
        die("Error process command line options (%s)\n", message);
      }

      if (strcmp(option_name, "allow-weak-motifs") == 0) {
        options->allow_weak_motifs = true;
      } else if (strcmp(option_name, "bg-file") == 0) {
        options->bg_filename = option_value;
      } else if (strcmp(option_name, "blocksize") == 0) {
        options->max_chars = atoi(option_value);
      } else if (strcmp(option_name, "e-thresh") == 0) {
        options->e_thresh = atof(option_value);
      } else if (strcmp(option_name, "eg-cost") == 0) {
        options->egcost_set = true;
        options->egcost = atof(option_value);
      } else if (strcmp(option_name, "fancy") == 0) {
        options->print_fancy = true;
      } else if (strcmp(option_name, "gap-extend") == 0) {
        options->gap_extend_set = true;
        options->gap_extend = atof(option_value);
      } else if (strcmp(option_name, "gap-open") == 0) {
        options->gap_open_set = true;
        options->gap_open = atof(option_value);
      } else if (strcmp(option_name, "gff") == 0) {
        options->gff_filename = option_value;
      } else if (strcmp(option_name, "max-gap") == 0) {
        options->max_gap_set = true;
        options->max_gap = atoi(option_value);
        // use of max-gap option forces --zselo
        options->zero_spacer_emit_lo = true;
      } else if (strcmp(option_name, "maxseqs") == 0) {
        options->max_seqs = atoi(option_value);
      } else if (strcmp(option_name, "min-score") == 0) {
        options->min_score_set = true;
        options->dp_thresh = atof(option_value);
      } else if (strcmp(option_name, "motif-scoring") == 0) {
        options->motif_scoring = true;
      } else if (strcmp(option_name, "noheader") == 0) {
        options->print_header = false;
      } else if (strcmp(option_name, "noparams") == 0) {
        options->print_params = false;
      } else if (strcmp(option_name, "nosort") == 0) {
        options->sort_output = false;
      } else if (strcmp(option_name, "notime") == 0) {
        options->print_time = false;
      } else if (strcmp(option_name, "p-thresh") == 0) {
        options->use_pvalues = true;
        // setting p-threshold forces motif_scoring.
        options->motif_scoring = true;
        options->motif_pthresh = atof(option_value);
      } else if (strcmp(option_name, "pam") == 0) {
        options->pam_distance_set = true;
        options->pam_distance = atoi(option_value);
      } else if (strcmp(option_name, "progress") == 0) {
        options->progress_every = atoi(option_value);
      } else if (strcmp(option_name, "pseudo-weight") == 0) {
        options->beta_set = true;
        options->beta = atof(option_value);
      } else if (strcmp(option_name, "quiet") == 0) {
        options->print_header = options->print_params = options->print_time = false;
        verbosity = QUIET_VERBOSE;
      } else if (strcmp(option_name, "score-file") == 0) {
        options->sc_filename = option_value;
      } else if (strcmp(option_name, "use-obs-gc-synth") == 0) {
        options->use_obs_gc = true;
      } else if (strcmp(option_name, "synth") == 0) {
        options->use_synth = true;
      } else if (strcmp(option_name, "text") == 0) {
        options->text_only = true;
      } else if (strcmp(option_name, "verbosity") == 0) {
        verbosity = (VERBOSE_T)atoi(option_value);
      } else if (strcmp(option_name, "width") == 0) {
        options->output_width = atoi(option_value);
      } else if (strcmp(option_name, "zselo") == 0) {
        options->zero_spacer_emit_lo = true;
      } 
    }

    // Read the two required arguments.
    if (option_index + 2 != argc) {
      fprintf(stderr, "%s", usage);
      exit(1);
    }
    options->hmm_filename = argv[option_index];
    options->seq_filename = argv[option_index+1];
     
    /***********************************************
     * Verify the command line.
     ***********************************************/

    // Make sure we got the required files. 
    if (options->hmm_filename == NULL) {
      die("No HMM file given.\n");
    }
    if (options->seq_filename == NULL) {
      die("No sequence file given.\n");
    }

    // Check p-threshold is in range [0<p<=1].
    if (options->use_pvalues && (options->motif_pthresh <= 0 || options->motif_pthresh > 1)) {
      die("You may only specify p-thresh in the range [0<p<=1].\n");
    }

    if (options->max_gap_set) {

      // When max-gap option is given, don't allow the gap-open,
      // gap-extend or min-score options.
      if(
        options->gap_open_set
        || options->gap_extend_set
        || options->min_score_set
      ) {
        die("You may not use the gap-extend, gap-open, or min-score options\n"
            "in combination with the max-gap option.\n");
      } 

      // Check that length is legal.
      if (options->max_gap < 0) {
        die("You may not specify a negative max-gap length.\n");
      }

    }

    // Check that min-score is positive.
    if (options->min_score_set && options->dp_thresh <= 0.0) {
      die("You may only specify a positive min-score.\n");
    }

    // FIXME: both-strands not implemented
    if (options->both_strands) {
      die("Sorry, -both-strands not yet implemented.");
    }

    // Check egcost ok.
    if (options->egcost_set && options->min_score_set) {
      die("You may not use both the egcost and min-score options.\n");
    }

    if (options->egcost_set && ! options->max_gap_set) {
      die("If you use the egcost option, you must also use the\n"
          "max-gap options.\n"
      );
    }

    if (options->egcost_set && options->egcost <= 0.0) {
      die("The egcost must be positive.\n");
    }

    // Check that -bg-file given if -synth given
    if (options->use_synth == true && options->bg_filename == NULL) {
      die("You may only use -synth if you specify -bg-file.\n");
    }

}

int main(int argc, char * argv[]) {

  MHMMSCAN_OPTIONS_T options;
  FILE* seq_file;

  // Set option defaults. 
  init_mhmmscan_options(&options);

  // Parse the command line into mhmmscan options.
  process_command_line(argc, argv, &options);

  // Use the sequence file name to get sequence file handle.
  if (open_file(options.seq_filename, "r", true, "sequence", "sequences", &seq_file) 
    == 0) {
    exit(1);
  }

  // Read the HMM. 
  MHMM_T* hmm;
  read_mhmm(options.hmm_filename, &hmm);

  // Scan the sequences with the HMM.
  mhmmscan(argc, argv, &options, hmm, seq_file);

  return(0);
}
#endif // MAIN
