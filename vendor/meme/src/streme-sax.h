#ifndef STREME_SAX_H
#define STREME_SAX_H

#include <libxml/parser.h>

#include <stdarg.h>

/*****************************************************************************
 * Method for obtaining negative dataset enumeration
 ****************************************************************************/
enum streme_neg {
  STREME_NEG_SHUFFLED, // apply a shuffle to the positive set
  STREME_NEG_FILE,     // use a supplied file for the negative dataset
  STREME_NEG_NONE      // CD only: use no negative dataset
};
typedef enum streme_neg STREME_NEG_EN;

/*****************************************************************************
 * Strands
 ****************************************************************************/
enum streme_strands {
  STREME_STRANDS_BOTH, 		// use both strands
  STREME_STRANDS_FORWARD,	// given strand
  STREME_STRANDS_NONE		// single strand alphabet
};
typedef enum streme_strands STREME_STRANDS_EN;

/*****************************************************************************
 * Objective function 
 ****************************************************************************/
enum streme_objfun {
  STREME_DIFFERENTIAL_ENRICHMENT,
  STREME_CENTRAL_DISTANCE
};
typedef enum streme_objfun STREME_OBJFUN_EN;

/*****************************************************************************
 * Test 
 ****************************************************************************/
enum streme_test {
  STREME_FISHER_EXACT_TEST,
  STREME_BINOMIAL_TEST,
  STREME_CUMULATIVE_BATES_DISTRIBUTION
};
typedef enum streme_test STREME_TEST_EN;

/*****************************************************************************
 * Background sources enumeration
 ****************************************************************************/
enum streme_bg {
  STREME_BG_FROM_DATASET, // use negative dataset frequencies
  STREME_BG_FROM_FILE  // use pre-calculated frequencies from a file
};
typedef enum streme_bg STREME_BG_EN;

struct streme_io_xml_callbacks {

/*****************************************************************************
 * Error message handler
 ****************************************************************************/
  void (*error)(void * ctx, const char *format, va_list args);

/*****************************************************************************
 * streme
 * contains: model, motifs, reason_for_stopping, run_time
 *
 *  release         the release date.
 *  version         the program version.
 ****************************************************************************/
  void (*start_streme)(void *ctx, int major_version, int minor_version, 
      int patch_version, char *release_date);

/*****************************************************************************
 * /streme
 ****************************************************************************/
  void (*end_streme)(void *ctx);

/*****************************************************************************
 * STREME > model
 * contains: command_line, train_positives, train_negatives, test_positives,
 *          test_negatives, alphabet, strands, sequence_db, background_frequencies,
 *          stop, objfun, srch, test, minw, maxw, kmer, hofract, neval, nref,
 *          niter, mmr, seed, host
 ****************************************************************************/
  void (*start_model)(void *ctx);

/*****************************************************************************
 * STREME > /model
 ****************************************************************************/
  void (*end_model)(void *ctx);

/*****************************************************************************
 * STREME > model > command_line
 * the command-line used to run the program.
 ****************************************************************************/
  void (*handle_command_line)(void *ctx, char *command_line);

/*****************************************************************************
 * STREME > model > train_positives
 *
 *  count           the number of sequences in the positive input set
 *  positions       the number of characters in the positive training set
 *  file            the file containing the positive input set
 ****************************************************************************/
  void (*handle_train_positives) (void *ctx, long sequence_count, long sequence_positions,
      char *positives_file_path);

/*****************************************************************************
 * STREME > model > train_negatives
 *
 *  count           the number of sequences in the negative dataset
 *  positions       the number of characters in the negative training set
 *  from            the source of the negative dataset (eg shuffled positives)
 *  file            the file containing the negative dataset (optional)
 ****************************************************************************/
  void (*handle_train_negatives) (void *ctx, long sequence_count, long sequence_positions,
      STREME_NEG_EN negatives_source, char *negatives_file_path);

/*****************************************************************************
 * STREME > model > test_positives
 *
 *  count           the number of sequences in the positive training set
 *  positions       the number of characters in the positive training set
 ****************************************************************************/
  void (*handle_test_positives) (void *ctx, long sequence_count, long sequence_positions);

/*****************************************************************************
 * STREME > model > test_negatives
 *
 *  count           the number of sequences in the negative training set
 *  positions       the number of characters in the negative training set
 ****************************************************************************/
  void (*handle_test_negatives) (void *ctx, long sequence_count, long sequence_positions);

/*****************************************************************************
 * STREME > model > sequence_db
 *
 *  This is a dummy entry with no contents.
 ****************************************************************************/
  void (*handle_sequence_db) (void *ctx, int dummy);

/*****************************************************************************
 * STREME > model > alphabet
 *
 *  name            is the name of the alphabet
 *  extends_flags   is the alphabet being extended or 0
 ****************************************************************************/
  void (*start_alphabet) (void *ctx, char *name, int extends_flags);

/*****************************************************************************
 * STREME > model > /alphabet
 *
 *  name            is the name of the alphabet
 ****************************************************************************/
  void (*end_alphabet) (void *ctx);

/*****************************************************************************
 * STREME > model > alphabet > letter
 *
 *  id              the identifier used to represent the symbol as XML attributes
 *  symbol          the symbol
 *  aliases         additional symbols that resolve to the same meaning
 *  complement      the complement symbol in a double stranded alphabet
 *  equals          the comprising symbols
 *  name            the name of the negative dataset
 *  colour          the colour in RGB triplet (first 3 bytes)
 ****************************************************************************/
  void (*handle_alphabet_letter)(void *ctx, char *id, char symbol,
      char *aliases, char complement, char *equals, char *name, int color);

/*****************************************************************************
 * STREME > model > strands
 *
 *  strands         how is streme handling the strands
 ****************************************************************************/
  void (*handle_strands) (void *ctx, STREME_STRANDS_EN strands);

/*****************************************************************************
 * STREME > model > background
 *
 *  nfreqs          is the count of frequencies (same as alphabet core size)
 *  freqs           is the array of frequencies (ordered as in alphabet core)
 *  from            from the positive dataset or a background file
 ****************************************************************************/
  void (*handle_background) (void *ctx, int nfreqs, double* freqs, 
      STREME_BG_EN background_source);

  void (*start_background_frequencies)(void *, char *); // source of data
  void (*end_background_frequencies)(void *);
  //start background frequencies
  void (*start_bf_alphabet_array)(void *);
  void (*end_bf_alphabet_array)(void *);
  //start alphabet_array
  void (*handle_bf_aa_value)(void *, char *, double);
  //end alphabet_array
  //end background frequencies

/*****************************************************************************
 * STREME > model > stop
 *
 *  pvt         the maximum p-value threshold for motifs
 *  nmotifs     the maximum number of motifs to find
 *  time        the maximum time to run (seconds)
 ****************************************************************************/
  void (*handle_stop) (void *ctx, double pvt, int nmotifs, double time);

/*****************************************************************************
 * STREME > model > objfun
 ****************************************************************************/
  void (*handle_objfun) (void *ctx, STREME_OBJFUN_EN objfun);

/*****************************************************************************
 * STREME > model > test 
 ****************************************************************************/
  void (*handle_test) (void *ctx, STREME_TEST_EN test);

/*****************************************************************************
 * STREME > model > minw
 ****************************************************************************/
  void (*handle_minw) (void *ctx, int minw);

/*****************************************************************************
 * STREME > model > maxw
 ****************************************************************************/
  void (*handle_maxw) (void *ctx, int maxw);

/*****************************************************************************
 * STREME > model > kmer 
 ****************************************************************************/
  void (*handle_kmer) (void *ctx, int kmer);

/*****************************************************************************
 * STREME > model > hofract
 ****************************************************************************/
  void (*handle_hofract) (void *ctx, double hofract);

/*****************************************************************************
 * STREME > model > neval 
 ****************************************************************************/
  void (*handle_neval) (void *ctx, int neval);

/*****************************************************************************
 * STREME > model > nref 
 ****************************************************************************/
  void (*handle_nref) (void *ctx, int nref);

/*****************************************************************************
 * STREME > model > niter 
 ****************************************************************************/
  void (*handle_niter) (void *ctx, int niter);

/*****************************************************************************
 * STREME > model > patience
 ****************************************************************************/
  void (*handle_patience) (void *ctx, long patience);

/*****************************************************************************
 * STREME > model > seed 
 ****************************************************************************/
  void (*handle_seed) (void *ctx, long seed);

/*****************************************************************************
 * STREME > model > useer
 ****************************************************************************/
  void (*handle_useer) (void *ctx, int useer);

/*****************************************************************************
 * STREME > model > minscore
 ****************************************************************************/
  void (*handle_minscore) (void *ctx, long minscore);

/*****************************************************************************
 * STREME > model > ignore_depth
 ****************************************************************************/
  void (*handle_ignore_depth) (void *ctx, long ignore_depth);

/*****************************************************************************
 * STREME > model > nsubsets 
 ****************************************************************************/
  void (*handle_nsubsets) (void *ctx, long nsubsets);

/*****************************************************************************
 * STREME > model > min_pal_ratio
 ****************************************************************************/
  void (*handle_min_pal_ratio) (void *ctx, double min_pal_ratio);

/*****************************************************************************
 * STREME > model > max_pal_ed
 ****************************************************************************/
  void (*handle_max_pal_ed) (void *ctx, double max_pal_ed);

/*****************************************************************************
 * STREME > model > cand 
 ****************************************************************************/
  void (*handle_cand) (void *ctx, int cand);

/*****************************************************************************
 * STREME > model > experimental
 ****************************************************************************/
  void (*handle_experimental) (void *ctx, int experimental);

/*****************************************************************************
 * STREME > model > totallength
 ****************************************************************************/
  void (*handle_totallength) (void *ctx, long totallength);

/*****************************************************************************
 * STREME > model > align
 ****************************************************************************/
  void (*handle_align) (void *ctx, long align);

/*****************************************************************************
 * STREME > model > host 
 ****************************************************************************/
  void (*handle_host) (void *ctx, char *host);

/*****************************************************************************
 * STREME > model > description
 ****************************************************************************/
  void (*handle_description)(void *ctx, char *experiment_description);

/*****************************************************************************
 * STREME > motifs
 * contains: motif*
 ****************************************************************************/
  void (*start_motifs)(void *ctx);

/*****************************************************************************
 * STREME > /motifs
 * contains: motif*
 ****************************************************************************/
  void (*end_motifs)(void *ctx);

/*****************************************************************************
 * STREME > motifs > motif
 *
 *  id                the identifier used by STREME
 *  alt               alternate identifier used by STREME
 *  width             the width of the motif
 *  initial_width     the width of the prefix seed for the motif
 *  seed              the seed word for the motif
 *  score_threshold   the log-odds score threshold for calling sites
 *  train_pos_count   positive training sites
 *  train_neg_count   negative training sites
 *  train_log_pvalue  log10 p-value on training set
 *  train_pvalue      p-value on training set
 *  train_dtc         CD only: average site distance to center in training set
 *  train_bernoulli   Binomial Test only: Bernoulli prob. in training set
 *  test_pos_count    positive test sites
 *  test_neg_count    negative test sites
 *  test_log_pvalue   log10 p-value on test set
 *  test_pvalue       p-value on test set
 *  test_dtc          CD only: average site distance to center in test set
 *  test_bernoulli    Binomial Test only: Bernoulli prob. in test set
 *  is_palindromic    is the motif a perfect palindrome?
 *  elapsed time      time elapsed until motif found (seconds)
 ****************************************************************************/
  void (*start_motif)(void *ctx, 
    char *id, char *alt, int width, int initial_width, char *seed, double score_threshold,
    long train_pos_count, long train_neg_count, double train_log_pvalue, char *train_pvalue, double train_dtc, double train_bernoulli,
    long test_pos_count, long test_neg_count, double test_log_pvalue, char *test_pvalue, double test_dtc, double test_bernoulli,
    char *is_palindromic, double elapsed_time);

/*****************************************************************************
 * STREME > motifs > /motif
 * contains: pos+, match*
 ****************************************************************************/
  void (*end_motif)(void *ctx);

/*****************************************************************************
 * STREME > motifs > motif > pos
 *
 *  position          index of the motif position
 *  nfreqs            number of frequencies
 *  freqs             array of frequencies
 ****************************************************************************/
  void (*handle_pos)(void *ctx, int position, int nfreqs, double* freqs);

/*****************************************************************************
 * STREME > reason_for_stopping
 ****************************************************************************/
  void (*handle_reason_for_stopping)(void *ctx, char *reason_for_stopping);

/*****************************************************************************
 * STREME > run_time
 *
 *  cpu             the time the program was running on a cpu.
 ****************************************************************************/
  void (*handle_run_time)(void *ctx, double cpu_time);

};
typedef struct streme_io_xml_callbacks STREME_IO_XML_CALLBACKS_T;

/*****************************************************************************
 * Register handlers on the xmlSAXHandler structure
 ****************************************************************************/
void register_streme_io_xml_sax_handlers(xmlSAXHandler *handler);

/*****************************************************************************
 * Creates the data to be passed to the SAX parser
 ****************************************************************************/
void* create_streme_io_xml_sax_context(void *user_data, STREME_IO_XML_CALLBACKS_T *callbacks);

/*****************************************************************************
 * Destroys the data used by the SAX parser
 ****************************************************************************/
void destroy_streme_io_xml_sax_context(void *context);

#endif
