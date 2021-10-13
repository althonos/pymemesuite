/********************************************************************
 * FILE: ame.c
 * AUTHOR: Robert McLeay & Timothy L. Bailey
 * CREATE DATE: 19/08/2008
 * PROJECT: MEME suite
 * COPYRIGHT: 2008-2017, Robert McLeay & Timothy L. Bailey
 *
 * AME seeks to assist in
 * determining whether a given transcription factor regulates a set
 * of genes that have been found by an experimental method that
 * provides ranks, such as microarray or ChIP-chip.
 ********************************************************************/
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>
#include <libgen.h>
#include "alphabet.h"
#include "ame.h"
#include "array.h"
#include "config.h"
#include "fasta-io.h"
#include "fisher_exact.h"
#include "io.h"
#include "matrix.h"
#include "motif-db.h"
#include "motif-in.h"
#include "motif_regexp.h"
#include "projrel.h"
#include "ranksum_test.h"
#include "regress.h"
#include "simple-getopt.h"
#include "spearman-rank-correlation.h"
#include "string-list.h"
#include "utils.h"
#include "red-black-tree.h"

VERBOSE_T verbosity = NORMAL_VERBOSE;

#define MIN_SEQS 2		// Require at least this many sequences in primary and control sets to prevent crashing.  Change min_seqs to match in website/js/ame.js if this changes.
#define MIN_PEARSON_SAMPLES 30	// Don't trust p-value estimate below 30 samples.
#define NO_WARNING 0
#define DUPLICATE_WARNING 1
#define MARGIN_WARNING 2

// Globals.
bool first_sig_motif = true;	// first significant motif

// Options struct.
typedef struct  {
  char* alph_file; // file containing xalph
  ALPH_T *alphabet; // the alphabet
  char* bg_source; //file to get base freqs from
  char* outputdir; // where to send outputs
  ARRAYLST_T *motif_sources; //filenames of the motif library
  bool scan_both_strands; // scan both strands
  bool scan_separately; // scan strands separately
  char* seq_source; //(primary) input sequences in fasta format
  char* control_filename; //control input sequences in fasta format
  int kmer;		// preserve counts of words of this size when shuffling
  int seed;		// random number seed
  char* commandline; // command line with path stripped from program name
  float pseudocount; //add to the motif frequency counts
  int scoring; //AVG_ODDS or MAX_ODDS
  int verbose;
  int rs_method; //QUICK_RS or BETTER_RS
  int positive_list; //POS_FASTA or POS_PWM
  int pvalue_method; //RANKSUM_METHOD or FISHER_METHOD
  double fisher_fasta_threshold;
  //double pvalue_threshold; //threshold for the mhg test with pwms.
  double hit_lo_fraction;		// fraction of maximum log-odds for strong site
  double evalue_report_threshold; //threshold for reporting a motif in output.
  bool log_fscores;
  bool log_pwmscores;
  bool linreg_switchxy;
  bool clobber; // true if we can replace existing output directory
  bool text_only; // true if we only are output the TSV results to stdout
  bool noseq;	// suppress output of sequence TSV file
  char* pearson_dump_dir;
  int fix_partition;
  ARRAYLST_T *include_patterns; // Wildcard patterns of motif names to include
  ARRAYLST_T *exclude_patterns; // Wildcard patterns of motif names to exclude
  // derived from command line:
  FILE *tsv_output;
  FILE *seq_output;
  HTMLWR_T *html_output;
  JSONWR_T *json_output;
} AME_OPTION_T;

typedef struct {
  int db_idx;
  MOTIF_T *motif;
  int pos;
  int neg;
  int tp;
  int fp;
  double p_thresh;
  double tp_thresh;
  double log_pleft;
  double log_pright;
  double log_pboth;
  double u;
  double rho;
  int num_tests;
  double log_corrected_pvalue;
  double log_evalue;
  bool is_significant;
} AME_RESULT_T;

typedef struct {
  int seq_idx;
  double f_rank;
  double f_score;
  double pwm_rank;
  double pwm_score;
  double rand;
} AME_RANK_T;

/*
 * Global Variables
 */
static AME_OPTION_T default_options;
static AME_OPTION_T options;
static time_t  t0; /* measuring time */
static clock_t c0; /* measuring cpu_time */
static double* factorial_array;
#define DEFAULTDIR "ame_out"
static const char *TSV_FILENAME = "ame.tsv";
static const char *SEQ_FILENAME = "sequences.tsv";
static const char *HTML_FILENAME = "ame.html";
static const char *TEMPLATE_FILENAME = "ame_template.html";
const char *PVALUE_METHOD_NAMES[] = {"Wilcoxon rank-sum test", "Fisher's exact test", "multihg", 
  "long_multihg", "Pearson's correlation coefficient", "Spearman's correlation coefficient"};
const char *SCORING_METHOD_NAMES[] = {"avg_odds", "max_odds", "sum_odds", 
  "total_hits"};

/*************************************************************************
 * Initializations of parameters to defaults
 *************************************************************************/

void ame_set_defaults() {
  default_options.alphabet = NULL;
  default_options.alph_file = NULL;
  default_options.bg_source = NULL;
  default_options.outputdir = DEFAULTDIR;
  default_options.motif_sources = arraylst_create();
  default_options.scan_both_strands = true; // AME scans both strands.
  default_options.scan_separately = false;
  default_options.seq_source = NULL;
  default_options.control_filename = NULL;
  default_options.kmer = 2;
  default_options.seed = 1;
  default_options.commandline = NULL;
  default_options.pseudocount = 0.1;
  default_options.scoring = AVG_ODDS; 
  default_options.verbose = NORMAL_VERBOSE;
  default_options.rs_method = QUICK_RS;
  default_options.positive_list = POS_FASTA;
  default_options.pvalue_method = FISHER_METHOD; 
  default_options.fisher_fasta_threshold = 1e-3;
  //default_options.pvalue_threshold = 2e-4;
  default_options.hit_lo_fraction = 0.25;
  default_options.evalue_report_threshold = 10.0;
  default_options.log_pwmscores = false;
  default_options.log_fscores = false;
  default_options.linreg_switchxy = true;
  default_options.clobber = true;
  default_options.text_only = false;
  default_options.noseq = false;
  default_options.fix_partition = -1; //i.e. disabled.
  default_options.include_patterns = arraylst_create();
  default_options.exclude_patterns = arraylst_create();
  default_options.tsv_output = stdout;
  default_options.seq_output = NULL;
  default_options.html_output = NULL;
  default_options.json_output = NULL;
  // Now copy the defaults into the real options
  memcpy(&options, &default_options, sizeof(AME_OPTION_T));
} // ame_set_defaults

/***********************************************************************
* Print the contents or a result instance for debug.
***********************************************************************/
void print_rsr(AME_RESULT_T *rsr){
  fprintf(stderr, "db_idx %d motif %s pos %d neg %d tp %d fp %d\n", 
    rsr->db_idx, get_motif_st_id(rsr->motif), rsr->pos, rsr->neg, rsr->tp, rsr->fp);
  fprintf(stderr, "p_thresh %g tp_thresh %g log_pleft %g log_pright %g log_pboth %g u %g rho %g\n",
    rsr->p_thresh, rsr->tp_thresh, rsr->log_pleft, rsr->log_pright, rsr->log_pboth, rsr->u, rsr->rho);
}

/***********************************************************************
* Create a new instance of a result, with a given motif and initial values
* if a non-null pointer is passed in, use its value, otherwise pass overwrite
* the location pointed to by new_result; in either case return new_result
***********************************************************************/
AME_RESULT_T *init_result(AME_RESULT_T *new_result,
  int db_idx,
  MOTIF_T *motif,
  int pos,
  int neg,
  int tp,
  int fp,
  double p_thresh,
  double tp_thresh,
  double log_pleft,
  double log_pright,
  double log_pboth,
  double u,
  double rho,
  int num_tests
) {
  if (!new_result) new_result = mm_malloc(sizeof(AME_RESULT_T));
  new_result->db_idx = db_idx;
  new_result->motif = motif;
  new_result->pos = pos;
  new_result->neg = neg;
  new_result->tp = tp;
  new_result->fp = fp;
  new_result->p_thresh = p_thresh;
  new_result->tp_thresh = tp_thresh;
  new_result->log_pleft = log_pleft;
  new_result->log_pright = log_pright;
  new_result->log_pboth = log_pboth;
  new_result->u = u;
  new_result->rho = rho;
  new_result->num_tests = num_tests;
  return new_result;
}

/***********************************************************************
 Free memory allocated in options processing
 ***********************************************************************/
static void cleanup_options() {
  alph_release(options.alphabet);
  arraylst_destroy(NULL, options.include_patterns);
  arraylst_destroy(NULL, options.exclude_patterns);
  arraylst_destroy(NULL, options.motif_sources);
}

/*************************************************************************
 * Setup the JSON writer and output a lot of pre-calculation data
 *************************************************************************/
static void ame_start_output(
  int argc, 
  char** argv, 
  ARRAY_T* bg_freqs, 
  SEQ_T** sequences, 
  int pos_num_seqs, 
  int num_skipped, 
  int neg_num_seqs, 
  int neg_num_skipped, 
  bool fasta_scores,
  ARRAYLST_T* dbs
) {
  int i;
  MOTIF_DB_T* db;

  // setup the html output writer 
  if ((options.html_output = htmlwr_create(get_meme_data_dir(), TEMPLATE_FILENAME, false))) {
    htmlwr_set_dest_name(options.html_output, options.outputdir, HTML_FILENAME);
    htmlwr_replace(options.html_output, "ame_data.js", "data");
    options.json_output = htmlwr_output(options.html_output);
    if (options.json_output == NULL) die("Template does not contain data section.\n");
  } else {
    DEBUG_MSG(QUIET_VERBOSE, "Failed to open HTML template file.\n");
    options.html_output = NULL;
    options.json_output = NULL;
  }

  // write out some information
  jsonwr_str_prop(options.json_output, "version", VERSION);
  jsonwr_str_prop(options.json_output, "revision", REVISION);
  jsonwr_str_prop(options.json_output, "release", ARCHIVE_DATE);
  jsonwr_str_prop(options.json_output, "program", "AME");
  //jsonwr_str_array_prop(options.json_output, "cmd", argv, argc);
  jsonwr_args_prop(options.json_output, "cmd", argc, argv);
  // options
  jsonwr_property(options.json_output, "options");
  jsonwr_start_object_value(options.json_output);
  jsonwr_str_prop(options.json_output, "scoring", SCORING_METHOD_NAMES[options.scoring]);
  //jsonwr_str_prop(options.json_output, "rs_method", (options.rs_method ? "better" : "quick"));
  jsonwr_str_prop(options.json_output, "positive_list", (options.positive_list == POS_FASTA) ? "FASTA" : "MOTIF");
  jsonwr_str_prop(options.json_output, "pvalue_method", PVALUE_METHOD_NAMES[options.pvalue_method]);
  jsonwr_dbl_prop(options.json_output, "pseudocount", options.pseudocount);
  jsonwr_dbl_prop(options.json_output, "fisher_fasta_threshold", options.fisher_fasta_threshold);
  jsonwr_dbl_prop(options.json_output, "hit_lo_fraction", options.hit_lo_fraction);
  jsonwr_dbl_prop(options.json_output, "evalue_report_threshold", options.evalue_report_threshold);
  jsonwr_bool_prop(options.json_output, "noseq", options.noseq); 
  jsonwr_bool_prop(options.json_output, "log_fscores", options.log_fscores);
  jsonwr_bool_prop(options.json_output, "log_pwmscores", options.log_pwmscores);
  jsonwr_bool_prop(options.json_output, "linreg_switchxy", options.linreg_switchxy);
  jsonwr_bool_prop(options.json_output, "fix_partition", (options.control_filename || options.fix_partition > 0));

  if (options.control_filename) {
    jsonwr_str_prop(options.json_output, "partition", "number of primary sequences");
  } else if (options.fix_partition > 0) {
    jsonwr_lng_prop(options.json_output, "partition", options.fix_partition);
  }
  jsonwr_dbl_prop(options.json_output, "kmer", options.kmer); 
  jsonwr_end_object_value(options.json_output);

  // output the motif dbs
  jsonwr_property(options.json_output, "motif_dbs");
  jsonwr_start_array_value(options.json_output);
  for (i = 0; i < arraylst_size(options.motif_sources); i++) {
    db = arraylst_get(i, dbs);
    jsonwr_start_object_value(options.json_output);
    jsonwr_str_prop(options.json_output, "source", db->source);
    jsonwr_lng_prop(options.json_output, "count", arraylst_size(db->motifs));
    jsonwr_end_object_value(options.json_output);
  }
  // finish writing motif dbs
  jsonwr_end_array_value(options.json_output);

  // Write alphabet to json data.
  jsonwr_property(options.json_output, "alphabet");
  alph_print_json(options.alphabet, options.json_output);

  // Write the background model to json data.
  jsonwr_property(options.json_output, "background");
  jsonwr_start_object_value(options.json_output);
  jsonwr_str_prop(options.json_output, "source", options.bg_source);
  if (!strcmp(options.bg_source, "--sequences--")) jsonwr_str_prop(options.json_output, "file", options.seq_source);
  jsonwr_property(options.json_output, "frequencies");
  jsonwr_start_array_value(options.json_output);
  for (i = 0; i < alph_size_core(options.alphabet); i++) {
    jsonwr_dbl_value(options.json_output, get_array_item(i, bg_freqs));
  }
  jsonwr_end_array_value(options.json_output);
  jsonwr_end_object_value(options.json_output);

  // Output the sequence databases
  jsonwr_property(options.json_output, "sequence_db");
  jsonwr_start_object_value(options.json_output);
  jsonwr_str_prop(options.json_output, "source", options.seq_source);
  jsonwr_lng_prop(options.json_output, "count", pos_num_seqs);
  jsonwr_lng_prop(options.json_output, "skipped", num_skipped);
  jsonwr_bool_prop(options.json_output, "fasta_scores", fasta_scores);
  jsonwr_end_object_value(options.json_output);
  if (options.control_filename) {
    jsonwr_property(options.json_output, "control_db");
    jsonwr_start_object_value(options.json_output);
    jsonwr_str_prop(options.json_output, "source", options.control_filename);
    jsonwr_lng_prop(options.json_output, "count", neg_num_seqs);
    jsonwr_lng_prop(options.json_output, "skipped", neg_num_skipped);
    jsonwr_end_object_value(options.json_output);
  }
} // ame_start_output

/*****************************************************************************
* Get the FASTA scores.
* Check if all sequences have FASTA scores in their ID lines.
* Set FASTA scores equal to input order if not all sequences have them.
*****************************************************************************/
ARRAY_T *read_fasta_scores(
  int num_seqs,			// the number of sequences
  SEQ_T **sequences, 		// the sequences
  bool *fasta_scores		// OUT set to true if FASTA scores were given
) {
  int i;

  bool scores_found = false;
  int n_found = 0;
  ARRAY_T *seq_fscores = allocate_array(num_seqs);

  // Read in the FASTA scores.
  for (i=0; i<num_seqs; i++) {
    char *endptr;
    double fscore = strtod(get_seq_description(sequences[i]), &endptr);
    if (endptr != get_seq_description(sequences[i])) n_found++;
    set_array_item(i, fscore, seq_fscores);
  }

  // Check.
  if (n_found == 0) {
    DEBUG_MSG(NORMAL_VERBOSE, "No FASTA scores were read from sequence ID lines.\n");
  } else if (n_found < num_seqs) {
    DEBUG_FMT(NORMAL_VERBOSE, "Only %d of %d sequences contained FASTA scores in their ID lines.\n", n_found, num_seqs);
  } else {
    DEBUG_MSG(NORMAL_VERBOSE, "FASTA scores were read from the sequence ID lines.\n");
    scores_found = (n_found == num_seqs);
  }

  // See if scores were required.
  if (!scores_found) {
    if (PEARSON_METHOD == options.pvalue_method) {
      die("You must specify scores in the FASTA IDs of your sequences with --method pearson.\n");
    }
    DEBUG_MSG(NORMAL_VERBOSE, "Setting FASTA scores equal to the input order of the sequence.\n");
    for (i=0; i<num_seqs; i++) set_array_item(i, i+1, seq_fscores);
  }

  *fasta_scores = scores_found;
  return(seq_fscores);
} // read_fasta_scores

/*************************************************************************
 * Read all the sequences into an array of SEQ_T *.
 *************************************************************************/
static void ame_read_sequences(
  int max_width,		// ignore sequences shorter than this
  ALPH_T *alph,			// sequence alphabet
  char *seq_source,		// sequence file name
  char *type,			// for error messages
  int *num_seqs,		// IN/OUT the total number of sequences (all files so far)
  int *num_skipped,		// IN/OUT the total number of sequences skipped (this file)
  double *avg_length,		// IN/OUT average length of sequences
  SEQ_T ***sequences,		// IN/OUT sequences read (all files so far)
  bool get_fasta_scores,	// get seq_fscores and set fasta_scores
  bool *fasta_scores,		// OUT set true/false depending on if FASTA scores were provided in sequence IDs 
  ARRAY_T **seq_fscores		// OUT the FASTA scores (or input order ranks) of the (primary) sequences
) {

  // Open the sequence file.
  FILE *seq_file = NULL;
  if (open_file(seq_source, "r", false, "FASTA", "sequences", &seq_file) == 0) {
    die("Couldn't open the %ssequence file '%s'.\n", type, seq_source);
  }

  // Read in the sequences.
  int old_num_seqs = *num_seqs;
  read_many_fastas(alph, seq_file, MAX_SEQ_LENGTH, num_seqs, sequences);

  // Remove the new sequences that are too short.
  int i, j;
  double total_length = 0;
  for (i=j=old_num_seqs; i < *num_seqs; i++) {
    int length = get_seq_length((*sequences)[i]);
    if (length < max_width) {
      char *name = get_seq_name((*sequences)[i]);
      DEBUG_FMT(NORMAL_VERBOSE, "Removing %ssequence '%s' from file '%s' because it is shorter (%d) than widest motif (%d).\n", 
	type, name, seq_source, length, max_width);
      free_seq((*sequences)[i]);
    } else {
      total_length += length;
      (*sequences)[j++] = (*sequences)[i];
    }
  }

  if (j == old_num_seqs) die("All %ssequences are shorter than the longest motif (%d).\n", type, max_width);

  if (j != *num_seqs) {
    DEBUG_FMT(NORMAL_VERBOSE, "Warning: Removed %d %ssequence(s) that were shorter than the longest motif (%d).\n", 
      *num_seqs-j, type, max_width);
  }
  *num_skipped = *num_seqs - j;	// number of sequences skipped (from this file)
  *num_seqs = j;		// number of sequences saved (from all files)
  *avg_length = ((old_num_seqs*(*avg_length))+total_length) / *num_seqs;
  if (get_fasta_scores) *seq_fscores = read_fasta_scores(*num_seqs, *sequences, fasta_scores);

} // ame_read_sequences

/*************************************************************************
 * Send line to tsv motif output.
 *************************************************************************/
void tsvwr (char* format, ...) {
  va_list arg_ptrs;
  va_start(arg_ptrs, format);
  vfprintf(options.tsv_output, format, arg_ptrs);
  va_end(arg_ptrs);
} // tsvwr

/*************************************************************************
 * Send line to sequence output.
 *************************************************************************/
void seqwr (char* format, ...) {
  va_list arg_ptrs;
  va_start(arg_ptrs, format);
  vfprintf(options.seq_output, format, arg_ptrs);
  va_end(arg_ptrs);
} // seqwr

/*****************************************************************************
* Output the TP and FP sequences in TSV format.
*****************************************************************************/
void output_tp_fp_sequences(
  bool first_motif,		// true if first motif ever
  char *db_source,		// filename of the motif DB
  MOTIF_T *motif,		// the motif
  AME_RESULT_T *result,		// the motif's result
  int num_seqs,			// number of sequences in rankings
  SEQ_T **sequences,		// the sequences
  AME_RANK_T **rankings		// sequence rankings
) {
  int i;

  bool pos_fasta = (POS_FASTA == options.positive_list);

  // Open the file and print the header line.
  if (first_motif) {
    const char *MOTIF_DB = "motif_DB";
    const char *MOTIF_ID = "motif_ID";
    const char *MOTIF_ALT_ID = "motif_ALT_ID";
    const char *SEQ_ID = "seq_ID";
    const char *FASTA_SCORE = "FASTA_score"; 
    const char *PWM_SCORE = "PWM_score";
    const char *CLASS = "class";

    // Create the list of headers for the given method.
    ARRAYLST_T *headers = arraylst_create();
    arraylst_add((char *)MOTIF_DB, headers);
    arraylst_add((char *)MOTIF_ID, headers);
    arraylst_add((char *)MOTIF_ALT_ID, headers);
    arraylst_add((char *)SEQ_ID, headers);
    if (pos_fasta) {
      arraylst_add((char *)FASTA_SCORE, headers);
      arraylst_add((char *)PWM_SCORE, headers);
    } else {
      arraylst_add((char *)PWM_SCORE, headers);
      arraylst_add((char *)FASTA_SCORE, headers);
    }
    arraylst_add((char *)CLASS, headers);

    // Print the header line.
    for (i=0; i<arraylst_size(headers); i++) {
      seqwr(i==0?"%s":"\t%s", arraylst_get(i, headers));
    }
    seqwr("\n");
    arraylst_destroy(NULL, headers);
  }

  // Print the true and false positives;
  int n_found = 0;
  for (i=0; i<num_seqs; i++) {
    bool label_pos;
    bool class_pos;
    double p_score, tp_score;
    if (pos_fasta) {
      p_score = rankings[i]->f_score;
      tp_score = rankings[i]->pwm_score;
      label_pos = (p_score <= result->p_thresh);
      class_pos = (tp_score >= result->tp_thresh);
    } else {
      p_score = rankings[i]->pwm_score;
      tp_score = rankings[i]->f_score;
      label_pos = (p_score >= result->p_thresh);
      class_pos = (tp_score <= result->tp_thresh);
    }
    // Print only classified positives (tp and fp);
    if (class_pos) {
      n_found++;
      if (!first_motif && n_found==1) seqwr("#\n");
      seqwr("%s\t%s\t%s\t%s\t%.6g\t%.6g\t%s\n",
	db_source,
	get_motif_id(result->motif),
	get_motif_id2(result->motif) ? get_motif_id2(result->motif) : "",
	get_seq_name(sequences[rankings[i]->seq_idx]),
	p_score,
	tp_score,
	((label_pos && class_pos) ? "tp" : "fp")
      );
    }
  } // seqs
} // output_tp_fp_sequences

/*****************************************************************************
 * Write a motif result to JSON.
 *****************************************************************************/
static void output_json_result(
  AME_RESULT_T* result 
) {
  int i, j, len, alen;
  MATRIX_T *freqs;
  ARRAY_T *row;
  jsonwr_start_object_value(options.json_output);
  jsonwr_lng_prop(options.json_output, "db", result->db_idx);
  jsonwr_str_prop(options.json_output, "id", get_motif_id(result->motif));
  if (get_motif_id2(result->motif)[0] != '\0')
    jsonwr_str_prop(options.json_output, "alt", get_motif_id2(result->motif));
  jsonwr_lng_prop(options.json_output, "len", get_motif_length(result->motif));
  jsonwr_dbl_prop(options.json_output, "motif_evalue", get_motif_evalue(result->motif));
  jsonwr_dbl_prop(options.json_output, "motif_nsites", get_motif_nsites(result->motif));
  if (has_motif_url(result->motif)) 
    jsonwr_str_prop(options.json_output, "url", get_motif_url(result->motif));
  alen = alph_size_core(get_motif_alph(result->motif));
  len = get_motif_length(result->motif);
  freqs = get_motif_freqs(result->motif);
  jsonwr_property(options.json_output, "pwm");
  jsonwr_start_array_value(options.json_output);
  for (i = 0; i < len; i++) {
    row = get_matrix_row(i, freqs);
    jsonwr_start_array_value(options.json_output);
    for (j = 0; j < alen; j++) {
      jsonwr_dbl_value(options.json_output, get_array_item(j, row));
    }
    jsonwr_end_array_value(options.json_output);
  }
  jsonwr_end_array_value(options.json_output);
    
  // Stuff for partitioning or for optimizing over PWM in Fisher method
  jsonwr_lng_prop(options.json_output, "pos", result->pos);
  jsonwr_lng_prop(options.json_output, "neg", result->neg);
  jsonwr_lng_prop(options.json_output, "tp", result->tp);
  jsonwr_lng_prop(options.json_output, "fp", result->fp);
  jsonwr_dbl_prop(options.json_output, "positive_thresh", result->p_thresh);
  jsonwr_dbl_prop(options.json_output, "true_positive_thresh", result->tp_thresh);

  // Stuff for all methods
  jsonwr_dbl_prop(options.json_output, "number_of_tests", result->num_tests);
  char *corrected_pvalue = print_log_value(NULL, result->log_corrected_pvalue, 2);
  char *evalue = print_log_value(NULL, result->log_evalue, 2);
  jsonwr_str_prop(options.json_output, "corrected_pvalue", corrected_pvalue);
  jsonwr_str_prop(options.json_output, "evalue", evalue);
  myfree(corrected_pvalue);
  myfree(evalue);

  // Extra stuff for PEARSON_METHOD
  if (options.pvalue_method == PEARSON_METHOD) {
    jsonwr_dbl_prop(options.json_output, "pearsons_rho", result->rho);
    jsonwr_dbl_prop(options.json_output, "m", result->log_pleft);
    jsonwr_dbl_prop(options.json_output, "b", result->log_pboth);
    jsonwr_dbl_prop(options.json_output, "mean_square_error", result->u);
  } 

  // Extra stuff for SPEARMANS_METHOD
  if (options.pvalue_method == SPEARMAN_METHOD) {
    jsonwr_dbl_prop(options.json_output, "spearmans_rho", result->rho);
  }

  jsonwr_end_object_value(options.json_output);
} // output_json_result

/*************************************************************************
 * Close the output files.
 *************************************************************************/
static void final_print_results () {
  fclose(options.tsv_output);
  if (options.html_output) {
    if (htmlwr_output(options.html_output) != NULL) {
      die("Expected only one replacement variable in template.\n");
    }
    htmlwr_destroy(options.html_output);
    options.html_output = NULL;
  }
} // final_print_results

/*************************************************************************
 * End AME.
 *************************************************************************/
void ame_terminate(
    AME_OPTION_T *options,
    int status
  ) {
  /* Display time of execution */
  if (verbosity >= HIGH_VERBOSE) {
    time_t t1 = time(NULL);
    clock_t c1 = clock();
    fprintf (stderr, "Elapsed wall clock time: %ld seconds\n", (long)(t1 - t0));
    fprintf (stderr, "Elapsed CPU time:        %g seconds\n", 
      (float)(c1 - c0) / CLOCKS_PER_SEC);
  }
  // don't risk closing the files if something went wrong: 
  // should be more specific on when not to do this
  if (!status && !options->text_only && options->seq_source) final_print_results();
  exit(status);
} // ame_terminate

/*
 * COMPARISON METHODS FOLLOW
 */

/*************************************************************************
 * Compare two log evalues, increasing.
 * Break ties using motif idx, then db idx.
 *************************************************************************/
int ame_compare_log_evalues(const void *a, const void *b)
{
  AME_RESULT_T r1 = **(AME_RESULT_T**)a;
  AME_RESULT_T r2 = **(AME_RESULT_T**)b;
  if (r2.log_evalue - r1.log_evalue < 0.0) {
    return 1;
  } else if (r2.log_evalue - r1.log_evalue > 0.0){
    return -1;
  } else {
    if (r1.motif->idx < r2.motif->idx) {
      return 1;
    } else if (r1.motif->idx > r2.motif->idx) {
      return -1;
    } else {
      return(r1.db_idx - r2.db_idx);
    }
  }
} // ame_compare_log_evalues

/*************************************************************************
 * Compare two mean squared error, increasing.
 *************************************************************************/
int ame_compare_mse (const void *a, const void *b)
{
  AME_RESULT_T r1 = **(AME_RESULT_T**)a;
  AME_RESULT_T r2 = **(AME_RESULT_T**)b;
  if (r2.log_pboth - r1.log_pboth < 0.0) {
    return 1;
  } else if (r2.log_pboth - r1.log_pboth > 0.0){
    return -1;
  } else {
    return 0;
  }
} // ame_compare_mse

/*************************************************************************
 * Compare two doubles, increasing.
 *************************************************************************/
int ame_compare_doubles (const void *a, const void *b)
{
  if (*(const double *)a - *(const double *)b < 0) {
    return 1;
  } else if (*(const double *)a - *(const double *)b > 0) {
    return -1;
  } else {
    return 0;
  }
} // ame_compare_doubles

/*************************************************************************
 * Compare two f_scores, increasing.  Break ties conserving sequence order.
 *************************************************************************/
int ame_compare_ranks_fasta_score (const void *a, const void *b) {
  AME_RANK_T r1 = **(AME_RANK_T**)a;
  AME_RANK_T r2 = **(AME_RANK_T**)b;

  if (r2.f_score - r1.f_score < 0.0) {
    return 1;
  } else if (r2.f_score - r1.f_score > 0.0){
    return -1;
  } else {
    // keeps sequence order if all scores are the same
    // so that input order of sequences matters
    return(r2.seq_idx - r1.seq_idx);
  }
} // ame_compare_ranks_fasta_score

/*************************************************************************
 * Compare two pwm_scores, decreasing.  Break ties randomly.
 *************************************************************************/
int ame_compare_ranks_pwm_score_random (const void *a, const void *b) {
  AME_RANK_T r1 = **(AME_RANK_T**)a;
  AME_RANK_T r2 = **(AME_RANK_T**)b;

  if (r1.pwm_score - r2.pwm_score < 0.0) {
    return 1;
  } else if (r1.pwm_score - r2.pwm_score > 0.0) {
    return -1;
  } else if (r2.rand - r1.rand < 0.0) {
    return 1;
  } else if (r2.rand - r1.rand > 0.0) {
    return -1;
  } else {
    // Just in case 2 random numbers are the same!
    return(r2.seq_idx - r1.seq_idx);
  }
} // ame_compare_ranks_pwm_score_random

/*************************************************************************
 * Compare two pwm_scores, decreasing.  Break ties on FASTA ranks, decreasing
 * so that runs of equal PWM have negatives first.
 *************************************************************************/
int ame_compare_ranks_pwm_score_fasta (const void *a, const void *b) {
  AME_RANK_T r1 = **(AME_RANK_T**)a;
  AME_RANK_T r2 = **(AME_RANK_T**)b;

  if (r1.pwm_score - r2.pwm_score < 0.0) {
    return 1;
  } else if (r1.pwm_score - r2.pwm_score > 0.0){
    return -1;
  } else {
    // breaks ties putting *larger* FASTA ranks first
    // so that runs of equal PWM have negatives first
    return(r2.f_rank - r1.f_rank);
  }
} // ame_compare_ranks_pwm_score_fasta

/*****************************************************************************
 * Print a usage message with an optional reason for failure.
 *****************************************************************************/
void ame_usage(char *format, ...) {
  va_list argp;
  char *usage = 
    "Usage: ame [options] <sequence file> <motif file>+\n"
    "     <sequence file>           file of sequences in FASTA format\n"
    "     <motif file>+             file(s) of motifs in MEME format\n"
    "\n"
    "     --o <output dir>          output directory; default: ame_out\n"
    "     --oc <output dir>         overwrite output; default: ame_out\n"
    "     --text                    output TSV format to stdout; overrides --o and --oc;\n"
    "                               default: create HTML and TSV files in <output_dir>\n"
    "     --control <control file>  control sequences in FASTA format or the keyword\n"
    "                               '--shuffle--' to use shuffled versions of the\n"
    "                               primary sequences\n"
    "     --kmer <k>                preserve k-mer frequencies when shuffling;\n"
    "                               default: %d\n"
    "     --seed <s>                random number seed (integer); default: %d\n"
    "     --method [fisher|3dmhg|4dmhg|ranksum|pearson|spearman]\n"
    "                               statistical test; default: fisher\n"
    "     --scoring [avg|max|sum|totalhits]\n"
    "                               sequence scoring method; default: avg\n"
    "     --hit-lo-fraction <fraction>\n"
    "                               fraction of maximum log-odds for a hit;\n"
    "                               default: %g\n"
    "     --evalue-report-threshold <ev>\n"
    "                               motif significance reporting threshold;\n"
    "                               default: %g\n"
    "     --fasta-threshold <ft>    maximum FASTA score for sequence to be positive\n"
    "                               (requires --poslist pwm); default: %g\n"
    "     --fix-partition <int>     number of sequences in positive partition;\n"
    "     --poslist [fasta|pwm]     partition on affinity (fasta) or motif (pwm)\n"
    "                               scores; default: fasta\n"
    "     --log-fscores             use log of FASTA scores (pearson) or log of\n"
    "                               ranks (spearman)\n"
    "     --log-pwmscores           use log of log of PWM scores\n"
    "     --linreg-switchxy         switch roles of X=FASTA scores and Y=PWM scores\n"
    "     --xalph <alph file>       motifs will be converted to this custom alphabet\n"
    "     --bfile <bfile>           background model file; default: motif file freqs\n"
    "                               default: unconstrained partition maximization\n"
    "     --motif-pseudo <pc>       pseudocount for creating PWMs from motifs;\n"
    "                               default: %g\n"
    "     --inc <pattern>           name pattern to select as motif; may be\n"
    "                               repeated; default: all motifs are used\n"
    "     --exc <pattern>           name pattern to exclude as motif; may be\n"
    "                               repeated; default: all motifs are used\n"
    "     --verbose [1|2|3|4|5]     controls program verbosity (5=maximum verbosity);\n"
    "                               default: %d\n"
    "     --help                    print this message and exit\n"
    "     --version                 print the version and exit\n"
    "\n";

  if (format) {
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, usage, 
    default_options.hit_lo_fraction, 
    default_options.evalue_report_threshold, 
    default_options.fisher_fasta_threshold,
    default_options.pseudocount,
    default_options.kmer,
    default_options.seed,
    default_options.verbose
  );

  fflush(stderr);
  exit(EXIT_FAILURE);
} // ame_usage

/*****************************************************************************
* Get and process the command line options.
*****************************************************************************/
void ame_getopt(int argc, char *argv[]) {
  const int num_options = 27; // change this if the number of options changes
  cmdoption const motif_scan_options[] = {
    {"xalph", REQUIRED_VALUE},
    {"control", REQUIRED_VALUE},
    {"kmer", REQUIRED_VALUE},
    {"seed", REQUIRED_VALUE},
    {"bfile", REQUIRED_VALUE},
    {"bgfile", REQUIRED_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"text", NO_VALUE},
    {"motif-pseudo", REQUIRED_VALUE},
    {"motif", REQUIRED_VALUE},
    {"inc", REQUIRED_VALUE},
    {"exc", REQUIRED_VALUE},
    {"pseudocount", REQUIRED_VALUE},
    {"fasta-threshold", REQUIRED_VALUE},
    {"hit-lo-fraction", REQUIRED_VALUE},
    {"evalue-report-threshold", REQUIRED_VALUE},
    {"scoring", REQUIRED_VALUE},
    {"method", REQUIRED_VALUE},
    {"poslist", REQUIRED_VALUE},
    {"verbose", REQUIRED_VALUE},
    {"noseq", NO_VALUE},
    {"linreg-switchxy", NO_VALUE},
    {"log-fscores", NO_VALUE},
    {"log-pwmscores", NO_VALUE},
    {"fix-partition", REQUIRED_VALUE},
    {"h", NO_VALUE},
    {"help", NO_VALUE},
    {"version", NO_VALUE}
  };

  int option_index = 0;
  char* option_name = NULL;
  char* option_value = NULL;
  const char * message = NULL;
  bool bad_options = false;
  int i;

  if (simple_setopt(argc, argv, num_options, motif_scan_options) != NO_ERROR) {
    die("Error processing command line options: option name too long.\n");
  }

  /*
   * Now parse the command line options
   */
  //simple_getopt will return 0 when there are no more options to parse
  while(simple_getopt(&option_name, &option_value, &option_index) > 0) {
    if (options.verbose >= HIGHER_VERBOSE) {
      fprintf(stderr, "AME Option: %s - %s\n", option_name, option_value);
    }
    if (strcmp(option_name, "xalph") == 0) {
      options.alph_file = option_value;
    } else if (strcmp(option_name, "control") == 0) {
      options.control_filename = option_value;
    } else if (strcmp(option_name, "kmer") == 0) {
      options.kmer = atoi(option_value);
    } else if (strcmp(option_name, "seed") == 0) {
      options.seed = atoi(option_value);
    } else if (strcmp(option_name, "bfile") == 0 || strcmp(option_name, "bgfile") == 0) {
      options.bg_source = option_value;
    } else if (strcmp(option_name, "o") == 0) {
      options.outputdir = option_value;
      options.clobber = false;
    } else if (strcmp(option_name, "oc") == 0) {
      options.outputdir = option_value;
      options.clobber = true;
    } else if (strcmp(option_name, "text") == 0) {
      options.text_only = true;
    } else if (strcmp(option_name, "motif-pseudo") == 0
      || strcmp(option_name, "pseudocount") == 0) {
      options.pseudocount = atof(option_value);
    } else if (strcmp(option_name, "inc") == 0 || strcmp(option_name, "motif") == 0) {
     arraylst_add(option_value, options.include_patterns);
    } else if (strcmp(option_name, "exc") == 0) {
     arraylst_add(option_value, options.exclude_patterns);
    } else if (strcmp(option_name, "fasta-threshold") == 0) {
      options.fisher_fasta_threshold = atof(option_value);
    } else if (strcmp(option_name, "hit-lo-fraction") == 0) {
      options.hit_lo_fraction= atof(option_value);
    } else if (strcmp(option_name, "evalue-report-threshold") == 0) {
      options.evalue_report_threshold = atof(option_value);
    } else if (strcmp(option_name, "scoring") == 0) {
      if (strcmp(option_value,"avg")==0) {
        options.scoring = AVG_ODDS;
        DEBUG_MSG(NORMAL_VERBOSE, "Using average odds scoring.\n");
      } else if (strcmp(option_value,"sum") == 0) {
        options.scoring = SUM_ODDS;
        DEBUG_MSG(NORMAL_VERBOSE, "Using sum of odds scoring.\n");
      } else if (strcmp(option_value,"max") == 0) {
        options.scoring = MAX_ODDS;
        DEBUG_MSG(NORMAL_VERBOSE, "Using maximum odds scoring.\n");
      } else if (strcmp(option_value,"totalhits") == 0) {
        options.scoring = TOTAL_HITS;
        DEBUG_MSG(NORMAL_VERBOSE, "Using total hits scoring.\n");
      } else {
        ame_usage(NULL);
        ame_terminate(&options, 1);
      }
    } else if (strcmp(option_name, "poslist") == 0) {
      if (strcmp(option_value,"fasta") == 0) {
        options.positive_list = POS_FASTA;
      } else if (strcmp(option_value,"pwm") == 0) {
        options.positive_list = POS_PWM;
      } else {
        ame_usage(NULL);
        ame_terminate(&options, 1);
      }
    } else if (strcmp(option_name, "method") == 0) {
      if (strcmp(option_value,"ranksum")==0) {
        options.pvalue_method = RANKSUM_METHOD;
      } else if (strcmp(option_value,"fisher") == 0) {
        options.pvalue_method = FISHER_METHOD;
      } else if (strcmp(option_value,"3dmhg") == 0) {
        options.pvalue_method = MULTIHG_METHOD;
      } else if (strcmp(option_value,"4dmhg") == 0) {
        options.pvalue_method = LONG_MULTIHG_METHOD;
      } else if (strcmp(option_value,"pearson") == 0 || strcmp(option_value,"linreg") == 0) {
        options.pvalue_method = PEARSON_METHOD;
      } else if (strcmp(option_value,"spearman") == 0) {
        options.pvalue_method = SPEARMAN_METHOD;
      } else {
        ame_usage("Unrecognized --method value: %s\n", option_value);
        ame_terminate(&options, 1);
      }
    } else if (strcmp(option_name, "verbose") == 0) {
      options.verbose = atoi(option_value);
      verbosity = options.verbose;
      if (options.verbose <= INVALID_VERBOSE || options.verbose > DUMP_VERBOSE) {
        ame_usage("You must specify a value for --verbose in the range [1..5].\n");
        ame_terminate(&options, 1);
      }
    } else if (strcmp(option_name, "noseq") == 0) {
      options.noseq = true; 		// suppress sequence TSV file
    } else if (strcmp(option_name, "log-fscores") == 0) {
      options.log_fscores = true;
    } else if (strcmp(option_name, "log-pwmscores") == 0) {
      options.log_pwmscores = true;
    } else if (strcmp(option_name, "fix-partition") == 0) {
      options.fix_partition = atoi(option_value);
    } else if (strcmp(option_name, "linreg-switchxy") == 0) {
      options.linreg_switchxy = false;
    } else if (strcmp(option_name, "help") == 0 || strcmp(option_name, "h") == 0) {
      ame_usage(NULL);
      ame_terminate(&options, 0);
    } else if (strcmp(option_name, "version") == 0) {
      fprintf(stdout, VERSION "\n");
      exit(EXIT_SUCCESS);
    } else {
      ame_usage("Error: %s is not a recognised option.\n", option_name);
      ame_terminate(&options, 1);
    }

    option_index++;
  }

  if (options.pvalue_method == PEARSON_METHOD || options.pvalue_method == SPEARMAN_METHOD) {
    if (options.linreg_switchxy) {
      DEBUG_MSG(NORMAL_VERBOSE, "In LR/Spearman mode, x=PWM score, y=FASTA score\n");
    } else {
      //Standard mode
      DEBUG_MSG(NORMAL_VERBOSE, "In LR/Spearman mode, x=FASTA score, y=PWM score\n");
    }
  }

  // Must have sequence and motif files
  if (argc < option_index + 2) {
    ame_usage("Error: Must specify a sequence file and motif file(s).\n");
  }
  options.seq_source = argv[option_index++];

  // Record the motif file names
  for (; option_index < argc; option_index++) {
    arraylst_add(argv[option_index], options.motif_sources);
    DEBUG_FMT(NORMAL_VERBOSE, "Added %s to motif_sources which now has %d file names.\n", argv[option_index], arraylst_size(options.motif_sources));
    DEBUG_FMT(NORMAL_VERBOSE, "Motif file name is %s.\n", (char *) arraylst_get(0, options.motif_sources));
  }

  // --text implies --noseq
  if (options.text_only) options.noseq = true;

  /* Now validate the options.
   *
   * Illegal combinations are:
   *   - (control or spearman) and fix-partition
   *   - control and (spearman or pearson)
   *   - sequence bg and no sequence
   *   - no motif
   *   - no sequences
   *   - each file exists.
   */

  if (options.control_filename) {
    if (options.pvalue_method == PEARSON_METHOD) {
      fprintf(stderr, "Error: You may not specify '--method pearson' and '--control'.\n");
      bad_options = true;
    }
    if (options.pvalue_method == SPEARMAN_METHOD) {
      fprintf(stderr, "Error: You may not specify '--method spearman' and '--control'.\n");
      bad_options = true;
    }
    if (options.pvalue_method == MULTIHG_METHOD) {
      fprintf(stderr, "Error: You may not specify '--method 3dmhg' and '--control'.\n");
      bad_options = true;
    }
    if (options.pvalue_method == LONG_MULTIHG_METHOD) {
      fprintf(stderr, "Error: You may not specify '--method 4dmhg' and '--control'.\n");
      bad_options = true;
    }
    if (options.fix_partition > 0) {
      fprintf(stderr, "Error: You may not specify '--fix-partition' and '--control'.\n");
      bad_options = true;
    }
    if (options.positive_list == POS_PWM) {
      fprintf(stderr, "Error: You may not specify '--poslist pwm' and '--control'.\n");
      bad_options = true;
    }
  }

  if (options.fisher_fasta_threshold != default_options.fisher_fasta_threshold) {
    if (options.pvalue_method != FISHER_METHOD && options.positive_list != POS_PWM) { 
      fprintf(stderr, "Error: '--fasta-threshold' requires '--method fisher' and '--poslist pwm'.\n");
      bad_options = true;
    }
  }

  if (options.positive_list == POS_PWM) {
    if (options.fix_partition > 0) {
      fprintf(stderr, "Error: You may not specify '--fix-partition' and '--poslist pwm'.\n");
      bad_options = true;
    }
  } 

  if (options.motif_sources == NULL) {
    fprintf(stderr, "Error: Motif file not specified.\n");
    bad_options = true;
  } else {
    int i;
    int number_motif_files = arraylst_size(options.motif_sources);
    for (i = 0; i < number_motif_files; i++) {
      if (!file_exists(arraylst_get(i, options.motif_sources))) {
        fprintf(stderr, "Error: Specified motif '%s' file does not exist.\n", (char *) arraylst_get(i, options.motif_sources));
        bad_options = true;
      }
    }
  }

  if (options.seq_source == NULL) {
    fprintf(stderr, "Error: Sequence file not specified.\n");
    bad_options = true;
  } else if (!file_exists(options.seq_source)) {
    fprintf(stderr, "Error: Specified sequence file '%s' does not exist.\n", options.seq_source);
    bad_options = true;
  }

  if (options.scoring != TOTAL_HITS) {
    if ((options.pvalue_method == MULTIHG_METHOD || options.pvalue_method == LONG_MULTIHG_METHOD ) ) {
      fprintf(stderr, "Error: You must use '--scoring totalhits' with '--method 3dmhg' or '--method 4dmhg'.\n");
      bad_options = true;
    }
  }

  if (options.log_pwmscores && options.pvalue_method != PEARSON_METHOD) {
    fprintf(stderr, "Error: You may only specify '--log-pwmscores' with '--method pearson'.\n");
    bad_options = true;
  }

  if (options.log_fscores && options.pvalue_method != PEARSON_METHOD) {
    fprintf(stderr, "Error: You may only specify '--log-fscores' with '--method pearson'.\n");
    bad_options = true;
  }

  if (bad_options) {
    ame_usage("Type 'ame --h' to see the complete command usage message.\n");
    ame_terminate(&options, 1);
  }

  // make enough space for all the command line options, with one space between each
  int line_length = 0;
  for (i = 0; i < argc; i++)
    line_length += strlen(i == 0 ? basename(argv[0]) : argv[i]);
  // add on argc to allow one char per word for separating space + terminal '\0'
  options.commandline = (char*) mm_malloc(sizeof(char)*(line_length+argc));
  int nextpos = 0;
  for (i = 0; i < argc; i++) {
    // been here before? put in a space before adding the next word
    if (nextpos) {
      options.commandline[nextpos] = ' ';
      nextpos ++;
    }
    char * nextword = i == 0 ? basename(argv[0]) : argv[i];
    strcpy(&options.commandline[nextpos], nextword);
    nextpos += strlen (nextword);
  }

  // Sequences are only output if the method is FISHER.
  if (FISHER_METHOD != options.pvalue_method) options.noseq = true;

  // Open the output files unless --text.
  if (! options.text_only) {
    if (create_output_directory(options.outputdir, options.clobber, 
       options.verbose > QUIET_VERBOSE)) { // only warn in higher verbose modes
      die("Failed to create output directory `%s' or already exists.\n", options.outputdir);
    }
    char *path = make_path_to_file(options.outputdir, TSV_FILENAME);
    options.tsv_output = fopen(path, "w"); //FIXME check for errors: MEME doesn't either and we at least know we have a good directory
    myfree(path);
    if (! options.noseq) {
      path = make_path_to_file(options.outputdir, SEQ_FILENAME);
      options.seq_output = fopen(path, "w"); //FIXME check for errors: MEME doesn't either and we at least know we have a good directory
    } 
    myfree(path);
  }

} // ame_getopt

/*************************************************************************
 * Correct a p-value for multiple tests.
 * If log_p is too small, the power forumula rounds to 0.
 *************************************************************************/
double inline ame_bonferroni_correction(double log_p, double numtests) {
  return LOGEV(log(numtests), log_p);
} // ame_bonferroni_correction

/*************************************************************************
 * Read in the motif databases.
 *************************************************************************/
ARRAYLST_T *ame_load_motifs_and_background(
  bool separate_namespaces,	// keep motif DB namespaces separate
  bool xalph,			// convert motifs to alphabet specified in options
  ARRAY_T **background,		// OUT the background
  int *max_width		// OUT maximum motif width
) 
{
  int i, num_motifs = 0;
  ARRAYLST_T *dbs = arraylst_create_sized(arraylst_size(options.motif_sources));
  RBTREE_T *motif_names = rbtree_create(rbtree_strcmp, NULL, NULL, NULL, NULL);

  ALPH_T *alph = options.alphabet;
  assert(alph != NULL);
  bool use_rc = alph_has_complement(alph);

  bool stdin_used = false;	// set if path is "-"
  *max_width = 0;		// maximum motif width
  for (i = 0; i < arraylst_size(options.motif_sources); i++) {

    // Get the name of the next motif db.
    char *motif_source = (char*) arraylst_get(i, options.motif_sources);
    DEBUG_FMT(NORMAL_VERBOSE, "Loading motifs from file '%s'\n", motif_source);

    // Load the motifs from this file.
    MOTIF_DB_T *db = read_motifs_and_background(
      i, 			// id of DB file
      motif_source,		// motif file name (or special word)
      "Query motifs",		// type of database for error messages
      NULL,			// get one motif by name
      NULL,			// get one motif by index
      options.include_patterns,	// get set of motifs by name (or NULL)
      options.exclude_patterns,	// exclude this set of motifs by name (or NULL)
      true,			// allow motifs with zero probability entries
      false,			// don't create RC copies, appended
      options.pseudocount,	// multiply times background model
      false,			// set_trim
      0,			// trim_bit_threshold
      &(options.bg_source),	// background file; may be changed
      true, 			// make bg symmetrical if alph complementable
      background,		// will be set if id==0
      options.seq_source,	// sequence file name
      alph,			// sequence alphabet
      xalph,			// set motif conversion alphabet
      false,			// don't remove extension from name
      false,			// don't remove ".meme" extension from name
      false,			// don't replace underscores in name
      &stdin_used		// IN/OUT check and set if path is "-"
    );

    // Get maximum motif width and check for motif name uniqueness across all DBs.
    int warn_type = NO_WARNING;
    int i;
    for (i=0; i<arraylst_size(db->motifs); i++) {
      MOTIF_T *motif = arraylst_get(i, db->motifs);
      if (! separate_namespaces) {
	bool created;
	RBNODE_T *node = rbtree_lookup(motif_names, get_motif_id(motif), true, &created);
	if (!created) {
	  clump_motif_db_warning(&warn_type, DUPLICATE_WARNING, "Warning: The following "
	    "duplicate motifs in '%s' were excluded:\n  ", motif_source, get_motif_id(motif));
	  destroy_motif(motif);
	  db->loaded--;
	  db->excluded++;
	  continue;
	}
      } // namespaces
      // Get width of widest motif.
      int w = get_motif_length(motif);
      if (w > *max_width) *max_width = w;
    }

    // Add DBs to the list of DBs
    arraylst_add(db, dbs);

    // Number of motifs found.
    num_motifs += arraylst_size(db->motifs);

    if (warn_type && verbosity >= NORMAL_VERBOSE) fprintf(stderr, "\n");
  } // motif_files

  // Check that we found suitable motifs.
  if (num_motifs == 0) die("No acceptable motifs found.");

  // Cleanup
  rbtree_destroy(motif_names);

  return(dbs);
} // ame_load_motifs_and_background

/*************************************************************************
 * Make a PSSM for AME.
 *************************************************************************/
static PSSM_T* ame_make_pssm(
  bool use_log_odds,		// type of PSSM to make; true->log-odds, false->_odds
  ARRAY_T* bg_freqs,		// background model
  MOTIF_T *motif		// motif (contains frequencies)
) {
  PSSM_T *pssm;
  // Build PSSM for motif and tables for p-value calculation.
  pssm = build_motif_pssm(
    motif,		// motif frequencies p_ia (IN)
    bg_freqs,		// background frequencies b_a for pssm (IN)
    bg_freqs,		// background frequencies b_a for p-values (IN)
    NULL,		// Distribution of priors. May be NULL (IN)
    0,			// Scale factor for non-specific priors. Unused if prior_dist is NULL.
    100,		// range of scaled scores is [0..w*range]
    0,			// no GC bins
    ! use_log_odds	// no_log
  );
  return pssm;
} // ame_make_pssm

/*************************************************************************
 * Calculate the odds score for each motif-sized window at each
 * site in the sequence using the given nucleotide frequencies.
 *
 * This function is a lightweight version based on the one contained in
 * motiph-scoring. Several calculations that are unnecessary for GOMo
 * have been removed in order to speed up the process
 *************************************************************************/
static double ame_score_sequence(
  SEQ_T *seq,		// sequence to scan (IN)
  MOTIF_T *motif,	// motif (IN)
  PSSM_T *pssm,		// PSSM (may be odds or scaled log-odds (IN)
  int method,		// method used for scoring (IN)
  double threshold,	// threshold to use in TOTAL_HITS mode
  ARRAY_T *bg_freqs	// background model
) {
  assert(seq != NULL);
  assert(motif != NULL);

  ALPH_T *alphabet = get_motif_alph(motif);
  int asize = get_pssm_alphsize(pssm);
  int w = get_pssm_w(pssm);

  char* raw_seq = get_raw_sequence(seq);
  int seq_length = get_seq_length(seq);

  int max_index = seq_length - get_motif_length(motif);
  if (max_index < 0) max_index = 0;
  double *scores = (double*) mm_malloc(sizeof(double)*max_index);

  // For each site in the sequence
  int seq_index;
  for (seq_index = 0; seq_index < max_index; seq_index++) {
    double score = 1;

    // For each site in the motif window
    int motif_position, aindex;
    for (motif_position = 0; motif_position < w; motif_position++) {
      char c = raw_seq[seq_index + motif_position];

      // Skip gaps and ambiguous characters.
      if (c == '-' || c == '.') break;
      aindex = alph_index(alphabet, c);
      if (aindex >= asize) break;
      score *= get_pssm_score(motif_position, aindex, pssm);
    }
    scores[seq_index] = score;
  }

  // return odds as requested (MAX or AVG scoring)
  double requested_odds = 0.0;
  if (method == AVG_ODDS || method == SUM_ODDS) {
    for (seq_index = 0; seq_index < max_index; seq_index++) {
      requested_odds += scores[seq_index];
    }
    if (method == AVG_ODDS) requested_odds /= max_index + 1;
  } else if (method == MAX_ODDS) {
    for (seq_index = 0; seq_index < max_index; seq_index++) {
      if (scores[seq_index] > requested_odds) requested_odds = scores[seq_index];
    }
  } else if (method == TOTAL_HITS) {
    for (seq_index = 0; seq_index < max_index; seq_index++) {
      if (scores[seq_index] >= threshold) requested_odds++; 	// It's a hit.
      DEBUG_FMT(HIGHER_VERBOSE, "Window Data: %s\t%s\t%i\t%g\t%g\n",
      get_seq_name(seq), get_motif_id(motif), seq_index, scores[seq_index], threshold);
    }
  } else {
    die("Unknown value for --scoring in ame_score_sequence.");
  }

  myfree(scores);

  return(requested_odds);
} // ame_score_sequence

/*************************************************************************
 * Scan a sequence with the given motif (and its RC).
 *************************************************************************/
double ame_sequence_scan(
    SEQ_T* sequence,    // the sequence to scan INPUT
    MOTIF_T* motif,     // the motif to scan with INPUT
    MOTIF_T* rc_motif,	// the RC motif
    PSSM_T* pssm,	// the PSSM
    PSSM_T* rc_pssm,	// the RC PSSM
    int scoring,        // the scoring function to apply AVG_ODDS, MAX_ODDS or TOTAL_HITS
    bool scan_both_strands, // should we scan with both motifs and combine scores
    double threshold,   // threshold to use in TOTAL_HITS mode with a PWM
    ARRAY_T* bg_freqs   //background model
) {
  // Score the forward strand.
  double score = ame_score_sequence(
    sequence,
    motif,
    pssm,
    scoring,
    threshold,
    bg_freqs
  );
  // Score the reverse strand.
  if (scan_both_strands) {
    double rc_score = ame_score_sequence(
      sequence,
      rc_motif,
      rc_pssm,
      scoring,
      threshold,
      bg_freqs
    );
    if (scoring == AVG_ODDS){
      score = (score + rc_score)/2.0;
    } else if (scoring == MAX_ODDS){
      score = max(score, rc_score);
    } else if (scoring == SUM_ODDS){
      score = score + rc_score;
    } else if (scoring == TOTAL_HITS) {
      score = score + rc_score;
    } else {
      die("Unknown value for --scoring in ame_sequence_scan.");
    }
  } // both strands

  return(score);
} // ame_sequence_scan

/*************************************************************************
 * Scan the sequences with the given motif.  Returns the list of scores.
 *************************************************************************/
double *ame_scan_sequences(
  double min_score_thresh, // the hit threshold
  int num_seqs,		// the number of sequences
  SEQ_T **seq_list,	// the sequences
  ARRAY_T *bg_freqs,	// background model
  int motif_idx,	// index of motif (0..)
  MOTIF_T *motif,	// the motif
  PSSM_T *pssm,		// the PSSM
  MOTIF_T *rc_motif,	// the reverse complement motif (or NULL)
  PSSM_T *rc_pssm 	// the reverse complement PSSM (or NULL)
) {
  int i;
  bool use_rc = rc_motif != NULL;
  //double scaled_pvalue_threshold = options.pvalue_threshold;
  double *scores = mm_malloc(sizeof(double) * num_seqs);  // list of sequence scores

  // Score each sequence.
  for (i=0; i<num_seqs; i++) {
    SEQ_T *sequence = seq_list[i];
    int seq_len = get_seq_length(sequence);

    // Progress counter.
    if ((i+1)%1000==0 || i+1==num_seqs) DEBUG_FMT(NORMAL_VERBOSE, "MOTIF: %d SEQ: %d/%d\r", motif_idx+1, i+1, num_seqs);

    // Score the sequence.
    scores[i] = ame_sequence_scan(
      sequence,		// the sequence to scan INPUT
      motif,		// the motif to scan with INPUT
      rc_motif,		// the RC motif
      pssm,		// the PSSM
      rc_pssm,		// the RC PSSM
      options.scoring,	// the scoring function to apply AVG_ODDS, MAX_ODDS or TOTAL_HITS
      use_rc, 		// Should we scan with both motifs and combine scores
      min_score_thresh,	// threshold to use in TOTAL_HITS mode with a PWM
      bg_freqs		//background model
    );

  } // sequence

  if (verbosity >= NORMAL_VERBOSE && motif_idx==0) fprintf(stderr, "\n");

  return(scores);
} // ame_scan_sequences

/*************************************************************************
 * Pearson CC or Spearman CC.
 *************************************************************************/
AME_RESULT_T *ame_do_pearson_test(
  AME_RANK_T** rankings, 
  int num_seqs,
  int db_id,
  MOTIF_T *motif,
  bool use_ranks
) {
  // Assorted vars
  int i;

  // Vars for the regression
  double *x;
  double *y;

  // Vars for scoring
  AME_RESULT_T *lowest_motif_result;

  // Allocate memory or set initial values.
  // Allocate space, as a ptr to this will go in the array later.
  lowest_motif_result = mm_malloc(sizeof(AME_RESULT_T)); 
  //that's why we don't free it in this loop.
  x = mm_malloc(sizeof(double)*num_seqs);
  y = mm_malloc(sizeof(double)*num_seqs);

  // Now we need to copy the scores into two double arrays
  // Use LOG macro so that log(0) 'works'
  // Use times -x so pearson_correlation will give significance level
  // for positive correlation coefficient.
  for (i=0; i < num_seqs; i++) {
    if (options.log_fscores == true) {
      x[i] = -LOG(rankings[i]->f_score);
    } else {
      x[i] = use_ranks ? rankings[i]->f_rank : -rankings[i]->f_score;
    }

    if (options.log_pwmscores == true) {
      y[i] = LOG(rankings[i]->pwm_score);
    } else {
      y[i] = use_ranks ? rankings[i]->pwm_rank : rankings[i]->pwm_score;
    }

    if (options.verbose >= HIGHER_VERBOSE) {
      fprintf(stderr, "Rank %i: LR F-Score %.3g (%.3g) LR motif-Score: %.3g "
	"(%.3g)\n", i, x[i], rankings[i]->f_score, y[i], 
	rankings[i]->pwm_score);
    }
  }

  // We start with a minimum of three sequences so that the data
  // is over-described.
  int min = 3;
  int max = num_seqs;
  if (options.fix_partition > 0) min = options.fix_partition;
  if (min < 3) min = 3;
  if (options.fix_partition > 0) max = min;
  if (max > num_seqs) max = num_seqs;
  double lowest_mse = 10000;
  double lowest_log_pv = 1;
  int num_tests = 0;			// number of statistical tests
  for (i=min; i <= max; i++) {
    double m = 0;
    double b = 0;
    double rho = 0;
    double log_pv = 0;
    double mse = 0;
    if (options.linreg_switchxy) {
      if (use_ranks) {
        rho = pearson_correlation(i, y, x, NULL, NULL, NULL, &log_pv, false, NULL);
      } else {
        rho = pearson_correlation(i, y, x, &m, &b, &mse, &log_pv, false, NULL);
      }
    } else {
      if (use_ranks) {
	rho = pearson_correlation(i, x, y, NULL, NULL, NULL, &log_pv, false, NULL);
      } else {
	rho = pearson_correlation(i, x, y, &m, &b, &mse, &log_pv, false, NULL);
      }
    }

    // fix NANs
    if (isnan(mse) || isnan(rho) || log_pv > 0) {
      mse = 9999;	// smaller than initial value to make sure a value is stored in motif name
      log_pv = 0;	// smaller ""
      rho = 0;
      m = b = 0;
    }

    // Undo multiplication of -1 * X.
    m = -m;

    if (options.verbose >= HIGH_VERBOSE) {
      fprintf(stderr, "%s p-value of motif '%s' top %i seqs: ",
        use_ranks ? "Spearman CC" : "LinReg", get_motif_st_id(motif), i);
      print_log_value(stderr, log_pv, 2);
      fprintf(stderr, " (m: %.3g b: %.3g rho: %.3g mse: %.3g)\n", m, b, rho, mse);
    }

    // Add to our motif list if best so far.
    double p_thresh = (POS_FASTA == options.positive_list) ? rankings[i-1]->f_score : rankings[i-1]->pwm_score;
    if (i < MIN_PEARSON_SAMPLES) {	// not enough samples for p-values yet; optimize MSE
      if (mse < lowest_mse) {
	lowest_mse = mse;
        init_result(
          lowest_motif_result, 
          db_id, 
          motif, 
	  i, // pos
	  0, // neg; non-positives ignored
	  -1, // tp
	  -1, // fp
	  p_thresh, // p_thresh
	  0, // tp_thresh
	  m, // log_pleft: m (slope)
	  0, // log_pright: log_pv = log(1)
	  b, // log_pboth: b (intercept)
	  mse, // u: mse
	  rho,  // rho 
	  0 // num_tests
	);
      }
    } else {			// enough samples; optimize significance
      num_tests++;
      if (log_pv < lowest_log_pv) {
	lowest_log_pv = log_pv;
        init_result(
          lowest_motif_result, 
          db_id,
          motif,
	  i, // pos
	  0, // neg; non-positives ignored
	  -1, // tp
	  -1, // fp
	  p_thresh, // p_thresh
	  0, // tp_thresh
	  m, // log_pleft: m (slope)
	  log_pv, // log_pright: log_pv
	  b, // log_pboth: b (intercept)
	  mse, // u: mse
	  rho,  // rho 
	  0 // num_tests
	);
      }
    }
  }

  // Free the arrays.
  free(x);
  free(y);

  // Record the number of tests.
  lowest_motif_result->num_tests = max(1, num_tests);
  return lowest_motif_result;
} // ame_do_pearson_test

/*************************************************************************
 * Ranksum test.
 *************************************************************************/
AME_RESULT_T *ame_do_ranksum_test(
  AME_RANK_T** rankings, 
  int num_seqs,
  double *scores,			// the sequence scores
  int db_id,
  MOTIF_T *motif
) {
  int n,na; // for rs stats test.
  double ta_obs; // for rs stats test.
  int i;
  double lowest_log_pval;
  RSR_T *r;
  AME_RESULT_T *lowest_motif_result;

  n = num_seqs;
  na = ta_obs = 0;
  lowest_log_pval = 0; 	// 0 is largest possible log_pval
  lowest_motif_result = mm_malloc(sizeof(AME_RESULT_T)); 

  // The top of the rankings list are considered "positives".
  int num_tests = 0;			// number of statistical tests
  for (i=0; i < num_seqs; i++) {
     /*
     * We need to keep track of:
     *      int n,          number of samples
            int na,         number of positives
            double ta_obs   sum of ranks of positives
     */
    na++;
    // Use the rank according to the "other" score (the one not used for sorting).
    if (POS_PWM == options.positive_list) {
      ta_obs += rankings[i]->f_rank;
    } else {
      ta_obs += rankings[i]->pwm_rank;
    }

    // Compute p-value if at allowed partition point.
    if (options.fix_partition <= 0 || i == options.fix_partition-1) {
      if (QUICK_RS == options.rs_method) {
	// Note: 1.0 is needed to force float arithmetic to avoid overflow
	r = ranksum_from_stats(n, n - na, n*(n+1.0)/2 - ta_obs);
      } else {
	r = ranksum_from_sets(scores+i+1, num_seqs-(i+1), scores, i+1);
      }

      if (options.verbose >= HIGH_VERBOSE) {
        double ta_neg = n*(n+1.0)/2 - ta_obs;
	fprintf(stderr, "Ranksum p-values of motif '%s' top %i seqs "
          "(left,right,twotailed): %g %g %g U-value: %.4g Ranksum pos %.4g neg %.4g\n",
	  get_motif_st_id(motif), i+1,
	  exp(RSR_get_log_p_left(r)), exp(RSR_get_log_p_right(r)), exp(RSR_get_log_p_twotailed(r)), RSR_get_u(r), ta_obs, ta_neg);
      }

      // Add to our motif list
      num_tests++;
      double p_thresh = (POS_FASTA == options.positive_list) ? rankings[i]->f_score : rankings[i]->pwm_score;
      if (lowest_log_pval >= RSR_get_log_p_right(r)) {
	lowest_log_pval = RSR_get_log_p_right(r);
	init_result(
          lowest_motif_result, 
          db_id,
          motif,
	  i+1, // pos
	  0, // neg; non-positives ignored
	  -1, // tp
	  -1, // fp
	  p_thresh, // p_thresh
	  0, // tp_thresh
	  RSR_get_log_p_left(r), // log_pleft
	  RSR_get_log_p_right(r), // log_pright
	  RSR_get_log_p_twotailed(r), // log_pboth
	  RSR_get_u(r), // u
	  -2, // rho 
	  0 // num_tests
        );
      }
      destroy_rsr(r);
      // Exit loop if fixed partition
      if (options.fix_partition > 0) break;  
    } // possible partition

  } // split

  // Record the number of tests.
  lowest_motif_result->num_tests = max(1, num_tests);

  return lowest_motif_result;
} // ame_do_ranksum_test

/*************************************************************************
 * Fisher's exact test
 *************************************************************************/
AME_RESULT_T *ame_do_fisher_test(
  AME_RANK_T** rankings, 
  int num_seqs,
  int db_id,
  double min_pwm_thresh,
  MOTIF_T *motif
) {
  int i;
  double lowest_log_pval;
  double log_p=0,log_pleft=0,log_pright=0,log_pboth=0;
  RSR_T* r;
  AME_RESULT_T* lowest_motif_result;
  int tp,fn,fp,tn;
  double p_thresh = BIG, tp_thresh;

  // Allocate space, as a ptr to this will go in the array later.
  // That's why we don't free it in this loop.
  lowest_motif_result = mm_malloc(sizeof(AME_RESULT_T)); 

  //Used for testing so that we can be sure that we've initialised this variable
  lowest_motif_result->u = 0;

  lowest_log_pval = 100;		// A large number.

  int num_tests = 0;			// number of statistical tests
  if (options.fix_partition > 0) {
    // Fixed-partition mode.
    // Positive sequences are those with low FASTA rank.
    // Sequences are sorted by PWM score.
    // Find best motif threshold.
    int pos = options.fix_partition;
    int neg = num_seqs - pos;
    tp = fp = 0; 		// counts
    int split = options.fix_partition;	// number of positive sequences
    for (i=0; i<num_seqs; i++) {
      // Count positive and negative sequences with scores above the motif threshold
      if (rankings[i]->f_rank <= split) {
        tp++;
      } else {
        fp++;
      }
      // Add to our motif list if most significant so far.
      // --poslist is not allowed with --fix-partition or control
      // so p_thresh is on FASTA and tp_thresh is on PWM scores.
      if (rankings[i]->f_rank == split) p_thresh = rankings[i]->f_score; // save true p_thresh
      tp_thresh = rankings[i]->pwm_score;
      if (i>0 && tp_thresh < min_pwm_thresh) continue;
      log_pright = getLogFETPvalue(tp, pos, fp, neg, false);
      num_tests++;
      if (lowest_log_pval >= log_pright) {
	lowest_log_pval = log_pright;
	init_result(
          lowest_motif_result, 
          db_id,
          motif,
	  pos,	// pos
	  neg,	// neg
	  tp,	// tp
	  fp, 	// fp
	  p_thresh, // p_thresh
	  tp_thresh, // tp_thresh
	  0, // log_pleft
          lowest_log_pval, // log_pright
          0, // log_pboth 
          -1, // u
          -2,  // rho
	  0 // num_tests
        );
      }
      if (options.verbose >= HIGH_VERBOSE) {
	fprintf(stderr, "M3: %s Threshold: %i tp: %i fn: %i tn: %i fp: %i P ",
	  get_motif_st_id(motif), i+1, tp, pos-tp, neg-fp, fp);
	print_log_value(stderr, log_pright, 2);
        fprintf(stderr, " lowest P ");
	print_log_value(stderr, lowest_log_pval, 2);
	fprintf(stderr, "\n");
      }
    } // sequence

  } else {
    // Partition-maximization mode.
    // Loop over possible partitions.
    // Sequences are sorted by FASTA or by PWM score.
    tp = fn = fp = tn = 0;	// counts
    int split;
    for (split=1; split<=num_seqs-1; split++) {	// rank of split
      if (split==1) {			// initialize first time thru
	for (i=0; i<num_seqs; i++) {
	  if (i < split) {
	    // Positive sequence
	    if (POS_FASTA == options.positive_list) {
	      (rankings[i]->pwm_score >= min_pwm_thresh) ? tp++ : fn++;
	    } else {
	      (rankings[i]->f_score <= options.fisher_fasta_threshold) ? tp++ : fn++;
	    }
	  } else {
	    // Negative sequence
	    if (POS_FASTA == options.positive_list) {
	      (rankings[i]->pwm_score >= min_pwm_thresh) ? fp++ : tn++;
	    } else {
	      (rankings[i]->f_score <= options.fisher_fasta_threshold) ? fp++ : tn++;
	    }
	  }
	} // sequence
      } else {			// use DP to compute next values
	i = split - 1;		// update values for this sequence
	if (POS_FASTA == options.positive_list) {
	  // It used to be a negative sequence:
	  (rankings[i]->pwm_score >= min_pwm_thresh) ? fp-- : tn--;
	  // It now is a positive sequence:
	  (rankings[i]->pwm_score >= min_pwm_thresh) ? tp++ : fn++;
	} else {
	  // It used to be a negative sequence:
	  (rankings[i]->f_score <= options.fisher_fasta_threshold) ? fp-- : tn--;
	  // It now is a positive sequence:
	  (rankings[i]->f_score <= options.fisher_fasta_threshold) ? tp++ : fn++;
	}
      }

      // Add to our motif list if most signficant so far.
      if (POS_FASTA == options.positive_list) {
        p_thresh = rankings[split-1]->f_score;
        tp_thresh = min_pwm_thresh;
      } else {
        p_thresh = rankings[split-1]->pwm_score;
        tp_thresh = options.fisher_fasta_threshold;
      }
      log_pright = getLogFETPvalue(tp, tp+fn, fp, tn+fp, false);
      num_tests++;
      if (lowest_log_pval >= log_pright) {
	lowest_log_pval = log_pright;
	init_result(
          lowest_motif_result, 
          db_id,
          motif,
	  split, // pos
	  num_seqs - split, // neg
	  tp, // tp
	  fp, // fp
	  p_thresh, // p_thresh
	  tp_thresh, // tp_thresh
	  0,
          lowest_log_pval, // log_pright
	  0,
          -1, // u
          -2,  // rho
	  0 // num_tests
        );
      }
      if (options.verbose >= HIGH_VERBOSE) {
	fprintf(stderr, "M3: %s Threshold: %i tp: %i fn: %i tn: %i fp: %i P ",
	  get_motif_st_id(motif), split, tp, fn, tn, fp);
	print_log_value(stderr, log_pright, 2);
	fprintf(stderr, "\n");
      }
    } // split
  }

  // Record the number of tests.
  lowest_motif_result->num_tests = max(1, num_tests);
  // Record the correct value of p_thresh if in fixed partition mode.
  if (options.fix_partition > 0) lowest_motif_result->p_thresh = p_thresh;
  assert(lowest_motif_result != NULL);
  assert(lowest_motif_result->u != 0);
  return lowest_motif_result;
} // ame_do_fisher_test

/*************************************************************************
 * 3- and 4-D multiHG test.
 *************************************************************************/
AME_RESULT_T *ame_do_multihg_test(
  AME_RANK_T** rankings, 
  int num_seqs,
  int db_id,
  MOTIF_T *motif
) {
  int i,j,k;
  double lowest_log_pval;
  double log_p;
  RSR_T* r;
  AME_RESULT_T* lowest_motif_result;
  int b0,b1,b2,b3;
  int n,N;
  int i0,i1,i2,i3,B0,B1,B2,B3;

  lowest_log_pval = 0; //a large number.
  lowest_motif_result = mm_malloc(sizeof(AME_RESULT_T)); //allocate space, as a ptr to this will go in the array later
  //that's why we don't free it in this loop.

  //Used for testing so that we can be sure that we're initialised this variable
  lowest_motif_result->u = 0;

  int num_tests = 0;			// number of statistical tests

  //Get per-class totals
  B0 = B1 = B2 = B3 = 0;
  for (i=0; i < num_seqs; i++) {
    if (rankings[i]->pwm_score == 0) {
      B0++;
    } else if (rankings[i]->pwm_score == 1) {
      B1++;
    } else if (rankings[i]->pwm_score == 2) {
      B2++;
    } else {
      B3++;
    }
  }
  // If we are only 3d rather than 4d...
  if (options.pvalue_method == MULTIHG_METHOD) {
    B2+=B3;
  }

  /*
   * Do it for a fixed partition only.
   */
  if (options.fix_partition > 0) {
    int i = options.fix_partition - 1;
    /* Need to do count table:
     *
     *     a0 - above threshold, 0 hits
     *     a1 - above threshold, 1 motif hit
     *     a2 - above threshold, >=2 motif hits
     *
     */
    b0 = b1 = b2 = b3 = 0;
    for (j=0;j<=i; j++) {
      if (rankings[j]->pwm_score == 0) {
        b0++;
      } else if (rankings[j]->pwm_score == 1) {
        b1++;
      } else if (rankings[j]->pwm_score == 2) {
        b2++;
      } else {
        b3++;
      }
    }

    if (options.pvalue_method == MULTIHG_METHOD) {

      b2 += b3;
      /* Formula details:
       *
       * Please see page 31 of Robert's lab book.
       *
       *
       * The notation is the same as Eden et al. 2007.
       *
       * If we count two or more. (3>2, so we add it to B2 etc), since this
       * piece of code has a max dimension of 3.
       *
       */
      n = i+1; //n is the threshold
      N = num_seqs; //total set size;
      log_p = LOGZERO;

      for (i0=b0; i0<=B0 && i0<=n; i0++) {
        for (i1=b1; i1<=B1 && i1<=n-i0; i1++) {
          for (i2=b2; i2<=B2 && i2<=n-(i0+i1); i2++) {
            log_p = LOG_SUM(log_p, //We're in log space, remember.
                (
                 (
                  factorial_array[n] - (factorial_array[i0] + 
                    factorial_array[i1] + factorial_array[i2])
                 ) + (
                   (factorial_array[N-n]) - (factorial_array[B0-i0] + 
                     factorial_array[B1-i1] + factorial_array[B2-i2])
                   )
                ) - (
                  (factorial_array[N] - (factorial_array[B0] + 
                      factorial_array[B1] + factorial_array[B2]))
                  )
                );
          }
        }
      }

    } else { // We're going to do the complicated method with four loops.

      /* Formula details:
       *
       * Please see page 31 of Robert's lab book.
       *
       *
       * The notation is the same as Eden et al. 2007.
       *
       *
       */
      n = i+1; //n is the threshold
      N = num_seqs; //total set size;
      log_p = LOGZERO;

      for (i0=b0; i0<=B0 && i0<=n; i0++) {
        for (i1=b1; i1<=B1 && i1<=n-i0; i1++) {
          for (i2=b2; i2<=B2 && i2<=n-(i0+i1); i2++) {
            for (i3=b3; i3<=B3 && i3<=n-(i0+i1+i2); i3++) {
              log_p = LOG_SUM(log_p, //We're in log space, remember.
                  (
                   (
                    factorial_array[n] - (factorial_array[i0] + 
                      factorial_array[i1] + factorial_array[i2] + 
                      factorial_array[i3])
                   ) + (
                     (factorial_array[N-n]) - (factorial_array[B0-i0] + 
                       factorial_array[B1-i1] + factorial_array[B2-i2] + 
                       factorial_array[B3-i3])
                     )
                  ) - (
                    (factorial_array[N] - (factorial_array[B0] + 
                            factorial_array[B1] + factorial_array[B2]) + 
                     factorial_array[B3])
                    )
                  );
            }
          }
        }
      }
    }

    //Add to our motif list
    if (options.verbose >= HIGH_VERBOSE) {
      fprintf(stderr, "Motif: %s Threshold: %i P-value: ", get_motif_st_id(motif), i);
      print_log_value(stderr, log_p, 2);
      fprintf(stderr, "\n");
    }
    double p_thresh = (POS_FASTA == options.positive_list) ? rankings[i]->f_score : rankings[i]->pwm_score;
    lowest_log_pval = log_p;
    num_tests++;
    init_result(
      lowest_motif_result, 
      db_id,
      motif,
      i+1, // pos
      0, // neg; non-positives ignored
      -1, // tp
      -1, // fp
      p_thresh, // p_thresh
      0, 	// tp_thresh
      lowest_motif_result->log_pleft, // log_pleft
      log_p, // log_pright
      lowest_motif_result->log_pboth, // log_pboth
      -1, // u
      -2, // rho
      0 // num_tests
      );
  } else {

    //TODO: Make this not n^2
    for (i=0; i < num_seqs; i++) {
      /* Need to do count table:
       *
       *     a0 - above threshold, 0 hits
       *     a1 - above threshold, 1 motif hit
       *     a2 - above threshold, >=2 motif hits
       *
       */
      b0 = b1 = b2 = b3 = 0;
      for (j=0;j<=i; j++) {
        if (rankings[j]->pwm_score == 0) {
          b0++;
        } else if (rankings[j]->pwm_score == 1) {
          b1++;
        } else if (rankings[j]->pwm_score == 2) {
          b2++;
        } else {
          b3++;
        }
      }

      if (options.pvalue_method == MULTIHG_METHOD) {

        b2 += b3;
        /* Formula details:
         *
         * Please see page 31 of Robert's lab book.
         *
         *
         * The notation is the same as Eden et al. 2007.
         *
         * If we count two or more. (3>2, so we add it to B2 etc), since this
         * piece of code has a max dimension of 3.
         *
         */
        n = i+1; //n is the threshold
        N = num_seqs; //total set size;
        log_p = LOGZERO;

        for (i0=b0; i0<=B0 && i0<=n; i0++) {
          for (i1=b1; i1<=B1 && i1<=n-i0; i1++) {
            for (i2=b2; i2<=B2 && i2<=n-(i0+i1); i2++) {
              log_p = LOG_SUM(log_p, //We're in log space, remember.
                  (
                   (
                    factorial_array[n] - (factorial_array[i0] + 
                      factorial_array[i1] + factorial_array[i2])
                   ) + (
                     (factorial_array[N-n]) - (factorial_array[B0-i0] + 
                       factorial_array[B1-i1] + factorial_array[B2-i2])
                     )
                  ) - (
                    (factorial_array[N] - (factorial_array[B0] + 
                        factorial_array[B1] + factorial_array[B2]))
                    )
                  );
            }
          }
        }

      } else { // We're going to do the complicated method with four loops.

        /* Formula details:
         *
         * Please see page 31 of Robert's lab book.
         *
         *
         * The notation is the same as Eden et al. 2007.
         *
         *
         */
        n = i+1; //n is the threshold
        N = num_seqs; //total set size;
        log_p = LOGZERO;

        for (i0=b0; i0<=B0 && i0<=n; i0++) {
          for (i1=b1; i1<=B1 && i1<=n-i0; i1++) {
            for (i2=b2; i2<=B2 && i2<=n-(i0+i1); i2++) {
              for (i3=b3; i3<=B3 && i3<=n-(i0+i1+i2); i3++) {
                log_p = LOG_SUM(log_p, //We're in log space, remember.
                    (
                     (
                      factorial_array[n] - (factorial_array[i0] + 
                        factorial_array[i1] + factorial_array[i2] + 
                        factorial_array[i3])
                     ) + (
                       (factorial_array[N-n]) - (factorial_array[B0-i0] + 
                         factorial_array[B1-i1] + factorial_array[B2-i2] + 
                         factorial_array[B3-i3])
                       )
                    ) - (
                      (factorial_array[N] - (factorial_array[B0] + 
                              factorial_array[B1] + factorial_array[B2]) + 
                       factorial_array[B3])
                      )
                    );
              }
            }
          }
        }

      }

      //TODO :fix below here.
      //Add to our motif list
      if (options.verbose >= HIGH_VERBOSE) {
        fprintf(stderr, "Motif: %s Threshold: %i P-value: ", get_motif_st_id(motif), i);
        print_log_value(stderr, log_p, 2);
        fprintf(stderr, "\n"); 
      }
      double p_thresh = (POS_FASTA == options.positive_list) ? rankings[i]->f_score : rankings[i]->pwm_score;
      num_tests++;
      if (lowest_log_pval >= log_p) {
        lowest_log_pval = log_p;
        init_result(
          lowest_motif_result, 
          db_id,
	  motif,
	  i+1,			// pos
	  0,		 	// neg; non-positives are ignored
	  -1, 			// tp
	  -1, 			// fp
	  p_thresh,		// p_thresh
	  0,			// tp_thresh
	  lowest_motif_result->log_pleft, // log_pleft
	  log_p, // log_pright
	  lowest_motif_result->log_pboth, // log_pboth
          -1, // u
	  -2, // rho
           0 // num_tests
        );
      }
    }
  }

  // Record the number of tests.
  lowest_motif_result->num_tests = max(1, num_tests);
  assert(lowest_motif_result != NULL);
  assert(lowest_motif_result->u != 0);
  return lowest_motif_result;
} // ame_do_multihg_test

/*****************************************************************************
* Assign ranks to a set of scores in a rankings array.
*****************************************************************************/
void ame_assign_ranks(
  AME_RANK_T **rankings,	// array of rankings
  int num_seqs,			// the number of sequences (hence, rankings)
  int type,			// 1: FASTA 2: PWM
  bool handle_ties,		// assign ties the average rank if true
  bool first			// print status messages if true
) {
  int i, j;
  int prev_i;
  double prev_score;
  double tied_rank_sum;
  double avg_tied_rank;

  if (handle_ties) {
    // Handle tied scores by assigning average rank.
    if (type == 1) {	// FASTA
      if (first) DEBUG_MSG(NORMAL_VERBOSE, "Assigning tied FASTA scores the average rank.\n");
      prev_i = 0;
      prev_score = rankings[0]->f_score;
      tied_rank_sum = 1;
      for (i=1; i<num_seqs; i++) {
	if (prev_score == rankings[i]->f_score) {
	  tied_rank_sum += (i+1);
	} else {
	  double avg_tied_rank = tied_rank_sum / (i-prev_i);
	  for (j=prev_i; j<i; j++) rankings[j]->f_rank = avg_tied_rank;
	  prev_i = i;
	  prev_score = rankings[i]->f_score;
	  tied_rank_sum = (i+1);
	}
      }
      avg_tied_rank = tied_rank_sum / (i-prev_i);
      for (j=prev_i; j<i; j++) rankings[j]->f_rank = avg_tied_rank;
     } else {		// PWM
      if (first) DEBUG_MSG(NORMAL_VERBOSE, "Assigning tied PWM scores the average rank.\n");
      prev_i = 0;
      prev_score = rankings[0]->pwm_score;
      tied_rank_sum = 1;
      for (i=1; i<num_seqs; i++) {
	if (prev_score == rankings[i]->pwm_score) {
	  tied_rank_sum += (i+1);
	} else {
	  double avg_tied_rank = tied_rank_sum / (i-prev_i);
	  for (j=prev_i; j<i; j++) rankings[j]->pwm_rank = avg_tied_rank;
	  prev_i = i;
	  prev_score = rankings[i]->pwm_score;
	  tied_rank_sum = (i+1);
	}
      }
      avg_tied_rank = tied_rank_sum / (i-prev_i);
      for (j=prev_i; j<i; j++) rankings[j]->pwm_rank = avg_tied_rank;
    }
  } else {
    // Don't handle ties.
    if (type == 1) {	// FASTA
      for (i=0; i<num_seqs; i++) rankings[i]->f_rank = (i+1);
    } else {		// PWM
      for (i=0; i<num_seqs; i++) {
        rankings[i]->pwm_rank = (i+1);
      }
    }
  }
} // ame_assign ranks

/*****************************************************************************
* Get the enrichment score for a motif.
*****************************************************************************/
AME_RESULT_T *ame_get_motif_score(
  MOTIF_T *motif,		// the motif 
  int id,			// index of motif
  int db_id,			// index of DB
  char *db_source,		// filename of the motif DB
  int num_seqs,			// number of sequences
  SEQ_T **sequences,		// the sequences
  double *scores,		// sequence PWM scores for this motif
  double min_score_thresh,	// minumum score threshold for motif
  int num_motifs,		// the number of motifs (for E-value)
  bool fasta_scores,		// FASTA scores given?
  ARRAY_T *seq_fscores		// FASTA scores for sequences
) {
  int i;
  double highest_score = 0; 	// We use this to make sure that a PWM has scored on at least 1 seq.
				// Note that both odds scores and totalhits are >= 0.
  bool first = (id==0 && db_id==0);	// first motif ever?
  AME_RESULT_T *rsr = NULL;
  int handle_ties = (SPEARMAN_METHOD == options.pvalue_method || RANKSUM_METHOD == options.pvalue_method);
  int type;

  // Create the list of ranks.
  AME_RANK_T **rankings = mm_malloc(sizeof(AME_RANK_T*)*num_seqs);
  for (i=0; i<num_seqs; i++) {
    if (scores[i] > highest_score) highest_score = scores[i]; // keep track of the PWM high score.
    rankings[i] = mm_malloc(sizeof(AME_RANK_T));
    rankings[i]->seq_idx = i;
    rankings[i]->pwm_score = scores[i];
    // Set FASTA score to input order if no scores were given or they were all the same.
    rankings[i]->f_score = fasta_scores ? get_array_item(i, seq_fscores) : i+1;
    rankings[i]->f_rank = (i+1);	// initialize f-ranks to input sequence order
    rankings[i]->rand = drand_mt();	// random number
  }

  // Get the ranking according to the FASTA score if they were given.
  if (fasta_scores) {
    if (first) DEBUG_MSG(NORMAL_VERBOSE, "Sorting sequences by FASTA score to get FASTA ranks; breaking ties to conserve input order.\n");
    qsort(rankings, num_seqs, sizeof(AME_RANK_T*), ame_compare_ranks_fasta_score);
    type = 1;	// FASTA
    ame_assign_ranks(rankings, num_seqs, type, handle_ties, first);
  }

  // Get the ranking according to PWM score.
  // Break ties by putting *larger* f-ranks first so negatives come first.
  if (first) DEBUG_MSG(NORMAL_VERBOSE, 
    "Sorting sequences by sequence PWM score to get PWM ranks; breaking ties to put negatives first.\n");
  qsort (rankings, num_seqs, sizeof(AME_RANK_T*), ame_compare_ranks_pwm_score_fasta);
  type = 2;	// PWM
  ame_assign_ranks(rankings, num_seqs, type, handle_ties, first);

  // Re-sort via FASTA score if necessary.
  // Don't resort by FASTA if FISHER and fixed partition since we are
  // going to optimize over PWM score.
  if (POS_PWM == options.positive_list || 
    (POS_FASTA == options.positive_list && 
    (options.pvalue_method == FISHER_METHOD && (options.control_filename || options.fix_partition > 0)))) {
    // Don't re-sort via FASTA score.
    if (id == 0) DEBUG_MSG(NORMAL_VERBOSE, "Leaving sequences sorted by PWM score.\n");
  } else {
    // Re-sort via FASTA score.
    if (fasta_scores) {
      if (id == 0) DEBUG_MSG(NORMAL_VERBOSE, "Resorting sequences by FASTA score.\n");
    } else {
      if (id == 0) DEBUG_MSG(NORMAL_VERBOSE, "Resorting sequences to their original input order.\n");
    }
    qsort(rankings, num_seqs, sizeof(AME_RANK_T*), ame_compare_ranks_fasta_score);
  }

  // Debugging code.
  if (options.verbose >= HIGHER_VERBOSE) {
    fprintf(stderr, "\n");
    for (i=0; i < num_seqs; i++) {
      fprintf(stderr, "M2: %s - Seq: %s Rankings[%i] -\tpwm: %.8f\tprank: %f\tf: %.8f\tfrank: %f\n", 
	get_motif_st_id(motif), get_seq_name(sequences[rankings[i]->seq_idx]), i,
	rankings[i]->pwm_score, rankings[i]->pwm_rank, rankings[i]->f_score, rankings[i]->f_rank);
    }
  }

  // Describe the type of optimization.
  if (MULTIHG_METHOD != options.pvalue_method && LONG_MULTIHG_METHOD != options.pvalue_method) {
    if (id == 0) {
      if (FISHER_METHOD == options.pvalue_method
	&& (options.fix_partition > 0 || options.control_filename) ) {
	  DEBUG_MSG(NORMAL_VERBOSE, "Optimizing over sequence PWM score threshold.\n");
      } else if (options.fix_partition <= 0) {
	if (POS_FASTA == options.positive_list) {
	  if (fasta_scores) {
	    DEBUG_MSG(NORMAL_VERBOSE, "Performing partition maximization over FASTA score thresholds.\n");
	  } else {
	    DEBUG_MSG(NORMAL_VERBOSE, "Performing partition maximization over the sequence input order.\n");
	  }
	} else {
	  DEBUG_MSG(NORMAL_VERBOSE, "Performing partition maximization over sequence PWM score thresholds.\n"); 
	}
      }
    }
  }

  // Apply the selected method.
  if (highest_score > 0) {
    if (FISHER_METHOD == options.pvalue_method) {
      rsr = ame_do_fisher_test(rankings, num_seqs, db_id, min_score_thresh, motif);
    } else if (RANKSUM_METHOD == options.pvalue_method) {
      rsr = ame_do_ranksum_test(rankings, num_seqs, scores, db_id, motif);
    } else if (PEARSON_METHOD == options.pvalue_method || SPEARMAN_METHOD == options.pvalue_method) {
      bool use_ranks = SPEARMAN_METHOD == options.pvalue_method;
      rsr = ame_do_pearson_test(rankings, num_seqs, db_id, motif, use_ranks);
    } else if (MULTIHG_METHOD == options.pvalue_method || LONG_MULTIHG_METHOD == options.pvalue_method) {
      rsr = ame_do_multihg_test(rankings, num_seqs, db_id, motif);
    }
    // Set the adjusted p-value, E-value and if result is printable.
    rsr->log_corrected_pvalue = ame_bonferroni_correction(rsr->log_pright, rsr->num_tests);
    rsr->log_evalue = rsr->log_corrected_pvalue + log(num_motifs);
    rsr->is_significant = (rsr->log_evalue <= log(options.evalue_report_threshold));
    // Print the TP and FP sequences to the sequence TSV file.
    if (!options.noseq && rsr->is_significant) {
      output_tp_fp_sequences(
	first_sig_motif,// first motif ever
	db_source,	// filename of motif DB
	motif,		// the motif
	rsr,		// the motif's result
	num_seqs,	// number of sequences in rankings
	sequences,	// the sequences
	rankings	// sequence rankings
      ); 
      first_sig_motif = false;
    }
  } else {
    // If no sequence has scored at all, then we give a null result.
    rsr = init_result(
      NULL, 
      db_id,
      motif,
      num_seqs, // pos
      0, // neg
      0, // tp
      0, // fp
      -1, // p_thresh
      -1, // tp_thresh
      0, // log_pleft
      0, // log_pright: log_pv = log(1)
      0, // log_pboth
      0, // u
      0,  // rho
      1 // num_tests
    );
    rsr->log_corrected_pvalue = 0;
    rsr->log_evalue = 0;
    rsr->is_significant = false;
  }

  // Free up some space - TODO: Move this into a more appropriate place
  for (i=0; i<num_seqs; i++) free(rankings[i]);
  free(rankings);

  return(rsr);
} // ame_get_motif_score

/*****************************************************************************
* Output a motif result in TSV format.
*****************************************************************************/
void output_tsv_result(
  AME_RESULT_T *result,
  int rank,			// rank of result
  char *db_source		// name of motif DB
) {
  const char *RANK = "rank";
  const char *DB = "motif_DB";
  const char *ID = "motif_ID";
  const char *ID2 = "motif_alt_ID";
  const char *CONS = "consensus";
  const char *PVALUE = "p-value"; 
  const char *ADJ_PVALUE = "adj_p-value";
  const char *EVALUE = "E-value";
  const char *NUM_TESTS = "tests";
  const char *POS = "pos";
  const char *NEG = "neg";
  const char *TP = "TP";
  const char *TP_PCT = "%TP";
  const char *FP = "FP";
  const char *FP_PCT = "%FP";
  const char *FASTA = "FASTA_max";
  const char *PWM = "PWM_min";
  const char *PLEFT = "pleft";
  const char *PRIGHT = "pright";
  const char *PBOTH = "pboth";
  const char *U = "U";
  const char *ADJ_PLEFT = "adj_pleft";
  const char *ADJ_PRIGHT = "adj_pright";
  const char *ADJ_PBOTH = "adj_pboth";
  const char *PRHO = "Pearson_CC";
  const char *SRHO = "Spearmans_CC";
  const char *MSE = "mean_squared_error";
  const char *M = "slope";
  const char *B = "intercept";

  // Create the list of headers for the given method.
  ARRAYLST_T *headers = arraylst_create();
  arraylst_add((char *)RANK, headers);
  arraylst_add((char *)DB, headers);
  arraylst_add((char *)ID, headers);
  arraylst_add((char *)ID2, headers);
  arraylst_add((char *)CONS, headers);
  arraylst_add((char *)PVALUE, headers);
  arraylst_add((char *)ADJ_PVALUE, headers);
  arraylst_add((char *)EVALUE, headers);
  arraylst_add((char *)NUM_TESTS, headers);
  arraylst_add((POS_FASTA == options.positive_list) ? (char *)FASTA : (char *)PWM, headers);
  arraylst_add((char *)POS,  headers);
  arraylst_add((char *)NEG,  headers);
  if (FISHER_METHOD == options.pvalue_method) {
    arraylst_add((POS_FASTA == options.positive_list) ? (char *)PWM : (char *)FASTA, headers);
    arraylst_add((char *)TP,  headers);
    arraylst_add((char *)TP_PCT,  headers);
    arraylst_add((char *)FP,  headers);
    arraylst_add((char *)FP_PCT,  headers);
  } else if (RANKSUM_METHOD == options.pvalue_method) {
    arraylst_add((char *)U, headers); 
    arraylst_add((char *)PLEFT,headers);
    arraylst_add((char *)PRIGHT, headers);
    arraylst_add((char *)PBOTH, headers);
    arraylst_add((char *)ADJ_PLEFT, headers);
    arraylst_add((char *)ADJ_PRIGHT,headers);
    arraylst_add((char *)ADJ_PBOTH, headers);
  } else if (PEARSON_METHOD == options.pvalue_method) {
    arraylst_add((char *)PRHO, headers);
    arraylst_add((char *)MSE, headers);	
    arraylst_add((char *)M, headers);
    arraylst_add((char *)B, headers);
  } else if (SPEARMAN_METHOD == options.pvalue_method) {
    arraylst_add((char *)SRHO, headers);
  }

  // Print the header line.
  int i;
  if (rank==1) {
    for (i=0; i<arraylst_size(headers); i++) {
      tsvwr(i==0?"%s":"\t%s", arraylst_get(i, headers));
    }
    tsvwr("\n");
  }
  arraylst_destroy(NULL, headers);

  // Start result line.
  tsvwr("%i\t%s\t%s\t%s\t%s\t",
    rank,
    db_source,
    get_motif_id(result->motif),
    get_motif_id2(result->motif),
    get_motif_consensus(result->motif)
  );
  print_log_value(options.tsv_output, result->log_pright, 2);
  tsvwr("\t");
  print_log_value(options.tsv_output, result->log_corrected_pvalue, 2);
  tsvwr("\t");
  print_log_value(options.tsv_output, result->log_evalue, 2);
  tsvwr("\t%i", result->num_tests);
  if (options.control_filename) {
    tsvwr("\t%i", result->pos);
  } else {
    tsvwr("\t%.3g", result->p_thresh);
  }
  tsvwr("\t%i", result->pos);
  tsvwr("\t%i", result->neg);

  // Finish with extra columns for some methods.
  if (FISHER_METHOD == options.pvalue_method) {
    tsvwr("\t%.3g", result->tp_thresh);
    tsvwr("\t%i", result->tp);
    tsvwr("\t%.2f", result->pos > 0 ? (100.0 * result->tp)/result->pos : 0);
    tsvwr("\t%i", result->fp);
    tsvwr("\t%.2f", result->neg > 0 ? (100.0 * result->fp)/result->neg : 0);
  } else if (RANKSUM_METHOD == options.pvalue_method) {
    tsvwr("\t%.4g\t", result->u);
    print_log_value(options.tsv_output, result->log_pleft, 2);
    tsvwr("\t");
    print_log_value(options.tsv_output, result->log_pright, 2);
    tsvwr("\t");
    print_log_value(options.tsv_output, result->log_pboth, 2);
    tsvwr("\t");
    print_log_value(options.tsv_output, ame_bonferroni_correction(result->log_pleft, result->num_tests), 2);
    tsvwr("\t");
    print_log_value(options.tsv_output, ame_bonferroni_correction(result->log_pright, result->num_tests), 2); 
    tsvwr("\t");
    print_log_value(options.tsv_output, ame_bonferroni_correction(result->log_pboth, result->num_tests), 2); 
  } else if (PEARSON_METHOD == options.pvalue_method) {
    tsvwr("\t%.3g\t%.3g\t%.3g\t%.3g",
      result->rho, 		// correlation coefficient
      result->u,		// mse
      result->log_pleft,	// m
      result->log_pboth  	// b
    );
  } else if (SPEARMAN_METHOD == options.pvalue_method) {
    tsvwr("\t%.3g",
      result->rho  		// correlation coefficient
    );
  }
  tsvwr("\n");
} // output_tsv_result

/*************************************************************************
 * Finish outputting the results in HTML and TSV.
 *************************************************************************/
void ame_finish_output(
  AME_OPTION_T *options,		// the options
  ARRAYLST_T *dbs,		// the motif DBs
  int num_motifs,		// number of motifs with results
  AME_RESULT_T **rsr,		// results for each motif
  int num_seqs,			// the number of sequences
  double min_sequence_score,
  double max_sequence_score,
  double max_median_sequence_score
) {
  int i;

  // Output the min, max and max-median scores.
  if (!options->text_only) {
    jsonwr_dbl_prop(options->json_output, "min_sequence_score", min_sequence_score);
    jsonwr_dbl_prop(options->json_output, "max_sequence_score", max_sequence_score);
    jsonwr_dbl_prop(options->json_output, "max_median_sequence_score", max_median_sequence_score);
  }

  // Sort the results.
  qsort(rsr, num_motifs, sizeof(AME_RESULT_T*), ame_compare_log_evalues);

  // Output the motif results
  if (!options->text_only) {
    jsonwr_property(options->json_output, "motifs");
    jsonwr_start_array_value(options->json_output);
  }
  for (i = 0; i < num_motifs; ++i) {
    AME_RESULT_T *result = rsr[i];
    MOTIF_DB_T *db = arraylst_get(result->db_idx, dbs);
    if (result->is_significant) {
      if (!options->text_only) output_json_result(result);
      output_tsv_result(result, i+1, db->source);
    }
  }
  if (!options->text_only) {
    jsonwr_end_array_value(options->json_output);
  }

  // Finish the TSV output
  char *version_message = "# AME (Analysis of Motif Enrichment): Version " VERSION " compiled on " __DATE__ " at " __TIME__ "\n";
  //  "# Copyright (c) Robert McLeay <r.mcleay@imb.uq.edu.au> & Timothy L. Bailey <tlawbailey@gmail>, 2009-2017.\n";
  tsvwr("\n%s", version_message);
  tsvwr("# The format of this file is described at %s/%s.\n", SITE_URL, "doc/ame-output-format.html");
  tsvwr("# %s\n", options->commandline);

} // ame_finish_output

/*************************************************************************
 * Compute the enrichment of each motif.
 *************************************************************************/
AME_RESULT_T **compute_motif_enrichment(
  ARRAYLST_T *dbs,		// the motif DBs
  ARRAY_T *background,		// the motif background
  int num_seqs, 		// the number of sequences
  double avg_length,		// average length of sequences
  SEQ_T **sequences,		// the sequences
  bool fasta_scores, 		// FASTA scores given? 
  ARRAY_T *seq_fscores,		// FASTA scores for sequences
  int *n_motifs, 		// OUT total number of motifs in all DBs
  double *min_sequence_score,	// OUT minimum score over all motifs and sequences
  double *max_sequence_score,	// OUT maximum score ovar all motifs and sequences
  double *max_median_sequence_score	// OUT maximum over motifs of the median score
) {
  int i, db_i, motif_i;
  MOTIF_DB_T *db = NULL;
  MOTIF_T *motif = NULL, *rc_motif = NULL;
  PSSM_T *pssm = NULL, *rc_pssm = NULL;
  //bool use_log_odds = (TOTAL_HITS == options.scoring);
  bool use_log_odds = false;
  bool use_rc = alph_has_complement(options.alphabet);
  *min_sequence_score = +1e300;
  *max_sequence_score = -1e300;
  *max_median_sequence_score = -1e300;
  // Get the total number of motifs in all databases.
  int total_motifs = 0;		// total number of motifs in all DBs
  for (db_i = 0; db_i < arraylst_size(options.motif_sources); db_i++) {
    db = arraylst_get(db_i, dbs);
    total_motifs += arraylst_size(db->motifs);	// number of motifs in this DB
  }
  // Create the array to hold results per motif.
  AME_RESULT_T **rsr = (AME_RESULT_T **)mm_malloc(sizeof(AME_RESULT_T *) * total_motifs);
  int id = 0;			// overall index of motif
  // Compute enrichment of each motif.
  for (db_i = 0, i = 1; db_i < arraylst_size(options.motif_sources); db_i++) {
    db = arraylst_get(db_i, dbs);
    int num_motifs = arraylst_size(db->motifs);	// number of motifs in this DB
    for (motif_i = 0; motif_i < num_motifs; motif_i++, id++) {
      motif = (MOTIF_T *) arraylst_get(motif_i, db->motifs);
      // Create the PSSM.
      pssm = ame_make_pssm(use_log_odds, background, motif);
      // Set the minimum score threshold.
      double min_score_thresh = 1;
      int n_sites = (avg_length - pssm->w + 1) * (use_rc ? 2 : 1);
      // A strong site has odds score corresponding to hit_lo_fraction * the maximum log-odds score possible.
      double epsilon = 1e-15;		// because of roundoff
      double strong_site_score = (1-epsilon)*exp(options.hit_lo_fraction * log(pssm->nolog_max_score));
      if (MAX_ODDS == options.scoring) {
	// sequence has a strong site
	min_score_thresh = strong_site_score;
      } else if (SUM_ODDS == options.scoring) {
	// sequence has a strong site and rest of sequence averages odds=1 per site
	min_score_thresh = strong_site_score + n_sites - 1;
      } else if (AVG_ODDS == options.scoring) {
	// sequence has a strong site and rest of sequence averages odds=1 per site
	min_score_thresh = (strong_site_score + n_sites - 1) / n_sites;
      } else if (TOTAL_HITS == options.scoring) {
	// sequence has at least one site
	min_score_thresh = strong_site_score;
      }
      // If required, do the same for the reverse complement motif.
      if (use_rc) {
        rc_motif = dup_rc_motif(motif);
        rc_pssm = ame_make_pssm(use_log_odds, background, rc_motif);
      }
      // Scan with the motif to determine sequence scores.
      double *scores = ame_scan_sequences(
        min_score_thresh,// the hit threshold
	num_seqs,	// the number of sequences
	sequences,	// the sequences
	background,	// background model
        id,		// overall index of motif
	motif,		// the motif
	pssm,		// the PSSM
	rc_motif,	// the reverse complement motif (or NULL)
	rc_pssm 	// the reverse complement PSSM (or NULL)
      );
      // Get the statistical score for this motif.
      rsr[id] = ame_get_motif_score(
	motif,		// the motif
        id,		// overall index of motif
	db_i,		// index of motif DB
        db->source,	// filename of the motif DB
	num_seqs,	// number of sequences
        sequences,	// the sequences
	scores,		// sequence PWM scores for this motif
	options.scoring==TOTAL_HITS ? 1 : min_score_thresh, // minumum pwm threshold for motif
        total_motifs,	// the number of motifs (for E-value)
	fasta_scores,	// FASTA scores given?
	seq_fscores	// FASTA scores for sequences
      );
      // Update the minimum and maximum sequence scores over all sequences and motifs
      // and compute the median for this motif and save the largest median.
      qsort(scores, num_seqs, sizeof(double), ame_compare_doubles);
      if (scores[num_seqs-1] < *min_sequence_score) *min_sequence_score = scores[num_seqs-1];
      if (scores[0] > *max_sequence_score) *max_sequence_score = scores[0];
      double median = scores[num_seqs/2];
      if (median > *max_median_sequence_score) *max_median_sequence_score = median;
      //fprintf(stdout, "min %f max %f median %f max median %f\n", *min_sequence_score,*max_sequence_score,  median, *max_median_sequence_score);
      // Free memory associated with this motif.
      free_pssm(pssm);
      if (rc_pssm) free_pssm(rc_pssm);
      if (rc_motif) destroy_motif(rc_motif);
      myfree(scores);
    } // motif
  } // motif DB

  DEBUG_MSG(NORMAL_VERBOSE, "\n");

  assert(total_motifs == id);
  *n_motifs = total_motifs;
  return(rsr);
} // compute_motif_enrichment

/*************************************************************************
 * Entry point for AME
 *************************************************************************/
int main(int argc, char *argv[]) {
  int i;
  ARRAY_T *background = NULL;		// the motif background
  ARRAYLST_T *dbs = NULL;		// the motif DBs

  ame_set_defaults(); 		// Set default cmd line options.
  ame_getopt(argc, argv); 	// Get command line options.

  // set random number generators
  srand_mt(options.seed);
  set_randfunc((randfunc_t) random_mt); // for ushuffle

  /* Record the execution start and end times */
  if (verbosity >= HIGH_VERBOSE) {
    t0 = time(NULL);
    c0 = clock();
  }

  DEBUG_FMT(NORMAL_VERBOSE, "E-value threshold for reporting results: %g\n", options.evalue_report_threshold);

  // Set the alphabet "pseudo-option".
  bool xalph = (options.alph_file != NULL);
  options.alphabet = NULL;
  if (xalph) {
    options.alphabet = alph_load(options.alph_file, true);
    if (options.alphabet == NULL) exit(EXIT_FAILURE);
  }

  // Read all the alphabets and make sure they are the same.
  DEBUG_FMT(NORMAL_VERBOSE, "Checking alphabets in %d motif files.\n", arraylst_size(options.motif_sources));
  read_motif_alphabets(options.motif_sources, xalph, &(options.alphabet));

  if (!alph_has_complement(options.alphabet)) {
    if (options.scan_separately)
      ame_usage("The option --sep cannot be used with an unstranded alphabet.");
    options.scan_both_strands = false;
  }

  // Load the motifs and the background model.
  int max_width;
  dbs = ame_load_motifs_and_background(
    true,		// keep motif DB namespaces separate
    xalph,		// convert motifs to the user-specified alphabet
    &background,	// OUT background model
    &max_width		// OUT maximum motif width
  );

  // Load the primary sequences.
  SEQ_T **sequences = NULL;
  char *type = options.control_filename ? "primary " : "single set of ";
  DEBUG_FMT(NORMAL_VERBOSE, "Loading %ssequences.\n", type);
  int num_seqs = 0, num_skipped = 0;
  double avg_length = 0;	// average sequence length
  bool get_fasta_scores = (options.control_filename == NULL);
  bool fasta_scores = false;	// true if sequence ID lines have scores in them
  ARRAY_T *seq_fscores;		// the FASTA scores (or input order ranks) of the (primary) sequences
  ame_read_sequences(max_width, options.alphabet, options.seq_source, type,
    &num_seqs, &num_skipped, &avg_length, &sequences, get_fasta_scores, &fasta_scores, &seq_fscores);
  int pos_num_seqs = num_seqs;

  // Append the control sequences.
  int neg_num_seqs = 0;
  int neg_num_skipped = 0;
  if (options.control_filename) {
    options.fix_partition = num_seqs;
    if (strcmp(options.control_filename, "--shuffle--") == 0) {
      DEBUG_FMT(NORMAL_VERBOSE, "Creating control sequences by shuffling input sequences preserving %i-mers.\n", options.kmer);
      int i, j;
      int n_copies = (1000/num_seqs)+1;
      sequences = (SEQ_T**)mm_realloc(sequences, sizeof(SEQ_T*) * (n_copies+1) * num_seqs);
      for (i=1; i<=n_copies; i++) {
        int start = i * num_seqs;
	for (j=0; j<num_seqs; j++) {
	  sequences[start+j] = shuffle_seq(sequences[j], options.kmer, i);
	}
      }
      neg_num_seqs = num_seqs * n_copies;
      num_seqs += neg_num_seqs;
    } else {
      DEBUG_MSG(NORMAL_VERBOSE, "Loading control sequences.\n");
      get_fasta_scores = false;		// Never read FASTA scores from control sequences.
      ame_read_sequences(max_width, options.alphabet, options.control_filename, "control ", 
	&num_seqs, &neg_num_skipped, &avg_length, &sequences, get_fasta_scores, &fasta_scores, &seq_fscores);
      neg_num_seqs = num_seqs - options.fix_partition;
    }
  } // control

  if (num_seqs < MIN_SEQS) {
    die("You must provide a total of at least %d (valid) primary and control sequences.\n"
        "(%d primary and %d control sequences were removed.)\n", MIN_SEQS, num_skipped, neg_num_skipped);
  }

  // Check that fix_partition value is legal.
  if (options.fix_partition > num_seqs) {
    die("'--fix-partition' value cannot be larger than the total number of input sequences.\n");
  }
  if (options.control_filename) {
    DEBUG_FMT(NORMAL_VERBOSE, "Not in partition maximization mode. "
      "Fixing partition at the number of primary sequences (%d).\n", options.fix_partition);
  } else if (options.fix_partition > 0) {
    DEBUG_FMT(NORMAL_VERBOSE, "Not in partition maximization mode. "
      "Fixing partition at %i.\n", options.fix_partition);
  } else {
    DEBUG_MSG(NORMAL_VERBOSE, "In partition maximization mode.\n");
  }

  // Prepare factorial tables if necessary.
  if (options.pvalue_method == MULTIHG_METHOD || options.pvalue_method == LONG_MULTIHG_METHOD) {
    fisher_exact_init(num_seqs); // Generate the table of log factorials
    factorial_array = fisher_exact_get_log_factorials();
  }

  // Start the output of HTML.
  if (! options.text_only) {
    ame_start_output(
      argc, 
      argv, 
      background,
      sequences, 
      pos_num_seqs, 
      num_skipped, 
      neg_num_seqs, 
      neg_num_skipped, 
      fasta_scores,
      dbs
    );
  }

  // Compute enrichment of each motif.
  int total_motifs = 1;
  double min_sequence_score = 0;
  double max_sequence_score = 0;
  double max_median_sequence_score = 0;
  AME_RESULT_T **rsr = compute_motif_enrichment(
    dbs,
    background,		// the motif background
    num_seqs, 		// the number of sequences
    avg_length,		// average length of sequences
    sequences,		// the sequences
    fasta_scores, 	// FASTA scores given? 
    seq_fscores,	// FASTA scores for sequences
    &total_motifs, 	// OUT total number of motifs in all DBs
    &min_sequence_score, // minimum score over all motifs and sequences
    &max_sequence_score, // maximum score ovar all motifs and sequences
    &max_median_sequence_score // maximum over motifs of the median score
  );

  // Output results to HTML and TSV.
  ame_finish_output(
    &options,
    dbs,		// the motif DBs
    total_motifs,	// number of motifs with results
    rsr,		// results for each motif
    num_seqs,		// the number of sequences
    min_sequence_score,
    max_sequence_score,
    max_median_sequence_score
  );

  // Free the motif databases.
  arraylst_destroy(destroy_motif_db, dbs);

  for (i=0; i<num_seqs; i++) free_seq(sequences[i]);

  cleanup_options();

  ame_terminate(&options, 0); 	//Successfully end.

  return(0); 		//shuts up a warning.
} // main
