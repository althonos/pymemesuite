#ifndef STREME_UTILS_H
#define STREME_UTILS_H

#define _GNU_SOURCE // for asprintf
#include <stdio.h>
#include <stdarg.h> // for va_start, va_end
#include <getopt.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <libgen.h> // for basename
#include "config.h" // for VERSION
#include "alphabet.h"
#include "banner.h"
#include "utils.h"
#include "bates.h"
#include "binomial.h"
#include "fisher_exact.h"
#include "string-list.h"
#include "io.h"
#include "heap.h"
#include "ushuffle.h"
#include "st_types.h"
#include "st_symboldef.h"
#include "st_streetyp.h"
#include "st_spacedef.h"
#include "st_streeacc.h"
#include "st_streemac.h"
#include "st_megabytes.h"
#include "st_streehuge.h"
#include "st_arraydef.h"
#include "st_multidef.h"
#include "inputmultiseq.h"
#include "macros.h"
#include "matrix.h"
#include "motif.h"	// for motif2consensus
#include "motif-db.h"
#include "xml-out.h"
#include "projrel.h"	// for ARCHIVE_DATE

// Imported from st_construct.c
Sint constructstree(Suffixtree *stree, SYMBOL *text, Uint textlen);
void freestree(Suffixtree *stree);
// Imported from st_mapfile.c
void mmcheckspaceleak(void);
#define PRINTTIME(text) \
{ \
  DEBUG_FMT(NORMAL_VERBOSE, "# %s\t%.2f seconds \t(cumulative %.2f seconds)\n", (text), mytime(1)/1E6, mytime(0)/1E6); \
}

// Maximum allowed width of motifs.
#define MAX_WIDTH 30

// Minimum allowed width of motifs (do not make less than 3 or bad things happen!)
#define MIN_WIDTH 3

// Minimum width of seeds.
#define MIN_SEED_WIDTH 3

// Minimum number of input sequences.
#define MIN_SEQUENCES 2

// Prior on counts when initializing model.
#define SBM_SEED_PRIOR 0.01
#define REF_COUNT_PRIOR 0.1
#define NESTED_COUNT_PRIOR 0.1 

// Minimum size of hold-out set.
#define MIN_HO_SIZE 5

// Prior for creating higher-order background model.
#define MEME_BG_PRIOR 1

// Maximum allowed order.
#define MAX_ORDER 5

// Size for switching from list to bit table.
#define LISTSIZE 200

// Size of memory resize chunk for saving all sites.
#define RCHUNK 10000

// Size of p-value cache in doubles.
#define MAX_CACHE_SIZE 10000000

#define FLOAT_EPS 1e-10

// Print the current memory usage.
// Activate with: make CFLAGS=" -O3 -DMEMORY_MONITOR" <program>
#ifdef MEMORY_MONITOR
  #define PRINT_MEMORY_USAGE() \
   if (verbosity >= NORMAL_VERBOSE) { \
     char line[1000]; \
     snprintf(line, 1000, "ps -o pid,vsz,rss,%%mem | awk '$1 == %d {print \"## MEMORY_USAGE: \", \"vsz:\", $2, \"rss:\", $3, \"%%mem:\", $4; exit(0)};'", getpid()); \
     system(line); \
   }
#else
  #define PRINT_MEMORY_USAGE()
#endif

// Objective function type.
typedef enum {DE, CD, NO_OBJFUN} OBJFUN_T;

// Saving matches type.
typedef enum {NONE, ZOOPS, PASSING} SAVE_MATCHES_T;

// Site positional distribution alignment type.
typedef enum {LEFT, CENTER, RIGHT} SEQ_ALIGN_T;

// Threshold type.
typedef enum {EVALUE, QVALUE, PVALUE} THRESH_T;

DEXTERN(OBJFUN_T, DEFAULT_OBJFUN, DE);
DEXTERN(int, DEFAULT_ALPHABET_TYPE, Dna);
DEXTERN(int, DEFAULT_DNA_ORDER, 2);
DEXTERN(int, DEFAULT_RNA_ORDER, 2);
DEXTERN(int, DEFAULT_PROT_ORDER, 0);
DEXTERN(int, DEFAULT_CUSTOM_ORDER, 0);
DEXTERN(int, DEFAULT_MINWIDTH, 8);
DEXTERN(int, DEFAULT_MAXWIDTH, 15);
DEXTERN(int, DEFAULT_NEVAL, 25);
DEXTERN(int, DEFAULT_NREF, 4);
DEXTERN(int, DEFAULT_NITER, 20);
DEXTERN(int, DEFAULT_NMOTIFS, 5);		// use pvt by default; kicks in if pvt == 0
DEXTERN(int, DEFAULT_TIME, 0);			// must be >= 0; 0->no limit
DEXTERN(double, DEFAULT_EVT, 0.05);
DEXTERN(double, DEFAULT_PVT, 0.05);
DEXTERN(int, DEFAULT_PATIENCE, 3);		// quit after this many consecutive non-significant motifs 
DEXTERN(double, DEFAULT_MIN_PAL_RATIO, 0.85);
DEXTERN(double, DEFAULT_MAX_PAL_ED, 5);
DEXTERN(double, DEFAULT_HOFRACT, 0.10);
DEXTERN(double, DEFAULT_MINSCORE, 0);		// minimum score during approximate matching
DEXTERN(double, DEFAULT_OPTSCORE, 0);		// minimum score during motif p-value optimization
DEXTERN(int, DEFAULT_NSUBSETS, 1);
DEXTERN(int, DEFAULT_IGNORE_DEPTH, 5);		// Don't enforce score cutoff below this depth
DEXTERN(int, DEFAULT_SEED, 0);
DEXTERN(double, DEFAULT_PSEUDOCOUNT, 0);
DEXTERN(double, DEFAULT_LOG_PVALUE, 0);
DEXTERN(int, DEFAULT_CCUT, 0);
DEXTERN(int, DEFAULT_VERBOSITY, NORMAL_VERBOSE);
DEXTERN(VERBOSE_T, verbosity, NORMAL_VERBOSE);	// for MEME routines

// The alphabets.
DEXTERN(char*, DNA, "ACGT");
DEXTERN(char*, RNA, "ACGU");
DEXTERN(char*, PROTEIN, "ACDEFGHIKLMNPQRSTVWY");

// Functions for converting MEME-style alphabet to st-style.
// Convert alphabet letter or its complement to an index in [0,alen) 
// and convert a letter to its complement, or an index to its complement.
#define MAX_ALENGTH 256
EXTERN Uchar ALPH[MAX_ALENGTH+1];
EXTERN Uchar COMPALPH[MAX_ALENGTH+1];
EXTERN int ALPHINDEX[MAX_ALENGTH+1];
EXTERN int COMPINDEX[MAX_ALENGTH+1];
#define I2A(i) ALPH[(i)]
#define A2I(a) ALPHINDEX[(Uint) (a)]
#define C2I(a) COMPINDEX[(Uint) (a)]
#define COMP(a) COMPALPH[(Uint) (a)]
#define I2CI(i) C2I(I2A(i))

//
// Convert a model to a palindrome.
//
#define PALINDROMIZE(model) \
{ \
  int i, j; \
  int w = (model)->width; \
  int alen = (model)->alen; \
  for (j=0; j<=w/2; j++) { \
    for (i=0; i<=alen; i++) { \
      (model)->probs[i][j] = \
        (model)->probs[I2CI(i)][w-j-1] = \
           ((model)->probs[i][j] + (model)->probs[I2CI(i)][w-j-1])/2.0; \
    } \
  } \
}

//
// Determine whether the given sequence is a positive or negative,
// and whether it is the input sequence or the reverse complement of the input sequence.
// Returns the number of the input sequence corresponding to seqno (input_seqno) and
// the sequence number of the reverse complement of the given sequence (rc_seqno),
// which may be either the input strand or the (added) rc strand.
#define GET_SEQUENCE_NUMBER_AND_TYPE(multiseq, seqno, on_revcomp, is_negative, input_seqno, rc_seqno) \
{ \
  Uint npos = (multiseq)->npos; \
  Uint nneg = (multiseq)->nneg; \
  Uint ntot = npos + nneg; \
  on_revcomp = ((multiseq)->do_rc && (seqno) >= ntot); \
  /* Get the input sequence number and the reverse-complement sequence number */ \
  if (on_revcomp) { \
    rc_seqno = input_seqno = seqno - ntot; \
  } else { \
    rc_seqno = seqno + ntot; \
  } \
  is_negative = ((input_seqno) >= npos); \
} 

// Determine if a sequence number is the reverse complement
// of an input sequence.
#define IS_SEQUENCE_ON_REVCOMP(multiseq, seqno) \
  ((multiseq)->do_rc && (seqno) >= ((multiseq)->npos + (multiseq)->nneg)) 

// Initialize a model PROB and PSSM matrices to a consensus sequence.
#define INIT_PSSM_TO_CONSENSUS(model, consensus, prior, background, order) \
  { \
    int i, j; \
    for (j=0; j<(model)->width; j++) { \
      for (i=0; i<(model)->alen; i++) { \
	  (model)->probs[i][j] = ((consensus[j]) == I2A(i)) ? ((1+(prior)*(background)[i]))/(1+(prior)) : (prior)*(background[i])/(1+(prior)); \
      } \
    } \
    INIT_PSSM_FROM_PROBS((model), (background), (order)); \
    (model)->seed = (consensus); \
  }

// Create a PSSM matrix from the PROB and BACKGROUND matrices.
// If the background order > 0, the odds denominators are all set to 1.0.
#define INIT_PSSM_FROM_PROBS(model, background, order) \
  { \
    int i, j; \
    (model)->min_possible_score = 0; \
    for (j=0; j<w; j++) { \
      double min_possible_score_j = 0; \
      for (i=0; i<(model)->alen; i++) { \
	(model)->pssm[i][j] = ( (order) == 0 ? log((model)->probs[i][j] / (background)[i]) : log((model)->probs[i][j]) ) / log(2); \
	if ((model)->pssm[i][j] < min_possible_score_j) min_possible_score_j = (model)->pssm[i][j]; \
      } \
      (model)->min_possible_score += min_possible_score_j; \
    } \
  }

// Print the MODEL matrices.
#define PRINT_MODEL(model, text, type) \
  { \
    int i, j; \
    fprintf(stderr, "%s MODEL:\n", (text)); \
    for (i=0; i<(model)->alen; i++) { \
      for (j=0; j<(model)->width; j++) { \
        if ((type) == 0) { \
	  fprintf(stderr, " %.3f", (model)->counts[i][j]); \
        } else if ((type) == 1) { \
	  fprintf(stderr, " %.3f", (model)->probs[i][j]); \
        } else if ((type) == 2) { \
	  fprintf(stderr, " %.3f", (model)->pssm[i][j]); \
        } \
      } \
      fprintf(stderr, "\n"); \
    } \
  }

// Print the MOTIF in MEME format.
//  1) If have_holdout, prints the hold-out p-value in the "P=" field (the "E=" field in MEME output),
//  2) else, prints the training set p-value in the "S=" field (the "E=" field in MEME output).
// The "nsites=" field is set to the SUM of the positive test and hold-out matches.
// Uses the single-letter consensus for the motif name, prepending "text".
#define PRINT_MOTIF(model, motifno, index, text, stream, have_holdout) \
  { \
    int _i, _j; \
    double m1, e1, prec=1; \
    char *pv_str; \
    if ((have_holdout)) { \
      if ((index) == -1) { \
        exp10_logx((model)->test_log_evalue/log(10), m1, e1, prec); \
        pv_str = "E="; \
      } else { /* cand */ \
	exp10_logx((model)->test_log_pvalue/log(10), m1, e1, prec); \
        pv_str = "P="; \
      } \
    } else { \
      exp10_logx((model)->train_log_pvalue/log(10), m1, e1, prec); \
      pv_str = "S="; \
    } \
    if (*(text) != '\0') { \
      fprintf((stream), "MOTIF %s%d-%g-%s\n", (text), (motifno), (double) (index), (model)->consensus); \
    } else { \
      fprintf((stream), "MOTIF %d-%s STREME-%d\n", (motifno), (model)->consensus, (motifno)); \
    } \
    fprintf((stream), "letter-probability matrix: alength= %d w= %d nsites= %d %s %3.1fe%+04.0f\n", \
      (model)->alen, (model)->width, (model)->train_pos_count + (model)->test_pos_count, pv_str, m1, e1); \
    for (_j=0; _j<(model)->width; _j++) { \
      for (_i=0; _i<(model)->alen; _i++) { \
        fprintf((stream), " %.6f", (model)->probs[_i][_j]); \
      } \
      fprintf((stream), "\n"); \
    } \
    fprintf((stream), "\n"); \
    fflush(stream); \
  }

//
// Print runtime statistics and information.
//
#define MAX_HOSTNAME 255
#define PRINT_STATS(stream, commandline) \
  { \
    char hostname[MAX_HOSTNAME]; \
    int result = gethostname(hostname, MAX_HOSTNAME); \
    PSTARS(stream); \
    fprintf((stream), "COMMAND:\t%s\n", (commandline)); \
    PSTARS(stream); \
    fprintf((stream), "CPU:\t\t%s\n", (result<0) ? "unknown" : hostname); \
    PSTARS(stream); \
    fprintf((stream), "FINALTIME:\t%.2f seconds\n", mytime(0)/1E6); \
    PSTARS(stream); \
  }

/*
 get the exponent and mantissa of large numbers expressed as log10(x)
 for printing with prec digits after the decimal point in mantissa
*/
#define exp10_logx(logx, m, e, prec) { \
(e) = floor(logx); \
(m) = pow(10.0, (logx) - (e)); \
  if (m+(.5*pow(10,-prec)) >= 10) { (m) = 1; (e) += 1;} \
}

// 
// Get the corresponding head in the reverse-complement of the sequence.
//
#define get_rc_head(rc_head, multiseq, seqno, pos, node_maxdepth) \
{ \
  int seqstart = (multiseq)->seqstarts[seqno]; \
  int seqlen = (multiseq)->seqlengths[seqno]; \
  int totallength = ((multiseq)->totallength/2); \
  BOOL on_revcomp = IS_SEQUENCE_ON_REVCOMP((multiseq), (seqno)); \
  int rc_seqstart = (on_revcomp) ? seqstart - totallength - 1 : seqstart + totallength + 1; \
  (rc_head) = (multiseq)->sequence + rc_seqstart + (seqlen - (pos) - (node_maxdepth)); \
}

//
// Compute the p-value of a set of sites and round it to RNDDIG places.
//
#define GET_PVALUE(log_pvalue, objfun, use_binomial, pos_count, neg_count, npos, nneg, bernoulli, dtc_sum, w, len, text, cache, cache_length) \
  { \
    if ((objfun) == DE) { \
      Uint tot_count = (pos_count) + (neg_count); \
      /* Initialize the sub_cache for this total count. */ \
      if (tot_count < (cache_length) && (cache)[tot_count] == NULL) { \
        (cache)[tot_count] = (double*) malloc((tot_count+1)*sizeof(double)); \
        Uint i; \
        for (i=0; i<=tot_count; i++) (cache)[tot_count][i] = 2; \
      } \
      /* Get the hashed p-value if possible */ \
      if (tot_count < (cache_length) && (cache)[tot_count][(pos_count)] <= 1) { \
        (log_pvalue) = (cache)[tot_count][(pos_count)]; \
      } else { \
	(log_pvalue) = \
          (use_binomial) ? \
	    /* DE: Binomial Test */ \
	    ((pos_count) == 0 ? 0 : log_betai((pos_count), (neg_count) + 1, (bernoulli))) : \
	    /* DE: Fisher Exact Test */ \
	    getLogFETPvalue((pos_count), (npos), (neg_count), (nneg), False); \
      } \
      if (tot_count < (cache_length) && (cache)[tot_count][(pos_count)] > 1) { \
        (cache)[tot_count][(pos_count)] = (log_pvalue); \
      } \
    } else if ((objfun) == CD) { \
      (log_pvalue) = \
        /* CD: Cumulative Bates distribution */ \
	(pos_count > 0) ? get_dtc_log_pv((dtc_sum)/(pos_count), (pos_count), (w), (len)) : 0; \
    } else { \
      fprintf(stderr, "ERROR: unknown objfun in %s\n", text); \
      exit(EXIT_FAILURE); \
    } \
    RND((log_pvalue), 8, (log_pvalue)); \
    (log_pvalue) = MIN(0, (log_pvalue)); \
  }

//
// Set the score of a node for the given weight from the srch_state data.
//
#define SETNODESCORE(srch_state, nodestdptr, w, score) \
  { \
    if ((nodestdptr)->scores == NULL) { \
      int minw = (nodestdptr)->min_width; \
      int maxw = MIN((nodestdptr)->max_width, (srch_state)->maxwidth); \
      (nodestdptr)->scores = (double *) malloc((maxw - minw + 1) * sizeof(double)); \
    } \
    (nodestdptr)->scores[w-minw] = (score); \
  }

//
// Compute the log conditional background probability (lcbp) 
// of the current letter given the Markov background model 
// and subtract it from llr.
//
#define SUBLCBP(word, w, alen, lcbp, bg_order, llr) \
  { \
    int _i, index = 0, a2n = 1, offset = 0; \
    int w0 = (w)>(bg_order) ? (bg_order+1) : (w); \
    for (_i=0; _i<w0; _i++) { \
      if (_i > 0) { \
        a2n *= (alen); \
        offset += a2n; \
        index *= (alen); \
      } \
      index += A2I((word)[(w)-w0+_i]); \
    } \
    index += offset; \
    (llr) -= (lcbp)[index]; \
  }

// Motif Site
typedef struct {
  Uint seqno;		// input sequence number of site
  Uint pos;		// position of site in input sequence
  Uchar strand;		// '-': site is on negative (not input) strand
  BOOL is_positive;	// true if sequence is positive
  double score;		// score of site
  Uint pos_count;	// number of positive sequences with this site
  Uint neg_count;	// number of negative sequences with this site
  Uint pos_zoops;	// number of positive sequences where this is the best site
  Uint neg_zoops;	// number of negative sequences where this is the best site
  double log_pvalue;	// log p-value of site based on all occurrences (not just ZOOPS)
  int dtc;		// average distance between site centers and sequence centers
} Site;

// Sequence with a passing site.
typedef struct {
  int dbno;		// sequence database: 0,1
  int seqno;		// sequence number within the database
  Uchar *desc;		// Sequence ID
  double score;		// Best score of sequence
  BOOL is_tp;		// True if primary sequence, False if control sequence.
} Passing_seq;

// Motif MODEL.
// MAX_CUST_ALEN = 10 digits + 52 letters + 4 other characters (Resize(name, i+RCHUNK, char) = 66
#define MAX_CUST_ALEN 66
typedef struct {
  ALPH_T *alph;					// the model's alphabet
  Uint alen;					// rows
  Uint width;					// columns
  Uint initial_width;				// initial width of seed node
  Uchar *seed;					// seed string
  double probs[MAX_CUST_ALEN][MAX_WIDTH+2];	// PSPM
  double pssm[MAX_CUST_ALEN][MAX_WIDTH+2];	// log-odds PSSM
  double counts[MAX_CUST_ALEN][MAX_WIDTH+2];	// counts
  double min_possible_score;			// the minimum possible score of the PSSM
  // training set
  Uint train_pos_count;				// best number of positives
  Uint train_neg_count;				// best number of negatives
  double train_dtc;				// mean distance from site center to sequence center
  double score_threshold;			// best threshold for separating pos and neg sites
  double train_log_pvalue;			// best log p-value on training set
  double train_ratio;				// ratio of positive to negative counts on training set
  int n_eff_tests;				// number of successful multiple tests (if no hold-out set)
  MOTIF_T *motif;				// the MEME-formatted motif entry
  int db_idx;					// index of DB containing motif
  // hold-out set
  Uint test_pos_count;				// best number of hold-out positives
  Uint test_neg_count;				// best number of hold-out negatives 
  double test_dtc;				// best hold-out mean distance from site center to sequence center
  double test_log_pvalue;			// best log p-value on hold-out set
  double test_log_evalue; 			// best log E-value on hold-out set
  double test_log_qvalue;	 		// log q-value on hold-out set (used by SEA only)
  double test_ratio;				// ratio of positive to negative counts on hold-out set
  // sites
  Site *matches;				// positive and negative matches to erase (malloced per model)
  Uint nmatches;				// length of matches list
  int total_sites;				// total number positive sequences with PASSING sites
  int *site_distr;				// counts of PASSING sites; sequences are aligned to center of longest sequence
  int max_sites;				// maximum number of PASSING sites in any positive sequence
  int *site_hist;				// number of PASSING positive sequences with each number of PASSING sites <= max_sites
  // matching sequences
  Passing_seq *passing_seqs;			// sequences with a passing site
  int npassing;					// number of sequences with a passing site
  // other
  BOOL is_palindromic;				// model is palindromic
  double ed;					// palindrome edit distance
  int iter;					// number of refinement iterations done
  double elapsed_time;				// elapsed time (seconds) up until now
  char *consensus;				// single-letter consensus of model (malloced per model)
} Model;

// To hold current depth and number of errors during MODEL approximate matching.
typedef struct {
  Uint depth;
  double errors;
  double score;
} Modelstate;
DECLAREARRAYSTRUCT(Modelstate);

// Data needed for adding the sequence counts to the suffix tree.
typedef struct {
  Suffixtree *stree;	// the suffix tree
  Multiseq *multiseq;	// the input target sequences
  ArrayBref mother;	// the stack of mother nodes
  ArrayBref ancestor;	// the stack of not-too-deep ancestor nodes
  Uint mindepth;	// minimum node depth to store sequence counts for
  Uint maxdepth;	// maximum node depth to store sequence counts for
  ArrayUint valid_node_indices;// indices of valid nodes in stree->subtreedata
  OBJFUN_T objfun;	// objective function
  double **pv_cache;	// cache for p-values, indexed by pos_count + neg_count
  Uint cache_length;	// length of pv_cache array
} Countstate;

// Data needed for searching suffix tree for approximate matches.
typedef struct {
  Suffixtree *stree;	// the suffix tree
  Reference *rootref;	// the root
  Multiseq *multiseq;	// the input training sequences
  Multiseq *test_multiseq;// the input hold-out sequences
  ArrayBref mother;	// the stack of mother nodes
  Uint mindepth;	// minimum node depth to get approximate matches for
  Uint maxdepth;	// maximum node depth to get approximate matches for
  Uint maxwidth;	// maximum width of motifs being searched for
  double mmr;		// allowed mismatches per motif position [0,..1)
  double minscore;	// allowed minimum score for score-based matching
  int ignore_depth;	// ignore minimum score up to this depth for score-based matching
  Model *model;		// model being matched
  ArrayModelstate modelstate; // stack for depth and errors/score on path prior to current node
  ArraySubtreedataptr matching_nodes[MAX_WIDTH+3];	// node pointers for matching words of different widths
} Searchstate;

// Result data from score-based matching, indexed by motif width.
typedef struct {
  double *log_pvalues;
  double *score_thresholds;
  int *pos_counts;
  int *neg_counts;
  double *dtc_sum;
} Sbmdata;

// Structure for holding an evaluated seed.
typedef struct {
  Subtreedata *seed;
  int width; 
  double mmr;
  double log_pvalue;
} Evaluatedseed;
DECLAREARRAYSTRUCT(Evaluatedseed);

// Program options.
DECLAREARRAYSTRUCT(double);
typedef struct options {
  char *output_dirname;			// directory for STREME output
  BOOL allow_clobber;			// allow overwrite of directory
  BOOL text_only;			// output text only; don't create directory
  char *posfile;			// positive sequences filename
  char *negfile;			// negative sequences filename
  int totallength;			// truncate each sequence set (+/-) to this length
  OBJFUN_T objfun;			// objective function
  int order;				// use this order background model and sequence shuffling
  char *bfile;                          // name of MEME-style background model file
  ALPHABET_T alphabet_type;		// alphabet type (Dna, Rna, Protein, Custom)
  ALPH_T *alph;				// alphabet for dataset 
  char *alph_file;			// custom alphabet file
  int minwidth;				// minimum motif width
  int maxwidth;				// maximum motif width
  int width;				// set min and max width to this
  int neval;				// number of initial seeds of each width to evaluate
  int nref;				// number of evaluated seeds to refine
  int niter;				// maximum number of refinement iterations/seed
  //double pvt;				// p-value threshold for motifs
  double thresh;			// significance threshold for reporting enriched motifs
  THRESH_T thresh_type;			// statistic for significance threshold (E-, q-, or p-value)
  int patience;				// quit after this many consecutive non-significant motifs
  int nmotifs;				// number of motifs to find
  double time;				// try to quit before <time> CPU seconds consumed
  double hofract;                       // fraction of sequences in hold-out set
  BOOL usepv;				// choose best model based on p-value
  BOOL cand;                            // print candidate motifs of each width
  double minscore;			// minimum score for score-based matching
  int ignore_depth;			// ignore minimum score up to this depth
  int nsubsets;				// number of nested subsets
  int seed;				// random seed for choosing hold-out set
  double min_pal_ratio;			// minimum log_pvalue ratio for choosing palindromic model
  double max_pal_ed;			// maximum Euclidean distance for choosing palindromic model
  char *description;			// value of --desc or the contents of file named dfile
  char *dfile;				// name of file containing description
  BOOL experimental;			// set true if compiled with EXP defined
  char *align_txt;			// align sequences left/center/right for distr plots
  SEQ_ALIGN_T align;			// align sequences left/center/right for distr plots
//
// Derived options:
  int kmer;				// use this for sequence shuffling
  FILE *text_output;
  char *commandline;
  double **pv_cache;			// cache for p-values, indexed by pos_count + neg_count
  int cache_length;			// length of pv_cache array
//
// Options for SEA only:
  FILE *sequences_tsv_output;		// name of TSV file for the passing sequences of each significant motif
  BOOL noseqs;				// don't output passing sequences if true
  ARRAYLST_T *motif_sources;            // filenames of the motif library
  ARRAYLST_T *include_patterns;         // Wildcard patterns of motif names to include
  ARRAYLST_T *exclude_patterns;         // Wildcard patterns of motif names to exclude
  double pseudocount;                   // add to the motif frequency counts
// Options for fasta-holdout-set only:
  int ccut;				// centrally trim sequences longer than this
} STREME_OPTIONS_T;

void shuffle_multiseq(
  Multiseq *multiseq,                   // multiseq to shuffle
  int kmer,                             // preserve frequencies of k-mers
  Uchar fixed                           // the character to fix locations of
);

Multiseq *read_pos_neg_seqs(
  STREME_OPTIONS_T *options,
  BOOL do_rc,                           // append reverse complement of pos+neg sequences
                                        // negative sequences by shuffling positive sequences
  BOOL allow_ambigs,                    // don't convert ambiguous characters to SEPARATOR
  BOOL no_trim,                         // don't trim sequences to average length
  BOOL set_back,                        // set the background in the multiseq objects
  ARRAY_T *background,                  // IN MEME-style background model
  Multiseq **test_multiseq_ptr          // OUT hold-out set
);

void set_multiseq_background(
  STREME_OPTIONS_T *options,
  Multiseq *multiseq,
  ALPH_T *alph,
  BOOL do_rc,
  ARRAY_T *back                 // MEME-style background model (or NULL)
);

void set_logcumback(
  STREME_OPTIONS_T *options,
  Multiseq *multiseq
);

void append_to_multiseq(
  Multiseq *multiseq,           // old sequences
  Multiseq *multiseq2,          // sequences to append
  BOOL no_trim,
  BOOL is_ho
);

void append_rc_to_multiseq (Multiseq *multiseq);

void score_model_pssm(
  STREME_OPTIONS_T *options,            // STREME options
  Multiseq *multiseq,                   // the positive and negative sequences
  Model *model,                         // the model
  BOOL find_score_threshold,            // find the best score threshold and set in model
  BOOL set_pvalue,                      // compute the enrichment p-value of the model
  SAVE_MATCHES_T save_matches,          // ZOOPS: save best match (sorted) in every sequence (positive or negative) for nesting
                                        // PASSING: save all matches with score >= model->score_threshold for erasing
  BOOL is_ho,                           // data is the hold-out set
  BOOL use_cache                        // use the p-value cache
);

void get_passing_sequences(
  STREME_OPTIONS_T *options,    // STREME options
  Model *model,                 // the model, including its PASSING sites (see score_model_pssm)
  Multiseq *multiseq,           // the sequences
  BOOL append,                  // append passing sites to model if true
  BOOL sort                     // sort the passing sequences by score
);

void get_site_distr(
  Model *model,                 // The model, including its PASSING sites (see score_model_pssm)
  Multiseq *multiseq,           // The sequences
  SEQ_ALIGN_T align             // Align sequences left, center or right.
);

int compare_site_score(
  const void *v1,
  const void *v2
);

int compare_sequence_score(
  const void *v1,
  const void *v2
);

Uint initialize_st_alphabet(
  ALPH_T *alph                  // the MEME-style alphabet
);

char *get_single_letter_consensus(
  Model *model
);

//
// Read the contents of a description file.
//
char *get_description_file(char *filename);

// Imported from st_mapfile.c
void mmcheckspaceleak(void);

#endif
