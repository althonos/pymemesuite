/***********************************************************************
*                                                                      *
* MEME                                                                 *
* Copyright 1994-2017, The Regents of the University of California     *
* Author: Timothy L. Bailey                                            *
*                                                                      *
***********************************************************************/
/**********************************************************************/
/*
  Display routines for the results of MEME
*/
/**********************************************************************/
// 7-23-99 tlb; replace nsites() calls with model->nsites 
// 7-16-99 tlb; move SITE to meme.h, get_sites to meme_util.c 
// 7-14-99 tlb; change simplified prob. matrix to be observed frequencies,
//              consensus to be based on observed frequencies
// 4-7-99 tlb; fix bug in get_sites: setting of best site in oops and zoops 

#include <assert.h>
#include "display.h"
#include "meme-dtd.h"
#include "mast.h"
#include "projrel.h"
#include "motif.h"
#include "motif-spec.h"
#include "config.h"
#include "ceqlogo.h"
#include "dir.h"
#include "alphabet.h"
#include <string.h>
#include "string-builder.h"
#ifdef UNIX
#include <unistd.h>
#endif
#include "xml-out.h"

#ifndef PARALLEL
#define mpMyID() 0
#endif

static char *yesno[] = {"no", "yes"};
static char *stars = NULL;
static LO **los = NULL;		/* logodds structure for each motif */
static double **pv = NULL;	/* p-value tables for each motif */

// number of different integral score values for computing p-values 
#define PVAL_RANGE 100

// stepsize and size of smoothing window for get_q 
#define NSTEPS 100
// #define WINDOW (NSTEPS/10)
#define WINDOW 0

// round double to integer; round up if half way between 
// only works for positive numbers 
#ifndef NINT
  #define NINT(x) ((int) ((x)+0.5))
#endif

// maximum number of levels in consensus sequence 
#define MAXDEPTH ((int) (1.0/MINCONS))

// encode and decode a (seq #, offset) pair as one integer 
// this will bomb if seq. length > MAX_SEQ_LEN or
// if number of sequences * MAX_SEQ_LEN exceeds size of int
#define MAX_SEQ_LEN 10000
#define ENCODE(I,J) (I) * MAX_SEQ_LEN + (J);
#define DECODE(N, SEQ, OFF) {\
  (OFF) = (N) % MAX_SEQ_LEN;\
  (SEQ) = ((N) - (OFF)) / MAX_SEQ_LEN;\
}

// distance to indent start of RE histogram, consensus and simplified motif 
#define IND 13
#define IND2 6

// record describing accuracy of motif 
typedef struct {
  double thresh; // optimal threshold 
  int err; // classification errors using thresh 
  double roc; // ROC 
} ACCURACY;

// sortable sequence score record 
typedef struct {
  double score; // sequence score 
  int class; // class of sequence,  1=pos, 0=neg 
  char *id;
} SORTED_SCORE;

// sortable letter value record 
typedef struct {
  char letter;
  double value;
} LETTER_VALUE;

// local functions 
static void print_sites(
  DATASET *dataset, // the dataset 
  MODEL *model,     // the model 
  char *consensus,  // single-letter consensus of motif
  int format,       // 0=BLOCKS; 1=FASTA 
  char *com,        // comment to append 
  FILE *outfile     // where to print 
);
static void print_log_odds(
  int imotif,         // motif number 
  char *consensus,    // single-letter consensus
  DATASET *dataset,   // the dataset 
  int w,              // width of motif 
  double **logodds,   // log-odds matrix 
  double bayes,       // threshold 
  double log_ev,      // log E-value of motif 
  FILE *outfile       // where to print 
);
static void print_entropy(
  bool logo,       // true: prints the logo, false: otherwise 
  MODEL *model,       // the model 
  DATASET *dataset,   // the dataset 
  char *str_space,    // space for printing strand direction 
  FILE *outfile       // stream for output 
);
static void print_logo(
  MODEL *model,       // the model 
  DATASET *dataset    // the dataset 
);
static void print_candidates(
  CANDIDATE *candidates, // list of candidate models IN 
  DATASET *dataset,      // the dataset IN 
  int max_w,             // maximum width for motifs IN 
  FILE *outfile          // stream for output IN 
);
static double meme_score_sequence(
  uint8_t *eseq, // integer-coded sequence to score 
  int length, // length of the sequence 
  int w, // width of motif 
  double **logodds1, // log-odds matrix: LOG2(m_ij/b_j) 
  double **logodds2 // log-odds matrix: LOG2(m_ij/n_ij) 
);
static int s_compare(
  const void *v1,
  const void *v2
);
static int lv_compare(
  const void *v1,
  const void *v2
);
static LO *create_lo(
  ALPH_T *alph, // alphabet
  int imotif, // load number of motif 
  int w, // width of motif 
  double **logodds, // single-letter logodds matrix 
  double threshold // Bayes optimal threshold 
);
static void score_sites(
  DATASET *dataset, // the dataset 
  MODEL *model, // the model 
  LO *lo, // LO structure 
  double *pv // p-values for scores of this motif 
);
static void print_meme_header_xml(FILE *outfile);
static void print_meme_training_set_xml(
  DATASET *dataset, // the dataset 
  DATASET *neg_dataset, // the control dataset IN 
  FILE* outfile     // file for output 
);
static void print_meme_model_xml(
  MODEL *model,           // the model 
  DATASET *dataset,       // the dataset 
  char *stopping_reason,  // LO structure 
  FILE* outfile           // file for output 
);
static void print_meme_motifs_xml(
  MODEL *model,                   // the model IN 
  DATASET *dataset,               // the dataset 
  int nmotifs,                    // number of motifs 
  MOTIF_SUMMARY *motif_summaries, // List of final motif properties 
  FILE* outfile                   // output file 
);
static void print_meme_pssm_xml(
  double **logodds, // pointer to matrix of log-odds scores 
  ALPH_T *alph, // alphabet
  int width,       // width of the motif 
  FILE* outfile    // pointer to output file 
);
static void print_meme_psfm_xml(
  THETA theta,  // pointer to matrix of frequencies 
  ALPH_T *alph, // alphabet
  int width,    // width of the motif 
  FILE* outfile // pointer to output file 
);
static void print_meme_regular_expression_xml(
  char* regexp,     // regular expression  
  FILE* outfile     // pointer to output file 
);

static void print_meme_contributing_sites_xml(
  MODEL *model,				// the model 
  MOTIF_SUMMARY *motif_summary,         // summary of the motif 
  DATASET* dataset,                     // the dataset 
  FILE*  outfile                        // pointer to output file 
);

static void print_meme_scanned_sites_xml(
  MODEL *model,       // the model 
  DATASET *dataset,   // the dataset 
  int nmotifs,        // number of motifs 
  double **pv,        // p-value tables for each motif 
  FILE *outfile
);

static void print_site_diagrams(
  DATASET *dataset,   // the dataset 
  MODEL *model,       // the model 
  int nmotifs,        // number of motifs in los 
  char *consensus,    // single-letter consensus for motif
  FILE *outfile       // where to print 
);

void get_aligned_sequence_parts(
  ALPH_T *alph,
  int motif_start,
  int motif_width,
  bool ic,
  int lseq,
  char *seq,
  char *pre,
  char *site,
  char *post
);
static void align_sites(
  DATASET *dataset,     // the dataset 
  MODEL *model,         // the model 
  char *consensus,    // single-letter consensus
  LO *lo,               // LO structure 
  double *pv,           // pvalues for scores of this motif 
  FILE *outfile         // stream for output 
);
static void print_block_diagrams(
  MODEL *model,       // the model 
  DATASET *dataset,   // the dataset 
  int nmotifs,        // number of motifs 
  double **pv,        // p-value tables for each motif 
  FILE *outfile
);
static char *get_regexp(
  MODEL *model,       // the model 
  DATASET *dataset,   // the dataset 
  bool prosite        // true: [AK]-E-[EG]. false: [AK]E[EG] 
);

/**********************************************************************/
/*
	get_total_res

	Get the total size of the dataset ignoring noise sequences.
*/
/**********************************************************************/
static int get_total_res(
  DATASET *dataset    // the dataset 
) 
{
  int i;
  // Get size of dataset if NZ is objective function.
  int total_res = 0;
  if (dataset->objfun == NZ) {
    int n_samples = dataset->n_group[0];
    for(i = 0; i < n_samples; i++) {
      SAMPLE *s = dataset->input_order[i];
      total_res += s->length;
    }
  } else {
    total_res = dataset->total_res;
  }
  return(total_res);
} // get_total_res 

/**********************************************************************/
/*
  Record the results of EM
  Doesn't support negative datasets.
*/
/**********************************************************************/
void record_results(
  DATASET *dataset,              // the dataset IN 
  MODEL *model,                  // the model IN 
  MOTIF_SUMMARY *motif_summaries // summaries of final motifs IN 
)
{
  int i, j;
  int w = model->w;                     // width of last component 
  int nsites_dis = model->nsites_dis;   // # of sites after discretiz.*/
  double m1, e1, m2, e2;      // for printing significance 
  THETA theta = model->theta;
  THETA obs = model->obs;
  double lambda = model->lambda;
  double **logodds;
  ARRAY_T *back = dataset->back;
  int imotif = model->imotif-1;     // index of motif 
  double thresh;                    // Bayes optimal threshold 
  ALPH_T *alph = dataset->alph;

  // get p-value and E-value of motif 
  exp10_logx(model->logpv/log(10.0), m1, e1, 1);
  exp10_logx(model->logev/log(10.0), m2, e2, 1);

  // Record the results for the model as a whole 
  motif_summaries[imotif].width = w;
  motif_summaries[imotif].num_sites = nsites_dis;
  motif_summaries[imotif].ic = model->ic;
  motif_summaries[imotif].re = w * model->rel;
  motif_summaries[imotif].llr = model->llr;
  motif_summaries[imotif].p_value_mant = m1;
  motif_summaries[imotif].p_value_exp = e1;
  motif_summaries[imotif].e_value_mant = m2;
  motif_summaries[imotif].e_value_exp = e2;
  motif_summaries[imotif].sites = mm_malloc(model->nsites_dis * sizeof(p_prob));
  memcpy(motif_summaries[imotif].sites, model->maxima, model->nsites_dis * sizeof(p_prob));

  // make the log-odds matrices 
  logodds = make_log_odds(theta, NULL, back, 0, w, alph_size_core(alph));
  // calculate the optimal threshold (min classification error or Bayes' 
  thresh = LOG2((1-lambda)/lambda); // Bayes' threshold 
  motif_summaries[imotif].bayes = thresh;
  create_2array(motif_summaries[imotif].pssm, double, w + 1, alph_size_wild(alph));
  for (i = 0; i < w; i++) {
    for (j = 0; j < alph_size_core(alph); j++) {
      motif_summaries[imotif].pssm[i][j] = logodds(i, j);
    }
  }
  // destroy logodds matrix
 for (i = 0; i < w; i++) free(logodds[i]);
  free(logodds);

  create_2array(motif_summaries[imotif].psfm, double, w + 1, alph_size_wild(alph));
  for (i = 0; i < w; i++) {
    for (j = 0; j < alph_size_core(alph); j++) {
      motif_summaries[imotif].psfm[i][j] = obs(i, j);
    }
  }
  motif_summaries[imotif].regexp = get_regexp(model, dataset, false);
  motif_summaries[imotif].consensus = get_single_letter_consensus(model, dataset->alph);

  // Record elapsed time 
  motif_summaries[imotif].elapsed_time = mytime(0)/1E6;

} // Record results 

/**********************************************************************/
/*
  Print the results of EM
*/
/**********************************************************************/
void print_results(
  DATASET *dataset,      // the dataset IN 
  MODEL *model,          // the model 
  CANDIDATE *candidates, // candidate models found IN 
  FILE* outfile          // file for text output IN 
)
{
  char *regexp;
  int i;
  int max_w = model->max_w;               // maximum width tried 
  int nstrands = model->invcomp ? 2 : 1;  // # of strands to use 
  int w = model->w;                       // width of last component 
  int nsites_dis = model->nsites_dis;     // # of sites after discretize
  double m1, e1, m2, e2;                  // for printing significance 
  char *str_space = (nstrands == 1) ? "" : "       ";
  THETA theta = model->theta;
  THETA obs = model->obs;
  double lambda = model->lambda;
  double **logodds;
  char *cons = model->cons;
  ARRAY_T *back = dataset->back;
  double thresh;        // Bayes optimal threshold 
  int imotif = model->imotif-1;
  ARRAY_T *rounded_back = NULL;
  rounded_back = allocate_array(alph_size_core(dataset->alph));

  // Get a single-letter consensus for the motif
  char *consensus = get_single_letter_consensus(model, dataset->alph);

  // get entropy and E-value of motif 
  calc_entropy(model, dataset);

  // Retrieve the array of sites predicted in the model:
  P_PROB pred_sites = model->maxima;
  int n_psites = model->nsites_dis;

  // If requested, print the final MEME site predictions (for the "best"
  // starting point after EM has been completed):
  if (dataset->print_pred) {
    fprintf(stdout, "PREDICTED SITES AFTER MEME COMPLETION MOTIF %i-%s\n",
            model->imotif, consensus);
    print_site_array(pred_sites, n_psites, stdout, model->w, dataset);
    double sig = model->logev;
    fprintf(stdout, "FINAL MODEL SIGNIFICANCE = %f\n", sig);
  }

  // get p-value and E-value of motif 
  exp10_logx(model->logpv/log(10.0), m1, e1, 1);
  exp10_logx(model->logev/log(10.0), m2, e2, 1);

  // print the significant models 
  if (VERBOSE) {
    print_candidates(candidates, dataset, max_w, outfile);
  }

  // print the results for the model as a whole 
  fprintf(
    outfile,
    "\n\n%s\nMOTIF %s MEME-%d\twidth = %3d  sites = %3d  ",
    stars, consensus, model->imotif, w, nsites_dis
  );

  if (dataset->objfun == Classic || dataset->objfun == NC) {
    // E-value independent of number of starts
    fprintf(
      outfile,
      "llr = %.0f  E-value = %3.1fe%+04.0f\n%s\n",
      model->llr, m2, e2, stars
    );
  } else {			// E-value depends on number of starts
    fprintf(
      outfile,
      "llr = %.0f  p-value = %3.1fe%+04.0f  E-value = %3.1fe%+04.0f\n%s\n",
      model->llr, m1, e1, m2, e2, stars
    );
  }

  if (VERBOSE) {
    fprintf(
      outfile,
      "p-value = %3.1fe%+04.0f   E-value = %3.1fe%+04.0f\n%s\n",
      m1, e1, m2, e2, stars
    );
  }

  // print results for motif 
  if (VERBOSE) {
    char *cons0 = model->cons0;
    fprintf(outfile, "\n(best) %s --> %s\n", cons0, cons);
    fprintf(
      outfile,
      "(best) w %3d nsites %5.1f lambda %8.7f RE/col %6.3f\n",
      w, lambda*wps(dataset, w, 0), lambda, model->rel
    );
    fprintf(outfile, "\n(best) RE %6.3f\n\n", w * model->rel);
  }

  /*
    print motif description
  */
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fputc('\n', outfile);
  fprintf(outfile, "\tMotif %s MEME-%d Description\n", consensus, model->imotif);
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fputc('\n', outfile);

  // make the log-odds matrices 
  // calculate the optimal threshold (min classification error or Bayes' 
  logodds = make_log_odds(theta, NULL, back, 0, w, alph_size_core(dataset->alph));
  thresh = LOG2((1-lambda)/lambda); // Bayes' threshold 
  print_theta(
    0, NULL, 2, model->nsites_dis, obs, w, model->logev,
    str_space, dataset, outfile
  );
  print_entropy(true, model, dataset, str_space, outfile);
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fprintf(outfile, "\n\n");

  // create an LO structure and store it in the local array 
  Resize(los, imotif+1, LO*);
  los[imotif] = create_lo(dataset->alph, model->imotif, w, logodds, thresh);

  // round background so p-values are same as for MAST
  for (i = 0; i < alph_size_core(dataset->alph); i++) {
    double value, rounded_value;
    value = get_array_item(i, back);
    RND(value, 3, rounded_value);
    set_array_item(i, rounded_value, rounded_back);
  }

  // create a table of p-values and store it in the array 
  // pv[imotif] = calc_cdf(los[imotif], PVAL_RANGE, dataset->back);
  Resize(pv, imotif+1, double*);
  pv[imotif] = calc_pssm_cdf(
    los[imotif]->w, alph_size_core(los[imotif]->alph), PVAL_RANGE,
    los[imotif]->logodds, rounded_back
  );

  // score the sites and sort by position p-value 
  score_sites(dataset, model, los[imotif], pv[imotif]);

  // print alignment of the sites 
  if (dataset->n_samples <= dataset->brief) align_sites(dataset, model, consensus, los[imotif], pv[imotif], outfile);

  // print diagrams of the sites 
  if (dataset->n_samples <= dataset->brief) print_site_diagrams(dataset, model, model->imotif, consensus, outfile);

  // print the sites "making up" the model 
#ifdef PARALLEL
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: at print_print_sites\n", mpMyID()); fflush(stderr);
#endif
#endif
  if (dataset->n_samples <= dataset->brief) print_sites(dataset, model, consensus, PRINT_FASTA, "", outfile);
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: past print_print_sites\n", mpMyID()); fflush(stderr);
#endif

  // print the logodds matrix 
  print_log_odds(
    model->imotif, consensus, dataset, w, logodds, thresh, model->logev, outfile
  );
  // destroy logodds matrix
  for (i = 0; i < w; i++) free(logodds[i]);
  free(logodds);

  // print the observed frequency matrix 
  print_theta(
    model->imotif, consensus, 1, model->nsites_dis, obs, w, model->logev,
    str_space, dataset, outfile
  );

  // print a regular expression corresponding to motif
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fputc('\n', outfile);
  fprintf(outfile, "\tMotif %s MEME-%d regular expression\n", consensus, model->imotif);
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fputc('\n', outfile);
  regexp = get_regexp(model, dataset, false);
  fprintf(outfile, "%s\n", regexp);
  free(regexp);
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fprintf(outfile, "\n\n");

  // display elapsed time 
  fprintf(outfile, "\n\n\nTime %5.2f secs.\n\n", mytime(0)/1E6);

  // print line of stars 
  fprintf(outfile, "%s\n", stars);

  // flush 
  fflush(outfile);

  free_array(rounded_back);
  myfree(consensus);

} // print_results 

/**********************************************************************/
/*
  create_lo

  Create an an LO structure from the logodds matrix;
  include 'X' in the alphabet.
  PSSM is scaled by 100 rounded to the nearest integer,
  and then rescaled to PVAL_RANGE.
*/
/**********************************************************************/
static LO *create_lo(
  ALPH_T *alph,                 // alphabet
  int imotif,			// index of motif 
  int w,			// width of motif 
  double **logodds,		// single-letter logodds matrix 
  double threshold		// Bayes optimal threshold 
)
{
  int i, j, len, tmp;

  // create a logodds structure 
  LO *lo = NULL;      // the LO structure 
  lo = mm_malloc(sizeof(LO));
  memset(lo, 0, sizeof(LO));

  // initialize it 
  lo->alph = alph_hold(alph);
  lo->alen = alph_size_wild(alph); // add 'X' column
  lo->w = lo->ws = w;
  lo->imotif = imotif;
  lo->thresh = threshold;
  //convert the index into a string
  len = 1, tmp = imotif;
  while (tmp >= 1) { 
    tmp /= 10;
    ++len;
  }
  lo->meme_name = mm_malloc(sizeof(char) * len);
  snprintf(lo->meme_name, len, "%d", imotif);

  // make a copy of the logodds matrix and scale it to [0..range] 
  create_2array(lo->logodds, double, w, lo->alen);
  // make like when we print it out!
  for (i = 0; i < w; i++) {
    for (j = 0; j < lo->alen; j++) {
      lo->logodds(i,j) = NINT(100 * logodds(i,j));
    }
  }
  scale_lo(&lo, 1, PVAL_RANGE); // scale 
  make_double_lo(&lo, 1);       // make a double-letter logodds matrix 

  return(lo);
} // create_lo 

/**********************************************************************/
/*
  print_block_diagrams

  Tile the dataset sequences with the motifs in los[] and print
  the block diagrams with the p-value of the product of p-values.
*/
/**********************************************************************/
static void print_block_diagrams(
  MODEL *model,       // the model 
  DATASET *dataset,   // the dataset 
  int nmotifs,        // number of motifs 
  double **pv,        // p-value tables for each motif 
  FILE *outfile
)
{
  int i;
  STYPE stype = alph_has_complement(dataset->alph) ? (model->invcomp ? Combine : Norc) : Unstranded;
  XLATE_T *xlate = NULL;        // don't translate DNA 
  bool best_motifs = false;     // use all motifs 
  double m_thresh = 1e-4;       // show motifs over p-value 0.0001 
  double w_thresh = m_thresh;   // show strong motifs only 
  bool use_seq_p = false;       // use postion p-values 
  char *f = "%-*s%s %8s  %s\n"; // format 

  for (i=0; i<PAGEWIDTH; i++) { fputc('-', outfile); } fputc('\n', outfile);
  fprintf(outfile, "\tCombined block diagrams:");
  fprintf(outfile, " non-overlapping sites with p-value < %6.4f\n", m_thresh);
  for (i=0; i<PAGEWIDTH; i++) { fputc('-', outfile); } fputc('\n', outfile);
  fprintf(outfile, f, MSN, "SEQUENCE NAME", "", "COMBINED P-VALUE", "MOTIF DIAGRAM");
  fprintf(outfile, f, MSN, "-------------", "", "----------------", "-------------");

  // Handle noise samples in NZ mode.
  int n = (dataset->objfun == NZ) ? dataset->n_group[0] : dataset->n_samples;
  for (i=0; i<n; i++) {
    SAMPLE *s = dataset->input_order[i];
    char *name = s->sample_name;
    int lseq = s->length;
    char *sequence = s->seq;
    int col, offset;
    TILING tiling = score_tile_diagram(NULL, dataset->alph, xlate, sequence, lseq, los, nmotifs,
      stype, false, best_motifs, true, pv, m_thresh, w_thresh,
      use_seq_p, false, NULL);
    fprintf(outfile, "%-*.*s %16.2e  ", MSN, MSN, name, tiling.pv);
    col = 43; offset = 0;
    while (tiling.diagram[offset] != '\0') {
      int j, max, last;
      max = offset + (PAGEWIDTH - col - 1);
      for (j = offset, last = offset; j < max; j++) {
        if (tiling.diagram[j] == '\0') {
          last = j;
          break;
        }
        if (tiling.diagram[j] == '_') last = j+1;
      }
      if (tiling.diagram[last] != '\0') {
        fprintf(outfile, "%.*s\\\n", (last-offset), tiling.diagram+offset); col = 0;
        fputs("    ", outfile); col += 4;
      } else {
        fprintf(outfile, "%s\n", tiling.diagram+offset); col = 0;
      }
      offset = last;
    }
    free_tiling(tiling);
  } // sequence 

  // print a final line of hyphens 
  for (i=0; i<PAGEWIDTH; i++) { fputc('-', outfile); } fprintf(outfile, "\n\n");

} // print_block_diagrams 

/**********************************************************************/
/*
  print_meme_scanned_sites_xml

  Tile the dataset sequences with the motifs in los[] and print
  the XML for the motif occurence strings.
*/
/**********************************************************************/
static void print_meme_scanned_sites_xml(
  MODEL *model,     // the model 
  DATASET *dataset, // the dataset 
  int nmotifs,      // number of motifs 
  double **pv,      // p-value tables for each motif 
  FILE *outfile     // write XML to this file handle 
)
{
  int i_seq;
  STYPE stype = alph_has_complement(dataset->alph) ? (model->invcomp ? Combine : Norc) : Unstranded;
  XLATE_T *xlate = NULL;      // don't translate DNA 
  bool best_motifs = false;   // use all motifs 
  double m_thresh = 1e-4;     // show motifs over p-value 0.0001 
  double w_thresh = m_thresh; // show strong motifs only 
  bool use_seq_p = false;     // use postion p-values 

  fprintf(
    outfile,
    "<scanned_sites_summary p_thresh=\"%.4f\">\n",
    m_thresh
  );
  // Handle noise samples in NZ mode.
  int n = (dataset->objfun == NZ) ? dataset->n_group[0] : dataset->n_samples;
  for (i_seq=0; i_seq<n; i_seq++) {
    SAMPLE *s = dataset->input_order[i_seq];
    char *name = s->sample_name;
    int lseq = s->length;
    char *sequence = s->seq;
    TILING tiling = score_tile_diagram(NULL, dataset->alph, xlate,
      sequence, lseq, los, nmotifs,
      stype, false, best_motifs, true, pv, m_thresh, w_thresh,
      use_seq_p, false, NULL
    );
    // Count the number of sites
    int i_pos = 0;
    int num_sites = 0;
    for (i_pos = 0; i_pos < lseq; i_pos++) {
      if (tiling.hits[i_pos] != 0) {
        num_sites++;
      }
    }
    // Print opening scanned sites tag
    // contains info about sequence.
    fprintf(
      outfile, 
      "<scanned_sites sequence_id=\"sequence_%d\""
      " pvalue=\"%.2e\" num_sites=\"%d\">",
      i_seq,
      tiling.pv,
      num_sites
    );
    for (i_pos = 0; i_pos < lseq; i_pos++) {
      // Print a scanned_site element for each site.
      // contains info about motif occurence.
      int motif_index = tiling.hits[i_pos]; 
      if (motif_index != 0) {
        char* strand = "none";
        if (stype != Unstranded) {
          strand = motif_index > 0 ? "plus" : "minus";
        }
        fprintf(
          outfile, 
          "<scanned_site motif_id=\"motif_%d\" strand=\"%s\""
          " position=\"%d\" pvalue=\"%.2e\"/>\n",
          abs(tiling.hits[i_pos]),
          strand,
          i_pos,
          tiling.pvalues[i_pos]
        );
      }
    }
    // Print closing scanned sites tag
    fprintf( outfile, "</scanned_sites>\n");
    free_tiling(tiling);
  } // sequence 
  fprintf(
    outfile,
    "</scanned_sites_summary>\n"
  );

} // print_meme_scanned_sites_xml 

/**********************************************************************/
/*
  print_log_odds

  Print the log-odds matrix
*/
/**********************************************************************/
static void print_log_odds(
  int imotif,       // motif number 
  char *consensus,  // single-letter consensus
  DATASET *dataset, // the dataset 
  int w,            // width of motif 
  double **logodds, // log-odds matrix 
  double bayes,     // threshold 
  double log_ev,    // log E-value of motif 
  FILE *outfile     // output file 
)
{
  int i, j;
  int n = wps(dataset, w, 0);		// weighted possible sites 
  char *type = "";			// type of matrix 
  double m1, e1;			// for printing significance 

#ifdef PARALLEL
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: at print_print_logodds\n", mpMyID()); fflush(stderr);
#endif
#endif

  // get E-value of motif 
  exp10_logx(log_ev/log(10.0), m1, e1, 1);

  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n\tMotif %s MEME-%d position-specific scoring matrix\n", consensus, imotif);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n");
  fprintf(outfile,
    "log-odds matrix: alength= %d w= %d n= %d bayes= %g E= %3.1fe%+04.0f %s\n",
    alph_size_core(dataset->alph), w, n, bayes, m1, e1, type);

  for (i=0; i < w; i++) {         // site position 
    for (j=0; j < alph_size_core(dataset->alph); j++) { // letter 
      fprintf(outfile, "%6d ", NINT(100*logodds(i,j)));
    }
    fprintf(outfile, "\n");
  }
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n\n");

#ifdef PARALLEL
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: leaving print_print_logodds\n", mpMyID()); fflush(stderr);
#endif
#endif
} // print_log_odds 

/**********************************************************************/
/*
  print_entropy

  Displays the relative entropy of each column of the motif
  as a bar graph and as a logo.
*/
/**********************************************************************/
static void print_entropy(
  bool logo,          // true: prints the logo, false: otherwise 
  MODEL *model,       // the model 
  DATASET *dataset,   // the dataset 
  char *str_space,    // space for printing strand direction 
  FILE *outfile       // stream for output 
)
{
  int i, j;
  int w = model->w;                     // width of motif 
  THETA obs = model->obs;               // observed frequencies 
  double *rentropy = model->rentropy;   // RE of each column 
  double re = w *model->rel;            // motif relative entropy 
  ARRAY_T *back = dataset->back;        // background model 
  char restring[15];                    // print string for re 
  char *consensus;                      // consensus strings 
  double min_freq;                      // minimum background freq 
  double maxre;                         // maximum relative entropy 
  int nsteps;                           // number of steps histogram 

  // get minimum background frequency and maximum relative entropy 
  for (i = 0, min_freq = 1; i < alph_size_core(dataset->alph); i++) {
    if (get_array_item(i, back) < min_freq) min_freq = get_array_item(i, back);
  }
  maxre = -LOG2(min_freq);              // maximum relative entropy 

  // create string containing RE 
  sprintf(restring, "(%.1f bits)", re);

  // print the relative entropy of each column as a bar graph 
  nsteps = 10;
  for (i=0; i<nsteps; i++) {
    double level = maxre - (i * maxre/nsteps);
    fprintf(outfile, (i==0 ? "%*.*s %*.1f " : "%-*.*s %*.1f "), IND, IND,
      (i==0 ? "bits" : i==4 ? "Relative" : i==5 ? "Entropy" :
        i==6 ? restring : ""), IND2, level);
    for (j=0; j<w; j++) {
      if (NINT(nsteps * rentropy[j] / maxre) >= nsteps-i) {
        fputc('*', outfile);
      } else {
        fputc(' ', outfile);
      }
    }
    fputc('\n', outfile);
  }
  fprintf(outfile, "%-*.*s %*.1f ", IND, IND, "", IND2,0.0);
  for (i=0; i<w; i++) fputc('-', outfile);
  fprintf(outfile, "\n\n");
  // get and print the consensus sequences 
  consensus = get_consensus(obs, w, dataset, MAXDEPTH, MINCONS);
  for (i=0; i < MAXDEPTH && i < alph_size_core(dataset->alph); i++) { // print next levels of consensus 
    fprintf(outfile, "%-*.*s %*.0s %*.*s\n", IND, IND,
      (i==0 ? "Multilevel" : i == 1 ? "consensus" : i == 2 ? "sequence" : ""),
      IND2, "", w, w, consensus+(i*w));
  }
  // free up space 
  myfree(consensus);

  // Prints a logo in EPS and PNG format to two files in the output directory 
  if(logo) print_logo(model, dataset);

} // print_entropy 

/**********************************************************************/
/*
  print_logo

  Print the logo of a motif
*/
/**********************************************************************/
static void print_logo(
  MODEL *model,       // the model 
  DATASET *dataset    // the dataset 
)
{
  char* logodir = dataset->output_directory;
  int both_strands = (alph_has_complement(dataset->alph) && model->invcomp);

  if(logodir != NULL) {
    // convert theta to motif 
    MOTIF_T motif;
    memset(&motif, 0, sizeof(MOTIF_T));
    set_motif_id("0", 1, &motif);
    set_motif_strand('+', &motif);
    motif.num_sites  = model->nsites_obs;
    motif.freqs      = convert_matrix(model->obs, model->w, alph_size_core(dataset->alph)); 
    motif.length     = model->w;
    motif.alph       = dataset->alph; 
    motif.flags      = (both_strands ? MOTIF_BOTH_STRANDS : 0);
    motif.evalue     = 0.0;
    motif.complexity = 0.0;
    motif.trim_left = 0;
    motif.trim_right = 0;

    // create the output path
    char *path = NULL;
    Resize(path, strlen(logodir)+29, char);// room for "/logo_ssc<16digts>\0"

    // create logo without small sample correction
    sprintf(path, "%s/logo%d", logodir, model->imotif);
    CL_create1(
      &motif,         // first motif
      false,          // no error bars
      false,          // ssc
      "MEME (no SSC)",// program name
      path,           // output file path
      true,           // output eps
      true            // output png
    );

    if (alph_has_complement(dataset->alph)) {
      reverse_complement_motif(&motif);
      sprintf(path, "%s/logo_rc%d", logodir, model->imotif);
      CL_create1(
        &motif,         // first motif
        false,          // no error bars
        false,          // ssc
        "MEME (no SSC)",// program name
        path,           // output file path
        true,           // output eps
        true            // output png
      );
    }

    myfree(path);
    free_matrix(motif.freqs);
  }
} // print_logo

/**********************************************************************/
/*
  print_theta

    format=1    floating point; pos x letter
    format=2    1 digit; letter x pos

  Print the probability array.
*/
/**********************************************************************/
void print_theta(
  int imotif,       // number of motif 
  char *consensus,  // single-letter consensus of motif
  int format,       // 1 = floating point
                    // 2 = integer 
  int nsites,       // number of sites (discrete) 
  THETA theta,      // theta 
  int w,            // width of motif 
  double log_ev,    // log motif E-value 
  char *str_space,  // space for printing strand direction 
  DATASET *dataset, // the dataset 
  FILE *outfile     // file to print to 
)
{
  int i, j;

  if (format == 1) {
    double e1, m1;
    exp10_logx(log_ev/log(10.0), m1, e1, 1);
    for (i=0; i<PAGEWIDTH; i++) {
      fputc('-', outfile);
    }
    if (consensus) {
      fprintf(outfile, "\n\tMotif %s MEME-%d position-specific probability matrix\n", consensus, imotif);
    } else {
      fprintf(outfile, "\n\tMotif %d position-specific probability matrix\n", imotif);
    }
    for (i=0; i<PAGEWIDTH; i++) {
      fputc('-', outfile);
    }
    fprintf(
      outfile,
      "\nletter-probability matrix: alength= %d w= %d "
      "nsites= %d E= %3.1fe%+04.0f ",
      alph_size_core(dataset->alph), w, nsites, m1, e1
    );
    fprintf(outfile, "\n");
    for (i=0; i < w; i++) {
      for (j=0; j < alph_size_core(dataset->alph); j++) {
        fprintf(outfile, "%9.6f ", theta(i, j));
      }
      fprintf(outfile, "\n");
    }
    for (i=0; i<PAGEWIDTH; i++) {
      fputc('-', outfile);
    }
    fprintf(outfile, "\n");

  } else if (format == 2) {
    // print theta: rows = letter; cols = position in motif; 1-digit integer 
    for (i=0; i < alph_size_core(dataset->alph); i++) {
      // print the letter 
      fprintf(
        outfile, "%-*.*s%*c  ", IND, IND,
        (i==0 ? "Simplified" : i==1 ? "pos.-specific" : i==2 ?
          "probability" : i==3 ? "matrix" : "" ), IND2, alph_char(dataset->alph, i)
      );
      for (j=0; j < w; j++) {
        int k = NINT(10.0 * theta(j,i)); // round to 1 digit 
        if (k == 0) {
          fprintf(outfile, ":");         // print 0 as colon 
        } else {
          fprintf(outfile, "%1x", k);    // print 1 digit 
        }

      }
      fprintf(outfile, "\n");
    }
  }

  fprintf(outfile, "\n");
}

/**********************************************************************/
/*
  print_zij
*/
/**********************************************************************/
void print_zij(
  DATASET *dataset,   // the dataset 
  MODEL *model        // the model 
)
{
  int i, j;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->input_order;
  FILE *out=stdout;

  fprintf(out, "z_ij: lambda=%f ll=%f\n", model->lambda, model->ll);
  for (i=0; i<n_samples; i++) {         // sequence 
    int lseq = samples[i]->length;
    Z_T *zi = samples[i]->z;		// zi[j], j in [-lseq...+lseq]
    int w = model->w;
    fprintf(out, ">%s\nz : ", samples[i]->sample_name);
    for (j=0; j<lseq-w+1; j++) {    	// position 
      int k = j+1;			// Z_i = k
      double z = model->invcomp ? MIN(1.0,Zi(-k)+Zi(k)) : Zi(k);
      int zij = NINT(10 * z);           // round z 
      fprintf(out, "%1x", zij);
    } // position
    // print s0 and s1 for backwards compatibility
    if (model->invcomp) {
      fprintf(out, "\ns0: ");
      for (j=0; j<lseq-w+1; j++) {    	// position 
	int k = j+1;			// Z_i = k
	double z = Zi(k);
	int zij = NINT(10 * z);  // round z 
	fprintf(out, "%1x", zij);
      } // position
      fprintf(out, "\ns1: ");
      for (j=0; j<lseq-w+1; j++) {    	// position 
	int k = j+1;			// Z_i = k
	double z = Zi(-k);
	int zij = NINT(10 * z);  // round z 
	fprintf(out, "%1x", zij);
      } // position
    }
    fprintf(out, "\n");
  } // sequence 
  printf("\n");
} // print_zij 

/**********************************************************************/
/*
  print_wij
*/
/**********************************************************************/
void print_wij(
  DATASET *dataset      // the dataset 
)
{
  int i,j;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->input_order;

  printf("w_ij:\n");
  for (i=0; i<n_samples; i++) {               // sequence 
    int len = samples[i]->length;
    WEIGHTS_T *weights = samples[i]->weights;
    printf(">%s\n", samples[i]->sample_name);
    for (j=0; j<len; j++) {                   // position 
      int w = NINT(10 * weights[j]);
      printf("%1x", w);
    }
    printf("\n");
  }
  printf("\n");
}

/**********************************************************************/
/*      get_consensus

  Get the consensus string from a motif.

  For each position, N consensus letters are found.
  If no letter has probability > min_prob,
        'x' is written for the first consensus
  letter and ' ' in the others.
        Otherwise, N letters are written in decreasing order of
  probability until one with min_prob is reached, and then ' ' is
  written for the remaining letters.
*/
/**********************************************************************/
char *get_consensus(
  THETA theta,      // motif theta 
  int w,            // width of motif 
  DATASET *dataset, // the dataset 
  int N,            // number of letters for each position 
  double min_prob   // minimum cumulative prob for N letters 
)
{
  int i, j, n;
  char *string = NULL;
  LETTER_VALUE *letterv = NULL;

  
  Resize(letterv, alph_size_core(dataset->alph), LETTER_VALUE);
  Resize(string, w*N+2, char);

  for (i=0; i < w; i++) {   // position in motif 
    int maxj[MAXDEPTH];     // array of max indices in Theta 

    // sort letters at position i, largest first 
    for (j = 0; j < alph_size_core(dataset->alph); j++) {
      letterv[j].letter = alph_char(dataset->alph, j);
      letterv[j].value = theta(i, j);
    }
    qsort(letterv, alph_size_core(dataset->alph), sizeof(LETTER_VALUE), lv_compare);

    // set up the consensus strings for position i 
    for (n = 0; n < N; n++) {                  // current depth 
      if (n < alph_size_core(dataset->alph)) {
        if (letterv[n].value < min_prob) {
          string[(n*w)+i] = (n==0 ? alph_wildcard(dataset->alph) : ' '); // below cutoff 
        } else {
          string[(n*w)+i] = letterv[n].letter; // set n'th consensus 
        }
      } else string[(n*w)+i] = ' ';      
    }
  }
  string[N*w] = '\0';     // terminate string 
  // cleanup
  free(letterv);

  return string;
} // get_consensus 


/**********************************************************************/
/*      get_regexp

  Get a regular expression using the same rules as for
  the multilevel consensus sequence.
*/
/**********************************************************************/
static char *get_regexp(
  MODEL *model,     	// the model 
  DATASET *dataset,   	// the dataset 
  bool prosite          // true: [AK]-E-[EG]. false: [AK]E[EG] 
)
{
  THETA obs = model->obs; // motif observed theta 
  int w = model->w;   	  // width of motif 
  int i, j;
  char *pcons = get_consensus(obs, w, dataset, MAXDEPTH, MINCONS);
  char *re = NULL;    	  // RE string to return 

  // regular expression string 
  Resize(re, w*(MAXDEPTH+6), char);

  // create the regular expression from the "packed" consensus 
  int pos = 0;
  for (i=0; i<w; i++) {
    // see if there is more than one letter in consensus for this column 
    if (i>0 && prosite) re[pos++] = '-';
    if (MAXDEPTH > 1 && pcons[w + i] != ' ') re[pos++] = '[';
    for (j=0; j<MAXDEPTH; j++) {  // copy consensus 
      char a = pcons[w*j + i];
      if (a == ' ') break;
      re[pos++] = a;
    }
    if (MAXDEPTH > 1 && pcons[w + i] != ' ') re[pos++] = ']';
  }
  if (prosite) re[pos++] = '.';
  re[pos++] = '\0';
  free(pcons);

  return re;
} // get_regexp 

/**********************************************************************/
/*
  print_candidates

  Print the candidate models found.

*/
/**********************************************************************/
static void print_candidates(
  CANDIDATE *candidates, // list of candidate models IN 
  DATASET *dataset,      // the dataset IN 
  int max_w,             // maximum width for motifs IN 
  FILE* outfile          // file for text output IN 
)
{
  int i, w;
  int hdr = 35;
  int tail = PAGEWIDTH - hdr;
  double m, e;           // for printing signficance 

  fprintf(outfile, "\nCandidate motifs:\n");
  fprintf(outfile, "width nsites  ll        signif     consensus\n");
  fprintf(outfile, "----- ------  --------- ---------- ---------\n");

  for (w=1; w<=max_w; w++) {

    if (candidates[w].sig > 1) continue;  // skip unused candidates 

    exp10_logx(candidates[w].sig/log(10.0), m, e, 3);

    fprintf(
      outfile,
      "%5d %6.1f %c%9.0f %5.3fe%+04.0f ",
      w,
      candidates[w].lambda * wps(dataset, w, 0),
      (candidates[w].pal ? 'P' : ' '),
      candidates[w].ll,
      m, e
    );
    fprintf(outfile, "%-*.*s\n", tail, tail, candidates[w].cons);
    for (i=tail; i < w; i+=tail) {
      fprintf(
        outfile, "%*.*s%-*.*s\n", hdr, hdr, "",
        tail, tail, candidates[w].cons+i
      );
    }
  }
} // print_candidates 

/**********************************************************************/
/*
  print_dataset_summary
*/
/**********************************************************************/
void print_dataset_summary (
  DATASET *dataset, // the dataset IN 
  FILE *outfile     // where to print IN 
)
{
  int i, pcol;
  int w = MSN + 15;
  char *datafile = dataset->datafile; // name of the training set file 
  char *negfile = dataset->negfile;   // name of negative example file 
  ALPH_T *alph = dataset->alph;       // alphabet of dataset

  // set up printing spacers 
  Resize(stars, PAGEWIDTH+1, char);
  for (i=0; i<PAGEWIDTH; i++) {
    stars[i] = '*';
  }
  stars[PAGEWIDTH] = '\0';

  // new-style alphabet declaration
  if (!alph_is_builtin(alph)) {
    fprintf(outfile, "%s\n", stars);
    if (mpMyID() == 0) alph_print_header(alph, outfile);
    fprintf(outfile, "%s\n", stars);
    if (mpMyID() == 0) alph_print(alph, false, outfile);
    fprintf(outfile, "%s\n\n", stars);
  }

  // announce the training set 
  fprintf(outfile, "%s\nTRAINING SET\n%s\n", stars, stars);

  // print name of file and alphabet 
  fprintf(outfile, "PRIMARY SEQUENCES= %s\n", datafile);

  // print name of negative dataset 
  if (negfile) {
    fprintf(outfile, "CONTROL SEQUENCES= %s\n", negfile);
  }

  // old-style alphabet declaration
  if (alph_is_builtin_dna(alph)) {
    fprintf(outfile, "ALPHABET= ACGT\n");
  } else if (alph_is_builtin_rna(alph)) {
    fprintf(outfile, "ALPHABET= ACGU\n");
  } else if (alph_is_builtin_protein(alph)) {
    fprintf(outfile, "ALPHABET= ACDEFGHIKLMNPQRSTVWY\n");
  }

  //
  //  print a table of sequence lengths
  //
  if (dataset->n_samples <= dataset->brief) {

    // print table header 
    for (pcol = w; pcol < 80; pcol += w) {
      fprintf(
	outfile,
	"%-*.*s %6s %6s%2s",
	MSN, MSN, "Sequence name", "Weight", "Length", " "
      );
    }
    fprintf(outfile, "\n");
    for (pcol = w; pcol < 80; pcol += w) {
      fprintf(
	outfile,
	"%-*.*s %6s %6s%2s",
	MSN, MSN, "-------------", "------", "------", " "
      );
    }
    fprintf(outfile, "\n");

    // print table columns 
    pcol = 0;
    // Handle noise samples in NZ mode.
    int n_samples = (dataset->objfun == NZ) ? dataset->n_group[0] : dataset->n_samples;
    for (i=0; i<n_samples; i++) {
      SAMPLE *sample = dataset->input_order[i];
      char *sample_name = sample->sample_name;
      double wgt = sample->sw;
      int lseq = sample->length;
      // print the sample name and its length 
      pcol += w;          // new line for print out? 
      if (pcol >= 80) {
	fprintf(outfile, "\n");
	pcol = w;
      }
      fprintf(
	outfile,
	"%-*.*s %6.4f %6d%2s",
	MSN, MSN, sample_name, wgt, lseq, " "
      );
    }
  } // sequence lengths

  // finish section 
  fprintf(outfile, "\n%s\n\n", stars);
} // print_dataset_summary 

/**********************************************************************/
/*
  print_command_summary

  Print the command line summary
*/
/**********************************************************************/
void print_command_summary(
  MODEL *model,     // the model IN 
  DATASET *dataset, // the dataset IN 
  FILE *outfile     // where to print IN 
)
{
  int i, pcol;
  char evt_string[20];

  if (dataset->evt == BIG) {
    strcpy(evt_string, "inf");
  } else {
    sprintf(evt_string, "%8g", dataset->evt);
  }

  fprintf(
    outfile,
    "%s\nCOMMAND LINE SUMMARY\n%s\n"
    "This information can also be useful in the event you wish to report a\n"
    "problem with the MEME software.\n\n"
    "command: %s\n\n"
    "model:  mod=      %8s    nmotifs=  %8d    evt=      %8s\n"
    "objective function:           em=       %s%s\n"
    "                              starts=   %s\n",
    stars, stars,
    dataset->command,
    dataset->mod, dataset->nmotifs, evt_string,
    (dataset->objfun == Classic) ? "E-value of product of p-values"
      : (dataset->objfun == NC) ? "E-value of log likelihood ratio"
	: (dataset->objfun == DE) ? "Differential Enrichment"
	  : (dataset->objfun == SE) ? "Selective"
	    : (dataset->objfun == CE) ? "Central Enrichment"
	      : (dataset->objfun == CD) ? "Central Enrichment: p-value of mean distance"
		: (dataset->objfun == NZ) ? "Noise-injected"
		  : "Unknown objective function",
    (dataset->test == mHG) ? " mHG"
      : (dataset->test == mBN) ? " mBN"
        : (dataset->test == mRS) ? " mRS"
          : "", 
    (dataset->objfun == Classic && !dataset->use_llr) ? "E-value of product of p-values"
      : (dataset->objfun == CE) ? "Central Enrichment binomial test"
        : (dataset->objfun == CD) ? "Mean distance of best site from center"
          : "log likelihood ratio (LLR)"
    );

  // Print how strands are handled.
  if (alph_has_complement(dataset->alph)) {
    fprintf(outfile, "strands: +");
    if (model->invcomp) {
      fprintf(outfile, " -");
    }
    fprintf(outfile, "\n");
  }

  // Print width range.
  fprintf(outfile, "width:  minw=     %8d    maxw=     %8d\n", 
    model->min_w, model->max_w);

  fprintf(outfile,
    "nsites: minsites= %8g    maxsites= %8g    wnsites=  %8g\n"
    "theta:  spmap=    %8s    spfuzz=   %8g\n", 
    dataset->min_nsites, dataset->max_nsites, dataset->wnsites,
    dataset->mapname, dataset->map_scale);

  // print global search parameters
  //fprintf(outfile,
  //  "global: substring=%8s    branching=%8s    wbranch=  %8s\n",
  //  (dataset->p_point->c > 0) ? "no" : "yes", 
  //  (dataset->branch_params->point_branch != NO_POINT_B) ? "yes" : "no",
  //  yesno[dataset->branch_params->w_branch]);

  // print branching parameters
  //if (dataset->branch_params->point_branch != NO_POINT_B){
  //  fprintf(outfile,
  //  "        bfactor=  %8i    heapsize= %8i\n",
  //  dataset->branch_params->bfactor,
  //   dataset->main_hs);
  //}

  // print local search parameters
  fprintf(outfile,
    "em:     prior=   %9s    b=        %8g    maxiter=  %8d\n"
    "        distance= %8g\n",
    dataset->priorname, dataset->beta, dataset->maxiter,
    dataset->distance);

  // Print trimming.
  if (dataset->ma_adj) {
    fprintf(
      outfile,
      "trim:   wg=       %8g    ws=       %8g    endgaps=  %8s\n",
      dataset->wg, dataset->ws, yesno[dataset->endgaps]
    );
  }

  // print properties of the dataset
  int n_samples = (dataset->objfun == NZ) ? dataset->n_group[0] : dataset->n_samples;
  fprintf(outfile,
    "data:   n=        %8d    N=        %8d\n", 
    get_total_res(dataset), n_samples);

  // print sampling and efficiency thresholds
  fprintf(outfile,
    "sample: seed=     %8d    hsfrac=   %8g\n"
    "        searchsize=%7d    norand=   %8s    csites=   %8d\n",
    dataset->seed, dataset->hsfrac,
    dataset->search_size, (dataset->no_rand ? "yes" : "no"), 
    dataset->classic_max_nsites
  );

  // print sequence priors
  if (dataset->plib_name) {
    fprintf(outfile, "Dirichlet mixture priors file: %s\n", dataset->plib_name);
  }

  // print dataset frequencies of letters in alphabet 
  fprintf(outfile, "Letter frequencies in dataset:\n");
  for (i = 0, pcol = 0; i < alph_size_core(dataset->alph); i++) {
    pcol += 8;          // start of next printed thing 
    if (pcol >= PAGEWIDTH) {pcol=8; fprintf(outfile, "\n");}
    fprintf(outfile, "%c %.3g ", alph_char(dataset->alph, i), dataset->res_freq[i]);
  }

  // print background frequencies of letters in alphabet 
  fprintf(
    outfile,
    "\nBackground letter frequencies (from file %s):\n",
    dataset->bfile ? dataset->bfile : "dataset with add-one prior applied"
  );
  for (i = 0, pcol = 0; i < alph_size_core(dataset->alph); i++) {
    pcol += 8;          // start of next printed thing 
    if (pcol >= PAGEWIDTH) {
      pcol=8;
      fprintf(outfile, "\n");
    }
    fprintf(outfile, "%c %.3g ", alph_char(dataset->alph, i), get_array_item(i, dataset->back));
  }
  // Print background model order
  fprintf(outfile, "\nBackground model order: %d\n", dataset->back_order);
  fprintf(outfile, "%s\n", stars);

} // print_command_summary 

/**********************************************************************/
/*
  meme_score_sequence

  Compute the sequence score for a motif.

  Returns the sequence score.
*/
/**********************************************************************/
static double meme_score_sequence(
  uint8_t *eseq,     // integer-coded sequence to score 
  int length,        // length of the sequence 
  int w,             // width of motif 
  double **logodds1, // log-odds matrix: LOG2(m_ij/b_j) 
  double **logodds2  // log-odds matrix: LOG2(m_ij/n_ij) 
)
{
  int i, j;
  double best = LITTLE;     // sequence score 
  double score, sc1, sc2;
  double loge2 = log(2);

  // score the sequence with motif 
  for (i=0; i <= length - w; i++) { // site start 
    // calculate score of subsequence 
    for (j=0, sc1=0, sc2=0; j<w; j++) { // position in sequence 
      sc1 += logodds1(j, eseq[i+j]);
      if (logodds2) sc2 += logodds2(j, eseq[i+j]);
    } // subsequence 
    score = logodds2 ? -LOGL_SUM(-sc1*loge2, -sc2*loge2)/loge2 : sc1;
    best = MAX(score, best);
  } // site start 

  return best;

} // meme_score_sequence 

/**********************************************************************/
/*
        s_compare

        Compare two scores in descending order.  Return <0 >0
        if the first is <, > the second.  If they are equal,
        resolves ties by returning <0 if the first has smaller class.
*/
/**********************************************************************/
static int s_compare(
  const void *v1,
  const void *v2
)
{
  const SORTED_SCORE * s1 = (const SORTED_SCORE *) v1;
  const SORTED_SCORE * s2 = (const SORTED_SCORE *) v2;
  double diff = s1->score - s2->score;
  if (diff == 0) diff = (double) (s1->class - s2->class);
  return ((diff > 0) ? -1 : ( (diff < 0) ? 1 : 0) );
} // s_compare 

/**********************************************************************/
/*
        lv_compare

        Compare two letter values in decending order
*/
/**********************************************************************/
static int lv_compare(
  const void *v1,
  const void *v2
)
{
  const LETTER_VALUE *lv1, *lv2;
  lv1 = (const LETTER_VALUE *) v1;
  lv2 = (const LETTER_VALUE *) v2;
  if (lv1->value > lv2->value) return -1;
  if (lv1->value < lv2->value) return 1;
  // value is equal
  if (lv1->letter < lv2->letter) return -1;
  if (lv1->letter > lv2->letter) return 1;
  // letter is equal
  return 0;
} // lv_compare 

/**********************************************************************/
/*
  print_sites

  Print the sites making up the model.
*/
/**********************************************************************/
static void print_sites(
  DATASET *dataset,     // the dataset 
  MODEL *model,		// the model 
  char *consensus, 	// single-letter consensus of motif
  int format,		// 0=BLOCKS; 1=FASTA 
  char *com,		// comment to append 
  FILE *outfile		// where to print 
)
{
  int i, j;
  int w = model->w;               // width of motif 
  P_PROB sites = model->maxima;   // sites "defining" model 
  int n = model->nsites_dis;      // number of sites 
  char *ftype = (format==0 ? "BLOCKS" : "FASTA");

  // print header 
  for (i=0; i<PAGEWIDTH; i++) { fputc('-', outfile); } fputc('\n', outfile);
  fprintf(outfile, "\tMotif %s MEME-%d in %s format%s\n", consensus, model->imotif, ftype, com);
  for (i=0; i<PAGEWIDTH; i++) { fputc('-', outfile); } fputc('\n', outfile);

  // print the sites 
  if (format == 0) fprintf(outfile, "BL   MOTIF %s width=%d seqs=%d\n",
    consensus, w, n);
  for (i=0; i<n; i++) {                                  // print each site 
    int seqno = sites[i].x;                              // sequence number 
    SAMPLE *s = dataset->input_order[seqno];             // sequence 
    bool ic = sites[i].ic;                               // strand direction 
    int y = sites[i].y;                                  // location of site 
    int off = ic ? s->length-w-y : y;                    // - strand offset from rgt. 
    uint8_t *res = ic ? s->resic+off : s->res+off;       // integer sequence 
    double weight = s->sw;                               // sequence weight 
    // double weight = sites[i].prob;

    // print sequence name and position of site 
    if (format == 0) {      // BLOCKS format 
      fprintf(outfile, "%-*.*s ( %4d) ", MSN, MSN, s->sample_name, y+1);
    } else {                // FASTA format 
      fprintf(outfile,">%-*.*s pos %4d\n", MSN, MSN, s->sample_name, y+1);
    }

    // print site 
   for (j=0; j<w; j++) { fputc(alph_char(dataset->alph, res[j]), outfile); }
    if (format == 0) {      // BLOCKS format 
      fprintf(outfile, "  %g ", weight);
    }
    fputc('\n', outfile);
  } // print each site 
  if (format == 0) {
    fprintf(outfile, "//\n\n");
  }
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n\n");

} // print_sites 

/**********************************************************************/
/*
  print_summary

  Print the summary of all the motifs.
*/
/**********************************************************************/
void print_summary(
  MODEL *model,     // the model IN 
  DATASET *dataset, // the dataset IN 
  int nmotifs,      // number of motifs IN 
  FILE *outfile     // where to print IN 
)
{
  // print the motif block diagrams using all the motifs 
  fprintf(outfile, "\n\n%s\nSUMMARY OF MOTIFS\n%s\n\n", stars, stars);
  print_block_diagrams(model, dataset, nmotifs, pv, outfile);
  fprintf(outfile, "%s\n\n", stars);
} // print_summary 

/**********************************************************************
  print_meme_file_xml

  Print MEME results in XML format. See DTD embeded in
  meme-dtd.h for description of MEME document.
 **********************************************************************/
void print_meme_file_xml(
  MODEL *model,                   // the model IN 
  DATASET *dataset,               // the dataset IN 
  DATASET *neg_dataset,  	  // the control dataset IN 
  int nmotifs,                    // number of motifs IN 
  MOTIF_SUMMARY *motif_summaries, // list of final motif properties IN 
  char *stopping_reason,          // description of reason for stopping IN 
  char* xml_filename              // full path to output file for xml IN 
)
{
  FILE* outfile = fopen(xml_filename, "w"); //FIXME CEG check for errors
  print_meme_header_xml(outfile);

  fputs("<MEME version=\"" VERSION "\" release=\"" ARCHIVE_DATE "\">\n", outfile);
  print_meme_training_set_xml(dataset, neg_dataset, outfile);
  print_meme_model_xml(model, dataset, stopping_reason, outfile);
  print_meme_motifs_xml(model, dataset, nmotifs, motif_summaries, outfile);
  if (nmotifs > 0 && (dataset->n_samples <= dataset->brief)) {
    print_meme_scanned_sites_xml(model, dataset, nmotifs, pv, outfile);
  }
  fprintf(outfile, "</MEME>\n");
  fclose(outfile);
} // print_meme_file_xml 

/**********************************************************************
  print_meme_header_xml

  Print the DTD and XML header for MEME
**********************************************************************/
static void print_meme_header_xml(FILE *outfile) {
  fputs(meme_dts, outfile);
  fputs("<!-- Begin document body -->\n", outfile);
}

/**********************************************************************
  sym_2_id

  Given a buffer return a string that can be used as the id for a symbol.
**********************************************************************/
static char* sym_2_id(char symbol, STR_T *buffer) {
  str_clear(buffer);
  if ((symbol >= 'A' && symbol <= 'Z') || (symbol >= 'a' && symbol <= 'z') || symbol == ':' || symbol == '_') {
    str_appendf(buffer, "%c", symbol); // these can be used as-is
  } else if ((symbol >= '0' && symbol <= '9') || symbol == '.' || symbol == '-') {
    str_appendf(buffer, "s%c", symbol); // these are valid but not NameStartChar
  } else {
    // any other symbols (eg '*') get encoded as hexadecimal
    str_appendf(buffer, "s%2X", (uint8_t)symbol);
  }
  return str_internal(buffer);
}

/**********************************************************************
  sym_equals

  Given a buffer return a string of all comprising symbols.
**********************************************************************/
static char* sym_equals(ALPH_T *alph, int index, STR_T *buffer) {
  int i;
  str_clear(buffer);
  for (i = 0; i < alph_ncomprise(alph, index); i++) {
    str_appendf(buffer, "%c", alph_char(alph, alph_comprise(alph, index, i)));
  }
  return str_internal(buffer);
}

/**********************************************************************
  print_meme_training_set_xml

  Print XML elements for the training set. See DTD embeded in
  print_meme_header_xml for description of the training set XML syntax.
**********************************************************************/
static void print_meme_training_set_xml(
  DATASET *dataset, // the dataset IN 
  DATASET *neg_dataset, // the control dataset IN 
  FILE* outfile     // file for output IN 
) {
  int i = 0;

  // allocate a buffer for create ids and escaping XML
  STR_T *b = str_create(10);

  // Print training set
  // Have to split print to keep from overwriting xmlify buffer.
  fprintf(
    outfile,
    "<training_set primary_sequences=\"%s\" primary_count=\"%d\" primary_positions=\"%d\" ", 
    xmlify(dataset->datafile, b, true),
    dataset->n_samples, dataset->total_res
  );
  fprintf(
    outfile,
    "control_sequences=\"%s\" control_count=\"%d\" control_positions=\"%d\">\n",
    xmlify(dataset->negfile, b, true),
    neg_dataset ? neg_dataset->n_samples : 0,
    neg_dataset ? neg_dataset->total_res : 0
  );

  alph_print_xml(dataset->alph, "alphabet", "", "", outfile);

  // Don't print sequences in -brief mode
  if (dataset->n_samples <= dataset->brief) {
    // Print sequences in training set
    // Handle noise samples in NZ mode.
    int n_samples = (dataset->objfun==NZ) ? dataset->n_group[0] : dataset->n_samples;
    for(i = 0; i < n_samples; i++) {
      SAMPLE *s = dataset->input_order[i];
      char *name = s->sample_name;
      fprintf(
	outfile,
	"<sequence "
	"id=\"sequence_%d\" "
	"name=\"%s\" "
	"length=\"%ld\" "
	"weight=\"%f\" "
	"/>\n",
	i,
	xmlify(name, b, true),
	dataset->input_order[i]->length,
	dataset->input_order[i]->sw
      );
    }
  } // not -brief

  fprintf(outfile, "<letter_frequencies>\n");
  fprintf(outfile, "<alphabet_array>\n");
  for(i = 0; i < alph_size_core(dataset->alph); i++) {
    fprintf(
      outfile,
      "<value letter_id=\"%s\">%.3g</value>\n",
	alph_xml_id(dataset->alph, i, b),
	dataset->res_freq[i]
      );
  }
  fprintf(outfile, "</alphabet_array>\n");
  fprintf(outfile, "</letter_frequencies>\n");

  fprintf(outfile, "</training_set>\n");
  str_destroy(b, false);
}

/**********************************************************************
  print_meme_model_xml

  Print XML elements for the model. See DTD embeded in
  print_meme_header_xml for the description of model XML syntax.
***********************************************************************/
static void print_meme_model_xml(
  MODEL *model,          // the model IN 
  DATASET *dataset,      // the dataset IN 
  char* stopping_reason, // reason for stopping IN 
  FILE* outfile          // output file IN 
) {
  char evt_string[20];
  STR_T *b = str_create(10);

  if (dataset->evt == BIG) {
    strcpy(evt_string, "inf");
  } else {
    sprintf(evt_string, "%g", dataset->evt);
  }

  char *hostname;
  hostname = HOSTNAME;
  if (hostname[0] == '\0') {
    hostname = "unknown";
  }

  fprintf(
    outfile,
    "<model>\n"
    "<command_line>%s</command_line>\n",
    xmlify(dataset->command, b, false)
  );
  fprintf(
    outfile,
    "<host>%s</host>\n"
    "<type>%s</type>\n"
    "<nmotifs>%d</nmotifs>\n"
    "<evalue_threshold>%s</evalue_threshold>\n"
    "<object_function>%s%s</object_function>\n"
    "<spfun>%s</spfun>\n"
    "<min_width>%d</min_width>\n"
    "<max_width>%d</max_width>\n",
    xmlify(hostname, b, false),
    dataset->mod, dataset->nmotifs, evt_string,
    (dataset->objfun == Classic) ? "E-value of product of p-values"
      : (dataset->objfun == NC) ? "E-value of log likelihood ratio"
	: (dataset->objfun == DE) ? "Differential Enrichment"
	  : (dataset->objfun == SE) ? "Selective"
	    : (dataset->objfun == CE) ? "Central Enrichment"
	      : (dataset->objfun == CD) ? "Central Enrichment: p-value of mean distance"
		: (dataset->objfun == NZ) ? "Noise-injected"
		  : "Unknown objective function",
    (dataset->test == mHG) ? " mHG"
      : (dataset->test == mBN) ? " mBN"
        : (dataset->test == mRS) ? " mRS"
          : "", 
    (dataset->objfun == Classic && !dataset->use_llr) ? "E-value of product of p-values"
      : (dataset->objfun == CE) ? "Central Enrichment binomial test"
        : (dataset->objfun == CD) ? "Mean distance of best site from sequence center"
          : "log likelihood ratio (LLR)",
    model->min_w, model->max_w
  );
  if (dataset->ma_adj) {
    fprintf(
      outfile,
      "<wg>%g</wg>\n"
      "<ws>%g</ws>\n"
      "<endgaps>%s</endgaps>\n",
      dataset->wg, dataset->ws, yesno[dataset->endgaps]);
  }
  fprintf(
    outfile,
    "<substring>%s</substring>\n",
    (dataset->p_point->c > 0) ? "no" : "yes");
  fprintf(
    outfile,
    "<minsites>%g</minsites>\n"
    "<maxsites>%g</maxsites>\n"
    "<wnsites>%g</wnsites>\n"
    "<spmap>%s</spmap>\n",
    dataset->min_nsites, dataset->max_nsites, dataset->wnsites,
    xmlify(dataset->mapname, b, false)
  );
  fprintf(
    outfile,
    "<spfuzz>%g</spfuzz>\n"
    "<prior>%s</prior>\n"
    "<beta>%g</beta>\n"
    "<maxiter>%d</maxiter>\n"
    "<distance>%g</distance>\n"
    "<num_positions>%d</num_positions>\n"
    "<seed>%d</seed>\n"
    "<hsfrac>%g</hsfrac>\n"
    "<searchsize>%d</searchsize>\n"
    "<maxsize>%d</maxsize>\n"
    "<norand>%s</norand>\n",
    dataset->map_scale,
    xmlify(dataset->priorname, b, false),
    dataset->beta, dataset->maxiter, dataset->distance,
    dataset->total_res, dataset->seed, 
    dataset->hsfrac, dataset->search_size, dataset->max_size,
    dataset->no_rand ? "yes" : "no"
  );
  if (dataset->objfun == Classic) fprintf(outfile, "<csites>%d</csites>\n", dataset->classic_max_nsites);
  fprintf(outfile, "<strands>%s</strands>\n", 
      (alph_has_complement(dataset->alph) ? (model->invcomp ? "both" : "forward") : "none"));
  fprintf(outfile, "<brief>%d</brief>\n", dataset->brief);

  // The position-specific priors file.
  fprintf(outfile, "<psp_file>");
  if (dataset->pspfile) {
    fprintf(outfile, "%s", xmlify(dataset->pspfile, b, false));
  }
  fprintf(outfile, "</psp_file>\n");

  // The model priors file.
  fprintf(outfile, "<priors_file>");
  if (dataset->plib_name) {
    fprintf(outfile, "%s", xmlify(dataset->plib_name, b, false));
  }
  fprintf(outfile, "</priors_file>\n");

  // Replace any newlines in stopping_reason with spaces.
  char *t = stopping_reason;
  while ((t = strchr(t, '\n')) != NULL) { *t = ' '; };
  fprintf(
    outfile,
    "<reason_for_stopping>%s</reason_for_stopping>\n",
    stopping_reason
  );
  char *bfile = dataset->bfile ? dataset->bfile : "--sequences--";
  fprintf(outfile, "<background_frequencies source=\"%s\" order=\"%d\">\n", 
      xmlify(bfile, b, true),
      dataset->back_order);
  fprintf(outfile, "<alphabet_array>\n");
  int i;
  for (i = 0; i < alph_size_core(dataset->alph); i++) {
    fprintf(
      outfile,
      "<value letter_id=\"%s\">%.3g</value>\n",
      alph_xml_id(dataset->alph, i, b),
      get_array_item(i, dataset->back)
    );
  }
  fprintf(outfile, "</alphabet_array>\n");
  fprintf(outfile, "</background_frequencies>\n");
  fprintf(outfile, "</model>\n");
  str_destroy(b, false);
}

/**********************************************************************
  print_meme_motifs_xml

  Print XML elements for the motifs. See DTD embeded in
  print_meme_header_xml for the description of model XML syntax.

***********************************************************************/
static void print_meme_motifs_xml(
  MODEL *model,                   // the model IN 
  DATASET *dataset,               // the dataset IN 
  int nmotifs,                    // number of motifs IN 
  MOTIF_SUMMARY *motif_summaries, // List of final motif properties IN 
  FILE* outfile                   // output file IN 
) {
  fprintf(outfile, "<motifs>\n");
  int i = 0;
  for (i = 0; i < nmotifs; i++) {
    fprintf(
      outfile,
      //"<motif id=\"motif_%d\" name=\"%d\" width=\"%d\" sites=\"%d\""
      "<motif id=\"motif_%d\" name=\"%s\" alt=\"MEME-%d\" width=\"%d\" sites=\"%d\""
      " ic=\"%.1f\" re=\"%.1f\""
      " llr=\"%.0f\" p_value=\"%3.1fe%+04.0f\" e_value=\"%3.1fe%+04.0f\" bayes_threshold=\"%g\""
      " elapsed_time=\"%f\">\n",
      i + 1,
      motif_summaries[i].consensus,	// TLB: new name for motifs
      i + 1,				// TLB: include alternate name 'MEME-i'
      motif_summaries[i].width,
      motif_summaries[i].num_sites,
      motif_summaries[i].ic,
      motif_summaries[i].re,
      motif_summaries[i].llr,
      motif_summaries[i].p_value_mant,
      motif_summaries[i].p_value_exp,
      motif_summaries[i].e_value_mant,
      motif_summaries[i].e_value_exp,
      motif_summaries[i].bayes,
      motif_summaries[i].elapsed_time
    );
    print_meme_pssm_xml(
      motif_summaries[i].pssm,
      dataset->alph,
      motif_summaries[i].width,
      outfile
    );
    print_meme_psfm_xml(
      motif_summaries[i].psfm,
      dataset->alph,
      motif_summaries[i].width,
      outfile
    );
    print_meme_regular_expression_xml(
      motif_summaries[i].regexp,
      outfile
    );
    print_meme_contributing_sites_xml(
      model,
      &(motif_summaries[i]),
      dataset,
      outfile
    );
    fprintf(outfile, "</motif>\n");
  }
  fprintf(outfile, "</motifs>\n");
}

/**********************************************************************
  print_meme_pssm_xml

  Print XML elements for the log-odds matrix. See DTD embeded in
  print_meme_header_xml for the description of model XML syntax.
**********************************************************************/
static void print_meme_pssm_xml(
  double **logodds,// pointer to matrix of log-odds scores IN 
  ALPH_T *alph,    // alphabet IN
  int width,       // width of the motif IN 
  FILE* outfile    // pointer to output file IN 
) {
  int i = 0;
  int j = 0;
  STR_T *b = str_create(3);

  fprintf(outfile, "<scores>\n<alphabet_matrix>\n");
  for (i = 0; i < width; i++) {   // site position 
    fprintf(outfile, "<alphabet_array>\n");
    for (j = 0; j < alph_size_core(alph); j++) { // letter 
      fprintf(outfile,
        "<value letter_id=\"%s\">%d</value>\n",
        alph_xml_id(alph, j, b),
        NINT(100*logodds(i,j))
      );
    }
    fprintf(outfile, "</alphabet_array>\n");
  }
  fprintf(outfile, "</alphabet_matrix>\n</scores>\n");
  str_destroy(b, false);
}

/**********************************************************************
  print_meme_psfm_xml

  Print XML elements describing the frequency matrix. See DTD embeded in
  print_meme_header_xml for the description of model XML syntax.
**********************************************************************/
static void print_meme_psfm_xml(
  THETA theta,  // pointer to matrix of frequencies 
  ALPH_T *alph, // alphabet
  int width,    // width of the motif 
  FILE* outfile // pointer to output file 
) {
  int i = 0;
  int j = 0;
  STR_T *b;

  b = str_create(3);
  fprintf(outfile, "<probabilities>\n<alphabet_matrix>\n");
  for (i=0; i < width; i++) { // site position 
    fprintf(outfile, "<alphabet_array>\n");
    for (j=0; j < alph_size_core(alph); j++) { // letter 
      fprintf(
        outfile,
        "<value letter_id=\"%s\">%f</value>\n",
        alph_xml_id(alph, j, b),
        theta(i, j)
      );
    }
    fprintf(outfile, "</alphabet_array>\n");
  }
  fprintf(outfile, "</alphabet_matrix>\n</probabilities>\n");
  str_destroy(b, false);
}

/**********************************************************************/
/*
  print_meme_regular_expression_xml

  Prints the XML element describing a motif's  regular expression.
  See DTD embeded in meme-dtd.h for the description of model
  XML syntax.
*/
/**********************************************************************/
static void print_meme_regular_expression_xml(
  char* regexp, // regular expression  
  FILE *outfile
) {
  STR_T *buffer;
  buffer = str_create(10);
  fprintf(
    outfile,
    "<regular_expression>\n%s\n</regular_expression>\n",
    xmlify(regexp, buffer, false)
  );
  str_destroy(buffer, false);
}

/**********************************************************************/
/*
  print_meme_contributing_sites_xml

  Print XML elements describing a motif's occurences. See DTD embeded in
  dtd.h for the description of model XML syntax.
*/
/**********************************************************************/
static void print_meme_contributing_sites_xml(
  MODEL *model,
  MOTIF_SUMMARY *motif_summary,
  DATASET *dataset,
  FILE* outfile
) {
  char *seq, *strand, site[MAXSITE+1], pre[10+1], post[10+1];
  int i, j, seqno, lseq, motif_start;
  bool ic;
  SAMPLE *s;
  ALPH_T *alph;
  STR_T *b;
  
  alph = dataset->alph;
  b = str_create(10);
  fprintf(outfile, "<contributing_sites>\n");
  for (i = 0; i < motif_summary->num_sites && (dataset->n_samples <= dataset->brief); i++) {
    seqno = motif_summary->sites[i].x;
    s = dataset->input_order[seqno];
    lseq = s->length; // length of sequence 
    seq = s->seq; // the ASCII sequence 
    ic = motif_summary->sites[i].ic;
    motif_start = motif_summary->sites[i].y;
    strand = "none";
    if (alph_has_complement(alph)) { 
      strand = ic ? "minus" : "plus";
    } 
    // get the aligned sequence parts 
    get_aligned_sequence_parts(alph, motif_start, motif_summary->width,
        ic, lseq, seq, pre, site, post);
    fprintf(outfile,
      "<contributing_site sequence_id=\"sequence_%d\" position=\"%d\""
      " strand=\"%s\" pvalue=\"%.2e\" >\n",
      seqno, motif_start, strand, motif_summary->sites[i].prob);
    fprintf(outfile, "<left_flank>%s</left_flank>\n<site>\n", xmlify(pre, b, false));
    for (j = 0; j < motif_summary->width; j++) {
      fprintf(outfile, "<letter_ref letter_id=\"%s\"/>\n",
          alph_xml_id(alph, alph_index(alph, site[j]), b));
    }
    fprintf(
      outfile,
      "</site>\n<right_flank>%s</right_flank>\n</contributing_site>\n",
      xmlify(post, b, false)
    );
  }
  fprintf(outfile, "</contributing_sites>\n");
  str_destroy(b, false);
}
/**********************************************************************/
/*
get_aligned_sequence_parts

Extract the sequence corresponding to the motif site and 10 bases
  on the left and right flanks.

*/
/**********************************************************************/
void get_aligned_sequence_parts(
  ALPH_T *alph,
  int motif_start,
  int motif_width,
  bool ic,
  int lseq,
  char *seq,
  char *pre,
  char *site,
  char *post
) {
  int i = 0;
  int ii = 0;
  if (!ic) { // + strand 
    // pre 
    for (i = motif_start - 10, ii = 0; i < motif_start; i++) {
      if (i<0) {
        continue;
      }
      pre[ii++] = seq[i];
    }
    pre[ii] = '\0';
    // site 
    for (i = motif_start, ii = 0; ii < motif_width; i++) {
      site[ii++] = seq[i];
    }
    site[ii] = '\0';
    // post 
    for (i = motif_start + motif_width, ii = 0; ii < 10 && i < lseq; i++) {
      post[ii++] = seq[i];
    }
    post[ii] = 0;
  }
  else { // - strand 
    // pre 
    for (i = motif_start + motif_width + 9, ii = 0; i >= motif_start + motif_width; i--) {
      if (i>=lseq) {
        continue;
      }
      pre[ii++] = comp_sym(alph, seq[i]);
    }
    pre[ii] = '\0';
    // site 
    for (i = motif_start + motif_width - 1, ii = 0; ii < motif_width; i--) {
      site[ii++] = comp_sym(alph, seq[i]);
    }
    site[ii] = '\0';
    // post 
    for (i = motif_start - 1, ii = 0; ii < 10 && i >= 0; i--) {
      post[ii++] = comp_sym(alph, seq[i]);
    }
    post[ii] = '\0';
  } // strand 
}

/**********************************************************************/
/*
  score_sites

  Score and get the pvalues of the sites in the model.
  Sort in order of increasing p-value.
*/
/**********************************************************************/
static void score_sites(
  DATASET *dataset,     // the dataset 
  MODEL *model,       // the model 
  LO *lo,       // LO structure 
  double *pv        // p-values for scores of this motif 
)
{
  int isite;
  P_PROB sites = model->maxima;   // sites "defining" model 
  int n = model->nsites_dis;    // number of sites 
  STYPE stype = (alph_has_complement(dataset->alph) ? (model->invcomp ? Combine : Norc) : Unstranded);
  SCORE **scores = NULL;    // the site scores 
  int old_seqno = -1;
  XLATE_T *xlate = NULL; // no translation;

  for (isite = 0; isite < n; isite++) {     // site 
    int seqno = sites[isite].x;             // sequence number 
    int y = sites[isite].y;                 // location of site 
    // sequence
    SAMPLE *s = dataset->input_order[seqno];
    int lseq = s->length;                   // length of sequence 
    char *seq = s->seq;                     // the ASCII sequence 
    double pvalue;                          // score p-value 

    // score the sequence if new 
    if (old_seqno != seqno) {
      if (old_seqno >= 0) free_2array(scores, 1);
      scores = score_sequence(dataset->alph, xlate, stype, false, seq, lseq, 1, &lo);
      old_seqno = seqno;
      s->minpv = 1.0;
    }

    pvalue = pv[(int) scores[0][y].score];  // p-value 

    // save MINUS the p-value in the .prob field of sites 
    sites[isite].prob = -pvalue;

    // update minimum p-value of sites 
    if (pvalue < s->minpv) s->minpv = pvalue;

  } // get p-values of sites 
  free_2array(scores, 1);               // free space 

  // sort the sites by p-value
  qsort((p_prob*)sites, n, sizeof(p_prob), pY_compare);

  // change sign of p-values back
  for (isite = 0; isite < n; isite++) sites[isite].prob *= -1;

} // score_sites 

/**********************************************************************/
/*
  print_site_diagrams

  Make block diagrams of the actual sites in the model
  and print them.
  Sequences are sorted by the minimum p-value of sites in them.
*/
/**********************************************************************/
static void print_site_diagrams(
  DATASET *dataset,   // the dataset 
  MODEL *model,       // the model 
  int nmotifs,        // number of motifs in los 
  char *consensus,    // single-letter consensus for motif
  FILE *outfile       // where to print 
)
{
  int i, j, isite;
  P_PROB sites = model->maxima;   // sites "defining" model 
  int n = model->nsites_dis;      // number of sites 
  int nseqs = dataset->objfun!=NZ ? dataset->n_samples : dataset->n_group[0]; // number of sequences in dataset 
  //bool dna = dataset->dna;       // dataset is DNA if true 
  // bool invcomp = model->invcomp; // use reverse complement strand, too 
  XLATE_T *xlate = NULL;       // not translating;
  bool best_motifs = false;      // use all sites 
  double m_thresh = 1;           // show all sites as strong 
  STYPE stype = (alph_has_complement(dataset->alph) ? (model->invcomp ? Combine : Norc) : Unstranded);
  int nseqs_with_sites;     // number of sequences with sites 
  int *seqlist = NULL;      // ordered list of sequences w/sites 
  int *hits = NULL;     // store hits 
  double *pvalues = NULL;   // store pvalues 
  char *f = "%-*s%s %8s  %s\n";   // format 

  // Create the list to contain sequences sorted by minimum p-value
  Resize(seqlist, nseqs, int);

  // Clear list of sites for each sequence
  for (i = 0; i < nseqs; i++) dataset->input_order[i]->nsites = 0;

  // Find which sequences have sites and create list of sites for each
  for (isite = 0, nseqs_with_sites = 0; isite < n; isite++) { // site 
    int seqno = sites[isite].x;   // sequence number 
    int y = sites[isite].y + 1;   // location of site (plus 1) 
    bool ic = sites[isite].ic;   // site on reverse complement strand 
    SAMPLE *s = dataset->input_order[seqno];// sequence 

    // record the sequence as containing a site 
    if (!s->nsites) seqlist[nseqs_with_sites++] = seqno;

    // store the site in its list of sites 
    Resize(s->sites, s->nsites+1, int);
    s->sites[(s->nsites)++] = ic ? -y : y;  // +/- site offset by 1 
  } // site 

  for (i = 0; i < PAGEWIDTH; i++) { fputc('-', outfile); } fputc('\n', outfile);
  fprintf(outfile, "\tMotif %s MEME-%d block diagrams\n", consensus, model->imotif);
  for (i = 0; i < PAGEWIDTH; i++) { fputc('-', outfile); } fputc('\n', outfile);
  fprintf(outfile, f, MSN, "SEQUENCE NAME", "", "POSITION P-VALUE", "MOTIF DIAGRAM");
  fprintf(outfile, f, MSN, "-------------", "", "----------------", "-------------");

  // create and print a block diagram for each sequence with sites
  for (i = 0; i < nseqs_with_sites; i++) {  // sequence 
    int seqno = seqlist[i];                 // current sequence 
    SAMPLE *s = dataset->input_order[seqno];// sequence 
    int lseq = s->length;                  // length of sequence 
    char *name = s->sample_name;           // name of sequence 
    double minpv = s->minpv;               // minimum p-value of sites 
    char hdr[MSN+20];                      // header 
    TILING tiling;                         // tiling struct 
    tiling.diagram = NULL;	           // silence compiler warnings

    // create storage for hits and pvalues and clear them
    Resize(hits, lseq, int);
    Resize(pvalues, lseq, double);
    for (j=0; j<lseq; j++) { hits[j] = 0; pvalues[j] = 0; }

    // copy hits from s->nsites into hits array 
    for (j=0; j<s->nsites; j++) {
      int y = abs(s->sites[j]) - 1; // position of site 
      int m = (s->sites[j] > 0) ? los[nmotifs-1]->imotif : -los[nmotifs-1]->imotif;
      hits[y] = m;      // +/- motif 
    }

    // put the hits in TILING struct 
    tiling.diagram = NULL;	// to prevent compiler warning
    tiling.hits = hits;
    tiling.pvalues = pvalues;

    // create the block diagram 
    tiling.diagram = create_diagram(xlate, stype, best_motifs, false,
      m_thresh, nmotifs, los, lseq, false, NULL, 0, 0, tiling);

    // print the diagram 
    sprintf(hdr, "%-*.*s %16.2g  ", MSN, MSN, name, minpv);
    print_diagram(tiling.diagram, hdr, outfile);

    myfree(tiling.diagram);   // release space 
  } // sequence 

  // print a final line of hyphens 
  for (i = 0; i < PAGEWIDTH; i++) { fputc('-', outfile); } fprintf(outfile, "\n\n");

  // cleanup
  free(hits);
  free(pvalues);

  myfree(seqlist);

} // print_site_diagrams 

/**********************************************************************/
/*
  align_sites

  Align all sites that make up the model.
*/
/**********************************************************************/
static void align_sites(
  DATASET *dataset,     // the dataset 
  MODEL *model,         // the model 
  char *consensus,      // single-letter consensus
  LO *lo,               // LO structure 
  double *pv,           // pvalues for scores of this motif 
  FILE *outfile         // stream for output 
)
{
  int i, ii, isite;
  int w = model->w;               // length of site 
  P_PROB sites = model->maxima;   // sites "defining" model 
  int n = model->nsites_dis;      // number of sites 
  bool invcomp = model->invcomp;  // use reverse complement strand, too 
  int imotif = lo->imotif;        // name of motif 
  char site[MAXSITE+1], pre[10+1], post[10+1];

  // print header 
  for (i=0; i<PAGEWIDTH; i++) { fputc('-', outfile); } fputc('\n', outfile);
  fprintf(outfile,
    "\tMotif %s MEME-%d sites sorted by position p-value\n", consensus, model->imotif);
  for (i=0; i<PAGEWIDTH; i++) { fputc('-', outfile); } fputc('\n', outfile);
  fprintf(outfile, "%-*.*s%s ", MSN, MSN, "Sequence name",
    invcomp ? " Strand" : "");
  fprintf(outfile, "%6s %9s %10s %*sSite%*s\n",
    "Start", "P-value", "", w/2 - 2, "", w - w/2 - 4, "");
  fprintf(outfile, "%-*.*s%s ", MSN, MSN, "-------------",
    invcomp ? " ------" : "");
  fprintf(outfile, "%6s %8s %10s ", "-----", "---------", "");
  for (i=0; i<w; i++) fputc('-', outfile);
  fputc('\n', outfile);

  // print sites that make up the model
  for (isite=0; isite<n; isite++) {         // site 
    int seqno = sites[isite].x;             // sequence number 
    int y = sites[isite].y;                 // location of site 
    bool ic = sites[isite].ic;              // strand direction 
    double pvalue = sites[isite].prob;      // position p-value 
    SAMPLE *s = dataset->input_order[seqno];// sequence 
    int lseq = s->length;                   // length of sequence 
    char *seq = s->seq;                     // the ASCII sequence 
    char *sample_name = s->sample_name;     // name of sample 

    // print name and strand 
    fprintf(outfile, "%-*.*s%s ", MSN, MSN, sample_name,
      invcomp ? (ic ? "     -" : "     +") : "");

    // print position and position p-value 
    fprintf(outfile, "%6d %9.2e", y+1, pvalue);

    // get the aligned sequence parts 
    if (!ic) {        // + strand 
      // pre 
      for (i=y-10, ii=0; i<y; i++) {
        if (i<0) continue;
        pre[ii++] = seq[i];
      }
      pre[ii] = '\0';
      // site 
      for (i=y, ii=0; ii<w; i++)  site[ii++] = seq[i];
      site[ii] = '\0';
      // post 
      for (i=y+w, ii=0; ii<10 && i<lseq; i++) post[ii++] = seq[i];
      post[ii] = 0;

    } else {        // - strand 
      // pre 
      for (i=y+w+9, ii=0; i>=y+w; i--) {
        if (i>=lseq) continue;
        pre[ii++] = comp_sym(dataset->alph,seq[i]);
      }
      pre[ii] = '\0';
      // site 
      for (i=y+w-1, ii=0; ii<w; i--) site[ii++] = comp_sym(dataset->alph, seq[i]);
      site[ii] = '\0';
      // post 
      for (i=y-1, ii=0; ii<10 && i>=0; i--) post[ii++] = comp_sym(dataset->alph, seq[i]);
      post[ii] = '\0';
    } // strand 

    // print the alignment 
    if (pre[0] == '\0') {     // print a dot in empty pre 
      fprintf(outfile, " %10s %-*s %-10s\n", ".", w, site, post);
    } else {
      fprintf(outfile, " %10s %-*s %-10s\n", pre, w, site, post);
    }

  } // site 

  // print line of hyphens 
  for (i=0; i<PAGEWIDTH; i++) { fputc('-', outfile); } fprintf(outfile, "\n\n");

} // align_sites 

/*
  Get a single-letter consensus from a model.
*/
extern char *get_single_letter_consensus(
  MODEL *model,
  ALPH_T *alphabet
) {
  // Create a motif from the model.
  MATRIX_T *freqs = convert_matrix(model->theta, model->w, model->alength);
  MOTIF_T *motif = allocate_motif("", "", alphabet, freqs, NULL);

  // Get a string representing the motif consensus.
  STR_T *id_buf = str_create(10);
  str_clear(id_buf);
  motif2consensus(motif, id_buf, true);
  char *motif_id = str_internal(id_buf);
  str_destroy(id_buf, true);

  // Free space.
  free_matrix(freqs);
  destroy_motif(motif);

  return(motif_id);
} // get_single_letter_consensus
