/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1995-2015, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/

 // 5-26-00 tlb; add cons0 to copy_model 
 // 5-23-00 tlb; fix finding of sites for OOPS to reflect scaling of zij 
 // 7-29-99 tlb; add nsites and nsites_obs to copy_model 
 // 7-27-99 tlb; add rentropy to MODEL; rem. p_point, nmotifs from create_model
 // 7-14-99 tlb; move get_sites and make_log_odds here from display.c 
 // 7-02-99 tlb; add logtheta to model, lambda_obs to copy_model 
 // 6-28-99 tlb; remove sigma from model, add psites 

#include "calculate_p_y.h"
#include "meme.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "seed.h"
#include "psp.h"

static int get_n_diffs_helper(
  char *str1, // A string in the aligned pair
  char *str2, // A string in the aligned pair
  int offset // Number of characters (>=0) str1 is shifted to the right of str2.
);


/**********************************************************************/
/*
  copy_theta

  Copy the theta matrix to another array.

*/
/**********************************************************************/
void copy_theta(
  THETA s, // source 
  THETA d, // destination 
  int w, // width of motif 
  int alength // length of alphabet 
)
{ 
  int i, j;
  for (i = 0; i < w; i++) { // col in motif 
    for (j = 0; j < alength; j++) { // row in motif 
        theta_ref(d, i, j) = theta_ref(s, i, j);
    }
  }
}


/**********************************************************************/
/*
        get_not_o

        Calculate the probability that each possible site
        start does not overlap a previously found site;
        this is taken to be the minimum of all the probabilities
        of a site starting within the possible site.

        Sets not_o and log_not_o in dataset.

*/
/**********************************************************************/
void get_not_o(
  DATASET *dataset, // the dataset
  int w // motif width
)
{
  int i,j,k;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;

  for (i=0; i<n_samples; i++){ // sequence 
    SAMPLE *s = samples[i];
    int lseq = s->length;
    WEIGHTS_T *weights = s->weights; // prb not in a previous site 
    double *not_o = s->not_o; // prb not overlapping a site 
    int *log_not_o = s->log_not_o; // log prb not overlapping a site 

    if (lseq < w) continue; // sequence to short for motif 
    skip_group_if_required(dataset, s, i);	// ignore holdout groups

    for (j=0; j <= lseq - w; j++) { // site start 
      not_o[j] = 1.0; // assume not overlapping 
      for (k=j; k < j+w; k++) { // position in sequence 
        if (weights[k] < not_o[j]) not_o[j] = weights[k];
      }
      log_not_o[j] = INT_LOG(not_o[j]);
    } // for j 

    for (j=lseq-w+1; j < lseq; j++) { // beyond possible site starts 
      not_o[j] = 1;
      log_not_o[j] = 0;
    }

  } // for i 
} // get_not_o 

/**********************************************************************/
/*
        create_model

        Create a model structure.
*/
/**********************************************************************/
MODEL *create_model(
  MOTYPE mtype, // type of model 
  bool invcomp, // use inv comp strand  
  int max_w, // maximum width of motif 
  ALPH_T *alph, // alphabet 
  OBJTYPE objfun // MEME objective function
)
{
  MODEL *model = (MODEL *) mm_malloc(sizeof(MODEL));

  model->mtype = mtype;
  model->invcomp = invcomp;
  model->alength = alph_size_core(alph);

 // add one extra column in case model is a palindrome 
  create_2array(model->theta, double, max_w+1, alph_size_wild(alph));
  create_2array(model->logtheta, double, max_w+1, alph_size_wild(alph));
  create_2array(model->logtheta_rc, double, max_w+1, alph_size_wild(alph));
  create_2array(model->obs, double, max_w+1, alph_size_core(alph));
  model->maxima = NULL;
  model->w = 0;
  model->max_w = max_w;
  model->cons[0] = 0;
  model->cons0[0] = 0;
  model->nsites_dis = MINSITES;

  model->iter = 0;
  model->objfun = objfun;

  return model;
}

/**********************************************************************/
/*
  copy_model

  Copy a model structure.
*/
/**********************************************************************/
void copy_model(
  MODEL *m1, // source 
  MODEL *m2, // destination 
  ALPH_T *alph// alphabet
)
{
  int i;

  m2->all_widths = m1->all_widths; // FIXME jj - this was left uncopied
  m2->mtype = m1->mtype;
  m2->min_w = m1->min_w;
  m2->max_w = m1->max_w;
  m2->pw = m1->pw;
  m2->psites = m1->psites;
  if (m1->maxima) { // copy maxima if they exist 
    Resize(m2->maxima, m1->nsites_dis+1, p_prob);
    bcopy((char *) m1->maxima, (char *) m2->maxima, 
      m1->nsites_dis*sizeof(p_prob));
  }
  m2->pal = m1->pal;
  m2->invcomp = m1->invcomp;
  m2->imotif = m1->imotif; 
  m2->w = m1->w;
  copy_theta(m1->theta, m2->theta, m1->w, alph_size_wild(alph));
  copy_theta(m1->obs, m2->obs, m1->w, alph_size_wild(alph));
  m2->lambda = m1->lambda;
  m2->lambda_obs = m1->lambda_obs;
  m2->nsites = m1->nsites;
  m2->nsites_obs = m1->nsites_obs;
  m2->nsites_dis = m1->nsites_dis;
  strcpy(m2->cons, m1->cons); 
  strcpy(m2->cons0, m1->cons0); 
  for (i=0; i<m1->w; i++) {
    m2->rentropy[i] = m1->rentropy[i];
  }
  m2->rel = m1->rel;
  m2->ic = m1->ic;
  m2->ll = m1->ll;
  m2->mll_0 = m1->mll_0;
  m2->mll_1 = m1->mll_1;
  m2->logpv = m1->logpv;
  m2->logev = m1->logev; 
  m2->llr = m1->llr;
  m2->site_threshold = m1->site_threshold;
  m2->ID = m1->ID; 
  m2->iter = m1->iter; 
  m2->iseq = m1->iseq;
  m2->ioff = m1->ioff;
  m2->objfun = m1->objfun;
} // copy_model 

/**********************************************************************/
/*
  free_model

  Free a model structure.
*/
/**********************************************************************/
void free_model(MODEL *model) {
  if (!model) return;
  free_2array(model->theta, model->max_w + 1);
  free_2array(model->logtheta, model->max_w+1);
  free_2array(model->logtheta_rc, model->max_w+1);
  free_2array(model->obs, model->max_w+1);
  free(model->maxima);
  free(model);
}

/**********************************************************************/
/*
        get_sites

        Get the principal sites contributing to a motif.

        OOPS            position with highest z_ij > 0.0 in each sequence
        ZOOPS           position with highest z_ij > 0.5 in each sequence
        TCM             all positions with z_ij > 0.5

        Assumes that z_ij is set for current motif.
        Returns a list consisting of sites and sets n, the number of sites.
*/
/**********************************************************************/
SITE *get_sites(
  DATASET *dataset, // the dataset 
  MODEL *model, // the model 
  int *n, // number of sites found 
  int *best_site // index of best site in array 
)
{
  int i, j;
  int nsites = 0; // number of sites 
  SITE *sites = NULL; // list of sites 
  int n_samples = dataset->n_samples; // number of sequences 
  SAMPLE **samples = dataset->samples; // the sequences 
  MOTYPE mtype = model->mtype; // type of model 
  int w = model->w; // width of motif 
  bool invcomp = model->invcomp; // using invcomp strand 
  double best_z = -1; // best overall zij 

  if (mtype==Oops || mtype==Zoops) { // (Z)OOPS model 

    for (i=0; i<n_samples; i++) { // sequence 
      SAMPLE *s = samples[i]; // sample 
      int lseq = s->length; // length of sequence 
      double max_zij = -1; // flag no z_ij found 
      int best_j = 0; // position of site (1..m) 
      Z_T *zi = s->z;                // z_i[j], j in [-lseq...+lseq]

      if (lseq < w) continue; // sequence to short for motif 

      // find maximum z_ij 
      for (j=0; j<lseq-w+1; j++) { // position 
        int k = j+1;                    // Z_i = k
        double zij = invcomp ? MIN(1.0,Zi(-k)+Zi(k)) : Zi(k);
        if (zij > max_zij) { // bigger found 
          max_zij = zij;
          best_j = j;
        }
      } // Z_i = j

      // record site 
      if ((max_zij && mtype == Oops) || max_zij > 0.5) { // site found 
        if (nsites % 100 == 0) Resize(sites, nsites+101, SITE);
        sites[nsites].seqno = i;
        sites[nsites].pos = best_j;
        sites[nsites].zij = max_zij;
        int best_k = best_j + 1;                // Z_i = best_k
        // on reverse complement strand?
        sites[nsites].invcomp = (invcomp & (Zi(-best_k)>Zi(best_k))); 
        if (max_zij > best_z) {
          *best_site = nsites;
          best_z = max_zij;
        }
        nsites++;
      } // site found 

    } // sequence 

  } else { // TCM model 

    for (i=0; i<n_samples; i++) { // sequence 
      SAMPLE *s = samples[i]; // sample 
      int lseq = s->length; // length of sequence 
      Z_T *zi = s->z;                // z_i[j], j in [-lseq...+lseq]

      if (lseq < w) continue; // sequence to short for motif 

      // find all z_ij > 0.5 
      for (j=0; j<lseq-w+1; j++) { // position 
        int k = j+1;                    // Z_i = k
        double zij = invcomp ? MIN(1.0,Zi(-k)+Zi(k)) : Zi(k); // z_ij 
        if (zij > 0.5) { // record site 
          if (nsites % 100 == 0) Resize(sites, nsites+101, SITE);
          sites[nsites].seqno = i;
          sites[nsites].pos = j;
          sites[nsites].zij = zij;
          // on reverse complement strand? 
          sites[nsites].invcomp = (invcomp && (Zi(-k)>Zi(k)));
          if (zij > best_z) {
            *best_site = nsites;
            best_z = zij;
          }
          nsites++;
        } // record site 
      } // position 
    } // sequence 

  } // TCM model 

 // return results 
  *n = nsites;
  return sites;
} // get_sites 

/**********************************************************************/
/*      
        make_log_odds

        Compute the log-odds matrix from the motif, negative motif and
        background letter frequencies.
                l = log(p/(q*n+(1-q)b))
        where   p = positive model,
                n = negative model,
                b = background model.

        Include 'X' in the logodds matrix.
        X is given a score equal to the weighted average score for the
        column.
*/
/**********************************************************************/
double** make_log_odds(
  THETA theta1, // positive motif theta 
  THETA theta0, // negative motif theta; use 0 if NULL 
  ARRAY_T *back, // background frequencies; use 0 if NULL 
  double q, // 0<=q<=1, mixing parameter 
  int w, // width of motif 
  int alength // length of alphabet 
)
{
  int i, j;
  double **logodds = NULL; // the logodds matrix 

  // compute the log-odds matrix 
  Resize(logodds, w, double *); // create rows 
  for (i=0; i<w; i++) { // site position 
    logodds[i] = NULL; // create columns 
    Resize(logodds[i], alength+1, double); // include 'X' 

    // calculate log-odds for this column 
    logodds(i,alength) = 0; // 'X' 
    for (j=0; j<alength; j++) { // letter 
      double p = theta1(i,j); // positive motif 
      double n; // negative motif 
      if (!theta0) { // no negative motif given 
        n = get_array_item(j, back); // background 
      } else if (!back) { // no background given 
        n = theta0(i,j); // negative motif 
      } else { // negative and background given 
        n = q*theta0(i,j) + (1-q) * get_array_item(j, back); // blend negative & background 
      }
      if (n==0) {
        logodds(i,j) = 0; // letter with zero prob 
      } else {
        logodds(i,j) = LOG2(p/n);
      }
      logodds(i,alength) += n * logodds(i,j); // 'X' 
    } // letter 
  } // column 

  return logodds;
} // make_log_odds 


/**
 * get_pred_sites
 *
 * Retrieve the sites predicted from the motif derived from the specified
 * seed. Do this according to the specified sequence model and number of
 * sites parameter.
 *
 * The calling function must pass in an array P_PROB with enough storage space
 * to store the maximum number of maxima, and the caller must also be
 * responsible for freeing that memory.
 *
 * \return The number of elements in the array of the most probable sites in the
 * dataset
 */
int get_pred_sites (
  P_PROB psites, // An array to contain the predicted sites
  MOTYPE mtype, // Model type
  int w, // Length of seed being evaluated
  char *seed, // ASCII version of seed being evaluated
  int *lmotif[MAXSITE], // Storage space for the motif model
  int lmap[MAXALPH][MAXALPH], // The sequence to theta log map
  DATASET *dataset, // The dataset of sequences
  bool ic // Use inverse complement
) {
  // Convert the seed into a motif using the specified mapping array:
  int seed_len; // Redundant, as w is known anyway
  uint8_t *e_seed = to_e_seed(dataset->alph, seed, &seed_len); // Integer encoding
  init_lmotif(w, e_seed, lmotif, lmap);

  // Evaluate the matches of that motif in the dataset
  // (ie initialise p(Y_ij | theta)):
  if (!ic) {
    get_pY(dataset, lmotif, w, 0);
  } else {
    get_pY(dataset, lmotif, w, 1);
    get_pY(dataset, lmotif, w, 2);
  }

  // Put highest pY into first scratch array if using both DNA strands:
  if (ic) combine_strands(dataset, w);

  // The maxima in pY indicate the sites that would be predicted by the seed:
  int n_maxima = get_max(mtype, dataset, false, w, 0, psites, ic, true);

  return n_maxima;
}


/**
 * print_site_array
 *
 * Prints the specified array of sites to the specified file handle.
 */
void print_site_array(
  P_PROB sites, // An array of sites to be printed
  int nsites, // Length of the array
  FILE *outfile, // The stream for output
  int w, // The size of each of the sites
  DATASET *dataset // Contains the sequences which contain the sites
) {
  // Print out the sites predicted by the seed:
  fprintf(outfile, "###########################\n");
  int site_idx;
  for (site_idx = 0; site_idx < nsites; site_idx++) {
    int seq_num = sites[site_idx].x;
    int site_loc = sites[site_idx].y;
    uint8_t *e_site = (dataset->samples[seq_num]->res)+site_loc;
    char *curr_site = to_str_seed(dataset->alph, e_site, w);
    SAMPLE *seq = dataset->samples[seq_num];
    char *seq_name = seq->sample_name;
    fprintf(stdout, "%7s %4i %s\n", seq_name, site_loc, curr_site);
  }
  fprintf(stdout, "---------------------------\n");
} // print_site_array


/**
 * get_first_siteloc
 *
 * Parse the specified sequence description, to retrieve the integer specifying
 * the location of the first site in the corresponding sequence.
 *
 * /return The location of the first site in the sequence, counting from 0
 */
int get_first_siteloc (
  char *descript // The description of the sequence; contains site info
) {
  // Ascii string containing position of site within sequence. Will be less
  // than 10**9, so string will be 10 or less characters total:
  char pos_tok[10];
  int descript_idx = 0;
  // NOTE: Assuming the correct format; errors in format are not dealt with.
  // This is OK as "get_first_siteloc" is currently only used when the
  // experimental option "planted_LLR" is used.
  // The integer before the first comma indicates the first (and only) site
  // position => Read the description as such:
  while ((descript[descript_idx] != '\0') &&
         (descript[descript_idx] != ',')){
    pos_tok[descript_idx] = descript[descript_idx];
    descript_idx++;
  }
  pos_tok[descript_idx] = '\0';

  if (strlen(pos_tok) > 0) {
    return atoi(pos_tok);
  } else {
    return -1;
  }
}


/**
 * get_n_strdiffs
 *
 * Calculate the number of differences between two strings under an alignment
 * specified by an "offset" value. The offset indicates the number of
 * characters by which the first character of str1 is shifted to the right of
 * str2. These numbers can be negative, indicating str1 is shifted to the left
 * of str2.
 *
 * This is a somewhat obscure function, since it represents the concept of
 * an ungapped alignment (of the two strings) via an offset value. However,
 * it may still be of general utility.
 *
 * \return Total number of differences between str1 and str2, given the
 * alignment specified by "offset". Each overhanging character by either string
 * is counted as a single difference.
 */ 
int get_n_strdiffs (
  char *str1, // A string in the aligned pair
  char *str2, // A string in the aligned pair
  int offset // Number of characters str1 is shifted to the right of str2
) {
  if (offset < 0) {
    return get_n_diffs_helper(str2, str1, (-offset));
  } else {
    return get_n_diffs_helper(str1, str2, offset);
  }
}


/**
 * get_n_diffs_helper
 *
 * Does the bulk of the work for the funcion get_n_strdiffs. Used for
 * simplification.
 *
 * \return Total number of diffs between str1 and str2 given the alignment
 * specified by "offset", which MUST BE POSITIVE.
 */
static int get_n_diffs_helper(
  char *str1, // A string in the aligned pair
  char *str2, // A string in the aligned pair
  int offset // Number of characters (>=0) str1 is shifted to the right of str2.
                     
) {
  assert(offset >= 0);

  int str1_len = strlen(str1);
  int str2_len = strlen(str2);

  // str1 can only be offset to the right of str2, hence the following
  // statement regarding the length of the (ungapped) alignment is true:
  int align_size = MAX(str2_len, (str1_len + offset));

  // Iterate throught the units in the alignment in order to calculate total
  // number of differences:
  int n_diffs = 0;
  int align_idx;
  for (align_idx = 0; align_idx < align_size; align_idx++) {
    char str1_char, str2_char;

    int str1_idx = align_idx - offset; // Current index into str1
    if ((str1_idx < 0) || (str1_idx >= str1_len)) {
      str1_char = '-'; // Current alignment index isn't inside str1 => "gap"
    } else {
      str1_char = str1[str1_idx];
    }

    int str2_idx = align_idx; // For readability's sake.
    if (str2_idx >= str2_len) {
      str2_char = '-'; // Current alignment index isn't inside str2 => "gap"
    } else {
      str2_char = str2[str2_idx];
    }

    if (str1_char != str2_char) {
      n_diffs++;
    }
  } // Indeces in the alignment

  return n_diffs;
}


/**
 * vector_subtract
 *
 * This function instatiates a column of integers with the difference between
 * two columns of integers.
 *  ith value of difference column =
 *  ith value of the first input column - ith value of second input column
 */
void vector_subtract(
  int *col_a, // Positive column
  int *col_b, // Negative column
  int *diff_col, // Column in which to store the differences
  int col_len // Length of all columns
) {
  int col_idx;
  for(col_idx = 0; col_idx < col_len; col_idx++) {
    diff_col[col_idx] = col_a[col_idx] - col_b[col_idx];
  }
}


/**
 * make_geometric_prog
 *
 * Generate a list of INTEGERS in a geometric progression, rounding values
 * to integer (using a cast to int) when they are not integers. The
 * last value in the progression is explicitly specified.
 *
 * \return the array of integers in the geometric progression
 */
int *make_geometric_prog (
  int min_val, // Minimum integer value in the geometric progression
  int max_val, // Maximum integer value in the geometric progression
  double factor, // Factor specifying rate of increase in progression
  int *n_vals // The final number of values in the progression - OUT
) {
  int curr_val;
  int *progression = NULL;
  *n_vals = 0;
  for (curr_val = min_val; curr_val < max_val;
       curr_val = (int)(curr_val*factor)) {
    (*n_vals)++; // Adding an extra value to the progression.
    Resize(progression, *n_vals, int);
    progression[*n_vals - 1] = curr_val;
  }
  // Add the terminal value:
  (*n_vals)++;
  Resize(progression, *n_vals, int);
  progression[*n_vals - 1] = max_val;

  return progression;
}

/**********************************************************************/
/*
        erase

        For all models:
          Reset the weights of the letters to probabilisticaly "erase"
          the ones which occur in sites already found.

*/
/**********************************************************************/
void erase(
  DATASET *dataset, // the dataset 
  MODEL *model // the model 
)
{
  int i, j, k;
  int n_samples = dataset->n_samples; // number of sequences 
  SAMPLE **samples = dataset->samples; // the sequences 
  int w = model->w; // width of motif 
  bool revcomp = model->invcomp;             // using both strands

  /*
    Set z from the maxima stored in the learned model.
  */
  set_z(model, dataset);

  /*
     z_ij is taken as the probability of a site occurring at i,j.
     The probability of a position being in a site is taken
     as the maximum of the z_ij for sites containing (overlapping) it.
     w_ij is set to 1-max(z_ij) times its previous value which
     reflects the independence assumption among motifs.
  */
  for (i=0; i<n_samples; i++)           { // sequence 
    WEIGHTS_T *weights = samples[i]->weights; // w_ij 
    int lseq = samples[i]->length; // seq length 
    Z_T *zi = samples[i]->z;                 // zi[j], j in [-lseq...+lseq]

    if (lseq < w) continue; // sample too short for motif 

    for (j=0; j<lseq; j++) { // position 
      double max_z = 0.0;
      // find largest probability that site overlaps this position 
      for (k=MAX(0,j-w+1); k<=j && k<lseq-w+1; k++) {
        int kk = k+1;                           // |Z_i| = kk
        double z = revcomp ? MIN(1.0,Zi(-kk)+Zi(kk)) : Zi(kk);
        max_z = MAX(max_z, z);
      }
      max_z = MIN(1.0, max_z); // fix roundoff errors 
      // update the probability that position not in a site 
      weights[j] *= 1.0 - max_z;
    }
  }

  if (PRINT_W) print_wij(dataset);
} // erase 

/**********************************************************************/
/*
  create_priors
*/
/**********************************************************************/
PRIORS *create_priors(
  PTYPE ptype, // type of prior to use
  // beta for dirichlet priors;
  // < 0 only returns alphabet
  double beta,
  DATASET *dataset, // the dataset
  char *plib_name // name of prior library
)
{
  int i;
  PRIORS *priors;
  priors = (PRIORS *) mm_malloc(sizeof(PRIORS));
  memset(priors, 0, sizeof(PRIORS));
  priors->ptype = ptype;

  // set up the prior counts
  switch (ptype) {
    case Addone: // add one prior
      for (i = 0; i < alph_size_core(dataset->alph); i++) {
        priors->prior_count[i] = 1.0;
      }
      break;
    case Dirichlet: // simple dirichlet prior
      for (i = 0; i < alph_size_core(dataset->alph); i++) {
        priors->prior_count[i] = beta * get_array_item(i, dataset->back);
      }
      break;
    case Dmix: // mixture of dirichlet's
    case Mega: // megaprior heuristic
    case MegaP: // mod. megaprior heuristic
    {
      priors->plib = read_PriorLib(plib_name, beta, dataset->alph);
      // get b=0 prior for modified mega prior heuristic
      if (ptype == MegaP || ptype == Mega) { // used adj freq with Mega
        double b = 0;
        priors->plib0 = read_PriorLib(plib_name, b, dataset->alph);
      }
      break;
    }
  }
  return priors;
} // create_priors

/***************************************************************************/
/*
        init_meme_background

        Read in the background Markov model (-bfile given), 
	or estimate a higher-order model from the dataset (-markov_order given), 
	or use the residue frequencies adjusted by add-one prior.

        Precalculate the log cumulative background probabilities.

        The log probability of any substring of any sequence will then
        be accessible via:
                Log_back(sample[i]->logcumback, j, w)
        where j is the start of the length-w substring of sequence i.
*/
/***************************************************************************/
void init_meme_background (
  char *bfile,			// background model file 
  bool rc,			// average reverse comps 
  DATASET *dataset,		// the dataset 
  char *alph_file,		// alphabet file
  ALPHABET_T alphabet_type,	// alphabet type
  int order,			// (maximum) order of Markov model to load or order to create
  char *seqfile,		// name of a sequence file (required)
  bool status			// print status message
)
{
  int i, j;
  SAMPLE **samples = dataset->samples; // the sequences 
  int n_samples = dataset->n_samples; // # of sequences in dataset 
  ARRAY_T *back; // freq. of tuples 

  // Get the background Markov model
  if (bfile) {				// -bfile given
    // read in background model
    // if order==-1, get the highest order available
    int old_order = order;
    back = load_markov_model(dataset->alph, &order, bfile);
    if (old_order != -1 && old_order != order) {
      fprintf(stderr, 
        "WARNING: Background order %d requested but the given file only contains order %d.\n"
        "         Setting -markov_order to %d\n",
	old_order, order, order);
    }
  } else if (order != -1) {		// -markov_order given
    back = get_markov_from_sequences(seqfile, &order, 0.1, dataset->alph, alph_file, alphabet_type, rc);
  } else {			// not -bfile and not -markov_order (old way)
    // use dataset frequencies + 1 as 0-order background model
    double *res_freq, wtr;
    int alen_core;
    res_freq = dataset->res_freq; // weighted residue freq 
    wtr = dataset->wgt_total_res; // weighted non-ambiguous letter count 
    alen_core = alph_size_core(dataset->alph);
    order = 0;
    back = allocate_array(alen_core);
    for (i = 0; i < alen_core; i++) {
      // adjust residue frequencies with add-one prior to avoid 0s
      // and round to avoid system dependencies
      double f;
      RND(((res_freq[i] * wtr) + 1) / (wtr + alen_core), 8, f);
      set_array_item(i, f, back);
    }
  }
  if (rc) average_rc_markov_model(dataset->alph, order, back);
  // add x-tuples to model
  extend_markov_model(dataset->alph, true, SUM_FREQS, back);
  // normalize for each prefix if the model is higher order
  for (i = 0; i < get_array_length(back); i += alph_size_wild(dataset->alph)) {
    normalize_subarray(i, alph_size_core(dataset->alph), 1e-7, back); 
    set_array_item(i + alph_wild(dataset->alph), 1.0, back);
  }
  // store background model
  dataset->back = back;
  dataset->back_order = order;
 // precalculate the log cumulative background probabilities for each seq 
  dataset->log_total_prob = 0;
  for (i = 0; i < n_samples; i++) { // sequence 
    SAMPLE *s = samples[i]; // sequence 
    // compute probabilites 
    dataset->log_total_prob += calculate_log_cumulative_background(
        dataset->alph, true, dataset->back_order, dataset->back, s->seq,
        s->logcumback);
  } // sequence 

  if (! NO_STATUS && status) fprintf(stderr, "BACKGROUND: using background model of order %d\n", order);

} // init_meme_background 

/**********************************************************************/
/*
        init_model

        Initialize a model from a starting point.

        Returns false if starting point was not valid.
*/
/**********************************************************************/
bool init_model(
  S_POINT *s_point, // the starting point 
  MODEL *model, // the model to intialize 
  DATASET *dataset, // the dataset 
  int imotif // motif number 
)
{
  int w0;

  // skip if no good starting point found for w0, nsites0 
  if (s_point->score == LITTLE) { return false; }

  // initialize the new motif 
  strcpy(model->cons0, s_point->cons0);
  if (VERBOSE) { fprintf(stderr, "\nStarting point: %s\n", model->cons0); }
  w0 = model->w = model->pw = s_point->w0;
  init_theta(model->theta, s_point->e_cons0, w0, dataset->map, dataset->alph);

  // initialize lambda 
  model->lambda = MIN(s_point->nsites0/wps(dataset, w0, 1), 1);
  model->pal = dataset->pal;

  // initialize prior estimate of number of sites 
  model->psites = s_point->nsites0;

  if (PRINTALL) {
    printf("component %2d: lambda= %8.6f ps= %8.0f\n",
      1, model->lambda, wps(dataset, w0, 1));
    print_theta(0, "", 2, model->nsites_dis, model->theta, model->w, 0, "",
      dataset, stdout);
  }

  // init motif number 
  model->imotif = imotif;
  model->iseq = s_point->iseq;
  model->ioff = s_point->ioff;

  return true;
} // init_model 

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
