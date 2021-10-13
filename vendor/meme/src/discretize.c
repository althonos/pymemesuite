/*#define DEBUG*/
/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 2000-2015, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/
#include "meme.h"
#include "dpalign.h"
#include "fisher_exact.h"
#include "ranksum_test.h"
#include "binomial.h"
#include "bates.h"

/* rounding stuff */
#define RNDEPS 1e-12

//FIXME:
//bool debug = true;
bool debug = false;

static int ma_adjust(
  P_PROB sites,                                 /* the sites */
  int nsites,                                   /* the number of sites */
  int w,                                        /* width of sites */
  int flank,                                    /* add flank cols on left+rgt */
  int min_w,		 			/* minimum width allowed */
  DATASET *dataset,  				/* the dataset */
  int *off 					/* best offset */
);

// For sorting column indices by column RE.
typedef struct pair {int ind; double re;} PAIR;
static int pairCompare(PAIR *e1, PAIR *e2) {
  // compares the RE of two pairs (used for ascending qsort order)
  return (e1->re < e2->re) ? -1 : (e1->re > e2->re) ? 1 : e1->ind - e2->ind;
}

// Keep around between invocations.
static THETA saved_theta = NULL;

/***********************************************************************/
/*
	print_maxima
*/
/***********************************************************************/
static void print_maxima(
  int n_maxima,
  MODEL *model
)
{
#ifdef PARALLEL
#undef printf
#endif
  int i;
  printf("\n");
  for (i=0; i<n_maxima; i++) {
    int x = model->maxima[i].x;
    int y = model->maxima[i].y;
    bool neg = model->maxima[i].negative;
    int rank = model->maxima[i].rank;
    double prob = model->maxima[i].prob;
    printf("max %d %15.12f x %d y %d set %s rank %d\n", 
      i+1, prob, x, y, neg ? "control" : "primary", rank);
  }
  fflush(stdout);
#ifdef PARALLEL
#define printf if (mpMyID() == 0) printf
#endif
} // print_maxima


/***********************************************************************/
/*
	shuffle_integers

	Returns a shuffled list of the integers from 0 to n-1.
*/
/***********************************************************************/
static int *shuffle_integers(
  int n
)
{
  int i, j;
  // Create the list of integers [0, n-1].
  int *list = NULL;
  Resize(list, n, int);
  for (i=0; i<n; i++) list[i] = i;

  // Fisher-Yates shuffle the list.
  for (i=n-1; i>0; i--) {
    j = drand_mt() * (i + 1);
    int tmp = list[i]; list[i] = list[j]; list[j] = tmp; // swap
  } 
  return(list);
} // shuffle_integers

/***********************************************************************/
/*
	int_cmp
*/
/***********************************************************************/
int int_cmp(const void *a, const void *b) 
{ 
  const int *ia = (const int *)a; // casting pointer types 
  const int *ib = (const int *)b;
  return *ia  - *ib; 
} 

/***********************************************************************/
/*
	get_best_nsites

	Find the value of nsites that maximizes the significance of the
	chosen objective function.

	If using a test, break ties using largest nsites, then largest LLR.

	Returns nsites; sets the column scores, best log p-value,
        best log e-value, its LLR and the best site threshold.
*/
/***********************************************************************/
int get_best_nsites(
  MODEL *model,					/* the model */
  DATASET *dataset,				/* the dataset */
  int min_nsites,				/* minimum sites */
  int max_nsites,				/* maximum sites */
  int w, 	 				/* width of motif; may differ in model */
  int n_pos_maxima,				/* number of primary maxima in model */
  int n_neg_maxima,				/* number of control maxima in model */
  P_PROB maxima,				// sorted maxima
  double *col_scores,				/* column scores */
  double *best_wN,				// best weighted number of sites
  double *best_log_pv,				/* best p-value */
  double *best_log_ev,				/* best E-value */
  double *best_llr,  				/* LLR of best p-value */
  double *best_threshold 			// site prob threshold 
)
{
  int i, j;
  int n_nsites; 				/* number of different nsites */
  int nsites;					/* number of sites */
  int best_nsites;				// best value of nsites
  S_POINT *s_points = NULL;			/* for use with align_top...*/
  OBJTYPE objfun = dataset->objfun;		// objective function
  bool use_ranksum = (dataset->test==mRS);
  bool use_binom = (dataset->test==mBN);

  // Return p-value = 1 if min_nsites too large.
  if (min_nsites > n_pos_maxima) {
    *best_log_pv = *best_log_ev = *best_llr = 0;
    return min_nsites;
  }

  /* initialize the s_points for all nsites in range [min_nsites,max_nsites] 
     scoring only positive maxima where the next ranked site was negative
     since those are at the end of a run of positives and should have
     the best significance. Note: this works for Classic etc. because
     all maxima ranks are set to -1.
  */
  max_nsites = MIN(n_pos_maxima, max_nsites);	
  n_nsites = MAX(1, max_nsites-min_nsites+1);
  Resize(s_points, n_nsites, S_POINT);
  for(nsites=min_nsites, j=0; nsites<=max_nsites; nsites++) {
    int rank = maxima[nsites-1].rank;	// maxima start at index 0
    int next_rank = nsites<n_pos_maxima ? maxima[nsites].rank : -1;
    if (nsites == max_nsites || next_rank != rank+1) {
      s_points[j].nsites0 = nsites;	/* number of sites */
      s_points[j].score = LITTLE;	/* no score yet */
      s_points[j].evaluate = true;	/* Evaluate at every s_point */
      j++;				// n_nsites 
    }
  }
  n_nsites = j;

  // NOTE: If not Classic, I think this call is only needed to get the LLR.
  // FIXME: Don't call align... if called from align...!
  if (objfun == Classic || objfun == NC) {
    /* get the probability that a site starting at position x_ij would
       NOT overlap a previously found motif; used in E_STEP.
    */
    //get_not_o(dataset, model->w);
    get_not_o(dataset, w);

    // TLB: Works better without the following line on Yeast examples.
    //add_psp_to_log_not_o(dataset, model->w, model->invcomp, model->mtype);

    /* 
       align the top nsites sorted subsequences and compute the 
       log_pop or LLR function on each alignment with nsites in [min_nsites,max_nsites]
    */
    bool orig_use_llr = dataset->use_llr;	// What to compute on alignment?
    dataset->use_llr = false;			// Causes log_pop to be used in Classic mode.
    (void) align_top_subsequences(model->mtype, w, dataset, 0, 0, 0, 0, n_nsites, 
      n_pos_maxima, 0, maxima, col_scores, s_points);
    dataset->use_llr = orig_use_llr;		// Restore setting.
  }

  // Set up parameters for CE (not really needed for CD).
  double ce_success_prob = 0;
  if (objfun == CE || objfun == CD) {
    int len = dataset->min_slength;
    ce_success_prob = 2 * dataset->ce_max_dist;	// Number of positions for motif center in central region.
    if (len+w % 2 == 0) ce_success_prob += 1;	// Centers can align perfectly gives an extra position.
    ce_success_prob /= len - w + 1;		// Divide by total possible positions for motif.
  }

  /* 
     Determine the significance of the score for each number of 
     sites and chose the number of sites to minimize it.
  */
  *best_llr = -BIG;
  *best_log_ev = BIG;
  best_nsites = 0;
  // Traverse the starting points in reverse order so we can take advantage of the
  // fact that the p-value will always be better for larger nsites at a given fraction
  // of observed successes (p_obs).  We can therefore skip evaluation if the p_obs
  // is no greater than that previously seen for a larger value of nsites for binomial test and FET.
  int slen = dataset->min_slength;		// All sequences should be same length for CE.
  double p_obs = 0;				// Fraction of successes at a given threshold
  double best_p_obs = -1;			// Best fraction of successes.
  double best_dtc = slen;			// Best distance to center.
  for (i=n_nsites-1; i>=0; i--) {		/* starting points */
    double score = s_points[i].score;		/* Classic: -log_pop; Otherwise: LLR */
    int N = s_points[i].nsites0;
    double wN = maxima[N-1].wN;
    double prob = maxima[N-1].prob;		// site probability score (threshold)
    double log_pv, log_ev;
    if (objfun==Classic || objfun==NC) {	// Classic or NC
      log_ev = get_log_sig(score, model->mtype, w, wN, N, model->invcomp, model->pal, dataset);
      RND(log_ev, RNDDIG, log_ev);
      log_pv = log_ev;				// Not real p-value.
      if (TRACE) {
	double m1, e1;
	exp10_logx(log_ev/log(10.0), m1, e1, 1);
	printf("get_best_nsites: score %f w %d wN %f N %d E-value %3.1fe%+04.0f\n", score, w, wN, N, m1, e1);
      }
    } else if (objfun == CE || objfun == CD) {
      if (objfun == CD) {			// CD
        // Get p-value of average distance less than or equal to that observed.
        // This optimization works because we are working backwards through ranks.
        double dtc = maxima[N-1].dist;
        log_pv = (dtc < best_dtc) ? get_dtc_log_pv(dtc, wN, w, slen) : 0;
      } else {					// CE
	// Binomial test on the number of sites in the central region.
	int n_trials = wN;				// Round down (weighted) rank.
	int n_successes = maxima[N-1].in_region;	// Round down (weighted) number of sites in central region.
	double p_success = ce_success_prob;
        p_obs = (double) n_successes/n_trials;
	log_pv = (p_obs > best_p_obs) ? log_betai(n_successes, n_trials - n_successes + 1, p_success) : 0;
      }
      RND(log_pv, RNDDIG, log_pv);
      log_ev = log_pv;				// Not real E-value.
    } else {					// DE, SE, NZ
      if (model->mtype == Oops || use_ranksum) {
        // OOPS model uses ranksum test
	int rank = maxima[N-1].rank;
        // This optimization works because we are working backwards through ranks.
        double ta_obs = maxima[N-1].ranksum;	// sum of positive ranks
        int n = n_pos_maxima + n_neg_maxima;
        int na = N;				// number of positives
        // Rank-sum test.
        RSR_T *r = ranksum_from_stats(n, na, ta_obs);
        log_pv = RSR_get_log_p_left(r);
	RND(log_pv, RNDDIG, log_pv);
	log_ev = log_pv;			// Not real E-value.
      } else if (use_binom) {
	int n_trials = maxima[N-1].rank;
	int n_successes = N;
        double p_success = ((double) n_pos_maxima)/(n_pos_maxima+n_neg_maxima);
        p_obs = (double) n_successes/n_trials;
	// Binomial test.
        // Compute p-value only if p_obs is higher than for more sites.
	log_pv = (p_obs > best_p_obs) ? log_betai(n_successes, n_trials - n_successes + 1, p_success) : 0;
	RND(log_pv, RNDDIG, log_pv);
	log_ev = log_pv;			// Not real E-value.
      } else {
	// non-OOPS models use mHG unless mRS or mBN requested
	int rank = maxima[N-1].rank;
	int pos = n_pos_maxima;
	int neg = n_neg_maxima;
	int pos_succ = N;
	int neg_succ = rank - N;
        double p_success = ((double) n_pos_maxima)/(n_pos_maxima+n_neg_maxima);
        p_obs = (double) pos_succ/rank;
	// Fisher Exact test.  Save time if it won't beat previous best threshold.
        // This optimization works because we are working backwards through ranks.
	log_pv = (p_obs > best_p_obs) ? getLogFETPvalue(pos_succ, pos, neg_succ, neg, false) : 0;
	RND(log_pv, RNDDIG, log_pv);
	log_ev = log_pv;			// Not real E-value.
      }
    }

    // Update the best observed success fraction.
    if (p_obs > best_p_obs) { best_p_obs = p_obs; } 

    if (TRACE) printf("get_best_nsites: w %d N %d wN %f log_ev %f \n", w, N, wN, log_ev);

    /* Save EV if best so far, breaking ties using LLR then nsites */
    if ( 
      (objfun==Classic 
        && RNDEPS < *best_log_ev - log_ev)
      || (objfun!=Classic
        && (
	  (log_ev < *best_log_ev)
	  || (log_ev == *best_log_ev && score > *best_llr )
	  || (log_ev == *best_log_ev && score == *best_llr && N > best_nsites) 
        )
      )
    ) {
      best_nsites = N;				/* number of sites */
      *best_wN = wN;
      *best_log_pv = log_pv;
      *best_log_ev = log_ev;
      *best_llr = score;
      *best_threshold = prob;
    }

  } /* nsites */

  myfree(s_points);

  return best_nsites;
} /* get_best_nsites */

/***********************************************************************/
/*
	set_z

	Set the z to 1/0 using list of sites (maxima).
*/
/***********************************************************************/
void set_z (
  MODEL *model,					/* the model */
  DATASET *dataset 				/* the dataset */
)
{
  int i, j;
  int nsites = model->nsites_dis;		/* new nsites */
  P_PROB maxima = model->maxima;		/* the maxima positions */
  int n_samples = dataset->n_samples;		/* number of sequences */
  SAMPLE **samples = dataset->samples;		/* sequences */
  bool invcomp = model->invcomp;

  /* set all z to 0 */
  for (i=0; i<n_samples; i++) {			/* sequence */
    SAMPLE *s = samples[i];			/* sample i */
    int lseq = s->length;			/* length of sequence */
    int min_j = invcomp ? -lseq : 0;		// minimum Z_i
    int max_j = lseq;				// maximum Z_i
    Z_T *zi = s->z;				// zi[j], j in [-lseq...+lseq] 
    for (j=min_j; j<=max_j; j++) {		// Z_i = j
      Zi(j) = 0;
    }
  }

  /* set z to 1 for selected sites */
  for (i=0; i<nsites; i++) {
    SAMPLE *s = samples[maxima[i].x];		/* sample */
    if (maxima[i].negative){printf("NSITE %d is a negative !\n", i); exit(1);}
    int y = maxima[i].y;			/* position of site */
    Z_T *zi = s->z;				// zi[j], j in [-lseq...+lseq] 
    int j = maxima[i].ic ? -(y+1) : y+1;	// value of Z_i
    Zi(j) = 1.0;
  }

} /* set_z */

/***********************************************************************/
/*
	set_pY

	Initialize pY from z with given offset.
*/
/***********************************************************************/
static void set_pY(
  int w, 				/* motif width */
  bool invcomp, 			/* use reverse complement strand, too */
  bool pal,				/* force palindrome */
  DATASET *dataset,			/* the dataset */
  int off				/* offset to shift motif */
)
{
  int i, j;
  int n_samples = dataset->n_samples;		/* number of sequences */
  SAMPLE **samples = dataset->samples;		/* sequences */

  /* 
    put integerized, weighted log z into pY array
  */
  for (i=0; i<n_samples; i++) {			/* sequence */
    SAMPLE *s = samples[i];			/* sequence i */
    int lseq = s->length;			/* sequence length */
    int last_j = lseq-w;			/* last start */
    Z_T *zi = s->z;				// zi[j], j in [-lseq...+lseq] 
    double sw = s->sw;				/* weight of sequence */
    int *pY = s->pY[0];				/* p(Y_j | theta_1) both */
    char *pYic = s->pYic;			/* site on - strand */

    if (lseq < w) continue;			/* sequence too short */
    skip_group_if_required(dataset, s, i);	// ignore holdout groups

    /* initialize pY and pYic */
    for (j=0; j<=last_j; j++) {
      pY[j] = INT_LOG(0.0);			/* z == 0 */
      pYic[j] = '\0';				/* site on + strand */
    }
    for (j=0; j<lseq; j++) {			/* site start */
      int jj = j + off;				/* new site start */
      int k = jj+1;				// Z_i = k

      if (jj<0 || jj>last_j) continue;

      /* no z available? */
      if (j > last_j) {				/* no z available */
        pY[jj] = 0;				/* no site */
        pYic[jj] = '\0';			/* strand doesn't matter */
        continue;
      } 

      /* not using inverse strand, too? */
      if (!invcomp) {
        pY[jj] = INT_LOG(sw * Zi(k));
        pYic[jj] = '\0';			/* site on + strand */
        continue;
      } 

      /* using inverse complement strand, too */
      if (pal) {				// use sum of Zi(-k)+Zi(k)
        pY[jj] = INT_LOG(sw * MIN(1.0,(Zi(-k)+Zi(k))));	// FIXME??
      } else if (Zi(-k) > Zi(k)) {		// - strand
        pY[jj] = INT_LOG(sw * Zi(-k));
      } else {					// + strand
        pY[jj] = INT_LOG(sw * Zi(k));
      }
      pYic[jj] = (Zi(-k) > Zi(k)) ? '\1' : '\0';	/* choose strand */

    } /* site start */

  } /* sequence */

} /* set_pY */

/***********************************************************************/
/*
	maxima_compare

	Compare maxima based on "prob".
        Return >0 if second has the higher probability, then
	>0 if second is a control and first is a positive, then
	>0 if second has larger X, then
	>0 if second has larger Y.
*/
/***********************************************************************/
static int maxima_compare(
  const void *v1,
  const void *v2
)
{
  const struct p_prob *s1 = (const struct p_prob *) v1;
  const struct p_prob *s2 = (const struct p_prob *) v2;

  if (s1->prob != s2->prob) {
    return(s2->prob > s1->prob ? +1 : -1);
  } else if (s1->negative != s2->negative) {
    return(s2->negative ? +1 : -1);
  } else if (s1->x != s2->x) {
    return(s2->x - s1->x);
  } else {
    return(s2->y - s1->y);
  }
} // maxima_compare

/***********************************************************************/
/*
	sort_maxima

	In order to sort the maxima, we assign their "prob"
	field either Z (classic) or the LLR of the site (otherwise).
  	Break ties to sort negatives first.  That way the last
  	positive tied with them will yield the best p-value.
  	This should make p-values (slightly) conservative if nsites is limited.

        Sets the ranks of the positive maxima and removes the control 
	maxima if there are any.
*/
/***********************************************************************/
static void sort_maxima(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the primary dataset */
  DATASET *control,			/* the control dataset */
  int n_pos_maxima,			/* number of possible primary sites */
  int n_neg_maxima			/* number of possible control sites */
)
{
  int i, j;
  bool invcomp = model->invcomp;     	/* reverse complement strand */
  int n_maxima = n_pos_maxima + n_neg_maxima;	// total possible real/control sites

  //
  // Set the "prob" field in the maxima.  Use:
  //   objfun == Classic: z_i of site
  //   otherwise: erased LLR of the site
  //
  for (i=0; i<n_maxima; i++) {
    int x = model->maxima[i].x;
    bool neg = model->maxima[i].negative;
    SAMPLE *s = neg ? control->samples[x] : dataset->samples[x];
    int y = model->maxima[i].y;
    double score;
//FIXME: Classic might be improved if it used erased LLR, too.
// TLB: I left this alone to keep Classic identical with v4.10.2
// if sequence order shuffling is not done.
    if (dataset->objfun == Classic) {
      Z_T *zi = s->z;	 			// Zi[j], j in [-lseq...+lseq]
      int j = y + 1;				// Z_i = j
      score = invcomp ? MIN(1.0,Zi(-j)+Zi(j)) : Zi(j);
    } else {
      // Use the (erased) LLR of sites for sorting maxima.
      int w = model->w;
      LCB_T *lcb = s->logcumback;		/* log cumulative bkg. probability */
      THETA logtheta1 = model->logtheta;	/* motif log(theta) */
      double *not_o = s->not_o;			/* Pr(V_ij = 1) */
      double init = -Log_back(lcb, y, w) + LOG(not_o[y]);
      double llr = init;
      for (j=0; j<w; j++) llr += logtheta1(j, (int) s->res[y+j]);
      // negative strand
      if (invcomp) {
        double llr_rc = init;
        THETA logtheta1_rc = model->logtheta_rc;/* motif log(theta) rev. comp. */
        for (j=0; j<w; j++) llr_rc += logtheta1_rc(j, (int) s->res[y+j]);
        llr = MAX(llr, llr_rc);
      }
      score = llr;
    }
    RND(score, 11, model->maxima[i].prob);
  }

  //
  // Sort the maxima by "prob" which contains Z or LLR of maxima.
  // Break ties to sort negatives first.  That way the last
  // positive tied with them will have the best p-value.
  // This should make p-values (slightly) conservative if nsites is limited.
  //
  qsort((char *) model->maxima, n_maxima, sizeof(p_prob), maxima_compare);

  //
  // Set the ranks of the positive maxima and remove the control maxima if any
  //
  if (n_neg_maxima > 0) {
    for (i=j=0; i<n_maxima; i++) {
      if (! model->maxima[i].negative) {
        int rank = i+1;
	model->maxima[i].rank = rank;
	model->maxima[i].ranksum = (j==0) ? rank : model->maxima[j-1].ranksum + rank;
	model->maxima[j++] = model->maxima[i];
      }
    }
  }

  //
  // Get the weighted number of sites for each maxima.
  // Also, get the CD and CE metrics if appropriate.
  //
  SAMPLE **samples = dataset->samples;  	/* the sequences */
  double dist = 0;				// average distance of sites from center
  double wN = 0;				// weighted number of sites
  double in_region = 0;				// weighted sum of sites in central region
  for (i=0; i<n_pos_maxima; i++) {
    //
    // Compute the weighted number of sites for each rank.
    int x = model->maxima[i].x;
    int y = model->maxima[i].y;
    SAMPLE *s = samples[x];				/* sequence */
    double sw = s->sw;					/* sequence weight */
    // FIXME: We are assumming that priors are always symmetrical here.
    double esw = sw * INT_DELOG(s->log_not_o[y]);     	// Pr(site not overlapped) * Pr(site)
    wN += esw;                                		/* total sequence wgt */
    model->maxima[i].wN = wN;				// weighted rank
    //
    // Set the average distance from the sequence center of the sites up to each rank for objfun CD.
    // Set the count of sites within the central region up to each rank for objfun CE.
    if (dataset->objfun == CE || dataset->objfun == CD) {
      int w = model->w;
      //int d = fabs((y+w/2.0) - (s->length/2.0)) + 0.5;// distance between site-center and sequence-center
      int d = fabs(y + (w - s->length)/2.0) + 0.5; 	// distance between site-center and sequence-center
							// the 0.5 makes the not-has-zero case round up, but
							// leaves the other unchanged
      dist += esw * d;					// (weighted) sum of distance between site-centers and sequence-centers
      model->maxima[i].dist = wN ? dist/wN : dist;	// average distance to center up to this rank for CD
      if (d <= dataset->ce_max_dist) in_region += esw;	// (weighted) number of sites in central region for CE
      model->maxima[i].in_region = in_region;		// (weighted) number of sites in central region for CE
    } // CE or CD

  } // maxima

  if (debug) print_maxima(n_pos_maxima, model);

} // sort_maxima


/***********************************************************************/
/*
	get_column_re_list

  Return a list of (col_num, col_re) ordered by col_re,
  where col_re is the relative entropy of the motif column.
*/
/***********************************************************************/
  PAIR *get_column_re_list(
    MODEL *model, 
    DATASET *dataset
)
{
  int i;

  // Get the relative entropy (RE) of each motif column.
  calc_entropy(model, dataset);
  PAIR *pairs = NULL;
  Resize(pairs, model->w, PAIR);

  // Set up the index/RE pairs for sorting.
  for (i=0; i< model->w; i++) {
    pairs[i].ind = i;
    pairs[i].re = model->rentropy[i];
  }

  // Sort column indices by increasing column RE.
  qsort((void *)pairs, model->w, sizeof(PAIR), (void *)pairCompare);

  return(pairs);

} // get_column_re_list

/***********************************************************************/
/*
	estep_maxima_nsites

	1) Do an E_STEP on the primary and control sets.
	2) Get the maxima for one or both datasets.
	3) Sort all maxima together.
	4) Find the value of nsites corresponding to score threshold
	if one given (except for OOPS model, of course).
	5) Call get_best_nsites to find optimum nsites and its p-value
	unless the objfun is Classic or NC.

	Returns the best number of sites, and sets the log(p-value)
	and LLR of that p-value and the site probability threshold.
*/
/***********************************************************************/
int estep_maxima_nsites(
  MODEL *model,				/* the model */
  DATASET *primary,			/* the primary dataset */
  DATASET *control,			/* the control dataset */
  double (*E_STEP)(MODEL *, DATASET *),	// E_STEP function
  bool primary_groups[3],		// which primary groups to include
  bool control_groups[3],		// which control groups to include
  int min_nsites,			/* minimum nsites */
  int max_nsites,		 	/* maximum nsites */
  double thresh_in,			// restrict to maxima
					// with prob >= thresh_in
					// set to -BIG to ignore
  double *log_pv,			/* log p-value of score */
  double *llr, 				// LLR of best p-value
  double *thresh_out  			// optimal site prob threshold
)
{
  int i, j;
  int nsites = 0;
  bool invcomp = model->invcomp;     	/* reverse complement strand */
  bool pal = model->pal;     		/* force DNA palindromes */
  MOTYPE mtype = model->mtype;		/* type of model */
  OBJTYPE objfun = primary->objfun;	// objective function
  double col_scores[MAXSITE];		/* column scores */
  int n_pos_maxima = 0;
  int n_neg_maxima = 0;
  double wN = 0;

  // Do an E_STEP and get the maxima on the primary dataset.
  set_seq_groups_to_include(primary, primary_groups);
  set_seq_regions_to_include(primary, true, true, true);	// Include all sequence regions
  get_not_o(primary, model->w);
  model->ll = E_STEP(model, primary);		// set z using model
  set_pY(model->w, invcomp, pal, primary, 0);
  n_pos_maxima = get_max(mtype, primary, false, model->w, 0, model->maxima, invcomp, false);
  restore_seq_groups_to_include(primary);
  restore_seq_regions_to_include(primary);

  // Do an E_STEP and get the maxima on the control dataset.
  if (control) {
    set_seq_groups_to_include(control, control_groups);
    set_seq_regions_to_include(control, true, true, true);	// Include all sequence regions
    get_not_o(control, model->w);		// Not set yet so initialize.
    E_STEP(model, control);
    set_pY(model->w, invcomp, pal, control, 0);
    n_neg_maxima = get_max(mtype, control, true, model->w, n_pos_maxima, model->maxima, invcomp, false);
    n_neg_maxima -= n_pos_maxima;
    restore_seq_groups_to_include(control);
    restore_seq_regions_to_include(control);
  }

  //
  // * Sort all maxima (in both datasets together) by site LLR (in Z).
  // 	 Break ties to sort negatives first. That way the last
  // 	 positive tied with them will yield the best p-value.
  // 	 This should make p-values (slightly) conservative if nsites is limited.
  // * Rank the positive maxima.
  // * Compute average distance to center at each rank for objfun==CD,
  // *   number of sites in central region for CE.
  // * Remove the control maxima.
  //
  sort_maxima(model, primary, control, n_pos_maxima, n_neg_maxima);

  //
  // Find the value of nsites corresponding to the score threshold.
  // (Don't do this for OOPS model).
  //
  if (mtype != Oops && thresh_in != -BIG) {
    // Binary search of maxima for maxima[i-1].prob < thresh_in <= maxima[i].prob
    int lo = MAX(0, primary->min_nsites-1);
    int hi = MIN(n_pos_maxima-1, primary->max_nsites-1);
    while (hi-lo > 1) {			// binary search
      int mid = (lo+hi)/2;		// midpoint
      if (model->maxima[mid].prob >= thresh_in) { lo = mid; } else { hi = mid; }
    }
    nsites = model->maxima[hi].prob >= thresh_in ? hi+1 : lo+1;
    max_nsites = MIN(max_nsites, nsites);
    max_nsites = MAX(min_nsites, max_nsites);
    min_nsites = max_nsites;		// force this number of sites in get_best_nsites()
  } else if (mtype == Oops) {
    nsites = min_nsites;
  }

  //
  // Get the best number of sites and get the p-value for non-Classic models.
  //
  if (! (objfun == Classic || objfun == NC) ) {
    double log_ev;
    nsites = get_best_nsites(model, primary, min_nsites, max_nsites,
      model->w, n_pos_maxima, n_neg_maxima, model->maxima, col_scores, &wN, log_pv, &log_ev, llr, thresh_out);
  }

  return(nsites);
} // estep_maxima_nsites

/***********************************************************************/
/*
	classic_get_width_nsites_pvalue

	Old documentation:
		1) get best nsites using E-value
		2) calculate p-value of each column of motif
		3) shorten using p-value
		4) get best nsites using E-value

	Adjust the width of the motif by
		1) using the multiple alignment trim procedure	
		2) optimizing E-value over all subsets of columns

	Updates the best motif information into the model.

	Returns the best starting offset.
*/
/***********************************************************************/
static int classic_get_width_nsites_pvalue(
  MODEL *model,				/* the model */
  DATASET *dataset, 			/* the dataset */
  double (*E_STEP)(MODEL *, DATASET *)	// E_STEP function
)
{
  int i, w, ini_w, best_w, min_w, max_w, ma_w, ma_off;
  int n_maxima; 				/* number of possible sites */
  MOTYPE mtype = model->mtype;			/* type of model */
  int nseqs = dataset->n_group[0];		// number of sequences for EM
  int csites = 					// maximum number of sites possible now
    (mtype==Tcm) ? dataset->classic_max_nsites : MIN(nseqs, dataset->classic_max_nsites);
  int min_nsites = MIN(csites, dataset->min_nsites);	/* minimum nsites */
  int max_nsites = MIN(csites, dataset->max_nsites);	/* maximum nsites */
  bool invcomp = model->invcomp;     	/* reverse complement strand */
  bool pal = model->pal;     		/* force DNA palindromes */
  double col_scores[MAXSITE];			/* column scores */
  double score;					// score of motif (-log_pop or LLR)
  double wN;					// weighted number of sites
  double log_pv, best_log_pv;			/* log p-value of motif */
  double log_ev;				/* log E-value of score */
  double llr;					// LLR of best E-value
  double thresh;				// site prob threshold
  OBJTYPE objfun = dataset->objfun;		// objective function

  //
  // initialize pY from Z
  //
  set_pY(model->w, invcomp, pal, dataset, 0);

  // Don't include the PSP probabilities in log_not_o (used in get_max())
  // because z already takes them into account.
  get_not_o(dataset, model->w);

  //
  // Get the maxima.
  //
  n_maxima = get_max(mtype, dataset, false, model->w, 0, model->maxima, invcomp, false);

  if (n_maxima < min_nsites) return 0; 	// Exit if not enough sites.

  //
  // Sort the maxima.
  //
  sort_maxima(model, dataset, NULL, n_maxima, 0);

  //
  // Set the minimum and maximum widths.
  //
  ini_w = best_w = model->w;			/* initial motif width */
  min_w = model->min_w;				/* minimum width */
  int best_off = 0;

  //
  // Get the best nsites using the full motif width.
  // 
  int best_nsites = (min_nsites==max_nsites) ? min_nsites :
    get_best_nsites(model, dataset, min_nsites, max_nsites,
      model->w, n_maxima, 0, model->maxima, col_scores, &wN, &log_pv, &log_ev, &llr, &thresh);

  //
  // Trim the alignment to include the minimum possible number of gaps
  //
  max_w = ma_w = ini_w;
  ma_off = 0;
  w = 0;
  if (dataset->ma_adj) {
    int flank = w/2;				/* amt. of flank for mult. a. */
    ma_w = ma_adjust(model->maxima, best_nsites, max_w, flank, min_w, dataset, &ma_off);
    max_w = ma_w;				/* new maximum width */
    /*
      update the maxima positions by shifting them
    */
    for (i=0; i<best_nsites; i++) {
      bool ic = model->maxima[i].ic;			/* on - strand */
      model->maxima[i].y += (ic ? ini_w-ma_w-ma_off : ma_off);
    }
  } /* ma_adj */

  //
  // Shorten the motif based on p-value of alignments.
  //
  /* 
    get the p-values of the columns given the best number of sites
    and the gap-trimmed width
  */
  (void) get_best_nsites(model, dataset, best_nsites, best_nsites,
    max_w, n_maxima, 0, model->maxima, col_scores, &wN, &log_pv, &log_ev, &llr, &thresh);

  /* 
    find subsequence of columns with best p-value
  */
  best_log_pv = log_ev;
  if (objfun == Classic) {
    best_log_pv = log_ev = BIG;
    for (w=min_w; w<=max_w; w++) {		/* width */
      int l_off;					/* left offset */

      for (l_off=0; l_off<=max_w-w; l_off++) {	/* left edge */

	/* get the product of column p-values */
	for (i=score=0; i<w; i++) score += col_scores[i+l_off];

	/* get the p-value of the pop */
	log_pv = get_log_sig(-score, mtype, w, best_nsites, 0, invcomp, pal, dataset);

	if (TRACE) 
	  printf(
	  "ini_w %d ma_w %d w %d ma_off %d off %d log_pv %f init cons %*.*s\n", 
	  ini_w, ma_w, w, ma_off, l_off, log_pv, w, w, model->cons0+l_off+ma_off);

	/* save if best so far: better log_pv */
	if (RNDEPS < best_log_pv - log_pv) {
	  if (TRACE) printf("better: w %d best log_pv %g\n", w, log_pv);
	  best_w = w;
	  best_off = l_off;
	  best_log_pv = log_pv;
	}
      } /* l_off */
    } /* w  */

    /*
      update the maxima positions by shifting them
    */
    for (i=0; i<best_nsites; i++) {
      bool ic = model->maxima[i].ic;			/* on - strand */
      model->maxima[i].y += (ic ? ma_w-best_w-best_off : best_off);
    }
  }
  
  /* 
    get the best number of sites for the shortened motif and the final E-value 
  */
  best_nsites =
    get_best_nsites(model, dataset, min_nsites, best_nsites,
      best_w, n_maxima, 0, model->maxima, col_scores, &wN, &log_pv, &log_ev, &llr, &thresh);

  /* 
    set the best motif info in the model
  */
  model->w = best_w;				/* best width */
  model->nsites_dis = MAX(MINSITES, best_nsites);	/* after discretization; not correct if samples held out */
  model->logpv = best_log_pv;			/* p-value */
  model->logev = log_ev;			/* E-value */
  model->llr = llr;				// LLR of best E-value
  model->site_threshold = thresh;		// score threshold

  return(best_off);
} // classic_get_width_nsites_pvalue

/***********************************************************************/
/*
	smhg_trim_motif

	Determine the set of columns to mask, masking
	some number of lowest relative entropy columns
	to the uniform distribution.
	Shift the motif (and shorten it) so that 
	there are no masked columns on the flanks.
        Use estep_maxima_best_nsites() with the *masked* model
	to determine nsites, LLR and pvalue.
	Replace model->theta with the same span, but using
	the unmasked versions of all columns.

	Sets model->w, ->nsites_dis, ->LLR, ->logpv and ->logev.

	Returns the best offset of the original motif, and
	sets site threshold and nsites.
*/
/***********************************************************************/
static int smhg_trim_motif(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the dataset */
  DATASET *neg_dataset,			/* the control dataset */
  double (*E_STEP)(MODEL *, DATASET *),	// E_STEP function
  double *threshold,			// best site threshold
  int *nsites				// best value of nsites
)
{
  int i, j;
  int min_w = model->min_w;			/* minimum width */
  bool invcomp = model->invcomp;     		/* reverse complement strand */
  bool pal = model->pal;     			/* force DNA palindromes */
  ALPH_T *alph = dataset->alph;			// alphabet
  GROUP_T primary_groups, control_groups;
  int min_nsites = dataset->min_nsites;		/* minimum nsites */
  int max_nsites = dataset->max_nsites;		/* maximum nsites */
  double log_pv;				/* log p-value of score */
  double llr;					// LLR of best E-value

  // Save a copy of the motif.
  THETA theta = model->theta;			//theta of motif
  if (saved_theta == NULL) {
    create_2array(saved_theta, double, model->max_w+1, alph_size_wild(alph));
  }
  copy_theta(theta, saved_theta, model->w, alph_size_wild(alph));

  //
  // Get a list of (col_num, col_re) ordered by col_re,
  // where col_re is the relative entropy of the motif column.
  //
  PAIR *pairs = get_column_re_list(model, dataset);

  //
  // Decide how many hi-RE columns to keep and hence the motif width.
  //
  int n_keep;
  // Tanaka 2014 code
  n_keep = model->w > 6 ? 6 + floor((model->w-6)/3) : model->w;
  //
  //n_keep = ceil(2 * sqrt(model->w));
  //n_keep = min_w + ceil(sqrt(model->w));
  //n_keep = ceil(sqrt(model->w));
  //n_keep = min_w + 2*ceil(sqrt(model->w));
  //n_keep = min_w + ceil(sqrt(model->w - min_w));
  //
  // Progression code:
  // Use sqrt(2) progression just like start widths in starts.c.
  // n_keep = min_w + floor((model->w - min_w)/3);
  // int prev_w = (int) (model->w/sqrt(2) + 0.5);	// prev w in progression
  // Make sure all smaller widths are reachable.
  // if (n_keep > prev_w + 1) n_keep = prev_w + 1;
  //
  // Make sure not keeping too many.
  if (n_keep > model->w) n_keep = model->w;
  // Make sure not keeping too few.
  if (n_keep < model->min_w) n_keep = model->min_w;	
  // How many columns to mask:
  int n_mask = model->w - n_keep;

  // Mask non-selected, low-RE columns.
  double f = 1.0/alph_size_core(alph);		// uniform frequency
  for (i=0; i<n_mask; i++) {
    int index = pairs[i].ind;
    for (j=0; j<alph_size_core(alph); j++) theta(index, j) = f;
    theta(index, j) = 0;			// not sure if necessary
  }

  // Get the start and end indices of the non-masked columns.
  int small_index = model->w - 1;
  int big_index = 0;
  for (i=model->w-1; i>=n_mask; i--) {
    int index = pairs[i].ind;
    if (index < small_index) small_index = index;
    if (index > big_index) big_index = index;
  }
  myfree(pairs);

  // Shift model and set new width (contains only non-masked columns).
  if (n_mask > 0) {
    int new_w = big_index - small_index + 1;
    // Increase width of motif to at least min_w, centering.
    if (new_w < min_w) {
      int left_margin = small_index;
      int right_margin = model->w - big_index - 1;
      while (new_w < min_w) {
	// Increase and center.
	if (left_margin >= right_margin) {
	  left_margin = --small_index;
	} else {
	  right_margin = model->w - (++big_index) - 1;
	}
	new_w++;
      }
    }
    if (small_index > 0) copy_theta(model->theta+small_index, model->theta, new_w, alph_size_wild(alph));
    model->w = new_w;
  }

  //
  // Use the mHG test to determine the best number of sites using the *masked* motif
  // and to get a (biased) p-value for the motif. 
  //
  double thresh_in = -BIG;
  double thresh_out;
  int best_nsites = estep_maxima_nsites(model, dataset, neg_dataset, E_STEP, 
    dataset->primary_groups.trim, dataset->control_groups.trim,
    min_nsites, max_nsites, thresh_in, &log_pv, &llr, &thresh_out);
  double best_thresh = thresh_out;

  //
  // Update the model using the unmasked span from the original theta matrix.
  //
  copy_theta(saved_theta+small_index, model->theta, model->w, alph_size_wild(alph));
  model->nsites_dis = (min_nsites==max_nsites) ? min_nsites : best_nsites;
  model->logpv = log_pv;			/* p-value */
  model->logev = log_pv;			/* E-value */
  model->llr = llr;				// LLR of best E-value

  *threshold = best_thresh;
  *nsites = best_nsites;
  return(small_index);
} // smhg_trim_motif


/***********************************************************************/
/*
	greedy_trim_motif

	Sort the columns in order of decreasing RE.
	Get score for each motif consisiting of the first i sorted columns
	such that the spanned width is within the allowed range.

	Sets model->theta to the selected range of columns
	containing their original values.
	Sets model->w to the new width.
        Sets model->logev to best log_pv.
	Sets model->llr to best llr.

	Returns the best offset and sets site threshold and final_nsites.
*/
/***********************************************************************/
static int greedy_trim_motif(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the dataset */
  DATASET *neg_dataset,			/* the control dataset */
  double (*E_STEP)(MODEL *, DATASET *),	// E_STEP function
  double *threshold,			// best site threshold (OUT)
  int *final_nsites			// best value of nsites
) {
  int i, w;
  int ini_w = model->w;				// initial motif width
  DATASET *control = NULL;
  OBJTYPE objfun = dataset->objfun;		// objective function
  bool pal = model->pal;     		/* force DNA palindromes */
  MOTYPE mtype = model->mtype;			/* type of model */
  int min_nsites = dataset->min_nsites;		/* minimum nsites */
  int max_nsites = dataset->max_nsites;		/* maximum nsites */
  ALPH_T *alph = dataset->alph;			// alphabet
  double log_pv;				/* log p-value of score */
  double llr;					// LLR of best p-value
  int nsites;
  double thresh_in = -BIG;			// site threshold (for input)
  double thresh_out;				// site threshold (for output)
  double scale=0;				// for adjusing min and max nsites

  // Save a copy of the motif.
  THETA theta = model->theta;			// theta of motif
  if (saved_theta == NULL) {
    create_2array(saved_theta, double, model->max_w+1, alph_size_wild(alph));
  }
  copy_theta(theta, saved_theta, model->w, alph_size_wild(alph));

  //
  // Get a list of (col_num, col_re) ordered by increasing col_re,
  // where col_re is the relative entropy of the motif column.
  //
  PAIR *pairs = get_column_re_list(model, dataset);

  //
  // Find the best motif using the objective function applied to each span of high-RE columns.
  //
  int best_start = 0;
  int best_end = 0;
  int best_nsites = 0;
  double best_log_pv = BIG;
  double best_llr = -BIG;
  double best_threshold = -BIG;
  int best_w = model->w;
  int maxsites=0, minsites=0;

  //
  // Set the parameters for scoring this motif.
  //
  if (objfun == NZ) {
    // Set parameters to score motif using mHG on primary group 0,2 vs primary group 1.
    // Note that group 0 are positive sequences, and groups 1 and 2 are negative sequences.
    control = dataset;
    scale = (double) (dataset->n_group[0] + dataset->n_group[2])/dataset->n_group[0];
    thresh_in = -BIG;
  } else if (objfun == DE) {
    // Set parameters to score motif using mHG on group 1 (holdout 1) vs. control group 0.
    control = neg_dataset;
    scale = (double) dataset->n_group[1]/dataset->n_samples;
    thresh_in = -BIG;
  } else if (objfun == CE || objfun == CD) {
    // Set parameters to score motif using CE or CD on group 1.
    control = NULL;
    scale = (double) dataset->n_group[1]/dataset->n_samples;
    thresh_in = -BIG;
  }

  // Adjust nsites range for sequences added to (NZ) or removed from (DE) primary sequences.
  minsites = MAX(2,scale*min_nsites);
  maxsites = MAX(2,scale*max_nsites);

  //
  // Add motif columns in order of decreasing relative entropy.
  // Stop when they span the maximum width.
  //
  int small_index, big_index;
  small_index = big_index = pairs[ini_w-1].ind;		// highest RE column
  for (i=ini_w-2; i>=0; i--) {
    // Get index of next highest RE column.
    int index = pairs[i].ind;
    if (index < small_index) {
      small_index = index;
    } else if (index > big_index) {
      big_index = index;
    } else {
      continue;				// no change in motif happened
    }
    int new_w = big_index - small_index + 1;

    // Motif width must be >= mininimum width
    if (new_w < model->min_w) continue;

    // Replace the motif with the new span of high-RE columns.
    if (small_index > 0) copy_theta(saved_theta+small_index, model->theta, new_w, alph_size_wild(alph));
    model->w = new_w;

    // Score motif using the parameters set above and find nsites and threshold.
    nsites = estep_maxima_nsites(model, dataset, control, E_STEP, 
      dataset->primary_groups.trim, dataset->control_groups.trim,
      minsites, maxsites, thresh_in, &log_pv, &llr, &thresh_out);

    //
    // Save motif offset, nsites and p-value information if best found so far.
    // Break ties using 1) width 2) LLR.
    //
    if (
      (log_pv < best_log_pv)
      || (log_pv == best_log_pv && new_w < best_w) 
      || (log_pv == best_log_pv && new_w == best_w && llr > best_llr)
    ) {
      best_start = small_index;
      best_end = big_index;
      best_w = best_end - best_start + 1; 
      best_nsites = nsites;
      best_log_pv = log_pv;
      best_llr = llr;
      best_threshold = thresh_out;
    }

    // Done if span is same as initial motif width;
    if (new_w == ini_w) break;
  }
  myfree(pairs);

  //
  // Extract the best columns from the motif and re-set its width.
  //
  model->w = best_w;
  if (best_start > 0) copy_theta(saved_theta+best_start, model->theta, best_w, alph_size_wild(alph));

  //
  // Store the SCORE and LLR of the best motif found above.
  //
  model->logev = best_log_pv;			// SCORE for comparing models
  model->llr = best_llr;			// LLR of best SCORE

  *threshold = best_threshold;
  *final_nsites = best_nsites;

  return(best_start);

} // greedy_trim_motif

/***********************************************************************/
/*
	get_width_nsites_pvalue

	Find the optimal number of sites and widths for the motif by
	optimizing:
		NZ: mHG (or mRS or mBN) on groups 0,2 vs. 1
		DE: mHG (or mRS or mBN) on group 1 vs. control
		LL: LLR on primary sequences (groups 0,1,2)
		CE/CD: CE/CD on group 1
	Note: In OOPS mode, the Wilcoxon rank-sum test replaces
	mHG.

	Widths are searched by using the motif columns
	spanned by the first 2, 3, ..., w highest-RE columns.
	The optimal site threshold is the one that gives the
	best value of the objective function (mHG, or LLR).
	Sets the best site (mHG, mRS, mBN or LLR) in model->logev, 
        sets the best LLR in model->llr, and the best nsites
	in model->nsites_dis. 

	These values are to be used for optimization over starting points.

	Find the final p-value for the motif using the optimal
	site threshold determined above and:
		NZ: Fisher exact (or rank-sum) test on group 0 vs. 2 
		DE: Fisher exact (or rank-sum) test on group 2 vs. control
		LL: Fisher exact (or rank-sum) test on primary vs. control
		CE/CD: CE/CD test on holdout group 2

	Sets model->w to the best width.
	Sets model->nsites_dis to the best value.
        Sets model->logpv to the best value.
        Sets model->logev to best log_pv.
	Sets model->llr to best llr.
        Sets model->site_threshold to best value.

	Returns the best starting offset.
*/
/***********************************************************************/
static int get_width_nsites_pvalue(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the dataset (including noise) */
  DATASET *neg_dataset,			// the control dataset
  double (*E_STEP)(MODEL *, DATASET *)	// E_STEP function
)
{
  int i, w;
  int ini_w = model->w;				// initial motif width
  DATASET *control = NULL;
  OBJTYPE objfun = dataset->objfun;		// objective function
  bool pal = model->pal;     		/* force DNA palindromes */
  MOTYPE mtype = model->mtype;			/* type of model */
  int min_nsites = dataset->min_nsites;		/* minimum nsites */
  int max_nsites = dataset->max_nsites;		/* maximum nsites */
  ALPH_T *alph = dataset->alph;			// alphabet
  double log_pv;				/* log p-value of score */
  double llr;					// LLR of best p-value
  int nsites;
  double thresh_in = -BIG;			// site threshold (for input)
  double thresh_out = -BIG;			// site threshold (for output)
  double scale=0;				// for adjusing min and max nsites
  double best_threshold = 0;			// best site threshold
  int best_nsites = 0;				// best value of nsites
  int best_start = 0;				// best offset into original motif

  //
  // Determine the best width and nsites for the motif.
  // Returns trimmed motif.
  //
  if (objfun == SE) {
    best_start = smhg_trim_motif(model, dataset, neg_dataset, E_STEP, &best_threshold, &best_nsites);
    // Done if not using DE to get an unbiased p-value for SE.
    if (dataset->hsfrac==0) { return(best_start); }
  } else {
    best_start = greedy_trim_motif(model, dataset, neg_dataset, E_STEP, &best_threshold, &best_nsites);
  }

  //
  // Set up parameters for estimating the final p-value of the motif.
  //
  if (objfun == NZ) {
    // Set up parameters to use the Fisher exact (or rank-sum) test with best_threshold on group 0 vs. 2.
    control = dataset;
    thresh_in = best_threshold;
    scale = 1;
  } else if (objfun == DE) {	// discriminative enrichment
    // Set up parameters to use the Fisher exact (or rank-sum) test on group 2 (holdout2) vs. control groups 1&2.
    control = neg_dataset;
    scale = (double) dataset->n_group[2]/dataset->n_samples;
    thresh_in = best_threshold;
  } else if (objfun == SE) {	// SmHG, SmBN, SmRS
    // Set up parameters to use the Fisher exact (or rank-sum) test on groups 1&2 (holdout1+holdout2) vs. control groups 1&2.
    control = neg_dataset;
    scale = (double) (dataset->n_group[1]+dataset->n_group[2])/dataset->n_samples;
    thresh_in = best_threshold;
    // Fake out estep_maxima_nsites
    dataset->objfun = DE;
  } else if (objfun == CE || objfun == CD) {	// central enrichment, central distance
    // Use CE or CD on group 2
    control = NULL;
    thresh_in = best_threshold;
    scale = (double) (dataset->n_group[2])/dataset->n_samples;
  }

  //
  // Get the final (unbiased) p-value using type of test determined by the values set above.
  //
  int minsites = MAX(2, scale*min_nsites);
  int maxsites = MAX(2, scale*max_nsites);
  nsites = estep_maxima_nsites(model, dataset, control, E_STEP, 
    dataset->primary_groups.pvalue, dataset->control_groups.pvalue, 
    minsites, maxsites, thresh_in, &log_pv, &llr, &thresh_out);
  model->logpv = log_pv;			/* FINAL p-value */

  if (debug) printf("p-value: ini_w %d w %d offset %d best_nsites %d new_nsites %d llr %g final logpv %g cons0 %s selection pv %g final pv %g\n", ini_w, model->w, best_start, best_nsites, nsites, llr, log_pv, model->cons0, exp(model->logev), exp(log_pv));

  //
  // Set up parameters to get the optimum value of nsites and to set the maxima for use by em().
  //
  if (objfun == NZ) {
    // Set up to use mHG on group 0 vs. groups 1,2. 
    control = dataset;
    thresh_in = -BIG;
  } else if (objfun == DE || objfun == SE) {
    // Set up to use mHG on groups 0,1,2 vs. control.
    control = neg_dataset;
    thresh_in = -BIG;
  } else if (objfun == CE || objfun == CD) {
    // Set up to use CE or CD on groups 0,1,2 (all primary sequences)
    control = NULL;
    thresh_in = -BIG;
  }

  //
  // Get the final nsites and set the maxima using the parameters set above.
  //
  minsites = min_nsites;
  maxsites = max_nsites;
  model->nsites_dis = estep_maxima_nsites(model, dataset, control, E_STEP, 
    dataset->primary_groups.nsites, dataset->control_groups.nsites,
    minsites, maxsites, thresh_in, &log_pv, &llr, &thresh_out);
  if (debug) printf("TUNED NSITES using ALL: ini_w %3d w %3d offset %4d selection p-value %.4g unbiased pv %.4g pv on full %.2g new llr %g selection nsites %3d full nsites %3d\n\n", ini_w, model->w, best_start, exp(model->logev), exp(model->logpv), exp(log_pv), llr, best_nsites, model->nsites_dis);

  //
  // Restore objfun (if needed).
  //
  dataset->objfun = objfun;

  // Return offset of final motif from start of original motif.
  return(best_start);
} // get_width_nsites_pvalue

/***********************************************************************/
/*
	discretize	

	Search over width and offset of motif to minimize the E-value
	of the chosen statistic.

	Sets z to 1.0 for sites, 0 for non-sites in dataset.
	Sets w, nsites_dis, maxima and sig in model.

	Returns the optimum number of sites.
*/
/***********************************************************************/
void discretize(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the dataset */
  DATASET *neg_dataset,			/* the control dataset */
  double (*E_STEP)(MODEL *, DATASET *)	// E_STEP function
)
{
  int i, j;
  MOTYPE mtype = model->mtype;			/* type of model */
  OBJTYPE objfun = dataset->objfun;		// objective function

  /* 
    set initial and minimum allowed widths
  */
  int ini_w = model->w;				/* initial motif width */
  int min_w = model->min_w;			/* minimum width */

  /* 
    create space for the maxima from both datasets
  */
  int n_maxima = (mtype==Tcm) ? ps(dataset, min_w) : dataset->n_samples;
  if (objfun == DE || objfun == SE) {
    n_maxima += (mtype==Tcm) ? ps(neg_dataset, min_w) : neg_dataset->n_samples;
  }
  Resize(model->maxima, n_maxima, p_prob);

  /* 
    Get best width and nsites for motif.
  */
  int best_off = 0;
  if (objfun==Classic || objfun==NC) {
    best_off = classic_get_width_nsites_pvalue(model, dataset, E_STEP);
  } else if (objfun==NZ) {
    DATASET *primary = dataset;
    DATASET *control = NULL;
    best_off = get_width_nsites_pvalue(model, primary, control, E_STEP);
  } else if (objfun==DE || objfun==SE) {
    DATASET *primary = dataset;
    DATASET *control = neg_dataset;
    best_off = get_width_nsites_pvalue(model, primary, control, E_STEP);
  } else if (objfun==CE || objfun==CD) {
    DATASET *primary = dataset;
    DATASET *control = NULL;
    best_off = get_width_nsites_pvalue(model, primary, control, E_STEP);
  } else {
    fprintf(stderr, "\nUnknown type of objective function.\n");
    exit(1);
  }
  
  if (!NO_STATUS) {
#if defined(PARALLEL) && defined(DEBUG_PARALLEL)
      printf("\nnode %d FINAL: ini_w %d w %d nsites %d selection %.17g pv %g llr %g\n", mpMyID(), ini_w, model->w, model->nsites_dis, exp(model->logev), exp(model->logpv), model->llr);
#else
      if (debug) printf("\nFINAL: ini_w %d w %d nsites %d selection %g p-value %g llr %g\n", ini_w, model->w, model->nsites_dis, exp(model->logev), exp(model->logpv), model->llr);
#endif
  }

  //
  // If some sequences were held out due to searchsize, put them back in
  // and get a better value of nsites using the score threshold found by em().
  //
  if (
    (objfun == Classic || objfun == NC)
      && (dataset->n_group[1] > 0 || dataset->n_group[2] > 0)
  ) {
    double dummy;
    double (*E_STEP)(MODEL *, DATASET *);        // expectation step
    E_STEP = model->mtype == Tcm ? tcm_e_step : e_step;
    int best_nsites = estep_maxima_nsites(model, dataset, NULL, E_STEP,
      dataset->primary_groups.nsites, dataset->control_groups.nsites,
      dataset->min_nsites, dataset->max_nsites, model->site_threshold, &dummy, &(model->llr), &dummy);
    model->nsites_dis = MAX(MINSITES, best_nsites);
  }

  //
  // Discretize the sites in Z from the maxima stored in the model.
  //
  set_z(model, dataset);

  if (TRACE) {
    int w = model->w;
    printf( 
      "ini_w %d w %d off %d nsites %d pv %9.2e EV %9.2e cons %*.*s %s\n",
      ini_w, w, best_off, model->nsites_dis, exp(model->logpv), 
      exp(model->logev), w, w, (model->cons0)+best_off, model->cons0);
  }

} /* discretize */

/**********************************************************************/
/*
	ma_adj

	Shorten a motif to the longest g-alignment of width at least
	min_w.  A g-alignment is an alignment with no more than g
	gapped sequences per column.  Values of g in [0..] are tried
	until an aligment of width min_w or greater is found.

	Returns best width and offset.
*/
/**********************************************************************/
static int ma_adjust(
  P_PROB sites,                                 /* the sites */
  int nsites,                                   /* the number of sites */ 
  int w,                                        /* width of sites */ 
  int flank,                                    /* add flank cols on left+rgt */
  int min_w,		 			/* minimum width allowed */
  DATASET *dataset,  				/* the dataset */
  int *off					/* best offset */
)
{
  char **ma;					/* multiple aligment */
  int left = MIN(flank, sites[0].y);		/* left edge after algnmnt. */
  int right = left+w-1;				/* right edge after algnmnt. */

  /* get the multiple alignment */
  ma = dp_multi_align(sites, nsites, w, flank, dataset);

  /* get longest g-alignment of width at least min_w */
  (void) g_align(ma, nsites, strlen(ma[0]), left, right, min_w, off, &w);

  /* free space */
  free_2array(ma, nsites);

  return w;					/* return the width */
} /* ma_adjust */

