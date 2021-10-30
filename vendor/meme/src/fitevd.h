#ifndef FITEVD_H

#include <stdio.h>
//#include <malloc.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include "mhmm-state.h"
#include "utils.h"

/******************************************************************************/
/*
  	constants
*/
/******************************************************************************/
#define BLKSIZ 1000
#define EPS1 1e-10		/* log-likelihood convergence criteria */
#define EPS2 1e-8		/* convergence criterion for first derivative */
#define MAXH 1000 		/* largest value of H to allow to prevent
				   problems with EPS2 */
#define MIN_SCORES 200		// Minimum number of scores for estimation.
//#define MAX_WIMP 0.9		// Maximimum allowed probability of a false "hit".
#define MAX_WIMP 1.0		// Maximimum allowed probability of a false "hit".
/* compute the p-value of value of n samples of an extreme value whose
   p-value is known
*/
#define EV(p, n) \
        (((p)*(n) < .001) ? (p)*(n) : 1.0 - pow((1.0 - (p)), (double)(n)))

//FIXME: TLB global for printing score sets
extern char *OUTPUT_DIRNAME;

/******************************************************************************/
/* 
  	data types
*/
/******************************************************************************/
// Type of score distribution to fit.
typedef enum {
  D_EVD,        // extreme value distribution
  D_EXP,        // 2-component: Exponential (mu_1), Gaussian(mu_2, sd_2); binned GC-content
  D_GCEXP       // 2-component: Exp(mu_1(g)=a+bg**d), Gaussian(mu_2, sd_2); (g=GC-content)
} DISTR_T;

/* score and associated info */
typedef struct score {
  double s;                             /* similarity score of two sequences */
  int t;				/* length of sequence */
  double gc;				// GC content around match.
  int nhits;				/* number of hits in match (repeated match) */
  int span;				/* length of match (repeated match) */
  bool ok_span;			// Span not too long?
  double pv;				/* p-value of score */
  char *id;				/* optional string identifier */
} SCORE;

/* set of score and length from a similarity search */
typedef struct score_set {
  int n;                                /* number of scores in set */
  int max_scores_saved;                 /* Counter for reservoir sampling */
  int num_scores_seen;                  /* Counter for reservoir sampling */
  SCORE *scores;                        /* array containing X = {(s_i,t_i)} */
  int max_t;                            /* maximum target length (t_i) */
  int q;                                /* length of query */
  bool negatives_only;		// All negatives?
  int max_gap;				// Maximum gap allowed.
  double total_length;   		// Total size of database.
  double min_e;				// Mininum E-value of scores.
  double avgw;				// Average motif width.
  double phit;				// Probability of a false hit.
  double wimp;				// Probability of false hit within max_gap.
  double egap; 	                	// E[gap length]
  double ehits; 	                // E[hits/match]
  double espan;				// E[length of match]
  double e_hit_score;			// E[score of single hit]
  double avg_min_pvalue;		// Geometric mean of min p-values for all motifs.
} SCORE_SET;

/* parameters of a distribution (exponential, gc-dependent-exponential or EVD) */
typedef struct evd {
  // Extreme value parameters.
  double lambda;                        /* parameters */
  double K;
  double H;
  int n;				/* number of scores used */
  // Other stuff for extreme values.
  int g;				/* length group number */
  double min_t;				/* min_t < t (or gc) <= max_t */
  double max_t;
  double mid_t;				/* midpoint of length range */
  // Exponential parameters.
  double mu1;				// Exponential mean.
  double mu2;				// Gaussian mean.
  double sigma2;			// Gaussian sd.	
  double c;				// Mixing parameter.
  // GC linear parameters.
  double a;                             // mu1(g) = a+b*g
  double b;                             // mu1(g) = a+b*g
  // Log likelihood (applies to both types of distributions).
  double L;                             /* log-likelihood */
} EVD;

/* set of distributions for different length ranges */
typedef struct evd_set {
  int n;				// Number of length ranges.
  double total_length;			// Database total length.
  int nscores;				// Number of scores.
  int outliers;				// # E-values < 1.
  int non_outliers;			// Total non-outliers among scores.
  int N;				// E = N p.
  DISTR_T dtype;			// Type of Score Distribution
  bool negatives_only;		// Distribution from synthetiic sequences.
  double min_e;				// Minimum Evalue.
  double sum_log_e;			// Sum log(Evalue < 1).
  char msg[300];			// Success/failure message.
  EVD *evds;				// Array of distributions.
} EVD_SET;

/* derivatives of extreme value distribution */
typedef struct deriv {
  double d1;				/* first derivative */
  double d2;				/* second derivative */
  double L;                             /* expected augmented log likelihood */
  double K;                             /* optimum value of K */
} DERIV;

/******************************************************************************/
/*
  	macros
*/
/******************************************************************************/

/******************************************************************************/
/*
	Exp_pvalue
	Get the p-value of a score using the exponential distribution.
*/
/******************************************************************************/
#define Exp_pvalue(evd, s) exp(-s/(evd).mu1)

/******************************************************************************/
/* 	
	Evd_set_evalue
	Get the E-value of a score.

		n * p-value(s, l, q)
*/
/******************************************************************************/
#define Evd_set_evalue(n, s, l, q, evd_set)		 		\
    ((n) * evd_set_pvalue((s), (l), (q), (evd_set)))

/******************************************************************************/
/* 	
	Push_score
	Push a score onto the end of a score set 
*/
/******************************************************************************/
#define Push_score(sc, set) {						\
  if ((set).n % BLKSIZ == 0)						\
    Resize((set).score, (set).n + BLKSIZ, SCORE);			\
  (set).score[(set).n++] = (sc);					\
}


/******************************************************************************/
/* 
	Score_bin
		s	maximum score

  	Get the number of bins to place scores into.
*/
/******************************************************************************/
#define Score_bin(s) ((int) (10*log(s)))

/******************************************************************************/
/*
	Evd_pvalue

	Get the pvalue of a score.  Returns pv, N and el.
*/
/******************************************************************************/
#define Evd_pvalue( 							\
  pv,				/* pvalue */				\
  N,				/* effective search space */		\
  el,				/* expected alignment length */		\
  evd, 				/* extreme value distribution */ 	\
  s, 				/* score */ 				\
  t, 				/* target length */ 			\
  q				/* query length */ 			\
) {									\
  Get_N((N), (el), (q), (t), (evd).K, (evd).H);				\
  (pv) = (evd).K * (N) * exp(-(evd).lambda * (s));			\
  if ((pv) > 0.01) (pv) = 1.0 - exp(-(pv));				\
} /* Evd_pvalue */

/******************************************************************************/
/*
	Get_N

  	Get the effective search space size.  Uses min(q,t)-1 if the 
  	expected alignment length is greater than min(q,t).
*/
/******************************************************************************/
#define Get_N(								\
  N,				/* effective search space */		\
  el,				/* expected alignment length */		\
  q,				/* query length */ 			\
  t, 				/* target length */ 			\
  K,				/* EVD parameter */			\
  H				/* EVD parameter */			\
) { 									\
  if (H) {								\
    (N) = (q) < (t) ? (q) : (t);					\
    (el) = log((K) * (q) * (t))/(H);					\
    if ((el) >= (N))  (el) = (N) - 1;					\
    (N) = ((q) - (el)) * ((t) - (el));					\
  } else {								\
    (N) = (q) * (t);							\
  }									\
} /* get_N */

/******************************************************************************/
/*
  	external prototypes
*/
/******************************************************************************/

/******************************************************************************/
/*
 *         get_n
 *
 *         Get the expected number of scores for use in converting
 *         p-values to E-values.
 *
 *         Returns n, the number of non-outliers.
 *
 */
/******************************************************************************/
extern int get_n(
  SCORE_SET score_set,		/* the set of scores and lens/gcs */
  EVD_SET evd_set		// Distribution.
);

/**********************************************************************/
/* 
	fit_score_distribution

  	Fit a distribution to a set of scores and lengths.
*/
/**********************************************************************/
extern EVD_SET fit_score_distribution(
  DISTR_T dtype,                        // Type of score distribution to fit.
  SCORE_SET score_set,                  /* the set of scores and lens */
  double H,                             /* initial estimate for H */
  int maxiter1,                         /* maximum iterations for ML */
  int maxiter2,                         /* maximum iterations for N-R */
  double eps1,                          /* error tolerance for ML */
  double eps2,                          /* error tolerance for N-R */
  int size,				/* size of length ranges */
  int min_scores,			// minimum number of scores per group
  double lspan 				/* minimum length ratio/group */
);

/**********************************************************************/
/* 

	evd_set_gbin

	Find EVD_SET the bin number of length/gc content

*/
/**********************************************************************/
extern int evd_set_bin(
  double t,					/* target length or GC-content */
  EVD_SET evd_set   /* set of EVDs */
);

/**********************************************************************/
/* 
	evd_set_pvalue

	Compute the p-value of (score, length) pair 
*/
/**********************************************************************/
extern double evd_set_pvalue(
  double s,					/* score */
  double t,					/* target length or GC-content */
  int q, 					/* query length */
  EVD_SET evd_set 				/* set of EVDs */
);

/**********************************************************************/
/*
        set_up_score_set

        Initialize a score_set to contain scores and lengths
        for use in computing the score distribution.
*/
/**********************************************************************/
extern SCORE_SET *set_up_score_set(
  PROB_T p_threshold,           // P-value threshold for motif hits.
  PROB_T dp_threshold,          // Score threshold for DP.
  int    max_gap,               // Maximum gap length to allow in matches.
  bool negatives_only,     // Negatives only in score set?
  MHMM_T* the_hmm               // The HMM itself.
);

#ifdef DEBUG_SCORE
//FIXME CEG
void print_score_set(DISTR_T dtype, SCORE_SET score_set, EVD evd);
#endif

#define FITEVD_H
#endif
