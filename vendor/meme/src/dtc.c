/************************************************************************
*       Copyright                                                       *
*       (1999-2016) The Regents of the University of California.        *
*       All Rights Reserved.                                            *
*       Author: Timothy L. Bailey
************************************************************************/

/***********************************************************************/
/*
        Distance to center statistics routines
*/
/***********************************************************************/

#ifdef MAIN
#define DEFINE_GLOBALS
#endif
#include <stdio.h>
#include "dtc.h"
#include "macros.h"
#include "general.h"
#include "io.h"
#include "logs.h"
#include "message.h"

/* probability distribution */ 
typedef struct distr {
  int maxd;		// maximum center-to-center distance 
  int maxn;		// maximum number of sequences
  double **pdf;		// array of (log) PDFs [nseqs][sum]
  double **cdf;		// array of (log) CDFs [nseqs][sum]
} DISTR;

//static int ndistrs = -1;        /* largest maxd for which distr known */
static int max_maxd = -1;
static int max_maxn = -1;
static DISTR *distrs = NULL;    /* array of distributions for different values of maxd */

static DISTR get_dtc_distr(
  double *dd,			// (log) distribution for distances [0, ..., m]
  int maxd,			// maximum distance
  int maxn			// maximum number of sequences
);

static double *cdf(
  double *d,                            /* integer valued distribution */
  int r                                 /* range [0..r] */
);

/************************************************************************/
/*
        init_dtc_pv_tables

        Initialize the central enrichment p-value tables.

*/
/************************************************************************/
void init_dtc_pv_tables(
  int minw,		// minimum motif width
  int maxw,		// maximum motif width
  int slen,		// sequence length
  int maxn		// maximum number of seqs
)
{
  int nseqs;		/* number of seqs */
  int w, d;

  if (!NO_STATUS)
    fprintf(stderr,
      "Initializing the central enrichment p-value tables for up to %d seqs...\n", maxn);

  // Get the minimum max-distance possible.  The 0.5 makes the not-has-zero case round up, but
  // leaves the other unchanged
  max_maxd = (slen - minw)/2.0 + 0.5;		// narrowest motif; GLOBAL
  int min_maxd = (slen - maxw)/2.0 + 0.5;	// widest motif
  max_maxn = maxn;

  // Create array of distribution objects.
  distrs = NULL;
  Resize(distrs, max_maxd+1, DISTR);

  // Fill in the distribution objects for each value of maxd.
  //for (maxd=min_maxd; maxd<=max_maxd; maxd++) {
  for (w=minw; w<=maxw; w++) {				// motif width
    int maxd = (slen - w)/2.0 + 0.5;			// Maximum distance for this width motif
    bool has_zero = (slen+w) % 2 ? false : true; 		// Is zero distance possible?
    double prob = has_zero ? log(2.0/(2*maxd + 1)) : log(1.0/maxd);	// prob. of any distance

    // Set up underlying distribution for this value of maxd.
    double *dd = NULL;
    Resize(dd, maxd+1, double);
    dd[0] = has_zero ? log(1.0/(2*maxd + 1)) : LOGZERO;		// Zero has 1/2 or no probability
    for (d=1; d<=maxd; d++) dd[d] = prob;

    // Get the PDFs and CDFs for different numbers of sequences for this value of maxd.
    distrs[maxd] = get_dtc_distr(dd, maxd, maxn);
  } // motif width

  if (!NO_STATUS)fprintf(stderr, "\nDone initializing.\n");

} /* init_dtc_pv_tables */

/******************************************************************************/
/*
        get_dtc_distr

        Compute the probability distributions for the total distance 
	between the center of a random site and the center of a sequence.

        Returns a distribution object with PDFs and CDFs.
                pdf[nseqs][sum] = log(Pr(dtc==sum))
                cdf[nseqs][sum] = log(Pr(dtc==sum))

*/
/******************************************************************************/
static DISTR get_dtc_distr(
  double *dd,			// (log) distribution for distances [0, ..., m]
  int maxd,			// maximum distance
  int maxn			// maximum number of sequences
)
{
  int nseqs;		// number sequences
  int sum;		// sum of distances
  int dist;		// center-to-center distance

  // Create the arrays of PDFs and CDFs for different maximum center-to-center distances.
  double **pdfs = NULL;
  Resize(pdfs, maxn+1, double *);
  double **cdfs = NULL;
  Resize(cdfs, maxn+1, double *);

  // Convolve the distribution N times.
  for (nseqs=1; nseqs<=maxn; nseqs++) {
    int old_max_sum = (nseqs-1) * maxd;		// range of distr for nseqs-1 sequences
    int new_max_sum = nseqs * maxd;		// range of distr for nseqs sequences
    pdfs[nseqs] = NULL; 
    Resize(pdfs[nseqs], new_max_sum+1, double);
    if (nseqs == 1) {
      // For one sequence, the distribution is the underlying one.
      for (dist=0; dist<=maxd; dist++) pdfs[nseqs][dist] = dd[dist];
    } else {
      // Set the new distribution to log(0).
      for (sum=0; sum<=new_max_sum; sum++) pdfs[nseqs][sum] = LOGZERO;
      // Now do the convolution.
      for (sum=0; sum<=old_max_sum; sum++) {
	for (dist=0; dist<=maxd; dist++) {
	  double prod = pdfs[nseqs-1][sum] + dd[dist];
	  pdfs[nseqs][sum+dist] = LOGL_SUM(pdfs[nseqs][sum+dist], prod);
	}
      }
    } // nseqs == 1?
    // Compute the CDF from the PDF.
    cdfs[nseqs] = cdf(pdfs[nseqs], new_max_sum);
  }

  // Create a distribution object.
  DISTR distr;
  distr.maxd = maxd;
  distr.maxn = maxn;
  distr.pdf = pdfs;
  distr.cdf = cdfs;

  return(distr);
} // get_dtc_distr

/******************************************************************************/
/*
        cdf

        Compute (log) CDF of an integer-valued (log) distribution.
        Smooths the CDF by linear interpolation so that adjacent positions
        in the table will have different values if possible.
*/
/******************************************************************************/
static double *cdf(
  double *d,                            /* integer valued distribution */
  int r                                 /* range [0..r] */
)
{
  double *cdf=NULL, slope=0;
  int I, i, j, k;

  Resize(cdf, r+1, double);
  cdf[r] = d[r];
  cdf[0] = d[0];
  for (I=1; I<=r; I++) {
    cdf[I] = LOGL_SUM(cdf[I-1], d[I]);
  }

  /* smooth CDF by linear interpolation in logs */
  for (i=0; i<=r; i++) {
    if (d[i] != LOGZERO) break;		// find first non-zero p
  }
  for (; i<=r; i=j) {
    for (j=i+1; j<=r && d[j]==LOGZERO; j++);	/* find next non-zero p */ 
    if (i!=j) slope = (cdf[j]-cdf[i])/(j-i);    /* slope */
    for (k=i+1; k<j; k++) cdf[k] = cdf[i] + (k-i)*slope;
  }

  return cdf;
} /* cdf */


/******************************************************************************/
/*
        get_dtc_pv

        Get the log p-value of the average, weighted distance to the sequence
	center of a randomly chosen site:
                dtc     	average, weighted distance to center
                n       	(weighted) number of sequences
                m	   	maximum possible distance to center
		has_zero	true if zero distance is possible

	If has_zero is true, then the underlying distribution is:
		[0, 1, 2, ..., m]  with p=[1/(2m+1), 2/(2m+1), ..., 2/(2m+1)],
	otherwise,
		[1, 2, ..., m]  with p=[1/m, ..., 1/m].

        Returns the log p-value of dtc.
*/
/******************************************************************************/
double get_dtc_pv(
  double dtc,			// average, weighted distance to center
  double n,			// (weighted) number of sequences
  int maxd,			// maximum possible distance to center
  bool has_zero		// true if zero distance is possible
)
{
  int i;
  double logpv;                         /* log pvalue */

  if (n<1) return 0.0;

  // Use Normal approximation to Bates distribution if n is too big.
  if (n > max_maxn) {
    int a = has_zero ? 0 : 1;				// Minimum/maximum distance distance to center
    int b = maxd;					// rounded up for non-zero case.
    double mu = (a + b)/2.0;				// mean of Bates distribution
    double sd = sqrt( ((b-a)*(b-a) ) / (12.0*n) );	// standard deviation of Bates distribution
    double z = (dtc - mu) / sd;				// Z-score of dtc
    double bates_pv = 0.5*erfc(-z/sqrt(2.0));		// Normal approximation to Pr(avg_dist <= x)
    return(bates_pv == 0 ? LOGZERO : log(bates_pv));	// return log p-value

  } else {
    // Use the precomputed (exact) table of p-values.
    /* return log of geometric mean of pv if n is not integral */
    double n0, n1;                        /* floor and ceil of n */
    if ( (n0=floor(n)) != (n1=ceil(n)) ) {
printf("n0 %d n1 %f\n", max_maxn, n);exit(0);
      return ( 
	(n1-n)*get_dtc_pv(dtc, n0, maxd, has_zero) +
	(n-n0)*get_dtc_pv(dtc, n1, maxd, has_zero)
      );
    }
    double I = dtc * n; 				// weighted sum of distances
    int N = (int) n;				// make number of sequences an integer
    int I0 = (int) I;					/* floor of position */
    int I1 = I0 + 1;					/* ceil. of position */
    if (I < 0) {					/* lower bound I */
      logpv = distrs[maxd].cdf[N][0];
    } else if (I0 >= N*maxd) {				/* upper bound I */
      logpv = distrs[maxd].cdf[N][N*maxd];
    } else {                                            /* lin. interpolate */
printf("dtc %f N % d maxd %d I %f I0 %d I1 %d N*maxd %d\n", dtc, N, maxd, I, I0, I1, N*maxd);
      logpv = distrs[maxd].cdf[N][I0] + (I-I0)*(distrs[maxd].cdf[N][I1]-distrs[maxd].cdf[N][I0]);
    }
    /* return log p-value */
    return(logpv);
  }
} /* get_dtc_pv */

#ifdef MAIN
/************************************************************************/
/*
        dtc <minw> <maxw> <slen> <maxn>

        Compute the probability distribution for the average distance 
	between a random site and a the center of <maxn> sequences of 
	length <slen>.
*/
/************************************************************************/
int main(
  int argc,
  char** argv
)
{
  int minw=1, maxw=1, slen=2, minn=1, maxn=1;
  double accuracy = -1;

  int i = 1;
  argv[0] = "dtc";
  DO_STANDARD_COMMAND_LINE(2,
    USAGE(<minw> <maxw> <slen> <minn> <maxn> [options]);
    USAGE(\n\t<minw>\tmininum motif width); 
    USAGE(\t<maxw>\tmaxinum motif width); 
    USAGE(\t<slen>\tlength of sequences); 
    USAGE(\t<minn>\tnumber of sequences); 
    USAGE(\t<maxn>\tnumber of sequences); 
    DATA_OPTN(1, minsum, <acc>, \tprint smallest sum where normal is accurate,
      accuracy = atof(_OPTION_)); 
    NON_SWITCH(1,\r,
      switch (i++) {
	case 1: minw = atoi(_OPTION_); break;
	case 2: maxw = atoi(_OPTION_); break;
	case 3: slen = atoi(_OPTION_); break;
	case 4: minn = atoi(_OPTION_); break;
	case 5: maxn = atoi(_OPTION_); break;
        default: COMMAND_LINE_ERROR;
      }
    );
    USAGE(\n\tCompute the probability distribution for the average);
    USAGE(\tdistance between the center of a random site and the)
    USAGE(\t center of a sequence of length <slen> for up to <nseq> sequences.);
    USAGE(\n\tCopyright);
    USAGE(\t(2016) The Regents of the University of California);
    USAGE(\tAll Rights Reserved.);
    USAGE(\tAuthor: Timothy L. Bailey);
  );
  if (i <= 5) { fprintf(stderr, "You need 5 command line arguments.\n"); exit(1);}

  double m1=0, e1=0, m2=0, e2=0;

  init_log();
  init_exp();

  // Initialize the distributions.
  //init_dtc_pv_tables(minw, maxw, slen, maxn);
  init_dtc_pv_tables(minw, maxw, slen, 5);


  // Print the distributions.
  int w;
  for (w=minw; w<=maxw; w++) {
    int maxd = (slen - w)/2.0 + 0.5;		// Maximum distance for this width motif, rounded up for non-zero case.
    bool has_zero = (slen + w) % 2 ? false : true;

    printf("#  w  slen  maxd   N   sum         p     cdf     AvgDist   Bates   Ratio\n");

    /* get distribution and print it */
    int n;
    for (n=minn; n<=maxn; n++) { 		// number of sequences
      int sum;
      double x=0, bates_pv=0, log_pv=0, ratio=0;
      for (sum=n*maxd; sum >= 0; sum--) { 
        x = ((double) sum)/n;			// average distance to center
        bates_pv = get_dtc_pv(x, n, maxd, has_zero);
        log_pv = distrs[maxd].cdf[n][sum];
        ratio = exp(bates_pv - log_pv);
        if (accuracy != -1 && fabs(ratio - 1) >= accuracy) break;
	if (distrs[maxd].pdf[n][sum] == LOGZERO) {
	  m1 = e1 = 0;
	} else {
	  exp10_logx(distrs[maxd].pdf[n][sum]/log(10.0), m1, e1, 1);
	}
	if (distrs[maxd].cdf[n][sum] == LOGZERO) {
	  m2 = e2 = 0;
	} else {
	  exp10_logx(log_pv/log(10.0), m2, e2, 1);
	}
        if (accuracy == -1 && log_pv != LOGZERO) 
	  printf("%4d %4d %4d %4d %4d    %3.1fe%+05.0f %3.1fe%+05.0f  %.2f %.1e   %.1g\n", 
  	  w, slen, maxd, n, sum, m1, e1, m2, e2, x, exp(bates_pv), ratio);
      } // sum
      if (accuracy != -1) 
        printf("%4d %4d %4d %4d %4d    %3.1fe%+05.0f %3.1fe%+05.0f    %.2f %.1e   %.1f\n", 
          w, slen, maxd, n, sum, m1, e1, m2, e2, x, exp(bates_pv), ratio);
    } // n
  } // w

  return 0;
}
#endif /* MAIN */
