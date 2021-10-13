/* user.h */
/*
	User settable parameters
*/

#ifndef user_h
#define user_h

#define MSN	 24 		/* maximum length of sample name */
				/* MSN + 40 < PAGEWIDTH (see meme.h) */
#define MAXALPH  28		/* maximum length of alphabet + 1 for 'X' */
#define MAXSITE 300		/* maximum length of a site */
#define MINSITES 2		/* minimum number of sites in valid motif */
#define LLR_RANGE 200		/* range of scaled LLR statistic */

#define MINCONS 0.2		/* Display 'X' as consensus if no letter f > */ 
#define LOGOHEIGHT 7.5		// height of sequence logo in cm.
#define MAXLOGOWIDTH 30		// maximum width of sequence logo in cm.

/* minimum allowable motif width before shortening; 
   never make less than 2 or will crash! */
#define MIN_W 8
/* maximum allowable length before shortening */
#define MAX_W 50

/* maximum number of MEME motifs to find if -evt set */
#define MAX_NMOTIFS 1000

/* default size of heap for branching search */
// Don't really need the heaps now that branching search is removed.
#define HSIZE 64
#define HS_DECREASE 2

/* default branching factor for branching search */
#define BFACTOR 3

/* Amount of error tolerated in probability column sums (they should sum to
   approximately 1): */
#define ERR_EPSILON 0.01

// Tomtom memory footprint is *cubic* in the query width.
// memory(w) = (w*(w+1)/2) * ((100*w)+1) * 8 * 2
// memory(100) = 0.8G
// memory(120) = 1.39G
// memory(150) = 2.8G
#define TOMTOM_MAX_QUERY_WIDTH 100

#endif
