/***********************************************************************
 * FILE: order.h
 * AUTHOR: William Stafford Noble
 * PROJECT: MHMM
 * COPYRIGHT: 2001-2008, WSN
 * DESCRIPTION: Data structure for representing the order and spacing
 *              of motifs in a linear model.
 ***********************************************************************/
#ifndef ORDER_H
#define ORDER_H

#include "motif.h"
#include "red-black-tree.h"

/*
 * The motifs are indexed from 1 to n.  In motif_spacers, the value at
 * position i is the length of the spacer following motif i, so
 * position 0 contains the length of the spacer before the first
 * motif.
 * */
typedef struct order ORDER_T;

/***********************************************************************
 * Create an order-and-spacing object from a user-specified string.
 *
 * The given string specifies the order and spacing of the motifs
 * within the model, and has the format "l=n=l=n=...=l=n=l", where "l"
 * is the length of a region between motifs, and "n" is a motif index.
 * Thus, for example, the string "34=3=17=2=5" specifies a two-motif
 * linear model, with motifs 3 and 2 separated by 17 letters and
 * flanked by 34 letters and 5 letters on the left and right.
 *
 * Returns the ORDER_T object, or NULL if the null string was given.
 ***********************************************************************/
ORDER_T* create_order
  (const char* order_string);

/***********************************************************************
 * Create an order-and-spacing object from a user-specified string.
 * <spacer>=<strand><motif ID>=<spacer>=<strand><motif ID>=<spacer>
 * eg
 * 5=+m1=12=-m1=15
 ***********************************************************************/
ORDER_T* create_order_from_ids
  (const char* order_string, RBTREE_T *ids2nums);

/***********************************************************************
 * Create an empty order-and-spacing object.
 ***********************************************************************/
ORDER_T* create_empty_order
  (int    num_occurs,
   double sequence_p);

/***********************************************************************
 * Add one spacer and the following motif to an order object.
 ***********************************************************************/
void add_occurrence
  (int    motif_i,
   int      spacer_length,
   ORDER_T* order);

/***********************************************************************
 * Determine whether a given motif order-and-spacing object contains a
 * given motif index.
 ***********************************************************************/
bool order_contains
  (int    motif_i,
   ORDER_T* order);

/***********************************************************************
 * Get the number of motif occurences 
 ***********************************************************************/
int get_order_occurs
  (ORDER_T* order);

/***********************************************************************
 * Get the length of the spacer after a given motif.  Motifs are
 * indexed from 1, so 0 gets the initial spacer.
 ***********************************************************************/
int get_spacer_length
  (ORDER_T* order,
   int      index);

/***********************************************************************
 * Get the motif number at a position.
 ***********************************************************************/
int get_order_motif
  (ORDER_T* order,
   int      index);

/***********************************************************************
 * Get and set the pvalue of a given order object.
 ***********************************************************************/
void set_pvalue
  (double pvalue,
   ORDER_T* order);

double get_pvalue
  (ORDER_T* order);

/***********************************************************************
 * Get the number of distinct motif occurrences.
 ***********************************************************************/
int get_num_distinct
  (ORDER_T* order);

/***********************************************************************
 * Re-order the motifs in the order in which they will appear in the model.
 *
 * The motifs are re-ordered according to the order requested in the
 * order object.
 *
 * order - The requested motif order and spacing.
 *
 * num_motifs - The number of motifs.
 *
 * motifs (input/output) - The motifs themselves, re-ordered according
 *        to order_string.
 ***********************************************************************/
/*
void reorder_motifs
  (ORDER_T* order_string,
   int*     num_motifs,
   MOTIF_T* motifs);
*/

/***********************************************************************
 * Print the order and spacing with hyphens between.
 ***********************************************************************/
void print_order_and_spacing
  (FILE*    outfile,
   ORDER_T* order_spacing);

/***********************************************************************
 * Free dynamic memory used by an order object.
 ***********************************************************************/
void free_order
  (ORDER_T* order);

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

