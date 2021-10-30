
#ifndef MOTIF_IN_H
#define MOTIF_IN_H

#include "motif.h"
#include "alphabet.h"
#include "motif-in-flags.h"
#include "array-list.h"


typedef struct mread MREAD_T;


/*
 *  Create a motif reader to read in a motif
 *  file of any MEME format. Note that the filename
 *  is not used for reading but is used for
 *  better error messages and hinting the
 *  file type.
 */
MREAD_T* mread_create(const char *filename, int options, bool symmetrical);

/*
 * Free all memory associated with the motif reader
 */
void mread_destroy(MREAD_T *mread);

/*
 * Update the readers with another chunk of the file.
 */
void mread_update(MREAD_T *mread, const char *buffer, size_t size, short end);

/*
 * Loads all the currently buffered motifs into a list.
 * If the file is set then this will read all the motifs in the
 * file into the list. If a list is not passed then
 * it will create a new one. 
 * returns the list.
 */
ARRAYLST_T* mread_load(MREAD_T *mread, ARRAYLST_T *motifs);

/*
 * Is there another motif to return
 */
bool mread_has_motif(MREAD_T *mread);

/*
 * Get the next motif. If there is no motif ready then
 * returns null.
 */
MOTIF_T* mread_next_motif(MREAD_T *mread);

/*
 * Convert the motifs to the given alphabet (or die trying).
 *
 * This cancels any previously set background option.
 */
void mread_set_conversion(MREAD_T *mread, ALPH_T *alph, const ARRAY_T *bg);

/*
 * Set the background to be used for pseudocounts.
 */
void mread_set_background(MREAD_T *mread, const ARRAY_T *bg, ALPH_T *alph);

/*
 * Set the background source to be used for pseudocounts.
 * The background is resolved by the following rules:
 * - If source is:
 *   null or "--nrdb--"       		use nrdb frequencies
 *   "--uniform--"            		use uniform frequencies
 *   "motif-file", --motif--, --query-- use frequencies in motif file
 *   file name                		read frequencies from bg file
 */
void mread_set_bg_source(MREAD_T *mread, const char *source, ALPH_T *alph);

/*
 * Set the pseudocount to be applied to the motifs.
 * Uses the motif background unless a background has been set.
 */
void mread_set_pseudocount(MREAD_T *mread, double pseudocount);

/*
 * Set the trimming threshold to be applied to the motifs
 */
void mread_set_trim(MREAD_T *mread, double trim_bits);

/* Get the background source. */
char* mread_get_other_bg_src(MREAD_T *mread);

/*
 * Get the alphabet. If the alphabet is unknown then returns
 * NULL.
 */
ALPH_T* mread_get_alphabet(MREAD_T *mread);

/*
 * Get the strands. 
 * DNA motifs have 1 or 2 strands,
 * otherwise return 0 (for both protein and unknown)
 */
int mread_get_strands(MREAD_T *mread);

/*
 * Get the background. If the background is unavailable
 * then returns null.
 */
ARRAY_T *mread_get_background(MREAD_T *mread);

/*
 * Get the motif occurrences. If the motif occurrences are unavailable
 * then returns null. This function only returns non-null once and
 * if called after that it will throw a fatal error.
 * In MEME text format the combined sites are at the
 * end of the file so you should read all the motifs before
 * calling this function.
 */
ARRAYLST_T *mread_get_motif_occurrences(MREAD_T *mread);

/*
 * Destroy the motif occurrences.
 */
void destroy_motif_occurrences(ARRAYLST_T *motif_occurrences);

#endif
