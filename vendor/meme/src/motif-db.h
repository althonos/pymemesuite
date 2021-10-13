#ifndef MOTIF_DB_H
#define MOTIF_DB_H
#include "array.h"
#include "alphabet.h"
#include "motif-in.h"

#define NO_WARNING 0
#define DUPLICATE_WARNING 1
#define MARGIN_WARNING 2

/**************************************************************************
 * Clump warning messages of the same type together.
 **************************************************************************/
inline void clump_motif_db_warning(int *current, int type, char *msg, char *file, char *motif_id);

/**************************************************************************
 * An object for storing information about a file containing a set of motifs
 **************************************************************************/
typedef struct {
  int id; 		// An id assigned to the database based on the order it was listed
  char *source;		// full path of motif file
  char *name;		// motif file name
  time_t last_mod;	// last modified date
  ARRAYLST_T *motifs;	// the motifs (may be NULL)
  int loaded;		// number of motifs loaded
  int excluded;		// number of motifs excluded
  // The following entries are used by Tomtom.
  int list_index;	// The index that in the list of motifs that this db starts at
  int list_entries;	//the number of entries in the list of motifs from this db
  RBTREE_T *matched_motifs;
} MOTIF_DB_T;

/**************************************************************************
 * Create a new motif database.
 * Makes a copy of the strings so there is no confusion on what can be
 * freed later.
 **************************************************************************/
MOTIF_DB_T* create_motif_db(
  int id,		// label to give DB
  char *path,		// path to motif file
  char *type,		// type of DB for error messages
  bool remove_ext,	// remove extension from name
  bool remove_meme_ext,	// remove ".meme" extension from name
  bool remove_underscores,	// remove underscores from name
  bool *stdin_used	// IN/OUT check and set if path is "-"
);

/**************************************************************************
 * Destroy the motif db.
 **************************************************************************/
void destroy_motif_db(void *motif_db);

/**************************************************************************
 * Destroy the motif db but not the motifs.
 **************************************************************************/
void destroy_motif_db_not_motifs(void *motif_db);

/******************************************************************
 *
 * This function reads in (selected) motifs and sets
 * the background model for them.
 *
******************************************************************/
MOTIF_DB_T *read_motifs_and_background(
  int db_index,			// (integer) name to give motif DB
  char *motif_source,		// name of motif file
  char *type,			// type of DB for error messages
  RBTREE_T *include_names,	// get only this set of motifs (by name, or NULL)
  RBTREE_T *include_indices,	// get only this set of motifs (by index, starts at 1, or <0)
  ARRAYLST_T *include_patterns,	// get only this set of motifs only (by filename pattern, or NULL)
				// Note: If include_names or include_indices contain entries,
				// then include_patterns will be ignored.
  ARRAYLST_T *exclude_patterns,	// exclude this set of motifs (by filename pattern, or NULL)
  bool allow_zeros,		// allow motifs with zero probability entries
  bool scan_separately,		// create RC copies of motifs, appended
  double pseudocount,		// total counts added to motif PSPM column
  bool set_trim,		// set trimming of motifs
  double trim_bit_threshold,	// trimming threshold (used by SpaMo)
  char **bg_source,		// (IN/OUT) name of background file, special string, 
				// or NULL -> from sequences
  bool bg_symmetrical,		// make background symmetrical if alphabet complementable
  ARRAY_T **bg,			// (IN/OUT) the background model
  char *seq_source,		// name of sequence file
  ALPH_T *alphabet, 		// sequence alphabet
  bool set_conv_alph,		// set conversion alphabet
  bool remove_ext,		// remove extension from name
  bool remove_meme_ext,		// remove ".meme" extension from name
  bool remove_underscores,	// remove underscores from name
  bool *stdin_used		// IN/OUT check and set if path is "-"
);

#endif
