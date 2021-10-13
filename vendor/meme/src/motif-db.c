#include <errno.h>
#include <sys/stat.h>
#include <time.h>
#include <fnmatch.h>
#include "motif-db.h"
#include "fasta-io.h"

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
) {
  MOTIF_DB_T *db;
  struct stat stbuf;
  time_t last_mod;
  char *name;

  assert(path != NULL);

  // Check that only one file uses stdin.
  if (strcmp(path, "-") == 0 && *stdin_used) {
      die("%s database file '%d' specifies standard input when it is "
	"already in use for another file.\n", type, id);
  }

  // Check that the path exists and get last modified time.
  if (strcmp(path, "-") == 0) {
    name = strdup("standard input");
    last_mod = time(NULL);
    *stdin_used = true;
  } else {
    // check if the path exists and determine last modified time
    errno = 0;
    if (stat(path, &stbuf) == -1) {
      if (errno == ENOENT) { // File not found
	die("Database file '%s' does not exist.\n", path);
      } else { // stat failed for some other reason
	die("Unable to check for status of file '%s'.\nError: %s.\n",
	  path, strerror(errno));
      }
    }
    last_mod = stbuf.st_mtime;
    // get the file name from the path
    name = file_name_from_path(
      path, 
      remove_ext,		// remove the extension
      remove_meme_ext, 		// remove the extension if == ".meme"
      remove_underscores	// replace underscores
    );
  }
  assert(name != NULL);

  db = mm_malloc(sizeof(MOTIF_DB_T));
  memset(db, 0, sizeof(MOTIF_DB_T));
  db->id = id;
  db->source = strdup(path);
  db->name = name;
  db->last_mod = last_mod;
  db->motifs = NULL;
  db->loaded = 0;
  db->excluded = 0;
  // The strings that Tomtom will be using are fixed size buffers of the MOTIF_T object
  // and they won't be deallocated in the life of the database entry so I don't need
  // a key copy or key free function.
  db->matched_motifs = rbtree_create(rbtree_strcmp, NULL, NULL, NULL, NULL);

  return db;
} // create_motif_db

/**************************************************************************
 * Used by destroy_motif_db and destroy_motif_db_not_motifs.
 **************************************************************************/
static void destroy_motif_db2(void *ptr, bool save_motifs) {
  MOTIF_DB_T *db = (MOTIF_DB_T *)ptr;
  free(db->source);
  free(db->name);
  if (save_motifs) {
    arraylst_destroy(NULL, db->motifs);
  } else {
    arraylst_destroy(destroy_motif, db->motifs);
  }
  if (db->matched_motifs) rbtree_destroy(db->matched_motifs);
  memset(db, 0, sizeof(MOTIF_DB_T));
  free(db);
} // destroy_motif_db

/**************************************************************************
 * Destroy the motif db.
 **************************************************************************/
void destroy_motif_db(void *ptr) {
  destroy_motif_db2(ptr, false);
}

/**************************************************************************
 * Destroy the motif db but not the motifs.
 **************************************************************************/
void destroy_motif_db_not_motifs(void *ptr) {
  destroy_motif_db2(ptr, true);
}

/******************************************************************
 *
 * This function tests if a given string matches any pattern
 * in a list of UNIX filename expressions.
 *
******************************************************************/
static inline bool pattern_list_match(
  char *motif_name,
  ARRAYLST_T *patterns
) {
  int i;
  for (i=0; i<arraylst_size(patterns); i++) {
    char *name_pattern = (char*)arraylst_get(i, patterns);
    if (fnmatch(name_pattern, motif_name, 0) == 0) return true;
  }
  return false;
}

/******************************************************************
 *
 * This function reads in (selected) motifs and sets
 * the background model for them.  It returns a motif database.
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
  bool scan_separately,		// append RC copies of motifs at end of motif list
  double pseudocount,		// total counts added to motif PSPM column
  bool set_trim,		// set trimming of motifs
  double trim_bit_threshold,	// trimming threshold (used by SpaMo)
  char **bg_source,		// (IN/OUT) name of background file, special string, 
				// or NULL -> from sequences
  bool bg_symmetrical,		// make background symmetrical if alphabet complementable
  ARRAY_T **bg,			// (IN/OUT) the background model
  char *seq_source,		// name of sequence file
  ALPH_T *alphabet, 		// sequence (or query motif) alphabet
  bool set_conv_alph,		// set conversion alphabet
  bool remove_ext,		// remove extension from name
  bool remove_meme_ext,		// remove ".meme" extension from name
  bool remove_underscores,	// remove underscores from name
  bool *stdin_used		// IN/OUT check and set if path is "-"
) {

  // Create the motif database struct to return.
  // Check that stdin is not used twice.
  MOTIF_DB_T *motif_db = create_motif_db(
    db_index, 
    motif_source, 
    type,		// type of DB for error messages
    remove_ext,
    remove_meme_ext,
    remove_underscores,
    stdin_used
  );

  // Create a motif reader.
  MREAD_T *mread = mread_create(motif_source, OPEN_MFILE, bg_symmetrical);

  // Initialize the pseudocount.  This just sets mread->pseudo_total.
  mread_set_pseudocount(mread, pseudocount);

  // Get the alphabet from the motifs if it is unknown.
  ALPH_T *db_alph = alph_hold(mread_get_alphabet(mread));
  if (alphabet == NULL) {			// no alphabet given
    alphabet = db_alph;				// use alphabet of motifs
  } else {
    if (! alph_equal(alphabet, db_alph)) {
      if (set_conv_alph) {
	// The alphabet does not match but we have been told to convert it so
	// we must check that such a thing is possible.
	switch(alph_core_subset(db_alph, alphabet)) {
	  case 0:
	    die("The motifs in '%s' have the '%s' alphabet which is not "
	      "a subset of the '%s' alphabet.\n", 
	      motif_source, alph_name(db_alph), alph_name(alphabet));
	    break;
	  case -1:
	    fprintf(stderr, "Warning: the alphabet expansion from '%s' to '%s'"
	      " requires changing complementation rules.\n",
	      alph_name(db_alph), alph_name(alphabet));
	    break;
	}
	mread_set_conversion(mread, alphabet, *bg);
      } else {
	die("Motif file '%s' uses a '%s' alphabet when a '%s' alphabet was expected.\n",
	  motif_source, alph_name(db_alph), alph_name(alphabet));
      }
    } else {		// alphabets are equal; clear set_conv_alph
      set_conv_alph = false;
    } 
  }
 
  // Get the background if unknown and initialize the motif reader with it.
  //if (get_background) {
  if (*bg == NULL) {		// background unknown
    if (*bg_source) {		// background from bgfile, or motifs, or uniform
      // Background from bgfile, or motifs, or uniform.
      // This frees mread->other_bg, other_bg_src, conv_alph; copies other_bg_src;
      // and then it sets other_bg_src=*bg_source and calls set_pseudo_bg(mread, NULL).
      // set_pseudo_bg(mread, NULL) calls get_alphabet(mread->formats->data) to get alph.
      //mread_set_bg_source(mread, *bg_source, db_alph);
      mread_set_bg_source(mread, *bg_source, alphabet);
      *bg = mread_get_background(mread);
      *bg_source = mread_get_other_bg_src(mread);
    } else {                            // background from sequences
      *bg = calc_bg_from_file(alphabet, seq_source, bg_symmetrical, false);
      mread_set_background(mread, *bg, alphabet);
      *bg_source = strdup("--sequences--");
    }
  } else if (*bg != NULL) {		// background known
    mread_set_background(mread, *bg, alphabet); // Use provided background model.
  }

  // Set the conversion again if needed (since it got cleared above).
  if (set_conv_alph) mread_set_conversion(mread, alphabet, *bg);

  // Activate motif trimming if requested.
  if (set_trim) mread_set_trim(mread, trim_bit_threshold);

  // Read the (selected) motifs.
  int motifs_read = 0;
  ARRAYLST_T *motifs = arraylst_create();
  bool check_include_names = include_names && rbtree_size(include_names) > 0;
  bool check_include_indices = include_indices && rbtree_size(include_indices) > 0;
  bool check_include_patterns = include_patterns && arraylst_size(include_patterns) > 0;
  bool check_exclude_patterns = exclude_patterns && arraylst_size(exclude_patterns) > 0;
  while (mread_has_motif(mread)) {
    motifs_read++;
    MOTIF_T *motif = mread_next_motif(mread);
    if (check_include_names || check_include_indices) {
      int motif_index = get_motif_idx(motif);

      if (include_names && rbtree_find(include_names, get_motif_id(motif))) {
	// keep motif as name is selected
      } else if (include_indices && rbtree_find(include_indices, &motif_index)) {
        // keep motif as index is selected
      } else {
        // not among the selected names or indices
        goto destroy_motif;
      }
    }
    if (check_include_patterns && !(check_include_names || check_include_indices)) {
      if (! pattern_list_match(get_motif_id(motif), include_patterns)) {
	DEBUG_FMT(NORMAL_VERBOSE, "Discarding motif '%s' in '%s' because it did not match any -inc pattern.\n",
	  get_motif_id(motif), motif_source);
        goto destroy_motif;
      }
    }
    if (check_exclude_patterns) {
      if (pattern_list_match(get_motif_id(motif), exclude_patterns)) {
        DEBUG_FMT(NORMAL_VERBOSE, "Discarding motif '%s' in '%s' it matched an -exc pattern.\n",
          get_motif_id(motif), motif_source);
        goto destroy_motif;
      }
    }

    // Got past all the exclusion tests, so save the motif.
    arraylst_add(motif, motifs);
    continue;
destroy_motif:
    // Motif was not selected or was excluded.
    destroy_motif(motif);
  } // motifs

  // Check that requested motifs found.
  int motifs_found = arraylst_size(motifs);
  if (motifs_found==0) {
    if (include_names && rbtree_size(include_names) == 1) {
      RBNODE_T *node = rbtree_first(include_names);
      char *motif_name = (char *)rbtree_value(node);
      die("Requested motif '%s' was not found in file '%s'.\n", motif_name, motif_source);
    } else if (include_indices && rbtree_size(include_indices) == 1) {
      RBNODE_T *node = rbtree_first(include_indices);
      int motif_index = *(int *)rbtree_value(node);
      die("Requested motif number %d  was not found in file '%s'.\n", motif_index, motif_source);
    } else if (check_include_patterns) {
      die("No motifs with names matching the %d given patterns were found in file '%s'.\n", 
        arraylst_size(include_patterns), motif_source);
    } else if (check_exclude_patterns) {
      fprintf(stderr, "Warning: Motif file '%s' contains no valid motifs or all were excluded.\n", motif_source);
    } else {
      fprintf(stderr, "Warning: Motif file '%s' contains no valid motifs.\n", motif_source);
    }
  } else if (include_patterns && motifs_found < arraylst_size(include_patterns)) {
    DEBUG_FMT(NORMAL_VERBOSE, "Only %d motif names in file '%s' matched the %d patterns you specified.\n",
      motifs_found, motif_source, arraylst_size(include_patterns));
  }

  // Remove duplicate ID motifs and motifs with zeros (if requested).
  if (motifs_found > 1) {
    int i;
    RBTREE_T *seen = rbtree_create(rbtree_strcmp, rbtree_strcpy, free, NULL, NULL);
    for (i = 0; i < motifs_found; i++) {
      MOTIF_T *motif = (MOTIF_T*)arraylst_get(i, motifs);
      bool exclude = false;
      if (!rbtree_make(seen, get_motif_id(motif), NULL)) {
        exclude = true;
	DEBUG_FMT(NORMAL_VERBOSE, "Discarding motif '%s' in file '%s' (non-unique ID).\n",
	  get_motif_id(motif), motif_source);
      } else if (!allow_zeros && has_motif_zeros(motif)) {
        exclude = true;
	DEBUG_FMT(NORMAL_VERBOSE, 
          "Discarding motif '%s' in file '%s' because it contains probabilities of zero\n"
          "in some columns; use option -motif-pseudo with a non-zero value.\n",
          get_motif_id(motif), motif_source);
      }
      if (exclude) {
	arraylst_remove(i--, motifs);
	motifs_found--; 
	destroy_motif(motif);
      }
    }
    rbtree_destroy(seen);
  }

  // Check that at least one motif was found.
  if (motifs_found == 0) {
    DEBUG_FMT(NORMAL_VERBOSE, 
     "Warning: Motif file '%s' contains no loadable motifs.\n", motif_source);
  }

  // Append reverse complement motifs to list if requested.
  if (scan_separately) add_reverse_complements(motifs);

  // Resize the array list.
  arraylst_fit(motifs);

  // Check that the motif alphabet is (now) equal to the sequence alphabet.
  if (motifs_found > 1) {
    if (!alph_equal(alphabet, get_motif_alph((MOTIF_T*)arraylst_peek(motifs)))) {
      die("Expected motifs in alphabet '%s' in file '%s'.\n", alph_name(alphabet),
        motif_source);
    }
  }

  // Return the background if it was not input.
  if (motifs_found > 1 && *bg == NULL) *bg = mread_get_background(mread);

  // Cleanup the motif reader.
  mread_destroy(mread);

  // Set the motifs and number of motifs found or excluded.
  motif_db->motifs = motifs;
  motif_db->list_entries = motifs_found;
  motif_db->loaded = motifs_found;
  motif_db->excluded = motifs_read - motifs_found;

  return(motif_db);
} // read_motifs_and_background

/**************************************************************************
 * Clumps together multiple motifs that have the same warning message.
 * Used by load_secondary_motifs.
 **************************************************************************/
inline void clump_motif_db_warning(int *current, int type, char *msg, char *file, char *motif_id) {
  if (*current != type) {//clump warnings of the same type
    if (verbosity >= NORMAL_VERBOSE) {
      //add a new line after a warning of a previous type
      if (*current)  fprintf(stderr, "\n");
      //output the message
      fprintf(stderr, msg, file);
    }
    *current = type;
  }
  if (verbosity >= NORMAL_VERBOSE) 
    fprintf(stderr, " %s", motif_id);
}
