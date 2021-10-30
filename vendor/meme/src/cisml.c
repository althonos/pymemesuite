/********************************************************************
 * FILE: cisml.c
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 9/25/2007
 * PROJECT: MEME suite
 * COPYRIGHT: 2007, UW
 ********************************************************************/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include "cisml.h"
#include "cisml-sax.h"
#include "dir.h"
#include "fasta-io.h"
#include "heap.h"
#include "io.h"
#include "string-builder.h"
#include "qvalue.h"
#include "xml-out.h"
#include "xml-util.h"
#include "cisml-dtd.h"
#include "utils.h"

#include "array.h"
const int PATTERN_INCREMENT = 5;
const int SEQUENCE_INCREMENT = 50;
const int ELEMENT_INCREMENT = 500;

// Define data structures related to CisML format.

struct multi_pattern {
  MULTI_PATTERN_MATCH_T *match; // May be NULL
  double *score; // May be NULL
  double *pvalue; // May be NULL
  int num_patterns; // Count of child patterns in this multi-pattern
  int num_allocated_patterns; // Size of pattern array
  PATTERN_T **patterns; // Array of child patterns in this multi-pattern
};

struct multi_pattern_match {
  char *clusterid; // Required
  char *seq_name; // Required
  char *sequence; // May be null
  int start; // Required
  int stop; // Required
  double evalue;
  double qvalue;
};

struct pattern {

  char *accession; // Required.
  char *name; // Required.

  double *pvalue; // May be NULL.
  double *score; // May be NULL.
  char *db; // May be NULL.
  char *lsid; // May be NULL.

  int num_allocated_sequences; // Size of child scanned-sequence array.
  int num_allocated_elements; // Size of matched-element array.
  int num_sequences; // Count of child scanned-sequences.
  long num_scanned_positions; // Count of all positions scanned for
                            // matched_element.
  int num_stored_matches; // Count of the matched-elements stored in the elements array.
  int max_stored_matches;   // maximum number of matches to store for this pattern.
  double max_pvalue_retained; // Largest pvalue of retained records.

  bool qvalues_computed; // Have q-values been calcuatled for matched-elements
  bool is_complete;     /// All matched elements have been added to pattern
  SCANNED_SEQUENCE_T **sequences; // Array of child scanned-sequence pointers.
  HEAP *element_heap; // Heap of matched elements ordered by ascending score.
  MATCHED_ELEMENT_T **elements; // Array of matched element pointers, ordered by p-value
};

struct scanned_sequence {

  char* accession; // Required.
  char* name; // Required.

  double *pvalue; // May be NULL.
  double *score; // May be NULL.
  int *length; // May be NULL.
  char* db; // May be NULL.
  char* lsid; // May be NULL.

  long num_scanned_positions; // Count of all positions scanned for
                            // matched_elements.
  int num_matched_elements; // Count of all matched_elements
                            // <= num_scanned_elements because of filtering.
  int num_allocated_elements; // Number of elements in elements array.

  MATCHED_ELEMENT_T **elements; // Array of all matched elements for this sequence
  PATTERN_T *parent_pattern; // Pointer to containing pattern.

};

struct matched_element {

  int start; // Required.
  int stop; // Required.

  double score;
  bool has_score;
  double pvalue;
  bool has_pvalue;
  double qvalue;
  bool has_qvalue;
  char* clusterid; // May be NULL.
  char* sequence; // May be NULL.
  char strand;

  SCANNED_SEQUENCE_T *parent_sequence; // Pointer to containing scanned-sequence.

};

struct cisml {

  char *program_name; // Required
  char *command_line; // Required
  char *pattern_file; // Required
  char *sequence_file; // Required

  char *background_file; // May return NULL
  double *pattern_pvalue_cutoff; // May return NULL
  double *sequence_pvalue_cutoff; // May return NULL
  double *site_pvalue_cutoff; // May return NULL
  double *site_qvalue_cutoff; // May return NULL
  int num_passing_cutoff; // Number of matches passing significance cuttoff
  char* sequence_filter; // May return NULL

  int num_allocated_multi_patterns; // Size of multi-pattern array.
  int num_allocated_patterns; //Size of patern array.
  int num_multi_patterns; // Count of multi-pattern objects.
  int num_patterns; // Count of pattern objects.

  MULTI_PATTERN_T **multi_patterns; // Array of pointers to multi-pattern objects.
  PATTERN_T **patterns; // Array of pointers to pattern objects.

};

struct cisml_match_it {
  CISML_T *cisml;
  int num_matches;            // Total number of matched elements in cisml
  int num_matches_returned;   // Number of matches returned by this iterator so far
  int *pattern_match_indices; // Current position in each pattern's array of matches
  int *pattern_match_limits;  // Max position in each pattern's array of matches
};

// Forward declarations
static int compare_matched_elements(void *p1, void *p2);
static void *copy_matched_element(void *p);
static void destroy_matched_element(void *p);
static void reduce_pattern_matched_elements(PATTERN_T *pattern);

/**********************************************************************
  allocate_cisml

  Constructor for the cisml data structure.
  Sets required fields to point to copies of the provided arguments.
  Other fields set to 0 or NULL.
**********************************************************************/
CISML_T *allocate_cisml(
  const char *program_name,
  const char *command_line,
  const char *pattern_file,
  const char *sequence_file
) {

  assert(program_name != NULL);
  assert(command_line != NULL);
  assert(pattern_file != NULL);
  assert(sequence_file != NULL);

  // Allocate memory and initialze fields
  CISML_T* cisml = mm_malloc(sizeof(CISML_T));
  cisml->program_name = strdup(program_name);
  cisml->command_line = strdup(command_line);
  cisml->pattern_file = strdup(pattern_file);
  cisml->sequence_file = strdup(sequence_file);
  cisml->background_file = NULL;
  cisml->pattern_pvalue_cutoff = NULL;
  cisml->sequence_pvalue_cutoff = NULL;
  cisml->site_pvalue_cutoff = NULL;
  cisml->site_qvalue_cutoff = NULL;
  cisml->num_passing_cutoff = 0;
  cisml->sequence_filter = NULL;
  cisml->num_multi_patterns = 0;
  cisml->num_allocated_multi_patterns = 0;
  cisml->multi_patterns = NULL;
  cisml->num_patterns = 0;
  cisml->num_allocated_patterns = 0;
  cisml->patterns = NULL;
  return cisml;
}

/**********************************************************************
  free_cisml

  Destructor for the cisml data structure.
**********************************************************************/
void free_cisml(CISML_T *cisml) {

  assert(cisml != NULL);

  while(cisml->num_multi_patterns > 0) {
    free_multi_pattern(cisml->multi_patterns[--cisml->num_multi_patterns]);
  }
  myfree(cisml->multi_patterns);
  cisml->num_allocated_multi_patterns = 0;

  while(cisml->num_patterns > 0) {
    free_pattern(cisml->patterns[--cisml->num_patterns]);
  }
  myfree(cisml->patterns);
  cisml->num_allocated_patterns = 0;

  myfree(cisml->sequence_filter);
  myfree(cisml->site_pvalue_cutoff);
  myfree(cisml->site_qvalue_cutoff);
  myfree(cisml->sequence_pvalue_cutoff);
  myfree(cisml->pattern_pvalue_cutoff);
  myfree(cisml->background_file);
  myfree(cisml->sequence_file);
  myfree(cisml->command_line);
  myfree(cisml->pattern_file);
  myfree(cisml->program_name);

  myfree(cisml);

}

/**********************************************************************
  get_cisml_program_name

  Gets the program_name member from a cisml object.
**********************************************************************/
char *get_cisml_program_name(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->program_name;
}

/**********************************************************************
  get_cisml_command_line

  Gets the command_line member in a cisml object.
**********************************************************************/
char *get_cisml_command_line(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->command_line;
}

/**********************************************************************
  get_cisml_pattern_file

  Gets the pattern_file member in a cisml object.
**********************************************************************/
char *get_cisml_pattern_file(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->pattern_file;
}

/**********************************************************************
  get_cisml_sequence_file

  Gets the sequence_file member from a cisml object.
**********************************************************************/
char *get_cisml_sequence_file(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->sequence_file;
}

/**********************************************************************
  set_cisml_background_file

  Sets the background_file member in a cisml object.
**********************************************************************/
void set_cisml_background_file(CISML_T *cisml, char *background_file) {
  assert(cisml != NULL);
  if (background_file == NULL) {
    if (cisml->background_file != NULL) {
      myfree(cisml->background_file);
    }
    cisml->background_file = NULL;
  }
  else {
    int new_length = strlen(background_file) + 1;
    int old_length = 0;
    if (cisml->background_file != NULL) {
      old_length = strlen(cisml->background_file) + 1;
    }
    if (old_length < new_length) {
      cisml->background_file = realloc(cisml->background_file, new_length);
    }
    strncpy(cisml->background_file, background_file, new_length);
  }
}

/**********************************************************************
  get_cisml_background_file

  Gets the background_file member from a cisml object.
  Return value may be NULL.
**********************************************************************/
char *get_cisml_background_file(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->background_file;
}

/**********************************************************************
  set_cisml_pattern_pvalue_cutoff

  Sets the pattern_pvalue_cutoff member in a cisml object.
**********************************************************************/
void set_cisml_pattern_pvalue_cutoff(
  CISML_T *cisml,
  double pattern_pvalue_cutoff
) {
  assert(cisml != NULL);
  if (cisml->pattern_pvalue_cutoff == NULL) {
    cisml->pattern_pvalue_cutoff = mm_malloc(sizeof(double));
  }
  *(cisml->pattern_pvalue_cutoff) = pattern_pvalue_cutoff;
}

/**********************************************************************
  clear_cisml_pattern_pvalue_cutoff

  Sets the pattern_pvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_pattern_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->pattern_pvalue_cutoff != NULL) {
    myfree(cisml->pattern_pvalue_cutoff);
  }
  cisml->pattern_pvalue_cutoff = NULL;
}

/**********************************************************************
  has_cisml_pattern_pvalue_cutoff

  Does a cisml object have a pattern_pvalue_cutoff?
**********************************************************************/
bool has_cisml_pattern_pvalue_cutoff(CISML_T *cisml) {
  return cisml->pattern_pvalue_cutoff != NULL ? true : false;

}

/**********************************************************************
  get_cisml_pattern_pvalue_cutoff

  Gets the pattern_pvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the pattern p-value cutoff has not been set.
**********************************************************************/
double get_cisml_pattern_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->pattern_pvalue_cutoff) {
    return *(cisml->pattern_pvalue_cutoff);
  }
  else {
    return 1.0;
  }
}

/**********************************************************************
  get_cisml_num_passing_cuttoff

  Returns the number of matches passing the significance cuttoff
**********************************************************************/
int get_cisml_num_passing_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->num_passing_cutoff;
}

/**********************************************************************
  set_cisml_sequence_pvalue_cutoff

  Sets the sequence_pvalue_cutoff member in a cisml object.
**********************************************************************/
void set_cisml_sequence_pvalue_cutoff(
  CISML_T *cisml,
  double sequence_pvalue_cutoff
) {
  assert(cisml != NULL);
  if (cisml->sequence_pvalue_cutoff == NULL) {
    cisml->sequence_pvalue_cutoff = mm_malloc(sizeof(double));
  }
  *(cisml->sequence_pvalue_cutoff) = sequence_pvalue_cutoff;
}

/**********************************************************************
  clear_cisml_sequence_pvalue_cutoff

  Sets the sequence_pvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_sequence_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->sequence_pvalue_cutoff != NULL) {
    myfree(cisml->sequence_pvalue_cutoff);
  }
  cisml->sequence_pvalue_cutoff = NULL;
}

/**********************************************************************
  has_cisml_sequence_pvalue_cutoff

  Does a cisml object have a sequence_pvalue_cutoff?
**********************************************************************/
bool has_cisml_sequence_pvalue_cutoff(CISML_T *cisml) {
  return cisml->sequence_pvalue_cutoff != NULL ? true : false;
}

/**********************************************************************
  get_cisml_sequence_pvalue_cutoff

  Gets the sequence_pvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the sequence p-value cutoff has not been set.
**********************************************************************/
double get_cisml_sequence_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->sequence_pvalue_cutoff) {
    return *(cisml->sequence_pvalue_cutoff);
  }
  else {
    return 1.0;
  }
}

/**********************************************************************
  set_cisml_site_pvalue_cutoff

  Sets the site_pvalue_cutoff member in a cisml object.
**********************************************************************/
void set_cisml_site_pvalue_cutoff(CISML_T *cisml, double site_pvalue_cutoff) {
  assert(cisml != NULL);
  if (cisml->site_pvalue_cutoff == NULL) {
    cisml->site_pvalue_cutoff = mm_malloc(sizeof(double));
  }
  *(cisml->site_pvalue_cutoff) = site_pvalue_cutoff;
}

/**********************************************************************
  clear_cisml_site_pvalue_cutoff

  Sets the site_pvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_site_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->site_pvalue_cutoff != NULL) {
    myfree(cisml->site_pvalue_cutoff);
  }
  cisml->site_pvalue_cutoff = NULL;
}

/**********************************************************************
  has_cisml_site_pvalue_cutoff

  Does a cisml object have a site_pvalue_cutoff?
**********************************************************************/
bool has_cisml_site_pvalue_cutoff(CISML_T *cisml) {
  return cisml->site_pvalue_cutoff != NULL ? true : false;
}

/**********************************************************************
  get_cisml_site_pvalue_cutoff

  Gets the site_pvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the site p-value cutoff has not been set.
**********************************************************************/
double get_cisml_site_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->site_pvalue_cutoff) {
    return *(cisml->site_pvalue_cutoff);
  }
  else {
    return 1.0;
  }
}

/**********************************************************************
  set_cisml_site_qvalue_cutoff

  Sets the site_qvalue_cutoff member in a cisml object.
**********************************************************************/
void set_cisml_site_qvalue_cutoff(CISML_T *cisml, double site_qvalue_cutoff) {
  assert(cisml != NULL);
  if (cisml->site_qvalue_cutoff == NULL) {
    cisml->site_qvalue_cutoff = mm_malloc(sizeof(double));
  }
  *(cisml->site_qvalue_cutoff) = site_qvalue_cutoff;
}

/**********************************************************************
  clear_cisml_site_qvalue_cutoff

  Sets the site_qvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_site_qvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->site_qvalue_cutoff != NULL) {
    myfree(cisml->site_qvalue_cutoff);
  }
  cisml->site_qvalue_cutoff = NULL;
}

/**********************************************************************
  has_cisml_site_qvalue_cutoff

  Does a cisml object have a site_qvalue_cutoff?
**********************************************************************/
bool has_cisml_site_qvalue_cutoff(CISML_T *cisml) {
  return cisml->site_qvalue_cutoff != NULL ? true : false;
}

/**********************************************************************
  get_cisml_site_qvalue_cutoff

  Gets the site_qvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the site q-value cutoff has not been set.
**********************************************************************/
double get_cisml_site_qvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->site_qvalue_cutoff) {
    return *(cisml->site_qvalue_cutoff);
  }
  else {
    return 1.0;
  }
}

/**********************************************************************
  get_cisml_sequence_filter

  Gets the sequence_filter member from a cisml object.
  May return NULL.
**********************************************************************/
char *get_cisml_sequence_filter(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->sequence_filter;
}

/**********************************************************************
  set_cisml_sequence_filter

  Sets the sequence_filter member in a cisml object.
**********************************************************************/
void set_cisml_sequence_filter(CISML_T *cisml, char *sequence_filter) {
  assert(cisml != NULL);
  if (sequence_filter == NULL) {
    if (cisml->sequence_filter != NULL) {
      myfree(cisml->sequence_filter);
    }
    cisml->sequence_filter = NULL;
  }
  else {
    int new_length = strlen(sequence_filter) + 1;
    int old_length = 0;
    if (cisml->sequence_filter != NULL) {
      old_length = strlen(cisml->sequence_filter) + 1;
    }
    if (old_length < new_length) {
      cisml->sequence_filter = realloc(cisml->sequence_filter, new_length);
    }
    strncpy(cisml->sequence_filter, sequence_filter, new_length);
  }
}

/**********************************************************************
  get_cisml_num_multi_patterns

  Gets the number of multi-patterns from a cisml object.
**********************************************************************/
int get_cisml_num_multi_patterns(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->num_multi_patterns;
}

/**********************************************************************
  get_cisml_multi_patterns

  Gets the array of pointers to multi_patterns from a cisml object.
  May return NULL.
**********************************************************************/
MULTI_PATTERN_T **get_cisml_multi_patterns(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->multi_patterns;
}

/**********************************************************************
  get_cisml_num_patterns

  Gets the number of patterns from a cisml object.
**********************************************************************/
int get_cisml_num_patterns(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->num_patterns;
}

/**********************************************************************
  get_cisml_patterns

  Gets the array of pointers to patterns from a cisml object.
  May return NULL.
**********************************************************************/
PATTERN_T **get_cisml_patterns(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->patterns;
}

/**********************************************************************
  get_cisml_num_stored_matches

  Gets the number of stored matches for all patterns in a cisml object.
**********************************************************************/
int get_cisml_num_stored_matches(CISML_T *cisml) {
  assert(cisml != NULL);
  int num_patterns = get_cisml_num_patterns(cisml);
  int total_stored_matches = 0;
  int i = 0;
  for(i = 0; i < num_patterns; ++i) {
    total_stored_matches += get_pattern_num_stored_matches(cisml->patterns[i]);
  }
  return total_stored_matches;
}

/**********************************************************************
  add_cisml_multi_pattern

  Adds a pattern to the array of pointers to patterns in a cisml object.
**********************************************************************/
void add_cisml_multi_pattern(CISML_T *cisml, MULTI_PATTERN_T* multi_pattern) {
  assert(cisml != NULL);
  if (multi_pattern == NULL) {
    return;
  }
  assert(cisml->num_multi_patterns <= cisml->num_allocated_multi_patterns);
  if (cisml->num_multi_patterns == cisml->num_allocated_multi_patterns) {
    cisml->num_allocated_multi_patterns += PATTERN_INCREMENT;
    cisml->multi_patterns = mm_realloc(
      cisml->multi_patterns,
      cisml->num_allocated_multi_patterns * sizeof(MULTI_PATTERN_T *)
    );
  }
  cisml->multi_patterns[cisml->num_multi_patterns] = multi_pattern;
  cisml->num_multi_patterns++;
}

/**********************************************************************
  add_cisml_pattern

  Adds a pattern to the array of pointers to patterns in a cisml object.
**********************************************************************/
void add_cisml_pattern(CISML_T *cisml, PATTERN_T* pattern) {
  assert(cisml != NULL);
  if (pattern == NULL) {
    return;
  }
  assert(cisml->num_patterns <= cisml->num_allocated_patterns);
  if (cisml->num_patterns == cisml->num_allocated_patterns) {
    cisml->num_allocated_patterns += PATTERN_INCREMENT;
    cisml->patterns = mm_realloc(
      cisml->patterns,
      cisml->num_allocated_patterns * sizeof(PATTERN_T *)
    );
  }
  cisml->patterns[cisml->num_patterns] = pattern;
  cisml->num_patterns++;
}

/**********************************************************************
  allocate_cisml_match_iterator

  Constructor for an iterator used to traverse all matched elements
  in in a cisml object. It stores the current index into the array of 
  matched elements for each of the patterns in the CisML object.
**********************************************************************/
CISML_MATCH_IT_T *allocate_cisml_match_iterator(CISML_T *cisml) {


  CISML_MATCH_IT_T *it = mm_malloc(sizeof(CISML_MATCH_IT_T));
  it->cisml = cisml;
  it->pattern_match_indices = calloc(cisml->num_patterns, sizeof(int));
  it->pattern_match_limits = calloc(cisml->num_patterns, sizeof(int));
  it->num_matches = 0;
  it->num_matches_returned = 0;

  int i;
  for (i = 0; i < cisml->num_patterns; ++i) {
    (it->pattern_match_limits)[i] = ((cisml->patterns)[i])->num_stored_matches;
    it->num_matches += (cisml->patterns[i])->num_stored_matches;
  }

  return it;

}

/**********************************************************************
  free_cisml_match_iterator

  Destructor for the CisML match iterator.
**********************************************************************/
void free_cisml_match_iterator(CISML_MATCH_IT_T *it) {
  myfree(it->pattern_match_indices);
  myfree(it->pattern_match_limits);
  myfree(it);
}

/**********************************************************************
  cisml_match_iterator_next

  Returns the next matched element from the CisML object.
  Matches are returned in order of ascending p-value.
***********************************************************************/
MATCHED_ELEMENT_T *cisml_match_iterator_next(CISML_MATCH_IT_T *it) {

  // Find the matched element in all the patterns with the minimum pvalue
  double min_pvalue = 1.1;

  if (it->num_matches_returned >= it->num_matches) {
    // No more matches in this cisml
    return NULL;
  }

  int pattern_index_of_min_pvalue = -1;
  MATCHED_ELEMENT_T *min_match = NULL;
  int i;
  for (i = 0; i < (it->cisml)->num_patterns; ++i) {
    if ((it->pattern_match_indices)[i] >= (it->pattern_match_limits)[i]) {
      // No matches left in this pattern
      continue;
    }
    int pattern_match_index = (it->pattern_match_indices)[i];
    MATCHED_ELEMENT_T *match = ((((it->cisml)->patterns[i])->elements)[pattern_match_index]);
    if (match->pvalue <= min_pvalue) {
      pattern_index_of_min_pvalue = i;
      min_match = match;
      min_pvalue = min_match->pvalue;
    }
  }
  ++(it->pattern_match_indices)[pattern_index_of_min_pvalue];
  ++(it->num_matches_returned);

  return min_match;
}

/**********************************************************************
  allocate_multi_pattern

  Constructor for the cisml multi_pattern data structure.
**********************************************************************/
MULTI_PATTERN_T *allocate_multi_pattern() {

  // Allocate memory and initialze fields
  MULTI_PATTERN_T *multi_pattern = mm_malloc(sizeof(MULTI_PATTERN_T));
  multi_pattern->match = NULL;
  multi_pattern->score = NULL;
  multi_pattern->pvalue = NULL;
  multi_pattern->num_patterns = 0;
  multi_pattern->num_allocated_patterns = 0;
  multi_pattern->patterns = NULL;
  return multi_pattern;
}

/**********************************************************************
  free_multi_pattern

  Destructor for the cisml multi_pattern data structure.
**********************************************************************/
void free_multi_pattern(MULTI_PATTERN_T *multi_pattern) {

  assert(multi_pattern != NULL);

  if (multi_pattern->match != NULL) {
    free_multi_pattern_match(multi_pattern->match);
  }

  while(multi_pattern->num_patterns > 0) {
    free_pattern(multi_pattern->patterns[--multi_pattern->num_patterns]);
  }
  myfree(multi_pattern->patterns);

  myfree(multi_pattern->score);
  myfree(multi_pattern->pvalue);

  myfree(multi_pattern);

}

/**********************************************************************
  set_multi_pattern_score

  Sets the score member in a cisml pattern object.
**********************************************************************/
void set_multi_pattern_score(MULTI_PATTERN_T *multi_pattern, double score) {
  assert(multi_pattern != NULL);
  if (multi_pattern->score == NULL) {
    multi_pattern->score = mm_malloc(sizeof(double));
  }
  *(multi_pattern->score) = score;
}

/**********************************************************************
  has_multi_pattern_score

  Does a multi_pattern object have a score?
**********************************************************************/
bool has_multi_pattern_score(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  return multi_pattern->score != NULL ? true : false;
}

/**********************************************************************
  clear_multi_pattern_score

  Sets the score member in a cisml multi_pattern object to null.
**********************************************************************/
void clear_multi_pattern_score(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  if (multi_pattern->score != NULL) {
    myfree(multi_pattern->score);
  }
  multi_pattern->score = NULL;
}

/**********************************************************************
  get_multi_pattern_score

  Gets the score member from a cisml multi_pattern object.
**********************************************************************/
double get_multi_pattern_score(MULTI_PATTERN_T* multi_pattern) {
  assert(multi_pattern != NULL);
  return *(multi_pattern->score);
}

/**********************************************************************
  set_multi_pattern_pvalue

  Sets the pvalue member in a cisml pattern object.
**********************************************************************/
void set_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern, double pvalue) {
  assert(multi_pattern != NULL);
  if (multi_pattern->pvalue == NULL) {
    multi_pattern->pvalue = mm_malloc(sizeof(double));
  }
  *(multi_pattern->pvalue) = pvalue;
}

/**********************************************************************
  has_multi_pattern_pvalue

  Does a multi_pattern object have a pvalue?
**********************************************************************/
bool has_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  return multi_pattern->pvalue != NULL ? true : false;
}

/**********************************************************************
  clear_multi_pattern_pvalue

  Sets the pvalue member in a cisml multi_pattern object to null.
**********************************************************************/
void clear_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  if (multi_pattern->pvalue != NULL) {
    myfree(multi_pattern->pvalue);
  }
  multi_pattern->pvalue = NULL;
}

/**********************************************************************
  get_multi_pattern_pvalue

  Gets the pvalue member from a cisml multi_pattern object.
**********************************************************************/
double get_multi_pattern_pvalue(MULTI_PATTERN_T* multi_pattern) {
  assert(multi_pattern != NULL);
  return *(multi_pattern->pvalue);
}

/**********************************************************************
  get_multi_pattern_num_patterns

  Gets the number of patterns from a cisml multi_pattern object.
**********************************************************************/
int get_multi_pattern_num_patterns(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  return multi_pattern->num_patterns;
}

/**********************************************************************
  get_multi_pattern_match

  Sets the match member in a cisml multi-pattern object.
**********************************************************************/
MULTI_PATTERN_MATCH_T *get_multi_pattern_match(
  MULTI_PATTERN_T *multi_pattern
) {
  assert(multi_pattern != NULL);
  return multi_pattern->match;
}

/**********************************************************************
  set_multi_pattern_match

  Sets the match member in a cisml multi-pattern object.
**********************************************************************/
void set_multi_pattern_match(
  MULTI_PATTERN_T *multi_pattern, 
  MULTI_PATTERN_MATCH_T *match
) {
  assert(multi_pattern != NULL);
  assert(match != NULL);
  multi_pattern->match = match;
}

/**********************************************************************
  get_multi_pattern_patterns

  Gets the array of pointers to patterns from a cisml object.
**********************************************************************/
PATTERN_T **get_multi_pattern_patterns(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  return multi_pattern->patterns;
}

/**********************************************************************
  add_multi_pattern_pattern

  Adds a pattern to the array of pointers to patterns in a multi_pattern
  object.
**********************************************************************/
void add_multi_pattern_pattern(
  MULTI_PATTERN_T *multi_pattern,
  PATTERN_T* pattern
) {
  assert(multi_pattern != NULL);
  if (pattern == NULL) {
    return;
  }
  assert(multi_pattern->num_patterns <= multi_pattern->num_allocated_patterns);
  if (multi_pattern->num_patterns == multi_pattern->num_allocated_patterns) {
    multi_pattern->num_allocated_patterns += PATTERN_INCREMENT;
    multi_pattern->patterns = mm_realloc(
      multi_pattern->patterns,
      multi_pattern->num_allocated_patterns * sizeof(PATTERN_T *)
    );
  }
  multi_pattern->patterns[multi_pattern->num_patterns] = pattern;
  multi_pattern->num_patterns++;
}

/**********************************************************************
  allocate_multi_pattern_match

  Constructor for the cisml multi_pattern_match data structure.
**********************************************************************/
MULTI_PATTERN_MATCH_T *allocate_multi_pattern_match(
    char *seq_name,
    char *seq,
    int start,
    int stop,
    double evalue 
) {

  assert(seq_name != NULL);
  assert(start >= 0);
  assert(stop >= 0);

  static int cluster_id = 0;

  // Allocate memory and initialze fields
  MULTI_PATTERN_MATCH_T *match = mm_malloc(sizeof(MULTI_PATTERN_MATCH_T));
  ++cluster_id;
  match->clusterid = mm_malloc(100);
  sprintf(match->clusterid, "cluster-%d", cluster_id);
  match->seq_name = strdup(seq_name);
  if (seq != NULL) {
    match->sequence = strdup(seq);
  }
  match->start = start;
  match->stop = stop;
  match->evalue = evalue;
  match->qvalue = NAN;

  return match;
}

/**********************************************************************
  free_multi_pattern_match

  Destructor for the cisml multi_pattern_match data structure.
**********************************************************************/
void free_multi_pattern_match(MULTI_PATTERN_MATCH_T *match) {

  assert(match != NULL);

  myfree(match->clusterid);
  myfree(match->seq_name);
  myfree(match->sequence);
  myfree(match);
}

/**********************************************************************
  set_multi_pattern_match_evalue

  Sets the evalue member in a cisml multi_pattern_match object.
**********************************************************************/
void set_multi_pattern_match_evalue(MULTI_PATTERN_MATCH_T *match, double evalue) {
  assert(match != NULL);
  match->evalue = evalue;
}

/**********************************************************************
  get_multi_pattern_match_pvalue

  Gets the evalue member from a cisml multi_pattern_match object.
**********************************************************************/
double get_multi_pattern_match_evalue(MULTI_PATTERN_MATCH_T* match) {
  assert(match != NULL);
  return match->evalue;
}

/**********************************************************************
  set_multi_pattern_match_qvalue

  Sets the qvalue member in a cisml multi_pattern_match object.
**********************************************************************/
void set_multi_pattern_match_qvalue(MULTI_PATTERN_MATCH_T *match, double qvalue) {
  assert(match != NULL);
  match->qvalue = qvalue;
}

/**********************************************************************
  get_multi_pattern_match_qvalue

  Gets the qvalue member from a cisml multi_pattern_match object.
**********************************************************************/
double get_multi_pattern_match_qvalue(MULTI_PATTERN_MATCH_T* match) {
  assert(match != NULL);
  return match->qvalue;
}

/**********************************************************************
  allocate_pattern

  Constructor for the cisml pattern data structure.
  Sets required fields to point to copies of the provided arguments.
  Other fields set to NULL.
**********************************************************************/
PATTERN_T *allocate_pattern(char *accession, char *name) {

  assert(accession != NULL);
  assert(name != NULL);

  // Allocate memory and initialze fields
  PATTERN_T *pattern = mm_malloc(sizeof(PATTERN_T));
  pattern->accession = NULL;
  pattern->name = NULL;
  pattern->pvalue = NULL;
  pattern->score = NULL;
  pattern->db = NULL;
  pattern->lsid = NULL;
  pattern->sequences = NULL;
  pattern->elements = NULL;

  pattern->num_allocated_sequences = 0;
  pattern->num_allocated_elements = 0;

  pattern->num_sequences = 0;
  pattern->num_scanned_positions = 0L;
  pattern->max_stored_matches = 100000;
  pattern->num_stored_matches = 0;
  pattern->max_pvalue_retained = 1.0;
  pattern->qvalues_computed = false;
  pattern->is_complete = false;

  // Set required fields
  pattern->accession = strdup(accession);
  pattern->name = strdup(name);
  pattern->element_heap = create_heap(
    pattern->max_stored_matches,
    compare_matched_elements,
    copy_matched_element,
    destroy_matched_element,
    NULL, // Key function
    NULL  // Print function
  );
  pattern->elements = NULL;

  return pattern;
}

/**********************************************************************
  free_pattern

  Destructor for the cisml pattern data structure.
**********************************************************************/
void free_pattern(PATTERN_T *pattern) {

  assert(pattern != NULL);

  myfree(pattern->lsid);
  myfree(pattern->db);
  myfree(pattern->score);
  myfree(pattern->pvalue);
  myfree(pattern->name);
  myfree(pattern->accession);

  while (pattern->num_sequences > 0) {
    free_scanned_sequence(pattern->sequences[--(pattern->num_sequences)]);
  }
  pattern->num_allocated_sequences = 0;
  pattern->num_sequences = 0;
  myfree(pattern->sequences)

  pattern->num_scanned_positions = 0L;
  destroy_heap(pattern->element_heap);
  // Free any elements stored in the elements array
  if (pattern->num_stored_matches > 0) {
    int i = 0;
    for (i = 0; i < pattern->num_stored_matches; i++) {
      free_matched_element(pattern->elements[i]);
    }
  }
  myfree(pattern->elements)
  myfree(pattern);

}

/**********************************************************************
  get_pattern_is_complete

  Returns a flag indicating whether or not all matched elements have
  been added to the pattern. If flag is true, the element heap is no
  longer available and all matched elements are stored in an array
  of matched element pointers, sorted by p-value.
**********************************************************************/
bool get_pattern_is_complete(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->is_complete;
}

/**********************************************************************
  set_pattern_is_complete

  Sets the flag indicating that all matched elements have
  been added to the pattern to true.  Moves  all matched elments out of
  the element heap into an array of matched element pointers, sorted by
  p-value. No further elements can be added to the pattern once
  this function has been called.
**********************************************************************/
void set_pattern_is_complete(PATTERN_T *pattern) {

  assert(pattern != NULL);
  assert(pattern->is_complete == false);

  pattern->is_complete = true;

  // Now that the pattern is complete, move elements from heap to array
  // sorted by p-value
  int num_elements = pattern->num_stored_matches;
  int i_element = 0;
  MATCHED_ELEMENT_T *element = NULL;
  pattern->elements = mm_malloc(sizeof(MATCHED_ELEMENT_T *) * num_elements);
  // Elements come off the heap in descending p-value order
  // Need to have array in ascending p-value order.
  for (i_element = num_elements - 1; i_element >= 0; --i_element) {
    element = (MATCHED_ELEMENT_T *) pop_heap_root(pattern->element_heap);
    pattern->elements[i_element] = element;
  }

  // Now that the pattern is complete, update the scanned sequeences
  // with the matched elements
  add_pattern_elements_to_scanned_seq(pattern);

}

/**********************************************************************
  get_pattern_max_pvalue_retained

  Returns the maximum p-value of the matched elements retained by
  the pattern.
**********************************************************************/
double get_pattern_max_pvalue_retained(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->max_pvalue_retained;
}

/**********************************************************************
  set_pattern_max_pvalue_retained

  Sets the maximum p-value of the matched elements retained by
  the pattern.
**********************************************************************/
void set_pattern_max_pvalue_retained(PATTERN_T *pattern, double max_pvalue) {
  assert(pattern != NULL);
  pattern->max_pvalue_retained = max_pvalue;
}

/**********************************************************************
  get_pattern_name

  Gets the pattern name member from a cisml pattern object.
  Caller should not free this string.
**********************************************************************/
char *get_pattern_name(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->name;
}

/**********************************************************************
  get_pattern_accession

  Gets the accession member from a cisml pattern object.
**********************************************************************/
char *get_pattern_accession(PATTERN_T* pattern) {
  assert(pattern != NULL);
  return pattern->accession;
}

/**********************************************************************
  set_pattern_pvalue

  Sets the pvalue member in a cisml pattern object.
**********************************************************************/
void set_pattern_pvalue(PATTERN_T *pattern, double pvalue) {
  assert(pattern != NULL);
  if (pattern->pvalue == NULL) {
    pattern->pvalue = mm_malloc(sizeof(double));
  }
  *(pattern->pvalue) = pvalue;
}

/**********************************************************************
  clear_pattern_pvalue

  Sets the pvalue member in a cisml pattern object to null.
**********************************************************************/
void clear_pattern_pvalue(PATTERN_T *pattern) {
  assert(pattern != NULL);
  if (pattern->pvalue != NULL) {
    myfree(pattern->pvalue);
  }
  pattern->pvalue = NULL;
}

/**********************************************************************
  has_pattern_pvalue

  Does a pattern object have a pvalue?
**********************************************************************/
bool has_pattern_pvalue(PATTERN_T *pattern) {
  return pattern->pvalue != NULL ? true : false;
}

/**********************************************************************
  get_pattern_pvalue

  Gets the pvalue member from a cisml pattern object.
**********************************************************************/
double get_pattern_pvalue(PATTERN_T* pattern) {
  assert(pattern != NULL);
  return *(pattern->pvalue);
}

/**********************************************************************
  has_pattern_qvalues

  Does the matched-elements for the pattern have qvalues available?
**********************************************************************/
bool has_pattern_qvalues(PATTERN_T *pattern) {
  return pattern->qvalues_computed;
}

/**********************************************************************
  set_pattern_score

  Sets the score member in a cisml pattern object.
**********************************************************************/
void set_pattern_score(PATTERN_T *pattern, double score) {
  assert(pattern != NULL);
  if (pattern->score == NULL) {
    pattern->score = mm_malloc(sizeof(double));
  }
  *(pattern->score) = score;
}

/**********************************************************************
  clear_pattern_score

  Sets the score member in a cisml pattern object to null.
**********************************************************************/
void clear_pattern_score(PATTERN_T *pattern) {
  assert(pattern != NULL);
  if (pattern->score != NULL) {
    myfree(pattern->score);
  }
  pattern->score = NULL;
}

/**********************************************************************
  has_pattern_score

  Does a pattern object have a score?
**********************************************************************/
bool has_pattern_score(PATTERN_T *pattern) {
  return pattern->score != NULL ? true : false;
}

/**********************************************************************
  get_pattern_score

  Gets the score member from a cisml pattern object.
**********************************************************************/
double get_pattern_score(PATTERN_T* pattern) {
  assert(pattern != NULL);
  return *(pattern->score);
}

/**********************************************************************
  set_pattern_db

  Sets the db member in a cisml pattern object.
**********************************************************************/
void set_pattern_db(PATTERN_T *pattern, char *db) {
  assert(pattern != NULL);
  if (db == NULL) {
    if (pattern->db != NULL) {
      myfree(pattern->db);
    }
    pattern->db = NULL;
  }
  else {
    int old_length = 0;
    int new_length = strlen(db) + 1;
    if (pattern->db != NULL) {
      old_length = strlen(pattern->db) + 1;
    }
    if (old_length < new_length) {
      pattern->db = mm_realloc(pattern->db, new_length);
    }
    strncpy(pattern->db, db, new_length);
  }
}

/**********************************************************************
  get_pattern_db

  Gets the db member from a cisml pattern object.
  May return NULL.
**********************************************************************/
char *get_pattern_db(PATTERN_T* pattern) {
  assert(pattern != NULL);
  return pattern->db;
}

/**********************************************************************
  set_pattern_lsid

  Sets the lsid member in a cisml pattern object.
**********************************************************************/
void set_pattern_lsid(PATTERN_T *pattern, char *lsid) {
  assert(pattern != NULL);
  if (lsid == NULL) {
    if (pattern->lsid != NULL) {
      myfree(pattern->lsid);
    }
    pattern->lsid = NULL;
  }
  else {
    int old_length = 0;
    int new_length = strlen(lsid) + 1;
    if (pattern->lsid != NULL) {
      old_length = strlen(pattern->lsid) + 1;
    }
    if (old_length < new_length) {
      pattern->lsid = mm_realloc(pattern->lsid, new_length);
    }
    strncpy(pattern->lsid, lsid, new_length);
  }
}

/**********************************************************************
  get_pattern_lsid

  Gets the lsid member from a cisml pattern object.
  May return NULL.
**********************************************************************/
char *get_pattern_lsid(PATTERN_T* pattern) {
  assert(pattern != NULL);
  return pattern->lsid;
}

/**********************************************************************
  get_pattern_num_scanned_sequences

  Gets the number of scanned_sequence objects in a cisml pattern object.
**********************************************************************/
int get_pattern_num_scanned_sequences(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->num_sequences;
}

/**********************************************************************
  get_pattern_max_stored_matches

  Gets the maximum number of matched elements that will be stored in a cisml
  pattern object.
**********************************************************************/
int get_pattern_max_stored_matches(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->max_stored_matches;
}

/**********************************************************************
  set_pattern_max_stored_matches

  Sets the maximum number of matched elements that will be stored in a cisml
  pattern object. This requires creating a new element heap using the new
  maximum number of elements, copying the elements from the existing heap,
  freeing the existing heap, and pointing the pattern to the new heap.

  It fails if the new max is smaller than the current number of nodes.

  Returns true if successful, false otherwise.
**********************************************************************/
bool set_pattern_max_stored_matches(PATTERN_T *pattern, int max) {

  assert(pattern != NULL);

  HEAP *current_heap = pattern->element_heap;
  if (pattern->max_stored_matches == max) {
    // No point in doing anything, the max isn't changing.
    return true;
  }
  else if (max > get_num_nodes(current_heap)) {
    HEAP *new_heap = create_heap(
      max,
      compare_matched_elements,
      copy_matched_element,
      destroy_matched_element,
      NULL, // Key function
      NULL  // Print function
    );
    MATCHED_ELEMENT_T *element = NULL;
    while ((element = pop_heap_root(current_heap))) {
      add_node_heap(new_heap, element);
    }
    destroy_heap(current_heap);
    pattern->element_heap = new_heap;
    pattern->max_stored_matches = max;
    return true;
  }
  else {
    // Can't make max size of heap smaller than the current number of nodes.
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(
        stderr,
        "Warning: The maximum size of the heap cannot be decreased.\n"
      );
    }
    return false;
  }
}

/**********************************************************************
  get_pattern_num_scanned_positions

  Gets the number of sites scanned with a cisml pattern object.
**********************************************************************/
long get_pattern_num_scanned_positions(PATTERN_T *pattern) {
  assert(pattern != NULL);
  assert(pattern->num_scanned_positions >= 0L);
  return pattern->num_scanned_positions;
}

/**********************************************************************
  get_pattern_num_stored_matches

  Gets the total number of matched element objects stored in a cisml
  pattern object.
**********************************************************************/
int get_pattern_num_stored_matches(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->num_stored_matches;
}

/**********************************************************************
  get_pattern_stored_matches

  Gets the array of matched element objects stored in a cisml
  pattern object.
**********************************************************************/
MATCHED_ELEMENT_T **get_pattern_stored_matches(PATTERN_T *pattern) {
  if (pattern->is_complete != true) {
    die("Attempt to retreive list of elements before pattern is complete.");
  }
  return pattern->elements;
}

/**********************************************************************
  add_pattern_matched_element

  Adds a pointer to a matched element to the heap of pointers to
  matched_elements in a cisml pattern object.

  Before adding to the heap, check the pvalue for the element.
  if it's greater the pvalues we've already thrown away, we refuse to
  add it. Also, check max_stored_matches. If we've
  hit the limit, purge the heap of the least significian matched elements.
  Set to false the Boolean that indicates that set of matched elements
  is complete.

  Returns true if the element was added, false otherwise
**********************************************************************/
bool add_pattern_matched_element(
  PATTERN_T *pattern,
  MATCHED_ELEMENT_T *element
) {

  assert(element != NULL);
  assert(pattern != NULL);

  if (pattern->is_complete == true) {
    // Don't add matched elements if pattern is marked as complete.
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(
        stderr,
        "Warning: trying to add matched elements to pattern marked as complete.\n"
      );
    }
    return false;
  }

  if (element->pvalue > pattern->max_pvalue_retained) {
    // Don't add the matched element if its pvalue is less
    // significant then pvalues we've already thrown away.
    return false;
  }

  if (pattern->max_stored_matches > 0
      && (long) pattern->max_stored_matches <= pattern->num_stored_matches) {
    // Max element storage has been reached. We have to drop
    // the least significant elements.
    reduce_pattern_matched_elements(pattern);
  }

  // The least significant pvalue may have changed
  // during the reduction.
  if (element->pvalue <= pattern->max_pvalue_retained) {
    // Add the element to the pattern.
    add_node_heap(pattern->element_heap, (void *) element);
    ++pattern->num_stored_matches;
    return true;
  }
  else {
    // Don't bother adding elements who's pvalue execeeds the
    // least sig. matched element retained.
    return false;
  }
}

/**********************************************************************
  add_pattern_matched_element_no_heap

  Adds a pointer to a matched element the array of
  matched_elements in a cisml pattern object.
  This is used when reading matches from an existing CisML file
  and there is no need to filter matches through the heap.

  Returns true if the element was added, false otherwise
**********************************************************************/
void add_pattern_matched_element_no_heap(
  PATTERN_T *pattern,
  MATCHED_ELEMENT_T *element
) {

  assert(element != NULL);
  assert(pattern != NULL);
  assert(pattern->num_stored_matches <= pattern->num_allocated_elements);

  if (pattern->num_stored_matches == pattern->num_allocated_elements) {
    pattern->num_allocated_elements += ELEMENT_INCREMENT;
    pattern->elements = mm_realloc(
      pattern->elements,
      pattern->num_allocated_elements * sizeof(MATCHED_ELEMENT_T *)
    );
  }
  pattern->elements[pattern->num_stored_matches] = element;
  pattern->num_stored_matches++;
}

/**********************************************************************
  add_pattern_elements_to_scanned_seq

  Updates each scanned sequence belonging to pattern with the
  matched_elements assocaited with that sequence.
  Should  not be called until pattern is complete.

**********************************************************************/
void add_pattern_elements_to_scanned_seq(PATTERN_T *pattern) {

  assert(pattern != NULL);
  assert(pattern->is_complete == true);
  assert(pattern->elements != NULL);

  int element_index = 0;
  for (element_index = 0; element_index < pattern->num_stored_matches; ++element_index) {
    MATCHED_ELEMENT_T *element = (pattern->elements)[element_index];
    SCANNED_SEQUENCE_T *seq = element->parent_sequence;
    add_scanned_sequence_matched_element(seq, element);
  }

}

/**********************************************************************
  reduce_pattern_matched_elements

  To conserve memory, remove the least significant matched_elements.
  We want to remove at least PERCENT_ELEMENT_DISCARD * num_stored_matches
  but we may remove more if there remain matched elments with a
  p-value equal to the delected elements.

**********************************************************************/
static void reduce_pattern_matched_elements(PATTERN_T *pattern) {

  assert(pattern != NULL);
  assert(pattern->is_complete == false);

  const float PERCENT_ELEMENT_DISCARD = 0.5;
  static bool have_moved_sequences = false;

  HEAP *heap = pattern->element_heap;

  // Delete PERCENT_ELEMENT_DISCARD with highest pvalue
  int num_elements_to_delete = PERCENT_ELEMENT_DISCARD * get_num_nodes(heap);
  if (verbosity > NORMAL_VERBOSE) {
          fprintf(
      stderr,
      "Deleting at least %d matched elements from pattern %s.\n",
      num_elements_to_delete,
      pattern->name
    );
  }
  // Delete least significant matched elements.
  double min_pvalue_discarded = 1.0;
  MATCHED_ELEMENT_T *victim = NULL;
  int deletion_count = 0;
  for (deletion_count = 0; deletion_count  < num_elements_to_delete; ++deletion_count) {
    victim = (MATCHED_ELEMENT_T *) pop_heap_root(heap);
    min_pvalue_discarded = victim->pvalue;
    --pattern->num_stored_matches;
    free_matched_element(victim);
  }

  // Keep deleting matched elements until we find an element more 
  // significant then the elements we've already deleted.
  while (((MATCHED_ELEMENT_T *) get_node(heap, 1))->pvalue 
         >= min_pvalue_discarded) {
    victim = (MATCHED_ELEMENT_T *) pop_heap_root(heap);
    assert(victim != NULL);
    --pattern->num_stored_matches;
    free_matched_element(victim);
    if (get_num_nodes(heap) == 0) {
      // All the matched elements have been deleted!
      break;
    }
  }

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(
      stderr, 
      "Warning: Reached max stored scores (%d).\n" 
      "Motif matches with p-value >= %3.2g have been "
      "dropped to reclaim memory.\n", 
      pattern->max_stored_matches,
      min_pvalue_discarded
    );
  }

  if (get_num_nodes(heap) > 0) {
    // Get the largest p-value retained from the top element of the heap.
    pattern->max_pvalue_retained = ((MATCHED_ELEMENT_T *) get_node(heap, 1))->pvalue;
  }
  else {
    // All items have been deleted!
    fprintf(
      stderr, 
      "Warning: there are no motif matches with p-value < %3.2g.\n"
      "Use --max-stored-scores to allocate more space for "
      "storing motif matches.\n", 
     min_pvalue_discarded 
    );
    // Set the largest p-value retained to something
    // slightly less the smallest p-value discarded.
    pattern->max_pvalue_retained 
      = get_next_smaller_double(min_pvalue_discarded);
  }

}

/**********************************************************************
  add_pattern_scanned_sequence

  Adds a pointer to a scanned_sequence to the array of pointers to
  scanned_sequences in a cisml pattern object.
**********************************************************************/
static void add_pattern_scanned_sequence(
  PATTERN_T *pattern,
  SCANNED_SEQUENCE_T *sequence
) {

  assert(pattern != NULL);
  assert(sequence != NULL);
  assert(pattern->num_sequences <= pattern->num_allocated_sequences);

  sequence->parent_pattern = pattern;

  if (pattern->num_sequences == pattern->num_allocated_sequences) {
    pattern->num_allocated_sequences += SEQUENCE_INCREMENT;
    pattern->sequences = mm_realloc(
      pattern->sequences,
      pattern->num_allocated_sequences * sizeof(SCANNED_SEQUENCE_T *)
    );
  }
  pattern->sequences[pattern->num_sequences] = sequence;
  pattern->num_sequences++;

}

/**********************************************************************
  get_pattern_scanned_sequences

  Gets the array of pointers to scanned_sequence objects in a cisml
  pattern object.  May return NULL.
**********************************************************************/
SCANNED_SEQUENCE_T **get_pattern_scanned_sequences(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->sequences;
}

/**********************************************************************
  allocate_scanned_sequence

  Constructor for the cisml scanned_sequence data structure.
  Sets required fields to point to copies of the provided arguments.
  Other fields set to NULL.
**********************************************************************/
SCANNED_SEQUENCE_T *allocate_scanned_sequence(
  char *accession,
  char *name,
  PATTERN_T *parent_pattern
) {

  assert(accession != NULL);
  assert(name != NULL);

  // Allocate memory and initialze fields
  SCANNED_SEQUENCE_T *scanned_sequence = mm_malloc(sizeof(SCANNED_SEQUENCE_T));
  scanned_sequence->accession = NULL;
  scanned_sequence->name = NULL;
  scanned_sequence->pvalue = NULL;
  scanned_sequence->score = NULL;
  scanned_sequence->length = NULL;
  scanned_sequence->db = NULL;
  scanned_sequence->lsid = NULL;
  scanned_sequence->num_scanned_positions = 0L;
  scanned_sequence->num_matched_elements = 0;
  scanned_sequence->num_allocated_elements = 0;

  // Set required fields
  scanned_sequence->accession = strdup(accession);
  scanned_sequence->name = strdup(name);
  add_pattern_scanned_sequence(parent_pattern, scanned_sequence);
  scanned_sequence->elements = NULL;

  return scanned_sequence;
}

/**********************************************************************
  free_scanned_sequence

  Destructor for the cisml scanned_sequence data structure.
**********************************************************************/
void free_scanned_sequence(SCANNED_SEQUENCE_T *scanned_sequence) {

  assert(scanned_sequence != NULL);

  myfree(scanned_sequence->lsid);
  myfree(scanned_sequence->db);
  myfree(scanned_sequence->length);
  myfree(scanned_sequence->score);
  myfree(scanned_sequence->pvalue);
  myfree(scanned_sequence->name);
  myfree(scanned_sequence->accession);
  myfree(scanned_sequence->elements); // Don't free individual elements
                                      // they belong to pattern.

  myfree(scanned_sequence);

}

/**********************************************************************
  get_scanned_sequence_accession

  Gets the accession member from a cisml scanned_sequence object.
**********************************************************************/
char *get_scanned_sequence_accession(SCANNED_SEQUENCE_T* scanned_sequence) {

  assert(scanned_sequence != NULL);

  return scanned_sequence->accession;

}

/**********************************************************************
  get_scanned_sequence_name

  Gets the program_name member from a cisml scanned_sequence object.
**********************************************************************/
char *get_scanned_sequence_name(SCANNED_SEQUENCE_T *scanned_sequence) {

  assert(scanned_sequence != NULL);

  return scanned_sequence->name;

}

/**********************************************************************
  get_scanned_sequence_parent

  Get a pointer to the pattern pointing to this scanned_sequence object.
**********************************************************************/
PATTERN_T *get_scanned_sequence_parent(SCANNED_SEQUENCE_T *scanned_sequence) {

  assert(scanned_sequence != NULL);

  return scanned_sequence->parent_pattern;

}

/**********************************************************************
  set_scanned_sequence_pvalue

  Sets the pvalue member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_pvalue(
  SCANNED_SEQUENCE_T *scanned_sequence,
  double pvalue
) {

  assert(scanned_sequence != NULL);

  if (scanned_sequence->pvalue == NULL) {
    scanned_sequence->pvalue = mm_malloc(sizeof(double));
  }
  *(scanned_sequence->pvalue) = pvalue;

}

/**********************************************************************
  clear_scanned_sequence_pvalue

  Sets the pvalue member in a cisml scanned_sequence object to null.
**********************************************************************/
void clear_scanned_sequence_pvalue(SCANNED_SEQUENCE_T *scanned_sequence) {

  assert(scanned_sequence != NULL);

  if (scanned_sequence->pvalue != NULL) {
    myfree(scanned_sequence->pvalue);
  }

  scanned_sequence->pvalue = NULL;

}

/**********************************************************************
  has_scanned_sequence_pvalue

  Does a scanned_sequence object have a pvalue?
**********************************************************************/
bool has_scanned_sequence_pvalue(SCANNED_SEQUENCE_T *scanned_sequence) {
  assert(scanned_sequence != NULL);
  return scanned_sequence->pvalue != NULL ? true : false;

}

/**********************************************************************
  get_scanned_sequence_pvalue

  Gets the pvalue member from a cisml scanned_sequence object.
**********************************************************************/
double get_scanned_sequence_pvalue(SCANNED_SEQUENCE_T* scanned_sequence) {
  assert(scanned_sequence != NULL);
  return *(scanned_sequence->pvalue);

}

/**********************************************************************
  set_scanned_sequence_score

  Sets the score member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_score(
  SCANNED_SEQUENCE_T *scanned_sequence,
  double score
) {
  assert(scanned_sequence != NULL);
  if (scanned_sequence->score == NULL) {
    scanned_sequence->score = mm_malloc(sizeof(double));
  }
  *(scanned_sequence->score) = score;
}

/**********************************************************************
  clear_scanned_sequence_score

  Sets the score member in a cisml scanned_sequence object to null.
**********************************************************************/
void clear_scanned_sequence_score(SCANNED_SEQUENCE_T *scanned_sequence) {
  assert(scanned_sequence != NULL);
  if (scanned_sequence->score != NULL) {
    myfree(scanned_sequence->score);
  }
  scanned_sequence->score = NULL;
}

/**********************************************************************
  has_scanned_sequence_score

  Does a scanned_sequence object have a score?
**********************************************************************/
bool has_scanned_sequence_score(SCANNED_SEQUENCE_T *scanned_sequence) {
  assert(scanned_sequence != NULL);
  return scanned_sequence->score != NULL ? true : false;
}

/**********************************************************************
  get_scanned_sequence_score

  Gets the score member from a cisml scanned_sequence object.
**********************************************************************/
double get_scanned_sequence_score(SCANNED_SEQUENCE_T* scanned_sequence) {
  assert(scanned_sequence != NULL);
  return *(scanned_sequence->score);
}

/**********************************************************************
  set_scanned_sequence_length

  Sets the length member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_length(
  SCANNED_SEQUENCE_T *scanned_sequence,
  int length
) {
  assert(scanned_sequence != NULL);
  if (scanned_sequence->length == NULL) {
    scanned_sequence->length = mm_malloc(sizeof(int));
  }
  *(scanned_sequence->length) = length;
}

/**********************************************************************
  clear_scanned_sequence_length

  Sets the length member in a cisml scanned_sequence object to null.
**********************************************************************/
void clear_scanned_sequence_length(SCANNED_SEQUENCE_T *scanned_sequence) {
  assert(scanned_sequence != NULL);
  if (scanned_sequence->length != NULL) {
    myfree(scanned_sequence->length);
  }
  scanned_sequence->length = NULL;
}

/**********************************************************************
  has_scanned_sequence_length

  Does a scanned_sequence object have a length?
**********************************************************************/
bool has_scanned_sequence_length(SCANNED_SEQUENCE_T *scanned_sequence) {
  return scanned_sequence->length != NULL ? true : false;
}

/**********************************************************************
  get_scanned_sequence_length

  Gets the length member from a cisml scanned_sequence object.
**********************************************************************/
int get_scanned_sequence_length(SCANNED_SEQUENCE_T* scanned_sequence) {
  assert(scanned_sequence != NULL);
  return *(scanned_sequence->length);
}

/**********************************************************************
  set_scanned_sequence_db

  Sets the db member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_db(SCANNED_SEQUENCE_T *scanned_sequence, char *db) {
  assert(scanned_sequence != NULL);
  if (db == NULL) {
    if (scanned_sequence->db != NULL) {
      myfree(scanned_sequence->db);
    }
    scanned_sequence->db = NULL;
  }
  else {
    int old_length = 0;
    int new_length = strlen(db) + 1;
    if (scanned_sequence->db != NULL) {
      old_length = strlen(scanned_sequence->db) + 1;
    }
    if (old_length < new_length) {
      scanned_sequence->db = mm_realloc(scanned_sequence->db, new_length);
    }
    strncpy(scanned_sequence->db, db, new_length);
  }
}

/**********************************************************************
  get_scanned_sequence_db

  Gets the db member from a cisml scanned_sequence object.
  May return NULL.
**********************************************************************/
char *get_scanned_sequence_db(SCANNED_SEQUENCE_T* scanned_sequence) {
  assert(scanned_sequence != NULL);
  return scanned_sequence->db;
}

/**********************************************************************
  set_scanned_sequence_lsid

  Sets the lsid member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_lsid(
  SCANNED_SEQUENCE_T *scanned_sequence,
  char *lsid
) {
  assert(scanned_sequence != NULL);
  if (lsid == NULL) {
    if (scanned_sequence->lsid != NULL) {
      myfree(scanned_sequence->lsid);
    }
    scanned_sequence->lsid = NULL;
  }
  else {
    int old_length = 0;
    int new_length = strlen(lsid) + 1;
    if (scanned_sequence->lsid != NULL) {
      old_length = strlen(scanned_sequence->lsid) + 1;
    }
    if (old_length < new_length) {
      scanned_sequence->lsid = mm_realloc(scanned_sequence->lsid, new_length);
    }
    strncpy(scanned_sequence->lsid, lsid, new_length);
  }
}

/**********************************************************************
  get_scanned_sequence_lsid

  Gets the lsid member from a cisml scanned_sequence object.
  May return NULL.
**********************************************************************/
char *get_scanned_sequence_lsid(SCANNED_SEQUENCE_T* scanned_sequence) {
  assert(scanned_sequence != NULL);
  return scanned_sequence->lsid;
}

/**********************************************************************
  add_scanned_sequence_matched_element

  Adds a matched element to the array of matched elements in
  the scanned sequence.
**********************************************************************/
void add_scanned_sequence_matched_element(
  SCANNED_SEQUENCE_T *sequence,
  MATCHED_ELEMENT_T *element
) {

  assert(sequence != NULL);
  assert(element != NULL);
  assert(sequence->num_matched_elements <= sequence->num_allocated_elements);

  if (sequence->num_matched_elements == sequence->num_allocated_elements) {
    if (sequence->num_allocated_elements == 0) {
      sequence->num_allocated_elements = 1;
    } else if (sequence->num_allocated_elements < ELEMENT_INCREMENT) {
      sequence->num_allocated_elements *= 2;
    } else {
      sequence->num_allocated_elements += ELEMENT_INCREMENT;
    }
    sequence->elements = mm_realloc(
      sequence->elements,
      sequence->num_allocated_elements * sizeof(MATCHED_ELEMENT_T *)
    );
  }
  sequence->elements[sequence->num_matched_elements] = element;
  sequence->num_matched_elements++;
}

/**********************************************************************
  compact_scanned_sequence_matched_elements

  Shrinks the array of matched elements to fit.
**********************************************************************/
void compact_scanned_sequence_matched_elements(
    SCANNED_SEQUENCE_T *sequence
) {
  if (sequence->num_allocated_elements > sequence->num_matched_elements) {
    sequence->elements = (MATCHED_ELEMENT_T**)mm_realloc(sequence->elements,
        sizeof(MATCHED_ELEMENT_T*) * sequence->num_matched_elements);
    sequence->num_allocated_elements = sequence->num_matched_elements;
  }
}

/**********************************************************************
  add_scanned_sequence_scanned_position

  Increments the count of scanned positions in a scanned_sequence.
**********************************************************************/
void add_scanned_sequence_scanned_position(SCANNED_SEQUENCE_T *sequence) {
  assert(sequence != NULL);
  sequence->num_scanned_positions++;
  sequence->parent_pattern->num_scanned_positions++;
}

/**********************************************************************
  get_scanned_sequence_num_matched_elements

  Gets the number of matched_element objects in a cisml
  scanned_sequence object.
**********************************************************************/
int get_scanned_sequence_num_matched_elements(SCANNED_SEQUENCE_T *sequence) {
  assert(sequence != NULL);
  return sequence->num_matched_elements;
}

/**********************************************************************
  get_scanned_sequence_num_scanned_positions

  Gets the number of positions in the scanned_sequence where we
  have scanned for a matched_element.
**********************************************************************/
long get_scanned_sequence_num_scanned_positions(SCANNED_SEQUENCE_T *sequence) {
  assert(sequence != NULL);
  return sequence->num_scanned_positions;
}

/**********************************************************************
  get_scanned_sequences_matched_elements

  Gets the array of pointers to matched_element objects in a cisml
  scanned_sequence object.  May return NULL.
**********************************************************************/
MATCHED_ELEMENT_T **get_scanned_sequence_matched_elements(
    SCANNED_SEQUENCE_T *sequence
) {
  assert(sequence != NULL);
  return sequence->elements;
}

/**********************************************************************
  allocate_matched_element

  Constructor for the cisml matched_element data structure.
  Sets required fields to the provided arguments.
  Other fields set to NULL.
**********************************************************************/
MATCHED_ELEMENT_T *allocate_matched_element(
  int start,
  int stop,
  SCANNED_SEQUENCE_T *parent_sequence
) {

  // Allocate memory and set required fields
  MATCHED_ELEMENT_T *element = mm_malloc(sizeof(MATCHED_ELEMENT_T));
  element->start = start;
  element->stop = stop;
  element->parent_sequence = parent_sequence;

  // Initialze optional fields
  element->score = 0.0;
  element->has_score = false;
  element->pvalue = 0.0;
  element->has_pvalue = false;
  element->qvalue = 0.0;
  element->has_qvalue = false;
  element->clusterid = NULL;
  element->sequence = NULL;
  element->strand = (start < stop) ? '+' : '-';

  return element;
}

/**********************************************************************
  allocate_matched_element_without_inversion

  Alternative Constructor for the cisml matched_element data structure.
  Sets required fields to the provided arguments.
  Other fields set to NULL.

  JH: I had to add this because the existing constructor inverts the
      DNA sequence if the start site is greater than the stop site. I
      could not simply copy an existing matched element's content and
      feed it to the constructor without the content being changed.

**********************************************************************/
MATCHED_ELEMENT_T *allocate_matched_element_without_inversion(
  int start,
  int stop,
  const char *seq,
  SCANNED_SEQUENCE_T *parent
) {

  // Allocate memory and set required fields
  MATCHED_ELEMENT_T *element = mm_malloc(sizeof(MATCHED_ELEMENT_T));
  element->start = start;
  element->stop = stop;
  int length = strlen(seq) + 1;
  element->sequence = mm_malloc(length * sizeof(char));
  strncpy(element->sequence, seq, length);

  element->parent_sequence = parent;

  // Initialze optional fields
  element->score = 0.0;
  element->has_score = false;
  element->pvalue = 0.0;
  element->has_pvalue = false;
  element->qvalue = 0.0;
  element->has_qvalue = false;
  element->clusterid = NULL;
  element->strand = '\0';

  return element;
}


/**********************************************************************
  allocate_matched_element_with_score

  Constructor for the cisml matched_element data structure.
  Sets required fields to the provided arguments.
  Sets score and pvalue.
  Other fields set to NULL.

**********************************************************************/
MATCHED_ELEMENT_T *allocate_matched_element_with_score(
  int start,
  int stop,
  double score,
  double pvalue,
  SCANNED_SEQUENCE_T *parent_sequence
) {

  // Allocate memory and set required fields
  MATCHED_ELEMENT_T *element = allocate_matched_element(start, stop, parent_sequence);

  set_matched_element_score(element, score);
  set_matched_element_pvalue(element, pvalue);

  return element;
}

/**********************************************************************
  free_matched_element

  Destructor for the cisml matched_element data structure.
**********************************************************************/
void free_matched_element(MATCHED_ELEMENT_T *element) {

  assert(element != NULL);

  myfree(element->clusterid);
  myfree(element->sequence);

  myfree(element);

}

/**********************************************************************
  get_matched_element_start

  Gets the start member from a cisml matched_element object.
**********************************************************************/
int get_matched_element_start(MATCHED_ELEMENT_T* element) {

  assert(element != NULL);

  return element->start;

}

/**********************************************************************
  set_matched_element_start

  Sets the start member from a cisml matched_element object.
**********************************************************************/
void set_matched_element_start(MATCHED_ELEMENT_T* matched_element, int newstart) {
  assert(matched_element != NULL);
  matched_element->start = newstart;
}

/**********************************************************************
  get_matched_element_stop

  Gets the stop member from a cisml matched_element object.
**********************************************************************/
int get_matched_element_stop(MATCHED_ELEMENT_T* element) {
  assert(element != NULL);
  return element->stop;
}

/**********************************************************************
  set_matched_element_stop

  sets the stop member from a cisml matched_element object.
**********************************************************************/
void set_matched_element_stop(MATCHED_ELEMENT_T* element, int newstop) {
  assert(element != NULL);
  element->stop = newstop;
}

/**********************************************************************
  set_matched_element_score

  Sets the score member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_score(
  MATCHED_ELEMENT_T *element,
  double score
) {
  assert(element != NULL);
  element->score = score;
  element->has_score = true;
}

/**********************************************************************
  has_matched_element_score

  Does a matched_element object have a score?
**********************************************************************/
bool has_matched_element_score(MATCHED_ELEMENT_T *element) {
  return element->has_score;
}

/**********************************************************************
  get_matched_element_score

  Gets the score member from a cisml matched_element object.
**********************************************************************/
double get_matched_element_score(MATCHED_ELEMENT_T* element) {
  assert(element != NULL);
  return element->score;
}

/**********************************************************************
  set_matched_element_pvalue

  Sets the pvalue member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_pvalue(
  MATCHED_ELEMENT_T *element,
  double pvalue
) {
  assert(element != NULL);
  element->pvalue = pvalue;
  element->has_pvalue = true;
}

/**********************************************************************
  has_matched_element_pvalue

  Does a matched_element object have a pvalue?
**********************************************************************/
bool has_matched_element_pvalue(MATCHED_ELEMENT_T *element) {
  return element->has_pvalue;
}

/**********************************************************************
  get_matched_element_pvalue

  Gets the pvalue member from a cisml matched_element object.
**********************************************************************/
double get_matched_element_pvalue(MATCHED_ELEMENT_T* element) {
  assert(element != NULL);
  return element->pvalue;
}

/**********************************************************************
  set_matched_element_qvalue

  Sets the qvalue member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_qvalue(
  MATCHED_ELEMENT_T *element,
  double qvalue
) {
  assert(element != NULL);
  element->qvalue = qvalue;
  element->has_qvalue = true;
}

/**********************************************************************
  has_matched_element_qvalue

  Does a matched_element object have a qvalue?
**********************************************************************/
bool has_matched_element_qvalue(MATCHED_ELEMENT_T *element) {
  return element->has_qvalue;
}

/**********************************************************************
  get_matched_element_qvalue

  Gets the qvalue member from a cisml matched_element object.
**********************************************************************/
double get_matched_element_qvalue(MATCHED_ELEMENT_T* element) {
  assert(element != NULL);
  return element->qvalue;
}

/**********************************************************************
  set_matched_element_clusterid

  Sets the clusterid member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_clusterid(
  MATCHED_ELEMENT_T *element,
  char *clusterid
) {

  assert(element != NULL);

  if (clusterid == NULL) {
    if (element->clusterid != NULL) {
      myfree(element->clusterid);
    }
    element->clusterid = NULL;
  }
  else {
    int old_length = 0;
    int new_length = strlen(clusterid) + 1;
    if (element->clusterid != NULL) {
      old_length = strlen(element->clusterid) + 1;
    }
    if (old_length < new_length) {
      element->clusterid =
        mm_realloc(element->clusterid, new_length);
    }
    strncpy(element->clusterid, clusterid, new_length);
  }

}

/**********************************************************************
  get_matched_element_clusterid

  Gets the clusterid member from a cisml matched_element object.
  May return NULL.
**********************************************************************/
char *get_matched_element_clusterid(MATCHED_ELEMENT_T* element) {

  assert(element != NULL);

  return element->clusterid;

}

/**********************************************************************
  get_matched_element_sequence

  Gets the sequence member from a cisml matched_element object.
  Caller should not free the string.
  May return NULL.
**********************************************************************/
const char *get_matched_element_sequence(MATCHED_ELEMENT_T* element) {

  assert(element != NULL);

  return (const char *) element->sequence;

}

/**********************************************************************
  set_matched_element_sequence

  Sets the sequence member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_sequence(MATCHED_ELEMENT_T* element, char *seq) {

  assert(element != NULL);

  if (element->sequence != NULL) {
    myfree(element->sequence);
  }
  element->sequence = strndup(seq, abs(element->stop - element->start) + 1);

}

/**********************************************************************
  set_matched_element_strand

  Sets the strand member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_strand(MATCHED_ELEMENT_T* element, char strand) {

  assert(element != NULL);

  element->strand = strand;

}

/**********************************************************************
  get_matched_element_strand

  Gets the strand member in a cisml matched_element object.
**********************************************************************/
char get_matched_element_strand(MATCHED_ELEMENT_T* element) {

  assert(element != NULL);

  return element->strand;

}

/**********************************************************************
  get_matched_element_scanned_seq

  Returns the parent_scanned_sequence member in a matched_element object.
**********************************************************************/
SCANNED_SEQUENCE_T *get_matched_element_scanned_seq(MATCHED_ELEMENT_T* element) {

  assert(element != NULL);

  return element->parent_sequence;

}

/**********************************************************************
  compare_matched_elements

  Compare two objects, treating them as matched elements
**********************************************************************/
static int compare_matched_elements(void *p1, void *p2) {

  MATCHED_ELEMENT_T *e1 = (MATCHED_ELEMENT_T *) p1;
  MATCHED_ELEMENT_T *e2 = (MATCHED_ELEMENT_T *) p2;

  if (e1->pvalue < e2->pvalue) {
    return 1;
  }
  else if (e1->pvalue > e2->pvalue) {
    return -1;
  }
  else {
    // If p-values are equal, compare staring postions
    // to break the time.
    return e1->start < e2->start ? 1: -1;
  }

}

/**********************************************************************
  copy_matched_element

  Copy an object, treating it as a matched element
  For use with HEAP.
**********************************************************************/
void *copy_matched_element(void *p) {

  MATCHED_ELEMENT_T *e = (MATCHED_ELEMENT_T *) p;
  MATCHED_ELEMENT_T *new_e = allocate_matched_element_with_score(
      e->start,
      e->stop,
      e->score,
      e->pvalue,
      e->parent_sequence
    );
  new_e->qvalue = e->qvalue;
  new_e->clusterid = e->clusterid;
  new_e->sequence = e->sequence;
  new_e->strand = e->strand;
  return new_e;

}

/**********************************************************************
  destroy_matched_element

  Destroy an object, treating it as a matched element
  For use with HEAP.
**********************************************************************/
static void destroy_matched_element(void *p) {

  MATCHED_ELEMENT_T *e = (MATCHED_ELEMENT_T *) p;
  free_matched_element(e);

}

/*********************************************************************
  print_cisml_matched_elements()

  Print the XML for the CisML matched elements under a scanned sequence.
**********************************************************************/
void print_cisml_matched_elements(
  CISML_T *cisml,
  FILE *out,
  int num_matched_elements,
  MATCHED_ELEMENT_T **elements
) {

  double qthresh = get_cisml_site_qvalue_cutoff(cisml);
  double pthresh = get_cisml_site_pvalue_cutoff(cisml);
  bool have_output_sequence_start = false;
  STR_T *buffer;
  buffer = str_create(10);

  // We're going to output at least one matched element
  // for this sequence. Output the scanned sequence
  // starting XML tag.
  if (!have_output_sequence_start) {
    have_output_sequence_start = true;
    // print_cisml_scanned_sequence_start(cisml, out, element->parent_sequence);
  }

  int i_element = 0;
  for(i_element = 0; i_element < num_matched_elements; i_element++) {

    MATCHED_ELEMENT_T *element = elements[i_element];

    if (element->pvalue > pthresh || element->qvalue > qthresh) {
      continue;
    }
    (cisml->num_passing_cutoff)++;

    fprintf(
      out,
      "<matched-element start=\"%d\" stop=\"%d\"",
      get_matched_element_start(element),
      get_matched_element_stop(element)
    );
    if (has_matched_element_score(element)) {
      double score = get_matched_element_score(element);
      fprintf(out, " score=\"%g\"", score);
    }
    if (has_matched_element_pvalue(element)) {
      double pvalue = get_matched_element_pvalue(element);
      fprintf(out, " pvalue=\"%.3g\"", pvalue);
    }
    char *clusterid = get_matched_element_clusterid(element);
    if (clusterid != NULL) {
      fprintf(out, " clusterid=\"%s\"", xmlify(clusterid, buffer, true));
    }
    fprintf(out, ">\n");
    const char *sequence = get_matched_element_sequence(element);
    if (sequence !=  NULL) {
      fprintf(out, "<sequence>%s</sequence>\n", xmlify(sequence, buffer, false));
    }
    if (has_matched_element_qvalue(element)) {
      double qvalue = get_matched_element_qvalue(element);
      fprintf(out, "<mem:qvalue>%.3g</mem:qvalue>\n", qvalue);
    }
    fputs("</matched-element>\n", out);

  }

  if (have_output_sequence_start) {
    // We output the scanned sequence start XML tag.
    // Output the scanned sequence closing tag.
    // print_cisml_scanned_sequence_end(out);
  }
  str_destroy(buffer, false);
}

/**********************************************************************
  print_cisml_scanned_sequence_start()

  Print the starting XML tag for a CisML scanned sequence element.
**********************************************************************/
void print_cisml_scanned_sequence_start(
  CISML_T *cisml,
  FILE *out,
  SCANNED_SEQUENCE_T *seq
) {
  STR_T *buffer;
  buffer = str_create(10);

    fprintf(
      out,
      "<scanned-sequence accession=\"%s\" ",
      xmlify(get_scanned_sequence_accession(seq), buffer, true)
    );
    fprintf(
      out,
      "name=\"%s\"",
      xmlify(get_scanned_sequence_name(seq), buffer, true)
    );

    if (has_scanned_sequence_score(seq)) {
      double score = get_scanned_sequence_score(seq);
      fprintf(out, " score=\"%g\"", score);
    }
    if (has_scanned_sequence_pvalue(seq)) {
      double pvalue = get_scanned_sequence_pvalue(seq);
      fprintf(out, " pvalue=\"%g\"", pvalue);
    }
    if (has_scanned_sequence_length(seq)) {
      int length = get_scanned_sequence_length(seq);
      fprintf(out, " length=\"%d\"", length);
    }
    char *db = get_scanned_sequence_db(seq);
    if (db != NULL) {
      fprintf(out, " db=\"%s\"", xmlify(db, buffer, true));
    }
    char *lsid = get_scanned_sequence_lsid(seq);
    if (lsid != NULL) {
      fprintf(out, " lsid=\"%s\"", xmlify(lsid, buffer, true));
    }
    fprintf(out, ">\n");
  str_destroy(buffer, false);
}

/**********************************************************************
  print_cisml_scanned_sequence_end()

  Print the ending XML tag for a CisML scanned sequence element.
**********************************************************************/
void print_cisml_scanned_sequence_end(FILE *out) {

    fputs("</scanned-sequence>\n", out);

}

/**********************************************************************
  print_cisml_scanned_sequences()

  Print the XML for a CisML scanned sequence element.
**********************************************************************/
void print_cisml_scanned_sequences(
  CISML_T *cisml,
  FILE *out,
  int num_seqs,
  SCANNED_SEQUENCE_T **sequences
) {

  int i_seq = 0;
  for(i_seq = 0; i_seq < num_seqs; i_seq++) {

    SCANNED_SEQUENCE_T *seq = sequences[i_seq];

    print_cisml_scanned_sequence_start(cisml, out, seq);

    // Skip sequences with no matched elements
    if (seq->num_matched_elements == 0) {
      print_cisml_scanned_sequence_end(out);
      continue;
    }

    // We don't know if any matched elements pass the
    // p-value or q-value  threshold, so the first
    // matched elemenet printed will trigger the output
    // for the scanned sequence tags.
    print_cisml_matched_elements(
      cisml,
      out,
      seq->num_matched_elements,
      seq->elements
    );

    print_cisml_scanned_sequence_end(out);

  }

}

/**********************************************************************
  print_cisml_start_pattern

  Print the starting tag for a CisML pattern
**********************************************************************/
void print_cisml_start_pattern(
  CISML_T *cisml,
  FILE *out,
  PATTERN_T *pattern
) {
  STR_T *buffer;
  buffer = str_create(10);

  fprintf(
    out,
    "<pattern accession=\"%s\"",
    xmlify(get_pattern_accession(pattern), buffer, true)
  );
  fprintf(
    out,
    " name=\"%s\"",
    xmlify(get_pattern_name(pattern), buffer, true)
  );
  if (has_pattern_score(pattern)) {
    double score = get_pattern_score(pattern);
    fprintf(out, " score=\"%g\"", score);
  }
  if (has_pattern_pvalue(pattern)) {
    double pvalue = get_pattern_pvalue(pattern);
    fprintf(out, " pvalue=\"%g\"", pvalue);
  }
  char *db = get_pattern_db(pattern);
  if (db != NULL) {
    fprintf(out, " db=\"%s\"", xmlify(db, buffer, true));
  }
  char *lsid = get_pattern_lsid(pattern);
  if (lsid != NULL) {
    fprintf(out, " lsid=\"%s\"", xmlify(lsid, buffer, true));
  }
  fputs(">\n", out);
  str_destroy(buffer, false);
}

/**********************************************************************
  print_cisml_end_pattern

  Print the ending tag for a CisML pattern
**********************************************************************/
void print_cisml_end_pattern(FILE *out) {
  fputs("</pattern>\n", out);
}

/**********************************************************************
  print_cisml_patterns

  Print pattern elements for CisML
**********************************************************************/
void print_cisml_patterns(
  CISML_T *cisml,
  FILE *out,
  int num_patterns,
  PATTERN_T **patterns
) {

  int i = 0;
  for(i = 0; i < num_patterns; i++) {
    int num_seq = 0;
    num_seq = get_pattern_num_scanned_sequences(patterns[i]);
    // only patterns with sequences can be printed to be conform with the DTD
    if (num_seq > 0) {
        print_cisml_start_pattern(cisml, out, patterns[i]);
        SCANNED_SEQUENCE_T **sequences
        = get_pattern_scanned_sequences(patterns[i]);
      print_cisml_scanned_sequences(cisml, out, num_seq, sequences);

      if (has_pattern_qvalues(patterns[i]) == true) {
          fputs("<mem:has-qvalues>yes</mem:has-qvalues>\n", out);
      }

      print_cisml_end_pattern(out);
    }
  }
}

/*********************************************************************
  print_cisml_multi_pattern_match()

  Print the XML for the multi-pattern match element under a scanned sequence.
**********************************************************************/
void print_cisml_multi_pattern_match(
  FILE *out,
  MULTI_PATTERN_MATCH_T *match
) {
  STR_T *buffer;
  buffer = str_create(10);

  fprintf(
    out,
    "<mem:match cluster-id=\"%s\" ",
    xmlify(match->clusterid, buffer, true)
  );
  fprintf(
    out,
    "seq-name=\"%s\" start=\"%d\" stop=\"%d\" "
    "evalue=\"%3.1g\" qvalue=\"%3.1g\">",
    xmlify(match->seq_name, buffer, true),
    match->start,
    match->stop,
    match->evalue,
    match->qvalue
  );
  if (match->sequence != NULL) {
    if (strlen(match->sequence) != (match->stop - match->start)) {
      fprintf(
        stderr, 
        "Match sequence fault! seq start = %d seq stop = %d, seq_len = %zd\n",
        match->start,
        match->stop,
        strlen(match->sequence)
      );
    }
    fprintf(out, "%s\n", xmlify(match->sequence, buffer, false));
  }
  fputs("</mem:match>\n", out);
  str_destroy(buffer, false);
}

/**********************************************************************
  print_cisml_multi_patterns

  Print multi_patterns element for CisML
**********************************************************************/
void print_cisml_multi_patterns(
  CISML_T *cisml,
  FILE *out,
  int num_multi_patterns,
  MULTI_PATTERN_T **multi_patterns
) {

  int i = 0;
  for(i = 0; i < num_multi_patterns; i++) {
    double *pvalue_cutoff = cisml->site_pvalue_cutoff;
    double *pvalue = multi_patterns[i]->pvalue;
    if (pvalue != NULL && pvalue_cutoff != NULL && *pvalue > *pvalue_cutoff) {
      continue;
    }
    fprintf(out, "<multi-pattern-scan");
    if (has_multi_pattern_score(multi_patterns[i])) {
      fprintf( out, " score=\"%g\"", *multi_patterns[i]->score);
    }
    if (has_multi_pattern_pvalue(multi_patterns[i])) {
      fprintf( out, " pvalue=\"%3.1g\"", *pvalue);
    }
    fprintf(out, ">\n");
    int num_patterns = get_multi_pattern_num_patterns(multi_patterns[i]);
    if (num_patterns > 0) {
      PATTERN_T **patterns = get_multi_pattern_patterns(multi_patterns[i]);
      print_cisml_patterns(cisml, out, num_patterns, patterns);
    }
    if (multi_patterns[i]->match != NULL) {
      print_cisml_multi_pattern_match(out, multi_patterns[i]->match);
    }
    fprintf(out, "</multi-pattern-scan>\n");
  }
}

/**********************************************************************
  print_cisml_parmeters

  Print parameters element for CisML
**********************************************************************/
void print_cisml_parameters(FILE *out, CISML_T *cisml) {
  STR_T *buffer;

  buffer = str_create(10);

  fprintf(out, "<parameters>\n");
  fprintf(
    out,
    "<command-line>%s</command-line>\n",
    xmlify(get_cisml_command_line(cisml), buffer, false)
  );
  fprintf(
    out,
    "<pattern-file>%s</pattern-file>\n",
    xmlify(get_cisml_pattern_file(cisml), buffer, false)
  );
  fprintf(
    out,
    "<sequence-file>%s</sequence-file>\n",
    xmlify(get_cisml_sequence_file(cisml), buffer, false)
  );
  char *background_file = get_cisml_background_file(cisml);
  if (background_file != NULL) {
    fprintf(
      out,
      "<background-seq-file>%s</background-seq-file>\n",
      xmlify(background_file, buffer, false)
    );
  }
  if (has_cisml_pattern_pvalue_cutoff(cisml)) {
    double pvalue = get_cisml_pattern_pvalue_cutoff(cisml);
    fprintf(
      out,
      "<pattern-pvalue-cutoff>%g</pattern-pvalue-cutoff>\n",
      pvalue
    );
  }
  if (has_cisml_sequence_pvalue_cutoff(cisml)) {
    double pvalue = get_cisml_sequence_pvalue_cutoff(cisml);
    fprintf(
      out,
      "<sequence-pvalue-cutoff>%g</sequence-pvalue-cutoff>\n",
      pvalue
    );
  }
  if (has_cisml_site_pvalue_cutoff(cisml)) {
    double pvalue = get_cisml_site_pvalue_cutoff(cisml);
    fprintf(
      out,
      "<site-pvalue-cutoff>%g</site-pvalue-cutoff>\n",
      pvalue
    );
  }
  char *filter = get_cisml_sequence_filter(cisml);
  if (filter != NULL) {
    fprintf(
      out,
      "<sequence-filtering on-off=\"on\" type=\"%s\" />\n",
      xmlify(filter, buffer, true)
    );
  }
  else {
    fprintf(
      out,
      "<sequence-filtering on-off=\"off\"/>\n"
    );
  }

  // Print MEME extensions for CisML


  fprintf(out, "</parameters>\n" );
  str_destroy(buffer, false);
}

/**********************************************************************
  print_cisml_xml_header

  Print the DTD and XML header for CisML
**********************************************************************/
static void print_cisml_xml_header(FILE* out, const char *stylesheet) {
  STR_T *buffer;
  buffer = str_create(10);
  if (stylesheet == NULL){
          fputs(cisml_dts,out);
  } else {
          fputs("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n", out);
          fputs("<?xml-stylesheet type=\"text/xsl\" ", out);
          fprintf(out, "href=\"%s\"?>\n", xmlify(stylesheet, buffer, true));
          fputs("<!-- Begin document body -->\n", out);
  }
  str_destroy(buffer, false);
}

/**********************************************************************
  print_cisml_start

  Print the opening section of the CisML XML
**********************************************************************/
void print_cisml_start(FILE* out, char *program_name, bool print_header,
                const char *stylesheet, bool print_namespace) {
  STR_T *buffer;
  buffer = str_create(10);

  if (print_header == true) {
    print_cisml_xml_header(out, NULL);
  }
  fputs("<cis-element-search\n", out);
  fputs("  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"", out);
  fputs("\n", out);
  fputs("  xsi:schemaLocation=", out);
  fputs("\"http://zlab.bu.edu/schema/cisml cisml.xsd\"\n", out);
  // FIXME the next line causes problems for the cisml reader
  if (print_namespace){
          fputs("  xmlns=\"http://zlab.bu.edu/schema/cisml\"\n", out);
  }
  fputs("  xmlns:mem=\"http://noble.gs.washington.edu/meme\"\n>\n", out);

  fprintf(
    out,
    "<program-name>%s</program-name>\n",
    xmlify(program_name, buffer, false)
  );
  str_destroy(buffer, 0);
}

/**********************************************************************
  print_cisml_end

  Print the closing section of the CisML XML
**********************************************************************/
void print_cisml_end(FILE* out) {
  fprintf(out, "</cis-element-search>\n");
}

/**********************************************************************
  print_cisml

  Print the cisml data structure as CisML XML
**********************************************************************/
void print_cisml(FILE* out, CISML_T *cisml, bool print_header,
                const char *stylesheet, bool print_namespace) {

  char *program_name = get_cisml_program_name(cisml);
  print_cisml_start(out, program_name, print_header, NULL, print_namespace);
  print_cisml_parameters(out, cisml);

  int num_multi_patterns = get_cisml_num_multi_patterns(cisml);
  if (num_multi_patterns > 0) {
    MULTI_PATTERN_T **multi_patterns = get_cisml_multi_patterns(cisml);
    print_cisml_multi_patterns(cisml, out, num_multi_patterns, multi_patterns);
  }

  int num_patterns = get_cisml_num_patterns(cisml);
  if (num_patterns > 0) {
    PATTERN_T **patterns = get_cisml_patterns(cisml);
    print_cisml_patterns(cisml, out, num_patterns, patterns);
  }

  print_cisml_end(out);

}

/*************************************************************************
 * Compare two matched-elements by pvalue for 'qsort'.
 *************************************************************************/
static int matched_elements_compare_by_pvalue
  (const void* elem1,
   const void* elem2)
{
  const double key1 = ((MATCHED_ELEMENT_T *) elem1)->pvalue;
  const double key2 = ((MATCHED_ELEMENT_T *) elem2)->pvalue;
  const int start1 = ((MATCHED_ELEMENT_T *)elem1)->start;
  const int start2 = ((MATCHED_ELEMENT_T *)elem2)->start;

  if (key1 < key2) {
    return(-1);
  } else if (key1 > key2) {
    return(1);
  } else if (start1 < start2) {
    return(-1);
  } else if (start1 > start2) {
    return(1);
  }
  return(0);
}

/*************************************************************************
 * sort_matched_elements
 *
 * Sort a an array of pointers to matched-elements sites by pvalue
 * or by sequence name and position.
 *************************************************************************/
void sort_matched_elements(
  bool sort_by_pvalue,
  int num_elements,
  MATCHED_ELEMENT_T **elements
) {

  assert(elements != NULL);

  // Tell the user what's up.
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Sorting %d matched elements ", num_elements);
    if (sort_by_pvalue) {
      fprintf(stderr, "by p-value.\n");
    } else {
      fprintf(stderr, "by sequence name and start position.\n");
    }
  }

  if (sort_by_pvalue) {
    qsort(
            (void *) elements,
            num_elements,
            sizeof(MATCHED_ELEMENT_T *),
            matched_elements_compare_by_pvalue
          );
  } else {
    /*
    qsort(
            (void *) elements,
            num_elements,
            sizeof(MATCHED_ELEMENT_T *),
            matched_elements_compare_by_position
          );
    */
  }
}

/*************************************************************************
 * Calculate the q-values for each matched-element in this pattern
 * from the p-values.
 *************************************************************************/
void pattern_calculate_qvalues(PATTERN_T *pattern, ARRAY_T *sampled_pvalues) {

  assert(pattern != NULL);
  assert(pattern->is_complete == true);

  int num_stored_matches = get_pattern_num_stored_matches(pattern);
  long num_scanned_positions = get_pattern_num_scanned_positions(pattern);
  if (verbosity >= HIGH_VERBOSE) {
    fprintf(stderr, "Num stored matches %d\n", num_stored_matches);
    fprintf(stderr, "Num scanned positions %ld\n", num_scanned_positions);
  }
  // Tell the user what's up.
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Computing q-values.\n");
  }

  if (num_stored_matches) {

    // Extract the p-values into an array.
    int i_element = 0;
    ARRAY_T* pvalues = allocate_array(num_stored_matches);
    MATCHED_ELEMENT_T *element = NULL;
    for (i_element = 0; i_element < num_stored_matches; i_element++) {
      element = pattern->elements[i_element];
      set_array_item(i_element, element->pvalue, pvalues);
    }

    // Convert them to q-values.
    compute_qvalues(
      false, // Don't stop with FDR
      true, // Try to esimate pi0
      NULL, // Don't store pi-zero in a file.
      NUM_BOOTSTRAPS,
      NUM_BOOTSTRAP_SAMPLES,
      NUM_LAMBDA,
      MAX_LAMBDA,
      num_scanned_positions,
      pvalues,
      sampled_pvalues
    );

    // Update the matched elements with the q-values.
    for (i_element = 0; i_element < num_stored_matches; i_element++) {
      set_matched_element_qvalue(
        pattern->elements[i_element],
        get_array_item(i_element, pvalues)
      );
    }

    free_array(pvalues);

    // Sort by sequence ID and position.
    // Since we are putting out XML there is no need to sort
    // matched-elements back into positon order.
  }

  pattern->qvalues_computed = true;

}

/*************************************************************************
 * Compare two matched-elements by pvalue for 'qsort'.
 *************************************************************************/
static int multi_pattern_compare_by_pvalue
  (const void* mp1,
   const void* mp2)
{
  const double key1 = *((*((MULTI_PATTERN_T **) mp1))->pvalue);
  const double key2 = *((*((MULTI_PATTERN_T **) mp2))->pvalue);

  if (key1 < key2) {
    return(-1);
  } 
  else if (key1 > key2) {
    return(1);
  }
  else {
    return(0);
  }
}

/*************************************************************************
 * sort_multi_patterns
 *
 * Sort a an array of pointers to multi-patters by pvalue
 *************************************************************************/
void sort_multi_patterns(
  int num_multi_patterns,
  MULTI_PATTERN_T **multi_patterns
) {

  assert(multi_patterns != NULL);

  // Tell the user what's up.
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(
      stderr, 
      "Sorting %d matches "
      "by p-value.\n", 
      num_multi_patterns
    );
  }

  qsort(
    (void *) multi_patterns,
    num_multi_patterns,
    sizeof(MULTI_PATTERN_T *),
    multi_pattern_compare_by_pvalue
  );

}

/*************************************************************************
 * Calculate the q-values corresponding from the the p-values of the
 * multi-patterns
 *************************************************************************/
void multi_pattern_calculate_qvalues(
  int num_multi_patterns,
  MULTI_PATTERN_T **multi_patterns, 
  ARRAY_T *sampled_pvalues
) {

  assert(multi_patterns != NULL);

  // Tell the user what's up.
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Computing q-values.\n");
  }

  // To calculate q-value pvalue must be sorted
  sort_multi_patterns(num_multi_patterns, multi_patterns);

  // Extract the p-values into an array.
  ARRAY_T* pvalues = allocate_array(num_multi_patterns);
  MULTI_PATTERN_T *mp = NULL;
  int i_mp = 0;
  for (i_mp = 0; i_mp < num_multi_patterns; i_mp++) {
    mp = multi_patterns[i_mp];
    set_array_item(i_mp, get_multi_pattern_pvalue(mp), pvalues);
  }

  // Convert them to q-values.
  compute_qvalues(
    false, // Don't stop with FDR
    true, // Try to esimate pi0
    NULL, // Don't store pi-zero in a file.
    NUM_BOOTSTRAPS,
    NUM_BOOTSTRAP_SAMPLES,
    NUM_LAMBDA,
    MAX_LAMBDA,
    num_multi_patterns,
    pvalues,
    sampled_pvalues
  );

  // Update the multi-patterns with the q-values.
  for (i_mp = 0; i_mp < num_multi_patterns; i_mp++) {
    set_multi_pattern_match_qvalue(
      multi_patterns[i_mp]->match,
      get_array_item(i_mp, pvalues)
    );
  }

  free_array(pvalues);

}

/**********************************************************************
 * This function saves CisML results as a set of files in a
 * directory. The file names are provided by the input parameters:
 *   cisml is a pointer to the ciml structure to be printed.
 *   ouput_dirname will be the name of the output directory
 *   xml_filename will be the CisML output
 *   html_filename will be the name of the HTML output
 *   text_filename will be the name of the plain text output
 *   gff_filename will be the name of the GFF output
 *   allow_clobber will determine whether or not existing files will
 *                 be overwritten.
 *   print_namespace will determine whether or not the standard name
 *                              space will be written
 *********************************************************************/
bool print_full_results(
  CISML_T *cisml,
  char *output_dirname,
  char *xml_filename,
  char *html_filename,
  char *text_filename,
  char *gff_filename,
  bool allow_clobber,
  bool print_namespace
) {

  bool success = true;

  // Create output directory
  if (create_output_directory(
       output_dirname,
       allow_clobber,
       false /* Don't print warning messages */
      )
    ) {
    // Failed to create output directory.
    die("Couldn't create output directory %s.\n", output_dirname);
  }
  // Create the paths to the output files
  const char* HTML_STYLESHEET = "cisml-to-html.xsl";
  const char* CSS_STYLESHEET = "cisml.css";
  const char* GFF_STYLESHEET = "cisml-to-gff.xsl";
  const char* TEXT_STYLESHEET = "cisml-to-text.xsl";
  char *html_stylesheet_path = make_path_to_file(get_meme_data_dir(), HTML_STYLESHEET);
  char *css_stylesheet_path = make_path_to_file(get_meme_data_dir(), CSS_STYLESHEET);
  char *text_stylesheet_path = make_path_to_file(get_meme_data_dir(), TEXT_STYLESHEET);
  char *gff_stylesheet_path = make_path_to_file(get_meme_data_dir(), GFF_STYLESHEET);
  char *xml_path = make_path_to_file(output_dirname, xml_filename);
  char *html_path = make_path_to_file(output_dirname, html_filename);
  char *text_path = make_path_to_file(output_dirname, text_filename);
  char *gff_path = make_path_to_file(output_dirname, gff_filename);
  char *html_stylesheet_copy_path = make_path_to_file(output_dirname, HTML_STYLESHEET);
  char *css_stylesheet_copy_path = make_path_to_file(output_dirname, CSS_STYLESHEET);
  FILE *xml_file = fopen(xml_path, "w");
  if (!xml_file) {
    die("Couldn't open file %s for output.\n", xml_path);
  }

  // Output XML
  print_cisml(xml_file, cisml, true, HTML_STYLESHEET,print_namespace);
  fclose(xml_file);

  // Output HTML
  success = print_xml_filename_to_filename_using_stylesheet(xml_path, html_stylesheet_path, html_path) && success;

  // Output text
  success = print_xml_filename_to_filename_using_stylesheet(xml_path, text_stylesheet_path, text_path) && success;

  // Output GFF
  success = print_xml_filename_to_filename_using_stylesheet(xml_path, gff_stylesheet_path, gff_path) && success;

  // Copy XML to HTML and CSS stylesheets to output directory
  copy_file(html_stylesheet_path, html_stylesheet_copy_path);
  copy_file(css_stylesheet_path, css_stylesheet_copy_path);

  myfree(html_stylesheet_path);
  myfree(html_stylesheet_copy_path);
  myfree(css_stylesheet_path);
  myfree(css_stylesheet_copy_path);
  myfree(text_stylesheet_path);
  myfree(gff_stylesheet_path);
  myfree(xml_path);
  myfree(html_path);
  myfree(text_path);
  myfree(gff_path);

  return(success);
}

/**********************************************************************
 * This function saves the CisML results as a plain text to stdout.
 *********************************************************************/
bool print_cisml_as_text(CISML_T *cisml) {

  const char* etc_dir = get_meme_data_dir();
  bool success = true;

  // Create temp file for CisML output.
  char tmp_filename[] = "CISMLXXXXXX";
  int fd = mkstemp(tmp_filename);
  if (fd == -1) {
    die("Couldn't create temporary file for text results\n");
  }
  FILE *xml_file = fdopen(fd, "w");
  if (!xml_file) {
    die("Couldn't open file %s for output.\n", "fimo.xml");
  }

  // Output CisML
  print_cisml(xml_file, cisml, false, NULL, true);

  fclose(xml_file);


  // Output text using stylesheet that converts CisML temp file to text.
  const char* TEXT_STYLESHEET = "cisml-to-text.xsl";
  char *text_stylesheet_path = make_path_to_file(etc_dir, TEXT_STYLESHEET);
  success = print_xml_filename_to_file_using_stylesheet(tmp_filename, text_stylesheet_path, stdout);


  myfree(text_stylesheet_path);
  close(fd);
  int result = remove(tmp_filename);
  if (result == -1) {
    fprintf(stderr, "Couldn't remove temporary file %s.\n", tmp_filename);
    success = false;
  }

  return(success);
}

/**********************************************************************
  Cisml Parser state
**********************************************************************/
typedef struct cismlp CISMLP_T;
struct cismlp {
  CISML_T *cisml;
  MULTI_PATTERN_T *multi;
  PATTERN_T *pattern;
  SCANNED_SEQUENCE_T *scan;
  MATCHED_ELEMENT_T *match;
  bool done;
};

/**********************************************************************
  Cisml Parser callbacks
**********************************************************************/

void cismlp_start_cisml(void* status) {
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  data->cisml = mm_malloc(sizeof(CISML_T));
  memset(data->cisml, 0, sizeof(CISML_T));
  data->done = false;
}
void cismlp_start_cis_element_search(void* status) {
  // <status>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
}
void cismlp_handle_program_name(void* status, char* program_name) {
  // <status> <program name>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  data->cisml->program_name = strdup(program_name);
}
void cismlp_start_parameters(void* status) {
  // <status>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
}
void cismlp_handle_command_line(void* status, char* command_line) {
  // <status> <command line>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  data->cisml->command_line = strdup(command_line);
}
void cismlp_handle_pattern_file(void* status, char* pattern_file) {
  // <status> <pattern file>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  data->cisml->pattern_file = strdup(pattern_file);
}
void cismlp_handle_sequence_file(void* status, char* sequence_file) {
  // <status> <sequence file>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  data->cisml->sequence_file = strdup(sequence_file);
}
// optional element
void cismlp_handle_background_seq_file(void* status, char* background_sequence_file) {
  // <status> <background sequence file>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  set_cisml_background_file(data->cisml, background_sequence_file);
}
// optional element
void cismlp_handle_pattern_pvalue_cutoff(void* status, double pattern_pvalue_cutoff) {
  // <status> <pattern pvalue cutoff>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  set_cisml_pattern_pvalue_cutoff(data->cisml, pattern_pvalue_cutoff);
}
// optional element
void cismlp_handle_sequence_pvalue_cutoff(void* status, double sequence_pvalue_cutoff) {
  // <status> <sequence pvalue cutoff>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  set_cisml_sequence_pvalue_cutoff(data->cisml, sequence_pvalue_cutoff);
}
// optional element
void cismlp_handle_site_pvalue_cutoff(void* status, double site_pvalue_cutoff) {
  // <status> <site pvalue cutoff>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  set_cisml_site_pvalue_cutoff(data->cisml, site_pvalue_cutoff);
}
void cismlp_handle_sequence_filtering(void* status, int filtering_used, char* filtering_program) {
  // <status> <filtering used?> [<filtering program>]    
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  if (filtering_used) set_cisml_sequence_filter(data->cisml, filtering_program != NULL ? filtering_program : "");
}
void cismlp_end_parameters(void* status) {
  // <status>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
}
void cismlp_start_multi_pattern_scan(void* status, double* pvalue, double* score) {
  // <status> [<pvalue>] [<score>]
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  // create the multi pattern scan but don't add it to the cisml yet
  assert(data->multi == NULL);
  data->multi = allocate_multi_pattern();
  if (pvalue != NULL) set_multi_pattern_pvalue(data->multi, *pvalue);
  if (score != NULL) set_multi_pattern_score(data->multi, *score);
}
void cismlp_start_pattern(void* status, char* accession, char* name, char* db,
    char* ls_id, double* pvalue, double* score) {
  // <status> <accession> <name> [<db>] [<lsId>] [<pvalue>] [<score>]
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  // create the pattern but don't add it to the cisml/multi yet
  assert(data->pattern == NULL);
  data->pattern = allocate_pattern(accession, name);
  set_pattern_db(data->pattern, db);
  set_pattern_lsid(data->pattern, ls_id);
  if (pvalue != NULL) set_pattern_pvalue(data->pattern, *pvalue);
  if (score != NULL) set_pattern_score(data->pattern, *score);
}
void cismlp_start_scanned_sequence(void* status, char* accession, char* name,
    char* db, char* ls_id, double* score, double* pvalue, long* length) {
  // <status> <accession> <name> [<db>] [<lsId>] [<score>] [<pvalue>] [<length>]
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  // create the scanned sequence but don't add it to the pattern yet
  assert(data->scan == NULL);
  assert(data->pattern != NULL);
  data->scan = allocate_scanned_sequence(accession, name, data->pattern);
  set_scanned_sequence_db(data->scan, db);
  set_scanned_sequence_lsid(data->scan, ls_id);
  if (score != NULL) set_scanned_sequence_score(data->scan, *score);
  if (pvalue != NULL) set_scanned_sequence_pvalue(data->scan, *pvalue);
  if (length != NULL) set_scanned_sequence_length(data->scan, (int)*length);
}
// read the matched element
void cismlp_start_matched_element(void* status, long start, long stop,
    double* score, double* pvalue, char* cluster_id) {
  // <status> <start> <stop> [<score>] [<pvalue>] [<clusterId>]
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  // create the matched element
  assert(data->match == NULL);
  assert(data->scan != NULL);
  data->match = allocate_matched_element((int)start, (int)stop, data->scan);
  if (score != NULL) set_matched_element_score(data->match, *score);
  if (pvalue != NULL) set_matched_element_pvalue(data->match, *pvalue);
  set_matched_element_clusterid(data->match, cluster_id);
}
void cismlp_handle_sequence(void* status, char* sequence) {
  // <status> <sequence>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  assert(data->match != NULL);
  set_matched_element_sequence(data->match, sequence);
}
void cismlp_end_matched_element(void* status) {
  // <status>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  assert(data->match != NULL);
  assert(data->scan != NULL);
  add_scanned_sequence_matched_element(data->scan, data->match);
  data->match = NULL;
}
void cismlp_end_scanned_sequence(void* status) {
  // <status>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  assert(data->scan != NULL);
  assert(data->pattern != NULL);
  // shrink the scan elements array to required size
  compact_scanned_sequence_matched_elements(data->scan);
  data->scan = NULL;
}
void cismlp_end_pattern(void* status) {
  // <status>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  assert(data->pattern != NULL);
  data->pattern->is_complete = true;
  if (data->multi != NULL) {
    add_multi_pattern_pattern(data->multi, data->pattern);
  } else {
    assert(data->cisml != NULL);
    add_cisml_pattern(data->cisml, data->pattern);
  }
  data->pattern = NULL;
}
void cismlp_end_multi_pattern_scan(void* status) {
  // <status>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  assert(data->multi != NULL);
  assert(data->cisml != NULL);
  add_cisml_multi_pattern(data->cisml, data->multi);
  data->multi = NULL;
}
void cismlp_end_cis_element_search(void* status) {
  // <status>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
}
void cismlp_end_cisml(void* status) {
  // <status>
  CISMLP_T *data;
  data = (CISMLP_T*)status;
  assert(data->cisml != NULL);
  assert(data->multi == NULL);
  assert(data->pattern == NULL);
  assert(data->scan == NULL);
  assert(data->match == NULL);
  data->done = true;
}

/**********************************************************************
  read_cisml

  Reads in a CisML XML file and create the cisml_t data structure
**********************************************************************/
CISML_T* read_cisml(const char* cisml_filename) {
  CISMLP_T data;
  CISML_CALLBACKS_T callbacks;
  // set data
  memset(&data, 0, sizeof(CISMLP_T));
  // set callbacks
  memset(&callbacks, 0, sizeof(CISML_CALLBACKS_T));
  callbacks.start_cisml = cismlp_start_cisml;
  callbacks.start_cis_element_search = cismlp_start_cis_element_search;
  callbacks.handle_program_name = cismlp_handle_program_name;
  callbacks.start_parameters = cismlp_start_parameters;
  callbacks.handle_command_line = cismlp_handle_command_line;
  callbacks.handle_pattern_file = cismlp_handle_pattern_file;
  callbacks.handle_sequence_file = cismlp_handle_sequence_file;
  callbacks.handle_background_seq_file = cismlp_handle_background_seq_file;
  callbacks.handle_pattern_pvalue_cutoff = cismlp_handle_pattern_pvalue_cutoff;
  callbacks.handle_sequence_pvalue_cutoff = cismlp_handle_sequence_pvalue_cutoff;
  callbacks.handle_site_pvalue_cutoff = cismlp_handle_site_pvalue_cutoff;
  callbacks.handle_sequence_filtering = cismlp_handle_sequence_filtering;
  callbacks.end_parameters = cismlp_end_parameters;
  callbacks.start_multi_pattern_scan = cismlp_start_multi_pattern_scan;
  callbacks.start_pattern = cismlp_start_pattern;
  callbacks.start_scanned_sequence = cismlp_start_scanned_sequence;
  callbacks.start_matched_element = cismlp_start_matched_element;
  callbacks.handle_sequence = cismlp_handle_sequence;
  callbacks.end_matched_element = cismlp_end_matched_element;
  callbacks.end_scanned_sequence = cismlp_end_scanned_sequence;
  callbacks.end_pattern = cismlp_end_pattern;
  callbacks.end_multi_pattern_scan = cismlp_end_multi_pattern_scan;
  callbacks.end_cis_element_search = cismlp_end_cis_element_search;
  callbacks.end_cisml = cismlp_end_cisml;
  // parse the file
  if (!parse_cisml(&callbacks, &data, cisml_filename)) {
    return NULL; //FIXME cleanup
  }
  // return the result
  return data.cisml;
}
