/****************************************************************************
 * FILE: bed-io.c
 * AUTHOR: Timothy L. Bailey
 * CREATE DATE: 3/28/19
 * COPYRIGHT: 2019 UNR
 ****************************************************************************/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "bed-io.h"
#include "io.h"
#include "utils.h"
#include "string-list.h"

#define NUM_ALLOC 1000 /* Allocate how many BED lines at once? */

//
// Read a BED file into memory.  
// Coordinates are assumed to be 0-based, half-open.
// Returns the number of entries or -1 if there was a problem.
//
int read_bed_file(
  FILE *bed_file,
  bool convert_coords,		// convert to 1-based, closed
  BED_TYPE type,		// type of BED file
  BED_ENTRY_T*** bed_entries
) {
  bool error = false;

  int required_fields = 7;
  switch (type) {
    case BED3:
      required_fields = 3;
      break;
    case BROADPEAK:
      required_fields = 7;
      break;
    case NARROWPEAK:
      required_fields = 7;
      break;
    default:
      die("Unknown BED file type."); 
  }

  // Allocate space for the BED entries.
  int num_allocated = NUM_ALLOC;
  *bed_entries = (BED_ENTRY_T**)mm_malloc(sizeof(BED_ENTRY_T*) * num_allocated);

  // Read in the BED file.
  char *line = NULL;
  int n_entries = 0;
  STRING_LIST_T *fields;
  while ((line = getline2(bed_file)) != NULL) {

    // Skip lines with an initial '#'. These are comments.
    if (line[0] == '#') {
      continue;
    }

    int len = strlen(line);
    if (line[len-1] == '\n') line[strlen(line)-1] = '\0';	// chomp()

    // Split the BED line fields on tab.
    fields = new_string_list_char_split('\t', line);
    int nfields = get_num_strings(fields);

    // Check that it had the right number of fields.
    if (nfields < required_fields) {
      error = true;
      fprintf(stderr, "BED file entry %d has only %d fields rather than the required %d.\n",
        n_entries+1, nfields, required_fields);
      fprintf(stderr, "Entry: '%s'\n", line);
      break;
    }

    // Create the BED entry.
    BED_ENTRY_T *bed_entry = mm_malloc(sizeof(BED_ENTRY_T));
    bed_entry->chrom = strdup(get_nth_string(0, fields));
    bed_entry->chrom_start = atol(get_nth_string(1, fields));
    bed_entry->chrom_end = atol(get_nth_string(2, fields));
    if (type != BED3) {
      bed_entry->signal_value = strtod(get_nth_string(6, fields), NULL);
    }

    // Check format.
    if (bed_entry->chrom_start < 0) {
      error = true;
      fprintf(stderr, "BED file entry %d is not in 0-based format: start %ld < 0\n",
        n_entries+1, bed_entry->chrom_start);
      fprintf(stderr, "Entry: '%s'\n", line);
      break;
    } else if (bed_entry->chrom_start >= bed_entry->chrom_end) {
      error = true;
      fprintf(stderr, "BED file entry %d is not in half-open format: start %ld >= end %ld\n",
        n_entries+1, bed_entry->chrom_start, bed_entry->chrom_end);
      fprintf(stderr, "Entry: '%s'\n", line);
      break;
    }

    // Convert to 1-based, closed format?
    if (convert_coords) bed_entry->chrom_start++;

    // free storage
    myfree(line);
    free_string_list(fields);

    // Save the BED entry, allocating more space if necessary.
    if (n_entries == num_allocated) {
      num_allocated += NUM_ALLOC;
      *bed_entries = (BED_ENTRY_T**)mm_realloc(*bed_entries, sizeof(BED_ENTRY_T*) * num_allocated);
    }
    (*bed_entries)[n_entries++] = bed_entry;
  } // read BED file

  // Return.
  if (error) { 
    myfree(line);
    free_string_list(fields);
    // Free space used by BED entries since there was an error;
    int i;
    for (i=0; i<n_entries; i++) free_bed_entry((*bed_entries)[i]);
    if (*bed_entries) myfree(*bed_entries);
    *bed_entries = NULL; 
    return -1;			// error occurred
  } else {
    return n_entries;		// success
  }
} // read_bed_file
