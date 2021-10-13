/****************************************************************************
 * FILE: gtf-io.c
 * AUTHOR: Timothy L. Bailey
 * CREATE DATE: 3/28/19
 * COPYRIGHT: 2019 UNR
 ****************************************************************************/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "gtf-io.h"
#include "io.h"
#include "utils.h"
#include "string-list.h"

#define NUM_ALLOC 1000 /* Allocate how many GTF lines at once? */

//
// Parse an ID or ID list attribute field.
// Strips off any quotation marks (and any extension).
// Sets the value if successful to the (first) ID.
// Returns true if successful, false otherwise.
// 
//	Field format:
// 	<name> ["]<ID>[.<extension>]["]
// 	<name> ["]<ID>.<extension>["]
//
static inline bool parse_id_attribute(
  char *field, 
  char *name, 
  int name_len,
  char **value
) {
  char *start=NULL, *end=NULL;

  if ((start=strstr(field, name)) == NULL) return false;

  start += name_len;
  for ( ;*start != '\0' && (*start == ' ' || *start == '"'); start++);
  if (*start == '\0') return false;
  //for (end=start; *end != '\0' && *end != '.' && *end != '"' && *end != ','; end++);
  for (end=start; *end != '\0' && *end != '"' && *end != ','; end++);
  int len = (int) (end - start);

  *value = strndup(start, len);
  return(true);
} // parse_id_attribute

//
// Parse a NAME attribute field.
// Strips off any quotation marks.
// Sets the value if successful.
// Returns true if successful, false otherwise.
// 
//	Field format:
// 	<name> ["]<NAME>]["]
//
static inline bool parse_name_attribute(
  char *field, 
  char *name, 
  int name_len,
  char **value
) {
  char *start=NULL, *end=NULL;

  if ((start=strstr(field, name)) == NULL) return false;

  start += name_len;
  for ( ;*start != '\0' && (*start == ' ' || *start == '"'); start++);
  if (*start == '\0') return false;
  for (end=start; *end != '\0' && *end != '"'; end++);
  int len = (int) (end - start);

  *value = strndup(start, len);
  return(true);
} // parse_name_attribute

//
// Parse an expression attribute field.
// Strips off any quotation marks.
// Sets the value if successful.
// Returns true if successful, false otherwise.
// 
//	Field format:
// 	<name>[12] ["]<value>["]
//
static inline bool parse_expr_attribute(
  char *field, 
  char *name, 
  int name_len,
  double *value
) {
  char *start=NULL, *end=NULL;

  if ((start=strstr(field, name)) == NULL) return false;

  start += name_len;
  if (*start == '1' || *start == '2') start++;		// <name>[12]
  for ( ;*start != '\0' && (*start == ' ' || *start == '"'); start++);
  if (*start == '\0') return false;

  *value = strtod(start, NULL);
  return(true);
} // parse_expr_attribute

//
// Read a GTF file into memory.  
// Coordinates are assumed to be 1-based, closed.
// Returns the number of entries or -1 if there was a problem.
//
int read_gtf_file(
  FILE *gtf_file,
  bool convert_coords,		// convert to 0-based, half-open format
  GTF_TYPE type,		// type of GTF file
  GTF_ENTRY_T*** gtf_entries
) {
  int i;
  bool error = false;

  // Set the number of matches expected for different attributes;
  int expected_gtf_fields = 9;		// GTF records have 9 fields
  int expected_gid_matches = 1;		// gene_id must be present
  int expected_tid_matches = 1;		// transcript_id must be present
  int expected_gnm_matches = 1;
  int expected_tty_matches = 1;
  int expected_expr_matches = 2;	// maximum number
  char *gene_id_tok = "gene_id";
  char *gene_name_tok = "gene_name";
  char *transcript_id_tok = "transcript_id";
  char *transcript_type_tok = "transcript_type";
  char *transcript_biotype_tok = "transcript_biotype";	// alternate for transcript_type
  char *expression_tok = "rpm";
  switch (type) {
    case GENCODE:
      expected_gnm_matches = -1;	// optional
      expected_expr_matches = 0;
      break;
    case CAGE:
      expected_gnm_matches = 0;
      expected_tty_matches = 0;
      transcript_id_tok = "trlist";
      break;
    case LONGPAP:
      expected_gnm_matches = 0;
      expected_tty_matches = 0;
      expression_tok = "PKM";
      break;
    default:
      die("Unknown GTF file type."); 
  }
  int gene_id_len = strlen(gene_id_tok);
  int transcript_id_len = strlen(transcript_id_tok);
  int gene_name_len = strlen(gene_name_tok);
  int transcript_type_len = strlen(transcript_type_tok);
  int transcript_biotype_len = strlen(transcript_biotype_tok);
  int expression_len = strlen(expression_tok);

  // Allocate space for the GTF entries.
  int num_allocated = NUM_ALLOC;
  *gtf_entries = (GTF_ENTRY_T**)mm_malloc(sizeof(GTF_ENTRY_T*) * num_allocated);

  // Read in the GTF file.
  char *line = NULL;
  int n_entries = 0;
  STRING_LIST_T *fields=NULL, *attribute_fields=NULL;
  while ((line = getline2(gtf_file)) != NULL) {
    int len = strlen(line);
    if (line[len-1] == '\n') line[strlen(line)-1] = '\0';	// chomp()

    // Split the GTF line fields on tab.
    fields = new_string_list_char_split('\t', line);
    int nfields = get_num_strings(fields);

    // Check that it had the right number of fields.
    if (nfields != expected_gtf_fields) {
      error = true;
      fprintf(stderr, "GTF file entry %d has %d fields rather than the required %d.\n", 
        n_entries+1, nfields, expected_gtf_fields);
      fprintf(stderr, "Entry: '%s'\n", line);
      break;
    }

    // Create the GTF entry and fill in from the non-attribute fields.
    GTF_ENTRY_T *gtf_entry = mm_malloc(sizeof(GTF_ENTRY_T));
    gtf_entry->seqname = strdup(get_nth_string(0, fields));
    gtf_entry->start = atol(get_nth_string(3, fields));
    gtf_entry->end = atol(get_nth_string(4, fields));
    gtf_entry->strand = get_nth_string(6, fields)[0];
    gtf_entry->gene_id = NULL;
    gtf_entry->gene_name = NULL;
    gtf_entry->transcript_id = NULL;
    gtf_entry->transcript_type = NULL;
    gtf_entry->expr = -1;

    // Check format.
    if (gtf_entry->start < 1) {
      error = true;
      fprintf(stderr, "GTF file entry %d is not in 1-based format: start %ld < 1\n",
	n_entries+1, gtf_entry->start);
      fprintf(stderr, "Entry: '%s'\n", line);
      break;
    } else if (gtf_entry->start > gtf_entry->end) {
      error = true;
      fprintf(stderr, "GTF file entry %d is not in closed format: start %ld > end %ld\n",
	n_entries+1, gtf_entry->start, gtf_entry->end);
      fprintf(stderr, "Entry: '%s'\n", line);
      break;
    }

    // Convert to 0-based, half-open format/
    if (convert_coords) gtf_entry->start--;

    // Split the attributes on semicolon
    char *attributes = get_nth_string(8, fields);
    attribute_fields = new_string_list_char_split(';', attributes);
    int nattr = get_num_strings(attribute_fields);

    // Parse the attributes.
    bool attr_error = false;
    int n_expr_found = 0;
    double expression[2];
    for (i=0; i<nattr; i++) {
      char *field = get_nth_string(i, attribute_fields);
      // Get the gene ID.
      if (expected_gid_matches != 0 && gtf_entry->gene_id == NULL) {
	if (parse_id_attribute(field, gene_id_tok, gene_id_len, &(gtf_entry->gene_id))) {
	  continue;
     	}
      }
      // Get the transcript ID.
      if (expected_tid_matches != 0 && gtf_entry->transcript_id == NULL) {
	if (parse_id_attribute(field, transcript_id_tok, transcript_id_len, &(gtf_entry->transcript_id))) {
	  continue;
     	}
      }
      // Get the gene name.
      if (expected_gnm_matches != 0 && gtf_entry->gene_name == NULL) {
	if (parse_name_attribute(field, gene_name_tok, gene_name_len, &(gtf_entry->gene_name))) {
	  continue;
        }
      }
      // Get the transcript type.  Allow alternate "transcript_biotype" for Ensembl format GTFs.
      if (expected_tty_matches != 0 && gtf_entry->transcript_type == NULL) {
	if (parse_name_attribute(field, transcript_type_tok, transcript_type_len, &(gtf_entry->transcript_type))) {
	  continue;
     	} else if (parse_name_attribute(field, transcript_biotype_tok, transcript_biotype_len, &(gtf_entry->transcript_type))) {
	  continue;
     	} 
      }
      // Get the expression.
      if (n_expr_found < expected_expr_matches) {
	if (parse_expr_attribute(field, expression_tok, expression_len, &expression[n_expr_found])) {
          n_expr_found++;
	  continue;
     	}
      }
    } // attributes

    // Average the (up to two) expression values.
    if (n_expr_found == 1) {
      gtf_entry->expr = expression[0];
    } else if (n_expr_found == 2) {
      gtf_entry->expr = (expression[0] + expression[1]) / 2.0;
    }

    //
    // Check that required attributes were set.
    //
    if (expected_gid_matches > 0 && gtf_entry->gene_id == NULL) {
      attr_error = error = true;
      fprintf(stderr, "Valid GTF gene ID attribute ('%s') was not found in line %d.\n", gene_id_tok, n_entries+1);
    }
    if (expected_tid_matches > 0 && gtf_entry->transcript_id == NULL) {
      attr_error = error = true;
      fprintf(stderr, "Valid GTF transcript ID attribute ('%s') was not found in line %d.\n", transcript_id_tok, n_entries+1);
    }
    if (expected_gnm_matches > 0 && gtf_entry->gene_name == NULL) {
      attr_error = error = true;
      fprintf(stderr, "Valid GTF gene name attribute ('%s') was not found in line %d.\n", gene_name_tok, n_entries+1);
    }
    if (expected_tty_matches > 0 && gtf_entry->transcript_type == NULL) {
      attr_error = error = true;
      fprintf(stderr, "Valid GTF transcript type attribute ('%s') was not found in line %d.\n", transcript_type_tok, n_entries+1);
    }
    if (expected_expr_matches > 0 && gtf_entry->expr == -1) {
      attr_error = error = true;
      fprintf(stderr, "Valid GTF expression attribute ('%s') was not found in line %d.\n", expression_tok, n_entries+1);
    }

    // Quit if there was a problem with the attributes.
    if (attr_error) {
      fprintf(stderr, "Attributes: '%s'\n", attributes);
      break;
    }

    // Save the GTF entry, allocating more space if necessary.
    if (n_entries == num_allocated) {
      num_allocated += NUM_ALLOC;
      *gtf_entries = (GTF_ENTRY_T**)mm_realloc(*gtf_entries, sizeof(GTF_ENTRY_T*) * num_allocated);
    }
    (*gtf_entries)[n_entries++] = gtf_entry;

    // free storage
    myfree(line);
    free_string_list(fields);
    free_string_list(attribute_fields);

  } // read GTF file

  // Return.
  if (error) { 
    myfree(line);
    free_string_list(fields);
    free_string_list(attribute_fields);
    // Free space used by GTF entries since there was an error;
    int i;
    for (i=0; i<n_entries; i++) free_gtf_entry((*gtf_entries)[i]);
    if (*gtf_entries) myfree(*gtf_entries);
    *gtf_entries = NULL; 
    return -1;			// error occurred
  } else {
    return n_entries;		// success
  }

} // read_gtf_file
