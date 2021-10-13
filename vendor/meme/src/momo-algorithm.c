/********************************************************************
 * MOMO Portal
 ********************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <regex.h>
#include <limits.h>
#include "display.h"
#include "dir.h"
#include "momo.h"
#include "momo-output.h"
#include "matrix.h"
#include "motif-in.h"
#include "motif-spec.h"
#include "seq-reader-from-fasta.h"
#include "binomial.h"
#include "meme.h"
#include "read_seq_file.h"
#include "momo-algorithm.h"
#include "momo-simple.h"
#include "momo-motifx.h"
#include "momo-modl.h"

/* local variables */
#define DATA_HASH_SIZE 100003

/**
 * The counts of matches given a MOMO motif ID (modified regular expression) and a set of sequences.
 *
**/
int get_counts_from_motifid(
  char *motifid,				// contains RE plus "-", "." and numbers
  ARRAYLST_T* seqs,				// sequences
  MOMO_OPTIONS_T* options			// momo options
) {
  #define NUM_SUBMATCHES 6
  #define ERROR_MESSAGE_SIZE 100
  #define BUFFER_SIZE 512
  char error_message[ERROR_MESSAGE_SIZE];
  int status = 0;
  static regex_t motif_regex;
  static regmatch_t matches[NUM_SUBMATCHES];
  int  i, j;

  // Create a regular expression string from the motifid.
  char *regexp_str = malloc(sizeof(char) * (strlen(motifid) + 1));
  for (i=j=0; i < strlen(motifid); i++) {
    // Check for single motif per mass (may be "__" or "_number").
    if (motifid[i] == '_' && (motifid[i+1] == '_' || (motifid[i+1]-'0' >= 0 && motifid[i+1]-'0' <= 9))) {
      // Convert '_' to '.' for (missing) central residue.
      regexp_str[j++] = '.';
    } else if (motifid[i] == '_' || (motifid[i]-'0' >= 0 && motifid[i]-'0' <= 9) || motifid[i] == '.') {
      // skip
    } else {
      // Convert 'x' to '.'.
      regexp_str[j++] = (motifid[i] == 'x') ? '.' : motifid[i];
    }
  }
  regexp_str[j] = '\0';
  //printf("motifid %s regexp_str %s\n", motifid, regexp_str);

  // Create a regular expression from the string.
  status = regcomp(
    &motif_regex,
    regexp_str,
    REG_EXTENDED
  );
  if (status != 0) {
    regerror(status, &motif_regex, error_message, ERROR_MESSAGE_SIZE);
    die(
      "Error while intitializing regular expression\n"
      "for counting occurrences in shuffled sequences: %s\n",
      error_message
    );
  }

  // Count the occurrences of the RE in the sequences.
  int size = arraylst_size(seqs);
  int n_matches = 0;
  for (i = 0; i < size; ++i) {
    char* curr_seq = get_raw_sequence((SEQ_T*) (options->eliminate_repeat_width ? hash_get_entry_value((HASH_TABLE_ENTRY*) arraylst_get(i, seqs)) : arraylst_get(i, seqs)));
    status = regexec(&motif_regex, curr_seq, NUM_SUBMATCHES, matches, 0);
    if (!status) n_matches++;
    if(status && status != REG_NOMATCH ){
      regerror(status, &motif_regex, error_message, 100);
      die("Error scanning sequence with regex: %s\n"
          "error message is %s\n",
          curr_seq,
          error_message
      );
    }
  }

  free(regexp_str);

  return(n_matches);
} // get_counts_from_motifid

/**
 * Using an array list of sequences (or hash table entries), creates and returns a count matrix.
 */
MATRIX_T* get_count_matrix(MATRIX_T* count,
                           ARRAYLST_T* sequences,
                           MOTIFX_STATUS_T** status,
                           MOMO_OPTIONS_T* options,
                           SUMMARY_T* summary) {
  int i;
  int j;
  
  bool sequences_are_hash_entries = options->eliminate_repeat_width > 0;
  const char* alph_letters = summary->alph_letters;
  
  if (!count) {
    count = allocate_matrix(options->width, strlen(alph_letters));
  }
  //init_matrix(1e-10, count);
  init_matrix(0.0, count);
  
  for (i = 0; i < arraylst_size(sequences); ++i) {
    if (!status || (*status)[i] == ACTIVE) {
      SEQ_T* seqobject = sequences_are_hash_entries ? hash_get_entry_value(arraylst_get(i, sequences)) : arraylst_get(i, sequences);
      char* raw_sequence = get_raw_sequence(seqobject);
      
      for (j = 0; j < options->width; j++) {
        if (strchr(alph_letters, raw_sequence[j])) {
          int aa_idx = strchr(alph_letters, raw_sequence[j]) - alph_letters;
          set_matrix_cell(j, aa_idx, get_matrix_cell(j, aa_idx, count) + 1.0, count);
        }
      }
    }
  }
  return count;
}

/**
 * Prints a matrix to the terminal. For debugging use.
 */
void print_matrix_to_terminal(MATRIX_T* matrix, MOMO_OPTIONS_T* options, SUMMARY_T* summary) {
  const char* alph_letters = summary->alph_letters;
  int i;
  int j;
  for (i = 0; i < options->width-1; i++) {
    for (j = 0; j < strlen(alph_letters); j++) {
      printf("%.1f\t", get_matrix_cell(i, j, matrix));
    }
    printf("\n");
  }
  printf("\n");
}


/**
 * For each mod inside the mod table that passes all the
 * filters specified by the user, create a motif. This also
 * sets the number of (passing) mods/modtypes.
 */
void create_motifs(MOMO_OPTIONS_T* options,
                          SUMMARY_T* summary
                          ) {
  
  ALPH_T* alph = summary->alph;
  const char* alph_letters = summary->alph_letters;
  int i;
  
  // Get the mod table
  HASH_TABLE mod_table = summary->mod_table;
  ARRAYLST_T * mod_table_keys = summary->mod_table_keys;
  
  // Initialize statistic counters.
  unsigned long num_mods = 0;
  unsigned long num_types = 0;
  unsigned long pass_mods = 0;
  unsigned long pass_types = 0;
  unsigned long bg_mods = 0;

  // For each mod inside the mod table
  for (i = 0; i < arraylst_size(mod_table_keys); i++) {
    // Get the mod entry
    HASH_TABLE_ENTRY * hash_entry = arraylst_get(i, mod_table_keys);
    MOD_INFO_T * mod_entry = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
    mod_entry->mod_name = hash_get_entry_key(hash_entry);
    
    // Count the number of mods and mod types and add to overall counter
    unsigned long num_occurrences = mod_entry->mod_occurrences;
    num_mods += num_occurrences;
    num_types++;
    
    // For each passing mod, we will add to passing mod and mod type counter and create a motif
    unsigned long fg_num_passing = arraylst_size(mod_entry->seq_list);
    unsigned long bg_num_passing = mod_entry->bg_seq_list ? arraylst_size(mod_entry->bg_seq_list) : ULONG_MAX;
    if (fg_num_passing >= options->min_occurrences && bg_num_passing >= options->min_occurrences) {
      // Increment counters
      pass_mods += fg_num_passing;
      pass_types++;
      bg_mods += bg_num_passing;
      
      // Create motif
      if (options->algorithm == Simple) {
        create_simple_motif(summary, options, mod_entry);
      } else if (options->algorithm == Motifx) {
        create_motifx_motifs(summary, options, mod_entry);
      } else if (options->algorithm == Modl) {
        create_modl_motifs(summary, options, mod_entry);
      }
    }
  }
  
  summary->num_mod = num_mods;
  summary->num_modtype = num_types;
  summary->num_mod_passing = pass_mods;
  summary->num_bg_mod = bg_mods;
  summary->num_modtype_passing = pass_types;
}
