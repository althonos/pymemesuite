/********************************************************************
 * MOMO Portal
 ********************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "display.h"
#include "dir.h"
#include "fasta-io.h"
#include "momo.h"
#include "array-list.h"
#include "io.h"
#include "linked-list.h"
#include "motif-in.h"
#include "motif-spec.h"
#include "seq-reader-from-fasta.h"
#include "read_seq_file.h"
#include "momo-input.h"
#include "prealigned-io.h"
#include "momo-algorithm.h"
#include "momo-modl.h"


/* local variables */
#define DATA_HASH_SIZE 100003
#define NUM_ALLOC 1000 /* Allocate how many sequences at once? */

/***********************************************************************
 Initialize mod info
 ***********************************************************************/
MOD_INFO_T* init_modinfo(MOMO_OPTIONS_T* options, SUMMARY_T* summary) {
  int i;
  const char* alph_letters = summary->alph_letters;
  
  MOD_INFO_T * modinfo = mm_malloc(sizeof(MOD_INFO_T));
  
  // Create and init an entry inside the mod table
  modinfo->mod_occurrences = 1;
  modinfo->seq_table = (options->eliminate_repeat_width) ? hash_create(DATA_HASH_SIZE, NULL) : NULL;
  modinfo->seq_list = arraylst_create();
  modinfo->motifinfos = arraylst_create();
  modinfo->amino_acids = (bool*) mm_malloc(sizeof(bool) * strlen(alph_letters));
  for (i = 0; i < strlen(alph_letters); ++i) {
    modinfo->amino_acids[i] = false;
  }
  modinfo->mod_name = NULL;
  modinfo->shuf_seq_table = NULL;
  modinfo->shuf_seq_list = NULL;
  modinfo->bg_seq_table = NULL;
  modinfo->bg_seq_list = NULL;
  
  modinfo->modl_ops = NULL;
  
  return modinfo;
}

/***********************************************************************
 Clean up mod info
 ***********************************************************************/
void cleanup_modinfo(
  MOMO_OPTIONS_T* options,
  MOD_INFO_T* modinfo) {
  int i;
  // Each MOD_INFO_T contains: motifinfos, mod_name, mod_occurrences, seq_list/seq_table, bg_sequences, amino_acids
  
  // Free motifinfos
  // Each MOTIF_INFO_T contains: motif, seqs (these point to the inside the seq_list/seq_table)
  for (i = 0; i < arraylst_size(modinfo->motifinfos); ++i) {
    MOTIF_INFO_T* motifinfo = arraylst_get(i, modinfo->motifinfos);
    destroy_motif(motifinfo->motif);
    arraylst_destroy(NULL, motifinfo->seqs);
    myfree(motifinfo);
  }
  arraylst_destroy(NULL, modinfo->motifinfos);
  
  // Free sequences from seq_list/seq_table
  for (i = 0; i < arraylst_size(modinfo->seq_list); ++i) {
    free_seq(modinfo->seq_table == NULL ? arraylst_get(i, modinfo->seq_list) : hash_get_entry_value(arraylst_get(i, modinfo->seq_list)));
  }

  arraylst_destroy(NULL, modinfo->seq_list);
  if (modinfo->seq_table != NULL) {
    hash_destroy(modinfo->seq_table);
  }
  
  // clean up bg sequences
  if (options->db_background && modinfo->bg_seq_list != NULL) {
    for (i=0; i < arraylst_size(modinfo->bg_seq_list); ++i) {
      free_seq(modinfo->bg_seq_table == NULL ? arraylst_get(i, modinfo->bg_seq_list) : hash_get_entry_value(arraylst_get(i, modinfo->bg_seq_list)));
    }
    arraylst_destroy(NULL, modinfo->bg_seq_list);
    if (modinfo->bg_seq_table != NULL) {
      hash_destroy(modinfo->bg_seq_table);
    }
  }
  myfree(modinfo->amino_acids);
  
  // clean up modl ops
  if (modinfo->modl_ops != NULL) {
    for (i = 0; i < arraylst_size(modinfo->modl_ops); ++i) {
      MODL_STEP_T* modl_step = arraylst_get(i, modinfo->modl_ops);
      cleanup_modl_step(modl_step);
    }
    arraylst_destroy(NULL, modinfo->modl_ops);
  }
  myfree(modinfo);
}

/**
 * If a protein database is provided, parse the sequences
 */
void read_protein_database_sequences(MOMO_OPTIONS_T* options,
                                     SUMMARY_T* summary,
                                     SEQ_T*** bg_sequences,
                                     int* num_bg_sequences) {
  char* protein_database_filename = options->protein_database_filename;
  ALPH_T* alph = summary->alph;
  
  FILE *fp = fopen(protein_database_filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot open the protein database file `%s'.\n", protein_database_filename);
    exit(1);
  }
  // Determine if file is FASTA format or not.
  int c;
  int previous_c = '\n';
  bool is_fasta = false;
  while ((c = getc(fp)) != EOF) {
    if (! (c == ' ' || c == '\t' || c == '\n')) {	// non-whitespace found
      if (c == '>' && previous_c == '\n') {		// FASTA ID line found
	is_fasta = true;
      } else if (c == '>') {
	die("Unknown file type: the '>' of FASTA ID lines must be the first character in a line.\n");
      }
      ungetc(c, fp);					// push back first non-whitespace
      break;
    }
    previous_c = c;
  }
  if (is_fasta) {
    read_many_fastas(alph, fp, MAX_SEQ, num_bg_sequences, bg_sequences);
  } else { // else Raw 
    read_many_line_delimited_sequences(alph, fp, MAX_SEQ, num_bg_sequences, options->width, bg_sequences);
  }
  fclose(fp);
  options->bg_filetype = is_fasta ? Fasta : Raw;
  if (! is_fasta) options->hash_fasta = false;
}

/**
 * Given a sequence and a start idx, we will change result to contain a peptide
 * of the length k. If the idx is is negative or greater than the length
 * of the sequence, we will replace with X.
 */
void set_kmer(char** result, char* sequence, int start_idx, int k) {
  int i;
  (*result)[k] = '\0';
  for (i = 0; i < k; ++i) {
    int curr_idx = start_idx + i; // start_idx + i
    if (curr_idx < 0 || curr_idx > strlen(sequence) - 1) {
      (*result)[i] = 'X';
    } else {
      (*result)[i] = sequence[curr_idx];
    }
  }
}

/**
 * Sets the O(1) lookup table. This hashes from every unique kmer within
 * the protein database to an arraylist of KMER_LOC_T objects.
 */
void create_hash_fasta_preprocess_table(MOMO_OPTIONS_T* options,
                                        SUMMARY_T* summary,
                                        SEQ_T** bg_sequences,
                                        int num_bg_sequences) {
  int i;
  int j;
  int k;
  
  if (options->hash_fasta && options->protein_database_filename) {
    HASH_TABLE hash_fasta_table = hash_create(DATA_HASH_SIZE, (void (*)(void*))arraylst_free);	// hash table of kmers
    int hash_fasta_width = options->hash_fasta_width;
    char* kmer = mm_malloc(hash_fasta_width + 1);
    double n_found = 0;
    double n_test = 0;
    for (i = 0; i < num_bg_sequences; i++) {
      char* raw_sequence = get_raw_sequence(bg_sequences[i]);
      int raw_sequence_length = strlen(raw_sequence);
      for (j = 0; j < raw_sequence_length - hash_fasta_width + 1; ++j) {
        n_test++;
        set_kmer(&kmer, raw_sequence, j, hash_fasta_width);
        KMER_LOC_T *kmerloc = mm_malloc(sizeof(KMER_LOC_T));
        kmerloc->seqidx = i;
        kmerloc->kmeridx = j;
        HASH_TABLE_ENTRY *hash_entry;
        char *dptr = 0;
        if (!(hash_entry = hash_lookup_str(kmer, hash_fasta_table))) {
          ARRAYLST_T *locations = arraylst_create();
          arraylst_add(kmerloc, locations);
          hash_insert_str_value(kmer, locations, hash_fasta_table);
        } else {
          n_found++;
          ARRAYLST_T *locations = hash_get_entry_value(hash_entry);
          arraylst_add(kmerloc, locations);
        }
      }
    }
    summary->hash_fasta_table = hash_fasta_table;
    
    // Clean up
    myfree(kmer);
  }
}

/**
 * Returns whether the given value in the filter field passes
 * the filter threshold specified by the user
 */
bool check_filter_threshold(char* value_str,
                                 MOMO_OPTIONS_T* options) {
  if (!options->filter) {
    return true;
  } else {
    FILTERTYPE_T filter_type = options->filter_type;
    double filter_threshold = options->filter_threshold;
    double value = atof(value_str);
    if (filter_type == Le) {
      return (value <= filter_threshold);
    } else if (filter_type == Lt) {
      return (value < filter_threshold);
    } else if (filter_type == Eq) {
      return (value == filter_threshold);
    } else if (filter_type == Gt) {
      return (value > filter_threshold);
    } else { // filter_type = ge
      return (value >= filter_threshold);
    }
  }
}

/**
 * Given a peptide (modless), we want to see if we can find the protein
 * it originated from. If not, result is not changed.
 */
void find_fasta( 
   SUMMARY_T* summary,
   MOMO_OPTIONS_T* options,
   SEQ_T **bg_sequences,
   int num_bg_sequences,
   char* modless,
   KMER_LOC_T* result 
) { 
  int i;
  int hash_fasta_width = options->hash_fasta_width;
  HASH_TABLE hash_fasta_table = summary->hash_fasta_table;
  
  // Use hash search if possible, otherwise linear search to find
  // the peptide in the protein database.
  if (options->hash_fasta && options->protein_database_filename && strlen(modless) >= hash_fasta_width) {
    // Use the O(1) lookup table to find the location of the peptide within the protein database
    // First, see if the prefix of the string is in the hash table.
    // Temporarily truncate the string at the hash width.
    char c_save = modless[hash_fasta_width];
    modless[hash_fasta_width] = '\0';
    HASH_TABLE_ENTRY *hash_entry = hash_lookup_str(modless, hash_fasta_table);
    modless[hash_fasta_width] = c_save;	// restore string
    if (hash_entry) {
      ARRAYLST_T *kmerloclist = (ARRAYLST_T*) hash_get_entry_value(hash_entry);
      // Then, see if the whole string matches one of the places beginning with its prefix.
      for (i = 0; i < arraylst_size(kmerloclist); ++i) {
	KMER_LOC_T *kmerloc = arraylst_get(i, kmerloclist);
	char* protein = get_raw_sequence(bg_sequences[kmerloc->seqidx]);
	char* peptide = protein + kmerloc->kmeridx;
	if (strncmp(modless, peptide, strlen(modless)) == 0) {
	  result->seqidx = kmerloc->seqidx;
            result->kmeridx = kmerloc->kmeridx;
            break;		// protein found
        }
      }
    }
  } else {
    // Use linear search to find the location of a peptide within
    // the protein database. Returns a kmer info object if found.
    for (i = 0; i < num_bg_sequences; i++) {
      char* raw_sequence = get_raw_sequence(bg_sequences[i]);
      char* found = strstr(raw_sequence, modless);
      if (found) {
	result->seqidx = i;
	result->kmeridx = found - raw_sequence;
	break;			// protein found
      }
    }
  }
}

/**
 * Adds a kmer to a sequence list (and table if eliminate repeats) can refer to either bg or fg sequences
 */
void add_kmer_to_sequence_list(HASH_TABLE* seq_table, ARRAYLST_T** seq_list, char* sequence, int idx, MOMO_OPTIONS_T* options
  ) {
  // Set up variables
  int motif_width = options->width;
  int eliminate_repeat_width = options->eliminate_repeat_width;
  
  // Create eliminate_repeat_str and motif_str. This is the mod
  // along with its flanking amino acids for the respective length.
  char* eliminate_repeat_str = mm_malloc(eliminate_repeat_width + 1);
  set_kmer(&eliminate_repeat_str, sequence, idx - eliminate_repeat_width/2, eliminate_repeat_width);
  char* motif_str = mm_malloc(motif_width + 1);
  set_kmer(&motif_str, sequence, idx - motif_width/2, motif_width);
  
  bool contains_unknowns = strchr(motif_str, 'X');
  
  // Store the sequence.
  if (!(options->remove_unknowns && contains_unknowns)) {
    SEQ_T* newsequence = allocate_seq("tempname", "tempdescription", 0, motif_str);
    if (eliminate_repeat_width) {
      if (!hash_lookup_str(eliminate_repeat_str, *seq_table)) {
        hash_insert_str_value(eliminate_repeat_str, newsequence, *seq_table);
        HASH_TABLE_ENTRY* hash_entry = hash_lookup_str(eliminate_repeat_str, *seq_table);
        arraylst_add(hash_entry, *seq_list);
      } else {
        free_seq(newsequence);
      }
    } else {
      arraylst_add(newsequence, *seq_list);
    }
  }
  
  // Clean Up
  myfree(eliminate_repeat_str);
  myfree(motif_str);
}

/**
 * Add a single mod with its flanking amino acids to the mod table
 */
void add_mod_to_mod_list(HASH_TABLE mod_table,
                    char* modless,
                    int mod_index,
                    char* mod_name,
                    char* protein_id,
                    char* scan,
                    SEQ_T** bg_sequences,
                    int num_bg_sequences,
                    MOMO_OPTIONS_T* options,
                    SUMMARY_T* summary,
                    bool passes_threshold) {
  int i;
  
  // Set up variables
  MOD_INFO_T * modinfo = hash_get_entry_value(hash_lookup_str(mod_name, mod_table));
  int motif_width = options->width;
  
  // For this mod, we will create a kmerinfo object.
  KMER_LOC_T* kmerloc = mm_malloc(sizeof(kmerloc));
  kmerloc->seqidx = -1;
  kmerloc->kmeridx = -1;
  find_fasta(summary, options, bg_sequences, num_bg_sequences, modless, kmerloc);
  
  char* protein = (kmerloc->seqidx >= 0) ? get_raw_sequence(bg_sequences[kmerloc->seqidx]) : modless;
  int start_of_mod = (kmerloc->kmeridx >= 0) ? kmerloc->kmeridx + mod_index : mod_index;
  
  if (passes_threshold) {
    add_kmer_to_sequence_list(&modinfo->seq_table, &modinfo->seq_list, protein, start_of_mod, options);
  }
  free(kmerloc);
}

/*
	Get the contents of the next non-empty line.
	(Empty means nothing but whitespace.)
	Return number of characters read, or -1 if EOF or error.
*/
static int getline_nonempty(
  char **linep,			// IN/OUT; pass pointer to NULL on first call;
				// used by getline; do not free
  size_t *linecapp,			// IN/OUT; used by getline
  FILE *stream			// IN; the open file pointer
) {
  int i;
  ssize_t n_read = 0;
  bool found = false;
  while (!found && (n_read = getline(linep, linecapp, stream)) != -1) {
    for (i=0; !found && i<n_read; i++) if (! isspace((*linep)[i])) found = true;
  }
  return(n_read);
} // getline_nonempty

/**
 * Get file type: PSM, FASTA or Raw 
 * Looks for 
 *   1) FASTA if first character of first non-empty line is ">"
 *   2) PSM if first non-empty line has >1 tab-delimited fields
 *   3) Raw, otherwise
*/
FILETYPE_T get_filetype(
  FILE *fp
) {
  FILETYPE_T filetype;
  // Find the first non-empty line.
  char *line = NULL;
  size_t len = 0;
  if (getline_nonempty(&line, &len, fp) != -1) {
    // Determine the file type.
    if (line[0] == '>') {
      filetype = Fasta;
    } else if (strstr(line, "\t") != NULL) {
      filetype = Psm;
    } else {
      filetype = Raw;
    }
  } else {
    filetype = Raw;
  }

  // Free the line string.
  //if (line) myfree(line);  // Acts like it wasn't malloc'ed.

  // Rewind the input.
  rewind(fp);
  
  return(filetype);   
} // get_filetype

/**
 * Parse a PSM header line. We want to locat the idx of filter_value, sequence, scan, and protein_id
 */
void parse_tsv_header(MOMO_OPTIONS_T* options,
                      int* index_of_filter,
                      int* index_of_sequence,
                      int* index_of_scan,
                      int* index_of_protein,
                      FILE* fp, char** header,
                      size_t* len) {
  int i = 0;
  if (getline_nonempty(header, len, fp) != -1) {
    *header = strsep(header, "\n");	// strip new line
    char *token;
    while((token = strsep(header, "\t")) != NULL) {
      // Trim leading & trailing whitespace from token.
      while (isspace(*token)) token++;
      int t_end = strlen(token)-1;
      while (t_end >= 0 && isspace(token[t_end])) t_end--;
      char t_save = token[t_end+1];
      token[t_end+1] = '\0';
      if (strcmp(token, options->sequence_column) == 0) {
	*index_of_sequence = i;
      }
      if (options->filter && strcmp(token, options->filter_field) == 0) {
	*index_of_filter = i;
      }
      if (strcmp(token, "scan") == 0) {
	*index_of_scan = i;
      }
      if (strcmp(token, "protein id") == 0) {
	*index_of_protein = i;
      }
      token[t_end+1] = t_save;
      i++;
    }
  }
  // ensure that input file has required columns
  if (*index_of_sequence == -1) {
    die("Could not find the modified peptide column named '%s' in the PSM file.", options->sequence_column);
  }
  if (options->filter && *index_of_filter == -1) {
    die("Could not find the filter column named '%s' in the PSM file.", options->filter_field);
  }
}

/**
 * Parse a tsv body line into filter_value, sequence, scan, and protein_id
 */
void parse_tsv_bodyline(MOMO_OPTIONS_T* options,
                        char** filter_value,
                        char** sequence,
                        char** scan,
                        char** protein_id,
                        int index_of_filter,
                        int index_of_sequence,
                        int index_of_scan,
                        int index_of_protein,
                        char* line) {
  char* token;
  int i = 0;
  while ((token = strsep(&line, "\t")) != NULL) {
    if (i == index_of_sequence) {
      *sequence = strdup(token);
    }
    if (i == index_of_scan) {
      *scan = strdup(token);
    }
    if (i == index_of_protein) {
      *protein_id = strdup(token);
    }
    if (options->filter && i == index_of_filter) {
      *filter_value = strdup(token);
    }
    i++;
  }
}

/**
 * Given a sequence, duplicate it into modless sequence and fill
 * an array with each position indicating whether or not it has a mod. 
 * If yes, then will contain the mod name. Otherwise, it is set to NULL.
 * The mod is always assume to be in the center of the sequence (position len/2).
 */
void get_prealigned_mod_array(ARRAYLST_T* mod_array, char** modless,  char* sequence, const char *alph_letters) {
  *modless = strdup(sequence);
  char* mod = mm_malloc(2);
  int center = strlen(sequence)/2;
  // Make sure center is a valid character.
  if (strchr(alph_letters, sequence[center])) { 
    mod[0] = sequence[center];
    mod[1] = '\0';
    arraylst_put(center, mod, mod_array);
  }
}

/**
 * Given a modified sequence, set the modless sequence and fill
 * an array with each position indicating whether or not it has a mod.
 * If yes, then will contain the mod name. Otherwise, it is set to NULL.
 */
void get_tsv_mod_array(ARRAYLST_T* mod_array, char** modless, char* sequence, const char *alph_letters) {
  size_t seq_len = strlen(sequence);
  *modless = mm_malloc(seq_len + 1);
  *modless[0] = '\0';
  
  char* sequence_mod_copy = strdup(sequence);
  char* sequence_mod_iterator = sequence_mod_copy;
  char* ptr1 = strchr(sequence_mod_iterator, '[');
  char* ptr2 = strchr(sequence_mod_iterator, ']');  
  
  // while mod has been found
  while (ptr1 && ptr2) {
    // find the length of the mod. For example, S[79.0] = 5 since we disregard the brackets.
    int mod_len = ptr2-ptr1;

    // make sure mod is valid letter
    if (ptr1 > sequence_mod_iterator && strchr(alph_letters, *(ptr1-1))) {
      // create a copy of the mod
      char* mod = mm_malloc(mod_len + 1);
      strncpy(mod, ptr1-1, 1);
      strncpy(mod+1, ptr1+1, mod_len - 1);
      mod[mod_len] = '\0';
      
      // update the modless sequence
      strncat(*modless, sequence_mod_iterator, ptr1-sequence_mod_iterator);
      
      // store the mod at its proper index.
      arraylst_put(strlen(*modless) - 1, mod, mod_array);
      
      //update
      sequence_mod_iterator = ptr2+1;
      ptr1 = strchr(sequence_mod_iterator, '[');
      ptr2 = strchr(sequence_mod_iterator, ']');
    }
  }
  // add the remaining amino acids to the end of modless
  strncat(*modless, sequence_mod_iterator, strlen(sequence_mod_iterator));
  
  myfree(sequence_mod_copy);
}

/**
 * Given a single phosphorylation file, add all the modifications and respective flanking
 * amino acid information into the mod table.
 */
void add_phospho_file_to_table(char* phospho_filename,
                                  MOMO_OPTIONS_T* options,
                                  SUMMARY_T* summary,
                                  SEQ_T** bg_sequences,
                                  int num_bg_sequences,
                                  SEQ_T*** fg_sequences,
                                  int* num_fg_sequences) {
  
  int i = 0, j = 0;
  const char* alph_letters = summary->alph_letters;
  ALPH_T* alph = summary->alph;

  // open and prepare file for reading
  FILE* fp;
  char* line = NULL;
  size_t len = 0;
  ssize_t read;
  fp = fopen(phospho_filename, "r");
  if (fp == NULL) {
    die("Failed to read filename: %s", phospho_filename);
  }

  // get the type of the input file
  options->filetype = get_filetype(fp);

  // parse header to find the indices of required data
  int index_of_filter = -1;           // the tab-delimited index of filter;
  int index_of_sequence = -1;         // the tab-delimited index of column sequence
  int index_of_scan = -1;             // the tab-delimited index of scan;
  int index_of_protein = -1;          // the tab-delimited index of protein id;

  // check for sequence_column if it is a PSM file, and parse the header.
  if (options->filetype == Psm) {
    if (options->sequence_column == NULL) {
      fprintf(stderr, "File \'%s\' is in PSM format so you must use option --psm-type or --sequence-column.\n", 
	phospho_filename);
      fprintf(stderr, "%s", options->usage);
      exit(EXIT_FAILURE);
    }
    parse_tsv_header(options, &index_of_filter, &index_of_sequence, &index_of_scan, &index_of_protein, fp, &line, &len);
  }
  
  HASH_TABLE mod_table = summary->mod_table;
  ARRAYLST_T * mod_table_keys = summary->mod_table_keys;
  
  // for each body line (PSM) inside the file:
  bool found_mod = false;
  while ((read = getline_nonempty(&line, &len, fp)) != -1) {
    // strip new line and make a copy
    line = strsep(&line, "\n");
    
    char* unprocessed_sequence = NULL;
    char* sequence = NULL;
    char* filter_value = NULL;
    
    // These are currently not used.
    char* scan = NULL;
    char* protein_id = NULL;
    
    if (options->filetype == Psm) {
      parse_tsv_bodyline(options, &filter_value, &unprocessed_sequence, &scan, &protein_id, index_of_filter, index_of_sequence, index_of_scan, index_of_protein, line);
    } else { // Raw/Fasta
      if (line == NULL || line[0] == '>') continue; 	// Skip FASTA ID lines.
      parse_prealigned_bodyline(&unprocessed_sequence, line, options->width);
    }
    
    if (options->filetype != Psm) {
      sequence = strdup(unprocessed_sequence);
    } else {
      // PSMs may appear in different formats, we will modify them to be in standard crux
      if (strlen(unprocessed_sequence) > 2 && unprocessed_sequence[1] == '.') {
        // sequences are in A.AS[79.9]A.- format, remove the dots and '-' from ends.
        sequence = mm_malloc(strlen(unprocessed_sequence) + 1);
        sequence[0] = '\0';
        if (unprocessed_sequence[0] != '-') {
          strncat(sequence, unprocessed_sequence, 1);
        }
        int pep_length = strlen(unprocessed_sequence) - 4;
        strncat(sequence, unprocessed_sequence + 2, pep_length);
        if (unprocessed_sequence[strlen(unprocessed_sequence) - 1] != '-') {
          strncat(sequence, unprocessed_sequence + strlen(unprocessed_sequence) - 1, 1);
        }
      } else {
        int count = 0;
        int plusminus_idx = -1;
        for (i = 0; i < strlen(unprocessed_sequence); ++i) {
          if (unprocessed_sequence[i] == '+' || unprocessed_sequence[i] == '-') {
            count++;
            if (plusminus_idx == -1) {
              plusminus_idx = i;
            }
          }
        }
        if (plusminus_idx != -1 && unprocessed_sequence[plusminus_idx-1] != '[') {
          // sequences are in AS+79.9A format, sequences are missing brackets, create bracket rep
          sequence = mm_malloc(strlen(unprocessed_sequence) + 2*count + 1);
          sequence[0] = '\0';
          for (i = 0; i < strlen(unprocessed_sequence); ++i) {
            if (unprocessed_sequence[i] == '+' || unprocessed_sequence[i] == '-') {
              strncat(sequence, "[", 1);
            }
            strncat(sequence, unprocessed_sequence + i, 1);
            if (!isalpha(unprocessed_sequence[i]) && isalpha(unprocessed_sequence[i+1])) {
              strncat(sequence, "]", 1);
            }
          }
        } else {
          // sequence is fine as it is
          sequence = strdup(unprocessed_sequence);
        }
      }
    }
    
    // check whether or not this PSM passes the filter threshold.
    bool passes_threshold = options->filetype == Raw || check_filter_threshold(filter_value, options);
    
    // NOTE: mod_name = X[mod_mass] when we create a mod per aminoacid + mod pair,
    // and simply [mod_mass] when we create a single mod per mass.
    
    // we will now create a modless version of sequence, while:
    // 1. update mod_table / mod_table_keys. This is a hash table
    // mapping from mod_name => MOD_INFO_T.
    // 2. create an arraylist: This stores the mod_name of mods that pass the threshold
    // at the index they are located inside the modless string.
    
    size_t seq_len = strlen(sequence);
    // This is supposed to store the locations of the mods within the sequence.
    ARRAYLST_T * mod_array = arraylst_create_sized(seq_len);
    for (i = 0; i < seq_len; i++) {
      arraylst_add(NULL, mod_array);
    }
    // modless sequence
    char* modless = NULL;
    
    // note that after we obtain the mods, we need to free them.
    int i_start, i_end;
    if (options->filetype == Psm) {
      // sequence contains multiple mods
      get_tsv_mod_array(mod_array, &modless, sequence, alph_letters);
      i_start = 0;
      i_end = seq_len;
    } else { // filetype = Raw or Fasta
      // sequence only contains a single mod
      get_prealigned_mod_array(mod_array, &modless, sequence, alph_letters);
      i_start = seq_len/2;
      i_end = i_start + 1;
    }

    // Add modless as a "sequence" to output if using the foreground as background.
    if (! options->db_background) {
      if (*num_fg_sequences % NUM_ALLOC == 0) {
	*fg_sequences = (*num_fg_sequences == 0) ?
	  (SEQ_T**)mm_malloc(sizeof(SEQ_T*) * NUM_ALLOC) :
	  (SEQ_T**)mm_realloc(*fg_sequences, sizeof(SEQ_T*) * ((*num_fg_sequences)+NUM_ALLOC));
      }
      (*fg_sequences)[*num_fg_sequences] = allocate_seq(NULL, NULL, 0, modless);
      (*num_fg_sequences)++;
    }

    // update the tables for each valid mod:
    for (i = i_start; i < i_end; i++) {
      char* mod = arraylst_get(i, mod_array);
      if (mod != NULL) {
        found_mod = true;
        char* mod_name = (options->single_motif_per_mass) ? mod + 1 : mod; // e.g. 79.9 or S79.9
        HASH_TABLE_ENTRY* hash_entry = hash_lookup_str(mod_name, mod_table);
        // update the mod_table for each mod.
        if (!hash_entry) {
          MOD_INFO_T* newmodinfo = init_modinfo(options, summary);
          
          // we are adding an entry to the hash table
          hash_insert_str_value(mod_name, (void *) newmodinfo, mod_table);
          hash_entry = hash_lookup_str(mod_name, mod_table);
          arraylst_add(hash_entry, mod_table_keys);
        } else {
          // Update the entry within the mod info
          MOD_INFO_T * currmodinfo = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
          currmodinfo->mod_occurrences++;
        }
        // update the list of amino acids associated with this mod
        MOD_INFO_T * currmodinfo = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
        currmodinfo->amino_acids[strchr(alph_letters, mod[0]) - alph_letters] = true;
        
        // find the flanks of this mod and add it to the list of sequences associated with the mod name
        add_mod_to_mod_list(mod_table, modless, i, mod_name, protein_id, scan, bg_sequences, num_bg_sequences, options, summary, passes_threshold);
      }
    }
    
    // clean up
    arraylst_destroy(free, mod_array);
    myfree(modless);
    
    myfree(unprocessed_sequence);
    myfree(sequence);
    myfree(filter_value);
    
    // These are currently not used.
    myfree(scan);
    myfree(protein_id);
    
  } // getline_nonempty
  
  // clean up
  myfree(line);
  fclose(fp);

  // Check that at least one modified peptide was found
  if (!found_mod && options->filetype == Psm) {
    fprintf(stderr, 
      "\nFATAL: No modified peptides were found in the PSM format (tab-delimited) file\n"
      "  \'%s\'.\n"
      "Check that the peptides are in an allowable format (with modification weights)\n"
      "or enter them in 'raw' or 'FASTA' format.\n"
      "Refer to %s/doc/psm-format.html for allowable formats.\n",
	phospho_filename, SITE_URL);
    exit(EXIT_FAILURE);
  }
}

/**
 * For each phosphorylation file, add all the modifications and respective flanking
 * amino acid information into the mod table.
 */
void add_phospho_files_to_table(MOMO_OPTIONS_T* options,
                                        SUMMARY_T* summary,
                                        SEQ_T** bg_sequences,
                                        int num_bg_sequences,
					SEQ_T*** fg_sequences,
                                        int* num_fg_sequences
                               ) {
  int i;
  *fg_sequences = NULL;
  *num_fg_sequences = 0;

  for (i = 0; i < arraylst_size(options->phospho_filenames); ++i) {
      add_phospho_file_to_table((char*) arraylst_get(i, options->phospho_filenames), options, summary, bg_sequences, num_bg_sequences, fg_sequences, num_fg_sequences);
  }
}

/**
  Add a shuffled copy of each foreground sequence (preserving central residue).
**/
void add_shuffled_sequences_to_table(
  MOMO_OPTIONS_T *options, 
  SUMMARY_T *summary,
  bool background			// use as bg_seq if true, otherwise as "shuf_seq"
) 
{
  if (options->algorithm == Motifx || options->algorithm == Modl) {
    int i, j, k;
    ARRAYLST_T *mod_table_keys = summary->mod_table_keys;
    const char *alph_letters = summary->alph_letters;
    char *noctr = mm_malloc(sizeof(char) * options->width);
    char *shuf1 = mm_malloc(sizeof(char) * options->width);
    char *shuf2 = mm_malloc(sizeof(char) * (options->width + 1));

    // for each mod
    for (i = 0; i < arraylst_size(mod_table_keys); ++i) {
      // initialize the shuffled table and list
      MOD_INFO_T *modinfo = (MOD_INFO_T*) hash_get_entry_value(arraylst_get(i, mod_table_keys));
      ARRAYLST_T *phospho_seqs = modinfo->seq_list;
      ARRAYLST_T *shuf_seq_list = arraylst_create();
      HASH_TABLE shuf_seq_table = (options->eliminate_repeat_width) ? hash_create(DATA_HASH_SIZE, NULL) : NULL;
      if (background) {
	modinfo->bg_seq_list = shuf_seq_list;
	modinfo->bg_seq_table = shuf_seq_table;
      } else {
	modinfo->shuf_seq_list = shuf_seq_list;
	modinfo->shuf_seq_table = shuf_seq_table;
      }

      // for each sequence
      for (j = 0; j < arraylst_size(phospho_seqs); ++j) {
        SEQ_T *active_sequence = 
          (options->eliminate_repeat_width) ? hash_get_entry_value(arraylst_get(j, phospho_seqs)) : arraylst_get(j, phospho_seqs);
        char *curr_seq = get_raw_sequence(active_sequence);
        // shuffle the sequence keeping the center untouched
        k = options->width/2;
        // remove the central position
        char save_char = curr_seq[k];
        strncpy(noctr, curr_seq, k);
        strcpy(noctr+k, curr_seq+k+1);
	// setup the globals in ushuffle
	ushuffle1(noctr, options->width-1, 1);
        // create a shuffled copy of the sequence from the globals in ushuffle
        ushuffle2(shuf1);
        shuf1[options->width-1] = '\0';
        // replace the central position
        strncpy(shuf2, shuf1, k);
        shuf2[k] = save_char;
        strcpy(shuf2+k+1, shuf1+k);
	add_kmer_to_sequence_list(&shuf_seq_table, &shuf_seq_list, shuf2, k, options);
      } // seq (j)
    } // mod (i)
    free(noctr);
    free(shuf1);
    free(shuf2);
  } // valid pgm
} // add_shuffled_sequences_to_table

/**
 * Given a list of sequences, set the background sequences for all modifications. For example,
 * if a mod is associated with amino acids S,T,Y, then find all the kmers centralized around
 * these letters and add these to the background sequences.
 */
void add_background_sequences_to_table(
  MOMO_OPTIONS_T* options, 
  SUMMARY_T* summary, 
  SEQ_T** protein_database_sequences, 
  int num_bg_sequences
) 
{
  if (options->algorithm == Motifx || options->algorithm == Modl) {
    int i, j, k;
    
    ARRAYLST_T* mod_table_keys = summary->mod_table_keys;
    const char* alph_letters = summary->alph_letters;
    
    // for each mod
    for (i = 0; i < arraylst_size(mod_table_keys); ++i) {
      // initialize the background table and list
      MOD_INFO_T* modinfo = (MOD_INFO_T*) hash_get_entry_value(arraylst_get(i, mod_table_keys));
      ARRAYLST_T* bg_seq_list = arraylst_create();
      HASH_TABLE bg_seq_table = (options->eliminate_repeat_width) ? hash_create(DATA_HASH_SIZE, NULL) : NULL;
      modinfo->bg_seq_list = bg_seq_list;
      modinfo->bg_seq_table = bg_seq_table;
      
      // amino acids for this mod
      bool* amino_acids = modinfo->amino_acids;

      // for each sequence
      for (j = 0; j < num_bg_sequences; ++j) {
        char* protein = get_raw_sequence(protein_database_sequences[j]);
        
        // if FASTA, we need to search the entire protein; k is the center
        //for (k = 0; k < strlen(protein); ++k) {
        for (k = options->width/2; k < strlen(protein) - options->width/2; ++k) {
          if (options->bg_filetype == Fasta || (options->bg_filetype == Raw && k == options->width/2)) {
            char* alph_letters_substring = strchr(alph_letters, protein[k]);
            if (alph_letters_substring) {
              int aa_idx = alph_letters_substring - alph_letters;
              if (amino_acids[aa_idx]) {
                add_kmer_to_sequence_list(&bg_seq_table, &bg_seq_list, protein, k, options);
              }
            }
          }
        }
      }
    }
  }
}
