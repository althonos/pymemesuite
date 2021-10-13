#ifndef MAKE_MOMO_MOTIFS_H
#define MAKE_MOMO_MOTIFS_H

#include "alphabet.h"
#include "array-list.h"
#include "projrel.h"
#include "utils.h"
#include "hash_table.h"

static const int MAX_ALPH_SIZE = 100;

// Structure for tracking momo command line parameters.
typedef enum filtertype {
  Le,                                 // <=
  Lt,                                   // <
  Eq,                                   // =
  Gt,                                   // >
  Ge                                  // >=
} FILTERTYPE_T;
static const char *filtertype_names[] = {"<=", "<", "=", ">", ">="};

// Structure for tracking momo command line algorithm.
typedef enum algorithm {
  No_algorithm,
  Unknown_algorithm,
  Simple,
  Motifx,
  Modl
} ALGORITHM_T;
static const char *algorithm_names[] ={"No algorithm", "Unknown algorithm", "simple", "motif-x", "MoDL"};

// Structure for tracking momo command line parameters.
typedef enum filetype {
  Psm,
  Fasta,
  Raw 
} FILETYPE_T;
static const char *filetype_names[] = {"PSM", "FASTA", "Raw"};

// Structure for tracking momo command line parameters.
typedef struct options {
  ALGORITHM_T algorithm;
  
  ARRAYLST_T* phospho_filenames;        // List of filenames containg PSMs.
  ARRAYLST_T* phospho_filesizes;        // List of file sizes of above.
  
  bool allow_clobber;			// Allow overwritting of files in output directory.
  bool filter;				// Filter?
  bool db_background;			// Use protein database as background (rather than shuffled foreground).
  bool hash_fasta;			// Use a O(1) lookup table on protein database to speed up search
  bool remove_unknowns;                 // Remove sequences with unknown flanks
  bool single_motif_per_mass; 		// Create a single motif per mass ?
  bool harvard;				// mimic original motif-x behavior by 
					// truncating binomial p-value at 1e-16 and subtracting count*1e-6
					// from log_p; this sorts ties by the number of matching peptides 
  bool printp;				// print normalized p-values for debuggin
  
  char* command_line;                   // Full command line
  char* filter_field;                   // Field to filter on
  char* html_path;                      // Path to MOMO HTML output file
  char* output_dirname;                 // Name of the output directory
  char* protein_database_filename;      // Name of file file containg protein fasta database.
  char* psm_type; 			// Type of PSM file
  char* sequence_column;                // Sequence column name in psm file
  char* text_path;                      // Path to MOMO motifs in meme format
  char* tsv_path;                       // Path to MOMO TSV output file
  char* paper;				// Name of paper to cite.
  char* link;				// Link of paper to cite.

  const char* HTML_FILENAME;            // Name of HTML output file.
  const char* TEMPLATE_FILENAME;        // Name of HTML template file.
  const char* TEXT_FILENAME;            // Name of motifs in MEME format output file.
  const char* TSV_FILENAME;             // Name of TSV output file.
  const char* usage;                    // Usage statment
  
  double filter_threshold;              // Threshold to Filter on
  double score_threshold;               // Motif-X Score Threshold
  
  FILETYPE_T filetype;			// Type of PSM input file
  FILETYPE_T bg_filetype;		// Type of protein database background file
  
  FILTERTYPE_T filter_type;             // Type to filter

  int seed;				// Random number seed
  int eliminate_repeat_width;           // Eliminate Repeat length
  int hash_fasta_width;                 // Kmer length used for O(1) lookup table
  int max_iterations;                   // Maximum number of iterations for modl
  int max_motifs;                       // Maximum number of motifs for modl
  int min_occurrences;                  // Number of Occurrences
  int max_no_decrease;       		// Number of iterations of no decrease in MDL before stopping algorithm
  int width;                            // Motif Width
  int verbosity;
} MOMO_OPTIONS_T;

// Structure for tracking summary of results
typedef struct summary {
  unsigned long num_mod;                          // Number of Modifications
  unsigned long num_modtype;                      // Number of Types of Modification
  unsigned long num_mod_passing;                  // Number of Modifications After Filtering
  unsigned long num_modtype_passing;              // Number of Types of Modifications After Filtering
  unsigned long num_bg_mod;			  // Number of Background Modifications
  
  HASH_TABLE mod_table;                           // Hashes mod name to a modinfo struct
  ARRAYLST_T * mod_table_keys;                    // Arraylist containing pointers to mod_table key entries
  
  HASH_TABLE hash_fasta_table;                    // O(1) lookup table for sequences in protein database
  
  ALPH_T * alph;                                  // Alphabet Used
  const char* alph_letters;
  ARRAY_T * bg_freqs;                              // Background frequencies (from protein database if provided)
} SUMMARY_T;

// Structure containing information for a given mod type
typedef struct modinfo {
  char *mod_name;                                 // Name of mod
  unsigned long mod_occurrences;                  // Number of times this mod occurs in the file
  ARRAYLST_T *seq_list;                           // If eliminate repeats, then this will be the keys
                                                  //   for seq_table. Otherwise, it is a list of sequences that pass filters
  HASH_TABLE seq_table;                           // If eliminate repeats, will hash from char* representation
                                                  //   for sequences of length eliminate_repeats that pass filters to respective
                                                  //   sequences.
  ARRAYLST_T *shuf_seq_list;			  // Shuffled foreground sequences.
  HASH_TABLE shuf_seq_table;			  // Shuffled foreground sequences.
  ARRAYLST_T *bg_seq_list;			  // From background file.
  HASH_TABLE bg_seq_table;			  // From background file.
  ARRAYLST_T *motifinfos;                         // the list of motifs for this mod.
  bool *amino_acids;                              // a boolean array for each amino acid. true if this mod is associated with the
                                                  // amino acid. false otherwise.
  ARRAYLST_T* modl_ops;                           // Only used for MoDL algorithm. Contains a list of MODL_STEP_T*.
} MOD_INFO_T;

typedef struct motifinfo {
  MOTIF_T* motif;                                 // Motif
  ARRAYLST_T* seqs;                               // Sequences that make up this motif
  double score;					  // motif score
  int n_tests;					  // number of tests made during search for motif
  int fg_matches;                                 // matches in remaining foreground sequences (harvard-style)
  int fg_size;                                    // possible sites in in remaining foreground sequences
  int bg_matches;                                 // matches in remaining background sequences
  int bg_size;                                    // possible sites in remaining background sequences
  int afg_matches;                                // matches in all foreground sequences
  int afg_size;                                   // possible sites in in all foreground sequences
  int abg_matches;                                // matches in all background sequences
  int abg_size;                                   // possible sites in in all background sequences
} MOTIF_INFO_T;

#endif
