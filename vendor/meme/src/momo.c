/********************************************************************
 * MoMo Portal
 ********************************************************************/

#define DEFINE_GLOBALS

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "display.h"
#include "dir.h"
#include "fasta-io.h"
#include "momo.h"
#include "momo-output.h"
#include "io.h"
#include "matrix.h"
#include "motif-in.h"
#include "motif-spec.h"
#include "seq-reader-from-fasta.h"
#include "simple-getopt.h"
#include "utils.h"
#include "binomial.h"
#include "meme.h"
#include "read_seq_file.h"
#include "momo-input.h"
#include "momo-algorithm.h"

char* program_name = "momo";

/* local variables */
#define DATA_HASH_SIZE 100003

// default verbosity
VERBOSE_T verbosity = NORMAL_VERBOSE;

/***********************************************************************
 Initialize options processing
 ***********************************************************************/
static MOMO_OPTIONS_T init_options(ALGORITHM_T algorithm) {
  /* Make sure various options are set to NULL or defaults. */
  MOMO_OPTIONS_T options;

  options.verbosity = 2; 

  options.command_line = NULL;
  options.html_path = NULL;
  options.algorithm = algorithm;

  options.output_dirname = "momo_out";
  options.allow_clobber = true;
  options.text_path = NULL;
  options.tsv_path = NULL;

  options.phospho_filenames = arraylst_create();
  options.phospho_filesizes = arraylst_create();
  options.psm_type = NULL;
  options.sequence_column = NULL;
  
  options.filter = false;
  options.filter_field = NULL;
  options.filter_threshold = 0.0;
  options.eliminate_repeat_width = -1;
  options.min_occurrences = (algorithm == Motifx) ? 20 : 5;
  options.remove_unknowns = (algorithm == Motifx) ? true : false;
  
  options.width = 7;
  options.single_motif_per_mass = false;

  options.protein_database_filename = NULL;
  options.db_background = false;
  options.seed = 0;
  options.hash_fasta = false;
  options.hash_fasta_width = 0;

  options.paper = "Alice Cheng, Charles Grant, Timothy L. Bailey and William Noble, \"MoMo: Discovery of statistically signficant post-translational modification motifs\", Bioinformatics, 35(16):2774-2782, 2018.";
  options.link = "https://doi.org/10.1093/bioinformatics/bty1058";

  // Motif-X options
  options.score_threshold = (algorithm == Motifx) ? 0.000001 : 0.0;
  options.harvard = false;

  // MoDL options
  options.max_motifs = 100;
  options.max_iterations = 50;
  options.max_no_decrease = 10;

  // Hidden options
  options.printp = false;

  // No longer user-settable.
  options.filetype = Psm;
  options.bg_filetype = Fasta;

  return options;
}

/***********************************************************************
 Initialize summmary
 ***********************************************************************/
static SUMMARY_T init_summary() {
  SUMMARY_T summary;
  summary.num_mod = 0;
  summary.num_modtype = 0;
  summary.num_mod_passing = 0;
  summary.num_modtype_passing = 0;
  
  summary.hash_fasta_table = NULL;
  HASH_TABLE mod_table = hash_create(DATA_HASH_SIZE, NULL);	// hash table of mods
  ARRAYLST_T * mod_table_keys = arraylst_create(); // hash table keys
  summary.mod_table = mod_table;
  summary.mod_table_keys = mod_table_keys;
  
  // Initalize the protein alphabet within summary
  ALPH_T* alph = alph_protein();
  summary.alph = alph;
  
  // Get list of amino acids
  STR_T* alph_letters_string = str_create(MAX_ALPH_SIZE);
  const char* alph_letters = alph_string(summary.alph, alph_letters_string);
  summary.alph_letters = alph_letters;
  
  // cleanup
  char* stored_string = str_destroy(alph_letters_string, true);
  return summary;
}

/***********************************************************************
 Free memory allocated in options processing
 ***********************************************************************/
static void cleanup_options(MOMO_OPTIONS_T* options) {
  myfree(options->command_line);
  myfree(options->html_path);
  myfree(options->text_path);
  myfree(options->tsv_path);
  if (options->filter) {
    myfree(options->filter_field);
  }
  arraylst_destroy(NULL, options->phospho_filenames);
  arraylst_destroy(NULL, options->phospho_filesizes);
}

/***********************************************************************
 Free memory allocated in summary
 ***********************************************************************/
static void cleanup_summary(
  MOMO_OPTIONS_T* options,
  SUMMARY_T* summary) {
  int i;

  // Each SUMMARY_T contains: alph, bg_freqs, mod_table/mod_table_keys, hash_fasta_table

  myfree(summary->alph_letters);
  
  // Clean up alph & bg_freqs
  alph_release(summary->alph);
  free_array(summary->bg_freqs);
  
  // Clean up mod_table/mod_table_keys
  for (i = 0; i < arraylst_size(summary->mod_table_keys); ++i) {
    HASH_TABLE_ENTRY * hash_entry = arraylst_get(i, summary->mod_table_keys);
    MOD_INFO_T * modinfo = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
    cleanup_modinfo(options, modinfo);
  }
  arraylst_destroy(NULL, summary->mod_table_keys);
  hash_destroy(summary->mod_table);
  
  // Clean up hash_fasta_table if used.
  if (summary->hash_fasta_table != NULL) hash_destroy(summary->hash_fasta_table);
}

static char* get_usage_message(
  ALGORITHM_T algorithm,
  MOMO_OPTIONS_T *options
) {

  if (algorithm == No_algorithm || algorithm == Unknown_algorithm) {
    return "Usage: momo <algorithm> [options] <arguments>\n"
      "\n"
      "MoMo supports the following algorithms:\n"
      "  simple\n"
      "  motifx\n"
      "  modl\n"
      "\n"
      "Options and arguments are specific to each command.\n\n"
      "Type \'momo <command>\' for details.\n\n";
  }

  char *msg1 = (char *) malloc(2000);
  char *msg2 = (char *) malloc(1000);
  char *msg3 = (char *) malloc(1000);

  // Common stuff.
  sprintf(msg1,
    "Usage: momo %s [options] <ptm file>+\n" 
    "\n"
    "   Options:\n"
    "     --o <output dir> (default: %s)\n"
    "     --oc <output dir> (default: %s)\n"
    "     --psm-type comet|ms-gf+|tide|percolator\n"
    "     --sequence-column [column name]\n"
    "     --width [positive odd integer] (default: %d)\n"
    "     --protein-database <protein sequence file> (default: %s)\n"
    "     --filter [field],lt|le|eq|ge|gt,[threshold] (default: no filter)\n"
    "     --remove-unknowns T|F (default: %c)\n"
    "     --eliminate-repeats [positive odd integer or 0 for no elimination] (default: width)\n"
    "     --min-occurrences [non-negative] (default: %d)\n"
    "     --single-motif-per-mass\n"
    "     --hash-fasta [positive integer or 0 for linear search] (default: %d)\n",
    (algorithm == Simple) ? "simple" : (algorithm == Motifx ? "motifx" : "modl"),
    options->output_dirname, options->output_dirname, 
    options->width,
    options->protein_database_filename ? options->protein_database_filename : "None",
    options->remove_unknowns ? 'T' : 'F', 
    options->min_occurrences,
    options->hash_fasta_width
  );

  if (algorithm == Simple) {
    msg2[0] = '\0';
  } else if (algorithm == Motifx) {
    sprintf(msg2, 
      "     --seed [non-negative integer]\n"
      "     --db-background\n"
      "     --score-threshold [positive value] (default: %.1e)\n"
      "     --harvard (default: compute binomial CDF correctly)\n",
      options->score_threshold
    );
  } else if (algorithm == Modl) {
    sprintf(msg2, 
      "     --seed [non-negative integer]\n"
      "     --db-background\n"
      "     --max-motifs [positive integer] (default: %d)\n"
      "     --max-iterations [positive integer] (default: %d)\n"
      "     --max-no-decrease [positive integer] (default: %d)\n",
      options->max_motifs,
      options->max_iterations,
      options->max_no_decrease
    );
  } else {
    msg2[0] = '\0';
  }

  // Common stuff.
  sprintf(msg3,
    "     --verbosity 1|2|3|4|5 (default: %d)\n"
    "     --version (print the version and exit)\n"
    "\n",
    options->verbosity 
  );

  char *msg = (char *) malloc((strlen(msg1)+strlen(msg2)+strlen(msg3)+1) * sizeof(char));
  strcpy(msg, msg1);
  strcat(msg, msg2);
  strcat(msg, msg3);
  myfree(msg1);
  myfree(msg3);
  myfree(msg3);
  return(msg);
} // get_usage_message

/***********************************************************************
 Process command line options
 ***********************************************************************/
static MOMO_OPTIONS_T process_momo_command_line(
                                              int argc,
                                              char* argv[]
                                              ) {
  
  ALGORITHM_T algorithm = No_algorithm;
  if (argc >= 2) {
    if (strcmp(argv[1], "simple") == 0) {
      algorithm = Simple;
    } else if (strcmp(argv[1], "motifx") == 0) {
      algorithm = Motifx;
    } else if (strcmp(argv[1], "modl") == 0) {
      algorithm = Modl;
    } else {
      algorithm = Unknown_algorithm;
    }
  }
  
  MOMO_OPTIONS_T options = init_options(algorithm);

  // Define the usage message.
  options.usage = get_usage_message(algorithm, &options);
  
  // Define command line options.
  cmdoption const momo_simple_options[] = {
    {"psm-type", REQUIRED_VALUE},
    {"eliminate-repeats", REQUIRED_VALUE},
    {"filter", REQUIRED_VALUE},
    {"hash-fasta", REQUIRED_VALUE},
    {"min-occurrences", REQUIRED_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"protein-database", REQUIRED_VALUE},
    {"remove-unknowns", REQUIRED_VALUE},
    {"sequence-column", REQUIRED_VALUE},
    {"single-motif-per-mass", NO_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"version", NO_VALUE},
    {"width", REQUIRED_VALUE},
  };
  cmdoption const momo_motifx_options[] = {
    {"psm-type", REQUIRED_VALUE},
    {"eliminate-repeats", REQUIRED_VALUE},
    {"filter", REQUIRED_VALUE},
    {"hash-fasta", REQUIRED_VALUE},
    {"min-occurrences", REQUIRED_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"seed", REQUIRED_VALUE},
    {"db-background", NO_VALUE},
    {"protein-database", REQUIRED_VALUE},
    {"remove-unknowns", REQUIRED_VALUE},
    {"score-threshold", REQUIRED_VALUE},
    {"sequence-column", REQUIRED_VALUE},
    {"single-motif-per-mass", NO_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"version", NO_VALUE},
    {"width", REQUIRED_VALUE},
    {"harvard", NO_VALUE},
    {"printp", NO_VALUE},
  };
  cmdoption const momo_modl_options[] = {
    {"psm-type", REQUIRED_VALUE},
    {"eliminate-repeats", REQUIRED_VALUE},
    {"filter", REQUIRED_VALUE},
    {"hash-fasta", REQUIRED_VALUE},
    {"min-occurrences", REQUIRED_VALUE},
    {"max-iterations", REQUIRED_VALUE},
    {"max-motifs", REQUIRED_VALUE},
    {"max-no-decrease", REQUIRED_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"seed", REQUIRED_VALUE},
    {"protein-database", REQUIRED_VALUE},
    {"db-background", NO_VALUE},
    {"remove-unknowns", REQUIRED_VALUE},
    {"sequence-column", REQUIRED_VALUE},
    {"single-motif-per-mass", NO_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"version", NO_VALUE},
    {"width", REQUIRED_VALUE},
    {"printp", NO_VALUE},
  };
  
  // Get the options for the corresponding algorithm
  int option_index = 0;
  
  if (algorithm == Simple) {
    const int num_options = sizeof(momo_simple_options) / sizeof(cmdoption);
    simple_setopt(argc - 1, argv + 1, num_options, momo_simple_options);
  } else if (algorithm == Modl) {
    const int num_options = sizeof(momo_modl_options) / sizeof(cmdoption);
    simple_setopt(argc - 1, argv + 1, num_options, momo_modl_options);
  } else { // algorithm == No_algorithm or algorithm == Motifx
    const int num_options = sizeof(momo_motifx_options) / sizeof(cmdoption);
    simple_setopt(argc - 1, argv + 1, num_options, momo_motifx_options);
  }
  
  // Parse the command line.
  while (true) {
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;
    
    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    }
    else if (c < 0) {
      (void) simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }
    else if (strcmp(option_name, "eliminate-repeats") == 0) {
      options.eliminate_repeat_width = atoi(option_value);
      if (options.eliminate_repeat_width > 0 && options.eliminate_repeat_width % 2 != 1) {
        die("--eliminate-repeats must be odd or 0.");
      }
    }
    else if (strcmp(option_name, "harvard") == 0) {
      options.harvard = true;
    }
    else if (strcmp(option_name, "filter") == 0){
      options.filter = true;
      char * value = strdup(option_value);
      char * pch = strtok(value, ",");
      options.filter_field = pch;
      pch = strtok(NULL, ",");
      if (pch == NULL) {
        die ("Error reading filter field.");
      } else if (strcmp(pch, "le") == 0) {
        options.filter_type = Le;
      } else if (strcmp(pch, "lt") == 0) {
        options.filter_type = Lt;
      } else if (strcmp(pch, "eq") == 0) {
        options.filter_type = Eq;
      } else if (strcmp(pch, "gt") == 0) {
        options.filter_type = Gt;
      } else if (strcmp(pch, "ge") == 0) {
        options.filter_type = Ge;
      } else {
        die ("Error reading filter comparison type.");
      }
      pch = strtok(NULL, ",");
      if (pch == NULL) {
        die ("Error reading filter threshold.");
      }
      options.filter_threshold = atof(pch);
    }
    else if (strcmp(option_name, "hash-fasta") == 0){
      options.hash_fasta_width = atoi(option_value);
      options.hash_fasta = (options.hash_fasta_width > 0) ? true : false;
    }
    else if (strcmp(option_name, "max-iterations") == 0){
      options.max_iterations = atoi(option_value);
    }
    else if (strcmp(option_name, "max-motifs") == 0){
      options.max_motifs = atoi(option_value);
    }
    else if (strcmp(option_name, "max-no-decrease") == 0){
      options.max_no_decrease = atoi(option_value);
    }
    else if (strcmp(option_name, "min-occurrences") == 0){
      options.min_occurrences = atoi(option_value);
      if (options.min_occurrences <= 0) die("Value of --min-occurrences must be greater than 0.");
    }
    else if (strcmp(option_name, "o") == 0){
      // Set output directory with no clobber
      options.output_dirname = option_value;
      options.allow_clobber = false;
    }
    else if (strcmp(option_name, "oc") == 0){
      // Set output directory with clobber
      options.output_dirname = option_value;
      options.allow_clobber = true;
    }
    else if (strcmp(option_name, "protein-database") == 0){
      options.protein_database_filename = option_value;
    }
    else if (strcmp(option_name, "db-background") == 0){
      options.db_background = true;
    }
    else if (strcmp(option_name, "seed") == 0){
      options.seed = atoi(option_value);
    }
    else if (strcmp(option_name, "score-threshold") == 0){
      options.score_threshold = atof(option_value);
    }
    else if (strcmp(option_name, "psm-type") == 0){
      options.psm_type = option_value;
    }
    else if (strcmp(option_name, "sequence-column") == 0){
      options.sequence_column = option_value;
    }
    else if (strcmp(option_name, "single-motif-per-mass") == 0){
      options.single_motif_per_mass = true;
    }
    else if (strcmp(option_name, "remove-unknowns") == 0){
      if (!(strcmp(option_value, "T") == 0 || strcmp(option_value, "F") == 0)) {
        die("Error reading remove-unknowns");
      }
      options.remove_unknowns = (strcmp(option_value, "T") == 0);
    }
    else if (strcmp(option_name, "verbosity") == 0){
      options.verbosity = atoi(option_value);
    }
    else if (strcmp(option_name, "version") == 0) {
      fprintf(stdout, VERSION "\n");
      exit(EXIT_SUCCESS);
    }
    else if (strcmp(option_name, "width") == 0){
      options.width = atoi(option_value);
      if (options.width % 2 != 1) {
        die("Width must be odd");
      }
    }
    else if (strcmp(option_name, "printp") == 0) {
      options.printp = true;
    }
  }
  option_index++;
  
  // Must have algorithm
  if (algorithm == No_algorithm) {
    fprintf(stderr, "\nYou must specify the algorithm!\n\n");
    fprintf(stderr, "%s", options.usage);
    exit(EXIT_FAILURE);
  } else if (algorithm == Unknown_algorithm) {
    fprintf(stderr, "\nYou must specify one of the known algorithms!\n\n");
    fprintf(stderr, "%s", options.usage);
    exit(EXIT_FAILURE);
  }
  
  // Must have peptide-spectrum-match (PSM) file name
  if (argc < option_index + 1) {
    fprintf(stderr, "\nYou must specify at least one PTM file!\n\n");
    fprintf(stderr, "%s", options.usage);
    exit(EXIT_FAILURE);
  }
  
  // Set the eliminate_repeat_width to width if it is not given. 
  if (options.eliminate_repeat_width == -1) {
    options.eliminate_repeat_width = options.width;
  }

  // Set the sequence_column if it is not given and
  // --psm-type is given.
  if (options.sequence_column == NULL) {
    if (options.psm_type != NULL) {
      if (strcmp(options.psm_type, "comet") == 0) {
        options.sequence_column = "modified sequence";
      } else if (strcmp(options.psm_type, "ms-gf+") == 0) {
        options.sequence_column = "Peptide";
      } else if (strcmp(options.psm_type, "tide") == 0) {
        options.sequence_column = "sequence";
      } else if (strcmp(options.psm_type, "percolator") == 0) {
        options.sequence_column = "sequence";
      } else {
	fprintf(stderr, "Unknown value for --psm-type %s\n", options.psm_type);
	fprintf(stderr, "%s", options.usage);
	exit(EXIT_FAILURE);
      }
    }
  }
  
  // --db-background requires --protein-database
  if (options.db_background && ! options.protein_database_filename) {
    die("You must use option --protein-database when you specify option --db-background!");
  }

  // set random number generators
  srand_mt(options.seed);
  set_randfunc((randfunc_t) random_mt); // for ushuffle 

  // Record the command line
  options.command_line = get_command_line(argc, argv);
  
  // Record the input file names
  while (option_index != argc) {
    arraylst_add(argv[option_index], options.phospho_filenames);
    option_index++;
  }
  
  // Set up path values for needed stylesheets and output files.
  options.HTML_FILENAME = "momo.html";
  options.TEMPLATE_FILENAME = "momo_template.html";
  options.TEXT_FILENAME = "momo.txt";
  options.TSV_FILENAME = "momo.tsv";
  
  options.html_path = make_path_to_file(options.output_dirname, options.HTML_FILENAME);
  options.text_path = make_path_to_file(options.output_dirname, options.TEXT_FILENAME);
  options.tsv_path = make_path_to_file(options.output_dirname, options.TSV_FILENAME);
  
  return options;
  
}

/**
 * Sets the background frequencies from the given sequences.
 */
static void initialize_background_frequencies(MOMO_OPTIONS_T* options,
                                            SUMMARY_T* summary,
                                            SEQ_T** bg_sequences,
                                            int num_bg_sequences) {
  int i, j;
  ARRAY_T* bg_freqs = NULL;
  const char* alph_letters = summary->alph_letters;
  
  bg_freqs = allocate_array(strlen(alph_letters));
  init_array(0.0, bg_freqs);
  
  for (i = 0; i < num_bg_sequences; ++i) {
    char* currseq = get_raw_sequence(bg_sequences[i]);
    for (j = 0; j < strlen(currseq); ++j) {
      if (strchr(alph_letters, currseq[j])) {
	int idx = strchr(alph_letters, currseq[j]) - alph_letters;
	if (idx >= 0) {
	  set_array_item(idx, get_array_item(idx, bg_freqs) + 1, bg_freqs);
	}
      }
    }
  }
  normalize(0.0, bg_freqs);
  summary->bg_freqs = bg_freqs;
}

/**
 * Analyzes the Peptide-Spectrum Match (PSM) files, create motifs, and returns a summary
 * on the results
 */
static SUMMARY_T get_summary(MOMO_OPTIONS_T* options) {
  int i;
  
  // Initialize the summary object
  SUMMARY_T summary = init_summary();

  // If a protein database is provided and in FASTA format, we will read it into
  // an array of sequences.
  SEQ_T** bg_sequences = NULL;
  int num_bg_sequences = 0;
  if (options->protein_database_filename) {
    read_protein_database_sequences(options, &summary, &bg_sequences, &num_bg_sequences);
    // See momo-input.c
    // If a protein database is provided and the user specified to use
    // an O(1) lookup table to speed up the process of finding a peptide
    // within the protein database, we initialize the O(1) lookup table.
    create_hash_fasta_preprocess_table(options, &summary, bg_sequences, num_bg_sequences);
  }
  
  // See momo-input.c
  // For each PSM file, we will add the information to a mod table. Each mod is hashed
  // to its mod entry.
  SEQ_T** fg_sequences = NULL;
  int num_fg_sequences = 0;
  add_phospho_files_to_table(options, &summary, bg_sequences, num_bg_sequences, &fg_sequences, &num_fg_sequences);

  // Free up the memory used by background sequences if they are no longer needed.
  if (! options->db_background && bg_sequences) {
    for (i = 0; i < num_bg_sequences; ++i) free_seq(bg_sequences[i]);
    myfree(bg_sequences);
    bg_sequences = NULL;
  }

  // Initialize background frequencies.
  if (options->db_background) {
    // Using the protein database for the background, so
    // initialize background frequencies from the protein database
    // and save the background peptides in the table.
    initialize_background_frequencies(options, &summary, bg_sequences, num_bg_sequences);
    add_background_sequences_to_table(options, &summary, bg_sequences, num_bg_sequences);
    //fprintf(stderr, "num_bg_sequences %d\n", num_bg_sequences);
  } else {
    // Using foreground for background frequencies.
    initialize_background_frequencies(options, &summary, fg_sequences, num_fg_sequences);
    // Create new (identical) set of shuffled foreground sequences to use as background peptides.
    srand_mt(options->seed);
    set_randfunc((randfunc_t) random_mt); // for ushuffle 
    add_shuffled_sequences_to_table(options, &summary, true);
  }

  // See momo-algorithm.c
  // Using the mod table, we will generate frequencies and create motifs for each mod.
  create_motifs(options, &summary);
  
  // Free up the memory used by background sequences.
  if (bg_sequences) {
    for (i = 0; i < num_bg_sequences; ++i) free_seq(bg_sequences[i]);
    myfree(bg_sequences);
  }

  // Free up the memory used by foreground sequences.
  if (fg_sequences) {
    for (i = 0; i < num_fg_sequences; ++i) free_seq(fg_sequences[i]);
    myfree(fg_sequences);
  }
  
  return summary;
}

/*************************************************************************
 * Entry point for momo
 *************************************************************************/
int main(int argc, char *argv[]) {
  // Start timing
  clock_t start = clock(), diff;
  
  // Get command line arguments
  MOMO_OPTIONS_T options = process_momo_command_line(argc, argv);
  
  // Create motifs and obtain summary
  SUMMARY_T summary = get_summary(&options);
  
  // Print results
  print_momo_results(argc, argv, options, summary);

  // Clean up.
  cleanup_summary(&options, &summary);
  cleanup_options(&options);
  
  // Print timing
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  //  printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
  
  return 0;
}
