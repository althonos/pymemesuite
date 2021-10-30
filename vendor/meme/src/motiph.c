/********************************************************************
 * FILE: motiph.c
 * AUTHOR: William Stafford Noble, Charles E. Grant, Timothy L. Bailey
 * CREATE DATE: 12/17/2004
 * PROJECT: EVOMCAST
 * COPYRIGHT: 2004, UW
 ********************************************************************/

#define DEFINE_GLOBALS
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "alignment.h"
#include "alphabet.h"
#include "cisml.h"
#include "dir.h"
#include "evomodel.h"
#include "fasta-io.h"
#include "mhmm-state.h"
#include "motif.h"
#include "motif-in.h"
#include "motiph-scoring.h"
#include "pssm.h"
#include "seq.h"
#include "simple-getopt.h"
#include "tree.h"
#include "object-list.h"

char* program_name = NULL;
VERBOSE_T verbosity = NORMAL_VERBOSE;

/*******************************************************************
  Print the column frequency distribution.
 ********************************************************************/
static void print_col_frequencies(
  ALPH_T *alph,
  ARRAY_T* alignment_column_freqs
)
{
  int i;
  int num_freqs = get_array_length(alignment_column_freqs);
  int asize = alph_size_core(alph);
  int num_leaves = NINT(log(num_freqs)/log(asize));
  char* alignment_col = mm_malloc((num_leaves + 1) * sizeof(char));
  for (i=0; i<num_freqs; i++) {
    unhash_alignment_col(
      alph,
      i,                              //col_index
      alignment_col,
      num_leaves
    );
    printf("%s %d %g\n", alignment_col, i+1,
      get_array_item(i, alignment_column_freqs));
  }
} // print_col_freqs

/*************************************************************************
 * Entry point for motiph
 *************************************************************************/
int main(int argc, char *argv[]) {

  bool allow_clobber = true;
  bool compute_qvalues = true;
  bool print_col_freqs = false;
  bool print_trimmed_tree = false;
  bool scan_both_strands = true;
  bool bls_both_strands = true;
  bool output_qthresh_set = false;
  bool text_only = false;
  bool use_file_list = false;
  bool use_halpern_bruno = false;
  bool use_BLS = false;
  int bls_distance = 20;
  int ref_seq_index = 0;  // Which sequence provides reference coords.
  int max_stored_scores = 100000;
  //FIXME: Eventually, above should be accessible via a command-line switch.
  char* command_line = NULL;
  char* bg_filename = NULL;
  char *output_dirname = "motiph_out";
  char* ustar_label = NULL;     // TLB; create uniform star tree
  double bg_rate = 1.0;
  double fg_rate = 1.0;
  double gap_cost = 0.0;
  double pseudocount = 0.1;
  double output_pthresh = 1e-4;
  double output_qthresh = 1.0;
  double purine_pyrimidine = 1.0; // r
  double transition_transversion = 0.5; // R
  GAP_SUPPORT_T gap_support = SKIP_GAPS;
  MODEL_TYPE_T model_type = F81_MODEL;
  COLUMN_FREQS_TYPE_T column_freqs_type = SIMULATED_COL_FREQS;
  STRING_LIST_T* selected_motifs = NULL;
  OBJECT_LIST_T* alignment_col_freqs_list = NULL;
  bool report_err = false;

  program_name = "motiph";

  /**********************************************
   * COMMAND LINE PROCESSING
   **********************************************/

  // Define command line options.
  // FIXME: Note that if you add or remove options you
  // must change n_options.
  int n_options = 29;
  cmdoption const motiph_options[] = {
    {"bg", REQUIRED_VALUE},
    {"bfile", REQUIRED_VALUE},
    {"bgfile", REQUIRED_VALUE},
    {"bls-dist", REQUIRED_VALUE},
    {"column-freqs", REQUIRED_VALUE},
    {"fg", REQUIRED_VALUE},
    {"flip", REQUIRED_VALUE}, //bls_both_strands
    {"gap", REQUIRED_VALUE},
    {"gap-cost", REQUIRED_VALUE},
    {"hb", NO_VALUE},
    {"list", NO_VALUE},
    {"max-stored-scores", REQUIRED_VALUE},
    {"model", REQUIRED_VALUE},
    {"motif", REQUIRED_VALUE},
    {"no-qvalue", NO_VALUE},
    {"norc", NO_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"output-pthresh", REQUIRED_VALUE},
    {"output-qthresh", REQUIRED_VALUE},
    {"print-col-freqs", NO_VALUE},
    {"print-trimmed-tree", NO_VALUE},
    {"pseudocount", REQUIRED_VALUE},
    {"pur-pyr", REQUIRED_VALUE},
    {"seed", REQUIRED_VALUE},
    {"text", NO_VALUE},
    {"transition-transversion", REQUIRED_VALUE},
    {"ustar", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE}
  };

  int option_index = 0;

  // Define the usage message.
  char *usage = 
  "Usage: motiph [options] <alignment> <tree> <motif>\n"
  "\n"
  "   Options:\n"

  // Evolutionary model parameters.
  "     --bg <float> (default=1.0)\n"
  "     --column-freqs simulated|empirical (default=simulated)\n"
  "     --fg <float> (default=1.0)\n"
  "     --gap skip | fixed | wildcard | minimum"
  " (default=skip)\n"
  "     --gap-cost <float> (default=0.0)\n"
  "     --hb (defalut false)\n"
  "     --model bls|single|average|jc|k2|f81|f84|hky|tn"
  " (default=f81)\n"
  "     --pur-pyr <float> (default=1.0)\n"
  "     --transition-transversion <float> (default=0.5)\n"
  "     --ustar <label>\n"      // TLB; create uniform star tree

  // Motif parameters.
  "     --motif <id> (default=all)\n"

  // Miscellaneous parameters
  "     --bfile <background> (default from motif file)\n"
  "     --list\n"
  "     --no-qvalue (default false)\n"
  "     --norc (default false)\n"
  "     --flip true|false (Allow BLS matches in reverse. Default: true)\n"
  "     --bls-dist <int> (Distance threshold for BLS. Default:20)\n"
  "     --max-stored-scores (Default 100,000)\n"
  "     --o <output dir> (default=motiph_out)\n"
  "     --oc <output dir> (default=motiph_out)\n" 
  "     --print-col-freqs\n"
  "     --print-trimmed-tree\n"
  "     --pseudocount <float> (default=0.1)\n"
  "     --output-pthresh <float> (default 1e-4)\n"
  "     --output-qthresh <float> (default 1.0)\n"
  "     --seed <long> (default from time())\n"
  "     --text (default false)\n"
  "     --verbosity [1|2|3|4] (default 2)\n"
  "\n";


  // Parse the command line.
  if (simple_setopt(argc, argv, n_options, motiph_options) != NO_ERROR) {
    die("Error processing command line options: option name too long.\n");
  }

  // Default seed the random number generator.
  long seed = time(0);

  while (true) {
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
      (void) simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "bg") == 0){
      bg_rate = atof(option_value);
    }
    else if (strcmp(option_name, "bfile") == 0 || strcmp(option_name, "bgfile") == 0){
      bg_filename = option_value;
    }     
    else if (strcmp(option_name, "bls-dist") == 0){
      bls_distance = atoi(option_value);
    }   
    //else if (strcmp(option_name, "bls") == 0){
    //  use_BLS = true;
    //} 
    else if (strcmp(option_name, "column-freqs") == 0) {
      if (strcmp(option_value, "simulated") == 0) {
        column_freqs_type = SIMULATED_COL_FREQS;
      } else if (strcmp(option_value, "empirical") == 0) {
        column_freqs_type = EMPIRICAL_COL_FREQS;
      } else {
        die("Unknown column-freqs type: %s\n", option_value);
      }
    }
    else if (strcmp(option_name, "fg") == 0){
      fg_rate = atof(option_value);
    }
    else if (strcmp(option_name, "flip") == 0) {
      if (strcmp(option_value, "false") == 0) {
        bls_both_strands = false;
      } else {
        bls_both_strands = true;
      }
    } 
    else if (strcmp(option_name, "gap") == 0){
      if (strcmp(option_value, "skip") == 0) {
        gap_support = SKIP_GAPS;
      } else if (strcmp(option_value, "fixed") == 0) {
        gap_support = FIXED_GAP_COST;
      } else if (strcmp(option_value, "wildcard") == 0) {
        gap_support = WILDCARD_GAP;
      } else if (strcmp(option_value, "minimum") == 0) {
        gap_support = MIN_GAPS;
      } else {
        die("Unknown gap handling method: %s\n", option_value);
      }
    }
    else if (strcmp(option_name, "gap-cost") == 0) {
      gap_cost = atof(option_value);
    }
    else if (strcmp(option_name, "hb") == 0){
        use_halpern_bruno = true;
    }
    else if (strcmp(option_name, "list") == 0){
      use_file_list = true;
    }
    else if (strcmp(option_name, "model") == 0) {
      if (strcmp(option_value, "jc") == 0) {
        model_type = JC_MODEL;
      } else if (strcmp(option_value, "k2") == 0) {
        model_type = K2_MODEL;
      } else if (strcmp(option_value, "f81") == 0) {
        model_type = F81_MODEL;
      } else if (strcmp(option_value, "f84") == 0) {
        model_type = F84_MODEL;
      } else if (strcmp(option_value, "hky") == 0) {
        model_type = HKY_MODEL;
      } else if (strcmp(option_value, "tn") == 0) {
        model_type = TAMURA_NEI_MODEL;
      } else if (strcmp(option_value, "single") == 0) {
        model_type = SINGLE_MODEL;

      } else if (strcmp(option_value, "bls") == 0) {
        model_type = SINGLE_MODEL;
        use_BLS = true;

      } else if (strcmp(option_value, "average") == 0) {
        model_type = AVERAGE_MODEL;
      } else {
        die("Unknown model: %s\n", option_value);
      }
    }
    else if (strcmp(option_name, "max-stored-scores") == 0) {
      // Use atof and cast to be able to read things like 1e8.
      max_stored_scores = (int)atof(option_value);
    }
    else if (strcmp(option_name, "motif") == 0){
      if (selected_motifs == NULL) {
        selected_motifs = new_string_list();
      }
      add_string(option_value, selected_motifs);
    }
    else if (strcmp(option_name, "no-qvalue") == 0) {
      compute_qvalues = false;
    }
    else if (strcmp(option_name, "norc") == 0) {
      scan_both_strands = false;
    }
    else if (strcmp(option_name, "o") == 0){
      // Set output directory with no clobber
      output_dirname = option_value;
      allow_clobber = false;
    }
    else if (strcmp(option_name, "oc") == 0){
      // Set output directory with clobber
      output_dirname = option_value;
      allow_clobber = true;
    }
    else if (strcmp(option_name, "output-pthresh") == 0){
      output_pthresh = atof(option_value);
      output_qthresh = 1.0;
    }
    else if (strcmp(option_name, "output-qthresh") == 0){
      output_qthresh = atof(option_value);
      output_pthresh = 1.0;
      output_qthresh_set = true;
    }
    else if (strcmp(option_name, "print-col-freqs") == 0){
      print_col_freqs = true;
    }

    else if (strcmp(option_name, "print-trimmed-tree") == 0){
      print_trimmed_tree = true;
    }
    else if (strcmp(option_name, "pseudocount") == 0){
      pseudocount = atof(option_value);
    }
    else if (strcmp(option_name, "pur-pyr") == 0){
      purine_pyrimidine = atof(option_value);
    }
    else if (strcmp(option_name, "seed") == 0){
      seed = atol(option_value);
    }
    else if (strcmp(option_name, "text") == 0){
      text_only = true;
    }
    else if (strcmp(option_name, "transition-transversion") == 0){
      transition_transversion = atof(option_value);
    }
    else if (strcmp(option_name, "ustar") == 0){        // TLB; create uniform star tree
      ustar_label = option_value;
    }
    else if (strcmp(option_name, "verbosity") == 0){
      verbosity = atoi(option_value);
    }
  }
 
  // Record the command line
  command_line = get_command_line(argc, argv);

  // Initialize random number generator
  srand_mt(seed);

  // Check that the qvalue options are consistent
  if (compute_qvalues == false && output_qthresh_set == true) {
    die("The --no-qvalue option cannot be used with the --output-qthresh options");
  }

  // Must have tree and alignment file names
  if (argc != option_index + 3) {
    fprintf(stderr, "%s", usage);
    exit(EXIT_FAILURE);
  }

  /****************************************************
   * Read the names of the alignments or sequence files
   ****************************************************/
  STRING_LIST_T* filenames = NULL;
  int num_filenames = 0;
  if (use_file_list) {
    filenames = read_string_list_from_file(argv[option_index]);
    num_filenames = get_num_strings(filenames);
  }
  else {
    filenames = new_string_list();
    add_string(argv[option_index], filenames);
    num_filenames = 1;
  }
  option_index++;

  /**********************************************
   * Read the phylogenetic tree.
   **********************************************/
  char* tree_filename = argv[option_index];
  option_index++;
  TREE_T* tree = read_tree_from_file(tree_filename);

  // Verify that the tree doesn't have duplicate species.
  STRING_LIST_T* tree_ids = make_leaf_list(tree);
  if (has_duplicates("Duplicate IDs in tree:", tree_ids)) {
    exit(1);
  }
  free_string_list(tree_ids);

  /**********************************************
   * Read the motifs.
   **********************************************/
  char* meme_filename = argv[option_index];
  ARRAYLST_T *motifs;
  ARRAY_T* bg_freqs = NULL;
  option_index++;

  MREAD_T *mread;
  ALPH_T *alph;

  // create motif parser
  mread = mread_create(meme_filename, OPEN_MFILE, scan_both_strands);
  mread_set_bg_source(mread, bg_filename, NULL);
  mread_set_pseudocount(mread, pseudocount);
  // read motifs
  motifs = mread_load(mread, NULL);
  alph = alph_hold(mread_get_alphabet(mread));
  bg_freqs = mread_get_background(mread);
  // destroy motif parser
  mread_destroy(mread);
  // check motifs
  if (arraylst_size(motifs) == 0) die("No motifs in motif file.");
  // TODO FIXME check that this is really required
  if (!alph_is_builtin_dna(alph)) die("Expected DNA motifs");
  assert(bg_freqs != NULL);

  if (scan_both_strands == true ) {
    calc_ambigs(alph, false, bg_freqs);
    // Make reverse complement motifs.
    add_reverse_complements(motifs);
  }

  // OK If we use the BLS Scan method we need to set up the inverse
  // motif right from the start so that the scans are done together.


  // TLB; need to resize bg_freqs array to 4 items
  // or copy array breaks in HB mode.  This throws away
  // the freqs for the ambiguous characters;
  // FIXME: fix this so that ambigs are allowed
  resize_array(bg_freqs, 4);

  // Only SIMULATED column frequencies are supported for SINGLE model.
  if (model_type == SINGLE_MODEL && column_freqs_type != SIMULATED_COL_FREQS) {
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Using simulated column frequencies.\n");
    }
    column_freqs_type = SIMULATED_COL_FREQS;
  }

  // Compute the column frequencies distributions from the input
  // alignments if necessary.  Remove all-gap sequences before
  // computing frequencies.
  if (column_freqs_type == EMPIRICAL_COL_FREQS) {
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(
        stderr,
        "Computing column frequencies from multiple alignment file(s).\n"
      );
    }
    alignment_col_freqs_list = get_alignment_column_freqs_list(alph, filenames, true);
  } else if (column_freqs_type == PRECOMPUTED_COL_FREQS) {
    // FIXME: add switch --column-freqs-file <file> and read that
    // file of column frequency distributions for different sets of species.
    die("PRECOMPUTED column frequencies not yet supported.\n");
  } else {
    // Using SIMULATED column frequencies.  Initialize list.
    alignment_col_freqs_list = new_object_list(equal_string_lists,
        (void*)copy_string_list,
        free_string_list,
        free_array);
  } // column freqs list

  // Create cisml data structure for recording results
  CISML_T *cisml = allocate_cisml("motiph", command_line, meme_filename, "clustal-w alignment");
  set_cisml_site_pvalue_cutoff(cisml, output_pthresh);
  set_cisml_site_qvalue_cutoff(cisml, output_qthresh);

  /**************************************************************
  * Score each of the alignments for each of the selected motifs.
  **************************************************************/
  int n_scored = 0;    // number of positions that were scored
  int motif_index;
  for (motif_index = 0; motif_index < arraylst_size(motifs); motif_index++) {

    MOTIF_T* motif = (MOTIF_T*)arraylst_get(motif_index, motifs);
    char* motif_id = (scan_both_strands ? get_motif_st_id(motif) : get_motif_id(motif));
    char* bare_motif_id = motif_id;

    // Create cisml pattern and add to cisml record
    PATTERN_T *pattern = allocate_pattern(motif_id, motif_id);
    set_pattern_max_stored_matches(pattern, max_stored_scores);
    add_cisml_pattern(cisml, pattern);

    // We may have specified on the command line that
    // only certain motifs were to be used.
    if (selected_motifs != NULL) {
      if (*bare_motif_id == '+' || *bare_motif_id == '-') {
        // The selected  motif id won't included a strand indicator.
        bare_motif_id++;
      }
      if (have_string(bare_motif_id, selected_motifs) == false) {
        continue;
      }
    }

    int window_size = get_motif_length(motif);
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(
        stderr,
        "Using motif %s of width %d.\n",
        motif_id, get_motif_length(motif)
      );
    }

    // Build an array of evolutionary models for each position in the motif.
    EVOMODEL_T** models = make_motif_models(
      motif,
      bg_freqs,
      model_type,
      fg_rate,
      bg_rate,
      purine_pyrimidine,
      transition_transversion,
      use_halpern_bruno
    );

    PSSM_T* pssm = NULL;
    ARRAY_T* alignment_col_freqs = NULL;
    OBJECT_LIST_T* pssm_list = new_object_list(equal_string_lists,
        (void*)copy_string_list,
        free_string_list,
        free_pssm);

    SCORED_SITES_T* site_list = new_scored_sites();

    /* -------------------------------------------------------------*/
    // IF we are using a BLS model and considering the inverse motif 
    // as a valid conserved instance (the flip flag) then we need to
    // build the models to scan the inverse motif at the same time.
    /* -------------------------------------------------------------*/
    MOTIF_T* inverse_motif = NULL;
    EVOMODEL_T** inverse_models = NULL;
    ARRAY_T* inverse_alignment_col_freqs = NULL;
    PSSM_T * inverse_pssm = NULL;
    OBJECT_LIST_T* inverse_pssm_list = NULL;

    if( use_BLS && bls_both_strands ) {
        // if the motif index is even then the inverse comes after,
        // otherwise it is before it.
        int inverse_index = motif_index + 1;
        if( (motif_index % 2) > 0)
          inverse_index = motif_index - 1;

        inverse_motif = (MOTIF_T*)arraylst_get(inverse_index, motifs);
        inverse_models = make_motif_models(
            inverse_motif, 
            bg_freqs,
            model_type,
            fg_rate, 
            bg_rate, 
            purine_pyrimidine, 
            transition_transversion, 
            use_halpern_bruno
        );

        inverse_pssm_list = new_object_list(equal_string_lists,
            (void*)copy_string_list,
            free_string_list,
            free_pssm);
    }

    // Consider each alignment in turn.
    int file_index;
    for(file_index = 0; file_index < num_filenames; file_index++) {

      int current_ref_seq_index = ref_seq_index;

      // Get the next alignment, removing any all-gap sequences and sorting
      // species alphabetically.
      char* filename = get_nth_string(file_index, filenames);
      ALIGNMENT_T* alignment = read_alignment_from_file(filename,
          true, // Sort?
          true, // Remove gaps?
          &current_ref_seq_index);
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(stderr, "Read alignment from %s.\n", filename);
      }

      // Make sure there are no duplicates.
      STRING_LIST_T* alignment_species = get_species_names(alignment);
      if (has_duplicates("Duplicate IDs in alignment:", alignment_species)) {
        exit(1);
      }

      // Create a scanned_sequence record and record it in pattern.
      SCANNED_SEQUENCE_T *scanned_seq =
        allocate_scanned_sequence(filename, filename, pattern);
      set_scanned_sequence_length(scanned_seq, get_alignment_length(alignment));

      // Trim the tree, eliminating species not in this alignment.
      TREE_T* trimmed_tree = trim_tree(true, tree, alignment_species);

      // Check that at least one species was found.
      if (trimmed_tree == NULL) {
        die("Your tree doesn't contain any of the species in your alignment.");
      }
      // Just print tree and exit?
      if (print_trimmed_tree) {
        printf("%s ", filename);
        write_tree(trimmed_tree, stdout);
        continue;
      }
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(stderr, "Trimmed tree: ");
        write_tree(trimmed_tree, stderr);
      }
      STRING_LIST_T* tree_species = make_leaf_list(trimmed_tree);
      sort_string_list(tree_species);   // keep species alphabetical

      // Trim the alignment, eliminating species not in this tree.
      ALIGNMENT_T* trimmed_alignment = remove_alignment_seqs(tree_species, alignment);
      free_alignment(alignment);
      free_string_list(tree_species);
      free_string_list(alignment_species);
      alignment_species = (get_species_names(trimmed_alignment));
      sort_string_list(alignment_species);      // sort species alphabetically
        
      // Try to retrieve the pre-computed PSSM from memory.
      pssm = (PSSM_T*)retrieve_object(alignment_species, pssm_list);
      if (use_BLS && bls_both_strands) {
        // Try to retrieve the pre-computed inverse PSSM from memory.
        inverse_pssm = (PSSM_T*)retrieve_object(alignment_species, inverse_pssm_list);
      }


      // JH - I am going to check whether this needs to be done now
      // so it gets done just the once.
      if (pssm == NULL || (use_BLS && bls_both_strands && inverse_pssm == NULL)) {
        // TLB; Convert the tree to a uniform star tree with
        // the target sequence at its center.
        if (ustar_label != NULL) {
          trimmed_tree = convert_to_uniform_star_tree(trimmed_tree, ustar_label);
          if (trimmed_tree == NULL) 
            die("Tree or alignment missing target %s\n", ustar_label);
          if (verbosity >= NORMAL_VERBOSE) {
            fprintf(
                stderr,
                "Target %s placed at center of uniform star tree:\n", 
                ustar_label
            );
            write_tree(trimmed_tree, stderr);
          }
        }
      }

      // Did we fail to retrieve the PSSM?
      if (pssm == NULL) {
        if (verbosity > NORMAL_VERBOSE) {
          fprintf(stderr, "Building PSSM for %s.\n",
          combine_string_list(alignment_species, "-"));
        }

        // Build PSSM *matrix* for motiph
        MATRIX_T* pssm_matrix = build_alignment_pssm_matrix(
          alph,
          alignment_species,
          window_size + 1,
          models,
          trimmed_tree,
          gap_support
        );

        // Get the column frequency distribution for computing p-values.
        alignment_col_freqs =
          retrieve_object(alignment_species, alignment_col_freqs_list);
        // Column frequencies not found?
        if (alignment_col_freqs == NULL) {
          if (column_freqs_type != SIMULATED_COL_FREQS) {
            die("No column frequency distribution was given for species: %s\n",
              combine_string_list(alignment_species, "-"));
          }
          // Get the simulated column frequencies from first row of PSSM.
          // They are based on the single-letter frequencies and the evo model.
          alignment_col_freqs = allocate_array(get_num_cols(pssm_matrix));
          copy_array(get_matrix_row(0, pssm_matrix), alignment_col_freqs);
          // Save the column frequencies for this species list for future use.
          store_object(
            (void*)alignment_col_freqs,
            (void*)alignment_species,
            0, // Score
            alignment_col_freqs_list
          );
          if (print_col_freqs && column_freqs_type == SIMULATED_COL_FREQS) {
            print_col_frequencies(alph, alignment_col_freqs);
            exit(1);
          }
        } else if (print_col_freqs) {
          print_col_frequencies(alph, alignment_col_freqs);
          exit(1);
        }

        // Remove first row from the PSSM.  It was saved above if needed.
        remove_matrix_row(0, pssm_matrix);

        // Create a PSSM object. Build tables to translate log-odds scores to p-values.
        pssm = build_matrix_pssm(
          alph,
          pssm_matrix, 
          alignment_col_freqs,
          NULL, // No prior dist
          1.0, // PSP paramater alpha not used
          PSSM_RANGE 
        );
        free_matrix(pssm_matrix);

        // Store them for later use.
        store_object((void*)pssm, (void*)alignment_species, 0, pssm_list);
      }

      // We got the PSSM, so get the other stuff that goes with it.
      else {
        if (verbosity > NORMAL_VERBOSE) {
          fprintf(stderr, "Retrieved PSSM for %s.\n",
          combine_string_list(alignment_species, "-"));
        }

        alignment_col_freqs = (ARRAY_T*) retrieve_object(
          alignment_species,
          alignment_col_freqs_list
        );
      }

      /* ----------------------------------------------------------------*/
      // We need to do all of the above again for inverse motif if we are
      // doing BLS on both strands
      /* ----------------------------------------------------------------*/
      if (use_BLS && bls_both_strands) {
        if (inverse_pssm == NULL) {
          if (verbosity > NORMAL_VERBOSE) {
            fprintf(stderr, "Building Inverse PSSM for %s.\n",
            combine_string_list(alignment_species, "-"));
          }

          // Build PSSM for Inverse Motif
          MATRIX_T* inverse_pssm_matrix = build_alignment_pssm_matrix(
            alph,
            alignment_species,
            window_size + 1, 
            inverse_models, 
            trimmed_tree, 
            gap_support
          );

          // Get the column frequency distribution for computing p-values.
          inverse_alignment_col_freqs = 
            retrieve_object(alignment_species, alignment_col_freqs_list);
          // Column frequencies not found?
          if (alignment_col_freqs == NULL) {
            if (column_freqs_type != SIMULATED_COL_FREQS) {
              die("No column frequency distribution was given for species: %s\n",
              combine_string_list(alignment_species, "-"));
            }
            // Get the simulated column frequencies from first row of PSSM.
            // They are based on the single-letter frequencies and the evo model.
            inverse_alignment_col_freqs = allocate_array(get_num_cols(inverse_pssm_matrix)); 
            copy_array(get_matrix_row(0, inverse_pssm_matrix), inverse_alignment_col_freqs);
            // Save the column frequencies for this species list for future use.
            store_object(
              (void*)inverse_alignment_col_freqs,
              (void*)alignment_species, 
              0, // Score
              alignment_col_freqs_list
            );
            if (print_col_freqs && column_freqs_type == SIMULATED_COL_FREQS) {
              print_col_frequencies(alph, inverse_alignment_col_freqs);
              exit(1);
            }
          } else if (print_col_freqs) {
            print_col_frequencies(alph, inverse_alignment_col_freqs);
            exit(1);
          }

          // Remove first row from the PSSM.  It was saved above if needed.
          remove_matrix_row(0, inverse_pssm_matrix);

          // Build tables to translate log-odds scores to p-values
          inverse_pssm = build_matrix_pssm(
            alph,
            inverse_pssm_matrix, 
            inverse_alignment_col_freqs,
            NULL, // No priors
            1.0, // PSP parameter alpha not used
            PSSM_RANGE 
          );
          free_matrix(inverse_pssm_matrix);

          // Store them for later use.
          store_object((void*)inverse_pssm, (void*)alignment_species, 0, inverse_pssm_list);

        } else {
          // We got the PSSM, so get the other stuff that goes with it.

          if (verbosity > NORMAL_VERBOSE) {
            fprintf(stderr, "Retrieved PSSM for %s.\n",
            combine_string_list(alignment_species, "-"));
          }

          inverse_alignment_col_freqs = (ARRAY_T*) retrieve_object(
              alignment_species,
              alignment_col_freqs_list
          );
        }
      }
      /* ------------------------------------------------------------------- */

      free_tree(true, trimmed_tree);

      // Free list of species
      free_string_list(alignment_species);

      // Build a table for converting the index into the alignment
      // to an index into ungapped reference sequence.
      int* coord_conv_table = make_alignment_to_seq_table(
          current_ref_seq_index,
          trimmed_alignment);

      if (use_BLS) {
        n_scored += bls_score_sequence_in_alignment(
            alph,
            current_ref_seq_index,
            trimmed_alignment, 
            (scan_both_strands ? get_motif_st_id(motif) : get_motif_id(motif)),
            tree, // tree
            window_size,
            NULL, // background frequencies
            models, 
            pssm,
            (scan_both_strands ? get_motif_st_id(inverse_motif) : get_motif_id(motif)),
            inverse_models,
            inverse_pssm,
            coord_conv_table,
            gap_support, 
            gap_cost,
            output_pthresh,
            bls_distance,
            scanned_seq
         );

      } else {
        n_scored += score_sequence_in_alignment(
            alph,
            current_ref_seq_index,
            trimmed_alignment, 
            (scan_both_strands ? get_motif_st_id(motif) : get_motif_id(motif)),
            NULL, // tree
            window_size,
            models, 
            NULL, // background frequencies
            pssm,
            coord_conv_table,
            gap_support, 
            gap_cost,
            output_pthresh,
            scanned_seq
         );
      }
      
      free_alignment(trimmed_alignment);
    } // score current alignment

    // Compute q-values, if requested.
    // Doesn't work for the BLS scores
    set_pattern_is_complete(pattern);
    if (compute_qvalues && !use_BLS) {
      pattern_calculate_qvalues(pattern, NULL);
    }

    if (models != NULL) {
      int model_index;
      int num_models = get_motif_length(motif) + 1;
      for (model_index = 0; model_index < num_models; model_index++) {
        free_model(models[model_index]);
      }
      myfree(models);
    }

    free_object_list(pssm_list);
    free_object_list(inverse_pssm_list);
    
  } // score using current motif

  // Write out results
  if (text_only) {
    print_cisml_as_text(cisml);
  }
  else {
    if (! print_full_results(
      cisml,
      output_dirname,
      "motiph.xml",
      "motiph.html",
      "motiph.txt",
      "motiph.gff",
      allow_clobber,
      true
    )) report_err = true;
  }

  /**********************************************
   * Clean up.
   **********************************************/
  free_motifs(motifs);
  free_object_list(alignment_col_freqs_list);
  free_array(bg_freqs);
  free_tree(true, tree);
  free_string_list(filenames);
  free_string_list(selected_motifs);
  free_cisml(cisml);

  //FIXME: (tlb) this is useful for checking p-value accuracy
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Total number of positions scored: %d\n", n_scored);
  }

  return(report_err ? 1 : 0);
} // main
