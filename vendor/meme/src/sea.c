#define DEFINE_GLOBALS
#include "streme-utils.h"
#include "qvalue.h"
#include "red-black-tree.h"
#include "html-monolith.h"

#define DEFAULT_SEA_HOFRACT 0.1
#define DEFAULT_SEA_ETHRESH 10.0
#define DEFAULT_SEA_PQTHRESH 0.05 
#define DEFAULT_SEA_PSEUDOCOUNT 0.1

static char *DEFAULT_OUTPUT_DIRNAME = "sea_out";
static const char *HTML_FILENAME = "sea.html";
static const char *TSV_FILENAME = "sea.tsv";
static const char *SEQUENCES_TSV_FILENAME = "sequences.tsv";
static const char *TEMPLATE_FILENAME = "sea_template.html";

/*****************************************************************************
 * Print a usage message with an optional reason for failure.
 *****************************************************************************/
static void usage(char *format, ...) {
  va_list argp;
  char *usage =
    "\n"
    "Usage: sea [options]\n"
    "\n"
    "   Options:\n"
    "     --p <filename>            primary (positive) sequence file name (required)\n"
    "     --m <filename>            motif file name (required, may be repeated)\n"
/*    "     --objfun de|cd            objective function (<objfun>)\n"
    "                                 de : Differential Enrichment\n"
    "                                 cd : Central Distance\n"
    "                               default: de\n"
*/
    "     --n <filename>            control (negative) sequence file name;\n"
    "                               defaults: if --n is not given, then SEA\n"
    "                               creates one by shuffling the primary sequences\n"
/*
    "                               creates control sequences as follows:\n"
    "                                 <objfun> = de, shuffle primary sequences\n"
    "                                          = cd, no control sequences allowed\n"
*/
    "     --order <m>               estimate an m-order background model\n"
    "                               and use an m-order shuffle if creating\n"
    "                               control sequences from primary sequences;\n"
    "                               default: %d (DNA), %d (RNA), %d (Protein), %d (custom)\n"
    "     --bfile <bfile>           use the background model contained in bfile instead\n"
    "                               of creating it from the control sequences;\n"
    "                               default: see --order\n"
    "     --o <output_dir>          output directory; default: '%s'\n"
    "     --oc <output_dir>         allow overwriting; default: '%s'\n"
    "     --text                    output text only; overrides --o and --oc;\n"
    "                               default: create text, HTML and TSV files in <output_dir>\n"
    "     --hofract <hofract>       fraction of sequences in hold-out set; default: %g\n"
    "     --seed <seed>             random seed for shuffling sequences;\n"
    "                               default: %d\n"
    "     --xalph <alph_file>       motifs will be converted to this custom alphabet\n"
    "                               converts to uppercase unless both cases in core\n"
    "     --motif-pseudo <pc>       pseudocount for creating PWMs from motifs;\n"
    "                               default: %g\n"
    "     --inc <pattern>           name pattern to select as motif; may be\n"
    "                               repeated; default: all motifs are used\n"
    "     --exc <pattern>           name pattern to exclude as motif; may be\n"
    "                               repeated; default: all motifs are used\n"
    "     --thresh <thresh>         significance threshold for reporting enriched motifs;\n"
    "                               default: E-value= %g (%g if --qvalue or --pvalue given)\n"
    "     --qvalue                  use q-value significance threshold; default: E-value\n"
    "     --pvalue                  use p-value significance threshold; default: E-value\n"
    "     --align left|center|right align sequences left/center/right for site\n"
    "                               positional distribution plots; default: center\n"
    "     --noseqs                  do not output matching sequences TSV file;\n"
    "                               default: output matching sequences\n"
    "     --desc <desc>             include this description text in HTML\n"
    "     --dfile <dfile>           include contents of this description file in HTML\n"
    "     --version                 print the program version and exit\n"
    "     --verbosity 1|2|3|4|5     level of diagnostic output (default: %d)\n"
    "\n"
    ;
  if (format) {
    fprintf(stderr, "\n");
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, usage, DEFAULT_DNA_ORDER, DEFAULT_RNA_ORDER, DEFAULT_PROT_ORDER, DEFAULT_CUSTOM_ORDER,
    DEFAULT_OUTPUT_DIRNAME, DEFAULT_OUTPUT_DIRNAME, DEFAULT_SEA_HOFRACT, DEFAULT_SEED,
    DEFAULT_SEA_PSEUDOCOUNT, DEFAULT_SEA_ETHRESH, DEFAULT_SEA_PQTHRESH,
    DEFAULT_VERBOSITY
  );
  fflush(stderr);
  exit(EXIT_FAILURE);
} // usage

// An easy way of assigning a number to each option.
// Be careful that there are not more than 62 options as otherwise the value
// of '?' will be used.
enum Opts {
  OPT_M, OPT_P, OPT_N, OPT_ORDER, OPT_O, OPT_OC, OPT_TEXT, //OPT_OBJFUN,
  OPT_XALPH, OPT_BFILE, OPT_PSEUDOCOUNT, OPT_INC, OPT_EXC,
  OPT_HOFRACT, OPT_SEED, OPT_THRESH, OPT_QVALUE, OPT_PVALUE, OPT_DISTR, OPT_ALIGN,
  OPT_NOSEQS, OPT_DESC, OPT_DFILE, OPT_VERBOSITY, OPT_VERSION
};

/***********************************************************************
 Process command line options and return the command line string.
 ***********************************************************************/
static char *process_command_line(int argc, char* argv[], STREME_OPTIONS_T *options) {
  int i;
  struct option streme_options[] = {
    {"m", required_argument, NULL, OPT_M},
    {"p", required_argument, NULL, OPT_P},
    {"n", required_argument, NULL, OPT_N},
    {"order", required_argument, NULL, OPT_ORDER},
    {"o", required_argument, NULL, OPT_O},
    {"oc", required_argument, NULL, OPT_OC},
    {"text", no_argument, NULL, OPT_TEXT},
    //{"objfun", required_argument, NULL, OPT_OBJFUN},
    {"xalph", required_argument, NULL, OPT_XALPH},
    {"hofract", required_argument, NULL, OPT_HOFRACT},
    {"seed", required_argument, NULL, OPT_SEED},
    {"bfile", required_argument, NULL, OPT_BFILE},
    {"motif-pseudo", required_argument, NULL, OPT_PSEUDOCOUNT},
    {"inc", required_argument, NULL, OPT_INC},
    {"motif", required_argument, NULL, OPT_INC},
    {"exc", required_argument, NULL, OPT_EXC},
    {"thresh", required_argument, NULL, OPT_THRESH},
    {"qvalue", no_argument, NULL, OPT_QVALUE},
    {"pvalue", no_argument, NULL, OPT_PVALUE},
    {"align", required_argument, NULL, OPT_ALIGN},
    {"noseqs", no_argument, NULL, OPT_NOSEQS},
    {"desc", required_argument, NULL, OPT_DESC},
    {"dfile", required_argument, NULL, OPT_DFILE},
    {"verbosity", required_argument, NULL, OPT_VERBOSITY},
    {"version", no_argument, NULL, OPT_VERSION},
    {NULL, 0, NULL, 0} //boundary indicator
  };

  // Set options to defaults.
  memset(options, 0, sizeof(STREME_OPTIONS_T));
  options->output_dirname = DEFAULT_OUTPUT_DIRNAME;
  options->allow_clobber = true;
  options->text_only = false;
  options->posfile = NULL;
  options->negfile = NULL;
  //options->objfun = DEFAULT_OBJFUN;
  options->order = -1;			// flags no option given
  options->alphabet_type = DEFAULT_ALPHABET_TYPE;
  options->alph = NULL;
  options->alph_file = NULL;
  options->minwidth = DEFAULT_MINWIDTH;
  options->maxwidth = DEFAULT_MAXWIDTH;
  options->hofract = DEFAULT_SEA_HOFRACT;
  options->seed = DEFAULT_SEED;
  options->motif_sources = arraylst_create();
  options->bfile = NULL;
  options->pseudocount = DEFAULT_SEA_PSEUDOCOUNT;
  options->include_patterns = arraylst_create();
  options->exclude_patterns = arraylst_create();
  options->thresh = DEFAULT_SEA_ETHRESH;
  options->thresh_type = EVALUE;
  options->align_txt = NULL;
  options->noseqs = false;
  options->description = NULL;
  options->dfile = NULL;

  // Process arguments.
  int longindex;
  while (1) {
    //int opt = getopt_long_only(argc, argv, "", streme_options, NULL);
    int opt = getopt_long_only(argc, argv, "", streme_options, &longindex);
    if (opt == -1) break;
    switch (opt) {
      case OPT_O: // Set output directory with no clobber
        options->output_dirname = optarg;
        options->allow_clobber = false;
        break;
      case OPT_OC: // Set output directory with clobber
        options->output_dirname = optarg;
        options->allow_clobber = true;
        break;
      case OPT_TEXT: // Output TSV only to standard output.
        options->text_only = true;
        break;
      case OPT_M:
        arraylst_add(optarg, options->motif_sources);
        break;
      case OPT_P:
        options->posfile = optarg;
        break;
      case OPT_N:
        options->negfile = optarg;
        break;
#ifdef EXP
      case OPT_OBJFUN:
        if (! strcmp(optarg, "de")) {
          options->objfun = DE;
        } else if (! strcmp(optarg, "cd")) {
          options->objfun = CD;
        } else {
          usage("Unknown value for --objfun (%s)\n", optarg);
        }
        break;
#endif
      case OPT_ORDER:
        options->order = atoi(optarg);
        break;
      case OPT_XALPH:
        options->alphabet_type = Custom;
        options->alph_file = optarg;
        break;
      case OPT_HOFRACT:
        options->hofract = atof(optarg);
        break;
      case OPT_SEED:
        options->seed = atoi(optarg);
        break;
      case OPT_BFILE:
        options->bfile = optarg;
        break;
      case OPT_PSEUDOCOUNT:
	options->pseudocount = atof(optarg);
        break;
      case OPT_INC:
        arraylst_add(optarg, options->include_patterns);
        break;
      case OPT_EXC:
        arraylst_add(optarg, options->exclude_patterns);
        break;
      case OPT_THRESH:
        options->thresh = atof(optarg);
        break;
      case OPT_QVALUE:
        options->thresh_type = QVALUE;
        break;
      case OPT_PVALUE:
        options->thresh_type = PVALUE;
        break;
      case OPT_ALIGN:
        options->align_txt = optarg;
        break;
      case OPT_NOSEQS:
        options->noseqs = true;
        break;
      case OPT_DESC:
        options->description = optarg;
        break;
      case OPT_DFILE:
        options->dfile = optarg;
        break;
      case OPT_VERBOSITY:
        verbosity = atoi(optarg);
        break;
      case OPT_VERSION:
        fprintf(stdout, VERSION "\n");
        exit(EXIT_SUCCESS);
        break;
      case '?':
        usage(NULL);
        break;
      default: // just in case we forget to handle a option
        die("Unhandled option %d", opt);
    }
  }

  // Check for primary sequences.
  if (options->posfile == NULL) {
    usage("You must supply a FASTA file with the primary sequences.");
  }

  // Check objective function
  if (options->objfun == CD && options->negfile != NULL) {
    DEBUG_FMT(QUIET_VERBOSE, "# Warning: ignoring control sequence file (%s) with --objfun cd.\n", options->negfile);
    options->negfile = NULL;
  }

  // Set the negfile = posfile if it is null and we need it.
  if (options->objfun != CD && options->negfile == NULL) {
    options->negfile = options->posfile;
  }

  // Check for motif files.
  if (arraylst_size(options->motif_sources) == 0) {
    usage("You must provide at least one motif file using --m.");
  }

  // Check that order is not too large.
  if (options->order > MAX_ORDER) {
    usage("The value of --order must be in the range [0,..,%d].  order = %d", MAX_ORDER, options->order);
  }

  // Check that hofract is not too large.
  if (options->hofract > 0.5) {
    usage("The value of --hofract must be in the range [0,..,0.5].  hofract = %d", options->hofract);
  }

  // Check the significance threshold.
  if (options->thresh_type == QVALUE || options->thresh_type == PVALUE) {
    if (options->thresh == DEFAULT_SEA_ETHRESH) options->thresh = DEFAULT_SEA_PQTHRESH;
    if (options->thresh <= 0 || options->thresh > 1) {
      usage("The value of --thresh must be in the range (0,..1] when you use --%cvalue. thresh = %g\n", 
        (options->thresh_type == QVALUE ? 'q' : 'p'), options->thresh);
    }
  } else {
    if (options->thresh <= 0) {
      usage("The value of --thresh must greater than 0. thresh = %g\n", options->thresh);
    }
  }

  // Set the alignment of sequences for site positional distribution plots.
  if (options->align_txt) {
    if (!strcmp(options->align_txt, "left")) {
      options->align = LEFT;
    } else if (!strcmp(options->align_txt, "center")) {
      options->align = CENTER;
    } else if (!strcmp(options->align_txt, "right")) {
      options->align = RIGHT;
    } else {
      usage("The value of --align must be 'left', 'center' or 'right'. align = %s\n", options->align_txt);
    }
  } else {
    options->align_txt = "center";
    options->align = CENTER;
  }
  
  // Get the contents of the description file.
  if (options->dfile) options->description = get_description_file(options->dfile);

  // Check the verbosity level.
  if (verbosity < 1 || verbosity > 5) {
    usage("The verbosity level must be in the range [1, ..., 5]. verbosity = %d", verbosity);
  }

  // make enough space for all the command line options, with one space between each
  // add on argc to allow one char per word for separating space + terminal '\0'
  size_t line_length = argc;
  for (i = 0; i < argc; i++) line_length += strlen(i == 0 ? basename(argv[0]) : argv[i]);
  char *commandline = (char*) malloc(sizeof(char)*line_length);
  int nextpos = 0;
  for (i = 0; i < argc; i++) {
    // been here before? put in a space before adding the next word
    if (nextpos) {
      commandline[nextpos] = ' ';
      nextpos++;
    }
    char *nextword = (i == 0) ? basename(argv[0]) : argv[i];
    strcpy(&commandline[nextpos], nextword);
    nextpos += strlen (nextword);
  }
  return(commandline);
} // process_command_line

//
// Free the storage and check for leaks.
//
void cleanup(
  STREME_OPTIONS_T *options,		// the program options
  Multiseq *test_multiseq
) {
  // Free the storage used.
  DEBUG_MSG(NORMAL_VERBOSE, "# Freeing storage...\n");

  // Free the things saved in options.
  alph_release(options->alph);
  arraylst_destroy(NULL, options->include_patterns);
  arraylst_destroy(NULL, options->exclude_patterns);
  arraylst_destroy(NULL, options->motif_sources);

  // free the sequences
  FREESPACE_TLB(test_multiseq->sequence);
  freemultiseq(test_multiseq);

  // Free the p-value cache.
  if (options->objfun == DE) {
    int i, j;
    for (i=0; i<options->cache_length; i++) {
      free(options->pv_cache[i]);
    }
    free(options->pv_cache);
  }

  // Check for leaks.
  mmcheckspaceleak();
} // cleanup

/*************************************************************************
 * Read in the motif databases.
 *************************************************************************/
ARRAYLST_T *sea_load_motifs_and_background(
  STREME_OPTIONS_T *options,
  bool separate_namespaces,     // keep motif DB namespaces separate
  ARRAY_T **back,               // OUT the background
  int *max_width                // OUT maximum motif width
)
{
  int i, num_motifs = 0;
  ARRAYLST_T *dbs = arraylst_create_sized(arraylst_size(options->motif_sources));
  RBTREE_T *motif_names = rbtree_create(rbtree_strcmp, NULL, NULL, NULL, NULL);

  // Set the alph "pseudo-option".
  bool xalph = (options->alph_file != NULL);
  options->alph = NULL;
  if (xalph) {
    options->alph = alph_load(options->alph_file, true);
    if (options->alph == NULL) exit(EXIT_FAILURE);
  }

  // Read all the alphabets and make sure they are the same.
  DEBUG_FMT(NORMAL_VERBOSE, "# Checking alphabets in %d motif files.\n", 
    arraylst_size(options->motif_sources));
  read_motif_alphabets(options->motif_sources, xalph, &(options->alph));

  ALPH_T *alph = options->alph;
  assert(alph != NULL);
  bool use_rc = alph_has_complement(alph);

  bool stdin_used = false;      // set if path is "-"
  *max_width = 0;               // maximum motif width
  char *save_bfile = options->bfile;	// in case it gets changed
  for (i = 0; i < arraylst_size(options->motif_sources); i++) {

    // Get the name of the next motif db.
    char *motif_source = (char*) arraylst_get(i, options->motif_sources);
    DEBUG_FMT(NORMAL_VERBOSE, "# Loading motifs from file '%s'\n", motif_source);

    // Load the motifs from this file.
    MOTIF_DB_T *db = read_motifs_and_background(
      i,                        // id of DB file
      motif_source,             // motif file name (or special word)
      "Query motifs",           // type of database for error messages
      NULL,                     // get one motif by name
      NULL,                     // get one motif by index
      options->include_patterns,// get set of motifs by name (or NULL)
      options->exclude_patterns,// exclude this set of motifs by name (or NULL)
      true,                     // allow motifs with zero probability entries
      false,                    // don't create RC copies, appended
      options->pseudocount,     // multiply times background model
      false,                    // set_trim
      0,                        // trim_bit_threshold
      &(options->bfile),        // background file; may be changed
      true,                     // make bg symmetrical if alph complementable
      back,                     // will be set if id==0
      options->posfile,         // sequence file name
      alph,                     // sequence alphabet
      xalph,                    // set motif conversion alphabet
      false,                    // don't remove extension from name
      false,                    // don't remove ".meme" extension from name
      false,                    // don't replace underscores in name
      &stdin_used               // IN/OUT check and set if path is "-"
    );

    // Get maximum motif width and check for motif name uniqueness across all DBs.
    int warn_type = NO_WARNING;
    int i;
    for (i=0; i<arraylst_size(db->motifs); i++) {
      MOTIF_T *motif = arraylst_get(i, db->motifs);
      if (! separate_namespaces) {
        bool created;
        RBNODE_T *node = rbtree_lookup(motif_names, get_motif_id(motif), true, &created);
        if (!created) {
          clump_motif_db_warning(&warn_type, DUPLICATE_WARNING, "# Warning: The following "
            "duplicate motifs in '%s' were excluded:\n  ", motif_source, get_motif_id(motif));
          destroy_motif(motif);
          db->loaded--;
          db->excluded++;
          continue;
        }
      } // namespaces
      // Get width of widest motif.
      int w = get_motif_length(motif);
      if (w > *max_width) *max_width = w;
    }

    // Add DBs to the list of DBs
    arraylst_add(db, dbs);

    // Number of motifs found.
    num_motifs += arraylst_size(db->motifs);

    if (warn_type && verbosity >= NORMAL_VERBOSE) fprintf(stderr, "\n");
  } // motif_files

  // Check that we found suitable motifs.
  if (num_motifs == 0) die("No acceptable motifs found.");

  // Delete the background unless --bfile was a special keyword (--uniform--, --motif--, etc.).
  bool special_keyword = strncmp(options->bfile, "--", 2) == 0;
  options->bfile = save_bfile;	// restore the original bfile 
  if (options->bfile == NULL || !special_keyword) {
    // bfile was not a special keyword.
    free_array(*back);
    *back = NULL;
  } else if (options->bfile && special_keyword) {
    // special keyword backgrounds are all 0-order
    if (options->order != 0) {
      DEBUG_FMT(NORMAL_VERBOSE, "# The specified --bfile (%s) is 0-order; setting --order to 0.\n", options->bfile);
      options->order = 0;
    }
  }

  // Cleanup
  rbtree_destroy(motif_names);

  return(dbs);
} // sea_load_motifs_and_background

/*************************************************************************
 * Convert motif to STREME model format.
 * Only sets the STREME PSPM matrix from the MEME motif frequency matrix.
 *************************************************************************/
Model *motif2streme(
  MOTIF_T *motif
) {
  int i, j;
  Model *model = (Model *)malloc(sizeof(Model));
  model->alph = get_motif_alph(motif);
  model->alen = alph_size_core(model->alph);
  model->width = get_motif_length(motif);
  model->npassing = 0;
  // PSPM
  for (i=0; i<model->alen; i++)
    for (j=0; j<model->width; j++)
      model->probs[i][j] = get_matrix_cell(j, i, get_motif_freqs(motif));
  model->score_threshold = 0;
  return(model);
} // motif2streme

//
// Compare function for sorting Model pointers
// with qsort in increasing order of test_log_pvalue.
// Return <0 >0
// if the second log_pvalue is <, > than the first.
// Break ties using the motif id, then db index.
//
static int compare_model_test_pvalue(
  const void *v1,
  const void *v2
)
{
  const Model *m1 = *(const Model **) v1;
  const Model *m2 = *(const Model **) v2;

  if (m1->test_log_pvalue > m2->test_log_pvalue) {
    return(+1);
  } else if (m1->test_log_pvalue < m2->test_log_pvalue) {
    return(-1);
  } else if (strcmp(m1->motif->id, m2->motif->id) > 1) {
    return(+1);
  } else if (strcmp(m1->motif->id, m2->motif->id) < 1) {
    return(-1);
  } else {
    return(m1->db_idx - m2->db_idx);
  }
} // compare_model_test_pvalue

//
// Compute q-values.
//
static void convert_sea_p_to_q(
  int num_motifs,		// the number of scored motifs
  Model **models 		// the scored motifs
) {
  ARRAY_T *pvalues;
  int i;

  // Extract all of the p-values.
  pvalues = allocate_array(num_motifs);
  for (i = 0; i < num_motifs; i++) {
    double pvalue = exp(models[i]->test_log_pvalue);
    set_array_item(i, pvalue, pvalues);
  }

  // Convert p-values to q-values in place.
  compute_qvalues(
    false, // Don't stop with FDR.
    true, // Estimate pi-zero.
    NULL, // Don't store pi-zero in a file.
    NUM_BOOTSTRAPS,
    NUM_BOOTSTRAP_SAMPLES,
    NUM_LAMBDA,
    MAX_LAMBDA,
    num_motifs,
    pvalues,
    NULL // No sampled p-values provided
  );

  // Set log q-values.
  // Replace with log E-value when q-value is 0 (since qvalue() works in non-log space).
  for (i = 0; i < num_motifs; i++) {
    Model *model = models[i];
    double qvalue = get_array_item(i, pvalues);
    model->test_log_qvalue = qvalue ? log(qvalue) : model->test_log_evalue;
  }

  // clean up
  free_array(pvalues);
} // convert_sea_p_to_q

//
// Output the site distribution and site count histogram in the positive sequences.
//
void output_json_site_distr(
  JSONWR_T *json_output,        // JSON output file pointer
  Model *model,                 // The model
  Multiseq *multiseq            // The sequences
) {
  int i;
  int w = model->width;         // width of motif
  int maxlen = multiseq->max_poslen;    // maximum length of positive sequences

  // Write the array of site counts.
  jsonwr_lng_prop(json_output, "total_sites", model->total_sites);
  jsonwr_property(json_output, "site_distr");
  jsonwr_start_array_value(json_output);
  for (i=0; i<maxlen-w; i++) {
    jsonwr_lng_value(json_output, model->site_distr[i]);
  }
  jsonwr_end_array_value(json_output);

  // Write the histogram of site counts.
  jsonwr_lng_prop(json_output, "max_sites", model->max_sites);
  jsonwr_property(json_output, "site_hist");
  jsonwr_start_array_value(json_output);
  for (i=0; i<=model->max_sites; i++) {
    jsonwr_lng_value(json_output, model->site_hist[i]);
  }
  jsonwr_end_array_value(json_output);
} // output_json_site_distr

//
// Output the SEQUENCES_TSV results.
//
static void print_passing_sequences_tsv(
  STREME_OPTIONS_T *options,    // SEA options
  ARRAYLST_T *dbs,              // the motif DBs
  Model *model,                 // the model including its passing sequences
  int imotif                    // motif number
) {
  int i;
  FILE *outfile = options->sequences_tsv_output;

  // Print the SEQUENCES_TSV header for first motif.
  if (imotif == 0) {
    fprintf(outfile,
      "%s\t"
      "%s\t"
      "%s\t"
      "%s\t"
      "%s\t"
      "%s\t"
      "%s\n",
      "motif_DB", "motif_ID", "motif_ALT_ID", "seq_ID", "seq_Score", "seq_Class", "is_holdout?"
    );
  } else {
    fprintf(outfile, "#\n");
  }

  // Print the sorted passing sequences.
  MOTIF_DB_T *db = arraylst_get(model->db_idx, dbs);
  MOTIF_T *motif = model->motif;
  Passing_seq *sequences = model->passing_seqs;
  int npassing = model->npassing;
  for (i=0; i<npassing; i++) {
    fprintf(outfile,
      "%s\t"    // Motif DB
      "%s\t"    // Motif ID
      "%s\t"    // Motif Alternate ID
      "%s\t"    // Sequence ID
      "%.2f\t"  // Score
      "%s\t"    // tf/fp
      "%d\n",   // holdout set?
      db->source, get_motif_id(motif), get_motif_id2(motif),
      sequences[i].desc, sequences[i].score, (sequences[i].is_tp ? "tp" : "fp"),
      sequences[i].dbno
    );
  } // passing sequence

} // print_passing_sequences_tsv

//
// Output the HTML, TSV and SEQUENCES_TSV results.
//
static void output_results(
  int argc,
  char **argv,
  STREME_OPTIONS_T *options,
  Multiseq *train_multiseq,	// training sequences
  Multiseq *test_multiseq,	// test sequences
  int num_motifs,		// the number of scored motifs
  Model **models,		// the scored motifs
  ARRAYLST_T *dbs		// the motif DBs
) {
  double *bg_freqs = test_multiseq->background;
  int pos_num_seqs = test_multiseq->npos + (train_multiseq ? train_multiseq->npos : 0);
  int neg_num_seqs = test_multiseq->nneg + (train_multiseq ? train_multiseq->nneg : 0);
  int max_pos_len = test_multiseq->max_poslen; 	// maximum length of test positive sequences
  bool use_binomial = test_multiseq->use_binomial; // Binomial test was used rather than Fisher
  bool is_holdout = (train_multiseq != NULL);	// there was a hold-out set
  bool do_rc = alph_has_complement(options->alph);
  HTMLWR_T *html_output=NULL;
  JSONWR_T *json_output=NULL;
  MOTIF_DB_T* db;
  int i, j;

  // Print the TSV header.
  fprintf(options->text_output, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\%s\t%s\t", 
    "RANK", "DB", "ID", "ALT_ID", "CONSENSUS", "TP", "TP%", "FP", "FP%", "ENR_RATIO", "SCORE_THR");
  fprintf(options->text_output, "%s\t%s\t%s\t%s\t%s\t%s\n",
     "PVALUE", "LOG_PVALUE", "EVALUE", "LOG_EVALUE", "QVALUE", "LOG_QVALUE");

  // setup the HTML output writer
  if (!options->text_only && (html_output = htmlwr_create(get_meme_data_dir(), TEMPLATE_FILENAME, false))) {
    htmlwr_set_dest_name(html_output, options->output_dirname, HTML_FILENAME);
    htmlwr_replace(html_output, "sea_data.js", "data");
    json_output = htmlwr_output(html_output);
    if (json_output == NULL) die("Template does not contain data section.\n");
  } else if (!options->text_only) {
    DEBUG_MSG(QUIET_VERBOSE, "# Failed to open HTML template file; no HTML output created.\n");
  }

  // Print the JSON initial section
  if (json_output) {
    // write out some information
    jsonwr_str_prop(json_output, "version", VERSION);
    jsonwr_str_prop(json_output, "revision", REVISION);
    jsonwr_str_prop(json_output, "release", ARCHIVE_DATE);
    jsonwr_str_prop(json_output, "program", "SEA");
    if (options->description) jsonwr_str_prop(json_output, "description", options->description);
    jsonwr_args_prop(json_output, "cmd", argc, argv);

    // options
    jsonwr_property(json_output, "options");
    jsonwr_start_object_value(json_output);
    jsonwr_str_prop(json_output, "objfun", options->objfun == DE ? "Differential Enrichment" : "Central Distance");
    jsonwr_str_prop(json_output, "test", options->objfun == DE ? 
      (use_binomial ? "Binomial Test" : "Fisher Exact Test")
      : "Central Distance");
    jsonwr_lng_prop(json_output, "order", options->order);
    jsonwr_lng_prop(json_output, "seed", options->seed);
    jsonwr_dbl_prop(json_output, "pseudocount", options->pseudocount);
    jsonwr_str_prop(json_output, "thresh_type", options->thresh_type == 
      EVALUE ?  "evalue" : (options->thresh_type == QVALUE ? "qvalue" : "pvalue")
    );
    jsonwr_dbl_prop(json_output, "thresh", options->thresh);
    jsonwr_dbl_prop(json_output, "hofract", options->hofract);
    jsonwr_bool_prop(json_output, "is_holdout", is_holdout);
    jsonwr_str_prop(json_output, "strands", (do_rc) ? "both" : "given");
    jsonwr_str_prop(json_output, "align", options->align_txt);
    jsonwr_str_prop(json_output, "noseqs", options->noseqs ? "true" : "false");
    jsonwr_end_object_value(json_output);

    // output the motif dbs
    jsonwr_property(json_output, "motif_dbs");
    jsonwr_start_array_value(json_output);
    for (i = 0; i < arraylst_size(options->motif_sources); i++) {
      db = arraylst_get(i, dbs);
      jsonwr_start_object_value(json_output);
      jsonwr_str_prop(json_output, "source", db->source);
      jsonwr_lng_prop(json_output, "count", arraylst_size(db->motifs));
      jsonwr_end_object_value(json_output);
    }
    jsonwr_end_array_value(json_output);

    // Write alphabet to json data.
    jsonwr_property(json_output, "alphabet");
    alph_print_json(options->alph, json_output);

    // Write the background model to json data.
    jsonwr_property(json_output, "background");
    jsonwr_start_object_value(json_output);
    jsonwr_str_prop(json_output, "source", options->bfile ? options->bfile : "--control--");
    if (options->bfile) jsonwr_str_prop(json_output, "file", options->negfile);
    jsonwr_property(json_output, "frequencies");
    jsonwr_start_array_value(json_output);
    for (i = 0; i < alph_size_core(options->alph); i++) {
      jsonwr_dbl_value(json_output, bg_freqs[i]);
    }
    jsonwr_end_array_value(json_output);
    jsonwr_end_object_value(json_output);

    // Output the sequence databases
    jsonwr_property(json_output, "sequence_db");
    jsonwr_start_object_value(json_output);
    jsonwr_str_prop(json_output, "source", options->posfile);
    jsonwr_lng_prop(json_output, "count", pos_num_seqs);
    jsonwr_lng_prop(json_output, "maxlen", max_pos_len);
    jsonwr_lng_prop(json_output, "holdout", (is_holdout ? train_multiseq->npos : 0));
    jsonwr_end_object_value(json_output);
    if (options->negfile) {
      jsonwr_property(json_output, "control_db");
      jsonwr_start_object_value(json_output);
      if (options->objfun != CD && options->negfile != options->posfile) {
        jsonwr_str_prop(json_output, "source", options->negfile);
      } else {
        jsonwr_str_prop(json_output, "from", "shuffled");
      }
      jsonwr_lng_prop(json_output, "count", neg_num_seqs);
      jsonwr_lng_prop(json_output, "holdout", (is_holdout ? train_multiseq->nneg : 0));
      jsonwr_end_object_value(json_output);
    }
  } // JSON initial section

  // Sort the results by test p-value.
  qsort(models, num_motifs, sizeof(Model *), compare_model_test_pvalue);

  // Compute q-values.
  convert_sea_p_to_q(num_motifs, models);

  // Print the significant results to JSON, TSV and SEQUENCES_TSV.
  double log_thresh = log(options->thresh);
  if (json_output) {
    jsonwr_property(json_output, "motifs");
    jsonwr_start_array_value(json_output);
  }
  int imotif;
  for (imotif=0; imotif<num_motifs; imotif++) {
    Model *model = models[imotif];
    MOTIF_T *motif = model->motif;
    MOTIF_DB_T *db = arraylst_get(model->db_idx, dbs);
    // Done if signficance threshold exceeded.
    double log_significance = options->thresh_type == PVALUE ? model->test_log_pvalue : 
      (options->thresh_type == QVALUE ? model->test_log_qvalue : model->test_log_evalue);
    if (log_significance > log_thresh) break;
    // Get the string versions of the p-values etc.
    int prec = 2;
    char *pvalue = print_log_value(NULL, model->test_log_pvalue, prec);
    char *evalue = print_log_value(NULL, model->test_log_evalue, prec);
    char *qvalue = print_log_value(NULL, model->test_log_qvalue, prec);
    double tp_percent = test_multiseq->npos ? 100 * ((double) model->test_pos_count/test_multiseq->npos) : 0;
    double fp_percent = test_multiseq->nneg ? 100 * ((double) model->test_neg_count/test_multiseq->nneg) : 0;
    // Output the motif to TSV.
    fprintf(options->text_output, 
      "%d\t"	// rank
      "%s\t"	// source
      "%s\t"	// ID
      "%s\t"	// ID2
      "%s\t"	// consensus
      "%d\t"	// tp
      "%.2f\t"	// tp%
      "%d\t"	// fp
      "%.2f\t"	// fp%
      "%.3g\t"	// enrichment ratio
      "%.2g\t"	// score threshold
      "%s\t"	// p-value
      "%.2f\t"	// log(p-value)
      "%s\t"	// E-value
      "%.2f\t"	// log(E-value)
      "%s\t"	// q-value
      "%.2f\n",	// log(q-value)
      imotif+1, db->source, get_motif_id(motif), get_motif_id2(motif), get_motif_consensus(motif),
      model->test_pos_count, tp_percent, model->test_neg_count, fp_percent, model->test_ratio, 
      model->score_threshold, pvalue, model->test_log_pvalue, evalue, model->test_log_evalue,
      qvalue, model->test_log_qvalue);

    // Output the motif to JSON.
    if (json_output) {
      jsonwr_start_object_value(json_output);
      jsonwr_lng_prop(json_output, "db", model->db_idx);
      jsonwr_str_prop(json_output, "id", get_motif_id(motif));
      if (get_motif_id2(motif)[0] != '\0')
	jsonwr_str_prop(json_output, "alt", get_motif_id2(motif));
      jsonwr_lng_prop(json_output, "len", get_motif_length(motif));
      jsonwr_str_prop(json_output, "consensus", get_single_letter_consensus(model));
      jsonwr_dbl_prop(json_output, "motif_evalue", get_motif_evalue(motif));
      jsonwr_dbl_prop(json_output, "motif_nsites", get_motif_nsites(motif));
      if (has_motif_url(motif)) jsonwr_str_prop(json_output, "url", get_motif_url(motif));
      int alen = alph_size_core(get_motif_alph(motif));
      int len = get_motif_length(motif);
      MATRIX_T *freqs = get_motif_freqs(motif);
      jsonwr_property(json_output, "pwm");
      jsonwr_start_array_value(json_output);
      for (i = 0; i < len; i++) {
	ARRAY_T *row = get_matrix_row(i, freqs);
	jsonwr_start_array_value(json_output);
	for (j = 0; j < alen; j++) {
	  jsonwr_dbl_value(json_output, get_array_item(j, row));
	}
	jsonwr_end_array_value(json_output);
      }
      jsonwr_end_array_value(json_output);
      jsonwr_str_prop(json_output, "pvalue", pvalue);
      jsonwr_dbl_prop(json_output, "log_pvalue", model->test_log_pvalue);
      jsonwr_str_prop(json_output, "evalue", evalue);
      jsonwr_dbl_prop(json_output, "log_evalue", model->test_log_evalue);
      jsonwr_str_prop(json_output, "qvalue", qvalue);
      jsonwr_dbl_prop(json_output, "log_qvalue", model->test_log_qvalue);
      jsonwr_lng_prop(json_output, "tp", model->test_pos_count);
      jsonwr_lng_prop(json_output, "fp", model->test_neg_count);
      jsonwr_dbl_prop(json_output, "enr_ratio", model->test_ratio);
      jsonwr_dbl_prop(json_output, "score_thresh", model->score_threshold);
      output_json_site_distr(json_output, model, test_multiseq);
      jsonwr_end_object_value(json_output);
      myfree(pvalue);
      myfree(evalue);
    } // json_out

    // Output the motif's passing sequences to SEQUENCES_TSV.
    if (! options->text_only && options->noseqs == false) {
      // Initialize the passing sequences from the main sequence set.
      get_passing_sequences(options, model, test_multiseq, false, (train_multiseq ? false : true));
      // Add the passing sequences from the hold-out sequence set.
      if (train_multiseq) {
	// Free the matches in the model.
	free(model->matches);
	// Score the second set of sequences getting just the ZOOPS sites (ignores score threshold).
	score_model_pssm(options, train_multiseq, model, false, false, ZOOPS, false, false);
	get_passing_sequences(options, model, train_multiseq, true, true);
      }
      print_passing_sequences_tsv(options, dbs, model, imotif);
    }
    
    // Free the matches in the model.
    free(model->matches);
    model->matches = NULL;
    model->nmatches = 0;
  } // motif

  // Finish the JSON output.
  if (json_output) jsonwr_end_array_value(json_output);
  if (json_output) jsonwr_end_object_value(json_output);

  // Finish the TSV output.
  char *version_message = "# SEA (Simple Enrichment Analysis): Version " VERSION " compiled on " __DATE__ " at " __TIME__ "\n";
  fprintf(options->text_output, "\n%s", version_message);
  fprintf(options->text_output, "# The format of this file is described at %s/%s.\n", SITE_URL, "doc/sea-output-format.html");
  fprintf(options->text_output, "# %s\n", options->commandline);

  // Finish the HTML output.
  if (html_output) {
    if (htmlwr_output(html_output) != NULL) {
      die("Expected only one replacement variable in template.\n");
    }
    htmlwr_destroy(html_output);
  }
} // output_results

/*
  Compute and print STREME objective function scores for each motif in the input.
*/
int main(int argc, char *argv[])
{
  int i, db_i, motif_i, num_motifs;
  STREME_OPTIONS_T options;
  ARRAY_T *back = NULL;                 // the MEME-style background model
  ARRAYLST_T *dbs = NULL;               // the motif DBs
  Model **models = NULL;

  // command line processing
  options.commandline = process_command_line(argc, argv, &options);

  srand_mt(options.seed);		// random seed for shuffling sequences
  set_randfunc((randfunc_t) random_mt); // for ushuffle
  INITPOWER2TABLE();			// needed by get_seqno_list()

  // Initialize the p-value lookup cache if we are using objective function DE.
  // It is indexed by the pos_count + neg_count, and each sub-cache
  // is indexed by pos_count.
  options.cache_length = sqrt(2*MAX_CACHE_SIZE);
  if (options.objfun == DE) options.pv_cache = (double**) calloc(options.cache_length, sizeof(double*));

  // Set up for output.
  if (options.text_only) {
    options.text_output = stdout;
  } else {
    // Create the output directory.
    if (create_output_directory(options.output_dirname, options.allow_clobber, false)) {
      exit(EXIT_FAILURE);
    }
    // Create the path of the TSV output file and open it.
    char *path = make_path_to_file(options.output_dirname, TSV_FILENAME);
    if ((options.text_output = fopen(path, "w")) == NULL) {
      fprintf(stderr, "ERROR: Unable to open TSV output file '%s' for writing.\n", path);
      exit(EXIT_FAILURE);
    }
    // Create the path of the SEQUENCES_TSV output file and open it.
    path = make_path_to_file(options.output_dirname, SEQUENCES_TSV_FILENAME);
    if ((options.sequences_tsv_output = fopen(path, "w")) == NULL) {
      fprintf(stderr, "ERROR: Unable to open SEQUENCES TSV output file '%s' for writing.\n", path);
      exit(EXIT_FAILURE);
    }
  }

  // Load the motifs and the background model.
  int max_width;
  dbs = sea_load_motifs_and_background(
    &options,
    true,               // keep motif DB namespaces separate
    &back,              // OUT MEME-style background model
    &max_width          // OUT maximum motif width
  );

  // Initialize the STREME-style alphabet.
  Uint alen = initialize_st_alphabet(options.alph);
  // Set the alphabet type.
  const char *aname = alph_name(options.alph);
  if (! strcmp(aname, "DNA")) {
    options.alphabet_type = Dna;
    if (options.order == -1) options.order = DEFAULT_DNA_ORDER;
  } else if (! strcmp(aname, "RNA")) {
    options.alphabet_type = Rna;
    if (options.order == -1) options.order = DEFAULT_RNA_ORDER;
  } else if (! strcmp(aname, "Protein")) {
    options.alphabet_type = Protein;
    if (options.order == -1) options.order = DEFAULT_PROT_ORDER;
  } else {
    options.alphabet_type = Custom;
    options.order = DEFAULT_CUSTOM_ORDER;
  }
  bool do_rc = alph_has_complement(options.alph);
  DEBUG_FMT(NORMAL_VERBOSE, "# Alphabet: %s\n", aname);

  // Input the sequences, converting non-core characters to SEPARATOR.
  // Note that this also estimates the background (if NULL) from the control sequences.
  // Note that we are swapping the training and hold-out sets here relative to STREME.
  options.nmotifs = 1;	// suppress warning message
  Multiseq *train_multiseq = NULL;
  Multiseq *test_multiseq = read_pos_neg_seqs(&options, do_rc, false, false, true, back, &train_multiseq);
  if (train_multiseq) set_logcumback(&options, train_multiseq);
  set_logcumback(&options, test_multiseq);

  // Count the number of motifs.
  num_motifs = 0;
  for (db_i=0; db_i < arraylst_size(options.motif_sources); db_i++) {
    MOTIF_DB_T *db = arraylst_get(db_i, dbs);
    num_motifs += arraylst_size(db->motifs); // number of motifs in this DB
  }

  // Score each motif.
  models = (Model **) malloc(sizeof(Model *) * num_motifs);
  for (db_i=0, i=0; db_i < arraylst_size(options.motif_sources); db_i++) {
    MOTIF_DB_T *db = arraylst_get(db_i, dbs);
    int nmotifs = arraylst_size(db->motifs); 	// number of motifs in this DB
    for (motif_i = 0; motif_i < nmotifs; motif_i++) {
      MOTIF_T *motif = (MOTIF_T *) arraylst_get(motif_i, db->motifs);
      // Create the PSSM.
      Model *model = motif2streme(motif);
      model->motif = motif;
      models[i++] = model;
      // DB index.
      model->db_idx = db_i;
      // Score the motif on the training set (small) to get score threshold.
      // Don't use the cache.
      if (train_multiseq) score_model_pssm(&options, train_multiseq, model, true, false, NONE, false, false);
      // Score the motif on the hold-out set (big) using the computed threshold if there is one.
      // Use the cache for speed.  Save all passing matches.
      score_model_pssm(&options, test_multiseq, model, !train_multiseq, true, PASSING, true, true);
      // Adjust for multiple tests if the score threshold was optimized on the test set.
      if (!train_multiseq) model->test_log_pvalue = LOGEV(log(model->n_eff_tests), model->test_log_pvalue);
      // Compute the E-value.
      model->test_log_evalue = log(num_motifs) + model->test_log_pvalue;
      // Get the total number of sequences with sites and the site count distribution.
      get_site_distr(model, test_multiseq, options.align);
      // Progress counter.
      DEBUG_FMT(NORMAL_VERBOSE, "DB: %d MOTIF: %d/%d\r", db_i+1, motif_i+1, nmotifs);
    }
  }

  // Output the HTML and TSV results.
  output_results(argc, argv, &options, train_multiseq, test_multiseq, num_motifs, models, dbs);

  // Clean up.
  for (i=0; i<num_motifs; i++) free(models[i]);
  cleanup(&options, test_multiseq);

  return(EXIT_SUCCESS);
} // main
