#include <assert.h>
#include <errno.h>
#include <fnmatch.h>
#include <getopt.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <libgen.h> // for basename
#include <time.h>
#include <unistd.h>
#include "matrix.h" // includes array.h with special defines set so must be first
#include "alphabet.h"
#include "array-list.h"
#include "config.h"
#include "dir.h"
#include "fasta-io.h"
#include "html-monolith.h"
#include "io.h"
#include "linked-list.h"
#include "macros.h"
#include "motif-in.h"
#include "motif.h"
#include "motif-db.h"
#include "projrel.h"
#include "spamo-matches.h"
#include "spamo-output.h"
#include "spamo-scan.h"
#include "red-black-tree.h"
#include "utils.h"
#include "xml-out.h"
#include "xml-util.h"

// default file names
#define TEMPLATE_FILENAME "spamo_template.html"
#define DATA_FILENAME "spamo_data.js"
#define HTML_FILENAME "spamo.html"
#define TEXT_FILENAME "spamo.tsv"

// default verbosity
VERBOSE_T verbosity = NORMAL_VERBOSE;

typedef struct SPAMO_OPTIONS {
  // pseudo random number generator seed
  unsigned int prng_seed;
  // program start time
  time_t start;
  // output directory
  char *outdir;
  // should existing output be overwritten?
  bool clobber;
  // should EPS files be created?
  bool create_eps;
  // should PNG files be created?
  bool create_png;
  // should images be made for redundant motifs?
  bool create_for_redundant;
  // should the matches be output for motifs declared significant?
  bool dump_seqs;
  // if the same file is specified for the primary motif and the secondary 
  // motif should the primary be included?
  bool keep_primary;
  // pseudocount added to motifs loaded
  double pseudocount;
  // motif background for distributing pseudocounts and creating pssms
  char *bg_filename;
  // maximum distance from primary motif to secondary motif and minimum 
  // distance to edge for primary motif.
  int margin;
  // bin size for calculations (1 is highly recommended)
  int bin;
  // range for significance testing
  int test_range;
  // bit threshold for trimming motifs
  double trim_bit_threshold;
  // fraction of bases needed to be equal before they are considered redundant
  double sequence_similarity;
  // minimum score to be accepted by the scanning code
  double score_threshold;
  // use only the best match to the secondary motif if true
  bool use_best_secondary;
  // pvalue threshold for spacings to be considered significant.
  double pvalue_cutoff;
  // Evalue threshold for reporting secondary motifs
  double motif_evalue_cutoff;
  // the minimum overlap in the best site
  int overlap;
  // the fraction of the intersection as relative to the smaller set
  double joint;
  // name of primary motif
  RBTREE_T *primary_name;
  // index of primary motif
  RBTREE_T *primary_index;
  // wildcard patterns for motif names to include
  ARRAYLST_T *include_patterns;
  // wildcard patterns for motif names to exclude
  ARRAYLST_T *exclude_patterns;
  // the sequences to be scanned
  char *sequences_file;
  // the primary motif file
  char *primary_file;
  // the files for the secondary motifs
  ARRAYLST_T *secondary_files;
  // should the alphabets of the secondary motifs be converted?
  bool xalph;
  // only output the spamo.tsv file not .html or .xml
  bool text_only;
  // Only dump sequences significant matches.
  bool dump_sig_only;
  // Control sampling for speedup of sequence similarity check.
  double odds;
} SPAMO_OPTIONS_T;

/**************************************************************************
 * Frees memory allocated to store the options.
 **************************************************************************/
static void cleanup_options(SPAMO_OPTIONS_T *options) {
  arraylst_destroy(NULL, options->include_patterns);
  arraylst_destroy(NULL, options->exclude_patterns);
  arraylst_destroy(NULL, options->secondary_files);
  if (options->bg_filename) free(options->bg_filename);
  if (options->primary_index) rbtree_destroy(options->primary_index);
  if (options->primary_name) rbtree_destroy(options->primary_name);
}

/**************************************************************************
 * Prints a usage message and exits. 
 * If given an error message it prints that first and will exit with
 * return code of EXIT_FAILURE.
 **************************************************************************/
static void usage(SPAMO_OPTIONS_T *options, char *format, ...) {
  va_list argp;

  char *usage = 
    "Usage:\n" 
    "    spamo [options] <sequences> <primary motif> <secondary motifs>+\n"
    "Options:\n"
    "    -o <directory>   create the directory and write output files in it;\n"
    "                       not compatible with -oc\n"
    "    -oc <directory>  create the directory but if it already exists allow\n"
    "                       overwriting the contents; default: spamo_out\n"
    "    -minscore <valu> the minimum score (bits) to accept as a motif match;\n"
    "                       default: 7\n"
    "                       Note: if <valu> is in the range [-1, 0), then the\n"
    "                       minimum score is set to: -<valu>*pwm_maximum_score\n"
    "    -margin <size>   edge margin excluded for primary motif matches and\n"
    "                       the maximum distance from the primary motif to the\n"
    "                       secondary motif; default: 150\n"
    "    -bin <size>      size of bins used for output; default: 1\n"
    "    -range <size>    the range from the primary to include in significance\n"
    "                       tests; default: 150\n"
    "    -usebestsec      use only the best match of the secondary motif;\n"
    "                       default: count all secondary matches above the\n"
    "                       score match threshold in the margins around the\n"
    "                       primary motif match\n"
    "    -shared <fract>  fraction of shared trimmed sequence content that\n"
    "                       is required to exclude the sequence as redundant;\n"
    "                       set <fract> to 0 to skip sequence redundancy check;\n"
    "                       default: 0.5\n"
    "    -odds <odds>     odds ratio used for speedup of the redundant\n"
    "                       sequence check; low values may cause some \n"
    "                       redundant sequences to be missed; set <odds>\n"
    "                       to 0 to do full-length check;\n"
    "                       default: 20\n"
    "    -cutoff <pvalue> cutoff for spacings to be considered significant;\n"
    "                       default: 0.05\n"
    "    -evalue <evalue> minimum E-value for printing secondary motifs;\n"
    "                       default: 10\n"
    "    -overlap <size>  number of bases that the most significant spacing\n"
    "                       must overlap before further redundancy testing is\n"
    "                       done; default: 2\n"
    "    -joint <fract>   fraction of sites making up the most significant\n"
    "                       spacing that must be in both sets for the less\n"
    "                       significant motif to be considered redundant;\n"
    "                       default: 0.5\n"
    "    -pseudo <count>  pseudocount added to loaded motifs;\n"
    "                       default: 0.1\n"
    "    -bgfile <file>   file containing background frequency information\n"
    "                       used to apply pseudocounts to motifs and create the\n"
    "                       scoring matrices; default: calculate from sequences\n"
    "    -xalph           Convert the alphabet of the secondary motif databases\n"
    "                       to the alphabet of the primary motif\n"
    "                       assuming the core symbols of the secondary motif\n"
    "                       alphabet are a subset; default: reject differences\n"
    "    -trim <bits>     trim the edges of motifs based on the information\n"
    "                       content; positions on the edges with information\n"
    "                       content less than bits will not be used in\n"
    "                       scanning\n"
    "                       default: 0.25\n"
    "    -numgen <seed>   specify the random seed for initializing the pseudo-\n"
    "                       random number generator used in breaking ties;\n"
    "                       the seed is included in the output so experiments\n"
    "                       can be repeated; special value 'time' seeds to the\n"
    "                       system clock; default: 1\n"
    "    -primary <name>  name of motif to select as the primary motif;\n"
    "                       overrides -primaryi\n"
    "    -primaryi <num>  index of motif to select as the primary motif\n"
    "                       counting from 1; default: 1\n"
    "    -keepprimary     by default all occurrences of the primary other\n"
    "                       than the chosen one are erased to reduce the\n"
    "                       number of motifs whose apparent enrichment is due\n"
    "                       to multiple occurrences of the primary motif. If\n"
    "                       the same file is specified for the primary and\n"
    "                       secondary motifs, then this option will override\n"
    "                       erasing and primary-primary spacings will be\n"
    "                       analyzed.\n"
    "    -inc <pattern>   name pattern to select as secondary motif; may be\n"
    "                       repeated; default: all motifs are used\n"
    "    -exc <pattern>   name pattern to exclude as secondary motif; may be\n"
    "                       repeated; default: all motifs are used\n"
    "    -text            output text file only (no HTML or XML)\n"
    "    -eps             output histograms in eps format; usable with -png\n"
    "    -png             output histograms in png format; usable with -eps\n"
    "    -dumpseqs        dump matching trimmed sequences to output files;\n"
    "                       matches are initially in sequence name order;\n"
    "                       see documentation for column descriptions\n"
    "    -dumpsigs        same as above except only secondary matches in\n"
    "                       in significant bins are dumped\n"
    "                       matches are initially in sequence name order;\n"
    "                       see documentation for column descriptions\n"
    "    -help            print this usage message\n"
    "    -verbosity <v>   set the verbosity of the output: 1 (quiet) - 5 (dump);\n"
    "                       default: 2 (normal)\n"
    "    -version         print the version and exit\n"
    "Description:\n"
    "    SpaMo looks for significant spacings between a primary motif and\n"
    "    each motif in a library of secondary motifs.\n";

  if (format) {
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fputs(usage, stderr);
    fflush(stderr);
  } else {
    puts(usage);
  }
  cleanup_options(options);
  if (format) exit(EXIT_FAILURE);
  exit(EXIT_SUCCESS);
}

// define constants for the arguments that would clash with others if given
// a single character code
#define OUTPUT 500
#define CLOBBER 501
#define CREATE_EPS 502
#define CREATE_PNG 503
#define PRIMARY_INDEX 504
#define PSEUDOCOUNT 505
#define ARG_BACKGROUND 506
#define MINSCORE 508
#define SHOW_VERSION 510
#define ODDS_RATIO 511
#define TEXT_ONLY 512

/**************************************************************************
 * Checks the consistency of arguments and loads them into the options 
 * structure for easy access.  Returns the command line.
 **************************************************************************/
static char *process_arguments(int argc, char **argv, SPAMO_OPTIONS_T *options) {
  int option_index = 1;
  bool bad_argument = false;
  char *endptr;
  // Note options are specified as follows:
  // <name> <has argument> <(not used)> <int to return>
  struct option spamo_options[] = {
    {"o", required_argument, NULL, OUTPUT},
    {"oc", required_argument, NULL, CLOBBER},
    {"eps", no_argument, NULL, CREATE_EPS},
    {"png", no_argument, NULL, CREATE_PNG},
    {"dumpseqs", no_argument, NULL, 'd'},
    {"dumpsigs", no_argument, NULL, 'z'},
    {"numgen", required_argument, NULL, 'n'},
    {"margin", required_argument, NULL, 'm'},
    {"minscore", required_argument, NULL, MINSCORE},
    {"usebestsec", no_argument, NULL, 'u'},
    {"bin", required_argument, NULL, 'b'},
    {"range", required_argument, NULL, 'r'},
    {"shared", required_argument, NULL, 's'},     //sequence similarity threshold
    {"cutoff", required_argument, NULL, 'c'},     //pvalue cutoff
    {"evalue", required_argument, NULL, 'e'},     //Evalue cutoff
    {"overlap", required_argument, NULL, 'o'},    //redundant site overlap
    {"joint", required_argument, NULL, 'j'},      //redundant joint fraction
    {"pseudocount", required_argument, NULL, PSEUDOCOUNT},
    {"bgfile", required_argument, NULL, ARG_BACKGROUND},
    {"trimmotifs", required_argument, NULL, 't'},
    {"primary", required_argument, NULL, 'p'},
    {"primaryi", required_argument, NULL, PRIMARY_INDEX},
    {"keepprimary", no_argument, NULL, 'k'},
    {"xalph", no_argument, NULL, 'a'},
    {"include", required_argument, NULL, 'i'},
    {"exclude", required_argument, NULL, 'x'},
    {"help", no_argument, NULL, 'h'},
    {"verbosity", required_argument, NULL, 'v'},
    {"version", no_argument, NULL, SHOW_VERSION},
    {"odds", required_argument, NULL, ODDS_RATIO},
    {"text", no_argument, NULL, TEXT_ONLY},
    {NULL, 0, NULL, 0} //boundary indicator
  };

  // set option defaults
  options->start = time(NULL);
  options->prng_seed = 1;
  options->clobber = true;
  options->outdir = "spamo_out";
  options->use_best_secondary = false;
  options->create_eps = false;
  options->create_png = false;
  options->create_for_redundant = true;
  options->dump_seqs = false;
  options->margin = 150;
  options->bin = 1;
  options->test_range = 150;
  options->trim_bit_threshold = 0.25;
  options->sequence_similarity = 0.5;
  options->score_threshold = 7;
  options->pvalue_cutoff = 0.05;
  options->dump_sig_only = false;
  options->motif_evalue_cutoff = 10.0;
  options->overlap = 2;
  options->joint = 0.5;
  options->primary_name = rbtree_create(rbtree_strcmp, NULL, NULL, NULL, NULL);
  options->primary_index = rbtree_create(rbtree_intcmp, rbtree_intcpy, free, NULL, NULL);
  options->keep_primary = false;
  options->include_patterns = arraylst_create();
  options->exclude_patterns = arraylst_create();
  options->sequences_file = NULL;
  options->primary_file = NULL;
  options->secondary_files = arraylst_create();
  options->odds = 20;
  options->text_only = false;
  options->xalph = false;

  //meme file loading defaults (same as FIMO)
  options->pseudocount = 0.1;
  options->bg_filename = NULL;

  // parse optional arguments
  int index;
  while (1) {
    int opt = getopt_long_only(argc, argv, "", spamo_options, NULL);
    if (opt == -1) break;
    switch (opt) {
      case OUTPUT:        //-o <dir>
        options->clobber = false;
      case CLOBBER:       //-oc <dir>
        options->outdir = optarg;
        break;
      case CREATE_EPS:    //-eps
        options->create_eps = true;
        break;
      case CREATE_PNG:    //-png
        options->create_png = true;
        break;
      case MINSCORE:    //-minscore
        options->score_threshold = strtod(optarg, NULL);
        break;
      case 'u':        //-usebestsec
        options->use_best_secondary = true;
        break;
      case 'd':           //-dumpseqs
        options->dump_seqs = true;
        break;
      case 'z':           //-dumpsigs
        options->dump_seqs = true;
        options->dump_sig_only = true;
        break;
      case 'n':           //-numgen <seed>|time
        if (strcmp(optarg, "time") == 0) {
          options->prng_seed = options->start;
        } else {
          options->prng_seed = strtoul(optarg, NULL, 10);
        }
        break;
      case 'm':           //-margin <num>
        options->margin = strtol(optarg, NULL, 10);
        break;
      case 'b':           //-bin <num>
        options->bin = strtol(optarg, NULL, 10);
        break;
      case 'r':           //-range <num>: pvalue testing range
        options->test_range = strtol(optarg, NULL, 10);
        break;
      case 's':           //-shared <fraction>: shared sequence
        options->sequence_similarity = strtod(optarg, NULL);
        break;
      case 'c':           //-cutoff <pvalue>: pvalue significance cutoff
        options->pvalue_cutoff = strtod(optarg, NULL);
        break;
      case 'e':           //-evalue <Evalue>: motif Evalue significance cutoff
        options->motif_evalue_cutoff = strtod(optarg, NULL);
        break;
      case 'o':           //-overlap <fraction>: minimum peak overlap
        options->overlap = (int)strtol(optarg, NULL, 10);
        break;
      case 'j':           //-joint <fraction>: minimum shared sites
        options->joint = strtod(optarg, NULL);
        break;
      case PSEUDOCOUNT:   //-pseudocount <fraction>: motif pseudocount
        options->pseudocount = strtod(optarg, NULL);
        if (options->pseudocount < 0) {
          usage(options, "Pseudocount must be positive but got \"%s\".", optarg);
        }
        break;
      case ARG_BACKGROUND:    //-bgfile <file>: background to distribute pseudocount
        options->bg_filename = optarg;
        break;
      case 't':           //-trimmotifs <bits>: bit threshold for trimming
        options->trim_bit_threshold = strtod(optarg, NULL);
        break;
      case 'p':           //-primary <name>: name of primary motif
        rbtree_make(options->primary_name, optarg, optarg);
        break;
      case PRIMARY_INDEX: //-primaryi <index>: index of primary motif
        index = strtol(optarg, &endptr, 10);
        if ((errno != 0 && index == 0) || *endptr != '\0') {
          usage(options, "Option -primaryi '%s' was not in integer.", optarg);
        } else if (index <= 0) {
          usage(options, "Option -primaryi '%s' cannot be zero or negative.", optarg);
        }
        index--;
        rbtree_make(options->primary_index, &index, &index);
        break;
      case 'k':           //-keepprimary: don't exclude the primary motif from the matched
        options->keep_primary = true;
        break;
      case 'a':           //-xalph
        options->xalph = true;
        break;
      case 'i':           //-inc <pattern>
        arraylst_add(optarg, options->include_patterns);
        break;
      case 'x':           //-exc <pattern>
        arraylst_add(optarg, options->exclude_patterns);
        break;
      case TEXT_ONLY:           //-text
        options->text_only = true;
        break;
      case 'h':           //-help
        usage(options, NULL);
        break;
      case 'v':           //-verbosity <1-5>
        verbosity = (int)strtol(optarg, NULL, 10);
        if (verbosity < 1 || verbosity > 5) {
          usage(options, "Verbosity must be between 1 and 5 inclusive. Got \"%s\".", optarg);
        }
        break;
      case SHOW_VERSION:
        fprintf(stdout, VERSION "\n");
        exit(EXIT_SUCCESS);
        break;
      case ODDS_RATIO:
        options->odds = strtof(optarg, NULL);
        if (options->pseudocount < 0) {
          usage(options, "Odds ratio must be positive but got \"%s\".", optarg);
        }
        break;
      case '?':           //unrecognised or ambiguous argument
        bad_argument = true;
    }
  }
  if (bad_argument) usage(options, "One or more unknown or ambiguous options were supplied.");
  option_index = optind;
  // parse required arguments
  // get sequences file
  if (option_index >= argc) usage(options, "No sequences file!");
  options->sequences_file = argv[option_index++];
  if (!file_exists(options->sequences_file))
    usage(options, "Sequences file \"%s\" does not exist!", options->sequences_file);

  if (rbtree_size(options->primary_name) + rbtree_size(options->primary_index) > 1)
    usage(options, "You may only specify up to one primary motif with -primary or -primaryi!\n");

  // free the index if a primary name was given
  if (rbtree_size(options->primary_name) > 0) {
    rbtree_destroy(options->primary_index);
    options->primary_index = NULL;
  }

  // get primary motif file
  if (option_index >= argc) usage(options, "No primary motif file!");
  options->primary_file = argv[option_index++];
  if (!file_exists(options->primary_file))
    usage(options, "Primary motif file \"%s\" does not exist!", options->primary_file);

  // get secondary motifs
  if (option_index >= argc) usage(options, "No secondary motif files!");
  while (option_index < argc) {
    // get secondary motif file
    if (!file_exists(argv[option_index]))
      usage(options, "Secondary motif file \"%s\" does not exist!", argv[option_index]);
    arraylst_add(argv[option_index++], options->secondary_files);
  }

  // make enough space for all the command line options, with one space between each
  int line_length = 0;
  int i;
  for (i = 0; i < argc; i++) line_length += strlen(i == 0 ? basename(argv[0]) : argv[i]);
  // add on argc to allow one char per word for separating space + terminal '\0'
  char *commandline = (char*) mm_malloc(sizeof(char)*(line_length+argc));
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
}

/**************************************************************************
 * Create the output directory.
 **************************************************************************/
static void create_spamo_output_directory(SPAMO_OPTIONS_T *options) {
  //make output directory
  if (create_output_directory(options->outdir, options->clobber, verbosity >= NORMAL_VERBOSE)) {
    die("Can not continue without an output directory.\n");
  }
}

/**************************************************************************
 * Loads the primary motif using the information specified in the options
 * and the pre-calculated background or from the options->bg_filename file.
 **************************************************************************/
static MOTIF_DB_T *load_primary_motif(
  SPAMO_OPTIONS_T *options, 
  ALPH_T *alph,
  ARRAY_T **background, 
  MOTIF_T **primary_motif, 
  bool *stdin_used
) {

  // Load the primary motif and the background.
  if (options->primary_name && rbtree_size(options->primary_name) > 0) {
    RBNODE_T *node = rbtree_first(options->primary_name);
    DEBUG_FMT(NORMAL_VERBOSE, "Loading Primary Motif (name=%s)\n", (char *)rbtree_value(node));
  } else if (options->primary_index && rbtree_size(options->primary_index) > 0) {
    RBNODE_T *node = rbtree_first(options->primary_index);
    DEBUG_FMT(NORMAL_VERBOSE, "Loading Primary Motif (index=%d)\n", *((int *)rbtree_value(node))+1);
  } else {
    DEBUG_MSG(NORMAL_VERBOSE, "Loading Primary Motif\n");
  }
  MOTIF_DB_T *primary_db = read_motifs_and_background(
    0,				// index to give DB
    options->primary_file,	// motif file name (or special word)
    "Primary motif",		// type of database for error messages
    options->primary_name,	// get one motif by name
    options->primary_index,	// get one motif by index
    NULL,			// get set of motifs by name
    NULL,			// exclude set of motifs by name
    true,			// allow motifs with zero probability entries
    false,			// create RC copies, appended
    options->pseudocount,	// multiply times background model
    true,			// set_trim
    options->trim_bit_threshold,// trim_bit_threshold
    &(options->bg_filename),	// background source; may be changed
    true, 			// make bg symmetrical if alph complementable
    background,			// will be set if id==0
    options->sequences_file,	// sequence file name
    alph,			// primary motif alphabet
    false,			// set motif conversion alphabet
    false,			// don't remove extension from name
    true,			// remove ".meme" extension from name
    false,			// don't replace underscores in name
    stdin_used			// IN/OUT check or set if path is "-"
  );

  // Get the primary motif from the database.
  *primary_motif = (MOTIF_T *)arraylst_get(0, primary_db->motifs);
  assert(*primary_motif);

  return(primary_db);
} // load_primary_motif

/**************************************************************************
 * Loads the sequences, finds the primary motif matches and trims
 * the sequences to only the portion surrounding the best primary match.
 * Also filters unmatched and redundant sequences.
 **************************************************************************/
static void calculate_trimmed_filtered_sequences(
  SPAMO_OPTIONS_T *options, 
  ALPH_T *alph,
  MOTIF_T *primary_motif, 
  ARRAY_T *background, 
  char **trimmed_sequence_data, 
  SEQUENCE_DB_T **sequence_db,
  RBTREE_T **sequences
) {
  RBTREE_T *seqs;
  RBNODE_T *node, *node2, *next;
  SEQUENCE_T *sequence, *sequence2;
  SEQ_T **all_sequences, *current;
  int count, i, j, trim_len, comp_len, match, match2, similarity, 
   pos_similarity, neg_similarity, name_len, duplicate_id_count,
   tooshort_count, nomatch_count, duplicate_count, width, margin;
  char *seq_dest, *p1, *p2, *n1, *name, *revcomp1, *revcomp2;
  FILE *fasta_fp;
  bool created, revcomp;
  
  revcomp = alph_has_complement(alph);

  *sequence_db = create_sequence_db(options->sequences_file);

  // calculate the trimmed length
  width = get_motif_trimmed_length(primary_motif);
  margin = options->margin;
  comp_len = 2 * margin;
  trim_len = 2 * margin + width;

  // Create reverse complement of sequence for similarity check.
  if (revcomp) {
    revcomp1 = (char*)mm_malloc(sizeof(char) * (trim_len+1));
    revcomp2 = (char*)mm_malloc(sizeof(char) * (trim_len+1));
  } else {
    revcomp1 = NULL;
    revcomp2 = NULL;
  }

  // Create a permutation of the integers [0..trim_len-1] for similarity check.
  // Note: the primary motif positions are not included for speed.
  int *perm = (int *)mm_malloc(sizeof(int) * comp_len); 
  for (i=0; i<comp_len; i++) perm[i] = (i < margin) ? i : i+width;
  SHUFFLE(perm, comp_len);

  // Calculate the maximum number of matches to achieve a given
  // odds-ratio for each number of comparisons given:
  //  1) the maximum desired similarity
  //  2) the background frequencies of the letters
  int *min_similarity = (int*) mm_malloc((comp_len+1) * sizeof(int));
  {
    // TODO check that this makes sense
    double bg_match_prob = 0;
    for (i = 0; i < alph_size_core(alph); i++) {
      double bgfreq = get_array_item(i, background);
      bg_match_prob += (bgfreq * bgfreq);
    }

    if (bg_match_prob >= options->sequence_similarity) 
      die("Value of -shared (%.2f) must be greater than expected similarity (%.2f).\n", options->sequence_similarity, bg_match_prob);

    // Lookup-table for the mininum number of matches required at each length to be too similar.
    // Value is the number of minimum number matching positions for the desired odds ratio.
    double p1 = bg_match_prob;
    double q1 = 1 - p1;
    double p2 = options->sequence_similarity;
    double q2 = 1 - p2;
    int r;
    double odds = options->odds; // Desired odds ratio.
    for (r=1; r<=comp_len; r++) {
      int t = (log(odds) - (r * log(q1/q2))) / log( (p1*q2)/(p2*q1) );
      double odds = pow(p1/p2, t) * pow(q1/q2, r-t);
      min_similarity[r] = (odds == 0 || (t == 0 && odds < odds)) ? -1 : t;
      //printf("r %d t %d %d odds_ratio %.3g\n", r, t, min_similarity[r], odds);
    }
  }

  // read the entire fasta file into memory...
  DEBUG_MSG(NORMAL_VERBOSE, "Loading All Sequences\n");
  if (!open_file(options->sequences_file, "r", false, "FASTA sequences", 
        "sequences", &fasta_fp)) exit(EXIT_FAILURE);
  count = 0;
  read_many_fastas(alph, fasta_fp, MAX_SEQ, &count, &all_sequences);
  fclose(fasta_fp);
  DEBUG_FMT(NORMAL_VERBOSE, "Loaded %d Sequences\n", count);
  (*sequence_db)->loaded = count;

  duplicate_id_count = 0;
  tooshort_count = 0;

  // create our set, ensuring no duplicate ids
  //seqs = rbtree_create(rbtree_strcmp, rbtree_strcpy, free, NULL, free);
  seqs = rbtree_create(rbtree_strcmp, rbtree_strcpy, free, NULL, free);
  for (i = 0; i < count; ++i) {
    current = all_sequences[i];

    // check that the sequence is long enough to have at least a single site for the primary
    if (get_seq_length(current) < trim_len) {
      tooshort_count++;
      continue;
    }

    node = rbtree_lookup(seqs, get_seq_name(current), true, &created);
    if (!created) {
      duplicate_id_count++;
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(stderr, "Skipping duplicate sequence identifier \"%s\" in FASTA file\n", 
        get_seq_name(current));
      }
      continue;
    } 
    sequence = (SEQUENCE_T*)mm_malloc(sizeof(SEQUENCE_T));
    sequence->idx = -1;
    sequence->contributes = false;
    sequence->index = -1;
    sequence->length = get_seq_length(current);
    sequence->data = get_raw_sequence(current);
    sequence->name = get_seq_name(current);
    sequence->primary_matches = NULL;
    sequence->trimmed_primary_matches = NULL;
    rbtree_set(seqs, node, sequence);
  }

  if (duplicate_id_count > 0) {
    DEBUG_FMT(QUIET_VERBOSE, "Warning: Eliminated %d sequences that had the "
      "same ID as another sequence.\nSequences must have unique IDs.\n", 
      duplicate_id_count);
  }
  (*sequence_db)->excluded_duplicate_id = duplicate_id_count;

  if (tooshort_count > 0) {
    DEBUG_FMT(QUIET_VERBOSE, "Warning: Eliminated %d sequences that were too "
      "short to be scanned with the current margin (%d) and primary motif "
      "length (%d). Sequences must be at least %d long.\n", tooshort_count, 
      margin, get_motif_trimmed_length(primary_motif), trim_len);
  }
  (*sequence_db)->excluded_tooshort = tooshort_count;

  // scan the primary motif
  DEBUG_MSG(NORMAL_VERBOSE, "Determining Best Primary Matches\n");
  scan_spamo_primary(margin, options->score_threshold, background, primary_motif, seqs, false);

  // remove sequences that don't have a match
  DEBUG_MSG(NORMAL_VERBOSE, "Eliminating Sequences Without Primary Matches\n");
  nomatch_count = 0;
  node = rbtree_first(seqs);
  while (node != NULL) {
    next = rbtree_next(node);
    sequence = (SEQUENCE_T*)rbtree_value(node);
    if (!(sequence->primary_matches)) {
      DEBUG_FMT(DUMP_VERBOSE, "Eliminating \"%s\": no primary match.\n", sequence->name);
      //rbtree_delete(seqs, node, NULL, NULL);
      SEQUENCE_T *del_seq;
      rbtree_delete(seqs, node, NULL, (void *) &del_seq);
      if (del_seq->primary_matches) free_array(del_seq->primary_matches);
      free(del_seq);
      nomatch_count++;
    }
    node = next;
  }
  DEBUG_FMT(NORMAL_VERBOSE, "Eliminated %d Sequences Without Primary Match\n", nomatch_count);
  (*sequence_db)->excluded_nomatch = nomatch_count;

  // Remove sequences that have ambigs in margins if we are removing
  // overly-similar sequences.  This is necessary since sequences
  // for the ambiguity heuristic to work correctly.
  (*sequence_db)->excluded_ambigs = 0;
  if (options->sequence_similarity > 0 && options->sequence_similarity < 1 ) {
    DEBUG_MSG(NORMAL_VERBOSE, "Eliminating Sequences With Ambiguous Characters in Margins\n");
    nomatch_count = 0;
    node = rbtree_first(seqs);
    while (node != NULL) {
      next = rbtree_next(node);
      sequence = (SEQUENCE_T*)rbtree_value(node);
      int match = get_array_item(0, sequence->primary_matches);
      if (match < 0) match = -match;
      p1 = sequence->data+(match - margin - 1);
      for (i=0; i<trim_len; i++) {
        if (alph_is_ambiguous(alph, p1[i])) {
          DEBUG_FMT(DUMP_VERBOSE, "Eliminating \"%s\": ambiguous characters in margins.\n", sequence->name);
          //rbtree_delete(seqs, node, NULL, NULL);
	  SEQUENCE_T *del_seq;
	  rbtree_delete(seqs, node, NULL, (void *) &del_seq);
	  if (del_seq->primary_matches) free_array(del_seq->primary_matches);
	  free(del_seq);
          nomatch_count++;
          break;
        }
      }
      node = next;
    }
    DEBUG_FMT(NORMAL_VERBOSE, "Eliminated %d Sequences With Ambiguous Characters in Margins\n", nomatch_count);
    (*sequence_db)->excluded_ambigs = nomatch_count;
  }

  // Remove sequences that are too similar to each other in the margin around the match
  if (options->sequence_similarity > 0 && options->sequence_similarity < 1 ) {
    DEBUG_MSG(NORMAL_VERBOSE, "Eliminating Similar Sequences\n");

    double seq_similarity = options->sequence_similarity;
    duplicate_count = 0;

    double j_sum = 0;
    double n_comp = 0;
    for (node = rbtree_first(seqs); node != NULL; node = rbtree_next(node)) {
      sequence = (SEQUENCE_T*)rbtree_value(node);

      // Get position of best match
      match = get_array_item(0, sequence->primary_matches);
      if (match < 0) match = -match;

      // Create reverse complement of sequence for similarity check.
      p1 = sequence->data+(match - margin - 1);
      if (revcomp) copy_string_with_rc(alph, p1, revcomp1, trim_len, true);

      // Compare with second sequence on both strands.
      for (node2=rbtree_next(node); node2 != NULL; node2=next) {
        next = rbtree_next(node2);
        sequence2 = (SEQUENCE_T*)rbtree_value(node2);
        match2 = get_array_item(0, sequence2->primary_matches);
        if (match2 < 0) match2 = -match2;

        // TLB: Check similarity with reverse complement of second sequence, too.
        p1 = sequence->data+(match - margin - 1);
        p2 = sequence2->data+(match2 - margin - 1);
        n1 = revcomp1;
        pos_similarity = neg_similarity = 0;
        int j;
        for (j = 0; j < comp_len; j++) {
          i = perm[j];                  // Compare positions in random order.
          if (p1[i] == p2[i]) {
            ++pos_similarity;
          } 
          if (revcomp && n1[i] == p2[i]) {
            ++neg_similarity;
          } 

          if (pos_similarity <= min_similarity[j+1] && neg_similarity < min_similarity[j+1]) break;
        } // compare sequences
        j_sum += j;
        n_comp++;
        similarity = pos_similarity > neg_similarity ? pos_similarity : neg_similarity;

        // Remove second sequence if too similar to either strand of first sequence.
        if (j == comp_len && (double)similarity/comp_len > seq_similarity) {
          if (verbosity == DUMP_VERBOSE && neg_similarity > pos_similarity) {   // create reverse complement string
            assert(revcomp);
            char *seq2 = sequence2->data+(match2 - margin - 1);
            copy_string_with_rc(alph, seq2, revcomp2, trim_len, true);
          }
          DEBUG_FMT(DUMP_VERBOSE, "Eliminating \"%s\": %d%% similarity to \"%s %s\".\n%*.*s\n%*.*s\n", 
              sequence2->name, 
              (int)(((double)similarity / comp_len) * 100 + 0.5), 
              neg_similarity > pos_similarity ? "(rc)" : "",
              sequence->name, 
              trim_len, trim_len, sequence->data+(match - margin - 1),
              trim_len, trim_len, neg_similarity > pos_similarity ? revcomp2 : sequence2->data+(match2 - margin - 1)
          );
          //rbtree_delete(seqs, node2, NULL, NULL);
	  SEQUENCE_T *del_seq;
	  rbtree_delete(seqs, node2, NULL, (void *) &del_seq);
	  if (del_seq->primary_matches) free_array(del_seq->primary_matches);
	  free(del_seq);
          duplicate_count++;
        }
      } // node2
    } // node1
    DEBUG_FMT(NORMAL_VERBOSE, "Eliminated %d Similar Sequences\n", duplicate_count);
    (*sequence_db)->excluded_similar = duplicate_count;
  } // eliminate similar sequences

  DEBUG_MSG(NORMAL_VERBOSE, "Trimming Sequences\n");
  //allocate space for the trimmed sequences (in sequential memory)
  *trimmed_sequence_data = mm_malloc(sizeof(char) * rbtree_size(seqs) * (trim_len + 1));
  //trim the sequences down to the margin around the match and allocate an index to each node in the tree
  seq_dest = *trimmed_sequence_data;
  for (node = rbtree_first(seqs), i = 0; node != NULL; node = rbtree_next(node), ++i) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    if (verbosity == DUMP_VERBOSE) fprintf(stdout, ">%s\n%s\n", sequence->name, sequence->data);
    //allocate an index
    sequence->index = i;
    //copy the sequence to trim
    match = get_array_item(0, sequence->primary_matches);
    if (match < 0) match = -match;
    memcpy(seq_dest, sequence->data+(match - margin - 1), sizeof(char) * trim_len);
    sequence->data = seq_dest;
    sequence->length = trim_len;
    seq_dest+=trim_len;
    seq_dest[0] = '\0';
    seq_dest++;
    //copy the sequence name
    name = sequence->name;
    name_len = strlen(name);
    sequence->name = mm_malloc(sizeof(char) * (name_len + 1));
    memcpy(sequence->name, name, sizeof(char) * name_len);
    sequence->name[name_len] = '\0';
  }

  //now we have the copied the name and data so we must change the free function
  rbtree_alter_value_free(seqs, destroy_sequence);

  // Erase any extra primary matches.
  int erased_count = 0;
  if (! options->keep_primary) {
    DEBUG_MSG(NORMAL_VERBOSE, "Erasing Extra Primary Matches\n");
    // scan the trimmed sequences with the primary motif 
    scan_spamo_primary(0, options->score_threshold, background, primary_motif, seqs, true);
    node = rbtree_first(seqs);
    while (node != NULL) {
      next = rbtree_next(node);
      sequence = (SEQUENCE_T*)rbtree_value(node);
      int n_primary_matches = get_array_length(sequence->trimmed_primary_matches);
      if (n_primary_matches > 1) {
        DEBUG_FMT(DUMP_VERBOSE, "Erasing %d extra matches in sequence \"%s\".\n", n_primary_matches, sequence->name);
        for (i = 0; i < n_primary_matches; i++) {
          int position = get_array_item(i, sequence->trimmed_primary_matches);
          position = (position < 0) ? -(position+1) : position-1; // positions are 1-relative 
          if (position == margin) continue; // Don't erase main primary match.
          // Erase extra primary match
          for (j = 0; j < width; j++, position++) sequence->data[position] = alph_wildcard(alph);
        }
        erased_count += n_primary_matches - 1;
      }
      node = next;
    }
    DEBUG_FMT(NORMAL_VERBOSE, "Erased %d Extra Primary Matches\n", erased_count);
  }
  (*sequence_db)->erased_primary_matches = erased_count;

  //deallocate the SEQ_T sequences
  DEBUG_MSG(NORMAL_VERBOSE, "Freeing Unneeded Sequences\n");
  for (i = 0; i < count; ++i) {
    current = all_sequences[i];
    free_seq(current); 
  }
  free(all_sequences);
  *sequences = seqs;
  free(perm);
  free(min_similarity);
  if (revcomp1)free(revcomp1);
  if (revcomp2)free(revcomp2);
}

/**************************************************************************
 * Loads the secondary motifs filtering them with the information in the
 * options structure.  Returns a list of MOTIF_DB_T.
 **************************************************************************/
static ARRAYLST_T *load_secondary_motifs(
  SPAMO_OPTIONS_T *options, 
  ALPH_T *alph,
  ARRAY_T *background, 
  RBTREE_T **secondary_motifs,	// IN/OUT
  bool *stdin_used		// IN/OUT check or set if path is "-"
) {
  int i;

  ARRAYLST_T *secondary_dbs = arraylst_create_sized(arraylst_size(options->secondary_files));
  *secondary_motifs = rbtree_create(secondary_key_compare, secondary_key_copy, free, NULL, NULL);
  RBTREE_T *secondary_names = rbtree_create(rbtree_strcmp, NULL, NULL, NULL, NULL);

  for (i = 0; i < arraylst_size(options->secondary_files); ++i) { //loop over the motif db files

    // Get the file name of the next motif db.
    char *secondary_file = (char*)arraylst_get(i, options->secondary_files);
    DEBUG_FMT(NORMAL_VERBOSE, "Loading Secondary Motifs from file '%s'\n", secondary_file);

    // Load the motifs from this file.
    MOTIF_DB_T *secondary_db = read_motifs_and_background(
      i+1, 			// start the id at 1 because we use 0 for the primary
      secondary_file,		// motif file name (or special word)
      "Secondary motif",	// type of database for error messages
      NULL,			// get one motif by name
      NULL,			// get one motif by index
      options->include_patterns,// get set of motifs by name (or NULL)
      options->exclude_patterns,// exclude this set of motifs by name (or NULL)
      true,			// allow motifs with zero probability entries
      false,			// create RC copies, appended
      options->pseudocount,	// multiply times background model
      true,			// set_trim
      options->trim_bit_threshold,// trim_bit_threshold
      NULL,			// background source
      true, 			// make bg symmetrical if alph complementable
      &background,		// will be set if id==0
      NULL,			// sequence file name
      alph,			// primary motif alphabet
      options->xalph,		// set motif conversion alphabet
      false,			// don't remove extension from name
      true,			// remove ".meme" extension from name
      false,			// don't replace underscores in name
      stdin_used		// IN/OUT check and set if path is "-"
    );

    // Post-process the motifs to make them the secondary motif type.
    // Reset the warning type for the warning clumping.
    int i;
    int warn_type = NO_WARNING;
    for (i=0; i<arraylst_size(secondary_db->motifs); i++) {
      MOTIF_T *motif = arraylst_get(i, secondary_db->motifs);
      // Check that the motif fits in the sequence margin.
      if ((options->margin - get_motif_trimmed_length(motif) + 1) < 1) {
        clump_motif_db_warning(&warn_type, MARGIN_WARNING, "Warning: The following "
	  "secondary motifs in '%s' were excluded as they didn't fit in "
	  "the margin:\n  ", secondary_file, get_motif_id(motif));
        //destroy_motif(motif);
        secondary_db->loaded--;
        secondary_db->excluded++;
        continue;
      }

      // Check for motif name uniqueness across all DBs.
      SECONDARY_KEY_T key;
      key.db_id = 0;
      key.motif_id = get_motif_id(motif);
      bool created;
      RBNODE_T *node = rbtree_lookup(secondary_names, get_motif_id(motif), true, &created);
      if (!created) {
        clump_motif_db_warning(&warn_type, DUPLICATE_WARNING, "Warning: The following "
          "duplicate secondary motifs in '%s' were excluded:\n  ", secondary_file, get_motif_id(motif));
        //destroy_motif(motif);
        secondary_db->loaded--;
        secondary_db->excluded++;
        continue;
      }
      key.db_id = secondary_db->id;
      key.motif_id = get_motif_id(motif);
      node = rbtree_lookup(*secondary_motifs, &key, true, &created);
      rbtree_set(*secondary_motifs, node, create_secondary_motif(options->margin, options->bin, secondary_db, motif));
    }
    arraylst_add(secondary_db, secondary_dbs);
    if (warn_type && verbosity >= NORMAL_VERBOSE) fprintf(stderr, "\n");
  } // secondary_files

  // Clean up.
  rbtree_destroy(secondary_names);

  // Check that we found suitable motifs.
  if (rbtree_size(*secondary_motifs) == 0) die("No acceptable secondary motifs found.");

  return(secondary_dbs);
}

/**************************************************************************
 * Loads the matches to the secondary motifs by
 * directly interfacing with FIMO to scan the sequences.
 **************************************************************************/
static void load_secondary_matches(
  SPAMO_OPTIONS_T *options, 
  RBTREE_T *sequences, 
  MOTIF_T *primary_motif, 
  RBTREE_T *secondary_motifs, 
  ARRAY_T *background
) {
  RBNODE_T *node;
  MOTIF_DB_T *db;
  SECONDARY_MOTIF_T *smotif;
  int i, count, total, report_step, *hits, hits_size, test_max;
  ARRAY_T **matches = NULL;

  test_max = (int)(options->test_range / options->bin) + (options->test_range % options->bin ? 1 : 0);
  int n_secondary_motifs = rbtree_size(secondary_motifs);

  DEBUG_MSG(NORMAL_VERBOSE, "Scanning For Secondary Matches\n");
  matches = (ARRAY_T**)mm_calloc(sizeof(ARRAY_T*), rbtree_size(sequences));
  hits_size = 4 * options->margin;
  hits = (int*)mm_malloc(sizeof(int) * hits_size);
  total = rbtree_size(secondary_motifs);
  report_step = (int)(((double)total / 100) + 0.5);
  if (report_step == 0) report_step = 1; // mod 0 causes Arithmetic exception

  for (node = rbtree_first(secondary_motifs), count = 0; node != NULL; node = rbtree_next(node), count++) {
    smotif = (SECONDARY_MOTIF_T*)rbtree_value(node);
    //if (verbosity >= NORMAL_VERBOSE) fprintf(stderr, "40.40%s\r", get_motif_id(smotif->motif));
    //scan
    scan_spamo_secondary(options->margin, options->score_threshold, options->use_best_secondary, background,
        smotif->motif, sequences, matches, hits, hits_size);
    //process
    process_matches(options->margin, options->bin, options->pvalue_cutoff, options->motif_evalue_cutoff,
        test_max, primary_motif, sequences, smotif, n_secondary_motifs, matches);
    //dump sequences (if requested)
    if (options->dump_seqs && smotif->sig_count > 0 && smotif->passes_evalue_cutoff) {
      output_sequence_matches(options->outdir, options->margin, 
        options->bin, options->pvalue_cutoff, options->dump_sig_only,
        sequences, primary_motif, smotif, matches);
    }
    if (count % report_step == 0) {
      DEBUG_FMT(NORMAL_VERBOSE, "\rScanned %d%% ", (int)(((double)(count + 1) / total) * 100 + 0.5));
    } else {
      //DEBUG_MSG(NORMAL_VERBOSE, ".");
    }
  }
  DEBUG_MSG(NORMAL_VERBOSE, "\n");
  for (i = 0; matches && i < rbtree_size(sequences); i++) free_array(matches[i]);
  free(matches);
  free(hits);
}

/**************************************************************************
 * Compares two secondary motifs by best pvalue
 **************************************************************************/
int secondary_motif_pvalue_comparator(void *v1, void *v2) {
  SECONDARY_MOTIF_T* smotif1 = (SECONDARY_MOTIF_T*)v1;
  SECONDARY_MOTIF_T* smotif2 = (SECONDARY_MOTIF_T*)v2;
  //if (smotif1->sigs->pvalue < smotif2->sigs->pvalue) {
  if (smotif1->min_pvalue < smotif2->min_pvalue) {
    return -1;
  //} else if (smotif1->sigs->pvalue == smotif2->sigs->pvalue) {
  } else if (smotif1->min_pvalue == smotif2->min_pvalue) {
    return 0;
  } else {
    return 1;
  }
}

/**************************************************************************
 * Returns a sorted array of the secondary motifs with a significant spacing, 
 * sorted by best spacing pvalue
 **************************************************************************/
LINKLST_T* sort_secondary_motifs(RBTREE_T *secondary_motifs) {
  RBNODE_T *node;
  LINKLST_T *list;
  SECONDARY_MOTIF_T *smotif;
  DEBUG_MSG(NORMAL_VERBOSE, "Sorting Significant Secondary Motifs\n");
  //copy entries with significant match into a list
  list = linklst_create();
  for (node = rbtree_first(secondary_motifs); node != NULL; node = rbtree_next(node)) {
    smotif = (SECONDARY_MOTIF_T*)rbtree_value(node);
    if (smotif->passes_evalue_cutoff) linklst_add(smotif, list);
  }
  //sort by pvalue
  linklst_sort(secondary_motif_pvalue_comparator, list);
  return list;
}

/**************************************************************************
 * Calculates the maximum value in a bin to be used in output histograms
 **************************************************************************/
int calculate_bin_max(LINKLST_T* sorted_secondary_motifs) {
  LL_LINK_T *node;
  SECONDARY_MOTIF_T *smotif;
  int binmax; 
  //calculate the maximum histogram count
  binmax = 0;
  for (node = linklst_first(sorted_secondary_motifs); node != NULL; node = linklst_next(node)) {
    smotif = (SECONDARY_MOTIF_T*)linklst_get(node);
    if (smotif->max_in_one_bin > binmax) {
      binmax = smotif->max_in_one_bin;
    }
  }
  return binmax;
}

/**************************************************************************
 * Calculate the maximum overlap of the two best peaks
 **************************************************************************/
static inline int peak_overlap(int bin_size, SECONDARY_MOTIF_T *smot1, SECONDARY_MOTIF_T *smot2) {
  int mot1_closest, mot1_furthest, mot2_closest, mot2_furthest, overlap;

  // Get distances (closest/furthest) of each of the two secondary motifs from primary motif.
  mot1_closest = smot1->sigs->bin * bin_size + 1;
  mot1_furthest = (smot1->sigs->bin + 1) * bin_size + get_motif_trimmed_length(smot1->motif) - 1;
  mot2_closest = smot2->sigs->bin * bin_size + 1;
  mot2_furthest = (smot2->sigs->bin + 1) * bin_size + get_motif_trimmed_length(smot2->motif) - 1;

  //check for overlap
  if (mot1_furthest < mot2_closest) return false; //no overlap
  if (mot2_furthest < mot1_closest) return false; //no overlap

  if (mot1_closest < mot2_closest) {
    overlap = mot1_furthest - mot2_closest + 1;
  } else {
    overlap = mot2_furthest - mot1_closest + 1;
  }
  return overlap;
}

/**************************************************************************
 * Count the intersection of two sorted sets of numbers.
 **************************************************************************/
static inline int set_intersect_count(int *list1, int count1, int *list2, int count2) {
  int same = 0;
  while (count1 > 0 && count2 > 0) {
    if (*list1 == *list2) {
      ++same;
      ++list1;
      --count1;
      ++list2;
      --count2;
    } else if (*list1 < *list2) {
      ++list1;
      --count1;
    } else {
      ++list2;
      --count2;
    }
  }
  return same;
}

/**************************************************************************
 * Determines if a motif has a result set so similar that it should be
 * declared redundant.
 **************************************************************************/
int is_redundant(SPAMO_OPTIONS_T *options, SECONDARY_MOTIF_T *best, SECONDARY_MOTIF_T *other) {
  int same, overlap, i;
  double joint;

  // Return if either motif has no significant spacing
  if (best->sigs == NULL || other->sigs == NULL) return false;

  // Return if either has seq_count set to zero.
  if (best->seq_count == 0 || other->seq_count == 0) return false;

  // calculate the overlap
  // check that they're on the same side
  int orient_best = best->sigs->orient;
  int orient_other = other->sigs->orient;
  if (!(
        (LEFT_SIDE(orient_best) && LEFT_SIDE(orient_other)) || 
        (RIGHT_SIDE(orient_best) && RIGHT_SIDE(orient_other))
      )) return false; //different side

  //check that the overlap is large enough
  overlap = peak_overlap(options->bin, best, other);
  if (overlap < options->overlap) return false;

  //check that the sequence set is similar enough
  same = set_intersect_count(best->seqs, best->seq_count, other->seqs, other->seq_count);
  if (best->seq_count < other->seq_count) {
    joint = (double)same / best->seq_count;
  } else {
    joint = (double)same / other->seq_count;
  }
  return joint > options->joint;
}

/**************************************************************************
 * Groups secondary motifs
 **************************************************************************/
void group_secondary_motifs(SPAMO_OPTIONS_T *options, LINKLST_T* sorted_secondary_motifs) {
  LL_LINK_T *best_node, *other_node, *next_node;
  SECONDARY_MOTIF_T *other;
  GROUPED_MOTIF_T *group;
  DEBUG_MSG(NORMAL_VERBOSE, "Grouping Redundant Secondary Motifs\n");
  for (best_node = linklst_first(sorted_secondary_motifs); 
      best_node != NULL; best_node = linklst_next(best_node)) {
    group = create_grouped_motif((SECONDARY_MOTIF_T*)linklst_get(best_node));
    other_node = linklst_next(best_node);
    while (other_node != NULL) {
      next_node = linklst_next(other_node);
      other = (SECONDARY_MOTIF_T*)linklst_get(other_node);
      if (is_redundant(options, group->best, other)) {
        linklst_remove(other_node, sorted_secondary_motifs);
        linklst_add(other, group->others);
      }
      other_node = next_node;
    }
    linklst_set(group, best_node);
  }
}

/**************************************************************************
 * Outputs TSV results.
 **************************************************************************/
void output_tsv(
  int argc, 
  char **argv, 
  char *commandline,
  SPAMO_OPTIONS_T *options, 
  int bin_max, 
  SEQUENCE_DB_T *sequence_db,
  MOTIF_DB_T *primary_db,
  MOTIF_T *primary_motif, 
  ARRAYLST_T *secondary_dbs,
  LINKLST_T *secondary_motifs, 
  int n_secondary_motifs
) {
  LL_LINK_T *node;
  GROUPED_MOTIF_T *gmotif;

  DEBUG_MSG(NORMAL_VERBOSE, "Outputting TSV\n");

  // open centrimo TSV file
  char *file_path = make_path_to_file(options->outdir, TEXT_FILENAME);
  FILE *tsv_file = fopen(file_path, "w");
  free(file_path);
  fputs(
    "prim_db\tprim_id\tprim_alt\tprim_cons\t"
    "sec_db\tsec_id\tsec_alt\tsec_cons\t"
    "trim_left\ttrim_right\t"
    "red_db\tred_id\tred_alt\t"
    "E-value\t"
    "gap\t"
    "orient\t"
    "count\t"
    "total\t"
    "adj_p-value\t"
    "p-value\n",
    tsv_file
  );

  for (node = linklst_first(secondary_motifs); node != NULL; node = linklst_next(node)) {
    gmotif = (GROUPED_MOTIF_T*)linklst_get(node);
    output_secondary_motif_tsv(tsv_file, primary_db, primary_motif, NULL, gmotif->best, n_secondary_motifs, gmotif->others);
  }

  // Finish the TSV output
  char *version_message = "# SpaMo (Spaced Motif Analysis Tool): Version " VERSION " compiled on " __DATE__ " at " __TIME__ "\n";
  fprintf(tsv_file, "\n%s", version_message);
  fprintf(tsv_file, "# The format of this file is described at %s/%s.\n", SITE_URL, "doc/spamo-output-format.html");
  fprintf(tsv_file, "# %s\n", commandline);
  fclose(tsv_file);
} // output_tsv

/**************************************************************************
 * Outputs JSON for the motif database
 **************************************************************************/
void output_motif_db_json(JSONWR_T *json, MOTIF_DB_T *db) {
  jsonwr_start_object_value(json);
  jsonwr_str_prop(json, "source", db->source);
  jsonwr_str_prop(json, "name", db->name);
  jsonwr_time_prop(json, "last_modified", &(db->last_mod));
  jsonwr_lng_prop(json, "loaded", db->loaded);
  jsonwr_lng_prop(json, "excluded", db->excluded);
  jsonwr_end_object_value(json);
}

/**************************************************************************
 * Outputs JSON for the motif pwm
 **************************************************************************/
void output_pwm_json(JSONWR_T *json, MOTIF_T *motif) {
  ALPH_T *alph;
  MATRIX_T *freqs;
  int i, j;
  alph = get_motif_alph(motif);
  freqs = get_motif_freqs(motif);
  jsonwr_start_array_value(json);
  for (i = 0; i < get_motif_length(motif); i++) {
    jsonwr_start_array_value(json);
    for (j = 0; j < alph_size_core(alph); j++) {
      jsonwr_dbl_value(json, get_matrix_cell(i, j, freqs));
    }
    jsonwr_end_array_value(json);
  }
  jsonwr_end_array_value(json);
}

/**************************************************************************
 * Outputs JSON for the motif
 **************************************************************************/
void output_motif_json(JSONWR_T *json, int db, MOTIF_T *motif) {
  char *alt, *url;
  double nsites, log_evalue;
  int i, j;
  alt = get_motif_id2(motif);
  nsites = get_motif_nsites(motif);
  log_evalue = get_motif_log_evalue(motif);
  url = get_motif_url(motif);
  jsonwr_start_object_value(json);
  jsonwr_lng_prop(json, "db", db);
  jsonwr_str_prop(json, "id", get_motif_id(motif));
  if (alt[0] != '\0') jsonwr_str_prop(json, "alt", alt);
  jsonwr_lng_prop(json, "len", get_motif_length(motif));
  jsonwr_lng_prop(json, "ltrim", get_motif_trim_left(motif));
  jsonwr_lng_prop(json, "rtrim", get_motif_trim_right(motif));
  if (nsites != 0) jsonwr_dbl_prop(json, "nsites", nsites);
  if (log_evalue > -HUGE_VAL) jsonwr_log10num_prop(json, "evalue", log_evalue, 1);
  if (url != NULL && url[0] != '\0') jsonwr_str_prop(json, "url", url);
  jsonwr_property(json, "pwm");
  output_pwm_json(json, motif);
  jsonwr_end_object_value(json);
}

/**************************************************************************
 * Outputs JSON describing the spacings between a primary and secondary motif.
 **************************************************************************/
void output_spacings_json(JSONWR_T *json, SECONDARY_MOTIF_T* smotif) {
  SPACING_T *spacings;
  SIGSPACE_T *sigspace;
  LINKLST_T *sigseqs;
  LL_LINK_T *node;
  int quad, i, nquads;
  bool revcomp;
  if (smotif->sig_count < 1) {
    fprintf(stderr, "WARNING: Skipping motif because it has no significant spacings.\n");
    fprintf(stderr, "This shouldn't happen so please report.\n");
    return;
  }
  revcomp = alph_has_complement(get_motif_alph(smotif->motif));
  nquads = (revcomp ? NQUADS : 2);
  jsonwr_start_object_value(json);
  jsonwr_lng_prop(json, "idx", smotif->idx);
  jsonwr_property(json, "counts");
  jsonwr_start_array_value(json);
  for (quad = 0; quad < nquads; quad++) {
    spacings = &(smotif->spacings[quad]);
    jsonwr_start_array_value(json);
    for (i = 0; i < spacings->bin_count; i++) {
      jsonwr_lng_value(json, spacings->count[i]);
    }
    jsonwr_end_array_value(json);
  }
  jsonwr_end_array_value(json);
  jsonwr_property(json, "spacings");
  jsonwr_start_array_value(json);
  for (i = 0; i < smotif->sig_count; i++) {
    sigspace = &(smotif->sigs[i]);
    jsonwr_start_object_value(json);
    jsonwr_lng_prop(json, "orient", sigspace->orient);
    jsonwr_lng_prop(json, "bin", sigspace->bin);
    jsonwr_dbl_prop(json, "pvalue", sigspace->pvalue);
    jsonwr_property(json, "inferred_pwm");
    output_pwm_json(json, sigspace->inferred_motif);
    jsonwr_property(json, "alignment_pwm");
    jsonwr_start_array_value(json);
    for (quad = 0; quad < nquads; quad++) {
      if (sigspace->alignment_motif[quad] != NULL) {
        output_pwm_json(json, sigspace->alignment_motif[quad]);
      } else {
        jsonwr_null_value(json);
      }
    }
    jsonwr_end_array_value(json);
    jsonwr_property(json, "seqs");
    jsonwr_start_array_value(json);
    sigseqs = smotif->spacings[sigspace->orient].sequences[sigspace->bin];
    for (node = linklst_first(sigseqs); node != NULL; node = linklst_next(node)) {
      jsonwr_lng_value(json, ((SEQUENCE_T*)linklst_get(node))->idx);
    }
    jsonwr_end_array_value(json);
    jsonwr_end_object_value(json);
  }
  jsonwr_end_array_value(json);
  jsonwr_end_object_value(json);
}

/**************************************************************************
 * Outputs html (containing JSON data) for the results
 **************************************************************************/
void output_html(
  int argc, 
  char **argv, 
  SPAMO_OPTIONS_T *options, 
  int bin_max, 
  ARRAY_T *background,
  SEQUENCE_DB_T *sequence_db,
  RBTREE_T *sequences,
  MOTIF_DB_T *primary_db,
  MOTIF_T *primary_motif, 
  ARRAYLST_T *secondary_dbs,
  RBTREE_T *all_secondary_motifs,
  LINKLST_T *secondary_motifs, 
  time_t start_time, 
  clock_t start_clock
) {
  HTMLWR_T *html;
  JSONWR_T *json;
  RBNODE_T *rbnode;
  LL_LINK_T *node;
  ALPH_T *alph;
  int i;
  time_t end_time;
  clock_t end_clock;
  // setup html monolith writer
  json = NULL;
  bool compress = true;
  if ((html = htmlwr_create(get_meme_data_dir(), TEMPLATE_FILENAME, compress))) {
    htmlwr_set_dest_name(html, options->outdir, HTML_FILENAME);
    htmlwr_replace(html, DATA_FILENAME, "data");
    json = htmlwr_output(html);
    if (json == NULL) die("Template does not contain data section.\n");
  } else {
    DEBUG_MSG(QUIET_VERBOSE, "Failed to open html template file.\n");
    return;
  }
  jsonwr_str_prop(json, "program", "SpaMo");
  jsonwr_str_prop(json, "version", VERSION);
  jsonwr_str_prop(json, "revision", REVISION);
  jsonwr_str_prop(json, "release", ARCHIVE_DATE);
  jsonwr_str_prop(json, "host", hostname());
  jsonwr_time_prop(json, "when", &(options->start));
  jsonwr_args_prop(json, "cmd", argc, argv);
  // output options
  jsonwr_property(json, "options");
  jsonwr_start_object_value(json);
  jsonwr_dbl_prop(json, "seq_min_hit_score", options->score_threshold);
  jsonwr_lng_prop(json, "margin", options->margin);
  jsonwr_lng_prop(json, "bin_size", options->bin);
  jsonwr_lng_prop(json, "bin_pvalue_calc_range", options->test_range);
  jsonwr_lng_prop(json, "use_best_secondary", options->use_best_secondary);
  jsonwr_dbl_prop(json, "seq_max_shared_fract", options->sequence_similarity);
  jsonwr_dbl_prop(json, "seq_odds_ratio", options->odds);
  jsonwr_dbl_prop(json, "bin_pvalue_cutoff", options->pvalue_cutoff);
  jsonwr_dbl_prop(json, "motif_evalue_cutoff", options->motif_evalue_cutoff);
  jsonwr_lng_prop(json, "redundant_overlap", options->overlap);
  jsonwr_dbl_prop(json, "redundant_joint", options->joint);
  jsonwr_dbl_prop(json, "motif_pseudocount", options->pseudocount);
  jsonwr_str_prop(json, "bgfile", options->bg_filename ? options->bg_filename : "");
  jsonwr_dbl_prop(json, "motif_trim", options->trim_bit_threshold);
  jsonwr_lng_prop(json, "xalph", options->xalph);
  jsonwr_dbl_prop(json, "seed", options->prng_seed);
  // Not an input 
  jsonwr_lng_prop(json, "bin_max", bin_max);
  //
  jsonwr_end_object_value(json);
  // output alphabet
  alph = get_motif_alph(primary_motif);
  jsonwr_property(json, "alphabet");
  alph_print_json(alph, json);
  // output background
  jsonwr_property(json, "background");
  jsonwr_start_array_value(json);
  for (i = 0; i < alph_size_core(alph); i++) {
    jsonwr_dbl_value(json, get_array_item(i, background));
  }
  jsonwr_end_array_value(json);
  // output sequence database
  jsonwr_property(json, "sequence_dbs");
  jsonwr_start_array_value(json); // write as list for easy future expansion
  jsonwr_start_object_value(json);
  jsonwr_str_prop(json, "source", sequence_db->source);
  jsonwr_str_prop(json, "name", sequence_db->name);
  jsonwr_time_prop(json, "last_modified", &(sequence_db->last_mod));
  jsonwr_lng_prop(json, "loaded", sequence_db->loaded);
  jsonwr_lng_prop(json, "excluded_too_short", sequence_db->excluded_tooshort);
  jsonwr_lng_prop(json, "excluded_no_match", sequence_db->excluded_nomatch);
  jsonwr_lng_prop(json, "excluded_ambigs", sequence_db->excluded_ambigs);
  jsonwr_lng_prop(json, "excluded_similar", sequence_db->excluded_similar);
  jsonwr_lng_prop(json, "erased_primary_matches", sequence_db->erased_primary_matches);
  jsonwr_end_object_value(json);
  jsonwr_end_array_value(json);
  // output primary motif database
  jsonwr_property(json, "primary_dbs");
  jsonwr_start_array_value(json); // write as a list for easy future expansion
  output_motif_db_json(json, primary_db);
  jsonwr_end_array_value(json);
  // output secondary motif databases
  jsonwr_property(json, "secondary_dbs");
  jsonwr_start_array_value(json);
  for (i = 0; i < arraylst_size(secondary_dbs); ++i) {
    output_motif_db_json(json, (MOTIF_DB_T*)arraylst_get(i, secondary_dbs));
  }
  jsonwr_end_array_value(json);
  // output secondary motifs
  jsonwr_property(json, "secondary_motifs");
  jsonwr_start_array_value(json);
  for (i = 0, rbnode = rbtree_first(all_secondary_motifs); rbnode != NULL; rbnode = rbtree_next(rbnode)) {
    SECONDARY_MOTIF_T *secondary_motif;
    secondary_motif = rbtree_value(rbnode);
    // skip motifs that don't pass the cutoff
    if (!secondary_motif->passes_evalue_cutoff) continue;
    // assign an index to the motif so we can refer to it later
    secondary_motif->idx = i++;
    // output the motif
    output_motif_json(json, secondary_motif->db->id - 1, secondary_motif->motif);
  }
  jsonwr_end_array_value(json);
  // output sequence names
  jsonwr_property(json, "sequence_names");
  jsonwr_start_array_value(json);
  for (i = 0, rbnode = rbtree_first(sequences); rbnode != NULL; rbnode = rbtree_next(rbnode)) {
    SEQUENCE_T *sequence;
    sequence = rbtree_value(rbnode);
    // skip sequences that did not contribute to a significant site
    if (!sequence->contributes) continue;
    // assign an index so we can refer to it later
    sequence->idx = i++;
    // output the sequence name
    jsonwr_str_value(json, sequence->name);
  }
  jsonwr_end_array_value(json);
  // loop over primary motifs (currently one one)
  jsonwr_property(json, "primaries");
  jsonwr_start_array_value(json);
  jsonwr_start_object_value(json);
  jsonwr_property(json, "motif");
  output_motif_json(json, 0, primary_motif);
  // loop over secondary motifs
  jsonwr_property(json, "secondaries");
  jsonwr_start_array_value(json);
  for (node = linklst_first(secondary_motifs); node != NULL; node = linklst_next(node)) {
    GROUPED_MOTIF_T *gmotif;
    LL_LINK_T *node2;
    gmotif = (GROUPED_MOTIF_T*)linklst_get(node);
    jsonwr_start_array_value(json);
    output_spacings_json(json, gmotif->best);
    for (node2 = linklst_first(gmotif->others); node2 != NULL; node2 = linklst_next(node2)) {
      output_spacings_json(json, (SECONDARY_MOTIF_T*)linklst_get(node2));
    }
    jsonwr_end_array_value(json); // cluster list
  }
  jsonwr_end_array_value(json); // secondaries list
  jsonwr_end_object_value(json); // primary object
  jsonwr_end_array_value(json); // primaries list
  // output runtime
  jsonwr_property(json, "run_time");
  jsonwr_start_object_value(json);
  time(&end_time);
  end_clock = clock();
  jsonwr_dbl_prop(json, "cpu", ((double)(end_clock - start_clock)) / CLOCKS_PER_SEC);
  jsonwr_dbl_prop(json, "real", difftime(end_time, start_time));
  jsonwr_end_object_value(json);
  // write rest of HTML file
  if (htmlwr_output(html) != NULL) {
    die("Found another JSON replacement!\n");
  }
  htmlwr_destroy(html);
}

/**************************************************************************
 * Gets the alphabet from the primary motif
 **************************************************************************/
ALPH_T *read_primary_alphabet(char *primary_file) {
  ALPH_T *alph;
  MREAD_T *mread;
  alph = NULL;
  mread = mread_create(primary_file, OPEN_MFILE, true);
  // get the alphabet from the motif file
  if (mread_has_motif(mread)) {
    alph = alph_hold(mread_get_alphabet(mread));
  } else {
    die("\"%s\" does not contain any motifs.\n", primary_file);
  }
  mread_destroy(mread);
  return alph;
}

/**************************************************************************
 * entry point for spamo
 **************************************************************************/
int main(int argc, char **argv) {
  SPAMO_OPTIONS_T options;
  ALPH_T *alph;
  char *trimmed_sequence_data;
  SEQUENCE_DB_T *sequence_db;
  RBTREE_T *sequences;
  MOTIF_DB_T *primary_db;
  MOTIF_T *primary_motif;
  ARRAY_T *background = NULL;
  ARRAYLST_T *secondary_dbs;
  RBTREE_T *secondary_motifs;
  LINKLST_T *sorted_secondary_motifs;
  int bin_max;
  time_t start_time;
  clock_t start_clock;
  time_t end_time;
  clock_t end_clock;
  bool stdin_used = false;

  time(&start_time);
  start_clock = clock();

  char *commandline = process_arguments(argc, argv, &options);

  alph = read_primary_alphabet(options.primary_file);

  // Check that all the motif files have the same alphabet or can
  // be converted to that of the primary, and get it.
  read_motif_alphabets(options.secondary_files, options.xalph, &alph);
  
  srand_mt(options.prng_seed);

  create_spamo_output_directory(&options);

  primary_db = load_primary_motif(&options, alph, &background, &primary_motif, &stdin_used);

  // Add the primary motif to the exclude list if requested.
  if (!options.keep_primary) arraylst_add(get_motif_id(primary_motif), options.exclude_patterns);

  calculate_trimmed_filtered_sequences(&options, alph, primary_motif, background, 
    &trimmed_sequence_data, &sequence_db, &sequences);

  secondary_dbs = load_secondary_motifs(&options, alph, background, &secondary_motifs, &stdin_used);

  load_secondary_matches(&options, sequences, primary_motif, secondary_motifs, background);

  sorted_secondary_motifs = sort_secondary_motifs(secondary_motifs);

  bin_max = calculate_bin_max(sorted_secondary_motifs);

  group_secondary_motifs(&options, sorted_secondary_motifs);

  output_tsv(argc, argv, commandline, &options, bin_max, sequence_db, primary_db, primary_motif, secondary_dbs, 
    sorted_secondary_motifs, rbtree_size(secondary_motifs));
  
  if (! options.text_only) output_html(argc, argv, &options, bin_max, background, sequence_db, 
    sequences, primary_db, primary_motif, secondary_dbs, secondary_motifs, sorted_secondary_motifs, 
    start_time, start_clock);

  //cleanup
  //free the groupings of secondary motifs (but not the secondary motifs)
  linklst_destroy_all(sorted_secondary_motifs, destroy_grouped_motif);
  //free secondary motifs
  rbtree_alter_value_free(secondary_motifs, destroy_secondary_motif);
  rbtree_destroy(secondary_motifs);
  arraylst_destroy(destroy_motif_db, secondary_dbs);
  //free sequences information (not including the actual trimmed sequences)
  rbtree_destroy(sequences);
  destroy_sequence_db(sequence_db);
  //free trimmed sequences
  free(trimmed_sequence_data);
  //don't free primary motif--its in primary_db
  //destroy_motif(primary_motif);
  destroy_motif_db(primary_db);
  free_array(background);
  alph_release(alph);
  //free the command line options information
  cleanup_options(&options);
  free(commandline);

  // print time to screen
  time(&end_time);
  end_clock = clock();
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "run_time cpu= %5.2f real= %5.2f\n", 
      ((double)(end_clock - start_clock)) / CLOCKS_PER_SEC,
      difftime(end_time, start_time));
  }

  return EXIT_SUCCESS;
}

