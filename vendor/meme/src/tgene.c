/********************************************************************
 * FILE: tgene.c
 * AUTHOR: Timothy Bailey
 * CREATE DATE: 26/03/2011
 * PROJECT: MEME suite
 * COPYRIGHT: 2019, UNR
 ********************************************************************/

#define DEFINE_GLOBALS
#define _GNU_SOURCE // for asprintf
#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <libgen.h> // for basename
#include <regex.h>
#include "config.h"
#include "dir.h"
#include "html-monolith.h"
#include "io.h"
#include "bed.h"
#include "bed-io.h"
#include "gtf-io.h"
#include "hash_table.h"
#include "macros.h"
#include "projrel.h"
#include "regress.h"
#include "string-list.h"
#include "utils.h"
#include "regex-utils.h"
#include "string-builder.h"
#include "qvalue.h"
#include "random-variate.h"
#include "mast.h"
#include "linked-list.h"

#include <sys/resource.h>
#include <stdarg.h>
#include <unistd.h>
#include <dirent.h>

#define RCHUNK 1000
#define HISTONE_FILE_NAME_FMT "^.*%s.*(broad|narrow)Peak$"
#define EXPRESSION_FILE_NAME_FMT "^.*%s.*[.]gtf$"
#define MAX_CHRM_LENGTH 1e9
// Number or random permutations for estimating null distribution
#define N_PERMS 10
// Mean of Gaussian noise is minimum/NOISE_FRACTION
#define NOISE_FRACTION 10.0

#define remove_trailing_slash(x) { \
  int len = strlen(x); \
  if ((x)[len-1] == '/') {(x)[len-1] = '\0';} \
}

char *out_dir = "tgene_out";
char *lfile = "";
char *afile = "";
char *tissue = "None";
char *tissues = "";
char *hrd = "";
char *histones = "";
char *mlds = "500000";
char *rna = "";
char *erd = "";
char *ttypes = "protein_coding,processed_transcript";
const double lecat = 0; 	// low expression correlation adjustment threshold
const int seed = 0; 
const double mpv = 0.05;

const char* TEMPLATE_FILENAME = "tgene_template.html";
const char* TSV_FILENAME = "tgene.tsv";		// Name of TSV tgene output
const char* HTML_FILENAME = "tgene.html";	// Name of HTML tgene output

VERBOSE_T verbosity = NORMAL_VERBOSE;

// Structure for tracking tgene command line parameters.
typedef struct options {
  char *locus_file; 		// name of file containing loci
  char *output_dirname; 	// name of the output directory
  bool allow_clobber; 		// allow overwriting of files in output directory
  char *annotation_file;	// name of the file with gene and start site annotation
  char *transcript_types;	// allowed types of transcripts 
  char *max_link_distances;	// comma-separated list of max distance between RE and target
  double max_pvalue; 		// maximum p-value for including non-CT, non-CL links
  char *tissues;		// comma-separated list of tissues
  char *histone_root;		// root directory for histone files
  char *histones;		// comma-separated list of histone mods
  char *rna_source;		// type of RNA expression data
  char *expression_root;	// root directory for expression files
  bool use_gene_ids;		// use gene_ids instead of transcript_ids
  bool no_closest_locus;	// include closest locus only if expression constraints satisfied
  bool no_closest_tss;		// include closest TSS only if expression constraints satisfied
  bool no_noise;		// do not add noise to expression and histone zeros
  double lecat;			// scale correlation if maximum expression of transcript < lecat
  int seed;			// random number seed
  char *desc; 			// description of job
  char *dfile;			// file containing description
  int n_tissues;		// not an option
  STRING_LIST_T *tissue_names;	// not an option
  int n_histones;		// not an option
  STRING_LIST_T *histone_names;	// not an option
  GTF_TYPE rna_source_type;	// not an option
  int *max_distances;           // not an option
  bool **valid_tissues;		// not an option; skip tissue j for histone i if valid_tissues[i][j] false
} TGENE_OPTIONS_T;

// Structure for holding one transcript.
# define get_transcript_chrom(ptr) ((ptr)->annotation->seqname)
# define get_transcript_start(ptr) ((ptr)->annotation->start)
# define get_transcript_end(ptr) ((ptr)->annotation->end)
# define get_transcript_id(ptr) ((ptr)->annotation->transcript_id)
typedef struct transcript {
  GTF_ENTRY_T *annotation;		// GTF entry from annotation file
  double *expr;				// expression values per tissue
  bool valid_type;			// transcript type is allowed
  double max_expr;			// maximum expression across tissues
} TRANSCRIPT_T;

// Structure for table of transcripts.
typedef struct transcript_table {
  int n_transcripts;			// number of transcripts
  int n_tissues;			// number of tissues
  int n_valid;				// number of valid transcripts for current histone
  GTF_ENTRY_T **annotations;		// array of GTF entry pointers
  TRANSCRIPT_T *transcripts;		// array of transcripts, sorted by coords
  bool have_tss_info;			// true if TSS IDs are present
  HASH_TABLE ht_chrom_names;		// known chromosome names
} TRANSCRIPT_TABLE_T;

// Structure for holding one locus.
# define get_locus_chrom(ptr) ((ptr)->reg_el->chrom)
# define get_locus_start(ptr) ((ptr)->reg_el->chrom_start)
# define get_locus_end(ptr) ((ptr)->reg_el->chrom_end)
typedef struct locus {
  BED_ENTRY_T *reg_el;			// BED entry from locus file
  double *histone_level;		// levels per tissue for current histone
  int used_in_link;			// number of links for current histone
  int num_closest_tss_links;		// number of closet TSS links for current histone
  double max_hist;			// maximum current histone level of locus across tissues
  double sd_hist;			// standard deviation of current histone level of locus across tissues
} LOCUS_T;

// Structure for table of loci.
typedef struct locus_table {
  int n_loci;				// number of loci
  int n_tissues;			// number of tissues
  BED_ENTRY_T **reg_els;		// array of BED entry pointers
  LOCUS_T *loci;			// array loci
} LOCUS_TABLE_T;

// Structure for holding links.
#define link_tss_id(ptr) (ptr)->transcript->annotation->transcript_id
#define link_tss_chrom(ptr) (ptr)->transcript->annotation->seqname
#define link_tss_start(ptr) (ptr)->transcript->annotation->start
#define link_tss_end(ptr) (ptr)->transcript->annotation->end
#define link_strand(ptr) (ptr)->transcript->annotation->strand
#define link_gene_id(ptr) (ptr)->transcript->annotation->gene_id
#define link_gene_name(ptr) (ptr)->transcript->annotation->gene_name
#define link_re_chrom(ptr) (ptr)->locus->reg_el->chrom
#define link_re_start(ptr) (ptr)->locus->reg_el->chrom_start
#define link_re_end(ptr) (ptr)->locus->reg_el->chrom_end
#define link_re_width(ptr) MAX(1, link_re_end(ptr) - link_re_start(ptr))
#define link_max_expr(ptr) (ptr)->transcript->max_expr
#define link_max_hist(ptr) (ptr)->locus->max_hist
#define link_sd_hist(ptr) (ptr)->locus->sd_hist
typedef struct link {
  TRANSCRIPT_T *transcript;
  LOCUS_T *locus;
  double correlation;
  double corr_pvalue;		// correlation p-value
  int distance;
  double distance_pvalue;	// Pr(d<=observed); assumes d ~ U[0, max_link_distance]
  double cnd_pvalue;		// p-value of product of correlation and distance p-values
  double qvalue;		// q-value of p-value (Distance or CnD)
  char closest_locus;
  char closest_tss;
  int i_hist;			// histone name index
  double *perm_correlation;
} LINK_T;

// Structure for table of links.
typedef struct link_table {
  int n_links;			// number of links
  LINK_T **links;		// array of link pointers
  bool have_tss_info;		// true if not just gene info
  int n_perms;			// number of null scores / link
} LINK_TABLE_T;

/*************************************************************************
 * Find index of word in list of words separated by whitespace, or return -1.
 *************************************************************************/
static int find_index_of_word_in_list(
  char *name,			// name to print if error
  char *word,			// exact word to match
  char *word_list,		// words separated by separator character
  char sep			// separator character
) {
  int i;

  if (! word_list) return(-1);

  STRING_LIST_T *strings = new_string_list_char_split(sep, word_list);
  int n = get_num_strings(strings);
  for (i=0; i<n; i++) {
    char *string = get_nth_string(i, strings);
    if (strcmp(string, word) == 0) break;
  }
  free_string_list(strings);

  return(i==n ? -1 : i);
} // find_index_of_word_in_list

/*************************************************************************
 * Get listing of files in directory.
 *************************************************************************/
static STRING_LIST_T *list_directory(
  char *dir 			// path to directory
) {
  STRING_LIST_T *listing = new_string_list();
  struct dirent *de;

  DIR *dr = opendir(dir);
  if (dr == NULL) die("Could not open directory %s.\n", dir);
  while ((de = readdir(dr)) != NULL) add_string(de->d_name, listing);
  closedir(dr);

  return(listing);
} // list directory

/*************************************************************************
 * Get paths of files of a given type for a given tissue.
*************************************************************************/
static STRING_LIST_T *get_file_names(
  char *root_dir,	// root directory
  char *tissue,		// name of folder
  char *format,		// format string for file name
  char *type 		// used by format string
) {
  // create space
  int dir_len = strlen(root_dir);
  int tis_len = strlen(tissue);
  int fmt_len = strlen(format);
  int typ_len = strlen(type);
  STR_T *dir = str_create(dir_len + tis_len + 100);
  STR_T *regex = str_create(fmt_len + typ_len + 100);

  // get directory listing for tissue
  str_setf(dir, "%s/%s/", root_dir, tissue);
  STRING_LIST_T *dir_listing = list_directory(str_internal(dir));
  str_setf(regex, format, type);
  STRING_LIST_T *file_names = get_matches_in_string_list(str_internal(regex), dir_listing);

  // prepend directory to all file names
  if (file_names != NULL) prepend_to_strings(str_internal(dir), file_names);

  // space
  free_string_list(dir_listing);
  str_destroy(dir, false);
  str_destroy(regex, false);

  return(file_names);
} // get_file_names

/*****************************************************************************
 * Print a usage message with an optional reason for failure.
 *****************************************************************************/
static void usage(bool msg, char *format, ...) {
  va_list argp;
  char *usage = 
    "\n"
    "Usage: tgene [options] <locus_file> <annotation_file>\n"
    "\n"
    "   Options:\n"
    "    --o <dir>                      output to the specified directory; default: %s\n"
    "    --oc <dir>                     output to the specified directory; default: %s\n"
    "    --transcript-types <ttypes>    comma-separated list of types of transcript to use from annotation file;\n"
    "                                          default: %s\n"
    "    --max-link-distances <mlds>    comma-separated list (no spaces) of maximum distances\n"
    "                                   between an RE and its target; default: %s\n"
    "                                          Note: only one allowed if no histone/expression data given\n"
    "                                          Note: there must be one distance for each histone name in <histones>\n"
    "    --max-pvalue <mp>              maximum p-value for including non-CT, non-CL links in output;\n"
    "                                          default: %.2f\n"
    "    --tissues <tissues>            comma-separated list (no spaces) of tissue names that are the\n"
    "                                   sources of the histone and expression data; default: None\n"
    "    --histone-root <hrd>           histone root directory; default: None\n"
    "                                          Note: histone files must be in subfolders '<hrd>/<t>'\n"
    "                                          where <t> is one of the tissue names in <tissues>\n"
    "    --histones <histones>          comma-separated list (no spaces) of histone names; default: None\n"
    "                                          Note: histone file names must match '<hrd>/<t>/*<hname>*[broad|narrow]Peak'\n"
    "                                          where <t> is one of the tissue names in <tissues>\n"
    "                                          and <hname> is one of the histone names listed in <histones>\n"
    "    --rna-source [Cage|LongPap]    type of expression data in expression files; default: None\n"
    "    --expression-root <erd>        expression root directory; default: None\n"
    "                                          Note: expression files must be in subfolders '<erd>/<t>'\n"
    "                                          where <t> is one of the tissue names in <tissues>, and have\n"
    "                                          extension '.gtf'\n"
    "    --use-gene-ids                 use the 'gene_id' field rather than 'transcript_id' field\n"
    "                                   to associate expression file and annotation file entries;\n"
    "                                          default: use 'transcript_id' field\n"
    "    --lecat <lecat>                scale correlation if maximum expression of transcript < <lecat>;\n"
    "                                          default: %.0f\n"
    "    --no-closest-locus             don't include closest locus for all targets\n"
    "                                   unless constraints are satisfied\n"
    "    --no-closest-tss               don't include closest TSS (target transcript) for all loci\n"
    "                                   unless constraints are satisfied\n"
    "    --no-noise                     do not add noise to expression and histone zeros\n"
    "    --seed <seed>                  seed for random number generator for noise and null model\n"
    "                                          default: %d\n"
    "    --desc <text>                  plain text description of the job\n"
    "    --dfile <file>                 file containing plain text description of the job\n"
    "    --verbosity 1|2|3|4|5          level of information messages output to terminal\n"
    "    --help                         display the usage message and exit\n"
    "    --version                      show version and exit\n"
    ;
  if (format) {
    fprintf(stderr, "\n");
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
  }
  if (msg) {
    fprintf(stderr, usage, out_dir, out_dir, ttypes, mlds, mpv, lecat, seed);
  } else {
    fprintf(stderr, "Type 'tgene --help' for more information.\n");
  }
  fflush(stderr);
  exit(EXIT_FAILURE);
} // usage

/***********************************************************************
 Process command line options and return the command line string.
 ***********************************************************************/
static char *process_command_line(int argc, char* argv[], TGENE_OPTIONS_T *options) {
  bool found_error = false;
  bool have_histone_root = false;
  bool have_expression_root = false;
  bool have_rna_source = false;
  // An easy way of assigning a number to each option.
  // Be careful that there are not more than 62 options as otherwise the value
  // of '?' will be used.
  enum Opts {
    OPT_O, OPT_OC,
    OPT_TRANSCRIPT_TYPES, 
    OPT_MAX_LINK_DISTANCES,
    OPT_MAX_PVALUE, 
    OPT_TISSUES, 
    OPT_HISTONE_ROOT, 
    OPT_HISTONES, 
    OPT_RNA_SOURCE, 
    OPT_EXPRESSION_ROOT, 
    OPT_USE_GENE_IDS, 
    OPT_LECAT,
    OPT_NO_CLOSEST_LOCUS, 
    OPT_NO_CLOSEST_TSS,
    OPT_NO_NOISE,
    OPT_SEED, 
    OPT_DESC, 
    OPT_DFILE,
    OPT_FDESC,
    OPT_VERBOSITY, 
    OPT_HELP, 
    OPT_VERSION
  };
  struct option tgene_options[] = {
    {"o", required_argument, NULL, OPT_O},
    {"oc", required_argument, NULL, OPT_OC},
    {"transcript-types", required_argument, NULL, OPT_TRANSCRIPT_TYPES},
    {"max-link-distances", required_argument, NULL, OPT_MAX_LINK_DISTANCES},
    {"max-pvalue", required_argument, NULL, OPT_MAX_PVALUE},
    {"tissues", required_argument, NULL, OPT_TISSUES},
    {"histone-root", required_argument, NULL, OPT_HISTONE_ROOT},
    {"histones", required_argument, NULL, OPT_HISTONES},
    {"rna-source", required_argument, NULL, OPT_RNA_SOURCE},
    {"expression-root", required_argument, NULL, OPT_EXPRESSION_ROOT},
    {"use-gene-ids", no_argument, NULL, OPT_USE_GENE_IDS},
    {"lecat", required_argument, NULL, OPT_LECAT},
    {"no-closest-locus", no_argument, NULL, OPT_NO_CLOSEST_LOCUS},
    {"no-closest-tss", no_argument, NULL, OPT_NO_CLOSEST_TSS},
    {"no-noise", no_argument, NULL, OPT_NO_NOISE},
    {"seed", required_argument, NULL, OPT_SEED},
    {"desc", required_argument, NULL, OPT_DESC},
    {"dfile", required_argument, NULL, OPT_DFILE},
    {"fdesc", required_argument, NULL, OPT_FDESC},
    {"verbosity", required_argument, NULL, OPT_VERBOSITY},
    {"help", no_argument, NULL, OPT_HELP},
    {"version", no_argument, NULL, OPT_VERSION},
    {NULL, 0, NULL, 0} //boundary indicator
  };
  // Make sure options are set to defaults.
  memset(options, 0, sizeof(TGENE_OPTIONS_T));
  options->output_dirname = out_dir;
  options->allow_clobber = true;
  options->locus_file = lfile;
  options->annotation_file = afile;
  options->transcript_types = ttypes;
  options->max_link_distances = mlds;
  options->max_pvalue = mpv;
  options->tissues = tissues;
  options->histone_root = hrd;
  options->histones = histones;
  options->rna_source = rna;
  options->expression_root = erd;
  options->use_gene_ids = false;
  options->lecat = lecat;
  options->no_closest_locus = false;
  options->no_closest_tss = false;
  options->no_noise = false;
  options->seed = seed;
  options->desc = NULL;
  options->dfile = NULL;
  verbosity = NORMAL_VERBOSE;

  // process arguments
  int i, j;
  while (1) {
    int opt = getopt_long_only(argc, argv, "", tgene_options, NULL);
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
      case OPT_TRANSCRIPT_TYPES:
        options->transcript_types = optarg;
        break;
      case OPT_MAX_LINK_DISTANCES:
        options->max_link_distances = optarg;
        break;
      case OPT_MAX_PVALUE:
	options->max_pvalue = atof(optarg);
        break;
      case OPT_TISSUES:
        options->tissues = optarg;
        break;
      case OPT_HISTONE_ROOT:
        options->histone_root = optarg;
        break;
      case OPT_HISTONES:
        options->histones = optarg;
        break;
      case OPT_RNA_SOURCE:
        options->rna_source = optarg;
        break;
      case OPT_EXPRESSION_ROOT:
        options->expression_root = optarg;
        break;
      case OPT_USE_GENE_IDS:
        options->use_gene_ids = true;
        break;
      case OPT_NO_CLOSEST_LOCUS:
	options->no_closest_locus = true;
        break;
      case OPT_NO_CLOSEST_TSS:
	options->no_closest_tss = true;
        break;
      case OPT_NO_NOISE:
        options->no_noise = true;
        break;
      case OPT_LECAT:
        options->lecat = atof(optarg);
        break;
      case OPT_SEED:
        options->seed = atoi(optarg);
        break;
      case OPT_DESC:
        options->desc = optarg;
        break;
      case OPT_DFILE:
      case OPT_FDESC:
        options->dfile = optarg;
        break;
      case OPT_VERBOSITY:
        verbosity = atoi(optarg);
        break;
      case OPT_HELP:
        usage(true, NULL);
        break;
      case OPT_VERSION:
        fprintf(stdout, VERSION "\n");
        exit(EXIT_SUCCESS);
        break;
      case '?':	// this is what getopt returns on failure
        usage(false, NULL);
        break;
      default: // just in case we forget to handle an option
        die("Unhandled option %d", opt);
    }
  }

  // Get the locus_file and annotation_file.
  if (optind < argc) options->locus_file = argv[optind++];
  if (optind < argc) options->annotation_file = argv[optind++];

  // locus_file
  if (strcmp(options->locus_file, "") == 0) {
    usage(false, "You must specify a BED file with chromosomal loci on the command line.\n"); 
  } else {
    if (! file_exists(options->locus_file)) {
      usage(false, "Locus file '%s' not found.\n"
        "You must specify a BED file with chromosomal loci on the command line.\n",
	options->locus_file);
    }
  }

  // annotation_file
  if (strcmp(options->annotation_file, "") == 0) {
    usage(false, "You must specify a valid annotation file in GTF format.\n");
  } else {
    if (! file_exists(options->annotation_file)) {
      usage(false, "Annotation file '%s' not found.\n"
        "You must specify a valid annotation file in GTF format.\n",
	options->annotation_file);
    }
  }

  // --transcript-types
  //	use defaults

  // --max-link-distances
  STRING_LIST_T *distances = new_string_list_char_split(',', options->max_link_distances);
  int num_distances = get_num_strings(distances);
  options->max_distances = (int *) mm_malloc(sizeof(int) * num_distances);
  for (i=0; i<num_distances; i++) {
    int distance = atoi(get_nth_string(i, distances));
    options->max_distances[i] = distance;
    if (distance <= 1) usage(false, "Illegal value for value-%d in --max-link-distances <mlds>: %s\n"
        "Must be > 0.", i+1, get_nth_string(i, distances));
  }
  free_string_list(distances);

  // --max-pvalue
  if (options->max_pvalue <= 0 || options->max_pvalue > 1) {
    usage(false, "Illegal value for --max-pvalue <mpv>: %f.  Must be >0 and <=1.\n",
      options->max_pvalue);
  }

  // --tissues
  options->n_tissues = 0;
  if (strcmp(options->tissues, "") != 0) {
    STRING_LIST_T *tissue_names = new_string_list_char_split(',', options->tissues);
    // Sort the tissue names so that input order won't matter.
    options->n_tissues = get_num_strings(tissue_names);
    sort_string_list(tissue_names);
    options->tissue_names = tissue_names;
  }

  // --histone-root
  if (strcmp(options->histone_root, "") != 0) {
    remove_trailing_slash(options->histone_root);
    if (!file_directory(options->histone_root)) {
      usage(false, "Histone root directory '%s' does not exist or is not a directory.\n", 
        options->histone_root);
    }
    have_histone_root = true;
  }

  // --histones
  options->n_histones = 0;
  if (strcmp(options->histones, "") != 0) {
    STRING_LIST_T *histone_names = new_string_list_char_split(',', options->histones);
    // Sort the histone names so that input order won't matter.
    options->n_histones = get_num_strings(histone_names);
    sort_string_list(histone_names);
    options->histone_names = histone_names;
  }

  // rna_source
  if (strcmp(options->rna_source, "") != 0) {
    int index = find_index_of_word_in_list("rna_source", options->rna_source, GTF_TYPE_STRINGS, ' ');
    if (index == -1 || index > 1) {		// only first two types allowed; see gtf-io.h
      usage(false, "You specified an invalid RNA source type on the command line (%s).\n", options->rna_source);
    }
    options->rna_source_type = (GTF_TYPE) index;		// strings and types must align!
    have_rna_source = true;
  }

  // --expression-root
  if (strcmp(options->expression_root, "") != 0) {
    remove_trailing_slash(options->expression_root);
    if (!file_directory(options->expression_root)) {
      usage(false, "Expression root directory '%s' does not exist or is not a directory.\n", options->expression_root);
    }
    have_expression_root = true;
  }

  // See if a tissue panel was specified correctly.
  if (! (options->n_tissues>0 || have_histone_root || options->n_histones>0 || have_expression_root || have_rna_source) ) {
    //
    // No tissue panel.  Check things.
    //
    // --max-link-distances
    if (num_distances > 1) {
      usage(false, "You may only specify one maximum link distance with --max-link-distances if there is no tissue panel.\n");
    }
  } else {
    //
    // Have a tissue panel. Check things.
    //
    // --max-link-distances
    if (num_distances != options->n_histones) {
      usage(false, "There must be the same number of (comma-separated) values in --max-link-distances as --histones:\n"
       "  '%s' (%d) vs '%s' (%d)\n.", options->max_link_distances, num_distances,
       options->histones, options->n_histones);
    }
    // --tissues
    if (options->n_tissues == 0) {
      usage(false, "You must specify a (comma-separated) list of tissues with --tissues <tissues>.\n");
    }
    // --histone-root
    if (! have_histone_root) {
      usage(false, "You must specify the root directory containing the histone modification files with --histone-root <hrd>.\n");
    }
    // --histones
    if (options->n_histones == 0) {
      usage(false, "You must specify a valid list of histone names with --histones <histones>.\n");
    }
    // --rna-source
    if (strcmp(options->rna_source, "Cage") != 0 && strcmp(options->rna_source, "LongPap") != 0) {
      usage(false, "You must specify a valid type or RNA expression data ('Cage' or 'LongPap') with --rna-source <type>.\n");
    }
    // --expression-root
    if (! have_expression_root) {
      usage(false, "You must specify the root directory containing the RNA expression files with --expression-root <erd>.\n");
    }

    //
    // Determine which tissues to skip for each histone; do so if expression or modification
    // data is missing for that tissue.
    //
    options->valid_tissues = options->n_histones==0 ? NULL : (bool **) mm_malloc(sizeof(bool *) * options->n_histones);
    int histone_root_len = strlen(options->histone_root);
    char *eq = "====================";
    for (i=0; i<options->n_histones; i++) {
      char *histone = get_nth_string(i, options->histone_names);
      int n_both = 0;
      int n_histone_found = 0;
      int n_expression_found = 0;
      options->valid_tissues[i] = (bool *) mm_malloc(sizeof(bool) * options->n_tissues);
      DEBUG_FMT(NORMAL_VERBOSE,
	"Checking for histone (%s) and expression data for the given tissues...\n\n"
	"%-20s %10s %10s %10s %10s\n" 
	"%-20.20s %10.10s %10.10s %10.10s %10.10s\n",
	histone, "TISSUE", "HISTONE", "EXPRESSION", "BOTH", "N_BOTH", eq, eq, eq, eq, eq);
      for (j=0; j<options->n_tissues; j++) {
	char *tissue = get_nth_string(j, options->tissue_names);

	// see if there is a histone file for the current tissue
	STRING_LIST_T *file_names = get_file_names(options->histone_root, tissue, HISTONE_FILE_NAME_FMT, histone);
	bool found_histone = (file_names && get_num_strings(file_names));
	free_string_list(file_names);

	// see if there is an expression file for the current tissue
	file_names = get_file_names(options->expression_root, tissue, EXPRESSION_FILE_NAME_FMT, options->rna_source);
	bool found_expression = (file_names && get_num_strings(file_names));
	free_string_list(file_names);

	bool found_both = (found_histone && found_expression);
	if (found_both) {
	  n_both++;
	  options->valid_tissues[i][j] = true;
	} else {
	  options->valid_tissues[i][j] = false;
	}
	if (found_histone) n_histone_found++;
	if (found_expression) n_expression_found++;
	DEBUG_FMT(NORMAL_VERBOSE, "%-20s %10s %10s %10s %10d\n",
	  tissue, found_histone ? "yes" : "no", 
	  found_expression ? "yes" : "no", found_both ? "yes" : "no", 
	  n_both);
      } // tissue
      DEBUG_MSG(NORMAL_VERBOSE, "\n");
      // check that there are enough tissues with both expression and histone data
      if (n_histone_found == 0) {
	fprintf(stderr,  
	    "No histone files found.\n"
	    "Based on your command line, histone files must be in the path\n\t'%s/<t>/*%s*[broad|narrow]Peak*'\n"
	    "where <t> is one of the tissues in the list\n\t'%s'.\n",
	    options->histone_root, histone, options->tissues
	  );
	found_error = true;
      }
      if (n_expression_found == 0) {
	fprintf(stderr,  
	  "No expression files found.\n"
	  "Based on your command line, expression files must be in the path\n\t'%s/<t>/*%s*.gtf',\n"
	  "where <t> is one of the tissues in the list\n\t'%s'.\n",
	  options->expression_root, options->rna_source,
	  options->tissues
	);
	found_error = true;
      } 
      if (found_error) {
	fprintf(stderr, "\n"); 
	usage(false, NULL);
      }
    } // histone
  } // Tissue panel specified.

  // return the command line in a string
  return(get_command_line(argc, argv));

} // process_command_line

/*************************************************************************
 * Replace commas in a string with blanks.
 *************************************************************************/
static char *remove_commas(
  char *string
) {
  char *new_string = strdup(string);
  char *c;
  for (c=new_string; *c != '\0'; c++) if (*c == ',') *c = ' ';
  return(new_string);
} // remove_commas

/*************************************************************************
 * Create a string with a locus.
 *************************************************************************/
static char *create_locus_string(
  char *chrom,
  long start,
  long end
)
{
  int len = strlen(chrom) + 1 + log(start+1)/log(10) + 1 + log(end+1)/log(10) + 2;
  char *locus = mm_malloc(len+1);
  snprintf(locus, len+1, "%s:%ld-%ld", chrom, start, end); 
  return(locus);
} // create_locus_string

/*************************************************************************
 * Output the HTML file.
 *************************************************************************/
static void output_html_file(
  TGENE_OPTIONS_T* options, 
  int argc, 
  char** argv, 
  LINK_TABLE_T *link_table
) {
  int i;
  HTMLWR_T *html;
  JSONWR_T *json;
  int n_histones = options->n_histones;

  // Setup HTML monolith writer.
  json = NULL;
  if ((html = htmlwr_create(get_meme_data_dir(), TEMPLATE_FILENAME, false))) {
    htmlwr_set_dest_name(html, options->output_dirname, HTML_FILENAME);
    htmlwr_replace(html, "tgene_data.js", "data");
    json = htmlwr_output(html);
    if (json == NULL) die("Template does not contain data section.\n");
  } else {
    DEBUG_MSG(QUIET_VERBOSE, "Failed to open html template file\n");
    return;
  }

  // Now output some JSON.
  // Output some top level variables.
  jsonwr_str_prop(json, "version", VERSION);
  jsonwr_str_prop(json, "revision", REVISION);
  jsonwr_str_prop(json, "release", ARCHIVE_DATE);
  jsonwr_str_prop(json, "program", "T-Gene");
  // Output options and command.
  jsonwr_args_prop(json, "cmd", argc, argv);
  jsonwr_property(json, "options");
  jsonwr_start_object_value(json);
  jsonwr_str_prop(json, "locus_file", options->locus_file);
  jsonwr_str_prop(json, "annotation_file", options->annotation_file);
  jsonwr_str_prop(json, "transcript_types", remove_commas(options->transcript_types));
  jsonwr_str_prop(json, "max_link_distances", remove_commas(options->max_link_distances));
  jsonwr_dbl_prop(json, "max_pvalue", options->max_pvalue);
  jsonwr_str_prop(json, "tissues", options->n_tissues ? combine_string_list(options->tissue_names, " ") : "");
  jsonwr_str_prop(json, "histone_root", options->histone_root);
  jsonwr_str_prop(json, "histones", options->n_histones ? combine_string_list(options->histone_names, " ") : "");
  jsonwr_str_prop(json, "rna_source", options->rna_source);
  jsonwr_str_prop(json, "expression_root", options->expression_root);
  jsonwr_str_prop(json, "use_gene_ids", (options->use_gene_ids ? "true" : "false"));
  jsonwr_lng_prop(json, "lecat", options->lecat);
  jsonwr_bool_prop(json, "inc_closest_locus", ! options->no_closest_locus);
  jsonwr_bool_prop(json, "inc_closest_tss", ! options->no_closest_tss);
  jsonwr_bool_prop(json, "noise", ! options->no_noise);
  jsonwr_lng_prop(json, "seed", options->seed);
  jsonwr_end_object_value(json);

  // Output description
  jsonwr_desc_prop(json, "job_description", options->dfile, options->desc);

  // Output general information.
  jsonwr_bool_prop(json, "have_tss_info", link_table->have_tss_info);
  jsonwr_lng_prop(json, "n_perms", N_PERMS);
  jsonwr_dbl_prop(json, "noise_fraction", NOISE_FRACTION);

  // Output Links.
  jsonwr_property(json, "regulatory_links");
  jsonwr_start_array_value(json);
  int n_links = link_table->n_links;
  int n_sig_links = 0;
  bool closest_locus = (options->no_closest_locus == false);
  bool closest_tss = (options->no_closest_tss == false);
  for (i=0; i<n_links; i++) {
    LINK_T *link = link_table->links[i];
    // Check if link passes p-value test unless we are outputing CT or CL links.
    if (
       (link->cnd_pvalue <= options->max_pvalue) ||		// passes p-value test
       (closest_locus && link->closest_locus == 'T') || 	// output a closest locus
       (closest_tss && link->closest_tss == 'T') 		// output a closest TSS 
    ) {
      if (link->cnd_pvalue <= options->max_pvalue) n_sig_links++;

      char *tss_locus = create_locus_string(link_tss_chrom(link), link_tss_start(link), link_tss_end(link));
      char *re_locus = create_locus_string(link_re_chrom(link), link_re_start(link), link_re_end(link));

      jsonwr_start_object_value(json);

      jsonwr_str_prop(json, "gene_id", link_gene_id(link));
      jsonwr_str_prop(json, "gene_name", link_gene_name(link));
      jsonwr_str_prop(json, "tss_id", link_tss_id(link));
      jsonwr_str_prop(json, "tss_locus", tss_locus);
      jsonwr_strn_prop(json, "strand", &link_strand(link), 1);
      if (n_histones > 0) jsonwr_dbl_prop(json, "max_expr", link_max_expr(link));
      jsonwr_str_prop(json, "re_locus", re_locus);
      if (n_histones > 0) jsonwr_dbl_prop(json, "max_hist", link_max_hist(link));
      jsonwr_lng_prop(json, "distance", link->distance);
      jsonwr_strn_prop(json, "closest_locus", &link->closest_locus, 1);
      jsonwr_strn_prop(json, "closest_tss", &link->closest_tss, 1);
      if (n_histones > 0) jsonwr_str_prop(json, "histone", get_nth_string(link->i_hist, options->histone_names));
      if (n_histones > 0) jsonwr_dbl_prop(json, "correlation", link->correlation);
      if (n_histones > 0) jsonwr_dbl_prop(json, "corr_pvalue", link->corr_pvalue);
      jsonwr_dbl_prop(json, "dist_pvalue", link->distance_pvalue);
      if (n_histones > 0) jsonwr_dbl_prop(json, "cnd_pvalue", link->cnd_pvalue);
      jsonwr_dbl_prop(json, "qvalue", link->qvalue);

      jsonwr_end_object_value(json);
      
      free(tss_locus);
      free(re_locus);
    }
  }

  DEBUG_FMT(NORMAL_VERBOSE, "Found %d links passing the p-value threshold (p-value <= %.2f) out of %d links (%.2f%%).\n",
    n_sig_links, options->max_pvalue, n_links, 100*(n_sig_links/(double)n_links));

  // Finish writing the links.
  jsonwr_end_array_value(json);

  // Finish writing the HTML file.
  if (html) {
    if (htmlwr_output(html) != NULL) {
      die("Found another JSON replacement!\n");
    }
    htmlwr_destroy(html);
  }

} // output_html_file

/*************************************************************************
 * Finish up JSON output and HTML output
 *************************************************************************/
static void end_json(HTMLWR_T* html, JSONWR_T* json) {
  // finish writing motifs
  if (json) jsonwr_end_array_value(json);
  // finish writing html file
  if (html) {
    if (htmlwr_output(html) != NULL) {
      die("Found another JSON replacement!\n");
    }
    htmlwr_destroy(html);
  }
} // end_json

/*************************************************************************
 * Compare two LOCUS_T objects by genomic coordinate.
 *************************************************************************/
static int compare_locus_coords(
 const void *p1,
 const void *p2
) {
  const LOCUS_T *item1 = (LOCUS_T *) p1;
  const LOCUS_T *item2 = (LOCUS_T *) p2;
  int cmp;

  // Compare chromosome names.
  cmp = strcmp(get_locus_chrom(item1), get_locus_chrom(item2));
  if (cmp != 0) return(cmp);

  // Compare locus starts.
  cmp = get_locus_start(item1) - get_locus_start(item2);
  if (cmp != 0) return(cmp);

  // Compare locus ends.
  cmp = get_locus_end(item1) - get_locus_end(item2);
  if (cmp != 0) return(cmp);

  die("Result of compare_locus_coords() was ambiguous.\n");
  return(cmp);

} // compare_locus_coords

/*************************************************************************
 * Compare two TRANSCRIPT_T objects by genomic coordinate.
 *************************************************************************/
static int compare_transcript_coords(
 const void *p1,
 const void *p2
) {
  const TRANSCRIPT_T *item1 = (TRANSCRIPT_T *) p1;
  const TRANSCRIPT_T *item2 = (TRANSCRIPT_T *) p2;
  int cmp;

  // Compare chromosome names.
  cmp = strcmp(get_transcript_chrom(item1), get_transcript_chrom(item2));
  if (cmp != 0) return(cmp);

  // Compare starts.
  cmp = get_transcript_start(item1) - get_transcript_start(item2);
  if (cmp != 0) return(cmp);

  // Compare ends.
  cmp = get_transcript_end(item1) - get_transcript_end(item2);
  if (cmp != 0) return(cmp);

  // Compare IDs.
  cmp = strcmp(get_transcript_id(item1), get_transcript_id(item2));
  if (cmp != 0) return(cmp);

  die("Result of compare_transcript_coords() was ambiguous.\n");

  return(cmp);
} // compare_transcript_coords

/*************************************************************************
 * Compare two LINK_T objects by fabs(correlation), then TSS_ID then RE_locus then histone.
 *************************************************************************/
static int compare_abs_correlation(
 const void *p1,
 const void *p2
) {
  const LINK_T *item1 = *(LINK_T **) p1;
  const LINK_T *item2 = *(LINK_T **) p2;
  int cmp;

  // Compare the |corelation| of two links.
  if (fabs(item1->correlation) > fabs(item2->correlation)) {
    return(-1);
  } else if (fabs(item1->correlation) < fabs(item2->correlation)) {
    return(+1);
  }

  // Compare the TSS IDs of two links.
  cmp = strcmp(link_tss_id(item1), link_tss_id(item2));
  if (cmp != 0) return(cmp);

  // Compare the locus chromosome names
  cmp = strcmp(link_re_chrom(item1), link_re_chrom(item2));
  if (cmp != 0) return(cmp);

  // Compare locus starts.
  cmp = link_re_start(item1) - link_re_start(item2);
  if (cmp != 0) return(cmp);

  // Compare locus ends.
  cmp = link_re_end(item1) - link_re_end(item2);
  if (cmp != 0) return(cmp);

  // Compare histone indices.
  cmp = item1->i_hist - item2->i_hist;
  if (cmp != 0) return(cmp);

  die("Result of compare_abs_correlation() was ambiguous.\n");

  return(cmp);
} // compare_abs_correlation

/*************************************************************************
 * Compare two LINK_T objects by TSS_ID
 * then RE_locus then histone.
 *************************************************************************/
inline static int compare_tss_ids(
 const void *p1,
 const void *p2
) {
  const LINK_T *item1 = *(LINK_T **) p1;
  const LINK_T *item2 = *(LINK_T **) p2;
  int cmp;

  // Compare the TSS IDs of two links.
  cmp = strcmp(link_tss_id(item1), link_tss_id(item2));
  if (cmp != 0) return(cmp);

  // Compare the locus chromosome names
  cmp = strcmp(link_re_chrom(item1), link_re_chrom(item2));
  if (cmp != 0) return(cmp);

  // Compare locus starts.
  cmp = link_re_start(item1) - link_re_start(item2);
  if (cmp != 0) return(cmp);

  // Compare locus ends.
  cmp = link_re_end(item1) - link_re_end(item2);
  if (cmp != 0) return(cmp);

  // Compare histone indices.
  cmp = item1->i_hist - item2->i_hist;

  return(cmp);
} // compare_tss_ids

/*************************************************************************
 * Compare two LINK_T objects by correlation p-value, 
 * then TSS_ID then RE_locus then histone.
 *************************************************************************/
static int compare_corr_pvalues(
 const void *p1,
 const void *p2
) {
  const LINK_T *item1 = *(LINK_T **) p1;
  const LINK_T *item2 = *(LINK_T **) p2;
  int cmp;

  // Compare the correlation p-values of two links.
  if (item1->corr_pvalue < item2->corr_pvalue) {
    return(-1);
  } else if (item1->corr_pvalue > item2->corr_pvalue) {
    return(+1);
  }

  // Compare the TSS IDs etc...
  cmp = compare_tss_ids(p1, p2);
  if (cmp != 0) return(cmp);

  die("Result of compare_corr_pvalues() was ambiguous.\n");

  return(cmp);
} // compare_corr_pvalues

/*************************************************************************
 * Compare two LINK_T objects by CnD p-value, 
 * then TSS_ID then RE_locus then histone.
 *************************************************************************/
static int compare_cnd_pvalues(
 const void *p1,
 const void *p2
) {
  const LINK_T *item1 = *(LINK_T **) p1;
  const LINK_T *item2 = *(LINK_T **) p2;
  int cmp;

  // Compare the CnD p-values of two links.
  if (item1->cnd_pvalue < item2->cnd_pvalue) {
    return(-1);
  } else if (item1->cnd_pvalue > item2->cnd_pvalue) {
    return(+1);
  }

  // Compare the TSS IDs etc...
  cmp = compare_tss_ids(p1, p2);
  if (cmp != 0) return(cmp);

  die("Result of compare_cnd_pvalues() was ambiguous.\n");

  return(cmp);
} // compare_cnd_pvalues

/*************************************************************************
 * Compute the distance between a LOCUS_T object and a TRANSCRIPT_T object.
 * Returns:
 *	-(m+1) 	if locus chromosome < transcript chromosome
 *	+(m+1) 	if locus chromosome > transcript chromosome
 *	0	if they overlap
 *	-d	if same chromosome and locus starts upstream of transcript
 *	+d 	if same chromosome and locus starts at or downstream of transcript
 *************************************************************************/
static long locus_to_trans_distance(
 const LOCUS_T *locus,
 const TRANSCRIPT_T *trans,
 const long m			// maximum distance
) {

  // Compare chromosome names.
  int cmp = strcmp(get_locus_chrom(locus), get_transcript_chrom(trans));
  if (cmp < 0) {
    return -(m+1);
  } else if (cmp > 0) {
    return +(m+1);
  }

  // Compute distance.
  long l_start = get_locus_start(locus);
  long t_start = get_transcript_start(trans);
  if (l_start < t_start) {			// locus is upstream of transcript
    return -MAX(0, t_start - get_locus_end(locus));
  } else {
    return MAX(0, l_start - get_transcript_end(trans));
  }

} // locus_to_trans_distance

/*************************************************************************
 * Determine if a LOCUS_T object and overlaps a BED_ENTRY object.
 * Returns:
 *	-2 	if locus chromosome < bed_entry chromosome
 *	+2 	if locus chromosome > bed_entry chromosome
 *	-1	if same chromosome and locus is upstream of bed_entry
 *	0	if same chromosome and locus overlaps bed_entry
 *	+1 	if same chromosome and locus is downstream of bed_entry
 *************************************************************************/
static int locus_to_bed_overlap(
 const LOCUS_T *locus,
 const BED_ENTRY_T *bed_entry
) {
  // Compare chromosome names.
  int cmp = strcmp(get_locus_chrom(locus), get_bed_chrom(bed_entry));
  if (cmp < 0) {
    return -2;
  } else if (cmp > 0) {
    return +2;
  }

  // Check overlap.
  long l_start = get_locus_start(locus);
  long b_start = get_bed_start(bed_entry);
  if (l_start < b_start) {				// locus starts upstream of bed_entry
    return get_locus_end(locus) < b_start ? -1 : 0;
  } else {
    return get_bed_end(bed_entry) < l_start ? +1 : 0;
  }
} // locus_to_bed_overlap

/*************************************************************************
 * Read the locus file and initialize the table containing the loci.
 * Does not fill in the table with histone level.
 *************************************************************************/
static LOCUS_TABLE_T *initialize_locus_table(
  TGENE_OPTIONS_T *options
) {
  int i;

  FILE *locus_file;
  if ((locus_file = fopen(options->locus_file, "r")) == NULL) {
    die("Unable to open locus file \"%s\" for reading.\n", options->locus_file);
  }

  // Read in the BED file of regulatory elements.
  DEBUG_FMT(NORMAL_VERBOSE, "Reading locus file %s of putative regulatory elements.\n", options->locus_file);
  BED_ENTRY_T **reg_els = NULL;
  // Convert coordinates to 1-based, closed.
  int n_loci = read_bed_file(locus_file, true, BED3, &reg_els);
  if (n_loci <= 0) {
    fprintf(stderr, "\nNo valid entries found in BED locus file '%s'.\n", options->locus_file);
    fprintf(stderr, "\nNo output was generated.\n");
    exit(1);
  }
  fclose(locus_file);
  DEBUG_FMT(NORMAL_VERBOSE, "Locus file contains %d entries.\n", n_loci);

  // Sort the regulatory elements by genomic start.
  qsort(reg_els, n_loci, sizeof(BED_ENTRY_T *), compare_bed_coords);

  //
  // Create a hash table of the loci and remove duplicates.
  //
  char *key;
  HASH_TABLE ht_loci = n_loci > 0 ? hash_create(n_loci, NULL) : NULL;
  int n_valid;
  for (i=n_valid=0; i<n_loci; i++) {
    BED_ENTRY_T *reg_el = reg_els[i];
    BED_ENTRY_T *prev_reg_el = (i==0) ? NULL : reg_els[i-1];
    int dummy=asprintf(&key, "%s:%ld-%ld", get_bed_chrom(reg_el), get_bed_start(reg_el), get_bed_end(reg_el));
    HASH_TABLE_ENTRY *hte = hash_lookup_str(key, ht_loci);
    if (hte) {			// key already in hash table
      DEBUG_FMT(NORMAL_VERBOSE, "Removing duplicate regulatory element locus: %s\n", key);
    } else {
      hash_insert_str_value(key, reg_el, ht_loci);	// value for this key is the reg_el
      if (n_valid < i) reg_els[n_valid] = reg_els[i];
      n_valid++;
    }
    free(key);
  } // reg_el

  hash_destroy(ht_loci);
  DEBUG_FMT(NORMAL_VERBOSE, "Locus table contains %d valid entries (%d duplicates removed).\n", n_valid, n_loci-n_valid);
  n_loci = n_valid;

  // Create array of loci objects (element and its histone state).
  int n_tissues = options->n_tissues;
  LOCUS_T *loci = mm_malloc(sizeof(LOCUS_T) * n_loci);
  for (i=0; i<n_loci; i++) {
    loci[i].reg_el = reg_els[i];		// pointer to the BED entry for this locus
    loci[i].histone_level = mm_calloc(n_tissues, sizeof(double));	// will hold the histone level per tissue
  }

  // Fill in the locus table and return it.
  LOCUS_TABLE_T *locus_table = mm_malloc(sizeof(LOCUS_TABLE_T));
  locus_table->n_loci = n_loci;
  locus_table->n_tissues = n_tissues;
  locus_table->reg_els = reg_els;
  locus_table->loci = loci;		// sorted by coord

  return(locus_table);

} // initialize_locus_table

/*************************************************************************
 * Create a link object, add it to the end of the link table and 
 * initialize it from a locus and a transcript with histone name and
 * link length (distance).
 * Do not free elements because they are all pointers to other object elements.
 *************************************************************************/
static LINK_T *create_link(
  LINK_TABLE_T *link_table,
  TRANSCRIPT_T *transcript,
  LOCUS_T *locus,
  int i_hist,
  int distance
) {
  // Add space to link table if necessary.
  if ((link_table->n_links % RCHUNK) == 0) {
    Resize(link_table->links, link_table->n_links+RCHUNK, LINK_T *);
  }
  // Create link and add it to the link table.
  LINK_T *link = (LINK_T *) mm_malloc(sizeof(LINK_T));
  link_table->links[link_table->n_links++] = link;

  // Initialize the link.
  link->transcript = transcript;
  link->locus = locus;
  link->correlation = 0;
  link->corr_pvalue = 1;
  link->distance_pvalue = 1;
  link->cnd_pvalue = 1;
  link->qvalue = 1;
  link->distance = link_strand(link) == '+' ? distance : -distance;
  link->closest_locus = 'F';
  link->closest_tss = 'F';
  link->i_hist = i_hist;
  link->perm_correlation = (double *) mm_malloc(sizeof(double) * link_table->n_perms);
  return(link);
} // create_link

/*************************************************************************
 * Finish the link table.
*************************************************************************/
static void finish_link_table(
  TGENE_OPTIONS_T* options, 
  LINK_TABLE_T *link_table
) {
  int i, j, i_hist;
  int n_links = link_table->n_links;
  int n_perms = link_table->n_perms;
  LINK_T **links = link_table->links;
  int n_histones = options->n_histones;
  HASH_TABLE_ENTRY *hte;

  // Check if any links were found.
  if (n_links == 0) {
    fprintf(stderr, "\nNo links were possible.\n");
    fprintf(stderr, "\nNo output was generated.\n");
    exit(1);
  }

  // Sort the link table by fabs(correlation) so that it
  // will be in the same order as the scores in convert_scores_to_pvalues. 
  if (n_links > 0) qsort(links, n_links, sizeof(LINK_T *), compare_abs_correlation);

  //
  // Get the p-values and the q-value of the link score.
  //
  ARRAY_T *scores = NULL;
  DEBUG_FMT(NORMAL_VERBOSE, "  Computing p-values and q-values for %d links\n", n_links);
  if (n_histones == 0) {
    // No tissue panel.
    for (i=0; i<n_links; i++) {
      LINK_T *link = links[i];
      // Compute distance "p-value" and use it for CnD p-value as well.
      int max_distance = options->max_distances[0];
      int d = abs(link->distance);		// absolute value of distance from edge of RE to TSS
      int w = link_re_width(link);		// width or RE
      link->correlation = 0;
      link->corr_pvalue = 1;
      link->distance_pvalue = MIN(1, (2.0*d + w) / (2.0*max_distance + w));
      link->cnd_pvalue = link->distance_pvalue;
    } // link
    DEBUG_FMT(NORMAL_VERBOSE, "Estimated %d Distance p-values.\n", n_links);
  } else {
    // Have a tissue panel.
    ARRAY_T *null_scores = NULL;
    for (i_hist=0; i_hist<n_histones; i_hist++) {
      scores = resize_array(scores, n_links);
      null_scores = resize_array(null_scores, n_links*n_perms);
      int n_scores = 0;
      int n_null = 0;
      char *histone = get_nth_string(i_hist, options->histone_names);
      // Copy the scores and null scores into the arrays for the current histone
      // where "score" is the absolute value of the correlation.
      for (i=0; i<n_links; i++) {
	LINK_T *link = links[i];
	if (link->i_hist == i_hist) {
	  set_array_item(n_scores++, fabs(link->correlation), scores);
	  for (j=0; j<n_perms; j++) {
	    set_array_item(n_null++, fabs(link->perm_correlation[j]), null_scores);
	  }
	}
      }
      DEBUG_FMT(NORMAL_VERBOSE, "    Computed %d link correlations for histone %s.\n", n_scores, histone);
      DEBUG_FMT(NORMAL_VERBOSE, "    Computed %d null link correlations for histone %s.\n", n_null, histone);

      scores = resize_array(scores, n_scores);
      null_scores = resize_array(null_scores, n_null);
      // Convert the scores to p-values using the null scores.
      // This will sort the score array, but it is already sorted so that is OK.
      convert_scores_to_pvalues(false, scores, null_scores);
      DEBUG_MSG(NORMAL_VERBOSE, "(Scores are the absolute values of the link correlations.)\n");
      // Save the p-value, distance_pvalue and cnd_pvalue for each link.
      int max_distance = options->max_distances[i_hist];
      int i_score = 0;
      for (i=0; i<n_links; i++) {
	LINK_T *link = links[i];
	if (link->i_hist == i_hist) {
	  link->corr_pvalue = get_array_item(i_score, scores);
	  // Compute distance "p-value" and the p-value of the product of p-values.
	  int d = abs(link->distance);		// absolute value of distance from edge of RE to TSS
	  int w = link_re_width(link);		// width or RE
	  link->distance_pvalue = MIN(1, (2.0*d + w) / (2.0*max_distance + w));
	  link->cnd_pvalue = qfast(2, link->corr_pvalue * link->distance_pvalue);
	  i_score++;
	}
      }
      DEBUG_FMT(NORMAL_VERBOSE, "Estimated %d p-values for histone %s.\n", i_score, histone);
    } // histone
  } // tissue panel

  // Convert the p-values to q-values.
  DEBUG_FMT(NORMAL_VERBOSE, "Converting %s p-values to q-values.\n", n_histones==0 ? "Distance" : "CnD");
  // Sort the link table by CnD p-value (will be = Distance p-values if no tissue panel).
  if (n_links > 0) qsort(links, n_links, sizeof(LINK_T *), compare_cnd_pvalues);
  // Copy the CnD p-values into the scores array (will be = Distance p-values if no tissue panel).
  scores = resize_array(scores, n_links);
  for (i=0; i<n_links; i++) {
    LINK_T *link = links[i];
    set_array_item(i, link->cnd_pvalue, scores);
  }
  compute_qvalues(
    false, 			// Don't stop with FDR.
    true, 			// Estimate pi_zero.
    NULL, 			// Don't store pi_zero in a file.
    NUM_BOOTSTRAPS,
    NUM_BOOTSTRAP_SAMPLES,
    NUM_LAMBDA,
    MAX_LAMBDA,
    n_links,
    scores,			// p-values
    NULL			// No sampled p-values provided
  );
  // Save the qvalues.
  for (i=0; i<n_links; i++) {
    LINK_T *link = links[i];
    link->qvalue = get_array_item(i, scores);
  }

} // finish_link_table
  
/*************************************************************************
 * Print the command line to a TSV file.
*************************************************************************/
static void print_command_line(
  FILE *outfile,
  char *commandline
) {
  fprintf(outfile, 
    "\n# T-Gene (Prediction of Target Genes): Version " VERSION " compiled on " __DATE__ " at " __TIME__ " \n"
    "# The format of this file is described in " SITE_URL "/doc/tgene-output-format.html.\n"
    "# %s\n", commandline);
} // print_command_line

/*************************************************************************
 * Print link table.
*************************************************************************/
static void output_links_tsv_file(
  TGENE_OPTIONS_T *options,
  char *commandline,
  LINK_TABLE_T *link_table
) {
  int i;
  int n_links = link_table->n_links;
  LINK_T **links = link_table->links;

  // Sort the links by score.
  if (n_links > 0) qsort(links, n_links, sizeof(LINK_T *), compare_corr_pvalues);

  // Open the file.
  int dir_len = strlen(options->output_dirname);
  STR_T *file_name = str_create(dir_len + 100);
  str_setf(file_name, "%s/%s", options->output_dirname, "links.tsv");
  FILE *outfile;
  if ((outfile = fopen(str_internal(file_name), "w")) == NULL) {
    die("Unable to open links TSV file \"%s\" for writing.\n", file_name);
  }
  str_destroy(file_name, false);

  // Print header.
  fprintf(outfile, 
    "Gene_ID\tGene_Name\t"
    "TSS_ID\tTSS_Locus\tStrand\t"
    "Max_Expr\t"
    "RE_Locus\t"
    "Max_Hist\t"
    "Distance\t"
    "Closest_Locus\t"
    "Closest_TSS\t"
    "Histone\t"
    "Correlation\t"
    "Correlation_P_Value\t"
    "Distance_P_Value\t"
    "CnD_P_Value\t"
    "Q_Value"
    "\n");

  // Print lines.
  bool closest_locus = (options->no_closest_locus == false);
  bool closest_tss = (options->no_closest_tss == false);
  for (i=0; i<n_links; i++) {
    LINK_T *link = links[i];
    // Check if link passes p-value test unless we are outputing CT or CL links.
    if (
       (link->cnd_pvalue <= options->max_pvalue) ||		// passes p-value test
       (closest_locus && link->closest_locus == 'T') || 	// output a closest locus
       (closest_tss && link->closest_tss == 'T') 		// output a closest TSS 
    ) {
      fprintf(outfile,
	"%s\t"		// gene_id
	"%s\t"		// gene_name
	"%s\t"		// tss_id
	"%s:%ld-%ld\t"	// tss_chrom:tss_start-tss_end
	"%c\t"		// strand
	"%0.2f\t"	// max_expr
	"%s:%ld-%ld\t"	// re_chrom:re_start-re_end
	"%0.2f\t"	// max_hist
	"%d\t"		// distance
	"%c\t"		// closest_locus
	"%c\t"		// closest_tss
	"%s\t"		// histone
	"%0.8f\t"	// correlation
	"%0.3e\t" 	// corr_pvalue
	"%0.3e\t"	// distance_pvalue
	"%0.3e\t"	// cnd_pvalue
	"%0.3e"		// qvalue
	"\n",
	link_gene_id(link), 
	link_gene_name(link),
	link_tss_id(link), 
	link_tss_chrom(link), link_tss_start(link), link_tss_end(link),  // tss_locus
	link_strand(link),
	link_max_expr(link),
	link_re_chrom(link), link_re_start(link), link_re_end(link),	// re_locus
	link_max_hist(link),
	link->distance,
	link->closest_locus, 
	link->closest_tss,
	options->n_histones == 0 ? "" :get_nth_string(link->i_hist, options->histone_names),	// histone
	link->correlation, 
	link->corr_pvalue,
	link->distance_pvalue,
	link->cnd_pvalue,
	link->qvalue
      );
    }
  } // link
  print_command_line(outfile, commandline);
  fclose(outfile);
} // output_links_tsv_file

/*************************************************************************
 * Add links for the current histone modification to the link table.
 *************************************************************************/
static int update_link_table(
  TGENE_OPTIONS_T *options,
  int i_hist,			// current histone index; -1 = distance only
  LOCUS_TABLE_T *locus_table,
  TRANSCRIPT_TABLE_T *transcript_table,
  LINK_TABLE_T *link_table
) {
  int i;
  char *key;
  int n_loci = locus_table->n_loci;
  LOCUS_T *loci = locus_table->loci;
  int n_transcripts = transcript_table->n_transcripts;
  TRANSCRIPT_T *transcripts = transcript_table->transcripts;
  char *histone_name = (i_hist == -1) ? "" : get_nth_string(i_hist, options->histone_names);
  int max_distance = options->max_distances[(i_hist == -1) ? 0 : i_hist];
  int first_link_index = link_table->n_links;	// save index of first added link
  TRANSCRIPT_T *transcript;
  LOCUS_T *locus;

  if (i_hist == -1) {
    DEBUG_FMT(NORMAL_VERBOSE, "  Distance: %d\n", max_distance);
  } else {
    DEBUG_FMT(NORMAL_VERBOSE, "  Histone: %s Distance: %d\n", histone_name, max_distance);
  }

  //
  // Find the closest loci for each TSS.  Save distances in a hash table.
  // Add links to table if they satisfy the distance criteria.
  // Add links to all closest loci unless --no-closest-loci given.
  //
  DEBUG_MSG(NORMAL_VERBOSE, "  Finding the closest loci to each TSS and creating links.\n");
  HASH_TABLE ht_closest_locus = hash_create(MAX(n_transcripts,1000), NULL);
  long *distance_to_locus = mm_malloc(n_transcripts * sizeof(long));
  long distance = 0;
  int save_locus_i = 0;
  int first_valid_dist_i = n_loci;
  int closest_locus_i = n_loci;
  long closest_locus_distance;
  int trans_i, locus_i;
  bool closest_locus_found;

  // Loop over transcripts.
  for (trans_i=locus_i=0; trans_i<n_transcripts; trans_i++) {
    TRANSCRIPT_T *transcript = &transcripts[trans_i];	// current transcript;

    // Loci for this TSS must come at or after the MIN(closest_locus_i, first_valid_dist_i).
    if (first_valid_dist_i == n_loci && closest_locus_i == n_loci) {
      locus_i = save_locus_i;
    } else {
      locus_i = MIN(closest_locus_i, first_valid_dist_i);
    }

    // Initialize the closest locus variables.
    save_locus_i = locus_i;	// current locus in case no closest/valid locus found
    first_valid_dist_i = n_loci;// records first locus that meets the distance constraint
    closest_locus_i = n_loci;	// index of closest locus (so far) to this TSS
    closest_locus_distance = -(MAX_CHRM_LENGTH+1);
    closest_locus_found = false;

    // Loop over remaining loci to record the closest locus index and valid links.
    for ( ; locus_i < n_loci; locus_i++) {
      // Get distance.
      locus = &loci[locus_i];
      distance = locus_to_trans_distance(locus, transcript, MAX_CHRM_LENGTH);

      // Skip if locus on prior chromosome.
      if (distance < -MAX_CHRM_LENGTH) continue;

      // Stop if locus on subsequent chromosome.
      if (distance > MAX_CHRM_LENGTH) break;

      // Stop if distance exceeds constraint and closest_locus_distance.
      if (distance > max_distance && distance > labs(closest_locus_distance)) break;

      // Record if this is the first locus meeting the distance criterion.
      if (first_valid_dist_i == n_loci && labs(distance) <= max_distance) first_valid_dist_i = locus_i;

      // Record if this is the closest locus so far for this TSS.
      if (labs(distance) < labs(closest_locus_distance)) {
        closest_locus_distance = distance;
        closest_locus_i = locus_i;
      }

      // Create link if it meets distance constraint.
      if (labs(distance) <= max_distance) {
        create_link(link_table, transcript, locus, i_hist, distance);
      }
    } // locus

    // Done with all TSSes if remaining loci are on prior chromosomes.
    if (closest_locus_i == n_loci && distance < MAX_CHRM_LENGTH) closest_locus_i = n_loci+1;

    // Was there a locus on the TSS chromosome?
    if (closest_locus_i < n_loci) {

      // Create an entry in the closest locus hash table; key is <transcript_id>.
      key = get_transcript_id(transcript);
      distance_to_locus[trans_i] = closest_locus_distance;
      hash_insert_str_value(key, &distance_to_locus[trans_i], ht_closest_locus);

      // Find and create links to all closest loci.
      int tied_i;
      for (tied_i=closest_locus_i; tied_i<n_loci; tied_i++) {
	locus = &loci[tied_i];
        distance = locus_to_trans_distance(locus, transcript, MAX_CHRM_LENGTH);
        // Stop if gone past tied loci.
        if (labs(distance) > labs(closest_locus_distance)) break;
        if (labs(distance) <= max_distance) {
          // Valid link; created above already.
        } else {
          // Add link unless --no-closest-locus was given.
          if (!options->no_closest_locus) {
	    create_link(link_table, transcript, locus, i_hist, closest_locus_distance);
          }
        }
      }
    } // closest_locus_i < n_loci

  } // transcript

  //
  // Create a hash table of the new links.
  // Key is <transcript_id>_<locus>.
  // This will allow us to check if a closest-tss link already exists.
  //
  int n_links = link_table->n_links;
  LINK_T **links = link_table->links;
  HASH_TABLE ht_links = hash_create(MAX(n_links,1000), NULL);
  for (i=first_link_index; i<n_links; i++) {
    LINK_T *link = links[i];
    int dummy=asprintf(&key, "%s_%s:%ld-%ld", link_tss_id(link), link_re_chrom(link), link_re_start(link), link_re_end(link));
    HASH_TABLE_ENTRY *hte = hash_lookup_str(key, ht_links);
    if (hte) {			// key already in hash table
      LINK_T *old_link = (LINK_T *) hash_get_entry_value(hte);
      fprintf(stderr, "Error: Duplicate link found: %s\n", key);
      fprintf(stderr, "\nNo output was generated.\n");
      exit(1);
    } else {
      hash_insert_str_value(key, link, ht_links);	// value for this key is the link
    }
    free(key);
  }

  //
  // Find the closest TSS(es) for each locus.  Save distances in a hash table.
  //
  DEBUG_MSG(NORMAL_VERBOSE, "  Finding the closest TSS(es) to each locus.\n");
  HASH_TABLE ht_closest_tss = n_loci > 0 ? hash_create(n_loci, NULL) : NULL;
  long *distance_to_tss = mm_malloc(n_loci * sizeof(long));
  long closest_trans_distance;
  int closest_trans_i = -1;
  for (locus_i=trans_i=0; locus_i<n_loci; locus_i++) {
    LOCUS_T *locus = &loci[locus_i];		// current locus
    locus->num_closest_tss_links = 0;		// number of closest TSS links to this locus

    // Back up TSS to previous best.
    if (closest_trans_i != -1) trans_i = closest_trans_i;

    // Initialize new closest TSS.
    closest_trans_distance = -(MAX_CHRM_LENGTH+1);
    closest_trans_i = -1;

    // Advance transcript index until on next chromosome
    // or until the absolute value of the distance increases.
    for ( ; trans_i < n_transcripts; trans_i++) {
      transcript = &transcripts[trans_i];
      distance = -locus_to_trans_distance(locus, transcript, MAX_CHRM_LENGTH);
      if (distance < -MAX_CHRM_LENGTH) continue;	// on prior chromosome
      if (distance > MAX_CHRM_LENGTH) break;		// on next chromosome
      if (labs(distance) > labs(closest_trans_distance)) break;	// distance increased 
      // Record if this is the closest TSS so far.
      if (labs(distance) < labs(closest_trans_distance)) {
        closest_trans_distance = distance;
        closest_trans_i = trans_i;
      }
    }

    // Done with all loci if remaining TSSs are on prior chromosomes.
    if (closest_trans_i == -1 && distance < MAX_CHRM_LENGTH) closest_trans_i = n_transcripts+1;

    // Create links for the closest TSSes if they aren't already there and
    // --no-closest-tss was not given.  Create distance entries in the closest TSS hash table.
    if (closest_trans_i != -1) {
      int tied_i;
      for (tied_i=closest_trans_i; tied_i<n_transcripts; tied_i++) {
	transcript = &transcripts[tied_i];
        distance = -locus_to_trans_distance(locus, transcript, MAX_CHRM_LENGTH);
        if (labs(distance) > labs(closest_trans_distance)) break;
        if (!options->no_closest_tss) {			// Include closest TSSes.
	  // Key is <transcript_id>_<locus>.
	  int dummy=asprintf(&key, "%s_%s:%ld-%ld", get_transcript_id(transcript), get_locus_chrom(locus), get_locus_start(locus), get_locus_end(locus));
	  HASH_TABLE_ENTRY *hte = hash_lookup_str(key, ht_links);
	  if (! hte) {					// Add new link.
	    LINK_T *link = create_link(link_table, transcript, locus, i_hist, closest_trans_distance);
	    hash_insert_str_value(key, link, ht_links);	// value for this key is the link
	  }
	  free(key);
        }
      }

      // Create an entry in the closest TSS hash table; key is <locus>.
      int dummy=asprintf(&key, "%s:%ld-%ld", get_locus_chrom(locus), get_locus_start(locus), get_locus_end(locus));
      distance_to_tss[locus_i] = closest_trans_distance;
      hash_insert_str_value(key, &distance_to_tss[locus_i], ht_closest_tss);
      free(key);
    }
  } // locus

  //
  // Mark all the links with closest_locus and closest_tss.
  //
  n_links = link_table->n_links;
  links = link_table->links;			// table pointer may have changed due to create_link above
  int link_i;
  HASH_TABLE_ENTRY *hte;
  for (link_i=0; link_i<n_links; link_i++) {
    LINK_T *link = links[link_i];
    // See if this is closest locus for TSS.
    key = link_tss_id(link);
    hte = hash_lookup_str(key, ht_closest_locus);
    if (hte) {
      long *shortest_distance_ptr = (long *) hash_get_entry_value(hte);
      // This link is to a closest Locus; label it.
      if (labs(link->distance) == labs(*shortest_distance_ptr)) {
        link->closest_locus = 'T';
      }
    } else {
      fprintf(stderr, "Error: Couldn't find TSS %s in the closest locus hash table.\n", key);
      exit(1);
    }
    // Get closest distance for this link's locus.
    int dummy=asprintf(&key, "%s:%ld-%ld", link_re_chrom(link), link_re_start(link), link_re_end(link));
    hte = hash_lookup_str(key, ht_closest_tss);
    if (hte) {
      long *shortest_distance_ptr = (long *) hash_get_entry_value(hte);
      // This link is to a closest TSS; label it.
      if (labs(link->distance) == labs(*shortest_distance_ptr)) {
        link->closest_tss = 'T';
        link->locus->num_closest_tss_links++;			// Number of links to closest TSSes to this locus
      }
    } else {
      fprintf(stderr, "Error: Couldn't find locus %s in the closest TSS hash table.\n", key);
      exit(1);
    }
    free(key);
  } // link

  DEBUG_FMT(NORMAL_VERBOSE, "  Created %d links.\n", n_links);

  // Free space.
  hash_destroy(ht_links);
  hash_destroy(ht_closest_tss);
  hash_destroy(ht_closest_locus);
  free(distance_to_tss);

  return(first_link_index);
} // update_link_table

/*************************************************************************
 * Get histone levels at loci used in links for the current histone modification.
 *************************************************************************/
static void get_histone_levels(
  TGENE_OPTIONS_T *options,
  int i_hist,			// current histone index
  LOCUS_TABLE_T *locus_table,
  LINK_TABLE_T *link_table,
  int first_link_index		// update from this link to last link
) {
  int i, j;
  char *histone_name = get_nth_string(i_hist, options->histone_names);
  int n_loci = locus_table->n_loci;
  int n_tissues = locus_table->n_tissues;
  LOCUS_T *loci = locus_table->loci;
  int n_links = link_table->n_links;
  LINK_T **links = link_table->links;

  // Mark which loci are used in links for this histone.
  for (i=0; i<n_loci; i++) {
    loci[i].used_in_link = 0;
  }
  for (i=first_link_index; i<n_links; i++) links[i]->locus->used_in_link++;

  // Set the histone levels to 0 for all loci used in links.
  for (i=0; i<n_loci; i++) {
    if (loci[i].used_in_link > 0) for (j=0; j<n_tissues; j++) loci[i].histone_level[j] = 0;
    loci[i].max_hist = 0;	// will hold maximum level for this histone
  }

  // Read in the BED files for different tissues for the current histone.
  // Update the histone levels.
  int tissue_i;
  DEBUG_FMT(NORMAL_VERBOSE, "  Creating the histone table for histone %s with %d loci.\n", 
    get_nth_string(i_hist, options->histone_names), n_loci); 
  DEBUG_FMT(NORMAL_VERBOSE, "    Reading in the %s histone files for %d tissues.\n", 
    get_nth_string(i_hist, options->histone_names), n_tissues); 
  for (tissue_i=0; tissue_i<n_tissues; tissue_i++) {		// tissue
    char *tissue = get_nth_string(tissue_i, options->tissue_names);
    // get listing of histone directory for current tissue
    STRING_LIST_T *file_names = get_file_names(options->histone_root, tissue, HISTONE_FILE_NAME_FMT, histone_name);
    if (file_names == NULL) {
      DEBUG_FMT(NORMAL_VERBOSE, "      Tissue: %s Histone file: %s Entries %d\n", tissue, "none", 0);
      continue;
    }
    int n_file_names = get_num_strings(file_names);
    for (j=0; j<n_file_names; j++) {				// histone file
      char *file_name = get_nth_string(j, file_names);
      // Read the histone file into an array of BED entries.
      FILE *histone_file;
      if ((histone_file = fopen(file_name, "r")) == NULL) {
	die("Unable to open histone file \"%s\" for reading.\n", file_name);
      }
      BED_ENTRY_T **bed = NULL;
      int n_bed = 0;
      BED_TYPE bed_type = BROADPEAK;
      if (strcmp("narrowPeak", file_name+strlen(file_name)-10) == 0) {
        bed_type = NARROWPEAK;
      }
      DEBUG_FMT(NORMAL_VERBOSE, "      Tissue: %s Histone file: %s ", tissue, file_name);
      // Convert coordinates to 1-based, closed.
      n_bed = read_bed_file(histone_file, true, bed_type, &bed);
      if (n_bed <= 0) {
        fprintf(stderr, "\nNo valid entries found in histone BED file '%s'.\n", file_name);
        fprintf(stderr, "\nNo output was generated.\n");
        exit(1);
      }
      fclose(histone_file);
      DEBUG_FMT(NORMAL_VERBOSE, "Entries: %d\n", n_bed);

      // Sort the histone bed entries by genomic start.
      if (n_bed > 0) qsort(bed, n_bed, sizeof(BED_ENTRY_T *), compare_bed_coords);

      // Update the used loci with this histone from the BED entries.
      // METHOD: 
      // Note: This assumes that loci are sorted by genomic start.
      // Note: This assumes that bed entries are sorted by genomic start.
      //   1) advance locus index to first locus that overlaps the current histone bed entry
      //   2) linear search from locus index forward until locus too far downstream is reached
      int locus_i = 0, bed_i = 0;
      for (bed_i=0; bed_i<n_bed; bed_i++) {
        BED_ENTRY_T *bed_entry = bed[bed_i];
	// Advance locus index until we find the first locus that
	// overlaps or is downstream of the bed entry.
	int overlap = 0;
	for ( ; locus_i < n_loci; locus_i++) {
          if (loci[locus_i].used_in_link > 0) {
	    overlap = locus_to_bed_overlap(&loci[locus_i], bed_entry);
	    if (overlap >= 0) break;
          }
	}
	if (locus_i == n_loci) locus_i--;

	// Do a linear search for loci stopping at the first locus 
        // that is downstream of the bed entry.
	int current_locus_i = locus_i;
	while (overlap <= 0) {
	  if (overlap == 0) {			// Found an overlap.
            loci[current_locus_i].histone_level[tissue_i] += get_bed_signal(bed_entry);
	  }
	  if (++current_locus_i == n_loci) {
	    break;
	  } else { 				// update overlap
            if (loci[current_locus_i].used_in_link > 0) {
	      overlap = locus_to_bed_overlap(&loci[current_locus_i], bed_entry);
            } else {
              overlap = -1;
            }
	  }
	} // overlap
        free_bed_entry(bed_entry);
      } // bed_entry
    } // histone file
    free_string_list(file_names);
    // Update the maximum histone level for each locus.
    if (! options->valid_tissues[i_hist][tissue_i]) continue;
    int locus_i;
    for (locus_i=0; locus_i<n_loci; locus_i++) {
      LOCUS_T *locus = &loci[locus_i];
      if (locus->used_in_link > 0) {
        if (locus->histone_level[tissue_i] > locus->max_hist) locus->max_hist = locus->histone_level[tissue_i];
      }
    }
  } // tissue
} // get_histone_levels

/*************************************************************************
 * Create table of transcripts sorted by genomic coordinates.
 *************************************************************************/
static TRANSCRIPT_TABLE_T *create_transcript_table(
  TGENE_OPTIONS_T *options,
  LOCUS_TABLE_T *locus_table
) {
  int i, j, tissue_i;
  int n_tissues = options->n_tissues;
  char *file_name = options->annotation_file;
  bool have_tss_info = false;

  // Read the annotation file.
  FILE *annotation_file;
  if ((annotation_file = fopen(file_name, "r")) == NULL) {
    die("Unable to open annotation file \"%s\" for reading.\n", options->annotation_file);
  }
  GTF_ENTRY_T **annotations = NULL;		// array of gtf entry pointers
  DEBUG_FMT(NORMAL_VERBOSE, "Reading the annotation file '%s'.\n", file_name); 
  // GENCODE GTF is in closed format.
  // RefSeq GTF is in closed format.
  int n_entries = read_gtf_file(annotation_file, false, GENCODE, &annotations);
  if (n_entries <= 0) {
    fprintf(stderr, "\nNo valid entries found in GTF annotation file '%s'.\n", options->annotation_file);
    fprintf(stderr, "\nNo output was generated.\n");
    exit(1);
  }
  fclose(annotation_file);
  DEBUG_FMT(NORMAL_VERBOSE, "The annotation file contained %d entries.\n", n_entries); 

  // Create a hash table of allowed transcript types to speed up checking transcript type.
  STRING_LIST_T *allowed_types_list = new_string_list_char_split(',', options->transcript_types);
  int num_allowed_types = get_num_strings(allowed_types_list);
  HASH_TABLE ht_allowed_types = num_allowed_types > 0 ? hash_create(num_allowed_types, NULL) : NULL;
  for (i=0; i<num_allowed_types; i++) {
    hash_insert_str(get_nth_string(i, allowed_types_list), ht_allowed_types);
  }
  free_string_list(allowed_types_list);
  
  // Create a hash table of the allowed chromosome names.
  HASH_TABLE ht_chrom_names = hash_create(1000, NULL);
  for (i=0; i<n_entries; i++) {
    GTF_ENTRY_T *gtf = annotations[i];
    HASH_TABLE_ENTRY *hte = hash_lookup_str(gtf->seqname, ht_chrom_names);
    if (! hte) hash_insert_str(gtf->seqname, ht_chrom_names);
  }

  // Check that the loci don't use any unknown chromosome names.
  int locus_i;
  LOCUS_T *loci = locus_table->loci;
  int n_loci = locus_table->n_loci;
  for (locus_i=0; locus_i<n_loci; locus_i++) {
    LOCUS_T *locus = &loci[locus_i];		// current locus
    HASH_TABLE_ENTRY *hte = hash_lookup_str(get_locus_chrom(locus), ht_chrom_names);
    if (! hte) {
      fprintf(stderr, "\nERROR: The locus file contains a chromosome name that is not in the annotation file: '%s'\n", get_locus_chrom(locus));
      fprintf(stderr, "\nThe chromosome names in the annotation file are:\n\n");
      STRING_LIST_T *keys = hash_get_keys(ht_chrom_names);
      sort_string_list(keys);
      write_string_list(", ", keys, stderr);
      hash_destroy(ht_chrom_names);
      free_string_list(keys);
      HASH_TABLE ht_chrom_names = hash_create(1000, NULL);
      for (locus_i=0; locus_i<n_loci; locus_i++) {
        locus = &loci[locus_i];		// current locus
	HASH_TABLE_ENTRY *hte = hash_lookup_str(get_locus_chrom(locus), ht_chrom_names);
	if (! hte) hash_insert_str(get_locus_chrom(locus), ht_chrom_names);
      }
      keys = hash_get_keys(ht_chrom_names);
      sort_string_list(keys);
      fprintf(stderr, "\nThe chromosome names in your locus file are:\n\n");
      write_string_list(", ", keys, stderr);
      fprintf(stderr, "\nMake sure that your locus file is for the same genome and uses\n"
        "the same chromosome nomenclature (e.g., UCSC or Ensembl) as the annotation file you are using.\n");
      hash_destroy(ht_chrom_names);
      free_string_list(keys);
      exit(1);
    }
  }

  // Create an array of transcript objects
  // and create a hash table to index it by transcript_id.
  // (Each entry in the hash table is a linked list of transcript ID pointers.)
  // Mark transcripts that are not for allowed transcript type(s).
  TRANSCRIPT_T *transcripts = mm_malloc(sizeof(TRANSCRIPT_T) * n_entries);
  HASH_TABLE ht_transcripts = n_entries > 0 ? hash_create(n_entries, NULL) : NULL;
  DEBUG_FMT(NORMAL_VERBOSE, "Creating the transcript table with %d entries.\n", n_entries); 
  bool non_unique_ids = false;
  for (i=0; i<n_entries; i++) {
    // Convert transcript to just the TSS.  
    if (annotations[i]->strand == '+') {
      annotations[i]->end = annotations[i]->start;
    } else {
      annotations[i]->start = annotations[i]->end;
    }
    // Create empty gene_name if necessary.
    if (annotations[i]->gene_name == NULL) annotations[i]->gene_name = mm_calloc(1, sizeof(char));

    // Find out if the GTF file contains different values for TSS and Gene ID.
    // This is useful to determine if the TSS column should be displayed.
    if (!have_tss_info) {
      char *tss_id = annotations[i]->transcript_id;
      char *gene_id = annotations[i]->gene_id;
      if (strcmp(tss_id, gene_id)) have_tss_info = true;
    }
    transcripts[i].annotation = annotations[i];	// pointer to the annotation entry for this transcript 
    transcripts[i].valid_type = hash_lookup_str(annotations[i]->transcript_type, ht_allowed_types);
    // will hold the expression per tissue
    transcripts[i].expr = n_tissues==0 ? NULL : (double *) mm_calloc(n_tissues, sizeof(double));
    transcripts[i].max_expr = 0;		// records the max expression of transcript across tissues
    
    // Make sure the transcript IDs are unique (unless -use-gene-ids given).
    // Store a list of transcript IDs associated with the gene_id if -use-gene-ids given.
    // Otherwise the list should just have one entry.
    char *id = (options->use_gene_ids) ? annotations[i]->gene_id : annotations[i]->transcript_id;
    HASH_TABLE_ENTRY *hte = hash_lookup_str(id, ht_transcripts);
    if (hte) {
      // non-unique allowed if using gene_ids
      if (options->use_gene_ids) {	
        // Add this transcript_id to the list for the gene_id.
        LINKLST_T *tr_list = (LINKLST_T *) hash_get_entry_value(hte);
        linklst_push(&transcripts[i], tr_list);
      } else {
        non_unique_ids = true;
	fprintf(stderr, "Error: Annotation file entry %d contains the non-unique transcript_id \"%s\".\n", i, id);
      }
    } else {
      // Create a list of transcript IDs (just one ID unless --use-gene-ids given).
      LINKLST_T *tr_list = linklst_create();
      linklst_push(&transcripts[i], tr_list);
      hash_insert_str_value(id, tr_list, ht_transcripts);
    }
  } // entry
  if (non_unique_ids) exit(1);

  // Read in the expression files for each tissue.
  if (n_tissues > 0) DEBUG_FMT(NORMAL_VERBOSE, "  Reading in the expression files for %d tissues.\n", n_tissues); 
  for (tissue_i=0; tissue_i<n_tissues; tissue_i++) {
    char *tissue = get_nth_string(tissue_i, options->tissue_names);

    // get expression file name for current tissue
    STRING_LIST_T *file_names = 
      get_file_names(options->expression_root, tissue, EXPRESSION_FILE_NAME_FMT, options->rna_source);
    if (file_names == NULL) {
      DEBUG_FMT(NORMAL_VERBOSE, "    Tissue: %s Expression file: %s Entries %d\n", tissue, "none", 0);
      continue;
    }
    // Only one expression file should be found in this tissue's directory.
    if (get_num_strings(file_names) != 1) {
      die("Found %d expression files matching '" EXPRESSION_FILE_NAME_FMT "' for tissue %s\n"
        "but there should only be 1 such file.\n", 
        get_num_strings(file_names), options->rna_source, tissue);
    }
    char *file_name = get_nth_string(0, file_names);

    // Read the expression file into an array of GTF entries.
    FILE *expression_file;
    if ((expression_file = fopen(file_name, "r")) == NULL) {
      die("Unable to open expression file \"%s\" for reading.\n", file_name);
    }
    GTF_ENTRY_T **expression = NULL;
    DEBUG_FMT(NORMAL_VERBOSE, "    Tissue: %s Expression file: %s ", tissue, file_name);
    // CAGE GTF is in closed format.  LongPap GTF is in closed format.
    int n_expr = read_gtf_file(expression_file, false, options->rna_source_type, &expression);
    if (n_expr <= 0) {
      fprintf(stderr, "\nNo valid entries found in GTF expression file '%s'.\n", file_name);
      fprintf(stderr, "\nNo output was generated.\n");
      exit(1);
    }
    fclose(expression_file);
    DEBUG_FMT(NORMAL_VERBOSE, "Entries: %d\n", n_expr);

    // Add the expression data for this tissue to each valid transcript.
    for (j=0; j<n_expr; j++) {

      GTF_ENTRY_T *entry = expression[j];
      char *id = entry->transcript_id;
      HASH_TABLE_ENTRY *hte = hash_lookup_str(id, ht_transcripts);
      if (! hte) {
	fprintf(stderr, "Error: Unknown transcript_id '%s' at line %d in expression file %s\n", id, j+1, file_name);
      } else {
        // Get the list of transcript_ids associated with this expression's ID (which
	// could be the same as the gene_id, in which case --use-gene-ids must be given).
        LINKLST_T *tr_list = (LINKLST_T *) hash_get_entry_value(hte);
        LL_LINK_T *ll_link = linklst_first(tr_list);
        // Update the expression for each transcript in the ID's list.
        for ( ; ll_link; ll_link = linklst_next(ll_link)) { 
          TRANSCRIPT_T *transcript = (TRANSCRIPT_T *) linklst_get(ll_link);
	  if (transcript->valid_type) {
	    if (transcript->expr == NULL) transcript->expr = (double *) mm_calloc(n_tissues, sizeof(double));
	    transcript->expr[tissue_i] = entry->expr;
	    if (entry->expr > transcript->max_expr) transcript->max_expr = entry->expr;
	  } // valid type
        } // transcripts associated with this ID (could be transcript_id or gene_id)
      } // known transcript
      // Free the expression entry.
      free_gtf_entry(entry);
    } // j (expression file entry)
    free_string_list(file_names);
  } // tissue_i (tissue)

  // Remove transcripts with invalid transcript types from array.
  DEBUG_MSG(NORMAL_VERBOSE, "  Removing transcripts for invalid transcript types from the table.\n");
  int n_valid = 0;
  for (i=0; i<n_entries; i++) {
    if (transcripts[i].valid_type) {
      if (n_valid < i) transcripts[n_valid] = transcripts[i]; 
      n_valid++;
    }
  }
  DEBUG_FMT(NORMAL_VERBOSE, "  Transcript table now contains %d transcripts (%d removed).\n", n_valid, n_entries-n_valid); 

  // Sort the transcripts by genomic coordinate.  Break ties by putting longest transcript first.
  if (n_valid > 0) qsort(transcripts, n_valid, sizeof(TRANSCRIPT_T), compare_transcript_coords);

  // Create a new hash table from transcript_id to (sorted) transcript.
  hash_destroy(ht_transcripts);
  
  // free stuff
  hash_destroy(ht_allowed_types);
  hash_destroy(ht_chrom_names);

  // Fill in the transcript table (n_transcripts X n_tissues) object and return it.
  TRANSCRIPT_TABLE_T *transcript_table = mm_malloc(sizeof(TRANSCRIPT_TABLE_T));
  transcript_table->n_transcripts = n_valid;
  transcript_table->n_tissues = n_tissues;
  transcript_table->annotations = annotations;
  transcript_table->transcripts = transcripts;
  transcript_table->have_tss_info = have_tss_info;
  return(transcript_table);
} // create_transcript_table

/*************************************************************************
 * Add random noise to transcripts.
 *************************************************************************/
void add_transcript_noise(
  TGENE_OPTIONS_T *options,
  TRANSCRIPT_TABLE_T *transcript_table
) {
  int tissue_i, trans_i;
  TRANSCRIPT_T *transcripts = transcript_table->transcripts;
  int n_transcripts = transcript_table->n_transcripts;
  int n_tissues = transcript_table->n_tissues;

  DEBUG_FMT(NORMAL_VERBOSE, "  Adding noise to the transcript table with %d entries.\n", n_transcripts); 
  // Get statistics on each tissue.
  double *min_expr = mm_calloc(n_tissues, sizeof(double));
  for (trans_i=0; trans_i<n_transcripts; trans_i++) {
    TRANSCRIPT_T *transcript = &transcripts[trans_i];	// current transcript;
    double *expr = transcript->expr;
    // Find the smallest non-zero expression for this tissue.
    for (tissue_i=0; tissue_i<n_tissues; tissue_i++) {
      double x = expr[tissue_i];
      double min_x = min_expr[tissue_i];
      if (x > 0 && (x < min_x || min_x == 0)) min_expr[tissue_i] = x ;
    }
  }

  // Set the noise mean to mu = minimum expression / NOISE_FRACTION, sigma = mu / 4.
  DEBUG_MSG(HIGH_VERBOSE, "    Minimum (non-zero) expression values for tissues: ");
  double *noise_mu = mm_malloc(n_tissues * sizeof(double));
  double *noise_sigma = mm_malloc(n_tissues * sizeof(double));
  for (tissue_i=0; tissue_i<n_tissues; tissue_i++) {
    DEBUG_FMT(HIGH_VERBOSE, " %.3g", min_expr[tissue_i]);
    // Add noise to zeros based on value of minimum expression of transcript.
    noise_mu[tissue_i] = min_expr[tissue_i] / NOISE_FRACTION;
    noise_sigma[tissue_i] = noise_mu[tissue_i] / 4.0;
  }
  DEBUG_MSG(HIGH_VERBOSE, "\n");

  // Add noise to zero expression only.
  srand_mt(options->seed);	// Make sure the noise is not affected by other things.
  for (tissue_i=0; tissue_i<n_tissues; tissue_i++) {
    if (noise_mu[tissue_i] == 0) continue;		// no expression data for tissue
    srand_mt(options->seed);	// Make sure the noise for this tissue is always the same.
    for (trans_i=0; trans_i<n_transcripts; trans_i++) {
      TRANSCRIPT_T *transcript = &transcripts[trans_i];	// current transcript;
      double *expr = transcript->expr;
      if (expr[tissue_i] != 0) continue;		// only add noise to 0's
      double noise = rand_gaussian(noise_mu[tissue_i], noise_sigma[tissue_i]);
      expr[tissue_i] += noise;
    } // transcript
  } // tissue

} // add_transcript_noise

/*************************************************************************
 * Output TSV file of expression for each transcript.
 *************************************************************************/
static void output_transcript_tsv_file(
  TGENE_OPTIONS_T *options,
  char *name,
  TRANSCRIPT_TABLE_T *transcript_table
) {
  int i;

  int dir_len = strlen(options->output_dirname);
  STR_T *file_name = str_create(dir_len + 100);
  str_setf(file_name, "%s/%s.tsv", options->output_dirname, name);
  FILE *tt_file;
  if ((tt_file = fopen(str_internal(file_name), "w")) == NULL) {
    die("Unable to open transcript table file \"%s\" for writing.\n", file_name);
  }
  str_destroy(file_name, false);

  fprintf(tt_file, "TSS_ID");
  for (i=0; i<transcript_table->n_tissues; i++) fprintf(tt_file, "\t%s", get_nth_string(i, options->tissue_names));
  fprintf(tt_file, "\n");
  for (i=0; i<transcript_table->n_transcripts; i++) {
    fprintf(tt_file, "%s", transcript_table->transcripts[i].annotation->transcript_id);
    int j;
    for (j=0; j<transcript_table->n_tissues; j++) {
      fprintf(tt_file, "\t%0.4f", transcript_table->transcripts[i].expr[j]);
    }
    fprintf(tt_file, "\n");
  }
  fclose(tt_file);
} // output_transcript_tsv_file

/*************************************************************************
 * Convert expression to log(1+expr).
 *************************************************************************/
static void log_convert_transcript_table(
  TGENE_OPTIONS_T *options,
  TRANSCRIPT_TABLE_T *transcript_table		// the transcripts
) {
  int i;
  for (i=0; i<transcript_table->n_transcripts; i++) {
    double *expr = transcript_table->transcripts[i].expr; 
    int j;
    for (j=0; j<transcript_table->n_tissues; j++) {
      expr[j] = log(1+expr[j]);    
    }
  }
} // log_convert_transcript_table

/*************************************************************************
 * Add random noise to histone levels.
 *************************************************************************/
void add_histone_noise(
  TGENE_OPTIONS_T *options,
  LOCUS_TABLE_T *locus_table
) {
  int tissue_i, locus_i;
  LOCUS_T *loci = locus_table->loci;
  int n_loci = locus_table->n_loci;
  int n_tissues = locus_table->n_tissues;
  srand_mt(options->seed);	// Make sure the noise is not affected by other things.

  DEBUG_FMT(NORMAL_VERBOSE, "    Adding noise to the locus table with %d entries.\n", n_loci); 
  // Get statistics on each tissue.
  double *min_level = mm_calloc(n_tissues, sizeof(double));
  for (locus_i=0; locus_i<n_loci; locus_i++) {
    LOCUS_T *locus = &loci[locus_i];		// current locus
    double *level = locus->histone_level;
    // Find the smallest non-zero histone level for this tissue.
    for (tissue_i=0; tissue_i<n_tissues; tissue_i++) {
      double x = level[tissue_i];
      double min_x = min_level[tissue_i];
      if (x > 0 && (x < min_x || min_x == 0)) min_level[tissue_i] = x ;
    }
  }

  // Set the noise mean to mu = minimum level / NOISE_FRACTION, sigma = mu / 4.
  DEBUG_MSG(HIGH_VERBOSE, "      Minimum (non-zero) histone levels for tissues: ");
  double *noise_mu = mm_malloc(n_tissues * sizeof(double));
  double *noise_sigma = mm_malloc(n_tissues * sizeof(double));
  for (tissue_i=0; tissue_i<n_tissues; tissue_i++) {
    DEBUG_FMT(HIGH_VERBOSE, " %.3g", min_level[tissue_i]);
    noise_mu[tissue_i] = min_level[tissue_i] / NOISE_FRACTION;
    noise_sigma[tissue_i] = noise_mu[tissue_i] / 4.0;
  }
  DEBUG_MSG(HIGH_VERBOSE, "\n");

  // Add noise to zero level only.
  for (tissue_i=0; tissue_i<n_tissues; tissue_i++) {
    if (noise_mu[tissue_i] == 0) continue;	// no histone for current tissue
    srand_mt(options->seed);			// Make sure the noise for this tissue is always the same.
    for (locus_i=0; locus_i<n_loci; locus_i++) {
      LOCUS_T *locus = &loci[locus_i];		// current locus
      double *level = locus->histone_level;
      if (level[tissue_i] != 0) continue;	// add noise to only to 0's
      double noise = rand_gaussian(noise_mu[tissue_i], noise_sigma[tissue_i]);
      level[tissue_i] += noise;
    } // locus
  } // tissue

} // add_histone_noise

/*************************************************************************
 * Convert histone level to log(1+histone_level).
 *************************************************************************/
static void log_convert_locus_table(
  TGENE_OPTIONS_T *options,
  LOCUS_TABLE_T *locus_table
) {
  int i;
  for (i=0; i<locus_table->n_loci; i++) {
    if (locus_table->loci[i].histone_level == NULL) continue;
    double *histone_level = locus_table->loci[i].histone_level; 
    int j;
    for (j=0; j<locus_table->n_tissues; j++) {
      histone_level[j] = log(1+histone_level[j]);    
    }
  }
} // log_convert_locus_table

/*************************************************************************
 * Output TSV file of histone levels for each locus.
 *************************************************************************/
static void output_histone_tsv_file(
  TGENE_OPTIONS_T *options,
  char *name,
  int i_hist,				// index of histone
  LOCUS_TABLE_T *locus_table
) {
  int i;
  char *histone_name = get_nth_string(i_hist, options->histone_names);

  int dir_len = strlen(options->output_dirname);
  STR_T *file_name = str_create(dir_len + 100);
  str_setf(file_name, "%s/%s.%s.tsv", options->output_dirname, name, histone_name);
  FILE *lt_file;
  if ((lt_file = fopen(str_internal(file_name), "w")) == NULL) {
    die("Unable to open histone level file \"%s\" for writing.\n", file_name);
  }
  str_destroy(file_name, false);

  fprintf(lt_file, "RE_Locus");
  for (i=0; i<locus_table->n_tissues; i++) fprintf(lt_file, "\t%s", get_nth_string(i, options->tissue_names));
  fprintf(lt_file, "\n");
  for (i=0; i<locus_table->n_loci; i++) {
    LOCUS_T *locus = &locus_table->loci[i];
    if (locus->used_in_link == 0) continue;
    fprintf(lt_file, "%s:%ld-%ld", locus->reg_el->chrom, locus->reg_el->chrom_start, locus->reg_el->chrom_end);
    int j;
    for (j=0; j<locus_table->n_tissues; j++) {
      fprintf(lt_file, "\t%0.4f", locus->histone_level[j]);
    }
    fprintf(lt_file, "\n");
  }
  fclose(lt_file);
} // output_histone_tsv_file

/*************************************************************************
 * Add links for the current histone modification to the link table.
 * Compute null scores using permuted expression columns.
 *************************************************************************/
static void get_correlations(
  TGENE_OPTIONS_T *options,
  int i_hist,			// index of current histone
  int first_link_index,		// current histone's links start here
  LINK_TABLE_T *link_table
) {
  int i, j, k, link_i, array_i;
  int n_links = link_table->n_links;
  LINK_T **links = link_table->links;
  int n_tissues = get_num_strings(options->tissue_names);
  int n_perms = link_table->n_perms;
  double lecat = options->lecat;	// low expression correlation adjustment threshold

  DEBUG_FMT(NORMAL_VERBOSE,
    "Getting the correlations for %d links for histone %s.\n", n_links, get_nth_string(i_hist, options->histone_names));

  // Get the indices of the valid tissues for this histone.
  int *valid_indices = (int *)mm_malloc(sizeof(int) * n_tissues);
  bool *is_valid_index = (bool *)mm_malloc(sizeof(bool) * n_tissues);
  int n_valid_tissues = 0;
  for (j=0; j<n_tissues; j++) {
    if (options->valid_tissues[i_hist][j]) {
      is_valid_index[j] = true;
      valid_indices[n_valid_tissues++] = j;
    } else {
      is_valid_index[j] = false;
    }
  }

  // Create a random ordering of the tissue indices to use for computing null distribution.
  double *perm_expr = (double *)mm_malloc(sizeof(double) * n_tissues);
  int **perm = (int **) mm_malloc(sizeof(int *) * n_perms);
  int *perm_indices = (int *)mm_malloc(sizeof(int) * n_valid_tissues);
  srand_mt(options->seed);	// Make sure the permutation is not affected by other things.
  for (i=0; i<n_perms; i++) {
    perm[i] = (int *)mm_malloc(sizeof(int) * n_tissues);
    //for (j=0; j<n_tissues; j++) (perm[i])[j] = j;
    //SHUFFLE(perm[i], n_tissues);
    for (j=0; j<n_valid_tissues; j++) perm_indices[j] = valid_indices[j];
    SHUFFLE(perm_indices, n_valid_tissues);
    for (j=k=0; j<n_tissues; j++) {
      (perm[i])[j] = is_valid_index[j] ? perm_indices[k++] : j;
    }
  }

  // Compute correlation of each link.
  array_i = 0;
  for (link_i=first_link_index; link_i<n_links; link_i++) {
    LINK_T *link = links[link_i];
    TRANSCRIPT_T *transcript = link->transcript;
    LOCUS_T *locus = link->locus;
    // Scale correlation if low expression correlation adjustment threshold
    // is set and the transcript has low expression.
    double scale_factor = 1;
    if (lecat > 0 && link_max_expr(link) < lecat) {
      scale_factor = link_max_expr(link) / lecat;
    }
    double rho = pearson_correlation(n_tissues, transcript->expr, locus->histone_level,
      NULL, NULL, NULL, NULL, true, options->valid_tissues[i_hist]);
      
    // Round results so order of tissues doesn't matter.
    RND(rho*scale_factor, 8, link->correlation);

    // Compute correlation with permuted order of expression tissues.
    for (i=0; i<n_perms; i++) {
      for (j=0; j<n_tissues; j++) perm_expr[j] = transcript->expr[(perm[i])[j]];
      rho = pearson_correlation(n_tissues, perm_expr, locus->histone_level,
        NULL, NULL, NULL, NULL, true, options->valid_tissues[i_hist]);
      RND(rho*scale_factor, 8, link->perm_correlation[i]);
    }
  } // link 

  // Free space.
  for (i=0; i<n_perms; i++) myfree(perm[i]);
  myfree(perm);
  myfree(perm_expr);
  myfree(perm_indices);
  myfree(valid_indices);
  myfree(is_valid_index);

} // get_correlations

/*************************************************************************
 * Entry point for tgene
 *************************************************************************/
int main(int argc, char *argv[]) {
  TGENE_OPTIONS_T options;
  HTMLWR_T *html;
  JSONWR_T *json;

  // command line processing
  char *commandline = process_command_line(argc, argv, &options);

  // Set random number generators.
  srand_mt(options.seed);

  // Create output directory
  if (create_output_directory(options.output_dirname, options.allow_clobber,
    (verbosity >= NORMAL_VERBOSE))) {
    die("Couldn't create output directory %s.\n", options.output_dirname);
  }

  // Read in the locus file and initialize the locus table (without histone levels).
  LOCUS_TABLE_T *locus_table = initialize_locus_table(&options);

  // Read in the annotation file and expression files to create the transcript table.
  TRANSCRIPT_TABLE_T *transcript_table = create_transcript_table(&options, locus_table);
  if (options.n_tissues > 0) {	// Tissue panel provided.
    // Print the transcript vs. tissue expression TSV file.
    output_transcript_tsv_file(&options, "TrExp", transcript_table);
    // Add noise to transcripts.
    if (! options.no_noise) add_transcript_noise(&options, transcript_table);
    // Print the transcript vs. tissue expression+noise TSV file.
    output_transcript_tsv_file(&options, "TrExp+noise", transcript_table);
    // Convert the expression levels to log(1+expr).
    log_convert_transcript_table(&options, transcript_table);
  }

  // Initialize the link table.
  LINK_TABLE_T link_table;
  link_table.n_links = 0;
  link_table.links = NULL;
  link_table.n_perms = N_PERMS;
  link_table.have_tss_info = transcript_table->have_tss_info;

  // Add links to link table for each histone; (needed because max distance can vary)
  int n_histones = options.n_histones;
  if (n_histones == 0) {
    DEBUG_MSG(NORMAL_VERBOSE, "Creating the link table based on distance alone.\n");
    // Initialize links.
    int first_link_index = update_link_table(&options, -1, locus_table, transcript_table, &link_table);
  } else {
    DEBUG_FMT(NORMAL_VERBOSE, "Creating the link table for %d histone type(s).\n", n_histones);
    int i_hist;
    for (i_hist=0; i_hist<n_histones; i_hist++) {
      // Initialize links for this histone.
      int first_link_index = update_link_table(&options, i_hist, locus_table, transcript_table, &link_table);
      // Compute histone levels for loci used in this histone's links.
      get_histone_levels(&options, i_hist, locus_table, &link_table, first_link_index);
      // Print histone TSV file.
      output_histone_tsv_file(&options, "HistLev", i_hist, locus_table);
      // Add noise to histone levels.
      if (! options.no_noise) add_histone_noise(&options, locus_table);
      // Print histone+noice TSV file.
      output_histone_tsv_file(&options, "HistLev+noise", i_hist, locus_table);
      // Convert the histone levels to log(1+histone_level).
      log_convert_locus_table(&options, locus_table);
      // Compute correlations and scores for this histone's links.
      get_correlations(&options, i_hist, first_link_index, &link_table);
    } // histone
  }

  // Finish the link table.
  finish_link_table(&options, &link_table);

  // Output the TSV and HTML files.
  output_links_tsv_file(&options, commandline, &link_table);		// sorts on gene score
  output_html_file(&options, argc, argv, &link_table);

  // Free stuff
  myfree(commandline);
  if (options.n_tissues > 0) free_string_list(options.tissue_names);
  if (options.n_histones > 0) free_string_list(options.histone_names);

  DEBUG_MSG(NORMAL_VERBOSE, "Program ends correctly.\n");
  return EXIT_SUCCESS;
} // main
