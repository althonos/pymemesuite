//
// fasta-holdout-set
//
#define DEFINE_GLOBALS
#include "streme-utils.h"

static char *DEFAULT_OUTPUT_DIRNAME = "fasta-holdout-set_out";

/*****************************************************************************
 * Print a usage message with an optional reason for failure.
 *****************************************************************************/
static void usage(char *format, ...) {
  va_list argp;
  char *usage =
    "\n"
    "Usage: fasta-holdout-set [options] --p <primary sequences>\n"
    "\n"
    "   Options:\n"
    "     --p <filename>           primary (positive) sequence file name (required)\n"
    "     --n <filename>           control (negative) sequence file name;\n"
    "     --order <m>              use an m-order shuffle if creating\n"
    "                              control sequences from primary sequences;\n"
    "                              default: %d (DNA), %d (RNA), %d (Protein), %d (custom)\n"
    "     --o <output_dir>         output directory; default: '%s'\n"
    "     --oc <output_dir>        allow overwriting; default: '%s'\n"
    "     --hofract <hofract>      fraction of sequences in hold-out set;\n"
    "                              default: %g\n"
    "     --seed <seed>            random seed for shuffling sequences;\n"
    "                              default: %d\n"
    "     --dna                    sequences use standard DNA alphabet (default)\n"
    "     --rna                    sequences use standard RNA alphabet\n"
    "     --protein                sequences use standard protein alphabet\n"
    "     --alph <alph_file>       sequences use alphabet defined in <alph_file>;\n"
    "                              converts to uppercase unless both cases in core\n"
    "     --ccut <size>            centrally trim primary sequences to this maximum size;\n"
    "                              default: 0 (no trimming)\n"
    "     --help                   print this message and exit\n"
    "     --verbosity 1|2|3|4|5    level of diagnostic output (default: %d)\n"
    "                              1: none 2: helpful 3: debug 4: tons 5: ludicrous\n"
    "\n"
    ;
  if (format) {
    fprintf(stderr, "\n");
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, usage, 
    DEFAULT_DNA_ORDER, DEFAULT_RNA_ORDER, DEFAULT_PROT_ORDER, DEFAULT_CUSTOM_ORDER,
    DEFAULT_OUTPUT_DIRNAME, DEFAULT_OUTPUT_DIRNAME, 
    DEFAULT_HOFRACT, DEFAULT_SEED, DEFAULT_VERBOSITY
  );
  fflush(stderr);
  exit(EXIT_FAILURE);
} // usage

// An easy way of assigning a number to each option.
// Be careful that there are not more than 62 options as otherwise the value
// of '?' will be used.
enum Opts {
  OPT_P, OPT_N, OPT_O, OPT_OC, OPT_ORDER, 
  OPT_DNA, OPT_RNA, OPT_PROTEIN, OPT_ALPH,
  OPT_HOFRACT, OPT_SEED, OPT_CCUT, OPT_HELP, OPT_VERBOSITY
};

/***********************************************************************
 Process command line options and return the command line string.
 ***********************************************************************/
static void process_command_line(int argc, char* argv[], STREME_OPTIONS_T *options) {
  int i, j;
  struct option streme_options[] = {
    {"p", required_argument, NULL, OPT_P},
    {"n", required_argument, NULL, OPT_N},
    {"order", required_argument, NULL, OPT_ORDER},
    {"o", required_argument, NULL, OPT_O},
    {"oc", required_argument, NULL, OPT_OC},
    {"dna", no_argument, NULL, OPT_DNA},
    {"rna", no_argument, NULL, OPT_RNA},
    {"protein", no_argument, NULL, OPT_PROTEIN},
    {"alph", required_argument, NULL, OPT_ALPH},
    {"hofract", required_argument, NULL, OPT_HOFRACT},
    {"seed", required_argument, NULL, OPT_SEED},
    {"ccut", required_argument, NULL, OPT_CCUT},
    {"help", no_argument, NULL, OPT_HELP},
    {"verbosity", required_argument, NULL, OPT_VERBOSITY},
    {NULL, 0, NULL, 0} //boundary indicator
  };

  // Set options to defaults.
  memset(options, 0, sizeof(STREME_OPTIONS_T));
  options->output_dirname = DEFAULT_OUTPUT_DIRNAME;
  options->allow_clobber = True;
  options->posfile = NULL;
  options->negfile = NULL;
  options->order = -1;                  // flags no option given
  options->alphabet_type = DEFAULT_ALPHABET_TYPE;
  options->alph = NULL;
  options->alph_file = NULL;
  options->hofract = DEFAULT_HOFRACT;
  options->seed = DEFAULT_SEED;
  options->ccut = DEFAULT_CCUT;

  // Process arguments.
  while (1) {
    int opt = getopt_long_only(argc, argv, "", streme_options, NULL);
    if (opt == -1) break;
    switch (opt) {
      case OPT_O: // Set output directory with no clobber
        options->output_dirname = optarg;
        options->allow_clobber = false;
        options->text_only = False;
        break;
      case OPT_OC: // Set output directory with clobber
        options->output_dirname = optarg;
        options->allow_clobber = True;
        options->text_only = False;
        break;
      case OPT_P:
        options->posfile = optarg;
        break;
      case OPT_N:
        options->negfile = optarg;
        break;
      case OPT_ORDER:
        options->order = atoi(optarg);
        break;
      case OPT_DNA:
        options->alphabet_type = Dna;
        break;
      case OPT_RNA:
        options->alphabet_type = Rna;
        break;
      case OPT_PROTEIN:
        options->alphabet_type = Protein;
        break;
      case OPT_ALPH:
        options->alphabet_type = Custom;
        options->alph_file = optarg;
        break;
      case OPT_HOFRACT:
        options->hofract = atof(optarg);
        break;
      case OPT_SEED:
        options->seed = atof(optarg);
        break;
      case OPT_CCUT:
        options->ccut = atoi(optarg);
        break;
      case OPT_VERBOSITY:
        verbosity = atoi(optarg);
        break;
      case OPT_HELP:
      case '?':
        usage(NULL);
        break;
      default: // just in case we forget to handle a option
        die("Unhandled option %d", opt);
    }
  }

  // Check that the input is valid.
  if (options->posfile == NULL) {
    usage("You must supply a FASTA file with the primary sequences.");
  }

  // Set the negfile = posfile if it is null and we need it.
  if (options->negfile == NULL) {
    options->negfile = options->posfile;
  }

  // Check that order is not too large.
  if (options->order > MAX_ORDER) {
    usage("The value of --order must be in the range [0,..,%d].  order = %d", MAX_ORDER, options->order);
  }

  // Set the MEME-style alphabet and set defaults if they were not given for --order.
  if (options->alphabet_type == Dna) {
    options->alph = alph_dna();
    if (options->order == -1) options->order = DEFAULT_DNA_ORDER;
  } else if (options->alphabet_type == Rna) {
    options->alph = alph_rna();
    if (options->order == -1) options->order = DEFAULT_RNA_ORDER;
  } else if (options->alphabet_type == Protein) {
    options->alph = alph_protein();
    if (options->order == -1) options->order = DEFAULT_PROT_ORDER;
  } else if (options->alphabet_type == Custom) {
    options->alph = alph_load(options->alph_file, True); // load custom alphabet
    // Die if a MEME-style alphabet did not load successfully.
    if (options->alph == NULL) exit(EXIT_FAILURE);
    if (options->order == -1) options->order = DEFAULT_CUSTOM_ORDER;
  }

  // Check the verbosity level.
  if (verbosity < 1 || verbosity > 5) {
    usage("The verbosity level must be in the range [1, ..., 5]. verbosity = %d", verbosity);
  }

} // process_command_line

//
// Print FASTA files.  
// Files will be empty if multiseq is NULL.
//
static void print_pos_neg_fasta(
  Multiseq *multiseq,		// multiseq object (or NULL)
  char *output_dirname,		// directory for output files
  char *pos_filename,		// file name for positive sequences
  char *neg_filename,		// file name for negative sequences
  BOOL shuffled_neg,		// negative file is shuffled sequences
  int ccut			// centrally trim positive sequences to this (maximum) length
) {
  int i, seqno;
  FILE *pos_fasta;
  FILE *neg_fasta;
  Uint npos = multiseq ? multiseq->npos : 0;
  Uint nneg = multiseq ? multiseq->nneg : 0;
  Uint ntot = npos + nneg;

  // Create the path of the positive FASTA output file and open it.
  char *path = make_path_to_file(output_dirname, pos_filename);
  if ((pos_fasta = fopen(path, "w")) == NULL) {
    fprintf(stderr, "ERROR: Unable to open primary sequence FASTA file '%s' for writing.\n", path);
    exit(EXIT_FAILURE);
  }
  myfree(path);

  // Create the path of the negative FASTA output file and open it.
  path = make_path_to_file(output_dirname, neg_filename);
  if ((neg_fasta = fopen(path, "w")) == NULL) {
    fprintf(stderr, "ERROR: Unable to open control sequence FASTA file '%s' for writing.\n", path);
    exit(EXIT_FAILURE);
  }
  myfree(path);

  // Output the sequences in the multiseq object.
  for (seqno=0; seqno<ntot; seqno++) {
    // Print the FASTA header
    Uchar *desc = DESCRIPTIONPTR(multiseq, seqno);
    fprintf(seqno < npos ? pos_fasta : neg_fasta, ">%s%s\n", 
      desc, (shuffled_neg && seqno >= npos) ? "_shuf" : "");
    // Print the FASTA sequence in lines of 50 characters.
    Uint seqstart = multiseq->seqstarts[seqno];
    Uint seqlen = multiseq->seqlengths[seqno];
    // Trim positive sequences to central region of size ccut.
    if (seqno < npos && ccut > 0 && seqlen > ccut) {
      seqstart += (seqlen - ccut) / 2;		// rounds down
      seqlen = ccut;
    }
    for (i=0; i<seqlen; i+=50) {
      int len = MIN(50, seqlen - i);
      fprintf(seqno < npos ? pos_fasta : neg_fasta, "%*.*s\n", len, len, multiseq->sequence + seqstart + i);
    }
  }
   
} // print_pos_neg_fasta

/*
  This is the main function.
*/
int main(int argc, char *argv[])
{
  STREME_OPTIONS_T options;
  char *train_pos = "train_pos.fa";
  char *train_neg = "train_neg.fa";
  char *holdout_pos = "holdout_pos.fa";
  char *holdout_neg = "holdout_neg.fa";

  // command line processing
  process_command_line(argc, argv, &options);

  srand_mt(options.seed);		// random seed for shuffling sequences
  set_randfunc((randfunc_t) random_mt); // for ushuffle

  // Create the output directory.
  if (create_output_directory(options.output_dirname, options.allow_clobber, False)) {
    exit(EXIT_FAILURE);
  }

  // Input the sequences, converting non-core characters to SEPARATOR.
  // Also reads or creates the background model.
  options.nmotifs = 1;				// suppress warning message
  options.objfun = NO_OBJFUN;			// create negatives if not given
  options.minwidth = 0;				// minimum allowed sequence length
  BOOL do_rc = False;				// don't create reverse complements
  BOOL no_trim = True;				// don't trim the sequences to average length
  BOOL set_back = False;			// don't set the background in the objects
  ARRAY_T *back = NULL;				// space for the MEME-style background
  Multiseq *holdout_multiseq = NULL;
  Multiseq *train_multiseq = read_pos_neg_seqs(&options, do_rc, True, no_trim, set_back, back, &holdout_multiseq);
  free_array(back);

  // 
  // Output the FASTA sequences.
  //
  BOOL shuffled_neg = options.negfile == options.posfile;
  print_pos_neg_fasta(train_multiseq, options.output_dirname, train_pos, train_neg, shuffled_neg, options.ccut);
  if (options.hofract>0) print_pos_neg_fasta(holdout_multiseq, options.output_dirname, holdout_pos, holdout_neg, shuffled_neg, options.ccut);

  return(EXIT_SUCCESS);
} // main
