//
// STREME -- Sensitive, Thorough, Rapid, Enriched Motif Elicitation
//
#define DEFINE_GLOBALS
#include "streme-utils.h"
#include "citation.js.h"

// Nasty global variable for compare_node_score().
// We'll need to change this for multithreading.
int CURRENT_WIDTH;

static char *DEFAULT_OUTPUT_DIRNAME = "streme_out";

/*****************************************************************************
 * Print a usage message with an optional reason for failure.
 *****************************************************************************/
static void usage(char *format, ...) {
  va_list argp;
  char *usage =
    "\n"
    "Usage: streme [options] --p <primary sequences>\n"
    "\n"
    "   Options:\n"
    "     --p <filename>            primary (positive) sequence file name (required)\n"
    "     --n <filename>            control (negative) sequence file name;\n"
    "                               default: if --n is not given, then STREME\n"
    "                               creates control sequences by shuffling each of\n"
    "                               the primary sequences preserving the positions\n"
    "                               of non-core characters and k-mer frequencies\n"
    "                               (see --order, below; ignored if --objfun cd given)\n"
    "     --order <m>               estimates an m-order background model for scoring\n"
    "                               sites and uses an m-order shuffle if creating\n"
    "                               control sequences from primary sequences;\n"
    "                               default: %d (DNA), %d (RNA), %d (Protein), %d (custom)\n"
    "     --kmer <m>                [deprecated: use --order instead]\n"
    "     --bfile <bfile>           use the background model contained in bfile instead\n"
    "                               of creating it from the control sequences;\n"
    "                               default: see --order\n"
    "     --objfun de|cd            objective function to optimize in motif discovery\n"
    "                                 de : Differential Enrichment\n"
    "                                 cd : Central Distance\n"
    "                               default: de\n"
    "     --o <output_dir>          output directory; default: '%s'\n"
    "     --oc <output_dir>         allow overwriting; default: '%s'\n"
    "     --text                    output text only; overrides --o and --oc;\n"
    "                               default: create text, HTML, TSV and XML files in <output_dir>\n"
    "     --dna                     sequences use standard DNA alphabet (default)\n"
    "     --rna                     sequences use standard RNA alphabet\n"
    "     --protein                 sequences use standard protein alphabet\n"
    "     --alph <alph_file>        sequences use alphabet defined in <alph_file>;\n"
    "                               converts to uppercase unless both cases in core\n"
    "     --minw <minwidth>         minimum width for motifs (must be >= %d); \n"
    "                               default: %d\n"
    "     --maxw <maxwidth>         maximum width for motifs (must be <= %d);\n"
    "                               default: %d\n"
    "     --w <w>                   sets <minwidth> and <maxwidth> to <w> (must be <= %d);\n"
    "                               default: see --minw and --maxw\n"
    "     --neval <neval>           evaluate <neval> seeds of each width;\n"
    "                               default: %d\n"
    "     --nref <nref>             refine <nref> evaluated seeds of each width;\n"
    "                               nref==0 means just evaluate single best seed;\n"
    "                               default: %d\n"
    "     --niter <niter>           iterate refinement at most <niter> times per seed;\n"
    "                               default: %d\n"
/*
    "     --pvt <pvt>               stop if hold-out set p-value greater than <pvt>\n"
    "                               (see --patience and --hofract, below);\n"
    "                               overrides --nmotifs;\n"
    "                               default: %g\n"
*/
    "     --thresh <thresh>         significance threshold for reporting enriched motifs;\n"
    "                               default: p-value= %g (%g if --evalue given)\n"
    "     --evalue                  use p-value significance threshold; default: p-value\n"
    "     --patience <patience>     quit after <patience> consecutive motifs exceed <thresh>;\n"
    "                               default: %d\n"
    "     --nmotifs <nmotifs>       stop if <nmotifs> motifs have been output;\n"
    "                               overrides --thresh if > 0;\n"
    "                               default: quit when new motif significance exceeds <thresh>\n"
    "     --time <t>                quit before <t> CPU seconds consumed;\n"
    "                               default: no time limit\n"
    "     --totallength <len>       truncate each sequence set to length <len>;\n"
    "                               default: 0 (do not truncate)\n"
    "     --hofract <hofract>       fraction of sequences in hold-out set;\n"
    "                               default: %g\n"
#ifdef EXP
    "     --max_pal_ed <mpe>        consider a palindromic model if the edit distance\n"
    "                               between the model and its reverse complement\n"
    "                               is less than <mpe>;\n" 
    "                               default: %f\n"
    "     --min_pal_ratio <mpr>     choose the palindromic model if the ratio of its\n"
    "                               log(pvalue) to the model's is greater than <mpr>;\n"
    "                               default: %f\n"
    "     --usepv                   choose best model using training set p-value;\n"
    "     --useer                   choose best model using training set enrichment ratio;\n"
    "                               default: choose best model using training set p-value\n"
    "     --minscore                minimum score to consider during approximate matching;\n"
    "                               default: %g\n"
    "     --ignore_depth            ignore minimum score below this depth;\n"
    "                               default: %d\n"
    "     --nsubsets                number of nested subsets to use in refinement;\n"
    "                               default: %d\n"
    "     --cand                    print all candidate motifs of each width\n"
    "                               to TEXT output only;\n"
    "                               default: print only final motifs\n"
#endif
    "     --seed <seed>             random seed for shuffling sequences;\n"
    "                               default: %d\n"
    "     --align left|center|right align sequences left/center/right for site\n"
    "                               positional distribution plots; default: center\n"
    "     --desc <desc>             include this description text in HTML\n"
    "     --dfile <dfile>           include contents of this description file in HTML,\n"
    "                               overrides --desc\n"
    "     --help                    print this message and exit\n"
    "     --version                 print the program version and exit\n"
    "     --verbosity 1|2|3|4|5     level of diagnostic output (default: %d)\n"
    "                               1: none 2: helpful 3: debug 4: tons 5: ludicrous\n"
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
    DEFAULT_OUTPUT_DIRNAME, DEFAULT_OUTPUT_DIRNAME, 
    MIN_WIDTH, DEFAULT_MINWIDTH, MAX_WIDTH, DEFAULT_MAXWIDTH, MAX_WIDTH,
    DEFAULT_NEVAL, DEFAULT_NREF, DEFAULT_NITER, DEFAULT_PVT, DEFAULT_EVT, DEFAULT_PATIENCE,
    DEFAULT_HOFRACT, 
#ifdef EXP
    DEFAULT_MAX_PAL_ED, DEFAULT_MIN_PAL_RATIO,
    DEFAULT_MINSCORE, DEFAULT_IGNORE_DEPTH, DEFAULT_NSUBSETS,
#endif
    DEFAULT_SEED, DEFAULT_VERBOSITY
  );
  fflush(stderr);
  exit(EXIT_FAILURE);
} // usage

// An easy way of assigning a number to each option.
// Be careful that there are not more than 62 options as otherwise the value
// of '?' will be used.
enum Opts {
  OPT_P, OPT_N, OPT_ORDER, OPT_KMER, OPT_BFILE, OPT_O, OPT_OC, OPT_TEXT, OPT_OBJFUN,
  OPT_DNA, OPT_RNA, OPT_PROTEIN, OPT_ALPH,
  OPT_MINWIDTH, OPT_MAXWIDTH, OPT_WIDTH, OPT_NEVAL, OPT_NREF, OPT_NITER, OPT_NMOTIFS, 
  OPT_TIME, OPT_TOTALLENGTH, OPT_PVT, OPT_THRESH, OPT_EVALUE, OPT_PATIENCE, OPT_MPR, OPT_MPE, OPT_HOFRACT, 
  OPT_USEER, OPT_USEPV, OPT_CAND, OPT_SEED, OPT_MINSCORE, OPT_IGNORE_DEPTH,
  OPT_NSUBSETS, OPT_ALIGN, OPT_DESC, OPT_DFILE, OPT_VERBOSITY, OPT_VERSION, OPT_HELP
};

/***********************************************************************
 Process command line options and return the command line string.
 ***********************************************************************/
static char *process_command_line(int argc, char* argv[], STREME_OPTIONS_T *options) {
  int i, j;
  struct option streme_options[] = {
    {"p", required_argument, NULL, OPT_P},
    {"n", required_argument, NULL, OPT_N},
    {"order", required_argument, NULL, OPT_ORDER},
    {"kmer", required_argument, NULL, OPT_KMER},
    {"bfile", required_argument, NULL, OPT_BFILE},
    {"o", required_argument, NULL, OPT_O},
    {"oc", required_argument, NULL, OPT_OC},
    {"text", no_argument, NULL, OPT_TEXT},
    {"objfun", required_argument, NULL, OPT_OBJFUN},
    {"dna", no_argument, NULL, OPT_DNA},
    {"rna", no_argument, NULL, OPT_RNA},
    {"protein", no_argument, NULL, OPT_PROTEIN},
    {"alph", required_argument, NULL, OPT_ALPH},
    {"minw", required_argument, NULL, OPT_MINWIDTH},
    {"maxw", required_argument, NULL, OPT_MAXWIDTH},
    {"w", required_argument, NULL, OPT_WIDTH},
    {"neval", required_argument, NULL, OPT_NEVAL},
    {"nref", required_argument, NULL, OPT_NREF},
    {"niter", required_argument, NULL, OPT_NITER},
    {"nmotifs", required_argument, NULL, OPT_NMOTIFS},
    {"time", required_argument, NULL, OPT_TIME},
    {"totallength", required_argument, NULL, OPT_TOTALLENGTH},
    {"pvt", required_argument, NULL, OPT_PVT},	// for backward compatibility
    {"thresh", required_argument, NULL, OPT_THRESH},
    {"evalue", no_argument, NULL, OPT_EVALUE},
    {"patience", required_argument, NULL, OPT_PATIENCE},
    {"mpr", required_argument, NULL, OPT_MPR},
    {"mpe", required_argument, NULL, OPT_MPE},
    {"hofract", required_argument, NULL, OPT_HOFRACT},
    {"useer", no_argument, NULL, OPT_USEER},
    {"usepv", no_argument, NULL, OPT_USEPV},
    {"cand", no_argument, NULL, OPT_CAND},
    {"seed", required_argument, NULL, OPT_SEED},
    {"minscore", required_argument, NULL, OPT_MINSCORE},
    {"ignore_depth", required_argument, NULL, OPT_IGNORE_DEPTH},
    {"nsubsets", required_argument, NULL, OPT_NSUBSETS},
    {"align", required_argument, NULL, OPT_ALIGN},
    {"desc", required_argument, NULL, OPT_DESC},
    {"dfile", required_argument, NULL, OPT_DFILE},
    {"version", no_argument, NULL, OPT_VERSION},
    {"help", no_argument, NULL, OPT_HELP},
    {"verbosity", required_argument, NULL, OPT_VERBOSITY},
    {NULL, 0, NULL, 0} //boundary indicator
  };

  // Set options to defaults.
  memset(options, 0, sizeof(STREME_OPTIONS_T));
  options->output_dirname = DEFAULT_OUTPUT_DIRNAME;
  options->allow_clobber = True;
  options->text_only = false;
  options->posfile = NULL;
  options->negfile = NULL;
  options->objfun = DEFAULT_OBJFUN;
  options->order = -1;			// flags no option given
  options->kmer = -1;			// flags no option given
  options->bfile = NULL;		// flags no option given
  options->alphabet_type = DEFAULT_ALPHABET_TYPE;
  options->alph = NULL;
  options->alph_file = NULL;
  options->minwidth = DEFAULT_MINWIDTH;
  options->maxwidth = DEFAULT_MAXWIDTH;
  options->width = 0;
  options->neval = DEFAULT_NEVAL;
  options->nref = DEFAULT_NREF;
  options->niter = DEFAULT_NITER;
  options->nmotifs = 0;			// use PVT by default
  options->thresh_type = PVALUE;
  options->thresh = DEFAULT_PVT;
  options->patience = DEFAULT_PATIENCE;
  options->min_pal_ratio = DEFAULT_MIN_PAL_RATIO;
  options->max_pal_ed = DEFAULT_MAX_PAL_ED;
  options->time = DEFAULT_TIME;
  options->totallength = 0;		// don't truncate input sets
  options->hofract = DEFAULT_HOFRACT;
  options->usepv = true;
  options->cand = false;
  options->minscore = DEFAULT_MINSCORE;
  options->ignore_depth = DEFAULT_IGNORE_DEPTH;
  options->nsubsets = DEFAULT_NSUBSETS;
  options->seed = DEFAULT_SEED;
  options->align_txt = NULL;
  options->description = NULL;
  options->dfile = NULL;
  options->text_output = NULL;
  options->commandline = NULL;
# ifdef EXP
  options->experimental = true;
#else
  options->experimental = false;
# endif

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
      case OPT_TEXT: // Output text only
        options->text_only = True;
        break;
      case OPT_P:
        options->posfile = optarg;
        break;
      case OPT_N:
        options->negfile = optarg;
        break;
      case OPT_OBJFUN:
        if (! strcmp(optarg, "de")) {
          options->objfun = DE;
        } else if (! strcmp(optarg, "cd")) {
          options->objfun = CD;
        } else {
          usage("Unknown value for --objfun (%s)\n", optarg);
        }
        break;
      case OPT_ORDER:
        options->order = atoi(optarg);
        options->kmer = options->order+1;
        break;
      case OPT_KMER:
        options->kmer = atoi(optarg);
        options->order = options->kmer - 1;
        DEBUG_FMT(QUIET_VERBOSE, "# Warning: --kmer is deprecated; setting --order to %d\n", options->order);
        break;
      case OPT_BFILE:
        options->bfile = optarg;
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
      case OPT_MINWIDTH:
        options->minwidth = atoi(optarg);
        break;
      case OPT_MAXWIDTH:
        options->maxwidth = atoi(optarg);
        break;
      case OPT_WIDTH:
        options->width = atoi(optarg);
        break;
      case OPT_NEVAL:
        options->neval = atoi(optarg);
        break;
      case OPT_NREF:
        options->nref = atoi(optarg);
        break;
      case OPT_NITER:
        options->niter  = atoi(optarg);
        break;
      case OPT_NMOTIFS:
        options->nmotifs = atoi(optarg);
        if (options->nmotifs > 0) options->thresh = DEFAULT_PVT;
        break;
      case OPT_PVT:		// for backwards compatibility
        options->thresh = atof(optarg);
        options->thresh_type = PVALUE;
        options->nmotifs = 0;
        break;
      case OPT_THRESH:
        options->thresh = atof(optarg);
        options->nmotifs = 0;
        break;
      case OPT_EVALUE:
        options->thresh_type = EVALUE;
        options->nmotifs = 0;
        break;
      case OPT_PATIENCE:
        options->patience = atoi(optarg);
        break;
      case OPT_TIME:
        options->time = (int) atof(optarg);
        break;
      case OPT_TOTALLENGTH:
        options->totallength = atof(optarg);
        break;
      case OPT_HOFRACT:
        options->hofract = atof(optarg);
        break;
      case OPT_MPR:
	options->min_pal_ratio = atof(optarg);
        break;
      case OPT_MPE:
	options->max_pal_ed = atof(optarg);
        break;
      case OPT_USEER:
        options->usepv = False;
        break;
      case OPT_USEPV:
        options->usepv = True;	// hidden option
        break;
      case OPT_MINSCORE:
        options->minscore = atof(optarg);
        break;
      case OPT_IGNORE_DEPTH:
        options->ignore_depth = atoi(optarg);
        break;
      case OPT_NSUBSETS:
        options->nsubsets = atoi(optarg);
        break;
      case OPT_CAND:
        options->cand = True;
        break;
      case OPT_SEED:
        options->seed = atoi(optarg);
        break;
      case OPT_ALIGN:
        options->align_txt = optarg;
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
      case OPT_HELP:
      case '?':
        usage(NULL);
        break;
      default: // just in case we forget to handle a option
        die("Unhandled option %d", opt);
    }
  }

  // Override --minw and --maxw if --w given
  if (options->width > 0) {
    options->minwidth = options->maxwidth = options->width;
  }

  // Check that the input is valid.
  if (options->posfile == NULL) {
    usage("You must supply a FASTA file with the primary sequences.");
  }
  if (options->minwidth < MIN_WIDTH) {
    usage("The minimum allowed motif width is %d.  minwidth = %d", MIN_WIDTH,
      options->minwidth);
  }
  if (options->maxwidth > MAX_WIDTH) {
    usage("The maximum allowed motif width is %d.  maxwidth = %d", MAX_WIDTH, 
      options->maxwidth);
  }
  if (options->maxwidth < options->minwidth) {
    usage("You must set minwidth < maxwidth.  minwidth = %d, maxwidth = %d\n",
      options->minwidth, options->maxwidth);
  }
  if (options->neval < 1) {
    usage("The value of --neval must be positive.  neval = %d",
      options->neval);
  }
  if (options->thresh_type == EVALUE) {
    // EVALUE
    if (options->thresh == DEFAULT_PVT) options->thresh = DEFAULT_EVT;
    if (options->thresh <= 0) {
      usage("The value of --thresh must greater than 0. thresh = %g\n", options->thresh);
    }
  } else {
    // PVALUE
    if (options->thresh <= 0 || options->thresh > 1) {
      usage("The value of --thresh must be in the range (0,..1] when you use --evalue. thresh = %g\n", options->thresh);
    }
  }
  if (options->nmotifs < 0) {
    usage("The value of --nmotifs may not be negative.  nmotifs = %f", options->nmotifs);
  }
  if (options->time < 0) {
    usage("The value of --time may not be negative.  time = %f", options->time);
  }
  if (options->totallength < 0) {
    usage("The value of --totallength may not be negative.  totallength = %d", options->totallength);
  }

  // Checks and warnings for CD objective function.
  if (options->objfun == CD) {
    if (options->negfile != NULL) {
      DEBUG_FMT(QUIET_VERBOSE, "# Warning: ignoring negative sequence file (%s) with --objfun cd.\n", options->negfile);
      options->negfile = NULL;
    }
  }

  // Set the negfile = posfile if it is null and we need it.
  if (options->objfun != CD && options->negfile == NULL) {
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
  options->kmer = options->order+1;

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

  // Make enough space for all the command line options, with one space between each
  int line_length = 2;		// include space for quotes around --desc text.
  char *pgm_name = basename(argv[0]);
  for (i = 0; i < argc; i++) line_length += strlen(i == 0 ? pgm_name : argv[i]);
  // add on argc to allow one char per word for separating space + terminal '\0'
  char *commandline = (char*) malloc(sizeof(char)*(line_length+argc));
  int nextpos = 0;
  for (i = 0; i < argc; i++) {
    // been here before? put in a space before adding the next word
    if (nextpos) {
      commandline[nextpos] = ' ';
      nextpos++;
    }
    char *nextword = (i == 0) ? pgm_name : argv[i];
    if (i > 0 && (strncmp(argv[i-1], "-desc", 5)==0 || strncmp(argv[i-1], "--desc", 6)==0)) {
      commandline[nextpos++] = '\'';
      strcpy(&commandline[nextpos], nextword);
      nextpos += strlen(nextword);
      commandline[nextpos++] = '\'';
    } else {
      strcpy(&commandline[nextpos], nextword);
      nextpos += strlen(nextword);
    }
  }

  return(commandline);
} // process_command_line

//
// Print the STREME output header.
//
void print_streme_header(
  ALPH_T *alph,
  BOOL do_rc,
  double *background,		// 0-order background
  FILE *stream
) {
  int i; 
  int alen = alph_size_core(alph);

  char *program = "STREME - Sensitive, Thorough, Rapid, Enriched Motif Elicitation"; 
  char *info = "For further information on how to interpret these results please access " SITE_URL ".\n" 
    "To get a copy of the MEME Suite software please access " SOURCE_URL ".\n"; 
  banner(program, info, STREME_CITE, stream); 

  if (alph_is_builtin_dna(alph)) { 
    fprintf(stream, "ALPHABET= %s\n", DNA); 
  } else if (alph_is_builtin_rna(alph)) { 
    fprintf(stream, "ALPHABET= %s\n", RNA); 
  } else if (alph_is_builtin_protein(alph)) {
    fprintf(stream, "ALPHABET= %s\n", PROTEIN);
  } else {
    PSTARS(stream);
    alph_print_header(alph, stream);
    PSTARS(stream);
    alph_print(alph, false, stream);
    PSTARS(stream);
  }
  fprintf((stream), "\n");
  fprintf((stream), "strands: %s\n\n", (do_rc) ? "+ -" : "+");
  fprintf((stream), "Background letter frequencies\n");
  for (i=0; i<alen; i++) {
    fprintf((stream), "%c %0.3g ", I2A(i), background[i]);
  }
  fprintf((stream), "\n\n");
  fflush(stream);
} // print_streme_header

//
// Output the SEQUENCES_TSV results.
//
static void print_streme_passing_sequences_tsv(
  STREME_OPTIONS_T *options,    // STREME options
  char *tsv_path,		// path to TSV output file
  Model **models,               // the discovered motifs
  BOOL have_holdout,		// true if there was a holdout set
  int n_out_motifs              // the number of motifs to output
) {
  int i;
  int imotif;			// motif number

  // Create the TSV output file.
  FILE *outfile;
  if ((outfile = fopen(tsv_path, "w")) == NULL) {
    DEBUG_FMT(QUIET_VERBOSE, "# Warning: Unable to open SEQUENCES TSV file \"%s\" for writing.\n", tsv_path);
    return;
  }

  // Print sequences for each model.
  for (imotif=0; imotif<n_out_motifs; imotif++) {
    Model *model = models[imotif];
    // Print the TSV header for first motif.
    if (imotif == 0) {
      fprintf(outfile,
	"%s\t"
	"%s\t"
	"%s\t"
	"%s\t"
	"%s\t"
	"%s\t"
	"%s\n",
	"motif_ID", "motif_ALT_ID", (have_holdout ? "motif_P-value" : "motif_Score"), "seq_ID", "seq_Score", "seq_Class", "is_holdout?"
      );
    } else {
      fprintf(outfile, "#\n");
    }

    // Print the sorted passing sequences.
    double m1, e1, prec=1;
    if (have_holdout) {
      exp10_logx(model->test_log_pvalue/log(10), m1, e1, prec); \
    } else {
      exp10_logx(model->train_log_pvalue/log(10), m1, e1, prec); \
    }
    MOTIF_T *motif = model->motif;
    Passing_seq *sequences = model->passing_seqs;
    int npassing = model->npassing;
    for (i=0; i<npassing; i++) {
      fprintf(outfile,
        "%d-%s\t"	// Motif ID
        "STREME-%d\t"	// Motif Alternate ID
	"%3.1fe%+04.0f\t"	// Motif P-value or Score
	"%s\t"		// seq_ID
	"%.2f\t"	// seq_Score
	"%s\t"		// seq_Class
	"%d\n",		// is_holdout?
	imotif+1,	// <count>-<consensus>
	model->consensus,		
	imotif+1,	// STREME-<count>
	m1, e1,
	sequences[i].desc, sequences[i].score, (sequences[i].is_tp ? "tp" : "fp"),
	sequences[i].dbno
      );
    } // passing sequence
  } // model
} // print_streme_passing_sequences_tsv

//
// Add counts to suffix tree: process a leaf.
//
Sint cnt_processleaf(Uint leafnumber, Bref lcpnode, void *info) {
  Countstate *cnt_state = (Countstate *) info;
  double avg_poslen = cnt_state->multiseq->avg_poslen;
  Branchinfo branchinfo;
  Uint nodeheadposition = leafnumber;

  // Get the mother of this leaf.
  Bref *motherptrptr = TOPARRAY(&(cnt_state->mother), Bref);
  Uint motherindex = BRADDR2NUM(cnt_state->stree, *motherptrptr);
  getbranchinfostree(cnt_state->stree, ACCESSDEPTH | ACCESSHEADPOS, &branchinfo, *motherptrptr);
  Uint motherdepth = branchinfo.depth;

  // Get the not-too-deep ancestor of this leaf.
  Bref *ancestorptrptr = TOPARRAY(&(cnt_state->ancestor), Bref);
  getbranchinfostree(cnt_state->stree, ACCESSHEADPOS, &branchinfo, *ancestorptrptr);
  Uint ancestordepth = branchinfo.depth;
  Uint ancestorheadposition = branchinfo.headposition;

  // Get the sequence corresponding to this leaf and set the bit in _seqno_bittab.
  // Keep track of the number of bits currently set.
  Uint seqno = -1;
  Uint pos = -1;
  Uint original_seqno = -1;
  BOOL on_revcomp = False;
  BOOL is_negative = False;
  // Ancestor is not the root and we have a valid leaf.
  if (ancestorptrptr != NULL && leafnumber < cnt_state->stree->textlen) {
    PairUint pp;
    // Convert leafnumber (position in multiseq) to position pair.
    if (pos2pospair(cnt_state->multiseq, &pp, leafnumber) >= 0) {	// bogus result if leafnumber == totallen
      seqno = pp.uint0; 
      pos = pp.uint1;
      Uint subtreeindex = cnt_state->stree->subtreedataindex[ancestorheadposition]; 
      Subtreedata *ancestorstdptr = PEEKARRAY(&(cnt_state->stree->subtreedata), subtreeindex, Subtreedata);

      // Store the positive or negative sequence number in the ancestor node subtree data.
      Uint dummy_seqno;
      original_seqno = seqno;
      GET_SEQUENCE_NUMBER_AND_TYPE(cnt_state->multiseq, seqno, on_revcomp, is_negative, seqno, dummy_seqno);
      Uint npos = cnt_state->multiseq->npos;
      Uint nneg = cnt_state->multiseq->nneg;
      if (! is_negative) {
        // positive sequence
        if (! ISIBITSET(ancestorstdptr->pos_seqno_bittab, seqno)) {
	  SETIBIT(ancestorstdptr->pos_seqno_bittab, seqno);
          if (cnt_state->objfun == CD) {
	    // Add the distance from site to sequence center to the ancestor's sum.
	    // Note that we are assumming all sequences have the same length.
	    // CD: note that only the one site's distance gets counted per sequence.
	    ancestorstdptr->dtc_sum += fabs(pos - avg_poslen/2);
	  }
        }
      } else {
        // negative sequence
        if (! ISIBITSET(ancestorstdptr->neg_seqno_bittab, seqno-npos)) {
	  SETIBIT(ancestorstdptr->neg_seqno_bittab, seqno-npos);
        }
      }

      // Store an entry in the subtreedata if this leaf should
      // be treated like a valid node.
      if (motherdepth < cnt_state->maxdepth) {
	Uchar *head = cnt_state->multiseq->sequence + nodeheadposition;
	// Get the minimum and maximum length of valid words ending in the leaf.
	Uint mindepth = motherdepth+1;
	Uint maxdepth, i;
        Uint totallength = cnt_state->multiseq->totallength;
	// See how deep the first separator is; that is the maxdepth.
	for (i=nodeheadposition, maxdepth=0; i<totallength && maxdepth<cnt_state->maxdepth; maxdepth++, i++) {
	  if (head[maxdepth] == SEPARATOR) break;
	}
	// Create a subtreedata node containing the data for this leaf if it is valid.
	if (maxdepth >= mindepth) {
	  Subtreedata std;
	  std.head = head;
	  std.min_width = mindepth;
	  std.max_width = maxdepth;
	  std.pos_count = is_negative ? 0 : 1;
	  std.neg_count = is_negative ? 1 : 0;
	  std.pos_seqno_bittab = NULL;
	  std.neg_seqno_bittab = NULL;
          std.pos_seqno_list = (Uint *) malloc(sizeof(Uint));
          std.neg_seqno_list = (Uint *) malloc(sizeof(Uint));
          std.pos_seqno_list[0] = is_negative ? -1 : seqno;
          std.neg_seqno_list[0] = is_negative ? seqno-npos : -1;
          std.dtc_sum = 0;
          if (verbosity >= HIGHER_VERBOSE) 
            fprintf(stderr, "cnt_processleaf: pos_count %d npos %d neg_count %d nneg %d bernoulli %f\n",
              std.pos_count, npos, std.neg_count, nneg, cnt_state->multiseq->bernoulli);
	  GET_PVALUE(
	    std.log_pvalue,
	    cnt_state->objfun, 
	    cnt_state->multiseq->use_binomial, 
            std.pos_count, 
	    std.neg_count, 
	    npos, 
	    nneg, 
	    cnt_state->multiseq->bernoulli,
	    std.dtc_sum, 
	    maxdepth, 
	    avg_poslen, 
	    "cnt_processleaf",
            cnt_state->pv_cache,
            cnt_state->cache_length
	  );
	  std.seqno = original_seqno;
	  std.seqpos = pos;
	  if (cnt_state->multiseq->do_rc) {
	    get_rc_head(std.rc_head, cnt_state->multiseq, original_seqno, pos, maxdepth);
	  } else {
	    std.rc_head = NULL;
	  }
          std.scores = NULL;
	  // Push the data for the leaf onto the subtree data stack.
	  PUSHARRAY(&(cnt_state->stree->subtreedata), Subtreedata, 128, std);
	  Uint subtreeindex = TOPINDEXARRAY(&(cnt_state->stree->subtreedata), Subtreedata);
	  // Save the index in the stack in the leaf data index at the leaf number.
	  cnt_state->stree->leafdataindex[leafnumber] = subtreeindex;
	} // valid leaf
 
      } // mother is not too deep

    } // pospair

  } // ancestor is not root

#ifdef PRINTREE
  // Print the branch with the path (label).
  Uint nodedepth = cnt_state->multiseq->totallength - leafnumber;
  Uint labelstart = nodeheadposition + motherdepth;
  Uint labellength = nodedepth - motherdepth;
  Uchar *label = cnt_state->multiseq->sequence + labelstart;
  fprintf(stderr, "B_%u --'%*.*s'--> L_%u", motherindex, labellength, labellength, label, leafnumber);
  fprintf(stderr, " seqno %d pos %d\n", (int) seqno, (int) pos);
#endif

  return(0);
} // cnt_processleaf

//
// Add counts to suffix tree: process a branch1.
//
BOOL cnt_processbranch1(Bref nodeptr, void *info)
{
  Countstate *cnt_state = (Countstate *) info;
  Uint npos = cnt_state->multiseq->npos;
  Uint nneg = cnt_state->multiseq->nneg;
  Branchinfo branchinfo;
  int i;

  // Get the depth and head position for the branch.
  getbranchinfostree(cnt_state->stree, ACCESSDEPTH | ACCESSHEADPOS, &branchinfo, nodeptr);
  Uint nodedepth = branchinfo.depth;
  Uint nodeheadposition = branchinfo.headposition;
  Uchar *head = cnt_state->multiseq->sequence + nodeheadposition;

  // Get the mother of this branch, the mother's index and the mother's info.
  Bref *motherptrptr = TOPARRAY(&(cnt_state->mother), Bref);
  Uint motherindex = BRADDR2NUM(cnt_state->stree, *motherptrptr);
  getbranchinfostree(cnt_state->stree, ACCESSDEPTH, &branchinfo, *motherptrptr);

  // Push this node on the mother stack.  This invalidates motherptrptr (potentially)
  // due to the use of realloc in PUSHARRAY().
  PUSHARRAY(&(cnt_state->mother), Bref, 128, nodeptr);

  // Store this node on the ancestor stack if it is not too deep.
  // The node is too deep to be an ancestor if:
  //   1) the depth of the mother is >= maxdepth
  //   2) or the path to the mother contains a separator
  // The top node in this stack is where we store counts and the bit table.
  Uint motherdepth = branchinfo.depth;
  BOOL node_too_deep = False;
  if (motherdepth >= cnt_state->maxdepth) {
    node_too_deep = True;
  } else { 
    // Look for separator leading to mother.
    Uchar *sequence = cnt_state->multiseq->sequence + nodeheadposition;
    for (i=motherdepth-1; i>=0; i--) {
      if (sequence[i] == SEPARATOR) {
        node_too_deep = True;
        break;
      }
    }
  }
  if (! node_too_deep) {
    PUSHARRAY(&(cnt_state->ancestor), Bref, 128, nodeptr);
    // Initialize the numbers of sequences in subtree and the bit table for this node.
    Subtreedata std;
    std.head = NULL;
    std.min_width = 0;
    std.max_width = 0;
    std.pos_count = 0;
    std.neg_count = 0;
    std.dtc_sum = 0;
    INITBITTAB_TLB(std.pos_seqno_bittab, npos);
    std.neg_seqno_bittab = NULL;
    if (nneg) INITBITTAB_TLB(std.neg_seqno_bittab, nneg);
    std.pos_seqno_list = NULL;
    std.neg_seqno_list = NULL;
    std.scores = NULL;
    // Push the subtree data for the node onto the tree's subtree data stack.
    PUSHARRAY(&(cnt_state->stree->subtreedata), Subtreedata, 128, std);
    // Save the index in the stack in the subtree data index at the node's head position.
    cnt_state->stree->subtreedataindex[nodeheadposition] = TOPINDEXARRAY(&(cnt_state->stree->subtreedata), Subtreedata);
  }

#ifdef PRINTTREE
  // Print the branch.
  // Get the ancestor at the top of the stack.
  Bref *ancestorptrptr = TOPARRAY(&(cnt_state->ancestor), Bref);
  Uint ancestorindex = BRADDR2NUM(cnt_state->stree, *ancestorptrptr);
  Uint labelstart = nodeheadposition + motherdepth;
  Uint labellength = nodedepth - motherdepth;
  Uchar *label = cnt_state->multiseq->sequence + labelstart;
  // Get the index of this node.
  Uint nodeindex = BRADDR2NUM(cnt_state->stree, nodeptr); 
  fprintf(stderr, "B_%u --> B_%u '%*.*s'", motherindex, nodeindex, nodedepth, nodedepth, head);
  fprintf(stderr, " ancestor %u\n", ancestorindex);
#endif

  return(True);
} // cnt_processbranch1

// 
// Add counts to suffix tree: process a branch2.
//
Sint cnt_processbranch2(Bref nodeptr, void *info)
{
  Countstate *cnt_state = (Countstate *) info;
  Multiseq *multiseq = cnt_state->multiseq;
  double avg_poslen = multiseq->avg_poslen;
  BOOL do_rc = multiseq->do_rc;
  Uint npos = multiseq->npos;
  Uint nneg = multiseq->nneg;
  BOOL use_binomial = multiseq->use_binomial;
  double bernoulli  = multiseq->bernoulli;
  Branchinfo branchinfo;
  Uint node_mindepth=0, node_maxdepth=0;
  int i, j, w;
  Uchar *rc_head = NULL;

  // Pop this node off the mother stack.
  (void) POPARRAY(&(cnt_state->mother), Bref);

  // Process this node if it is on the ancestor stack.
  if (*TOPARRAY(&(cnt_state->ancestor), Bref) == nodeptr) {
    // Get the node's head position and depth.
    Uint nodeindex = BRADDR2NUM(cnt_state->stree, nodeptr);
    getbranchinfostree(cnt_state->stree, ACCESSDEPTH | ACCESSHEADPOS, &branchinfo, nodeptr);
    Uint nodedepth = branchinfo.depth;
    Uint nodeheadposition = branchinfo.headposition;
    Uchar *head = multiseq->sequence + nodeheadposition;

    // Compute the numbers of positive and negative sequences in the subtree rooted at the node
    // by counting the numbers of bits that are set.
    Uint subtreeindex = cnt_state->stree->subtreedataindex[nodeheadposition]; 
    Subtreedata *nodestdptr = PEEKARRAY(&(cnt_state->stree->subtreedata), subtreeindex, Subtreedata);
    COUNTBITTAB(nodestdptr->pos_count, nodestdptr->pos_seqno_bittab, npos);
    if (nneg) COUNTBITTAB(nodestdptr->neg_count, nodestdptr->neg_seqno_bittab, nneg);

    // Get the mother of this node and set the node's minimum depth.
    Bref *motherptrptr = TOPARRAY(&(cnt_state->mother), Bref);
    Uint motherdepth = -1;
    if (motherptrptr) {				// Not the root
      getbranchinfostree(cnt_state->stree, ACCESSDEPTH, &branchinfo, *motherptrptr);
      motherdepth = branchinfo.depth;
    }
    node_mindepth = motherdepth + 1;

    //
    // Save the (valid) words corresponding to this node.
    // A node is valid if;
    //   1) it is not the root;
    //   2) it is not deeper than maxdepth
    //   3) no separators up to start of branch
    //   4) it is at least as deep as mindepth
    //
    BOOL node_is_valid = False;

    // Check that node is not root and node is not too deep.
    if (node_mindepth != 0 && node_mindepth <= cnt_state->maxdepth) {
      node_is_valid = True;

      // Check that there is no separator up to start of branch.
      for (i=0; i<=motherdepth && head[i] != SEPARATOR; i++) { /* empty */ }
      if (i <= motherdepth) node_is_valid = False;

      // Get the minimum and maximum depth for valid paths corresponding to this node.
      // Find the position of any separator on the branch label (up to maxdepth).
      for (j=motherdepth; j<nodedepth && j<cnt_state->maxdepth && head[j] != SEPARATOR; j++) { /* empty */ }

      // Set the node's maximum depth to just before the first separator
      // along the path from the head.
      node_maxdepth = j;
      if (node_mindepth > node_maxdepth) node_is_valid = False;

      nodestdptr->min_width = node_mindepth;
      nodestdptr->max_width = node_maxdepth;
      nodestdptr->head = head;

      // Get the sequence number and position of head of valid node and get the
      // corresponding head in the reverse-complement of the sequence if needed.
      if (node_is_valid) {
	PairUint pp;
	Uint seqno = 0, pos = 0;
	if (pos2pospair(multiseq, &pp, nodeheadposition) >= 0) {
	  seqno = pp.uint0; 
	  pos = pp.uint1;
	  if (do_rc) get_rc_head(rc_head, multiseq, seqno, pos, node_maxdepth);
	} else {
	  if (head[0] == SEPARATOR) {
	    if (do_rc) rc_head = head;
	  } else {
	    fprintf(stderr, "ERROR: Problem with pos2pospair in cnt_processbranch2; nodeheadposition %d.\n", nodeheadposition);
            exit(EXIT_FAILURE);
	  }
	}

	// Store the information in the subtree data struct.
	nodestdptr->seqno = seqno;
	nodestdptr->seqpos = pos;
	nodestdptr->rc_head = rc_head;
      } // node_is_valid

      // Compute the p-value of valid node and save it in the index if is not too shallow.
      if (node_is_valid && node_maxdepth >= cnt_state->mindepth) {
	if (verbosity >= HIGHER_VERBOSE) 
	  fprintf(stderr, "cnt_processbranch2: pos_count %d npos %d neg_count %d nneg %d bernoulli %f\n",
	    nodestdptr->pos_count, npos, nodestdptr->neg_count, nneg, cnt_state->multiseq->bernoulli);
        GET_PVALUE(
          nodestdptr->log_pvalue, 
          cnt_state->objfun, 
          use_binomial, 
          nodestdptr->pos_count, 
          nodestdptr->neg_count, 
          npos, 
          nneg, 
          bernoulli, 
          nodestdptr->dtc_sum, 
          node_maxdepth, 
          avg_poslen, 
          "cnt_processbranch2",
          cnt_state->pv_cache,
          cnt_state->cache_length
        );
	PUSHARRAY(&(cnt_state->valid_node_indices), Uint, 128, subtreeindex);
      }

#ifdef PRINTTREE
      if (node_is_valid) {
	// Print all the valid words corresponding to the node.
	// These are the valid words ending on the branch leading to it.
	Uint printdepth;
	for (printdepth=node_mindepth; printdepth<=node_maxdepth; printdepth++) {
	  BOOL print_it = False;
	  if (! do_rc) { 
	    print_it = True;
	    fprintf(stderr, "B_%u head '%*.*s' w %u", nodeindex, printdepth, printdepth, head, printdepth);
	  } else {
	    // Get the reverse complement of the head.
	    // Only print the node with the smaller head as they will have the same counts.
            Uchar *rc_word = rc_head + node_maxdepth - printdepth;
	    if (strncmp((char *)head, (char *)rc_word, printdepth) <= 0) {
	      print_it = True;
	      fprintf(stderr, "B_%u head '%*.*s' '%*.*s'", nodeindex, 
                printdepth, printdepth, head, printdepth, printdepth, rc_word);
	      fprintf(stderr, " w %u", printdepth);
	    }
	  }
	  if (print_it) {
	    fprintf(stderr, " headpos %u", nodeheadposition);
	    fprintf(stderr, " pos_count %d npos %u neg_count %d nneg %u log_pvalue %.2f",
	      nodestdptr->pos_count, npos, nodestdptr->neg_count, nneg, nodestdptr->log_pvalue);
	    fprintf(stderr, "\n");
	  }
	} // printdepth
      }
#endif

    } // Node is not root and not too deep.

    // Pop the node off the ancestor stack.
    (void) POPARRAY(&(cnt_state->ancestor), Bref);
    // Get this node's not-too-deep ancestor and update its sequence bit table.
    if (! EMPTYARRAY(&(cnt_state->ancestor), Bref)) {
      Bref ancestorptr = *TOPARRAY(&(cnt_state->ancestor), Bref);
      Uint ancestorindex = BRADDR2NUM(cnt_state->stree, ancestorptr);
      // Get the ancestor's head position.
      getbranchinfostree(cnt_state->stree, ACCESSHEADPOS, &branchinfo, ancestorptr);
      Uint ancestorheadposition = branchinfo.headposition;
      // Update the bit table of sequences used for this node's not-too-deep ancestor.
      Uint subtreeindex = cnt_state->stree->subtreedataindex[ancestorheadposition]; 
      Subtreedata *ancestorstdptr = PEEKARRAY(&(cnt_state->stree->subtreedata), subtreeindex, Subtreedata);
      double old_pos_matches=0, new_pos_matches=0;
      if (cnt_state->objfun == CD) COUNTBITTAB(old_pos_matches, ancestorstdptr->pos_seqno_bittab, npos);
      ORBITTABS(ancestorstdptr->pos_seqno_bittab, nodestdptr->pos_seqno_bittab, npos);
      if (cnt_state->objfun == CD) COUNTBITTAB(new_pos_matches, ancestorstdptr->pos_seqno_bittab, npos);
      if (nneg) ORBITTABS(ancestorstdptr->neg_seqno_bittab, nodestdptr->neg_seqno_bittab, nneg);

      // CD only: Add the sum of the distances from sites to the sequence centers, 
      // weighted by the number of sequences from the node that were new.
      if (cnt_state->objfun == CD && nodestdptr->pos_count > 0) {
        ancestorstdptr->dtc_sum += ((new_pos_matches-old_pos_matches)/nodestdptr->pos_count) * nodestdptr->dtc_sum;
      }

      // Free the node if it is not valid.
      if (!node_is_valid) {
	FREEBITTAB(nodestdptr->pos_seqno_bittab);
	if (nneg) FREEBITTAB(nodestdptr->neg_seqno_bittab);
	nodestdptr->pos_seqno_bittab = nodestdptr->neg_seqno_bittab = NULL;
      } else {
        // Replace the positive bit table with a list if it is small.
        if (nodestdptr->pos_count <= LISTSIZE) {
          nodestdptr->pos_seqno_list = (Uint *) malloc(nodestdptr->pos_count * sizeof(Uint));
          LISTBITTAB(nodestdptr->pos_seqno_list, nodestdptr->pos_seqno_bittab, npos);
	  FREEBITTAB(nodestdptr->pos_seqno_bittab);
        }
        // Replace the negative bit table with a list if it is small.
        if (nneg && nodestdptr->neg_count <= LISTSIZE) {
          nodestdptr->neg_seqno_list = (Uint *) malloc(nodestdptr->neg_count * sizeof(Uint));
          LISTBITTAB(nodestdptr->neg_seqno_list, nodestdptr->neg_seqno_bittab, nneg);
	  FREEBITTAB(nodestdptr->neg_seqno_bittab);
        }
      }
    } // update bit table
  } // process node on ancestor stack

  return(0);
} // cnt_processbranch2

//
// Process a leaf for finding score-based approximate matches by DFS.
//
Sint sbm_processleaf(Uint leafnumber, Bref lcpnode, void *info) 
{
  Searchstate *srch_state = (Searchstate *) info;
  int w;

  // Get the current depth and score.
  Modelstate *modelstateptr = TOPARRAY(&(srch_state->modelstate), Modelstate);
  Uint curr_depth = modelstateptr->depth;
  double curr_score = modelstateptr->score;

  // Ignore leaf if it is too deep.
  if (curr_depth >= srch_state->maxdepth) return(0);

  // Ignore leaf if minscore not reached.
  if (curr_depth > srch_state->ignore_depth && curr_score < srch_state->minscore) return(0);

  // Ignore leaf if it has no subtree data.
  Uint subtreeindex = srch_state->stree->leafdataindex[leafnumber];
  if (subtreeindex == 0) return(0);

  // Get the subtree data for the branch.
  Subtreedata *nodestdptr = PEEKARRAY(&(srch_state->stree->subtreedata), subtreeindex, Subtreedata);

  // Get the minimum and maximum widths for this node.
  int minw = nodestdptr->min_width;
  int maxw = MIN(nodestdptr->max_width, srch_state->maxdepth);

  // Push matching nodes for each legal width meeting score limit.
  int bg_order = srch_state->multiseq->bg_order;
  for (w=minw; w<=maxw; w++) {
    // Add the column score the current score. There cannot be a non-alphabetic character here.
    Uint aindex = A2I(nodestdptr->head[w-1]); 
    curr_score += srch_state->model->pssm[aindex][w-1];
    if (bg_order > 0) SUBLCBP(nodestdptr->head, w, srch_state->model->alen, srch_state->multiseq->lcbp, bg_order, curr_score);
    if (w >= srch_state->mindepth && w <= srch_state->maxdepth && curr_score >= srch_state->minscore) {
      SETNODESCORE(srch_state, nodestdptr, w, curr_score); 
      PUSHARRAY(&(srch_state->matching_nodes[w]), Subtreedataptr, 128, nodestdptr);
    }
  }

  return(0);
} // sbm_processleaf

// Process a branch1 for finding score-based approximate matches by DFS.
BOOL sbm_processbranch1(Bref nodeptr, void *info)
{
  Searchstate *srch_state = (Searchstate *) info;
  Branchinfo branchinfo;
  int w;

  // Get the current depth and score.
  Modelstate *modelstateptr = TOPARRAY(&(srch_state->modelstate), Modelstate);
  Uint curr_depth = modelstateptr->depth;
  double curr_score = modelstateptr->score;

  // Ignore subtree if node is too deep.
  if (curr_depth >= srch_state->maxdepth) return(False);

  // Ignore subtree if minscore not reached.
  if (curr_depth > srch_state->ignore_depth && curr_score < srch_state->minscore) return(0);

  // Get the head position for the branch.
  getbranchinfostree(srch_state->stree, ACCESSHEADPOS, &branchinfo, nodeptr);
  Uint nodeheadposition = branchinfo.headposition;
  Uchar *head = srch_state->multiseq->sequence + nodeheadposition;

  // Ignore subtree if head starts with a separator.
  if (head[0] == SEPARATOR) return(False);

  // Ignore subtree if it has no subtree data.
  Uint subtreeindex = srch_state->stree->subtreedataindex[nodeheadposition];
  if (subtreeindex == 0) return(False);

  // Get the subtree data for the branch.
  Subtreedata *nodestdptr = PEEKARRAY(&(srch_state->stree->subtreedata), subtreeindex, Subtreedata);

  // Set the minimum and maximum widths for the node.
  int minw = nodestdptr->min_width;
  int maxw = MIN(nodestdptr->max_width, srch_state->maxdepth);

  // Push the score of each word width on to the stack.
  Modelstate modelstate;
  modelstate.depth = curr_depth;
  modelstate.score = curr_score;
  int bg_order = srch_state->multiseq->bg_order;
  // Push one entry on the stack for each position in the branch.
  for (w=minw; w<=maxw; w++) {
    // Update the score for this width. There cannot be a non-alphabetic character here.
    Uint aindex = A2I(nodestdptr->head[w-1]);
    if (aindex == srch_state->model->alen) {
      fprintf(stderr, "Error in sbm_processbranch1: character '%c' is not in the alphabet.\n", nodestdptr->head[w-1]);
      exit(EXIT_FAILURE);
    }
    modelstate.score = (curr_score += srch_state->model->pssm[aindex][w-1]);
    if (bg_order > 0) {
      SUBLCBP(nodestdptr->head, w, srch_state->model->alen, srch_state->multiseq->lcbp, bg_order, curr_score);
      modelstate.score = curr_score;
    }
    modelstate.depth++;
    PUSHARRAY(&(srch_state->modelstate), Modelstate, 128, modelstate);
  }

  return(True);
} // sbm_processbranch1

// Process a branch2 for finding score-based approximate matches by DFS.
Sint sbm_processbranch2(Bref nodeptr, void *info)
{
  Searchstate *srch_state = (Searchstate *) info;
  Branchinfo branchinfo;
  int w;

  // Get the head position for the branch.
  getbranchinfostree(srch_state->stree, ACCESSHEADPOS, &branchinfo, nodeptr);
  Uint nodeheadposition = branchinfo.headposition;
  
  // Get the subtree data for the branch.
  Uint subtreeindex = srch_state->stree->subtreedataindex[nodeheadposition];
  Subtreedata *nodestdptr = PEEKARRAY(&(srch_state->stree->subtreedata), subtreeindex, Subtreedata);

  // Set the minimum and maximum widths for the node.
  int minw = nodestdptr->min_width;
  int maxw = MIN(nodestdptr->max_width, srch_state->maxdepth);

  // Pop the depth/score stack; 1 entry for each position in the node's branch.
  for (w=maxw; w>=minw; w--) {
    Modelstate *modelstate = POPARRAY(&(srch_state->modelstate), Modelstate);
    // Save the node in the matching_nodes[w] array if it matches the search word prefix of width w.
    if (w >= srch_state->mindepth && w <= srch_state->maxdepth && modelstate->score >= srch_state->minscore) {
      SETNODESCORE(srch_state, nodestdptr, w, modelstate->score); 
      PUSHARRAY(&(srch_state->matching_nodes[w]), Subtreedataptr, 128, nodestdptr);
    }
  }

  return(0);
} // sbm_processbranch2

// 
// Compare function for sorting Model pointers
// with qsort in increasing order of test_log_pvalue.
// Return <0 >0
// if the second log_pvalue is <, > than the first.
// Break ties using the train_log_pvalue, then consensus.
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
  } else if (m1->train_log_pvalue > m2->train_log_pvalue) {
    return(+1);
  } else if (m1->train_log_pvalue < m2->train_log_pvalue) {
    return(-1);
  } else {
    return(strncmp(m1->consensus, m2->consensus, MIN(m1->width, m2->width)));
  }

} // compare_model_test_pvalue

// 
// Compare function for sorting Subtree data pointers
// with qsort in increasing order of log_pvalue.
// Return <0 >0
// if the second log_pvalue is <, > than the first.
// Break ties by preferring wider, more positives.
//
static int compare_node_pvalue(
  const void *v1,
  const void *v2
)
{
  const Subtreedataptr *s1 = (const Subtreedataptr *) v1;
  const Subtreedataptr *s2 = (const Subtreedataptr *) v2;

  if ((*s1)->log_pvalue > (*s2)->log_pvalue) {
    return(+1);
  } else if ((*s1)->log_pvalue < (*s2)->log_pvalue) {
    return(-1);
  } else if ((*s1)->max_width < (*s2)->max_width) {
    return(+1);
  } else if ((*s1)->max_width > (*s2)->max_width) {
    return(-1);
  } else if ((*s1)->pos_count < (*s2)->pos_count) {
    return(+1);
  } else if ((*s1)->pos_count > (*s2)->pos_count) {
    return(-1);
  } else {
    return(strncmp((char *) (*s1)->head, (char *) (*s2)->head, (*s1)->max_width));
  }
} // compare_node_pvalue

// 
// Compare function for sorting Subtree data pointers
// with qsort in decreasing order of score.
// Return <0 >0
// if the second score is >, < than the first.
// Break ties by preferring lower p-values.
//
static int compare_node_score(
  const void *v1,
  const void *v2
)
{
  const Subtreedataptr *s1 = (const Subtreedataptr *) v1;
  const Subtreedataptr *s2 = (const Subtreedataptr *) v2;

  if ((*s1)->scores[CURRENT_WIDTH-(*s1)->min_width] < (*s2)->scores[CURRENT_WIDTH-(*s2)->min_width]) {
    return(+1);
  } else if ((*s1)->scores[CURRENT_WIDTH-(*s1)->min_width] > (*s2)->scores[CURRENT_WIDTH-(*s2)->min_width]) {
    return(-1);
  } else {
    return(compare_node_pvalue(s1, s2));
  }
} // compare_node_score

// Free the Sbmdata struct.
void free_sbmdata(
  Sbmdata *sbmdata
) {
  if (sbmdata) {
    free(sbmdata->log_pvalues);
    free(sbmdata->score_thresholds);
    free(sbmdata->pos_counts);
    free(sbmdata->neg_counts);
    if (sbmdata->dtc_sum) free(sbmdata->dtc_sum);
    free(sbmdata);
  }
} // free_sbmdata

//
// Used score-based matching to find all the ZOOPS sites 
// matching a model with a score of 0 or better.
// Determine the optimal score threshold for the sites.
// Compute the p-value of the model at each width.
//
Sbmdata *score_based_matching(
  STREME_OPTIONS_T *options,		// the program options
  Searchstate *srch_state,		// the struct
  BOOL save_sites			// save the sites in the model
) {
  int i, j, k, w;
  OBJFUN_T objfun = options->objfun;
  Suffixtree *stree = srch_state->stree;		// the suffix tree
  Model *model = srch_state->model;
  Uint minwidth = srch_state->mindepth;
  Uint maxwidth = srch_state->maxdepth;
  ArraySubtreedataptr *matching_nodes = srch_state->matching_nodes;
  Multiseq *multiseq = srch_state->multiseq;
  Reference *rootref = srch_state->rootref;	// the root
  Uint npos = multiseq->npos;
  Uint nneg = multiseq->nneg;
  Uint ntot = npos + nneg;
  double avg_poslen = multiseq->avg_poslen;
  double avg_neglen = multiseq->avg_neglen;
  BOOL use_binomial = multiseq->use_binomial;

  // Check for valid call.
  if (save_sites && minwidth != maxwidth) {
    fprintf(stderr, "ERROR: score_based_matching called illegally.\n");
    exit(1);
  }

  // Initialize the arrays of matching nodes for each width for the current model.
  for (w=minwidth; w<=maxwidth; w++) {
    INITARRAY(&(matching_nodes[w]), Subtreedataptr);
  }

  // Search for approximate matches to current model and save
  // them in srch_state.matching_nodes[w].
  depthfirststree(
    stree,
    rootref,
    &sbm_processleaf,
    &sbm_processbranch1,	// processbranch1
    &sbm_processbranch2,	// processbranch2
    NULL,			// stoptraversal
    NULL,			// stopinfo
    (void *)srch_state 		// holds state of search
  );

  // Create the array of sites array and struct to hold return values, indexed by word width.
  Uint overall_maxwidth = options->maxwidth+1;	// maximum depth ever searched (allow for palindrome expansion)
  Site **sites = (Site **) malloc((overall_maxwidth+1) * sizeof(Site *));
  int *nsites = (int *) malloc((overall_maxwidth+1) * sizeof(int));
  Sbmdata *sbmdata = (Sbmdata *) malloc(sizeof(Sbmdata));
  sbmdata->log_pvalues = (double *) malloc((overall_maxwidth+1) * sizeof(double));
  sbmdata->score_thresholds = (double *) malloc((overall_maxwidth+1) * sizeof(double));
  sbmdata->pos_counts = (int *) malloc((overall_maxwidth+1) * sizeof(int));
  sbmdata->neg_counts = (int *) malloc((overall_maxwidth+1) * sizeof(int));
  sbmdata->dtc_sum = (objfun == CD) ? (double *) malloc((overall_maxwidth+1) * sizeof(double)) : NULL;
  for(i=0; i<=overall_maxwidth; i++) {
    nsites[i] = 0;
    sbmdata->log_pvalues[i] = 1;
    sbmdata->score_thresholds[i] = -1;
    sbmdata->pos_counts[i] = sbmdata->neg_counts[i] = -1;
    if (sbmdata->dtc_sum) sbmdata->dtc_sum[i] = -1;
  }

  // Create bit tables to record which sequences are already used.
  Uint *pos_seqno_bittab=NULL, *neg_seqno_bittab=NULL;
  INITBITTAB_TLB(pos_seqno_bittab, npos);
  if (nneg) INITBITTAB_TLB(neg_seqno_bittab, nneg);

  //
  // For each legal width, identify the ZOOPS sites and compute the p-value.
  //
  for (w=minwidth; w<=maxwidth; w++) {
    Uint num_matching_nodes = TOPINDEXARRAY(&(matching_nodes[w]), Subtreedataptr) + 1;
    double pos_sites = npos * MAX(1, (avg_poslen - w + 1));
    double neg_sites = nneg * MAX(1, (avg_neglen - w + 1));
    double bernoulli = use_binomial ? pos_sites / (pos_sites + neg_sites) : -1;

    // Create array for the ZOOPS sites of width w.
    sites[w] = (Site *) malloc((ntot+1) * sizeof(Site));
    nsites[w] = 0;

    // Check that there are words of width w.
    if (num_matching_nodes == 0) {
      model->matches = NULL;
      model->nmatches = 0;
      free(sites[w]);
      continue;
    }

    // Clear the bit tables.
    CLEARBITTAB(pos_seqno_bittab, npos);
    if (nneg) CLEARBITTAB(neg_seqno_bittab, nneg);

    // Sort the matching word nodes by their scores.
    CURRENT_WIDTH = w;		// Nasty global.  Not thread safe.
    SORTARRAY(&(matching_nodes[w]), Subtreedataptr, compare_node_score);

    // Limit processing to one node per positive or negative sequence.
    // This prevents excessive run time when the motif is very non-specific
    // and has many times more matching nodes than there are sequences.
    if (num_matching_nodes > ntot) num_matching_nodes = ntot;

    //
    // Identify ZOOPS sites for the motif of this width.
    //
    int pos_matches=0, neg_matches=0, new_pos_matches=0, new_neg_matches=0;
    int added_pos_matches=0, added_neg_matches=0;
    // Loop over the matching nodes by decreasing score.
    for (i=0; i<num_matching_nodes; i++) {
      Subtreedataptr mnptr = *(PEEKARRAY(&(matching_nodes[w]), i, Subtreedataptr));
      Uchar *word = mnptr->head;

      // Get the counts of (new) sequences that are in this node
      // and update the used sequences.
      if (mnptr->pos_count) {
	if (mnptr->pos_seqno_bittab) {
	  ORBITTABS(pos_seqno_bittab, mnptr->pos_seqno_bittab, npos);
	} else {
	  for (j=0; j<mnptr->pos_count; j++) {
	    SETIBIT(pos_seqno_bittab, mnptr->pos_seqno_list[j]);
	  }
	}
      } // pos_count
      if (mnptr->neg_count) {
	if (mnptr->neg_seqno_bittab) {
	  ORBITTABS(neg_seqno_bittab, mnptr->neg_seqno_bittab, nneg);
	} else {
	  for (j=0; j<mnptr->neg_count; j++) {
	    SETIBIT(neg_seqno_bittab, mnptr->neg_seqno_list[j]);
	  }
	}
      } // neg_count
      COUNTBITTAB(new_pos_matches, pos_seqno_bittab, npos);
      added_pos_matches = new_pos_matches - pos_matches;
      if (nneg) {
        COUNTBITTAB(new_neg_matches, neg_seqno_bittab, nneg);
	added_neg_matches = new_neg_matches - neg_matches;
      }
      pos_matches += added_pos_matches;
      neg_matches += added_neg_matches;

      // Skip this node if it is not a ZOOPS site.
      if (added_pos_matches == 0 && added_neg_matches == 0) {
        continue;
      }

      // Add this site to the list for this width.
      Site *site = &(sites[w][nsites[w]++]);
      site->seqno = mnptr->seqno;
      site->pos = mnptr->seqpos;
      site->strand = '+';		  		// strand is always positive
      site->is_positive = (added_pos_matches > 0);	// at least one positive ZOOPS site
      site->score = mnptr->scores[w-mnptr->min_width];
      site->pos_count = mnptr->pos_count;
      site->neg_count = mnptr->neg_count;
      site->pos_zoops = added_pos_matches;
      site->neg_zoops = added_neg_matches;
      site->log_pvalue = mnptr->log_pvalue;
      site->dtc = (objfun == CD && site->pos_count > 0) ? mnptr->dtc_sum / site->pos_count : -1;
    } // matching node

    // Sort the sites by score, placing tied sites that have neg_zoops > 0 first.
    qsort(sites[w], nsites[w], sizeof(Site), compare_site_score);

    //
    // Get the optimum score threshold and compute the p-value for this width.
    //
    double log_pvalue, best_log_pvalue=0, best_score_threshold=0;
    double dtc_sum=0, best_dtc_sum=-1;
    int best_pos_matches=0, best_neg_matches=0;
    pos_matches = neg_matches = 0;
    double avg_pos_length = multiseq->avg_poslen;
    for (i=0; i<nsites[w]; i++) {
      Site *site = &(sites[w][i]);
      // Update the counts of matches with the number of sequences for which this is a ZOOPS site.
      pos_matches += site->pos_zoops;
      neg_matches += site->neg_zoops;
      // Compute the p-value and save the state if it is the best so far.
      log_pvalue = 0;
      if (objfun == CD) {
	// Add a the average dtc of the node multiplied by the number of ZOOPS sites
	// to the overall sum for the motif.
	if (site->is_positive) dtc_sum += site->pos_zoops * site->dtc;
	GET_PVALUE(log_pvalue, objfun, use_binomial, pos_matches, neg_matches, npos, nneg, bernoulli, dtc_sum, w, avg_poslen, 
	  "score_based_matching", options->pv_cache, options->cache_length);
      } else if (objfun == DE) {
	// DE:
        // Only compute the p-value when at the end of a run of
        // sites with the same score.
        // Since sites with positive sequences come last in runs,
        // p-value can't improve if last site in run is not positive.
	if ( site->is_positive &&
          (
	    i==nsites[w]-1 || // last site
	    fabs(site->score - (site+1)->score) > FLOAT_EPS
          )
	) {
	  GET_PVALUE(log_pvalue, objfun, use_binomial, pos_matches, neg_matches, npos, nneg, bernoulli, 0, 0, 0, 
	    "score_based_matching", options->pv_cache, options->cache_length);
	}
      } // DE

      // Update if the pvalue improved.
      if (log_pvalue < best_log_pvalue) {
        best_score_threshold = site->score;
        best_log_pvalue = log_pvalue;
        best_pos_matches = pos_matches;
        best_neg_matches = neg_matches;
        best_dtc_sum = dtc_sum;
      }

    } // nsites

    // Set the return values in the model.
    model->score_threshold = best_score_threshold - FLOAT_EPS;	// in case of rounding error later
    model->train_log_pvalue = best_log_pvalue;
    model->train_pos_count = best_pos_matches;
    model->train_neg_count = best_neg_matches;
    model->train_ratio = (objfun == DE) ? (model->train_pos_count+1.0)/(model->train_neg_count+1) : -1;
    model->train_dtc = (objfun == CD && best_pos_matches > 0) ? best_dtc_sum/best_pos_matches : -1;
    if (save_sites) {
      // Return the sites.
      model->matches = sites[w];
      model->nmatches = nsites[w];
    } else {
      // Free the sites of width w if we are not returning them.
      free(sites[w]);
      model->matches = NULL;
    }

    // Update the return data.
    sbmdata->log_pvalues[w] = best_log_pvalue;
    sbmdata->score_thresholds[w] = best_score_threshold;
    sbmdata->pos_counts[w] = best_pos_matches;
    sbmdata->neg_counts[w] = best_neg_matches;
    if (sbmdata->dtc_sum) sbmdata->dtc_sum[w] = best_dtc_sum;
  } // w

  // Free space.
  FREEBITTAB(pos_seqno_bittab);
  if (nneg) FREEBITTAB(neg_seqno_bittab);
  for (w=minwidth; w<=maxwidth; w++) FREEARRAY_TLB(&(matching_nodes[w]), Subtreedataptr);
  free(sites);
  free(nsites);

  return(sbmdata);
} // score_based_matching

// For each node up to maxdepth, add counts of positive and negative sequences below it,
// and the objective function p-value.  Return the valid nodes, sorted by log_pvalue.
// A node is valid if:
//   1) it is not the root;
//   2) it is not deeper than maxdepth
//   3) no separators up to start of branch
//   4) it is at least as deep as mindepth
//   Note: Some words in the width range may not have a valid node (they are at leaves).
//   Note: The count for such words will be 1.
//
ArraySubtreedataptr *add_counts_to_stree(
  STREME_OPTIONS_T *options,		// streme options
  Suffixtree *stree,			// the suffix tree
  Reference *rootref,			// the root
  Multiseq *multiseq 			// the positive and negative sequences
) {
  int i, cnt;
  Uint npos = multiseq->npos;
  Uint nneg = multiseq->nneg;

  DEBUG_MSG(NORMAL_VERBOSE, "# Adding sequence counts and p-values to nodes using depth first search...\n");

  // Do a depth-first search of the tree.
  Countstate cnt_state;
  cnt_state.stree = stree;
  cnt_state.multiseq = multiseq;
  cnt_state.mindepth = MIN_SEED_WIDTH;			// the minimum length of words to add counts for
  cnt_state.maxdepth = options->maxwidth+1;		// the maximum length of words to add counts for (+1 for palindromes)
  INITARRAY(&(cnt_state.valid_node_indices), Uint);
  cnt_state.objfun = options->objfun;			// objective function
  cnt_state.pv_cache = options->pv_cache;		// for efficiency
  cnt_state.cache_length = options->cache_length;

  // Initialize the stack that contains the mother of the current node.
  INITARRAY(&(cnt_state.mother), Bref);
  PUSHARRAY(&(cnt_state.mother), Bref, 128, rootref->address);

  // Initialize the stack that contains the closest ancestor that is not too deep.
  INITARRAY(&(cnt_state.ancestor), Bref);
  PUSHARRAY(&(cnt_state.ancestor), Bref, 128, rootref->address);

  // Initialize the numbers of sequences in subtree, 
  // and the bit table for the root.
  Subtreedata std;
  std.head = NULL;
  std.min_width = 0;
  std.max_width = 0;
  std.pos_count = 0;
  std.neg_count = 0;
  INITBITTAB_TLB(std.pos_seqno_bittab, npos);
  std.neg_seqno_bittab = NULL;
  if (nneg) INITBITTAB_TLB(std.neg_seqno_bittab, nneg);
  std.pos_seqno_list = NULL;
  std.neg_seqno_list = NULL;
  std.scores = NULL;
  // Push the subtree data for the root onto the tree's subtree data stack.
  INITARRAY(&(stree->subtreedata), Subtreedata);
  PUSHARRAY(&(stree->subtreedata), Subtreedata, 128, std);

  // Set up the index to subtree data; initialized to point to the root.
  stree->subtreedataindex = ALLOCSPACE_TLB(NULL,Uint,stree->nextfreeleafnum+1);
  stree->leafdataindex = ALLOCSPACE_TLB(NULL,Uint,stree->nextfreeleafnum+1);
  for (i=0; i<stree->nextfreeleafnum+1; i++) {
    stree->subtreedataindex[i] = stree->leafdataindex[i] = 0; // root
  }

  // Add the counts and p-values down to depth maxdepth.
  depthfirststree(
    stree,
    rootref,
    &cnt_processleaf,
    &cnt_processbranch1,// processbranch1
    &cnt_processbranch2,// processbranch2
    NULL,		// stoptraversal
    NULL,		// stopinfo
    (void *)&cnt_state 	// holds state of search
   );

  // Get the pointers to the valid nodes.
  Uint num_valid_nodes = TOPINDEXARRAY(&(cnt_state.valid_node_indices), Uint) + 1;
  ArraySubtreedataptr *valid_nodes = (ArraySubtreedataptr *) malloc(sizeof(ArraySubtreedataptr));
  INITARRAY(valid_nodes, Subtreedataptr);
  for (cnt = 0; cnt < num_valid_nodes; cnt++) {
    Uint subtreeindex = *(PEEKARRAY(&(cnt_state.valid_node_indices), cnt, Uint));
    Subtreedata *vnptr = PEEKARRAY(&(stree->subtreedata), subtreeindex, Subtreedata);
    PUSHARRAY(valid_nodes, Subtreedataptr, 128, vnptr);
  }
  PRINTTIME("ADDCOUNTS");

  // Sort the valid nodes by p-value.
  DEBUG_FMT(NORMAL_VERBOSE, "# Sorting %d initial seed nodes by p-value...\n", num_valid_nodes);
  SORTARRAY(valid_nodes, Subtreedataptr, compare_node_pvalue);
  PRINTTIME("SORTSEEDS");

  // free the count state
  FREEARRAY_TLB(&(cnt_state.mother), Bref);
  FREEARRAY_TLB(&(cnt_state.ancestor), Bref);
  FREEARRAY_TLB(&(cnt_state.valid_node_indices), Uint);

  return(valid_nodes);
} // add_counts_to_stree

//
// Compare two evaluated seed objects in increasing wgt_log_pvalue order.
// This function returns
//   1 if s1->wgt_log_pvalue > s2->wgt_log_pvalue,
//   -1 if s1->wgt_log_pvalue < s2->wgt_log_pvalue,
// Break ties by preferring wider motifs, wider initial widths.
//
int compare_evaluated_seeds(
  void *v1,
  void *v2
) {
  const Evaluatedseed *s1 = (const Evaluatedseed *) v1;
  const Evaluatedseed *s2 = (const Evaluatedseed *) v2;
  if (s1->log_pvalue > s2->log_pvalue) {
    return(+1);
  } else if (s1->log_pvalue < s2->log_pvalue) {
    return(-1);
  } else if (s1->width < s2->width) {
    return(+1);
  } else if (s1->width > s2->width) {
    return(-1);
  } else if (s1->seed->max_width < s2->seed->max_width) {
    return(+1);
  } else if (s1->seed->max_width > s2->seed->max_width) {
    return(-1);
  } else {
    return(strncmp((char *) s2->seed->head, (char *) s1->seed->head, s1->width));
  }
} // compare_evaluated_seeds

//
// Create a model that is the reverse complement of the given model.
// Only updates the PSPM (probs) matrix.
//
Model *create_model_rc(
  Model *model
) {
  int i, j;
  Model *new_model = (Model *) malloc(sizeof(Model));
  *new_model = *model;
  int alen = model->alen;
  int w = model->width;
  new_model->nmatches = 0;
  new_model->matches = NULL;
  new_model->consensus = NULL;

  for (i=0; i<alen; i++) {
    for (j=0; j<w; j++) {
      new_model->probs[I2CI(i)][w-j-1] = model->probs[i][j];
    }
  }

  return(new_model);
} // create_model_rc

//
// Create a shifted model.
//
Model *create_shifted_model(
  int shift,			// amount to shift the model
  Model *model,			// the model to shift
  Multiseq *multiseq 		// the sequences
) {
  int alen = model->alen;
  int w = model->width;
  int i, j;
  double *background = multiseq->background;	// 0-order background distribution

  // Create the shifted model.
  Model *shifted_model = (Model *) malloc(sizeof(Model));
  *shifted_model = *model;
  shifted_model->width = shift > 0 ? w + 2*shift : w + shift;

  // Set shifted model columns to 0-order background.
  for (i=0; i<alen; i++) for (j=0; j<shifted_model->width; j++) shifted_model->probs[i][j] = background[i];

  // Copy the (shifted) PSPM.
  for (i=0; i<alen; i++) {
    if (shift >= 0) {
      // Shift right.
      for (j=0; j<w; j++) shifted_model->probs[i][j+shift] = model->probs[i][j];
    } else {
      // Shift left.
      for (j=-shift; j<w; j++) shifted_model->probs[i][j+shift] = model->probs[i][j];
    }
  } 

  // Initialize the PSSM.
  INIT_PSSM_FROM_PROBS(shifted_model, multiseq->background, multiseq->bg_order);

  return(shifted_model);
} // create_shifted_model

//
// Create a palindromic model from a given model by aligning it
// with its reverse-complement.
//
// The best alignment of the model PSPM to its reverse complement is found.
// Best is based on the sum of the Euclidean distance between aligned columns, and
// between the background distribution and unaligned columns.
// The palindromic model is the average of the model PSPM and the (shifted) reverse
// complement model.
//
Model *create_palindromic_model(
  Model *model,			// the model to turn into a palindrome
  Multiseq *multiseq,		// the sequences
  BOOL trim			// trim result to original width
) {
  int alen = model->alen;
  int w = model->width;
  int i, j, dir, offset, shift;
  double *background = multiseq->background;	// 0-order background distribution

  // Create the reverse complement of the model.
  Model *rc_model = create_model_rc(model);

  if (verbosity >= HIGHER_VERBOSE) {
    PRINT_MODEL(model, "forward", 1)
    PRINT_MODEL(rc_model, "rc", 1)
  }

  // Find the best alignment of the RC to the original model.
  double best_ed = -1;
  int best_offset = 0;
  int best_dir = 0;
  for (dir=0; dir<2; dir++) {
    Model *mleft, *mright;
    if (dir == 0) {
      mleft = rc_model;
      mright = model;
    } else {
      mleft = model;
      mright = rc_model;
    }
    for (offset = dir ? 1 : 0; offset < w/2; offset++) {
      int pw = w + offset;
      double ed = 0;
      for (j=0; j<pw; j++) {
        double sumsq = 0;
        double diff;
        for (i=0; i<alen; i++) {
	  if (j < offset) {
	    diff = background[i] - mleft->probs[i][j];
	  } else if (j >= w) {
	    diff = background[i] - mright->probs[i][j-offset];
	  } else {
	    diff = mleft->probs[i][j] - mright->probs[i][j-offset];
	  }
          sumsq += diff * diff;
        } // letter
        ed += sqrt(sumsq);
      } // position
      if (best_ed < 0 || ed < best_ed) {
        best_ed = ed;
        best_offset = offset;
        best_dir = dir;
      }
    } // offset
  } // dir

  // Create the palindromic model.
  Model *pal_model = (Model *) malloc(sizeof(Model));
  *pal_model = *model;
  pal_model->nmatches = 0;
  pal_model->matches = NULL;		// Don't clobber this.
  pal_model->consensus = NULL;		// Don't clobber this.
  pal_model->train_log_pvalue = 0;	// Signal model not evaluated yet.
  // Set width of palindromic model.
  pal_model->width = w + best_offset;	
  // Zero palindromic model.
  for (i=0; i<alen; i++) for (j=0; j<pal_model->width; j++) pal_model->probs[i][j] = 0;
  // Add the original model probs (shifted) to the palindromic model.
  shift = (best_dir == 0) ? best_offset : 0;
  for (i=0; i<alen; i++) for (j=0; j<w; j++) pal_model->probs[i][j+shift] += model->probs[i][j];
  // Add the rc_model probs (shifted) to the palindromic model.
  shift = (best_dir == 0) ? 0 : best_offset;
  for (i=0; i<alen; i++) for (j=0; j<w; j++) pal_model->probs[i][j+shift] += rc_model->probs[i][j];
  // Compute the average probs.
  for (i=0; i<alen; i++) for (j=best_offset; j<w; j++) pal_model->probs[i][j] /= 2;

  // Trim the model to the original width if requested (old w+1 if offset is odd).
  if (trim && best_offset!=0) {
    // Compute trimmed width.
    int new_w = (best_offset % 2) ? w+1 : w;
    pal_model->width = new_w;
    // Shift model left by offset.
    for (j=0; j<new_w; j++) for (i=0; i<alen; i++) pal_model->probs[i][j] = pal_model->probs[i][j+best_offset/2];
  }

  if (verbosity >= HIGHER_VERBOSE) {
    PRINT_MODEL(pal_model, "pal", 1)
  }

  // Initialize the PSSM.
  INIT_PSSM_FROM_PROBS(pal_model, background, multiseq->bg_order);
  pal_model->is_palindromic = True;
  pal_model->ed = best_ed;

  // Free the RC model.
  free(rc_model);

  return(pal_model);
} // create_palindromic_model 

//
// Evaluate the initial seeds using score-based matching to find the best seed of each width.
// Updates an array containing the evaluated seeds.
//
void evaluate_initial_seeds(
  STREME_OPTIONS_T *options,		// the program options
  ArraySubtreedataptr *initial_seeds,	// the seeds to evaluate using score-based matching
  Searchstate *srch_state, 		// state variable for score-based matching
  int motifno,				// current motif number
  ArrayEvaluatedseed *evaluated_seeds	// IN/OUT; initialized here
) {
  int w, w0, cnt;
  Suffixtree *stree = srch_state->stree;		// the suffix tree
  Reference *rootref = srch_state->rootref; 		// the root
  Multiseq *multiseq = srch_state->multiseq;		// the training positive and negative sequences
  Uint minwidth = options->minwidth;			// the minimum length of words to add counts for
  Uint maxwidth = options->maxwidth;			// the maximum length of words to add counts for
  Model *model = srch_state->model;			// the model
  Uint npos = multiseq->npos;
  Sbmdata *sbmdata = NULL;				// initialized in score_based_matching()

  // Initialize the evaluated seeds array.
  INITARRAY(evaluated_seeds, Evaluatedseed);

  // Set the number of initial seeds remaining to evaluate for each width.
  // Note that seeds of widths as small as 1 are evaluated since we extend
  // them during the search.
  int n_done = 0;
  int n_initial_seeds[maxwidth+2];
  for (w = 1; w <= maxwidth+1; w++) n_initial_seeds[w] = options->neval;

  //
  // Loop over the initial seeds.
  // Use score-based approximate matching to compute the p-value of the 
  // each of the top N initial seeds of each width.  This allows a motif
  // with a strong signal in its left half to be detected here even though
  // the full-width seed gets a poor initial seed score.
  //
  Uint num_initial_seeds = TOPINDEXARRAY(initial_seeds, Subtreedataptr) + 1;
  DEBUG_FMT(NORMAL_VERBOSE, "#   Evaluating the top %d seeds (out of %d) of each width in range [%d, %d] ...\n", 
    options->neval, num_initial_seeds, MIN_SEED_WIDTH, maxwidth);
  int nseeds = options->neval * (maxwidth - MIN_SEED_WIDTH + 1);

  //
  // Evaluate the initial seeds using score-based matching.
  // Only neval seeds of each seed->max_width are evaluated.
  //
  for (cnt = 0; cnt < num_initial_seeds && n_done <= maxwidth-MIN_SEED_WIDTH; cnt++) {
    Subtreedata *seed = *(PEEKARRAY(initial_seeds, cnt, Subtreedataptr));

    // Skip the seeds (valid nodes from add_counts_to_stree) that are one wider
    // than the maximum width (added for palindromic +1 expansion).
    if (seed->max_width > maxwidth) continue;

    Uchar *seed_word = seed->head;
    if (verbosity >= HIGHER_VERBOSE) {
      w = seed->max_width;
      fprintf(stderr, "# Initial seed %d: %*.*s log_pvalue %.2f w0 %d n_done %d", 
	cnt, w,w,seed_word, seed->log_pvalue, w, n_done);
      int i;
      for (i=1; i<=maxwidth; i++) fprintf(stderr, " %2d", n_initial_seeds[i]);
      fprintf(stderr, "\n");
    }

    // Skip seed if we have evaluated this max_width enough.
    if (n_initial_seeds[seed->max_width] == 0) continue;

    // Process only lexicographically smaller seed if doing reverse complements
    // since it doesn't matter which of the two we evaluate--the process is symmetrical.
    if (multiseq->do_rc && strncmp((char *)seed->head, (char *)seed->rc_head, seed->max_width) > 0) continue; 
    
    // Try to extend the head to maxwidth or until the first separator.
    for (w0=seed->max_width; 
      w0<=maxwidth && 
      seed_word + w0 - multiseq->sequence <= multiseq->totallength && 
      seed_word[w0-1] != SEPARATOR; 
      w0++) /*empty*/;
    w0--;

    // Skip seed if it can't be extended to at least minwidth.
    if (w0 < minwidth) continue;

    // Only evaluate the neval initial seeds of each initial maximum width.
    if (--n_initial_seeds[seed->max_width] == 0) n_done++;

    // Finish initializing the search state.
    srch_state->mindepth = seed->min_width;	// minimum width to match; seeds can be smaller than options->minwidth
    srch_state->maxdepth = w0;			// maximum width to match

    // Initialize the MODEL to 0/1 consensus of SEED.
    model->width = w0;
    INIT_PSSM_TO_CONSENSUS(model, seed_word, SBM_SEED_PRIOR, multiseq->background, multiseq->bg_order);
    
    // Get the approximate matches and p-value for this initial seed.
    sbmdata = score_based_matching(options, srch_state, False);

    // Add the evaluated seeds to the array.
    for (w = minwidth; w <= srch_state->maxdepth; w++) {
      Evaluatedseed evaluated_seed;
      evaluated_seed.seed = seed;
      evaluated_seed.width = w;
      evaluated_seed.log_pvalue = sbmdata->log_pvalues[w];
      PUSHARRAY(evaluated_seeds, Evaluatedseed, 128, evaluated_seed);
    } // w

    // free space
    free_sbmdata(sbmdata);
  } // initial seeds

  DEBUG_FMT(NORMAL_VERBOSE, "%s", "\n");
  PRINTTIME("EVALUATESEEDS");
} // evaluate_initial_seeds

//
// Optimize the enrichment score of a model by testing models created from 
// nested subsets of its best sites.
//   1) Find the ZOOPS sites and enrichment p-value for the input model.
//   2) Until enrichment p-value does not improve:
//   3)   For new model created from the sites with score < Best-inc, Best-2*inc, Best-3*inc, ...:
//   4)     Find the ZOOPS sites and enrichment p-value using the new model.
//   5)     Save the new model if the p-value is better than any model's so far.
//   6)   Replace the input model with the best nested model, if any is better. 
//   7) Return model.
//
void nested_subset_enrichment(
  STREME_OPTIONS_T *options,		// STREME options
  Searchstate *srch_state,		// state variable for score-based matching
  int n_subsets				// number of nested subsets to use; 1 means use model->score_threshold
) {
  int i, j, i_iter, i_subset, i_match;
  int niter = options->niter;
  Suffixtree *stree = srch_state->stree;	// the suffix tree
  Model *model = srch_state->model;		// the motif model
  Multiseq *multiseq = srch_state->multiseq;	// the positive and negative sequences
  int alen = model->alen;
  int w = model->width;
  double *background = multiseq->background;	// 0-order background
  double m1, e1, prec=1;
  int best_i_subset = 0;

  // Check for valid call.
  if (w != srch_state->mindepth || w != srch_state->maxdepth) {
    fprintf(stderr, "ERROR: nested_subset_enrichment called illegally.\n");
    fprintf(stderr, "w %d mindepth %d maxdepth %d\n", w, srch_state->mindepth, srch_state->maxdepth);
    exit(1);
  }

  //
  // Loop until no better model is found.
  //
  for (i_iter=0; i_iter < niter; i_iter++) {
    // Get the ZOOPS sites and p-value for the initial model.
    Sbmdata* sbmdata = score_based_matching(options, srch_state, True);
    int nmatches = model->nmatches;	// number matching nodes (score > minscore) found for current model
    if (nmatches == 0) {
      break;
    }
    free_sbmdata(sbmdata);
    Site *matches = model->matches;			// the model's matches, sorted by score
    double best_score = matches[0].score;		// score of best match (could be negative sequence)
    double inc = best_score / (n_subsets+1.0);
    double best_threshold = model->score_threshold;

    //
    // Test models created from nested subsets of the model's sites.
    //
    // Initialize the nested model from the current model.
    Model nested_model = *model;
    if (verbosity >= HIGH_VERBOSE) {
      exp10_logx(nested_model.train_log_pvalue/log(10), m1, e1, prec);
      fprintf(stderr, "# Initial model pvalue %3.1fe%+04.0f npos %d nneg %d score_thresh %f\n",
	m1,e1, nested_model.train_pos_count, nested_model.train_neg_count, nested_model.score_threshold);
      PRINT_MODEL(&nested_model, "Initial", 1);
    }
    // Set the counts to the background counts.
    for (i=0; i<alen; i++) for (j=0; j<w; j++) nested_model.counts[i][j] = NESTED_COUNT_PRIOR * background[i];
    double total_wgt_pos_count = 0;
    int total_pos_count=0;
    Model best_nested_model;				// the best nested model found so far
    best_nested_model.train_log_pvalue = 1;		// signals invalid model
    i_match = 0;
    for (i_subset=1; i_subset<=n_subsets; i_subset++) {
      double threshold = n_subsets == 1 ? model->score_threshold : best_score - (i_subset * inc);
      if (n_subsets > 1) DEBUG_FMT(HIGH_VERBOSE, "# NESTED THRESHOLD %d: %f\n", i_subset, threshold);
      //
      // Add the next (positive) sites to the nested model's count matrix.
      //
      BOOL found_positive_match = False;
      for ( ; i_match<nmatches; i_match++) {
        Site *match = &(matches[i_match]);
        if (match->score < threshold) break;
        total_pos_count += match->pos_zoops; 	// The total positive ZOOPS count for seeing if we've got all sites.
        // Count to add: minus the node's log_pvalue scaled by the fraction of its sites that are ZOOPS sites.
        double pos_scale = match->pos_count > 0 ? (double) match->pos_zoops / match->pos_count : 0;
        double wgt_pos_count = -(match->log_pvalue) * pos_scale;
        total_wgt_pos_count += wgt_pos_count;	// Total weighted count for normalizing the counts matrix.
        if (match->is_positive && wgt_pos_count > 0) {
          found_positive_match = True;
          Uint pos = match->pos;
	  Uint seqstart = multiseq->seqstarts[match->seqno];
	  Uchar *word = multiseq->sequence+seqstart+pos;
	  for (j=0; j<w; j++) {
            if (word[j] == SEPARATOR) {
              fprintf(stderr, "ERROR: match contains the SEPARATOR (%*.*s).\n", w, w, word);
              exit(EXIT_FAILURE);
            }
	    Uint aindex = match->strand != '-' ? A2I(word[j]) : C2I(word[w-j-1]);
	    nested_model.counts[aindex][j] += wgt_pos_count;
	  } // col
        }
      } // i_match
   
      // Continue if not a single new positive site added to nested_model.
      if (!found_positive_match) continue;

      // Convert counts to probabilities.
      for (i=0; i<alen; i++) for (j=0; j<w; j++) 
	nested_model.probs[i][j] = nested_model.counts[i][j]/(total_wgt_pos_count + NESTED_COUNT_PRIOR);
      // Convert the nested model PROBS into a palindrome if requested.
      if (nested_model.is_palindromic) PALINDROMIZE(&nested_model);
      // Convert probabilities to PSSM.
      INIT_PSSM_FROM_PROBS(&nested_model, background, multiseq->bg_order);

      //
      // Score the sequences with the nested model and get its enrichment p-value.
      //
      srch_state->model = &nested_model;		// Nested model pointer passed via srch_state.
      Sbmdata* sbmdata = score_based_matching(options, srch_state, True);
      free_sbmdata(sbmdata);
      srch_state->model = model;			// Restore model pointer in srch_state.
      // Show result for debugging.
      if (verbosity >= HIGH_VERBOSE) {
        exp10_logx(nested_model.train_log_pvalue/log(10), m1, e1, prec); 
        fprintf(stderr, "# Nested model %d pvalue %3.1fe%+04.0f npos %d nneg %d score_thresh %f\n", 
          i_subset, m1,e1, nested_model.train_pos_count, nested_model.train_neg_count, nested_model.score_threshold);
        PRINT_MODEL(&nested_model, "Nested", 1);
      }
 
      // Free the nested model matches.
      free(nested_model.matches);
      nested_model.nmatches = 0;
      // Save the nested model if it is better than the model and the best nested model.
      if (nested_model.train_log_pvalue < model->train_log_pvalue && 
	nested_model.train_log_pvalue < best_nested_model.train_log_pvalue) {
	best_nested_model = nested_model;
        best_i_subset = i_subset;
        best_threshold = threshold;
      } // found improved model
      // Break if all ZOOPS sites already included.
      if (total_pos_count >= multiseq->npos) break;
    } // i_subset

    // Free the original model's sites.
    free(model->matches);
    model->matches = NULL;

    // Replace contents of model with contents of best nested model if it is better.
    if (best_nested_model.train_log_pvalue < model->train_log_pvalue) {
      *model = best_nested_model;
      exp10_logx(model->train_log_pvalue/log(10), m1, e1, prec);
      if (n_subsets == 1) {
	DEBUG_FMT(NORMAL_VERBOSE, "#     ITER %d pvalue %3.1fe%+04.0f pos %d neg %d score_threshold %.3f\n", 
	  i_iter+1, m1,e1, model->train_pos_count, model->train_neg_count, model->score_threshold);
      } else {
	DEBUG_FMT(NORMAL_VERBOSE, "#     ITER %d pvalue %3.1fe%+04.0f pos %d neg %d best_score %f best_threshold %.3f inc %f best_i_subset %d / %d\n", 
	  i_iter+1, m1,e1, model->train_pos_count, model->train_neg_count, best_score, best_threshold, inc, best_i_subset, n_subsets);
      }
    } else {
      // Did not find a better model.
      break;
    }
  } // i_iter

  // Set the model consensus sequence.
  model->consensus = get_single_letter_consensus(model);

  // Score the model on the hold-out sequences.
  if (srch_state->test_multiseq) {
    // Score the hold out sequences using the score_threshold in the model
    // to get an unbiased p-value.
    score_model_pssm(options, srch_state->test_multiseq, model, False, True, NONE, True, False);
  } else {
    // Set the p-value if there is no hold-out set.
    model->test_log_pvalue = DEFAULT_LOG_PVALUE;
  }

} // nested_subset_enrichment

// 
// Refine the final nref seeds of each width using score-based approximate matching.  
// Return the best overall model.
//
Model *refine_seeds(
  STREME_OPTIONS_T *options,		// the program options
  int motifno,				// number of current motif
  Searchstate *srch_state,		// state variable for score-based matching
  ArrayEvaluatedseed *evaluated_seeds	// array containing the evaluated initial seeds
) {
  int iter, w, iseed;
  Suffixtree *stree = srch_state->stree;		// the suffix tree
  Multiseq *multiseq = srch_state->multiseq;		// the training positive and negative sequences
  Multiseq *test_multiseq = srch_state->test_multiseq;	// the hold-out positive and negative sequences
  BOOL use_binomial = multiseq->use_binomial;		// use the binomial distribution
  Model *model = srch_state->model;			// the model
  Uint minwidth = options->minwidth;			// the minimum length of words to add counts for
  Uint maxwidth = options->maxwidth;			// the maximum length of words to add counts for
  int num_widths = maxwidth - minwidth + 1;		// number of widths being tried
  Uint npos = multiseq->npos;				// number of positive training sequences
  double best_log_pvalue = 1;				// enrichment p-value
  double best_ratio = 0;				// enrichment ratio
  double log_thresh = log(options->thresh);		// log of significance threshold
  Model *final_model;					// model to return
  final_model = (Model*) calloc(1, sizeof(Model));	// 0-filled
  final_model->test_log_pvalue = 1;			// not a reportable model yet
  final_model->width = minwidth;			// prefer wider models with same log_pvalue
  int nref = MAX(1, options->nref);			// in case options->nref 0 && == 0

  // Get the number of candidate seeds for refinement.
  int num_evaluated_seeds = TOPINDEXARRAY(evaluated_seeds, Evaluatedseed) + 1;
  if (options->nref > 0) {
    DEBUG_FMT(NORMAL_VERBOSE, "# Refining the %d best evaluated seeds of each width in range [%d, %d]\n"
      "#   using score-based refinement (%d evaluated seeds)...\n", 
      options->nref, minwidth, maxwidth, num_evaluated_seeds);
   } else {
    DEBUG_FMT(NORMAL_VERBOSE, "# Refining the single best evaluated seed with width in range [%d, %d]\n"
      "#   using score-based refinement (%d evaluated seeds)...\n", 
      minwidth, maxwidth, num_evaluated_seeds);
   }

  // Sort the evaluated initial seeds.
  SORTARRAY(evaluated_seeds, Evaluatedseed, compare_evaluated_seeds);
  // Print the initial seeds and their evaluation scores for debugging.
  if (verbosity >= HIGHER_VERBOSE) {
    for (iseed=0; iseed<num_evaluated_seeds; iseed++) {
      Evaluatedseed *evaluated_seed = PEEKARRAY(evaluated_seeds, iseed, Evaluatedseed);
      Subtreedata *seed = evaluated_seed->seed;
      w = evaluated_seed->width;
      DEBUG_FMT(NORMAL_VERBOSE, "# Evaluated seed %d: %*.*s log_pvalue %.2f width %d w0 %d\n", 
	iseed, w,w,seed->head, evaluated_seed->log_pvalue, evaluated_seed->width, evaluated_seed->seed->max_width);
    }
  }

  // Create a hash table to keep track of which words we have already evaluated.
  HASH_TABLE ht_words = hash_create(num_widths * nref, free);

  // Set the number of evaluated seeds remaining to refine.
  int n_evaluated_seeds[maxwidth+1];
  for (w = 1; w <= maxwidth; w++) {
    n_evaluated_seeds[w] = (w < minwidth) ? 0 : nref;
  }

  //
  // Refine nref of the best evaluated initial seeds of each width.
  //
  int num_refined_seeds, candno=0;
  int max_refined_seeds = num_widths * nref;
  BOOL first_motif = True;
  for (iseed=num_refined_seeds=0; iseed<num_evaluated_seeds && num_refined_seeds <= max_refined_seeds; iseed++) {

    // Get the next evaluated seed.
    Evaluatedseed *evaluated_seed = PEEKARRAY(evaluated_seeds, iseed, Evaluatedseed);
    Subtreedata *seed = evaluated_seed->seed;
    w = evaluated_seed->width;

    // See if the word for this seed has already been refined.
    Uchar *search_word = seed->head;
    char *key;
    int dummy = asprintf(&key, "%*.*s", w, w, search_word);
    HASH_TABLE_ENTRY *hte = hash_lookup_str(key, ht_words);
    if (hte) {			// search word already in hash table
      free(key);
      continue;
    } else {			// search word not in hash table
      // Only evaluate nref evaluated seeds of each width.
      if (--n_evaluated_seeds[w] < 0) {
        free(key);
        continue;
      }
      num_refined_seeds++;
      hash_insert_str(key, ht_words);
      free(key);
    }

    // Fix the model width.
    model->width = w;
    model->initial_width = seed->max_width;
    model->iter = 0;

    // Set the search_state seed, width range (single w).
    srch_state->mindepth = srch_state->maxdepth = w;

    // Initialize the PSSM to (0/1 + prior) consensus of SEED.
    INIT_PSSM_TO_CONSENSUS(model, search_word, REF_COUNT_PRIOR, multiseq->background, multiseq->bg_order);
    DEBUG_FMT(NORMAL_VERBOSE, "#   Refining evaluated seed %d: %*.*s log_pvalue %.2f width %d w0 %d\n", 
      num_refined_seeds, w,w,search_word, evaluated_seed->log_pvalue, evaluated_seed->width, evaluated_seed->seed->max_width);

    // 
    // Refine the model using nested subset enrichment.
    // Model is passed via srch_state and the best model found is returned there.
    nested_subset_enrichment(options, srch_state, options->nsubsets);

    // Print progress information.
    candno++;
    if (verbosity >= NORMAL_VERBOSE) {
      double m2, e2, m3, e3, m4, e4, prec=1;
      exp10_logx(model->train_log_pvalue/log(10), m2, e2, prec);
      exp10_logx(model->test_log_pvalue/log(10), m4, e4, prec);
      if (options->objfun == CD) {
        // CD
	fprintf(stderr, "#   CAND-%d w %d w0 %d '%*.*s' '%s' iter %d sd_pos %d sd_dtc %.1f sd_pv %.2f tr_pos %d tr_dtc %.1f tr_pv %3.1fe%+04.0f test_pos %d test_dtc %.1f test_pv %3.1fe%+04.0f\n",
	  motifno, w, seed->max_width, w, w, search_word, model->consensus, model->iter,
	  seed->pos_count, seed->dtc_sum/seed->pos_count, seed->log_pvalue,
          model->train_pos_count, model->train_dtc, m2, e2, 
	  model->test_pos_count, model->test_dtc, m4, e4
	);
      } else if (use_binomial) {
        // DE with binomial test
	fprintf(stderr, "#   CAND-%d-%d w %d w0 %d '%*.*s' '%s' sd_pos %d sd_neg %d sd_pv %.2f tr_pos %d tr_neg %d tr_bern %g tr_pv %3.1fe%+04.0f tr_rat %.1f test_pos %d test_neg %d test_bern %g test_pv %3.1fe%+04.0f test_rat %.1f\n",
	  motifno, candno, w, seed->max_width, w, w, search_word, model->consensus, 
	  seed->pos_count, seed->neg_count, seed->log_pvalue,
	  model->train_pos_count, model->train_neg_count, multiseq->bernoulli, m2, e2, model->train_ratio,
	  model->test_pos_count, model->test_neg_count, multiseq->bernoulli, m4, e4, model->test_ratio
	);
      } else {
        // DE with Fisher Exact test
//FIXME
	fprintf(stderr, "#   CAND-%d-%d w %d w0 %d '%*.*s' '%s' sd_pos %d sd_neg %d sd_pv %.2f tr_pos %d tr_neg %d tr_pv %3.1fe%+04.0f tr_rat %.1f test_pos %d test_neg %d test_pv %3.1fe%+04.0f test_rat %.1f score_threshold %g\n",
	  motifno, candno, w, seed->max_width, w, w, search_word, model->consensus,
	  seed->pos_count, seed->neg_count, seed->log_pvalue,
	  model->train_pos_count, model->train_neg_count, m2, e2, model->train_ratio,
	  model->test_pos_count, model->test_neg_count, m4, e4, model->test_ratio
//FIXME
, model->score_threshold
	);
      }
    } // print progress

    // Print candidate motifs to TEXT ONLY if requested.
    BOOL have_holdout = (test_multiseq != NULL);
    if (options->cand) PRINT_MOTIF(model, motifno, candno, "C", options->text_output, have_holdout);

    // Save this model if it is the best overall based on the training set.
    // Best model can be chosen either by training p-value or by training enrichment ratio.
    // If --nmotifs not given, peek at test_log_pvalue to try to find a model that can be reported.
    // A model is reportable if there is at least one positive training example and
    // it meets the nmotifs of p-value threshold.
    // convert E-value threshold to p-value threshold?
    double log_pvt = (options->thresh_type == EVALUE) ? log(options->thresh)-log(motifno) : log(options->thresh);
    BOOL current_reportable = (model->train_pos_count > 0 && (options->nmotifs > 0 || model->test_log_pvalue <= log_pvt));
    BOOL best_reportable = (final_model->train_pos_count > 0 && (options->nmotifs > 0 || final_model->test_log_pvalue <= log_pvt));
    if ( 
      first_motif ||					// save first motif
      (current_reportable && ! best_reportable) ||	// save first reportable motif
      (
	(current_reportable || ! best_reportable)  && 	// save best motif
	( 
          (options->usepv && 
            ( (model->train_log_pvalue < best_log_pvalue) ||
              (model->train_log_pvalue == best_log_pvalue && model->width > final_model->width)
            ) // prefer wider motifs
          ) ||
	  (! options->usepv && 
            ( (model->train_ratio > best_ratio) ||
              (model->train_ratio == best_ratio && model->width > final_model->width)
            ) // prefer wider motifs
          )
	)						// save better reportable motif
      )
    ) {
      best_ratio = model->train_ratio;
      best_log_pvalue = model->train_log_pvalue;
      *final_model = *model;
      DEBUG_FMT(HIGH_VERBOSE, "# New final model from seed %*.*s\n", w, w, search_word)
    } else {
      // Free malloced things in model.
      free(model->consensus);
    }
    first_motif = False;
  
    // Just refining a single seed, period?
    if (options->nref == 0) break;
    
  } // evaluated seed

  // See if we found a model at all. Return NULL if we didn't.
  if (final_model->train_pos_count == 0 || final_model->test_log_pvalue == 1) {
    free(final_model);
    return(NULL);
  }

  //
  // See if a palindromic model is better.
  //
  if (multiseq->do_rc) {
    Model *pal_model = create_palindromic_model(final_model, multiseq, True);
    srch_state->model = pal_model;
    // Set the search_state seed, width range (single w).
    srch_state->mindepth = srch_state->maxdepth = pal_model->width;
    nested_subset_enrichment(options, srch_state, options->nsubsets);
    double pal_ratio = pal_model->train_log_pvalue / final_model->train_log_pvalue;
    DEBUG_FMT(NORMAL_VERBOSE, "# Testing palindromic version of model: pal log_pvalue %.3f non-pal log_pvalue %.3f log_pvalue_ratio %.3f ED %.3f...\n",
      pal_model->train_log_pvalue, final_model->train_log_pvalue, pal_ratio, pal_model->ed);
    if (pal_ratio >= options->min_pal_ratio && pal_model->ed <= options->max_pal_ed) {
      DEBUG_FMT(NORMAL_VERBOSE, "# Using palindromic version of final model (%.3f >= %.3f AND %.3f <= %.3f)...\n", 
        pal_ratio, options->min_pal_ratio, pal_model->ed, options->max_pal_ed);
      free(final_model->consensus);
      free(final_model);
      final_model = pal_model;
    } else {
      DEBUG_MSG(NORMAL_VERBOSE, "# Using non-palindromic version of final model...\n");
      free(pal_model->consensus);
      free(pal_model);
    }
  } // test palindromic model

  //
  // Do one final refinement with ignore_depth set very high.
  // so that all sites are included in search (slow).
  //
  w = final_model->width;
  DEBUG_FMT(NORMAL_VERBOSE, "#  Refining final model: '%*.*s' '%s'\n", w, w, final_model->seed, final_model->consensus);
  free(final_model->consensus);			// will be replaced
  srch_state->ignore_depth = MAX_WIDTH+1;
  srch_state->mindepth = srch_state->maxdepth = w;
  srch_state->model = final_model;
  nested_subset_enrichment(options, srch_state, options->nsubsets);
  srch_state->ignore_depth = options->ignore_depth;

  // free space
  hash_destroy(ht_words);

  return(final_model);
} // refine_seeds

//
// 1) Evaluate the initial seeds using score-based matching to find the best seed of each width.
// 2) Refine the best seed of each width using iterative score-based matching.
// Returns (a pointer to) the model with the best p-value on the training set.
//
Model *find_best_model(
  STREME_OPTIONS_T *options,		// the program options
  int motifno,				// number of current motif
  Suffixtree *stree,			// the suffix tree
  Reference *rootref,			// the root
  Multiseq *multiseq,			// the training positive and negative sequences
  Multiseq *test_multiseq,		// the held-out positive and negative sequences
  ArraySubtreedataptr *initial_seeds	// the seeds to evaluate using approximate matching
) {
  int i, w, cnt;
  int maxwidth = options->maxwidth;
  Model *model = (Model*) calloc(1, sizeof(Model));	// 0-filled
  Model *final_model = NULL;

  // Initialize a model.
  model->alph = options->alph;
  model->alen = alph_size_core(options->alph);	// for convenience
  model->is_palindromic = False;
  model->train_pos_count = 0;
  model->train_neg_count = 0;
  model->train_dtc = -1;
  model->test_pos_count = 0;
  model->test_neg_count = 0;
  model->test_dtc = -1;
  model->seed = NULL;
  model->consensus = NULL;
  model->consensus = NULL;
  model->matches = NULL;
  model->nmatches = 0;

  // Initialize the search state.
  Searchstate srch_state;
  srch_state.stree = stree;
  srch_state.rootref = rootref;
  srch_state.multiseq = multiseq;
  srch_state.test_multiseq = test_multiseq;
  srch_state.mindepth = 0;
  srch_state.maxdepth = 0;
  srch_state.model = model;			// initial MODEL of word being matched
  for (w=0; w<=MAX_WIDTH+1; w++) INITARRAY(&(srch_state.matching_nodes[w]), Subtreedataptr);
  INITARRAY(&(srch_state.mother), Bref);
  // Initialize the stack for keeping the depth and the number of mismatches.
  INITARRAY(&(srch_state.modelstate), Modelstate);
  Modelstate modelstate;
  modelstate.depth = 0;
  modelstate.errors = 0;
  modelstate.score = 0;
  PUSHARRAY(&(srch_state.modelstate), Modelstate, 128, modelstate);
  // For score-based matching only.
  srch_state.minscore = options->minscore;
  srch_state.ignore_depth = options->ignore_depth;
  srch_state.maxwidth = options->maxwidth;	// maximum depth ever searched

  // Initialize the array of evaluated seeds.
  ArrayEvaluatedseed evaluated_seeds;
  INITARRAY(&evaluated_seeds, Evaluatedseed);

  // Evaluate initial seeds using score-based matching.
  evaluate_initial_seeds(options, initial_seeds, &srch_state, motifno, &evaluated_seeds);

  // Find the best model by refining the best seeds using nested subset refinement.
  final_model = refine_seeds(options, motifno, &srch_state, &evaluated_seeds);

  PRINTTIME("REFINESEEDS");

  if (final_model != NULL) {
    // Copy the seed into the model so it can be freed later.
    int w = final_model->width;
    final_model->seed = (Uchar *)strndup((char *)final_model->seed, w);
    // Set the seed last character to '*' if it is a separator or it runs off the end
    // of the last positive sequence due to the width increasing by +1 for palindrome.
    if (final_model->is_palindromic && final_model->seed[w-1] == SEPARATOR) final_model->seed[w-1] = '*';
    final_model->elapsed_time = mytime(0)/1E6;

    // Print progress: final model
    if (verbosity >= NORMAL_VERBOSE) {
      double m2, e2, m3, e3, m4, e4, prec=1;
      exp10_logx(final_model->train_log_pvalue/log(10), m2, e2, prec);
      exp10_logx(final_model->test_log_pvalue/log(10), m4, e4, prec);
      w = final_model->width;
      Uchar *seed = final_model->seed;
      if (options->objfun == CD) {
	fprintf(stderr, "# BEST-%d '%*.*s' '%s' w %d start_w %d pos_count %d dtc %.1f train_pvalue %3.1fe%+04.0f test_pos %d test_dtc %.1f test_pvalue %3.1fe%+04.0f score_thr %.2f\n\n", 
	  motifno, w, w, final_model->seed, final_model->consensus, w, final_model->initial_width,
	  final_model->train_pos_count, final_model->train_dtc,
	  m2, e2,
	  final_model->test_pos_count, final_model->test_dtc,
	  m4, e4, final_model->score_threshold
	);
      } else if (options->objfun == DE) {
	fprintf(stderr, "# BEST-%d '%*.*s' '%s' w %d start_w %d pos_count %d neg_count %d train_pvalue %3.1fe%+04.0f train_ratio %.1f test_pos %d test_neg %d test_pvalue %3.1fe%+04.0f test_ratio %.1f score_thr %.2f\n\n", 
	  motifno, w, w, final_model->seed, final_model->consensus, w, final_model->initial_width,
	  final_model->train_pos_count, final_model->train_neg_count, 
	  m2, e2, final_model->train_ratio, 
	  final_model->test_pos_count, final_model->test_neg_count,
	  m4, e4, final_model->test_ratio, final_model->score_threshold
	);
      } else {
	fprintf(stderr, "ERROR: Unknown objective function in find_best_model.\n");
	exit(EXIT_FAILURE);
      }
    } // print progress

  } // found a final model

  // Free storage.
  FREEARRAY_TLB(&(srch_state.modelstate), Modelstate);
  FREEARRAY_TLB(&evaluated_seeds, Evaluatedseed);

  // Free the scratch model.
  free(model);

  return(final_model);
} // find_best_model

//
// Erase the matches to the motif by setting them to SEPARATOR characters.
// Only the sites scoring better than the optimal score are erased.
// Sets the site_distr array and total_sites in the model if this is the training set.
// Returns the number of sites erased.
//
int erase_matches(
  STREME_OPTIONS_T *options,	// the options
  Model *model,			// the model
  Multiseq *multiseq,		// the dataset to erase sites from
  BOOL is_ho			// erasing the hold-out set
) {
  Sint n_pos_erased = 0;
  Sint n_neg_erased = 0;
  Uint w = model->width;
  int i, j;

  // Identify the sites to erase and store them in the model.
  score_model_pssm(options, multiseq, model, False, False, PASSING, is_ho, !is_ho);

  // Set the site distribution in the model using the PASSING sites in the training set.
  if (!is_ho) get_site_distr(model, multiseq, CENTER);

  // Add (or append if is_ho=true) the passing sequences to the model.
  get_passing_sequences(options, model, multiseq, is_ho, true);

  // Make PSSM using 0-order background for use in erasing.
  INIT_PSSM_FROM_PROBS(model, multiseq->background, 0);

  // Erase the matches by setting them to the SEPARATOR character.
  // Only erase positions in the matches that have positive scores according
  // to the 0-order PSSM to avoid erasing ends of words when the
  // motif is too wide.
  for (i=0; i<model->nmatches; i++) {
    Uint seqno = model->matches[i].seqno;
    Uint pos = model->matches[i].pos;
    Uint seqstart = multiseq->seqstarts[seqno];
    Uchar *word = multiseq->sequence+seqstart+pos;
    Uint seqlen, rc_seqstart;
    Uchar *rc_word = NULL;
    if (multiseq->do_rc) {
      seqlen = multiseq->seqlengths[seqno];
      rc_seqstart = seqstart + multiseq->totallength/2 + 1;
      rc_word = multiseq->sequence + rc_seqstart + (seqlen - pos) - w;
    }

    // Erase match (from both strands).
    BOOL erased = False;
    for (j=0; j<w; j++) {
      // Erase positions with positive scores.
      Uint aindex, aindex_rc;
      aindex = A2I(word[j]);
      if (multiseq->do_rc) aindex_rc = A2I(rc_word[w-j-1]);
      if (model->pssm[aindex][j] > 0 || (multiseq->do_rc && model->pssm[aindex_rc][w-j-1] > 0)) {
        erased = True;
	word[j] = SEPARATOR;
	if (multiseq->do_rc) rc_word[w-j-1] = SEPARATOR;
      }
    }
    // If nothing was erased, erase the middle of the match.
    if (! erased) {
      j = w/2;
      word[j] = SEPARATOR;
      if (multiseq->do_rc) rc_word[w-j-1] = SEPARATOR;
    }

    if (model->matches[i].is_positive) {
      n_pos_erased++;
    } else {
      n_neg_erased++;
    }
  }
  if (!is_ho) DEBUG_FMT(NORMAL_VERBOSE, "# Erasing %d positive and %d negative matches to the motif from seed: %*.*s\n", 
    n_pos_erased, n_neg_erased, w, w, model->seed);

  // Free the list of matching sites
  free(model->matches);
  model->matches = NULL;
  return(n_pos_erased + n_neg_erased);
} // erase_matches

void free_suffix_tree(
  Suffixtree *stree
) {
  // free the stree
  Subtreedata *stdptr;
  while ((stdptr = POPARRAY(&(stree->subtreedata), Subtreedata)) != NULL) {
    FREEBITTAB(stdptr->pos_seqno_bittab);
    if (stdptr->neg_seqno_bittab) FREEBITTAB(stdptr->neg_seqno_bittab);
    if (stdptr->pos_seqno_list) free(stdptr->pos_seqno_list);
    if (stdptr->neg_seqno_list) free(stdptr->neg_seqno_list);
   if (stdptr->scores) free(stdptr->scores);
  }
  FREEARRAY_TLB(&(stree->subtreedata), Subtreedata);
  freestree(stree);
} // free_suffix_tree

//
// Free the storage and check for leaks.
//
void cleanup(
  Multiseq *multiseq,
  Multiseq *test_multiseq
) {
  Uint ntot = multiseq->npos + multiseq->nneg;

  // Free the storage used.
  DEBUG_MSG(NORMAL_VERBOSE, "# Freeing storage...\n");

  // free the sequences
  FREESPACE_TLB(multiseq->sequence);
  freemultiseq(multiseq);
  if (test_multiseq) { 
    FREESPACE_TLB(test_multiseq->sequence);
    // Background was shared with multiseq and was free'd above.
    test_multiseq->background = test_multiseq->lcbp = NULL;
    freemultiseq(test_multiseq);
  }

  // Check for leaks.
  mmcheckspaceleak();
} // cleanup

/**********************************************************************
  print_streme_model_xml
***********************************************************************/
static void print_streme_model_xml(
  STREME_OPTIONS_T *options,	// STREME options
  Multiseq *multiseq,		// the training positive and negative sequences
  Multiseq *test_multiseq, 	// the held-out positive and negative sequences
  FILE *outfile          	// output file
) {
  int i;
  STR_T *b = str_create(10);

  // Get the name of the CPU.
  char *hostname = HOSTNAME;
  if (hostname[0] == '\0') hostname = "unknown";

  // Start the model section.
  fprintf(
    outfile,
    "  <model>\n"
    "    <command_line>%s</command_line>\n",
    xmlify(options->commandline, b, false)
  );

  // Print the sequence input information.
  fprintf(
    outfile,
    "    <train_positives count=\"%d\" positions=\"%d\" maxlen=\"%d\" file=\"%s\"/>\n",
    multiseq->npos, multiseq->pos_length + 1 - multiseq->npos, multiseq->max_poslen, 
    xmlify(options->posfile, b, True)
  );
  fprintf(
    outfile,
    "    <train_negatives count=\"%d\" positions=\"%d\" from=\"%s\"",
    (options->objfun==CD) ? 0 : multiseq->nneg, 
    (options->objfun==CD) ? 0 : multiseq->neg_length + 1 - multiseq->nneg,
    (options->objfun==CD) ? "none" : options->negfile != options->posfile ? "file" : "shuffled"
  );
  if (options->objfun != CD && options->negfile != options->posfile) {
    fprintf(outfile, " file=\"%s\"/>\n", xmlify(options->negfile, b, True));
  } else {
    fprintf(outfile, "/>\n");
  }
  fprintf(
    outfile,
    "    <test_positives count=\"%d\" positions=\"%d\"/>\n",
    test_multiseq ? test_multiseq->npos : 0, test_multiseq ? test_multiseq->pos_length + 1 - test_multiseq->npos : 0
  );
  fprintf(
    outfile,
    "    <test_negatives count=\"%d\" positions=\"%d\"/>\n",
    test_multiseq ? test_multiseq->nneg : 0, test_multiseq ? test_multiseq->neg_length + 1 - test_multiseq->nneg : 0
  );

  // Print the alphabet.
  alph_print_xml(multiseq->alph, "alphabet", "    ", "  ", outfile);
  fprintf(outfile, "    <strands>%s</strands>\n",
    (alph_has_complement(multiseq->alph) ? (multiseq->do_rc ? "both" : "forward") : "none"));

  // This just gives the raw letter frequencies in the control sequences.
  fprintf(outfile, "    <sequence_db");
  for(i = 0; i < alph_size_core(multiseq->alph); i++) {
    fprintf(outfile, " %s=\"%.3g\"", alph_xml_id(multiseq->alph, i, b), multiseq->freqs[i]);
  }
  fprintf(outfile, "/>\n");

  // Print the background model.
  char *bfile = options->bfile ? options->bfile : "--negatives--";
  fprintf(outfile, "    <background_frequencies source=\"%s\" order=\"%d\">\n", xmlify(bfile, b, true), options->order);
  fprintf(outfile, "      <alphabet_array>\n");
  for (i = 0; i < alph_size_core(multiseq->alph); i++) {
    fprintf(outfile, "        <value letter_id=\"%s\">%.3g</value>\n", alph_xml_id(multiseq->alph, i, b), multiseq->background[i]);
  }
  fprintf(outfile, "      </alphabet_array>\n");
  fprintf(outfile, "    </background_frequencies>\n");

  // Print the stopping criteria.
  fprintf(outfile, "    <stop");
  if (options->nmotifs > 0) {
    fprintf(outfile, " nmotifs=\"%d\"", options->nmotifs);
  } else {
    fprintf(outfile, " thresh_type=\"%s\"", options->thresh_type == EVALUE ? "evalue" : "pvalue");
    fprintf(outfile, " thresh=\"%g\"", options->thresh);
  }
  if (options->time > 0) {
    fprintf(outfile, " time=\"%g\"", options->time);
  }
  fprintf(outfile, "/>\n");

  // Print other input parameters
  fprintf(
    outfile,
    "    <objfun>%s</objfun>\n"
    "    <test>%s</test>\n"
    "    <minw>%d</minw>\n"
    "    <maxw>%d</maxw>\n"
    "    <kmer>%d</kmer>\n"
    "    <hofract>%g</hofract>\n"
    "    <neval>%d</neval>\n"
    "    <nref>%d</nref>\n"
    "    <niter>%d</niter>\n"
    "    <patience>%d</patience>\n"
    "    <seed>%d</seed>\n"
    "    <useer>%s</useer>\n"
    "    <minscore>%g</minscore>\n"
    "    <ignore_depth>%d</ignore_depth>\n"
    "    <nsubsets>%d</nsubsets>\n"
    "    <min_pal_ratio>%g</min_pal_ratio>\n"
    "    <max_pal_ed>%g</max_pal_ed>\n"
    "    <cand>%s</cand>\n"
    "    <experimental>%s</experimental>\n"
    "    <totallength>%d</totallength>\n"
    "    <align>%s</align>\n"
    "    <host>%s</host>\n",
    (options->objfun == DE) ? "Differential Enrichment"
      : (options->objfun == CD) ? "Central Distance"
        : "Unknown objective function",
    (options->objfun == DE && ! multiseq->use_binomial) ? "Fisher Exact Test"
     : (options->objfun == DE && multiseq->use_binomial) ? "Binomial Test"
       : (options->objfun == CD) ? "Cumulative Bates Distribution"
         : "Unknown test",
    options->minwidth, 
    options->maxwidth,
    options->kmer,
    options->hofract,
    options->neval,
    options->nref,
    options->niter,
    options->patience,
    options->seed,
    options->usepv ? "no" : "yes",
    options->minscore,
    options->ignore_depth,
    options->nsubsets,
    options->min_pal_ratio,
    options->max_pal_ed,
    options->cand ? "yes" : "no",
    options->experimental ? "yes" : "no",
    options->totallength,
    options->align_txt,
    xmlify(hostname, b, false)
  );
  if (options->description) { 
   fprintf(outfile, "    <description>%s</description>\n", options->description);
  }

  // Finish the model section
  fprintf(outfile, "  </model>\n");
  str_destroy(b, false);
} // print_streme_model_xml

/**********************************************************************
  print_streme_psfm_xml
**********************************************************************/
static void print_streme_psfm_xml(
  Model *model,		// the model
  FILE* outfile // pointer to output file
) {
  int i = 0;
  int j = 0;
  STR_T *b;
  int width = model->width;
  ALPH_T *alph = model->alph;

  b = str_create(3);
  for (i=0; i < width; i++) { // site position
    fprintf(outfile, "      <pos");
    for (j=0; j < alph_size_core(alph); j++) { // letter
      fprintf(
        outfile, " %s=\"%g\"", alph_xml_id(alph, j, b), model->probs[j][i]
      );
    }
    fprintf(outfile, "/>\n");
  }
  str_destroy(b, false);
} // print_streme_psfm_xml

/**********************************************************************
  print_streme_motifs_xml
***********************************************************************/
static void print_streme_motifs_xml(
  Model **models,		// the motifs
  int n_out_motifs,		// the number of motifs to output
  Multiseq *multiseq,		// the training positive and negative sequences
  Multiseq *test_multiseq,	// the test positive and negative sequences
  FILE *outfile			// output file
) {
  fprintf(outfile, "  <motifs>\n");
  int i, j;
  double m1, e1, m2, e2, m3, e3, prec=1;
  for (i = 0; i < n_out_motifs; i++) {
    int w = models[i]->width;
    exp10_logx(models[i]->train_log_pvalue/log(10), m1, e1, prec);
    exp10_logx(models[i]->test_log_pvalue/log(10), m2, e2, prec);
    exp10_logx(models[i]->test_log_evalue/log(10), m3, e3, prec);
    fprintf(
      outfile,
      "    <motif"
      " id=\"%d-%s\""
      " alt=\"STREME-%d\""
      " width=\"%d\""
      " initial_width=\"%d\""
      " seed=\"%s\""
      " score_threshold=\"%g\""
      " train_pos_count=\"%d\""
      " train_neg_count=\"%d\""
      " train_log_pvalue=\"%g\""
      " train_pvalue=\"%3.1fe%+04.0f\""
      " train_dtc=\"%.1f\""
      " train_bernoulli=\"%g\""
      " test_pos_count=\"%d\""
      " test_neg_count=\"%d\""
      " test_log_pvalue=\"%g\""
      " test_pvalue=\"%3.1fe%+04.0f\""
      " test_log_evalue=\"%g\""
      " test_evalue=\"%3.1fe%+04.0f\""
      " test_dtc=\"%.1f\""
      " test_bernoulli=\"%g\""
      " is_palindromic=\"%s\""
      " elapsed_time=\"%.1f\""
      " total_sites=\"%d\""
      " site_distr=\"",
      i + 1,				// <count>-<consensus>
      models[i]->consensus,		
      i + 1,				// STREME-<count>
      models[i]->width,
      models[i]->initial_width,
      models[i]->seed,
      models[i]->score_threshold,
      models[i]->train_pos_count,
      models[i]->train_neg_count,
      models[i]->train_log_pvalue/log(10),
      m1, e1,				// train_pvalue
      models[i]->train_dtc,
      multiseq->bernoulli,
      models[i]->test_pos_count,
      models[i]->test_neg_count,
      models[i]->test_log_pvalue/log(10),
      m2, e2,				// test_pvalue
      models[i]->test_log_evalue/log(10),
      m3, e3,				// test_evalue
      models[i]->test_dtc,
      test_multiseq ? test_multiseq->bernoulli : -1,
      multiseq->do_rc ? (models[i]->is_palindromic ? "yes" : "no") : "n/a",
      models[i]->elapsed_time,
      models[i]->total_sites
    );
    for (j=0; j<multiseq->max_poslen-w; j++) fprintf(outfile, " %d", models[i]->site_distr[j]);
    fprintf(outfile, "\" max_sites=\"%d\" site_hist=\"", models[i]->max_sites);
    for (j=0; j<=models[i]->max_sites; j++) fprintf(outfile, " %d", models[i]->site_hist[j]);
    fprintf(outfile, "\">\n");

    print_streme_psfm_xml(
      models[i],
      outfile
    );

    fprintf(outfile, "    </motif>\n");
  }
  fprintf(outfile, "  </motifs>\n");
} // print_streme_motifs_xml

//
// Print STREME results in TEXT format. 
//
static void print_streme_file_txt(
  STREME_OPTIONS_T *options,	// STREME options
  BOOL have_holdout,		// true if there is a hold-out set
  BOOL do_rc,			// doing reverse complement
  char *stopping_reason,	// reason search stopped
  double *background,		// 0-order background model
  Model **models,		// the discovered motifs
  int n_out_motifs		// the number of motifs to output
) {
  int i;

  // Print the header of the text output; already printed if --cand.
  if (! options->cand) print_streme_header(options->alph, do_rc, background, options->text_output);

  // Print the models.
  for (i=0; i<n_out_motifs; i++) PRINT_MOTIF(models[i], i+1, -1, "", options->text_output, have_holdout);

  // Print the stopping reason.
  PSTARS(options->text_output);
  fprintf(options->text_output, "%s\n", stopping_reason);
} // print_streme_file_txt

//
// Print positive sequences for each motif in TSV format.
//
static void print_streme_positive_sequences_tsv(
  STREME_OPTIONS_T *options,    // STREME options
  char *tsv_path,		// path to TSV output file
  Model **models,               // the discovered motifs
  int n_out_motifs              // the number of motifs to output
) {
  int imotif, imatch, ipassing;

  // Create the TSV output file.
  FILE *outfile;
  if ((outfile = fopen(tsv_path, "w")) == NULL) {
    DEBUG_FMT(QUIET_VERBOSE, "# Warning: Unable to open SEQUENCES TSV file \"%s\" for writing.\n", tsv_path);
    return;
  }

  // Reread the sequences because sites got erased.
  BOOL do_rc = alph_has_complement(options->alph);
  ARRAY_T *back = NULL;				// space for the MEME-style background
  Multiseq *test_multiseq = NULL;
  Multiseq *train_multiseq = read_pos_neg_seqs(options, do_rc, False, False, True, back, &test_multiseq);
  free_array(back);

  for (imotif=0; imotif<n_out_motifs; imotif++) {
    Model *model = models[imotif];
    double score_threshold = model->score_threshold;
    // Score the training set of sequences getting just the ZOOPS sites (ignores score threshold).
    score_model_pssm(options, train_multiseq, model, False, False, ZOOPS, False, False);
    int nmatches =  model->nmatches;
    int npassing = 0;
    int npos = train_multiseq->npos;
    int nneg = train_multiseq->nneg;
    int nseqs = npos + nneg;
    // Add each training sequence with a passing site to the list of passing sequences.
    Passing_seq *sequences = (Passing_seq *) malloc(nmatches * sizeof(Passing_seq));
    for (imatch=0; imatch<nmatches; imatch++) {
      double score = model->matches[imatch].score;
      if (score >= score_threshold) {
	int seqno = model->matches[imatch].seqno;
	sequences[npassing].dbno = 0;                 // present to allow unique sort order
	sequences[npassing].seqno = seqno;            // present to allow unique sort order
	sequences[npassing].desc = DESCRIPTIONPTR(train_multiseq, seqno);
	sequences[npassing].score = score;
	sequences[npassing].is_tp = (seqno < npos);
	npassing++;
      } // passing
    } // imatch

   // Add passing sites from the test set if there is one.
    if (test_multiseq) {
      // Free the matches in the model.
      free(model->matches);
      model->matches = NULL;
      model->nmatches = 0;
      // Score the test set of sequences getting just the ZOOPS sites (ignores score threshold).
      score_model_pssm(options, test_multiseq, model, False, False, ZOOPS, False, False);
      // Add the passing sequences to the array.
      int nmatches = model->nmatches;
      if (nmatches > 0) {
	npos = test_multiseq->npos;
	sequences = (Passing_seq *) realloc(sequences, (npassing+nmatches) * sizeof(Passing_seq));
	for (imatch=0; imatch<nmatches; imatch++) {
	  double score = model->matches[imatch].score;
	  if (score >= score_threshold) {
	    int seqno = model->matches[imatch].seqno;
	    sequences[npassing].dbno = 1;                 // present to allow unique sort order
	    sequences[npassing].seqno = seqno;            // present to allow unique sort order
	    sequences[npassing].desc = DESCRIPTIONPTR(test_multiseq, seqno);
	    sequences[npassing].score = score;
	    sequences[npassing].is_tp = (seqno < npos);
	    npassing++;
	  }
        } // imatch
      } // nmatches > 0
    } // test_multiseq != NULL

    // Sort the passing sequences.
    qsort(sequences, npassing, sizeof(Passing_seq), compare_sequence_score);

    // Print the TSV header for first motif.
    if (imotif == 0) {
      fprintf(outfile,
	"%s\t"
	"%s\t"
	"%s\t"
	"%s\t"
	"%s\t"
	"%s\t"
	"%s\n",
	"ID", "ALT_ID", (test_multiseq ? "P-value" : "Score"), "seq_ID", "seq_Score", "seq_Class", "is_holdout?"
      );
    } else {
      fprintf(outfile, "#\n");
    }

    // Print the sorted passing sequences.
    double m1, e1, prec=1;
    if (test_multiseq) {
      exp10_logx(model->test_log_pvalue/log(10), m1, e1, prec); \
    } else {
      exp10_logx(model->train_log_pvalue/log(10), m1, e1, prec); \
    }
    for (ipassing=0; ipassing<npassing; ipassing++) {
      fprintf(outfile,
	"%d-%s\t"		// Motif ID
	"STREME-%d\t"		// Motif Alternate ID
        "%3.1fe%+04.0f\t"	// Motif P-value or Score
	"%s\t"    // Sequence ID
	"%.2f\t"  // Score
	"%s\t"    // tf/fp
	"%d\n",   // holdout set?
	imotif+1, model->consensus, imotif+1, m1, e1,
	sequences[ipassing].desc, sequences[ipassing].score, (sequences[ipassing].is_tp ? "tp" : "fp"),
	sequences[ipassing].dbno
      );
    } // passing sequence

    // Free the passing sequences.
    free(sequences);
  } // motif
} // print_streme_positive_sequences_tsv

//
// Print STREME results in XML format. 
//
static void print_streme_file_xml(
  STREME_OPTIONS_T *options,	// STREME options
  char *xml_path,		// path of XML file to create
  char *stopping_reason,	// reason search stopped
  Multiseq *multiseq,		// the training positive and negative sequences
  Multiseq *test_multiseq, 	// the held-out positive and negative sequences
  Model **models,		// the discovered motifs
  int n_out_motifs		// the number of motifs to output
) {
  FILE *outfile;
  if ((outfile = fopen(xml_path, "w")) == NULL) {
    DEBUG_FMT(QUIET_VERBOSE, "# Warning: Unable to open XML file \"%s\" for writing.\n", xml_path);
    return;
  }
  fputs("<STREME version=\"" VERSION "\" release=\"" ARCHIVE_DATE "\">\n", outfile);
  print_streme_model_xml(options, multiseq, test_multiseq, outfile);
  print_streme_motifs_xml(models, n_out_motifs, multiseq, test_multiseq, outfile);
  // Replace any newlines in stopping_reason with spaces.
  char *t = stopping_reason;
  while ((t = strchr(t, '\n')) != NULL) { *t = ' '; };
  fprintf(outfile, "  <reason_for_stopping>%s</reason_for_stopping>\n", stopping_reason);
  fprintf(outfile, "  <run_time cpu=\"%.2f\"/>\n", mytime(0)/1E6);
  fprintf(outfile, "</STREME>\n");
  fclose(outfile);
} // print_streme_file_xml

//
// Print STREME results in HTML format.
//
static void print_streme_file_html(
  STREME_OPTIONS_T *options,	// STREME options
  char *xml_path, 		// path of XML file to process
  char *html_path 		// path of HTML file to create
) {
  // get the xml to html conversion script path
  char *prog = get_meme_libexec_file("streme_xml_to_html");
  if (prog != NULL) {
    // Create the conversion command.
    STR_T *cmd = str_create(0);
    str_append2(cmd, prog);
    str_append(cmd, " ", 1);
    str_append_path(cmd, 1, xml_path);
    str_append(cmd, " ", 1);
    str_append_path(cmd, 1, html_path);

    // Execute the command
    int ret = system(str_internal(cmd));
    if (!(WIFEXITED(ret) && WEXITSTATUS(ret) == 0)) {
      DEBUG_MSG(QUIET_VERBOSE, "# Warning: streme_xml_to_html exited abnormally and may "
        "have failed to create HTML output.\n");
    }

    // Cleanup.
    str_destroy(cmd, false);
    free(prog);
  } else {
    DEBUG_MSG(QUIET_VERBOSE, "# Warning: could not find streme_xml_to_html. "
      "The HTML output could not be created.\n");
  }
} // print_streme_file_html

/*
  This is the main function.
*/
int main(int argc, char *argv[])
{
  int i;
  STREME_OPTIONS_T options;
  Model **models = NULL;
  char *text_filename = "streme.txt";		// motifs in MEME format
  char *xml_filename = "streme.xml";		// interactive output
  char *html_filename = "streme.html";		// interactive output
  char *tsv_filename = "sequences.tsv";		// passing sequences
  char *xml_path=NULL, *html_path=NULL, *tsv_path=NULL;
  #define MAX_REASON_LENGTH 255
  char stopping_reason[MAX_REASON_LENGTH]; 

  // command line processing
  options.commandline = process_command_line(argc, argv, &options);

  // Initilize the clock for reporting times.
  mytime(0);

  srand_mt(options.seed);		// random seed for shuffling sequences
  set_randfunc((randfunc_t) random_mt); // for ushuffle
  INITPOWER2TABLE();			// needed by get_seqno_list()

  // Set up for output.
  if (options.text_only) {
    options.text_output = stdout;
  } else {
    // Create the output directory.
    if (create_output_directory(options.output_dirname, options.allow_clobber, False)) {
      exit(EXIT_FAILURE);
    }
    // Create the path of the TEXT output file and open it.
    char *path = make_path_to_file(options.output_dirname, text_filename);
    if ((options.text_output = fopen(path, "w")) == NULL) {
      fprintf(stderr, "ERROR: Unable to open MOTIFS IN MEME TEXT FORMAT output file '%s' for writing.\n", path);
      exit(EXIT_FAILURE);
    }
    myfree(path);
    // Create the paths for XML, HTML and TSV files.
    xml_path = make_path_to_file(options.output_dirname, xml_filename);
    html_path = make_path_to_file(options.output_dirname, html_filename);
    tsv_path = make_path_to_file(options.output_dirname, tsv_filename);
  }

  // Initialize the alphabet.
  Uint alen = initialize_st_alphabet(options.alph);
  BOOL do_rc = alph_has_complement(options.alph);

  // Input the sequences, converting non-core characters to SEPARATOR.
  // Also reads or creates the background model.
  ARRAY_T *back = NULL;				// space for the MEME-style background
  Multiseq *test_multiseq = NULL;
  Multiseq *multiseq = read_pos_neg_seqs(&options, do_rc, False, False, True, back, &test_multiseq);
  free_array(back);

  // Print the model header if we are outputting candidate motifs.
  if (options.cand) print_streme_header(options.alph, do_rc, multiseq->background, options.text_output);

  // Initialize the p-value lookup cache if we are using objective function DE.
  // It is indexed by the pos_count + neg_count, and each sub-cache
  // is indexed by pos_count.
  options.cache_length = sqrt(2*MAX_CACHE_SIZE);
  if (options.objfun == DE) options.pv_cache = (double**) calloc(options.cache_length, sizeof(double*));

  // Find motif, erase motif, output motif, repeat.
  Suffixtree stree;	// The suffix tree.
  int imotif, n_found_motifs = 0, n_nonsig_motifs = 0;
  for (imotif=0; options.nmotifs==0 || imotif<options.nmotifs; imotif++) {
    Reference rootref;
    // Construct the suffix tree with positive and negative sequences, and maybe reverse complements.
    DEBUG_MSG(NORMAL_VERBOSE, "# Building suffix tree of positive and negative sequences...\n");
    if (constructstree (&stree, multiseq->sequence, multiseq->totallength) != 0) { SIMPLESTANDARDMESSAGE; }
    rootref.toleaf = False;
    rootref.address = ROOT(&stree);
    PRINTTIME("SUFFIXTREE");

    // Add exact counts to stree to get the initial seeds.
    // Seeds of tiny (starting) width can be useful (MIN_SEED_WIDTH).
    ArraySubtreedataptr *initial_seeds = add_counts_to_stree(&options, &stree, &rootref, multiseq);

    // Evaluate initial seeds using score-based matching.
    // Refine the best seed(s) of each width to find the best MODEL.
    Model *model = find_best_model(&options, imotif+1, &stree, &rootref, multiseq, test_multiseq, initial_seeds);
    PRINTTIME("FINDBESTMODEL");
    FREEARRAY_TLB(initial_seeds, Subtreedataptr);

    // Report current memory usage.  Compile with -DMEMORY_MONITOR.
    PRINT_MEMORY_USAGE();

    // Check if we found a reportable model (>0 positive training sites).
    if (model == NULL) {
      snprintf(stopping_reason, MAX_REASON_LENGTH, "Stopped because the next motif had no positive sites with scores > %g.", DEFAULT_OPTSCORE);
      break;
    }

    // convert E-value threshold to p-value threshold?
    double log_pvt = (options.thresh_type == EVALUE) ? log(options.thresh)-log(imotif+1) : log(options.thresh);

    // Update the number of consecutive, non-significant motifs.
    n_nonsig_motifs = (model->test_log_pvalue <= log_pvt) ? 0 : n_nonsig_motifs + 1;

    // Erase the sites in the training set and the hold-out set.
    // This also sets the site_distr array of site counts in the model and total_sites.
    int n_erased = erase_matches(&options, model, multiseq, False);
    // Check that sites were erased to prevent endless loop.
    if (n_erased == 0) {
      snprintf(stopping_reason, MAX_REASON_LENGTH, "Stopped because the last motif found had no valid sites.");
      break;
    }
    if (test_multiseq) (void) erase_matches(&options, model, test_multiseq, True);

    // Save model (pointer) for output as TEXT, XML and HMTL.
    if (n_found_motifs % 10 == 0) Resize(models, n_found_motifs+10, Model *);
    models[n_found_motifs++] = model;

    // Stop if stopping criterion reached.
    if (options.nmotifs == 0 && model->test_log_pvalue > log_pvt) {
      if (n_nonsig_motifs >= options.patience) {
	snprintf(stopping_reason, MAX_REASON_LENGTH, "Stopped because %d consecutive motifs exceeded the %c-value threshold (%g).", 
          options.patience, (options.thresh_type == EVALUE ? 'E' : 'p'), options.thresh);
	break;
      }
    } else if (options.nmotifs > 0 && imotif+1 == options.nmotifs) {
      snprintf(stopping_reason, MAX_REASON_LENGTH, "Stopped because maximum number of motifs (%d) reached.", options.nmotifs);
      break;
    }

    // Check how much time we have left and stop if not enough for another motif.
    if (options.time > 0) {
      double avg_motif_time = (mytime(0)/1E6)/(imotif+1.0);
      if (options.time - mytime(0)/1E6 < 1.5 * avg_motif_time) {
        snprintf(stopping_reason, MAX_REASON_LENGTH, 
          "Stopped because STREME would probably exceed the allowed time (%.2f secs) before\nfinding the next motif.", options.time);
        break;
      }
    }

    // Free the tree.
    free_suffix_tree(&stree);

  } // imotif

  // Free the tree.
  free_suffix_tree(&stree);

  //
  // Output TEXT, XML and HTML.
  //
  // Sort the models by test_log_pvalue.
  qsort(models, n_found_motifs, sizeof(Model *), compare_model_test_pvalue);
  // Remove the non-significant models, saving at least patience of them, if not using nmotifs.
  int n_out_motifs = n_found_motifs;
  if (options.nmotifs == 0) {
    n_nonsig_motifs = 0;
    for (imotif=0; imotif<n_found_motifs; imotif++) {
      // Convert E-value threshold to p-value threshold?
      double log_pvt = (options.thresh_type == EVALUE) ? log(options.thresh)-log(imotif+1) : log(options.thresh);
      if (models[imotif]->test_log_pvalue > log_pvt) n_nonsig_motifs++;
      if (n_nonsig_motifs > options.patience) break;
    }
    n_out_motifs = imotif;
  }
  // Set the E-value in each output motif.
  for (imotif=0; imotif<n_out_motifs; imotif++) {
    models[imotif]->test_log_evalue = log(n_out_motifs) + models[imotif]->test_log_pvalue;
  }
  //
  // Output TEXT first since it sets the consensus, then XML and HTML.
  //
  print_streme_file_txt(&options, (test_multiseq != NULL), do_rc, stopping_reason, multiseq->background, models, n_out_motifs);
  if (! options.text_only) {
    print_streme_file_xml(&options, xml_path, stopping_reason, multiseq, test_multiseq, models, n_out_motifs);
    print_streme_file_html(&options, xml_path, html_path);
    print_streme_passing_sequences_tsv(&options, tsv_path, models, (test_multiseq != NULL), n_out_motifs);
    free(xml_path);
    free(html_path);
    free(tsv_path);
  }

  // Free the sequences.
  cleanup(multiseq, test_multiseq);

  // Clean up.
  for (i=0; i<n_found_motifs; i++) free(models[i]);
  free(models);

  // Print runtime statistics.
  PRINT_STATS(options.text_output, options.commandline);

  return(EXIT_SUCCESS);
} // main
