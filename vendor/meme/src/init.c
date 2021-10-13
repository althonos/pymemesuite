/***********************************************************************
*                                                                      *
* MEME                                                                 *
* Copyright 1994-2017, The Regents of the University of California     *
* Author: Timothy L. Bailey                                            *
*                                                                      *
***********************************************************************/
// init.c 
/*
  Initialize meme.
*/

#define ABS_MIN_W 2

#include "meme.h"
#include "general.h"
#include "banner.h"
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include "utils.h"
#include "io.h"
#include "psp.h"
#include "ushuffle.h"
#include "string-list.h"
#include "citation.js.h"

// priors 
#define PROTEIN_PLIB "prior30.plib"

#define ROUNDERROR (1E-12)
//FIXME: how large should these limits be?
#define DEMIN 2
#define NZMIN 2

#define MAX_MEME_SEQS_BRIEF 1000

// User input parameters 
static bool check_syntax = false; 	// exit after checking syntax if true 
static char *datafile = NULL;		// positive examples 
static char *sf = NULL;			// name to print for datafile 
static char *negfile = NULL;		// negative examples 
static char *obj = "classic";		// objective function 
static bool use_llr = false;		// use POP for starts in classic, not likelihood ratio
static char *stat = "";			// statistical test
static char *bfile = NULL;		// use default background Markov model file 
static int markov_order = -1;		// use highest-order background Markov model in bfile 
					// or 0-order model (add-1 prior) if no bfile given
static char *pspfile = NULL;		// use positional priors 
static bool psp2 = false;		// true if 2-stranded positional priors 
static char *default_output_dirname = "meme_out";  // default name of output directory 
static bool clobber = false;		// default is not to overwrite existing files 
static char *mod = "zoops";		// model type input string; default ZOOPS 
static ALPHABET_T alphabet_type = Protein;	// default alphabet IUPAC protein 1-letter
static char *alph_file = NULL;		// default use built-in alphabet
static bool revcomp = false;		// don't use reverse complement strand of DNA 
static int pal = 0;   			// = 0, no palindromes
           				// = 1, force DNA palindromes,
static bool ma_trim = true;		// trim width using multiple alignment method 
static double wg = 11;			// default gap weight 
static double ws = 1;			// default space weight 
static bool endgaps = true;		// count end gaps in multiple alignment 
static double distance = 1e-5;		// squared Euclidean distance for convergence 
static char *prior = NULL;		// prior type input string 
static double beta = -1;		// scale factor for prior; defaults differ 
static int nmotifs = 1;			// number of motifs to find 
static int maxiter = 50;		// max number iterations of EM on best start 
static double nsites = 0;		// try one value of nsites0 only if > 0 
static int min_nsites = 0;		// minimum nsites0 to try 
static int max_nsites = 0;		// maximum nsites0 to try 
static double wnsites = 0.8;		// weight on prior on nsites 
static int w = 0; 			// width of motifs 
static int min_w = MIN_W;		// minimum W0 to try 
static bool all_widths = false; 	// all widths between min and max 
static bool min_w_set; 			// if set on command line don't override 
static int max_w = MAX_W;		// maximum W0 to try 
static MAP_TYPE map_type; 		// type of sequence to theta mapping 
static char *mapname = NULL; 		// map type input string 
static double map_scale=-1;   		// scale of sequence to theta mapping:
          				//	Uni - size of add-n prior (n)
          				//	Pam - PAM distance (120)
           				//	Default set in init_em.
static STRING_LIST_T *spcons = NULL; 	// list of starting point consensus strings 
static int main_hs = HSIZE; 		// size of heap at "main" w values 
static double hs_decrease = HS_DECREASE; // Rate of decrease for heap size 
static bool x_branch = false; 		// Use x_branch regardless of seq model 
static bool no_x_branch = false; 	// Don't use x_branch, regardless of seq model
static bool w_branch = false; 		// Controls whether width branching occurs 
static bool print_heaps = false; 	// Print heaps after branching rounds
static bool print_pred = false; 	// Print out the predicted sites after each
                                   	// round of MEME search (eg subsequence, EM)
static int brief = 1000;      		// make output brief if > 1000 sequences
static int bfactor = BFACTOR; 		// branching factor for branching search 
static int search_size= 100000;		// search data size limit 
static int max_words = -1; 		// maximum number of seeds to test
static bool norand = false;		// do not randomize sequence order
static int classic_max_nsites = 1000;	// limit number of sites for Classic and NC
static int seed = 0;			// random number seed 
static double hsfrac = -1; 		// fraction of control sequences to use 
static double cefrac = -1; 		// fraction of sequence length for central region
static int kmer = 2; 			// size of kmer to preserve frequency of when shuffling 
static int max_size = 0;		// maximum allowed dataset size; ignore if 0
static double max_time = 0; 		// maximum allowed Wall time; ignore if 0 

/**********************************************************************/
/*
        init_meme

        Set up all the stuff for meme.
*/
/**********************************************************************/
void init_meme(
  int argc, // number of input arguments 
  char **argv, // input arguments 
  MODEL **model_p, // the model 
  MODEL **best_model_p, // the best model 
  MODEL **neg_model_p, // model of negative examples 
  DATASET **dataset_p, // the dataset 
  DATASET **neg_dataset_p, // dataset of negative examples 
  char *text_filename, // name of the text output file 
  char **output_dirname, // name of the output directory 
  FILE **text_output, // destination for text output
  ARRAYLST_T* seq_array,
  int width,
  bool eliminate_repeats
)
{
  int i, j, len, pos;
  char cc;
  OBJTYPE objfun = Classic; // type of objective function 
  TESTTYPE test = No_testtype; // use Fisher exact test
  MOTYPE mtype; // type of model 
  PTYPE ptype; // type of prior 
  PRIORS *priors; // the prior probabilities model 
  P_POINT *p_point; // previously learned starting points 
  ALPH_T *alph=NULL; // alphabet for dataset
  MODEL *model=NULL, *best_model=NULL, *neg_model=NULL;
  DATASET *dataset=NULL, *neg_dataset=NULL;
  int neg_n_samples=0;	// number of samples in negative dataset
  double evt = BIG; // no E-value cutoff 
  bool no_print = false; // turn off printing if parallel and not main 
  char *plib_name = NULL; // use default library 
  bool show_version = false;
  char *np = NULL;
  bool sampling = true;	// sampling is on by default
  bool mpi = false;	// not using MPI

  spcons = new_string_list();	// initialize the consensus sequence list

#ifdef PARALLEL
  // turn off printing if parallel and not the main processor 
  no_print = (mpMyID() != 0);
#endif

  *output_dirname = default_output_dirname;

  // get the command line arguments 
  i = 1;
#ifndef lint
  // print the command line 
  argv[0] = "";
  DO_STANDARD_COMMAND_LINE(1,
    USAGE(Usage:\tmeme\t<dataset> [optional arguments]\n);
    NON_SWITCH(1, <dataset> \t\tfile containing sequences in FASTA format\n,
      switch (i++) {
        case 1: datafile = _OPTION_; break;
        default: COMMAND_LINE_ERROR;
      });
     FLAG_OPTN(1, h, \t\t\tprint this message, USAGE_MESSAGE);
     FLAG_OPTN(2, h-exp, \t\tprint this message including experimental options, EXP_USAGE_MESSAGE);
     DATA_OPTN(1, o, <output dir>,
      \tname of directory for output files\n\t\t\t\twill not replace existing directory,
       *output_dirname = _OPTION_);
     DATA_OPTN(1, oc, <output dir>,
       \tname of directory for output files\n\t\t\t\twill replace existing directory,
       clobber=true; *output_dirname = _OPTION_);
     FLAG_OPTN(1, text, \t\t\toutput in text format (default is HTML),
       TEXT_ONLY = true);
     DATA_OPTN(1, objfun, classic|de|se|cd|ce, \tobjective function (default: classic),
       obj = _OPTION_);
     DATA_OPTN(2, objfun, classic|de|se|ce|cd|nc|nz, \tobjective function (default: classic),
       obj = _OPTION_);
     DATA_OPTN(1, test, mhg|mbn|mrs, \tstatistical test type (default: mhg),
       stat = _OPTION_);
     FLAG_OPTN(1, use_llr, \t\tuse LLR in search for starts in Classic mode,
       use_llr = true);
     FLAG_OPTN(2, use_markov, \t\tuse Markov LLR in search for Classic starts,
       use_llr = true); // Left in for backwards compatibility.
     DATA_OPTN(1, neg, <negdataset>,
       \tfile containing control sequences, negfile = _OPTION_);
     DATA_OPTN(1, shuf, <kmer>, \t\tpreserve frequencies of k-mers of size <kmer> \n\t\t\t\twhen shuffling (default: 2),
       kmer = atoi(_OPTION_));
     DATA_OPTN(1, hsfrac, <hsfrac>, \tfraction of primary sequences in holdout set \n\t\t\t\t(default: 0.5),
       hsfrac= atof(_OPTION_));
     DATA_OPTN(1, cefrac, <cefrac>, \tfraction sequence length for CE region \n\t\t\t\t(default: 0.25),
       cefrac= atof(_OPTION_));
     DATA_OPTN(1, searchsize, <ssize>, \tmaximum portion of primary dataset to use\n\t\t\t\tfor motif search (in characters),
       search_size = atoi(_OPTION_));
     DATA_OPTN(1, maxsize, <maxsize>, \tmaximum dataset size in characters, max_size = atoi(_OPTION_));
     FLAG_OPTN(1, norand, \t\tdo not randomize the order of the input \n\t\t\t\tsequences with -searchsize,
       norand = true);
     DATA_OPTN(1, csites, <csites>, \tmaximum number of sites for EM in Classic mode,
       classic_max_nsites = atoi(_OPTION_));
     DATA_OPTN(1, seed, <seed>, \t\trandom seed for shuffling and sampling,
       seed = atoi(_OPTION_));
     FLAG_OPTN(1, dna, \t\t\tsequences use DNA alphabet, alphabet_type = Dna);
     FLAG_OPTN(1, rna, \t\t\tsequences use RNA alphabet, alphabet_type = Rna);
     FLAG_OPTN(1, protein, \t\tsequences use protein alphabet, alphabet_type = Protein);
     DATA_OPTN(1, alph, <alph file>, \tsequences use custom alphabet,
       alph_file = _OPTION_; alphabet_type = Custom);
     FLAG_OPTN(1, revcomp, \t\tallow sites on + or - DNA strands,
       revcomp = true);
     FLAG_OPTN(1, pal, \t\t\tforce palindromes (requires -dna), pal = 1);
     DATA_OPTN(1, mod, oops|zoops|anr, \tdistribution of motifs, mod = _OPTION_);
     DATA_OPTN(1, nmotifs, <nmotifs>, \tmaximum number of motifs to find,
       nmotifs = atoi(_OPTION_));
     DATA_OPTN(1, evt, <ev>, \t\tstop if motif E-value greater than <evt>,
       evt = atof(_OPTION_));
     DATA_OPTN(1, time, <t>, \t\tquit before <t> seconds have elapsed,
       max_time = atof(_OPTION_));
     DATA_OPTN(1, nsites, <sites>, \tnumber of sites for each motif,
       nsites=atof(_OPTION_));
     DATA_OPTN(1, minsites, <minsites>,
       \tminimum number of sites for each motif, min_nsites=atoi(_OPTION_));
     DATA_OPTN(1,
       maxsites, <maxsites>, \tmaximum number of sites for each motif,
       max_nsites=atoi(_OPTION_));
     DATA_OPTN(1, wnsites, <wnsites>, \tweight on expected number of sites,
       wnsites=atof(_OPTION_));
     DATA_OPTN(1, w, <w>, \t\tmotif width, w = atoi(_OPTION_); min_w_set=true);
     DATA_OPTN(1, minw, <minw>, \t\tminimum motif width,
               min_w = atoi(_OPTION_); min_w_set=true);
     DATA_OPTN(1, maxw, <maxw>, \t\tmaximum motif width,
       max_w = atoi(_OPTION_));
     FLAG_OPTN(1, allw, \t\t\ttest starts of all widths from minw to maxw,
               all_widths=true);
     FLAG_OPTN(1, nomatrim,
       \t\tdo not adjust motif width using multiple\n\t\t\t\talignment,
       ma_trim = false);
     DATA_OPTN(1, wg, <wg>, \t\tgap opening cost for multiple alignments,
       wg=atof(_OPTION_));
     DATA_OPTN(1, ws, <ws>, \t\tgap extension cost for multiple alignments,
       ws=atof(_OPTION_));
     FLAG_OPTN(1, noendgaps, \t\tdo not count end gaps in multiple alignments,
       endgaps = false);
     DATA_OPTN(1, bfile, <bfile>, \tname of background Markov model file,
       bfile = _OPTION_);
     DATA_OPTN(1, markov_order, <order>, \t(maximum) order of Markov model to use or create,
       markov_order = atoi(_OPTION_));
     DATA_OPTN(1, psp, <pspfile>, \tname of positional priors file,
       pspfile = _OPTION_);
     DATA_OPTN(2, psp2, <pspfile>, \tname of 2-stranded positional priors file\n\t\t\t\t(requires -revcomp),
       pspfile = _OPTION_;psp2 = true);
     DATA_OPTN(1, maxiter, <maxiter>, \tmaximum EM iterations to run,
       maxiter = atoi(_OPTION_));
     DATA_OPTN(1, distance, <distance>, \tEM convergence criterion,
       distance = atof(_OPTION_));
     DATA_OPTN(1,
       prior, dirichlet|dmix|mega|megap|addone, \n\t\t\t\ttype of prior to use,
       prior = _OPTION_);
     DATA_OPTN(1, b, <b>, \t\tstrength of the prior, beta = atof(_OPTION_));
     DATA_OPTN(1, plib, <plib>, \t\tname of Dirichlet prior file,
       plib_name = strdup(_OPTION_));
     DATA_OPTN(1, spfuzz, <spfuzz>, \tfuzziness of sequence to theta mapping,
       map_scale = atof(_OPTION_));
     DATA_OPTN(1, spmap, uni|pam, \tstarting point seq to theta mapping type,
       mapname = _OPTION_);
     DATA_OPTN(1, cons, <cons>, \t\tconsensus sequence to start EM from,
       add_string(_OPTION_, spcons));
     DATA_OPTN(2, heapsize, <hs>,
               \tsize of heaps for widths where substring \n\t\t\t\tsearch occurs, 
               main_hs = atoi(_OPTION_));
     FLAG_OPTN(2, x_branch, \t\tperform x-branching, x_branch=true);
     FLAG_OPTN(2, no_x_branch, \t\tdo not perform x-branching, 
               no_x_branch=true);
     FLAG_OPTN(2, w_branch, \t\tperform width branching, w_branch=true);
     DATA_OPTN(2, bfactor, <bf>,
       \t\tbranching factor for branching search, bfactor = atoi(_OPTION_));
     FLAG_OPTN(2, print_pred, \t\tprint out the sites predicted by meme,
       print_pred = true);
     FLAG_OPTN(2, print_heaps, \t\tprint heaps after each branching round,
       print_heaps = true);
     DATA_OPTN(1, brief, <n>, \t\tomit sites and sequence tables in\n\t\t\t\toutput if more than <n> primary sequences,
       brief = atoi(_OPTION_));
     FLAG_OPTN(1, nostatus, \t\tdo not print progress reports to terminal,
       NO_STATUS = true);
     DATA_OPTN(1, p, <np>, \t\tuse parallel version with <np> processors, np = _OPTION_);
     DATA_OPTN(1, sf, <sf>, \t\tprint <sf> as name of sequence file, sf = _OPTION_);
     FLAG_OPTN(2, mpi, \t\tdo not use (set by exec_parallel), mpi = true);
     FLAG_OPTN(2, check_syntax, \t\tcheck input syntax and exit,
       check_syntax = true);
     FLAG_OPTN(1, V, \t\t\tverbose mode, VERBOSE = true);
     FLAG_OPTN(1, version, \t\tdisplay the version number and exit, show_version = true);
     FLAG_OPTN(2, trace, \t\ttrace starting points, TRACE = true);
     FLAG_OPTN(2, print_all, \t\tprint all debug information,
       PRINTALL = true);
     FLAG_OPTN(2, print_w, \t\tprint erasure matrix, PRINT_W = true);
     FLAG_OPTN(2, print_z, \t\tprint missing information matrix,
       PRINT_Z = true);
     FLAG_OPTN(2, print_ll, \t\tprint log-likelihood during EM,
       PRINT_LL = true);
     FLAG_OPTN(2, print_starts, \t\tprint starting points, 
       PRINT_STARTS = true);
  )
#endif

#ifndef PARALLEL
  if (np != NULL && !NO_STATUS) {
    fprintf(stderr, "Option '-p' given but Parallel MEME not configured. Running Serial MEME.\n");
  } 
#endif

  // exit if check_syntax is on 
  if (check_syntax) exit(0);

  if (show_version) {
    fprintf(stdout, VERSION "\n");
    exit(EXIT_SUCCESS);
  }

  // set random number generators 
  srand_mt(seed);
  set_randfunc((randfunc_t) random_mt); // for ushuffle

  if (TEXT_ONLY == true) {
    // Legacy: plain text output to standard out.
    *text_output = stdout;
  }
  else {
    if (!no_print) {
      // allow clobbering of the default output directory
      if (*output_dirname == default_output_dirname) { 
        clobber = true;
      } 
      if (create_output_directory(*output_dirname, clobber, !NO_STATUS)) {
        // Failed to create output directory.
        exit(1);
      }
      // Create the name of the text output file 
      // "<dir>/text_filename/" and open it for writing
      char *path = make_path_to_file(*output_dirname, text_filename);
      *text_output = fopen(path, "w"); //FIXME CEG check for errors
      myfree(path);
    }
  }
#ifdef PARALLEL
  // Send all text_output text from non-main node to bit bucket.
  if (mpMyID() != 0) *text_output = fopen("/dev/null", "w");
#endif

  // set all the print flags appropriately 
  if (PRINTALL) {
    PRINT_W = true;
    PRINT_Z = true;
    PRINT_LL = true;
    PRINT_STARTS = true;
  }

  // check that nmotifs >= 1 
  if (nmotifs < 1) {
    fprintf(stderr, "Option '-nmotifs' must be >= 1.\n");
    exit(1);
  }

  // set nmotifs to MAX_NMOTIFS if -evt was specified.
  if (evt != BIG) {
    nmotifs = MAX_NMOTIFS;
  }

  // get the objective function type
  if (!strcmp(obj, "classic")) {
    objfun = Classic;
  } else if (!strcmp(obj, "nc")) {
    objfun = NC;
  } else if (!strcmp(obj, "se") || !strcmp(obj, "smhg")) {
    objfun = SE;
  } else if (!strcmp(obj, "nz")) {
    objfun = NZ;
  } else if (!strcmp(obj, "de") || !strcmp(obj, "hs") || !strcmp(obj, "cv")) { // synonyms
    objfun = DE;
  } else if (!strcmp(obj, "ce")) {
    objfun = CE;
  } else if (!strcmp(obj, "cd")) {
    objfun = CD;
  } else {
    fprintf(stderr, "Unknown objective function type %s.\n", obj);
    exit(1);
  }

  // Get the statistical test type.
  if ((objfun == Classic || objfun == NC || objfun == CE || objfun == CD) && strcmp(stat, "")) {
    fprintf(stderr, "Ignoring option '-test' with option '-objfun %s'.\n", obj);
  }

  // Set default test type.
  if (!strcmp(stat, "") && (objfun == DE || objfun == NZ || objfun == SE)) {
    stat = "mhg";	// default test for DE, NZ, SE
  } else if (objfun == CE) { 
    stat = "mbn";	// only test for CE
  }
  if (!strcmp(stat, "")) {
    test = No_testtype;
  } else if (!strcmp(stat, "mhg")) {
    test = mHG;
  } else if (!strcmp(stat, "mbn") || !strcmp(stat, "msn")) {
    test = mBN;
  } else if (!strcmp(stat, "mrs")) {
    test = mRS;
  } else {
    fprintf(stderr, "Unknown statistical test type %s. \n", stat);
    exit(1);
  }

  // get the model type 
  if (!strcmp(mod, "anr") || !strcmp(mod, "tcm")) {
    mtype = Tcm;
    if (objfun == CE || objfun == CD) {
      fprintf (stderr, "Option '-mod %s' not allowed with option '-objfun %s'.\n", mod, obj);
      exit(1);
    }
  } else if (!strcmp(mod, "oops")) {
    mtype = Oops;
  } else if (!strcmp(mod, "zoops")) {
    mtype = Zoops;
  } else {
    mtype = Zoops; // prevent warning 
    fprintf(stderr, "Unknown model type '-mod %s'. \n", mod);
    exit(1);
  }

  // Check that model and objective function permitted for PSP.
  if (pspfile) {
    if (mtype == Tcm) {
      fprintf(stderr, "Option '-mod %s' not yet supported with option '-psp'.\n", mod);
      exit(1);
    } else if (objfun != Classic) {
      fprintf(stderr, "Option '-psp' only supported for '-objfun classic'.\n");
      exit(1);
    }
  }

  // turn off multiple alignment trimming unless Classic mode
  if (objfun != Classic) ma_trim = false;

  // hsfrac not given; set defaults
  if (hsfrac == -1) {
    if (objfun==DE || objfun==SE || objfun==CE || objfun==CD) {
      hsfrac = 0.5;
    } else if (objfun == NZ) {
      hsfrac = 1.0;
    } else if (objfun == Classic || objfun == NC) {
      hsfrac = 0;
    }
  } else {
    // Don't allow hsfrac if Classic
    if (objfun == Classic) {
      fprintf(stderr, "Option '-hsfrac' requires you to specify a non-classic objective function with '-objfun'.\n");
      exit(1);
    }
  }

  // check hsfrac 
  if (hsfrac == 0 && !(objfun == Classic || objfun == NC || objfun == SE)) {
    fprintf(stderr, "Option '-hsfrac' must be > 0 with option '-objfun %s'.\n", obj);
    exit(1);
  } else if (hsfrac < 0 && objfun == SE) {
    fprintf(stderr, "Option '-hsfrac' must be >= 0 with option '-objfun %s'.\n", obj);
    exit(1);
  } else if (hsfrac > 1) {
    fprintf(stderr, "Option '-hsfrac' must be <= 1 with option '-objfun %s'.\n", obj);
    exit(1);
  }

  // check the alphabet and set up default mappings and priors 
  switch(alphabet_type) {
    case Dna:
      alph = alph_dna(); // builtin DNA alphabet
      if (!mapname) mapname = "uni"; // uniform prior mapping 
      if (!prior) prior = "dirichlet"; // simple dirichlet prior 
      break;
    case Rna:
      alph = alph_rna(); // builtin RNA alphabet
      if (!mapname) mapname = "uni"; // uniform prior mapping 
      if (!prior) prior = "dirichlet"; // simple dirichlet prior 
      break;
    case Protein:
      alph = alph_protein(); // builtin Protein alphabet
      if (!mapname) mapname = "pam"; // PAM mapping 
      if (!prior) {
	switch (mtype) {
	  case Oops: prior = "dmix"; break;
	  case Zoops:
	  case Tcm: prior = "megap"; break;
	  default: prior = "dirichlet"; break;
	}
      }
      break;
    case Custom:
      alph = alph_load(alph_file, !no_print); // load custom alphabet
      if (alph == NULL) exit(EXIT_FAILURE);
      if (!mapname) mapname = "uni";// uniform prior mapping
      if (!prior) prior = "dirichlet"; // simple dirichlet prior
      break;
  } // alphabet_type

  // check that the alphabet can perform requested features
  if (!alph_has_complement(alph)) {
    if (revcomp) {
      fprintf(stderr,
        "You must use a complementable alphabet with option '-revcomp'.\n");
      exit(1);
    }
    if (pal) {
      fprintf(stderr, "You must use a complementable alphabet with option '-pal'.\n");
      exit(1);
    }
  }

  // find out type of prior 
  if (!strcmp(prior, "dirichlet")) {
    ptype = Dirichlet;
    if (beta < 0) beta = 0.01; // default b = 0.01 
  } else if (!strcmp(prior, "dmix")) {
    ptype = Dmix;
    if (beta < 0) beta = 0; // default b = 0 for dmix 
  } else if (!strcmp(prior, "megadmix") || !strcmp(prior, "mega")) {
    ptype = Mega; // use mega prior heuristic 
  } else if (!strcmp(prior, "megap")) {
    ptype = MegaP; // discretization uses b=0 
  } else if (!strcmp(prior, "addone")) {
    ptype = Addone;
    if (beta != -1) {
      fprintf(stderr, "You may not specify '-b' with option '-prior addone'.\n");
      exit(1);
    }
  } else {
    ptype = Dirichlet; // prevent warning 
    fprintf(stderr, "Unknown type of prior: %s!\n", prior);
    exit(1);
  }

  // Read the primary dataset and set up globals.
  dataset = (seq_array != NULL) ? create_meme_dataset_from_momo(seq_array, alph, width, eliminate_repeats) : read_seq_file(datafile, alph, revcomp, false);
  if (!dataset) exit(1);

  // read in the negative dataset if a negative file given 
  if (negfile) {
    if (objfun!=SE && objfun!=DE && objfun!=NZ) {
      fprintf(stderr, "You may not use option '-neg' with option '-objfun %s'.\n", obj);
      exit(1);
    }
    neg_dataset = read_seq_file(negfile, alph, revcomp, true);
    if (!neg_dataset) exit(1);
    neg_n_samples = neg_dataset->n_samples;
  }

  // Sort the samples alphabetically by their integer-encoding,
  // accounting for the possibility of two strands, then randomize the order.
  if (! norand) {
    shuffle_dataset_order(dataset);
    if (neg_dataset) shuffle_dataset_order(neg_dataset);
  }

  // Create a control dataset if none provided by reading
  // in a new copy of the primary dataset and shuffling preserving k-mers.
  if (!neg_dataset && (objfun==DE || objfun==NZ || objfun==SE)) {
    neg_dataset = read_seq_file(datafile, alph, revcomp, true);
    // Note: don't allow kmer<2 or runs of Ns will be broken up in repeatmasked
    // sequences, messing up the statistics.
    if (kmer < 1 || kmer > 6) {
      fprintf(stderr, "Option '-shuf' must be in the range [1,..,6].\n");
      exit(1);
    }
    // Shuffle the order of the sequences, then their letters.
    if (!NO_STATUS) fprintf(stderr, "Creating shuffled version of primary dataset (preserving %d-mer frequencies) as control...\n", kmer);
    shuffle_dataset_order(neg_dataset);
    shuffle_dataset_letters(neg_dataset, kmer, revcomp);
    neg_n_samples = dataset->n_samples;
  } else {
    kmer = -1;
  }

  // prevent too long jobs
  if (max_size > 0) {
    if (dataset->total_res > max_size) {
      fprintf(stderr, "The primary dataset is too large (%d > %d).  Rerun with larger -maxsize.\n", 
        dataset->total_res, max_size);
      exit(1);
    } else if (neg_dataset && neg_dataset->total_res > max_size) {
      fprintf(stderr, "The control dataset is too large (%d > %d).  Rerun with larger -maxsize.\n", 
        neg_dataset->total_res, max_size);
      exit(1);
    }
  }

  // Check that there are enough sequences.
  if (mtype != Tcm && dataset->n_samples < 2) {
    fprintf(stderr, "Your (primary) dataset must contain at least 2 sequences unless you choose the"
	"'Any Number of Repetitions' site distribution model'.\n");
    exit(1);
  }

  // If search_size is 0, set it and maxwords to size of dataset.
  if (search_size == 0) {
    search_size = dataset->total_res;
    max_words = dataset->total_res;
    sampling = false;
    if (!NO_STATUS) {
      if (mtype != Tcm) fprintf(stderr, "Turning off sequence sampling due to -searchsize 0.\n");
      fprintf(stderr, "Setting 'searchsize' to the actual size of primary sequences: %d\n", search_size);
    }
  } else {
    // Don't allow any sequences longer than search_size.
    if (dataset->max_slength > search_size) {
      fprintf(stderr, 
        "Your (primary) dataset must not contain sequences longer than set by -searchsize: %d.\n"
        "It contains sequence(s) of length %ld.\n", search_size, dataset->max_slength);
      exit(1);
    }
    search_size = MIN(dataset->total_res, search_size);
    // max_words scales as sqrt(search_size) so that overall runtime grows as n.
    max_words = MIN(dataset->total_res, 60000 + sqrt(search_size));
    max_words = MIN(max_words, search_size);
  }

  // Check that all sequences have same length if doing central enrichment and set
  // the central region distance threshold and prior probabilities.
  if (objfun == CE || objfun == CD) {
    if (cefrac == -1) {
       cefrac = 0.25;
    } else if (cefrac <= 0 || cefrac >= 1) {
      fprintf(stderr, "You must specify a value between 0 and 1 with option '-cefrac'.\n");
      exit(1);
    }
    if (dataset->min_slength != dataset->max_slength) {
      fprintf(stderr, "All sequences must have the same length with option '-objfun %s'.\n", obj);
      exit(1);
    }
    dataset->ce_frac = cefrac;
    dataset->ce_max_dist = 0.5*cefrac*dataset->min_slength;	// max possible distance to center
    if (!NO_STATUS) fprintf(stderr, "%s: cefrac %.2f length %ld central region %d\n", 
      objfun==CE ? "CE" : "CD", cefrac, dataset->min_slength, 2*dataset->ce_max_dist);
  } else if (cefrac != -1) {
    fprintf(stderr, "You may not use option '-cefrac' with option '-objfun %s'.\n", obj);
    exit(1);
  } else {
    cefrac = 1; // so search_size and highwater mark work
  } 

  // initialize the background model 
  init_meme_background(bfile, revcomp, dataset, alph_file, alphabet_type, markov_order, datafile, 1);
  // Note: use primary sequences to ensure that both background models are the
  // same or spurious motifs will be found.
  if (neg_dataset) init_meme_background(bfile, revcomp, neg_dataset, alph_file, alphabet_type, markov_order, datafile, 0);  // not negfile!

  // Find out the index of the last sequence where search_size is reached or exceeded.
  // It and all following sequences will be put in groups 1&2 for speed reasons.
  int search_size_nseqs = 0;					// move from here to groups 1&2
  int nseqs = dataset->n_samples;
  int slen = 0;
  if (search_size < dataset->total_res) {
    for (i=0; i<nseqs; i++) {
      slen += dataset->samples[i]->length * cefrac;		// adjust for removed flanks in CE mode
      // Get sequence index where search_size is reached or exceeded.
      if (!search_size_nseqs && slen >= search_size) search_size_nseqs = i+1;	// number of sequences in group 0
    }
  }
  if (search_size_nseqs == 0) search_size_nseqs = nseqs;

  // Limit the number of sequences in group 0 for Classic and NC with non-TCM models.
  if (mtype != Tcm && (objfun == Classic || objfun == NC) && search_size_nseqs > classic_max_nsites) {
    if (sampling) search_size_nseqs = classic_max_nsites;
  }

  // Check that there are at least 2 sequences unless TCM.
  if (mtype != Tcm && search_size_nseqs < 2) {
    fprintf(stderr, "Too few sequences will be searched (%d); increase '-searchsize'.\n", search_size_nseqs);
    exit(1);
  }

  // Assign the primary and control sequence groups to be used for each task.
  GROUP_T primary_groups = { .em={false,false,false}, .trim={false,false,false}, .pvalue={false,false,false}, .nsites={false,false,false} };
  GROUP_T control_groups = { .em={false,false,false}, .trim={false,false,false}, .pvalue={false,false,false}, .nsites={false,false,false} };
  char *group_string = "";
  if (objfun == Classic || objfun == NC) {
    group_string = "Starts/EM: p0; Trim: p0; pvalue: p0; nsites: p0,p1,p2";
    GROUP_T pgs = { .em={true,false,false}, .trim={true,false,false}, .pvalue={true,false,false}, .nsites={true,true,true} };
    GROUP_T cgs = { .em={false,false,false}, .trim={false,false,false}, .pvalue={false,false,false}, .nsites={false,false,false} };
    memcpy(&primary_groups, &pgs, sizeof(GROUP_T));
    memcpy(&control_groups, &cgs, sizeof(GROUP_T));
  } else if (objfun == DE) {
    group_string = "Starts/EM: p0 vs c0; Trim: p1 vs c1; pvalue: p2 vs c1,c2; nsites: p0,p1,p2 vs c0,c1,c2";
    GROUP_T pgs = { .em={true,false,false}, .trim={false,true,false}, .pvalue={false,false,true}, .nsites={true,true,true} };
    GROUP_T cgs = { .em={true,false,false}, .trim={true,false,false}, .pvalue={false,true,true}, .nsites={true,true,true} };
    memcpy(&primary_groups, &pgs, sizeof(GROUP_T));
    memcpy(&control_groups, &cgs, sizeof(GROUP_T));
  } else if (objfun == SE) {
    group_string = "Starts/EM: p0 vs c0; Trim: p0 vs c0; pvalue: p1,p2 vs c1,c2; nsites: p0,p1,p2 vs c0,c1,c2";
    GROUP_T pgs = { .em={true,false,false}, .trim={true,false,false}, .pvalue={false,true,true}, .nsites={true,true,true} };
    GROUP_T cgs = { .em={true,false,false}, .trim={true,false,false}, .pvalue={false,true,true}, .nsites={true,true,true} };
    memcpy(&primary_groups, &pgs, sizeof(GROUP_T));
    memcpy(&control_groups, &cgs, sizeof(GROUP_T));
  } else if (objfun == CE || objfun == CD) {
    group_string = "Starts/EM: p0; Trim: p1; pvalue: p2; nsites: p0,p1,p2";
    GROUP_T pgs = { .em={true,false,false}, .trim={false,true,false}, .pvalue={false,false,true}, .nsites={true,true,true} };
    GROUP_T cgs = { .em={false,false,false}, .trim={false,false,false}, .pvalue={false,false,false}, .nsites={false,false,false} };
    memcpy(&primary_groups, &pgs, sizeof(GROUP_T));
    memcpy(&control_groups, &cgs, sizeof(GROUP_T));
  } else if (objfun == NZ) {
    group_string = "Starts/EM: p0,p2; Trim: p0,p1 vs c1; pvalue: p1 vs c2; nsites: p0 vs c1,c2";
    GROUP_T pgs = { .em={true,false,true}, .trim={true,false,true}, .pvalue={true,false,false}, .nsites={true,false,false} };
    GROUP_T cgs = { .em={false,false,false}, .trim={false,true,false}, .pvalue={false,false,true}, .nsites={false,true,true} };
    memcpy(&primary_groups, &pgs, sizeof(GROUP_T));
    memcpy(&control_groups, &cgs, sizeof(GROUP_T));
  } else {
    die("Unknown objective function in init.c.\n");
  }

  // Put the last (after shuffling order) samples in the hold-out groups for 
  // for all objective functions except NZ.
  if (objfun != NZ) {				// objfun != NZ
    int nseqs = dataset->n_samples;		// number of primary sequences
    int c0 = MIN(search_size_nseqs, (1-hsfrac)*nseqs);
    //int c1 = (objfun == Classic) ? nseqs-c0 : (nseqs - c0)/2; // split control sequences 50:50
    int c1 = (objfun == Classic) ? nseqs-c0 : MIN(c0, (nseqs - c0)/2);	// same size as c0 for speed with DE
    int c2 = nseqs-c0-c1;
    dataset->n_group[0] = c0;
    dataset->n_group[1] = c1;
    dataset->n_group[2] = c2; 
    int h1 = c0;			// index of first sequence in group 1
    int h2 = h1 + c1;			// index of first sequence in group 2
    // Set the group of each sample.
    for (i=0; i<h1; i++) dataset->samples[i]->group = 0;		// To be safe.
    for (i=h1; i<h2; i++) dataset->samples[i]->group = 1;
    for (i=h2; i<nseqs; i++) dataset->samples[i]->group = 2;
    // Set the group last indices.
    dataset->group_last_idx[0] = h1-1;	// index of last sequence in group 0
    dataset->group_last_idx[1] = h2-1;	// index of last sequence in group 1
    dataset->group_last_idx[2] = nseqs-1; // index of last sequence in group 2
    if (!NO_STATUS) fprintf(stderr, "PRIMARY (%s): n %d p0 %d p1 %d p2 %d\n", obj, nseqs, c0, c1, c2);
    // Check that the groups are large enough
    if (hsfrac != 0) {
      if (nseqs < 3*DEMIN) {
	fprintf(stderr, 
	  "You can only use option '-objfun %s' if the dataset has at least %d sequences.\n", obj, 3*DEMIN);
	exit(1);
      }
      if (c1 < DEMIN || c2 < DEMIN) {
	double min_hsfrac = (2.0*DEMIN)/nseqs;
	fprintf(stderr, 
	  "You must specify '-hsfrac' at least %g with this dataset.\n", min_hsfrac);
	exit(1);
      }
    }
  } else {			// objfun == NZ
    //
    // Add the first hsfrac of the control sequences to the 
    // primary dataset.
    // Note: this should be done after setting max_nsites if we
    // don't want to bias the search.  On the other hand? 
    //
    int n = negfile ? neg_dataset->n_samples : dataset->n_samples;
    int nct = (hsfrac/2) * n; // split 50:50
    // Check that there are enough samples for NZ.
    if (!NO_STATUS) fprintf(stderr, "NZ: n %d holdout %d\n", n, nct);
    if (n < 2*NZMIN) {			// not enough sequences
      if (negfile) {
	fprintf(stderr, 
	  "You can only use '-objfun nz' if the control dataset has at least %d sequences.\n", 2*NZMIN);
      } else {
        fprintf(stderr, 
	  "You can only use '-objfun nz' if the dataset has at least %d sequences.\n", 2*NZMIN);
      }
      exit(1);
    }
    if (nct < NZMIN) {			// not enough control sequences
      double min_hsfrac = (2.0*NZMIN)/n;
      fprintf(stderr, 
        "You must specify '-hsfrac' at least %g with this dataset.\n", min_hsfrac);
      exit(1);
    }
    add_control_samples(dataset, neg_dataset, nct, nct, revcomp);
    neg_dataset = NULL; // Prevent it from ever being freed!
  }

  // Split (after shuffling order) the control samples into three groups
  // for all objective functions except NZ.
  if (neg_dataset && objfun != NZ) {		// objfun != NZ
    int nseqs = neg_dataset->n_samples;		// number of control sequences
    int c0 = MIN(search_size_nseqs, (1-hsfrac)*nseqs);
    int c1 = (nseqs - c0)/2;	// Split remaining control 50:50
    int c2 = nseqs-c0-c1;
    neg_dataset->n_group[0] = c0;
    neg_dataset->n_group[1] = c1;
    neg_dataset->n_group[2] = c2; 
    int h1 = c0;			// index of first sequence in group 1
    int h2 = h1 + c1;			// index of first sequence in group 2
    // Set the group of each sample.
    for (i=0; i<h1; i++) neg_dataset->samples[i]->group = 0;		// To be safe.
    for (i=h1; i<h2; i++) neg_dataset->samples[i]->group = 1;
    for (i=h2; i<nseqs; i++) neg_dataset->samples[i]->group = 2;
    // Set the group last indices.
    neg_dataset->group_last_idx[0] = h1-1;	// index of last sequence in group 0
    neg_dataset->group_last_idx[1] = h2-1;	// index of last sequence in group 1
    neg_dataset->group_last_idx[2] = nseqs-1;   // index of last sequence in group 2
    if (! NO_STATUS) fprintf(stderr, "CONTROL (%s): n %d c0 %d c1 %d c2 %d\n", obj, nseqs, c0, c1, c2);
  } // neg_dataset control groups
  
  if (! NO_STATUS && objfun != NZ) fprintf(stderr, "SEQUENCE GROUP USAGE-- %s\n", group_string);

  // Set the number of primary samples.
  int n_primary_samples = (objfun==NZ) ? dataset->n_group[0] : dataset->n_samples;

  // check we can actually use the dataset in this mode
  if (n_primary_samples == 1 && mtype != Tcm) {
    fprintf(stderr, "You must specify '-mod anr' to set the motif site model "
        "to 'Any Number of Repetitions [per sequence]' when only providing "
        "one sequence\n");
    exit(1);
  }

  // Initialize which sequence groups to use for starts and EM.
  set_seq_groups_to_include(dataset, primary_groups.em);
  if (objfun == DE || objfun == SE) set_seq_groups_to_include(neg_dataset, control_groups.em);

  // Initialize which sequence regions to use as source of starting words and for EM.
  if (objfun == CE || objfun == CD) {
    set_seq_regions_to_include(dataset, false, true, false);	// Only get starts and do EM in central region
    dataset->region_last_pos[0] = dataset->min_slength/2 - dataset->ce_max_dist - 1;	// left
    dataset->region_last_pos[1] = dataset->min_slength/2 + dataset->ce_max_dist - 1;	// center
    dataset->region_last_pos[2] = dataset->min_slength - 1;				// right
  } else {
    set_seq_regions_to_include(dataset, true, true, true);
  }

  // read in psp file if one given
  if (pspfile) {
    read_psp_file(pspfile, dataset, psp2, revcomp, mtype);

    // warn that we are using the W in the PSP file as minw 
    if (!min_w_set) {
      fprintf(stderr,"Setting minimum motif width to width of the prior in the PSP file: %d\n",
        dataset->psp_w);
      min_w = dataset->psp_w;
    }
  } else { // no PSP file 
    dataset->psp_w = min_w;
  }

  // Set the objective function.
  dataset->objfun = objfun;
  dataset->test = test;
  if (neg_dataset) neg_dataset->objfun = objfun;
  if (neg_dataset) neg_dataset->test = test;

  // create the priors 
  if (ptype == Dmix || ptype == Mega || ptype == MegaP) {
    // make the name of the prior library 
    if (!plib_name) {
      if (alph_is_builtin_protein(alph)) { // default mixture prior for proteins
        plib_name = make_path_to_file(get_meme_data_dir(), PROTEIN_PLIB);
      } else {
        fprintf(
          stderr, 
          "WARNING: When using DNA or a custom alphabet and specifiying a prior type of\n"
          "  'dmix', 'mega' or 'megap', a prior library must be provided.\n"
          "  No prior library was provided, so a simple Dirichlet prior will be used.\n"
        );
        prior = "dirichlet";
        ptype = Dirichlet;
        if (beta <= 0) beta = 0.01; // default b = 0.01 for simple Dirichlet
      }
    }
  }
  if ((ptype == Mega || ptype == MegaP) && beta == -1) {
    // tlb 5-9-97; wgt_total_res 
    //beta = 10.0 * dataset->wgt_total_res; // size of mega prior 
    beta = 5.0 * dataset->wgt_total_res; // size of mega prior 
  }
  priors = create_priors(ptype, beta, dataset, plib_name);

  // set number of occurrences of sites 
  if (nsites != 0) {
    if (mtype == Oops) {
      fprintf(stderr, "You may not specify option '-nsites' with option '-mod oops'.\n");
      exit(1);
    }
    min_nsites = max_nsites = nsites;
  }

  // set search range for nsites 
  if (mtype == Oops) {
    min_nsites = max_nsites = n_primary_samples;
  } else if (mtype == Zoops) {
    if (min_nsites > n_primary_samples) {
      fprintf(stderr, "Minimum number of sites is too large. Setting to 2.\n");
      min_nsites = 2;
    }
    if (!min_nsites) min_nsites = 2; // default 
    if (max_nsites > n_primary_samples) {
      fprintf(stderr, "Maximum number of sites is exceeded. Setting to %d.\n",
        n_primary_samples);
      max_nsites = n_primary_samples;
    }
    if (!max_nsites) max_nsites = n_primary_samples; // default 
  } else { // ANR/TCM model 
    if (min_nsites<2) min_nsites = 2; // default 
    if (!max_nsites) {
      if (objfun == Classic || objfun == NC) {
        //max_nsites = MIN(5*n_primary_samples, 50);
        max_nsites = MIN(5*n_primary_samples, classic_max_nsites);
      } else {
        max_nsites = 5*n_primary_samples;
      }
    }
  }

  // check that max number of sites >= min number of sites 
  if (min_nsites > max_nsites) {
    fprintf(
      stderr, 
      "The minimum number of sites is set to %d. "
      "It should be less than the max number of sites (%d).\n",
      min_nsites,
      max_nsites
    );
    exit(1);
  }
  // check that there are enough possible sites 
  if (min_nsites < 2) {
    fprintf(stderr, "You must specify a minimum of 2 or more sites.\n");
    exit(1);
  }
  if (max_nsites < 2) {
    fprintf(stderr, "It must be possible for at least 2 sites to fit.\n");
    exit(1);
  }

  // check weight on prior on nsites 
  if (wnsites >= 1 || wnsites < 0) {
    fprintf(stderr, "You must specify option '-wnsites' in the range [0..1).\n"); exit(1);
    exit(1);
  }

  // set up globals 
  if (w != 0) { // w specified; set min_w and max_w 
    max_w = min_w = w;
    if (! NO_STATUS) fprintf(stderr,"w set, setting max and min widths to %d\n",w);
  }

  // check that no sequence too short 
  if (dataset->min_slength < min_w) {
    fprintf(stderr,
     "All sequences must be at least %d characters long.  Set '-w' or '-minw' or remove ", min_w);
    fprintf(stderr, "shorter sequences\nand rerun.\n");
    exit(1);
  }

  // oops model: limit max_w to shortest seq 
  if (mtype == Oops && max_w > dataset->min_slength) {
    max_w = dataset->min_slength;
    fprintf(stderr,
      "Option '-maxw' is greater than the length of shortest sequence (%ld).", dataset->min_slength);
    fprintf(stderr, "  Setting '-maxw' to %d.\n", max_w);
  }
  // all models: limit max_w to longest seq 
  if (max_w > dataset->max_slength) {
    max_w = dataset->max_slength;
    fprintf(stderr,
      "Option '-maxw' is greather than the length of longest sequence (%ld).", dataset->max_slength);
    fprintf(stderr, "  Setting '-maxw' to %d.\n", max_w);
  }
  if (max_w > MAXSITE) {
    fprintf(stderr,
      "Option '-maxw' is larger than MAXSITE (%d).  Recompile with larger MAXSITE in user.h.\n", MAXSITE);
    exit(1);
  }
  if (max_w < 0) { // use default 
    max_w = MIN(MAXSITE, dataset->min_slength); // maximum W0 
  }

  // check that min_w <= max_w 
  if (min_w > max_w) {
    if (pspfile) {
      fprintf(stderr, "PSP file w = %d > maxw = %d. Respecify larger '-maxw'.\n",
              dataset->psp_w, max_w);
      exit(1);
    }
     fprintf(stderr, "minw > maxw.  Setting '-minw' to %d.\n", max_w);
     min_w = max_w;
  }

  // check that min_w and max_w are at least ABS_MIN_W 
  if (min_w < ABS_MIN_W) {
    fprintf(stderr,
      "Minimum width must be >= %d.  Respecify larger '-w' or '-minw'.\n",
      ABS_MIN_W);
    exit(1);
  } else if (max_w < ABS_MIN_W) {
    fprintf(stderr,
      "Maximum width must be >= %d.  Respecify larger '-w' or '-maxw'.\n",
      ABS_MIN_W);
    exit(1);
  }

  // must use TCM if only one sequence 
  if (mtype != Tcm && n_primary_samples==1) {
    fprintf(stderr,
      "You must specify '-mod anr' since your dataset contains only one sequence.\n"
    );
    fprintf(stderr,
      "Alternatively, you might wish to break your sequence into several sequences.\n"
    );
    exit(1);
  }

  // flag search for palindromes 
  dataset->pal = pal;
  if (neg_dataset) neg_dataset->pal = pal;

  // check that IUPAC alphabet if using PAM mapping 
  // mapname == "pam" && alphabet != PROTEIN0
  if (strcmp(mapname, "pam") == 0 && !(alph_is_builtin_protein(alph) || alph_is_builtin_dna(alph))) {
    fprintf(stderr,
     "Setting sequence-to-theta mapping type to 'uni' since the alphabet is "
     "not built-in and 'pam' is only supported for the built-in alphabets.\n");
    mapname = "uni";
  }
  // get the type of mapping between sequences and thetas 
  if (strcmp(mapname, "uni") == 0) {
    map_type = Uni;
    //if (map_scale == -1) map_scale = .5; // default add .5 
    int alen = alph_size_core(alph);
    if (map_scale == -1) map_scale = 2.0/alen; // default add 2/alen
  } else if (strcmp(mapname, "pam") == 0) {
    map_type = Pam;
    if (map_scale == -1) map_scale = 120; // default PAM 120 
  } else {
    fprintf(stderr, "Unknown mapping type '-spmap %s'. \n", mapname);
    exit(1);
  }

  // set up the sequence to theta mapping matrix for starts 
  dataset->map = init_map(map_type, map_scale, alph, dataset->back, false);
  dataset->lomap = init_map(map_type, map_scale, alph, dataset->back, true);

  // set up p_point with space for starting points specified by -cons
  int ncons = get_num_strings(spcons); 		// number of specified starts 
  p_point = (P_POINT *) mm_malloc(sizeof(P_POINT));
  p_point->c = ncons;
  if (ncons > 0) {
    p_point->w = NULL; Resize(p_point->w, ncons, int);
    p_point->nsites = NULL; Resize(p_point->nsites, ncons, double);
    p_point->e_cons0 = NULL; Resize(p_point->e_cons0, ncons, uint8_t*);
  }

  // setup user-specified start points 
  for (i = 0; i < p_point->c; i++) {
    char *cons = get_nth_string(i, spcons);
    int len = strlen(cons);					// length of starting point
    uint8_t *e_cons = NULL; Resize(e_cons, len, uint8_t);
    p_point->e_cons0[i] = e_cons;
    // encode as integer 
    for (j = 0; (cc = cons[j]) != '\0'; j++) {
      if (alph_is_known(alph, cc)) {
        // encode to core + wildcard
        e_cons[j] = alph_encodec(alph, cc);
      } else {
        fprintf(stderr, "Illegal letter %c in consensus string '-cons %s'.\n", cc, cons);
        exit(1);
      }
    }
    // set width to length of consensus
    p_point->w[i] = j;
    if (p_point->w[i] > max_w) {
      fprintf(stderr, "'-cons %s' requires -w or -maxw to be at least %d.\n", cons, p_point->w[i]);
      exit(1);
    }
  }

  // Setup heap size for storage of starting points:
  if (main_hs < 1) {
    fprintf(stderr, "Option '-heapsize' must be >= 1.\n");
    exit(1);
  } else {
    dataset->main_hs = main_hs;
    dataset->hs_decrease = hs_decrease;
  }

  // Setup a struct recording the desired parameters for branching search:
  BRANCH_PARAMS *branch_params = NULL;
  branch_params = mm_malloc(sizeof(BRANCH_PARAMS));

  // Setup branching factor for branching search:
  if (bfactor < 1) {
    fprintf(stderr, "Option '-bfactor' must be >= 1.\n");
    exit(1);
  } else {
    branch_params->bfactor = bfactor;
  }

  // Record whether the user wants w-branching carried out:
  branch_params->w_branch = w_branch;

  // Record whether the user wants x-branching carried out:
  if (x_branch){
    branch_params->point_branch = X_ONLY;
  } else { 
    branch_params->point_branch = NO_POINT_B;
  }

  if (false) {
    // This code can be used to set branching as the default 
    // for different alphabets/sequence models using the -x_branch
    // and -no_x_branch switches. Currently the default is set to
    // no branching (in previous if statement).
    /* Record whether the user wants x_branching carried out. The result
       is recorded in the "point_branch" attribute of branch_params. Note
       that the user is only able to specify x_branching or no x_branching.
       This is because we found ACGT (ie regular) branching to be of little
       benefit. */
    if (x_branch) {
      // User has specified they want x-branching... 
      if (no_x_branch) {
	// Invalid to specify x_branch and no_branch simultaneously:
	fprintf(stderr, "Options '-x_branch' and '-no_branch' cannot be specified"\
			"at the same time");
	exit(1);
      }
      else {
	branch_params->point_branch = X_ONLY;
      }
    } else if (no_x_branch) {
      branch_params->point_branch = NO_POINT_B;
    } else {
      // User did not specify x_branching => Decide based on sequence model:
      // NOTE: branching is the default for oops for DNA only
      if (mtype == Oops && alphabet_type == Dna) {
	branch_params->point_branch = X_ONLY; // Only x_branch under oops by default.
      } else {
	branch_params->point_branch = NO_POINT_B; // Only x_branch under oops by default.
      }
    } // Deciding branch_params->x_branch
  } // end if 

  // Store the branching parameters for future reference:
  dataset->branch_params = branch_params;

  // Print sites predicted by MEME:
  dataset->print_pred = print_pred;

  // Record whether heaps are to be printed:
  dataset->print_heaps = print_heaps;

  // Omit large tables of sites and sequences in .txt output?
  dataset->brief = brief;

  // make sure nmotifs is as large as the number of starting points 
  if (nmotifs < p_point->c) {
    nmotifs = p_point->c;
    fprintf(stderr, "Setting '-nmotifs' to %d.\n", nmotifs);
  }

  // prevent too long jobs 
  dataset->search_size = search_size;
  dataset->max_words = max_words;
  dataset->max_size = max_size;
  dataset->no_rand = norand;
  dataset->classic_max_nsites = (objfun == Classic || objfun == NC) ? classic_max_nsites : -1;

  // Set the high-water mark for restricting the total number of seeds tested.
  // Note: Does not include holdout sets in DE and SE modes.
  // Note: uses min_w as requested by user which may get reset below.
  int n_seed_seqs = (objfun==DE || objfun==SE) ? dataset->n_group[0] : dataset->n_samples;
  dataset->last_seed_seq = n_seed_seqs - 1;
  dataset->last_seed_pos = dataset->samples[n_seed_seqs-1]->length - min_w;
  if (dataset->total_res > max_words) {
    // Find where the last seed word of the minimum width would start.
    int total_res = 0;
    for (i=0; i<n_seed_seqs; i++) {
      int slen = dataset->samples[i]->length;
      if (total_res + slen - min_w + 1 >= max_words) {
        int start_of_central_region = (dataset->samples[i]->length - slen)/2;
        dataset->last_seed_seq = i;
        dataset->last_seed_pos = MAX(0, start_of_central_region + (max_words - total_res - min_w));
        break;
      }
      total_res += slen-min_w+1;
    }
    n_seed_seqs = dataset->last_seed_seq + 1;
  }
  if (!NO_STATUS) fprintf(stderr, "SEEDS: maxwords %d highwater mark: seq %d pos %d\n", max_words, n_seed_seqs, dataset->last_seed_pos);

  // Load balance the search for starting points in parallel mode.
  // Note: This must come after adding in control samples.
  balance_loop(dataset->samples, n_seed_seqs);

  // create the model 
  model = create_model(mtype, revcomp, max_w, alph, objfun);
  best_model = create_model(mtype, revcomp, max_w, alph, objfun);

  // initialize log and exp lookup tables 
  init_log();
  init_exp();

  // Initialize the probability tables for the objective function.  
  // Get for 2 and up because weighting may cause < min_nsites sites 
  if (objfun == Classic) {
    // Initialize the probability tables for the objective function.
    nsites = (mtype == Oops || mtype == Zoops) ? MIN(dataset->n_group[0], max_nsites) : max_nsites;
    init_llr_pv_tables(2, nsites, alph_size_core(alph), dataset->back, dataset->pal);
  }

  // set up scratch models 
  model->pal = pal; 
  model->min_w = min_w;
  model->max_w = max_w;
  model->all_widths = all_widths;
  copy_model(model, best_model, alph);

  // put meme parameters in dataset 
  dataset->use_llr = use_llr;
  dataset->priors = priors;
  dataset->p_point = p_point;
  dataset->wg = wg;
  dataset->ws = ws;
  dataset->endgaps = endgaps;
  dataset->min_nsites = min_nsites;
  dataset->max_nsites = max_nsites;
  dataset->wnsites = wnsites;
  dataset->ma_adj = ma_trim;
  dataset->distance = distance;
  dataset->nmotifs = nmotifs;
  dataset->maxiter = maxiter;
  dataset->evt = evt;
  dataset->mod = mod;
  dataset->mapname = mapname;
  dataset->map_scale = map_scale;
  dataset->priorname = prior;
  dataset->beta = beta;
  dataset->seed = seed;
  dataset->hsfrac = hsfrac;
  dataset->shuffle = kmer;
  dataset->motifs = NULL;
  // save name of psp file
  if (pspfile) {
    char *tmp, *base; // remove the directory, leaving just the file name
    for (tmp = base = pspfile; *tmp; tmp++) if (*tmp == '/') base = tmp + 1;
    dataset->pspfile = strdup(base); // copy the name
  } else {
    dataset->pspfile = NULL;
  }
  // save name of prior library 
  if (plib_name) {
    char *tmp, *base; // remove the directory, leaving just the file name
    for (tmp = base = plib_name; *tmp; tmp++) if (*tmp == '/') base = tmp + 1;
    dataset->plib_name = strdup(base); // copy the name
  } else {
    dataset->plib_name = NULL;
  }
  dataset->datafile = sf ? sf : datafile; // name to print 
  if (kmer == -1) {
    dataset->negfile = negfile ? negfile : "--none--";
    if (neg_dataset) neg_dataset->datafile = negfile;
  } else {
    char *tmp = NULL; 
    Resize(tmp, 60, char);
    snprintf(tmp, 60, "Primary sequences shuffled preserving %d-mers", kmer);
    dataset->negfile = tmp;
    if (neg_dataset) neg_dataset->datafile = tmp;
  }
  dataset->neg_n_samples = neg_n_samples;
  dataset->bfile = bfile;
  dataset->max_time = max_time;
  memcpy(&(dataset->primary_groups), &primary_groups, sizeof(GROUP_T));
  memcpy(&(dataset->control_groups), &control_groups, sizeof(GROUP_T));

  // save command line 
  dataset->command = NULL;
  argv[0] = "meme";
  /*for (i=pos=len=0; i<argc-2; i++) {*/ // don't save last arguments 
  // save all arguments; used to not save last 2; why??? 
  for (i=pos=len=0; i<argc; i++) {
    len += strlen(argv[i])+1; // +1 for space following 
    Resize(dataset->command, len+2, char); // +1 for null 
    strcpy((dataset->command)+pos, argv[i]);
    dataset->command[len-1] = ' ';
    dataset->command[len] = '\0';
    pos = len;
  }
  if (TEXT_ONLY) {
    dataset->output_directory = NULL;
  }
  else {
    dataset->output_directory = *output_dirname;
  }

  dataset->mpi = mpi;

  // set up return values 
  *model_p = model;
  *best_model_p = best_model;
  *neg_model_p = neg_model;
  *dataset_p = dataset;
  *neg_dataset_p = neg_dataset;

  // announce meme 
  char *program = "MEME - Motif discovery tool";
  char *info = "For further information on how to interpret these results please access " SITE_URL ".\n"
    "To get a copy of the MEME Suite software please access " SOURCE_URL ".\n";
  banner(program, info, MEME_CITE, *text_output);

  // cleanup 
  free(plib_name);

} // init_meme 
