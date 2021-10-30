/************************************************************************
*                                                                       *
*       MEME                                                            *
*       Author: Timothy L. Bailey                                       *
*                                                                       *
*       Copyright                                                       *
*       (1994 - 2018) The Regents of the University of California.      *
*       All Rights Reserved.                                            *
*                                                                       *
*       Permission to use, copy, modify, and distribute any part of     *
*       this software for educational, research and non-profit purposes,*
*       without fee, and without a written agreement is hereby granted, *
*       provided that the above copyright notice, this paragraph and    *
*       the following three paragraphs appear in all copies.            *
*                                                                       *
*       Those desiring to incorporate this software into commercial     *
*       products or use for commercial purposes should contact the      *
*       Technology Transfer Office, University of California, San Diego,*
*       9500 Gilman Drive, La Jolla, California, 92093-0910,            *
*       Ph: (619) 534 5815.                                             *
*                                                                       *
*       IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO     *
*       ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR         *
*       CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF   *
*       THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF CALIFORNIA  *
*       HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             *
*                                                                       *
*       THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE *
*       UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE          *
*       MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  *
*       THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND       *
*       EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, *
*       INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF        *
*       MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT    *
*       THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,           *
*       TRADEMARK OR OTHER RIGHTS.                                      *
*                                                                       *
************************************************************************/

 // 5-26-00 tlb; initialize model->cons0 in init_model 
 // 10-3-99 tlb; replace nrdb with back 
 // 7-29-99 tlb; set model->pal correctly in init_model 
 // 7-02-99 tlb; in erase normalize z_ij to 1.0 so erasing will be complete 
 // 7-02-99 tlb; in erase run just e_step of em 
 // 6-28-99 tlb; remove sigma from model 
 // 8-7-97 tlb; removed Mcm stuff from erase 

/*
        meme <datafile> [options]

        Reads in a sequence file (Pearson/FASTA format).
        Runs EM algorithm with selected values of W and lambda.

        <datafile>      set of samples: [>id sequence]+
*/

#ifdef MEME_MAIN
#define DEFINE_GLOBALS
#endif

#include <sys/types.h>
#include <sys/wait.h>

#include "dir.h"
#include "meme.h"
#include "calculate_p_y.h"
#include "seed_diffs.h"
#include "sp_matrix.h"
#include "utils.h"
#include "xml-util.h"
#include "string-builder.h"

#ifndef PARALLEL
#define mpMyID() 0
#endif

// external subroutines 
extern double sqrt(double x);

// local subroutines 
static bool save_candidate(
  MODEL *model, // final model 
  DATASET *dataset, // the dataset 
  S_POINT *s_point, // starting point 
  CANDIDATE *candidates, // candidate for best model of width 
  double best_sig, // best significance 
  int best_w, // best width
  double best_llr // best motif LLR
);

/*
 * Utility function for multi-process debugging.
 * When called it prints out the process ID and host and then waits for you to
 * attach the debugger. Once you have attached the debugger change the frame
 * to this function and run "set var i = 1" (or any number other than 0), then
 * "continue" to allow the program to run.
 */
static void await_debugger() {
  int i = 0;
  char hostname[256];
  char message[500];
  gethostname(hostname, sizeof(hostname));
  // work around the shadowing of printf and fprintf by macros
  snprintf(message, sizeof(message), "PID %d on %s ready for attach\n", getpid(), hostname);
  fputs(message, stdout);
  fflush(stdout);
  while (0 == i) sleep(5);
}

bool no_print = false; // turn off printing if parallel and not main 

// for getting rid of gcc's annoying warn_unused_result
static inline void ignore_value (int i) { (void) i; }
static inline void ignore_ptr (void* p) { (void) p; } 
/**********************************************************************/
/*
        main
*/
/**********************************************************************/
int meme_main(
  int argc,
  char *argv[],
  ARRAYLST_T* seq_array,
  int width,
  bool eliminate_repeats
)
{
  int i, imotif;
  DATASET *dataset; // the dataset 
  DATASET *neg_dataset; // the dataset of negative examples 
  MODEL *model; // the model 
  MODEL *neg_model; // the model of negatives 
  MODEL *best_model; // the best model found (so far) 
  CANDIDATE *candidates = NULL; // list of candidate models 
  MOTIF_SUMMARY *motif_summaries = NULL; // properties of final motifs 
  int n_starts = 0; // number of starting points 
  int nmotifs; // number of to find 
  double stop_time=0, m1time=0; // time when stopped, for motif 1 (secs) 
  int last_w; // keep track of change for PSP 
  bool output_error = false;

  char * text_filename = "meme.txt";
  char * xml_filename = "meme.xml";
  char * html_filename = "meme.html";
  FILE *text_output = NULL; // File for text output 

  int min_w, max_w;

#ifdef PARALLEL
  int start_start, incr;
  exec_parallel(argc, argv);     // Relaunch using MPI if -p <np> given.
  // Initialize MPI. 
  mpInit(&argc, &argv);
  // turn off printing if parallel and not the main processor 
  no_print = (mpMyID() != 0);
  incr = mpNodes();
  //await_debugger();
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "Running on %d nodes...\n", incr);
#endif 
#endif // PARALLEL 

#ifdef debug_ieee
  ieee_handler_setup("common");
#endif

  (void) mytime(0); 	// start elapsed time clock

  // initialize everything from the command line 
  char *output_dirname = NULL;
  init_meme(
    argc, argv, &model, &best_model, &neg_model,
    &dataset, &neg_dataset, text_filename, &output_dirname, &text_output, seq_array, width, eliminate_repeats
  );
  // Save the minimum and maximum widths.
  min_w = model->min_w;
  max_w = model->max_w;
  if (model->all_widths && !NO_STATUS)
    fprintf(stderr, "all widths from min to max\n");
  nmotifs = dataset->nmotifs; // number of motifs to find 
  Resize(candidates, model->max_w+2, CANDIDATE); // motif candidates; max_w+2 in case -pal given 
  Resize(motif_summaries, nmotifs, MOTIF_SUMMARY);

  char *hostname;
  hostname = HOSTNAME;
  if (hostname[0] == '\0') hostname = "unknown";

  if (!no_print && VERBOSE) fprintf(text_output, "\nCPU: %s\n\n", hostname);
  // print the command line 
  if (VERBOSE) {
    argv[0] = "meme";
    for (i=0; i<argc; i++) fprintf(text_output, "%s ", argv[i]);
    fprintf(text_output, "\n\n");
    fflush(text_output);
  }
  // describe the dataset 
  print_dataset_summary (dataset, text_output);

  // describe the command 
  print_command_summary(model, dataset, text_output);

  if (!NO_STATUS) {
    fprintf(stderr, "\nseqs=%6d, min_w=%4ld, max_w=%5ld, total_size=%9d\n",
      dataset->n_samples, dataset->min_slength, dataset->max_slength,
      dataset->total_res);
  }

  // Find a concept and erase it loop 
  OBJTYPE objfun = dataset->objfun;
  int ncons = dataset->p_point->c;			// number of consensus sequences provided by user
  for (imotif=1; imotif<=nmotifs; imotif++) {
    S_POINT *s_points = NULL; 		// array of starting points 
    int i_start; 			// index in s_points 
    uint8_t *e_cons = NULL;
    int best_w = 0; 			// best width 
    double best_sig; 			// motif significance 
    double best_llr;                    /* motif LLR */
    int iter = 0; 			// total number of EM iterations 
    model->min_w = min_w;		// Reset in case they were changed below.
    model->max_w = max_w;

    // Initialize from consensus?
    if (ncons >= imotif) {
      e_cons = dataset->p_point->e_cons0[imotif-1]; 	// consensus sequence 
      model->min_w = model->max_w = dataset->p_point->w[imotif-1];	
    }

    // 29-05-06: Also print motif number if printing of sites is requested:
    if (!NO_STATUS) {
      fprintf(stderr, "\nmotif=%d\n", imotif);
    }

    if (dataset->print_pred) {
      dataset->imotif = imotif;
    }

    // set up the array of starting points for EM 
    s_points = get_starts(dataset, neg_dataset, model, e_cons, &n_starts);

    // leave loop if no starts found 
    if (n_starts == 0) break;

    // tag each candidate width as unused; max_w+1 in case -pal given 
    for (i=0; i<=model->max_w+1; i++) {
      candidates[i].sig = BIG;
      candidates[i].s_point = NULL;
    }

    // run EM on each start and save best final model for each final width 
    best_sig = BIG; // best motif significance 
    best_llr = -BIG; // LLR of best motif

    // psp_w set to min_w if no PSP file
    last_w = dataset->psp_w;
#ifdef PARALLEL
    // Check whether to parallelize this loop. 
    start_start = mpMyID();
    // Make sure everybody has something to do. 
    if (start_start >= n_starts) start_start = n_starts-1;
    // Divide the various starting points among processors. 
    for (i_start=start_start; i_start < n_starts; i_start += incr) {
#else
    for (i_start=0; i_start<n_starts; i_start++) {
#endif // PARALLEL 
      S_POINT *s_point = s_points+i_start; // current starting point 
#ifdef DEBUG
      double s_point_time = mytime(0)/1E6;
#endif // DEBUG 

      // initialize the model from the starting point 
      if (! init_model(s_point, model, dataset, imotif)) continue;

      // Count iters per loop. 
      model->iter = 0;

      // don't optimize width if -cons given
      if (ncons >= imotif) model->min_w = model->max_w = dataset->p_point->w[imotif-1];	

      // Run EM from the starting model 
      em(model, dataset, neg_dataset, s_point);

      // Keep track of the total number of EM iterations. 
      iter += model->iter;

      // Adjust the p-value by the number of tests and store in E-value of model.
      if (objfun != Classic && objfun != NC) {
        model->logev = model->logpv + log(n_starts);		// n * p
        model->logpv = LOGEV(log(n_starts), model->logpv);	// 1-(1-p)**n
      }

      // If requested, print the final MEME site predictions (for the "best"
      // starting point after EM has been completed)...
      // Retrieve the array of sites predicted by the model:
      if (dataset->print_pred) {
        P_PROB pred_sites = model->maxima;
        int n_psites = model->nsites_dis;
        fprintf(stdout,
	  "\nPREDICTED SITES AFTER EM FROM STARTING POINT WITH W = %i AND"
	  " NSITES = %f:\n", s_point->w0, s_point->nsites0); //model->w, n_psites);
        print_site_array(pred_sites, n_psites, stdout, model->w, dataset);
        double sig = model->logev;
        fprintf(stdout, "MODEL SIGNIFICANCE = %f\n", sig);
      }

      /* store model as a candidate for best model;
         save model if best so far */
      if (save_candidate(model, dataset, s_point, candidates, best_sig, best_w, best_llr)) {
        best_w = model->w; // new best width 
        SWAP(MODEL*, model, best_model); // save model as best_model 
        best_sig = candidates[best_w].sig; // new best significance 
        best_llr = candidates[best_w].llr; // new best LLR
      }
    } // starting point loop 

    // found the best model 
    SWAP(MODEL*, model, best_model);

#ifdef PARALLEL
    // Copy the consensus sequence into the model. 
    store_consensus(model, candidates);

    // Do the reduction. 
    reduce_across_models(model, model->alength);
#ifdef DEBUG_PARALLEL
    fprintf(stderr, "%d: past reduction\n", mpMyID()); fflush(stderr);
#endif 
#endif // PARALLEL 

    // quit if model has too few sites 
    if (model->nsites_dis < MINSITES) break;
#ifdef DEBUG_PARALLEL
    fprintf(stderr, "%d: past few sites\n", mpMyID()); fflush(stderr);
#endif 

    /* quit if model fails E-value test */
    if (model->logev > log(dataset->evt)) break;
#ifdef DEBUG_PARALLEL
    fprintf(stderr, "%d: past E-value\n", mpMyID()); fflush(stderr);
#endif 

    // Store the total number of EM's in the model. 
    model->iter = iter;

    // ERASE the site and starts 
#ifdef DEBUG_PARALLEL
    fprintf(stderr, "%d: at erase\n", mpMyID()); fflush(stderr);
#endif 
    erase(dataset, model);
#ifdef DEBUG_PARALLEL
    fprintf(stderr, "%d: past erase\n", mpMyID()); fflush(stderr);
#endif 

    // Translate sequence indices in maxima to original order
    // so that output will preserve that.
    for (i=0; i<model->nsites_dis; i++) {
      model->maxima[i].x = dataset->samples[model->maxima[i].x]->orig_index;
    }

    // print results 
#ifdef DEBUG_PARALLEL
    fprintf(stderr, "%d: at print results\n", mpMyID()); fflush(stderr);
#endif 
    if (!no_print) {
      print_results(
        dataset,
        model,
        candidates,
        text_output
      );
    }

#ifdef DEBUG_PARALLEL
    fprintf(stderr, "%d: past print results\n", mpMyID()); fflush(stderr);
#endif 
    if (!no_print) record_results(dataset, model, motif_summaries);

    // stop if out of time 
    if (dataset->max_time && imotif<nmotifs) { // considering another motif 
      stop_time = mytime(0)/1E6; // current time 
#ifdef PARALLEL
      // use stop_time from process 0 to avoid it stopping while others continue
#ifdef DEBUG_PARALLEL
      fprintf(stderr, "%d: time broadcast size %lu\n", mpMyID(),sizeof(stop_time)); fflush(stderr);
#endif 
      mpBroadcast((void *)&stop_time, sizeof(stop_time), 0);
#ifdef DEBUG_PARALLEL
      fprintf(stderr, "%d: past time broadcast\n", mpMyID()); fflush(stderr); 
#endif 
#endif
      // record time if this is motif 1 
      if (imotif == 1) m1time = stop_time;
      if ((dataset->max_time - stop_time) < m1time) {
        ++imotif; // this motif OK 
        break;
      }
    } // check time 

    // The s_points for the current motif are no longer needed => delete them:
    int sp_idx;
    for (sp_idx = 0; sp_idx < n_starts; sp_idx++) {
      free_s_point(&(s_points[sp_idx]));
    }
    free(s_points);
    
#ifdef DEBUG_PARALLEL
    fprintf(stderr, "%d: past check time\n", mpMyID()); fflush(stderr);
#endif 
  } // nmotifs loop 
  --imotif; // number of motifs found 

  // print the motif block diagrams using all the motifs 
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: At print summary\n", mpMyID()); fflush(stderr);
#endif
  if (!no_print && imotif && (dataset->n_samples <= dataset->brief))
    print_summary(model, dataset, imotif, text_output);

  // Record the reason for stopping
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: At print stop reason\n", mpMyID()); fflush(stderr);
#endif
  #define MAX_REASON_LENGTH 255
  char stopping_reason[MAX_REASON_LENGTH];
  if (n_starts == 0) {
    snprintf(
      stopping_reason,
      MAX_REASON_LENGTH,
      "Stopped because couldn't find any more starting points for EM."
    );
  } else if (imotif == nmotifs) {
    snprintf(
      stopping_reason,
      MAX_REASON_LENGTH,
      "Stopped because requested number of motifs (%d) found.",
      nmotifs
    );
  } else if (dataset->max_time && (dataset->max_time - stop_time) < m1time) {
    snprintf(
      stopping_reason,
      MAX_REASON_LENGTH,
      "Stopped because MEME would probably exceed the allowed time (%.2f secs) before\nfinding the next motif.",
      dataset->max_time
    );
  } else if (model->nsites_dis < MINSITES) {
    snprintf(
      stopping_reason,
      MAX_REASON_LENGTH,
      "Stopped because next motif has fewer than %d sites.",
      MINSITES
    );
  } else {
    snprintf(
      stopping_reason,
      MAX_REASON_LENGTH,
      "Stopped because the next motif E-value > %8.2e.",
      dataset->evt
    );
  }
  // print reason for stopping 
  fprintf(text_output, "\n");
  PSTARS(text_output);
  fprintf(text_output, "%s\n", stopping_reason);
  PSTARS(text_output);
  fflush(text_output);

#ifdef UNIX
  if (!no_print) {
    #define MAX_HOSTNAME 255
    char hostname[MAX_HOSTNAME];
    int result = gethostname(hostname, MAX_HOSTNAME);
    if (result < 0) {
      // gethostname failed, but we don't really care why
      fprintf(text_output, "\nCPU: unknown\n\n");
    }
    else {
      fprintf(text_output, "\nCPU: %s\n\n", hostname);
    }
  }
#else
  if (!no_print) fprintf(text_output, "\nCPU: unknown\n\n");
#endif // UNIX 
  PSTARS(text_output);

  if (!NO_STATUS) fprintf(stderr, "\n");

  // Print results in XML and HTML format 
  if (!no_print && TEXT_ONLY == false) {
#ifdef DEBUG_PARALLEL
    fprintf(stderr, "%d: At make_path_to_file\n", mpMyID()); fflush(stderr);
#endif
    char * xml_path = make_path_to_file(output_dirname, xml_filename);
#ifdef DEBUG_PARALLEL
    fprintf(stderr, "%d: At print_meme_file_xml\n", mpMyID()); fflush(stderr);
#endif

    // Remove control samples from the primary dataset if needed.
    if (objfun==NZ) dataset->n_samples = dataset->n_group[0];

    print_meme_file_xml(
      model,
      dataset,
      neg_dataset,
      imotif,
      motif_summaries,
      stopping_reason,
      xml_path
    );
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: At free xml_path\n", mpMyID()); fflush(stderr);
#endif
    myfree(xml_path);

    // Free arrays in motif_summaries
    for (i = 0; i < imotif; i++) {
      free_2array(motif_summaries[i].pssm, motif_summaries[i].width + 1);
      free_2array(motif_summaries[i].psfm, motif_summaries[i].width + 1);
      free(motif_summaries[i].sites);
      free(motif_summaries[i].regexp);
      free(motif_summaries[i].consensus);
    }
    free(motif_summaries);

    // generate HTML
    char *prog;
    prog = get_meme_libexec_file("meme_xml_to_html");
    if (prog != NULL) {
      STR_T *cmd;
      int ret;
      cmd = str_create(0);
      str_append2(cmd, prog);
      str_append(cmd, " ", 1);
      str_append_path(cmd, 2, output_dirname, xml_filename);
      str_append(cmd, " ", 1);
      str_append_path(cmd, 2, output_dirname, html_filename);

      ret = system(str_internal(cmd));

      if (!(WIFEXITED(ret) && WEXITSTATUS(ret) == 0)) {
        output_error = true;
        fprintf(stderr, "Warning: meme_xml_to_html exited abnormally and may "
            "have failed to create HTML output.\n");
      }

      str_destroy(cmd, false);
      free(prog);
    } else {
      output_error = true;
      fprintf(stderr, "Warning: could not find meme_xml_to_html. "
          "The HTML output could not be created.\n");
    }
  } // !no_print

  // save mpi flag
  bool mpi = dataset->mpi;

  free_model(model);  
  free_model(neg_model);
  free_model(best_model);
  free(candidates);
  free_2array(dataset->map, alph_size_wild(dataset->alph));
  free_2array(dataset->lomap, alph_size_wild(dataset->alph));
  free(dataset->command);
  if (dataset->p_point->c > 0) {
    for (i = 0; i < dataset->p_point->c; i++) {
      myfree(dataset->p_point->e_cons0[i]);
    }
    myfree(dataset->p_point->w);
    myfree(dataset->p_point->nsites);
    myfree(dataset->p_point->e_cons0);
  }
  free(dataset->p_point);
  free_PriorLib(dataset->priors->plib);
  free(dataset->priors);
  free_array(dataset->back);
  free(dataset->branch_params);
  free(dataset->res_freq);
  for (i = 0; i < dataset->n_samples; i++) free_sample(dataset->samples[i]);
  free(dataset->samples);
  if (dataset->input_order != dataset->samples) free(dataset->input_order);
  free(dataset);

#ifdef PARALLEL
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: At mpFinalize\n", mpMyID()); fflush(stderr);
#endif
  if (mpi) mpFinalize();
#endif

  return (output_error ? 1 : 0);
} // meme_main

#ifdef MEME_MAIN
int main(
         int argc,
         char *argv[]
         )
{
  return meme_main(argc, argv, NULL, 0, false);
}
#endif

/**********************************************************************/
/*
        save_candidate

        Save the starting point and part of model if it is
        most significant model of its width so far:
                model->sig is smallest

        Returns true if the model is the best so far among all widths.
*/
/**********************************************************************/
static bool save_candidate(
  MODEL *model, // final model 
  DATASET *dataset, // the dataset 
  S_POINT *s_point, // starting point 
  CANDIDATE *candidates, // candidate for best model of width 
  double best_sig, // best motif significance 
  int best_w, // best motif width
  double best_llr // best motif LLR
)
{
  int w = model->w; // final motif w 
  bool classic = (dataset->objfun==Classic); // Classic meme

  /* objective function value */
  double sig = model->logev;
  double llr = model->llr;
  double m1, e1;                    // for printing significance 

  // get E-value of motif
  exp10_logx(sig/log(10.0), m1, e1, 1);

  /* print the results for this w0, nsites0 and THETA */
  if (PRINT_STARTS) {
#if defined(PARALLEL) && defined(DEBUG_PARALLEL)
    fprintf(stdout, "\nnode %d (start) %3d %6.1f %.*s --> %s ", mpMyID(),
#else
    printf("\n(start) %3d %6.1f %.*s --> %s ",
#endif
      s_point->w0, s_point->nsites0, s_point->w0, s_point->cons0, model->cons);
#if defined(PARALLEL) && defined(DEBUG_PARALLEL)
    fprintf(stdout, "w %3d nsites %4d sig %20.10g --> %3.1fe%+04.0f \n\n", 
      //w, model->nsites_dis, s_point->score, exp(sig));
      w, model->nsites_dis, s_point->score, m1, e1);
#else
    printf("w %3d nsites %4d sig %20.10g --> %3.1fe%+04.0f\n\n", 
      //w, model->nsites_dis, s_point->score, exp(sig));
      w, model->nsites_dis, s_point->score, m1, e1);
#endif
    fflush(stdout);
  }

  // save the results if best so far for this width.
  // Break ties on significance using LLR except in classic mode.
  if (
    (sig < candidates[w].sig)
    || (!classic && sig == candidates[w].sig && llr > candidates[w].llr) ) {
    candidates[w].s_point = s_point;
    candidates[w].w = w;
    candidates[w].pal = model->pal;
    candidates[w].invcomp = model->invcomp;
    candidates[w].lambda = model->lambda;
    strcpy(candidates[w].cons, model->cons);
    candidates[w].ic = model->ic;
    candidates[w].rel = model->rel;
    candidates[w].ll = model->ll;
    candidates[w].sig = sig;
    candidates[w].llr = llr;
  }

  // Return true if this model has the best significance,
  // breaking ties with 1) width and 2) LLR.
  // This rule must match that in max_packets().
  return (
    (sig < best_sig) 
    || (!classic 
      && (
	(sig == best_sig && w < best_w)
	|| (sig == best_sig && w == best_w && llr > best_llr)
      )
    )
  );
} // save candidates 
