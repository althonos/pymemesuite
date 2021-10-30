/* GLAM2 */
#include <assert.h>
#include <float.h>  /* DBL_MIN, DBL_EPSILON */
#include <sys/wait.h> /* WEXITSTATUS and associated macros */
#include "glam2_util.h"
#include "glam2_glam2.h"
#include "glam2_init.h"
#include "glam2_output.h"
#include "glam2_site_sample.h"
#include "glam2_column_sample.h"
#include "utils.h"

VERBOSE_T verbosity;		// needed by meme utilities

/* log probability of counts given probabilities */
double prob_vec_score(const prob_vec *p, const int *counts) {
  int i;
  double score = 0;
  for (i = 0; i < p->dim; ++i)
    score += counts[i] * p->log_probs[i];
  return score;
}

double bg_score(const glam2_scorer *s, const int *counts) {
  return prob_vec_score(&s->bg, counts);
}

double emission_score(glam2_scorer *s, const int *counts) {
  return dmix_score(&s->e_prior, counts);
}

double deletion_score(glam2_scorer *s, int delete_count, int no_delete_count) {
  return beta_score(&s->d_prior, delete_count, no_delete_count);
}

double insertion_score(glam2_scorer *s, int insert_count, int no_insert_count) {
  return beta_score(&s->i_prior, insert_count, no_insert_count);
}

double column_score(glam2_scorer *s, const glam2_col *col) {
  return emission_score(s, col->emission_counts)
    + deletion_score(s, col->delete_count, col->match_count)
    - bg_score(s, col->emission_counts);
}

/* Calculate the score of an alignment */
double aln_score(glam2_scorer *s, const glam2_aln *aln) {
  double score = 0;
  int i;

  for (i = 0; i < aln->width; ++i)
    score += column_score(s, &aln->cols[i]);

  for (i = 1; i < aln->width; ++i)
    score += insertion_score(s, aln->insert_counts[i-1], aln->aligned_seq);

  return score;
}

/* Calculate one sequence's contribution to the alignment score */
double marginal_score(glam2_scorer *s, glam2_aln *aln,
		      int seq, const fasta *f) {
  double score = aln->score;
  unalign(aln, seq, f);
  score -= aln_score(s, aln);
  realign(aln, seq, f);
  return score;
}

/* get a random starting alignment */
void start_aln(glam2_aln *aln, data *d) {
  int i;
#if 0
  aln->width = d->a.min_width;  /* ?? initial number of columns */
  aln->width = sqrt(d->a.max_width * d->a.min_width);  /* geometric mean */
#endif
  aln->width = d->a.init_width;
  aln_zero(aln);
  SHUFFLE(d->seq_order, aln->seq_num);
  for (i = 0; i < aln->seq_num; ++i)
    site_sample(aln, d->seq_order[i], d, 1);
  aln->score = aln_score(&d->scorer, aln);
}

void update_aln(glam2_aln *aln, data *d, const double temperature) {
  assert(aln->seq_num > 0);
  if (rand_dbl(d->a.column_sample_rate + 1) < 1) {
    const int seq_pick = rand_int(aln->seq_num);
    if (d->a.profile)
      fprintf(d->out, "site sample: seq=%d\n", seq_pick);
    site_sample(aln, seq_pick, d, temperature);
  } else {
    column_sample(aln, d, temperature);
  }
  aln->score = aln_score(&d->scorer, aln);
}

void optimise_aln(glam2_aln *best, data *d) {
  glam2_aln *aln = &d->aln;
  int no_improvement = 0;
  int i;

  aln_copy(aln, best);
  if (d->a.profile)
    fputs("Temperature, Columns, Sequences, Score:\n", d->out);

  for (i = 0; no_improvement < d->a.stop_after; ++i) {
    double temperature = d->a.temperature /
      xpow(d->a.cool, (double)i / d->a.stop_after);
    if (temperature < d->a.frozen)
      temperature = d->a.frozen;
    if (d->a.profile)
      fprintf(d->out, "%g\t%d\t%d\t%g\n",
	      temperature, aln->width, aln->aligned_seq, aln->score / xlog(2));
    /*
    print_aln(d->out, aln, d);
    */
    update_aln(aln, d, temperature);
    if (aln->score > best->score) {
      aln_copy(best, aln);
      no_improvement = 0;
    } else
      ++no_improvement;
  }

  if (d->a.profile)
    putc('\n', d->out);
  if (!d->a.quiet) fprintf(stderr, "%d iterations\n", i);
}

void print_misc_info(FILE *fp, const data *d) {
  const int alph_size = d->alph.size;
  int i;
  fprintf(fp, "Sequences: %d\n", d->seqs.seqnum);
  fprintf(fp, "Greatest sequence length: %d\n", d->seqs.maxlen);
  fputs("Residue counts: ", fp);
  for (i = 0; i <= alph_size; ++i)
    fprintf(fp, "%c=%d%c", d->alph.decode[i], d->scorer.bg.counts[i],
	    i < alph_size ? ' ' : '\n');
}

/* Alignment comparison function for sorting */
int aln_cmp(const void *a, const void *b) {
  const double x = ((const glam2_aln *)a)->score;
  const double y = ((const glam2_aln *)b)->score;
  return x < y ? +1 : x > y ? -1 : 0;
}

int main(int argc, char **argv) {
  data d;
  glam2_aln *alns;
  int r;
  bool output_error = false;

  prog_name = "glam2";  /* for error messages */
  getargs(&d.a, argc, argv);
  init(&d);

  fputs("GLAM2: Gapped Local Alignment of Motifs\nVersion "
#include "glam2_version.h"
	"\n\n", d.out);
  printargs(d.out, argc, argv);
  print_misc_info(d.out, &d);
  putc('\n', d.out);
  XMALLOC(alns, d.a.runs);

  for (r = 0; r < d.a.runs; ++r) {
    glam2_aln *aln = &alns[r];
    if (!d.a.quiet) {
      fprintf(stderr, "Run %d... ", r+1);
      fflush(stderr);
    }
    aln_init(aln, d.seqs.seqnum, d.a.max_width, d.alph.size);
    d.sm.underflow_flag = 1;  /* do we care about underflow in start_aln? */
    start_aln(aln, &d);
    optimise_aln(aln, &d);
    if (d.sm.underflow_flag < (d.a.algorithm == 2 ? DBL_EPSILON : DBL_MIN))
      fprintf(stderr, "%s: accuracy loss due to numeric underflow (%g)\nIf the alignment looks suspect, try rerunning with higher -u, or maybe lower -b\n", prog_name, d.sm.underflow_flag);
    if (d.a.profile)
      print_aln_info(d.out, aln, &d);
  }

  if (!d.a.quiet) putc('\n', stderr);

  SORT(alns, d.a.runs, aln_cmp);
  if (!d.a.profile)
    print_alns(d.out, alns, &d);

  xfclose(d.out);			// close text output file

  // Create the HTML output and MEME format output
  char *glam2html, *glam2psfm, *command;
  int glam2htmllen, glam2psfmlen, command_length, command_ret;
  // create the paths to the programs
  glam2html = get_meme_libexec_file("glam2html");
  glam2psfm = get_meme_libexec_file("glam2psfm");
  // allocate memory for the command
  glam2htmllen = (glam2html != NULL ? strlen(glam2html) : 0);
  glam2psfmlen = (glam2psfm != NULL ? strlen(glam2psfm) : 0);
  command_length = (glam2htmllen > glam2psfmlen ? glam2htmllen : glam2psfmlen) +
    strlen(d.txt_filename) + strlen(d.html_filename) + 50;
  command = mm_malloc(command_length);
  if (glam2html != NULL) {
    // run glam2html
    sprintf(command, "%s < %s > %s",  glam2html, d.txt_filename, d.html_filename);
    if ((command_ret = system(command)) != 0) {
      output_error = true;
      report_external_failure("glam2html", command_ret);
      fprintf(stderr, "Warning: failed to convert output to HTML!\n");
    }
  } else {
    output_error = true;
    fprintf(stderr, "Warning: could not find glam2html script! Failed to convert output to HTML!\n");
  }
  if (glam2psfm != NULL) {
    // run glam2psfm
    sprintf(command, "%s < %s > %s", glam2psfm, d.txt_filename, d.psfm_filename);
    if ((command_ret = system(command)) != 0) {
      output_error = true;
      report_external_failure("glam2psfm", command_ret);
      fprintf(stderr, "Warning: failed to convert output to MEME format motif!\n");
    }
  } else {
    output_error = true;
    fprintf(stderr, "Warning: could not find glam2psfm script! Failed to convert to MEME format motif!\n");
  }

  // Free everything.
  int i, j;
  free(command);
  free(glam2psfm);
  free(glam2html);
  alns_destroy(alns, d.a.runs, d.a.max_width);
  free(alns);

  free(d.txt_filename);
  free(d.html_filename);
  free(d.psfm_filename);
  free_mfasta(&(d.seqs));
  free(d.seq_order);
  alns_destroy(&(d.aln), 1, d.a.max_width);

  free2(d.sm.dp_mat, d.a.max_width+1);
  free2(d.sm.match_scores, d.a.max_width);
  free(d.sm.delete_scores);
  free(d.sm.insert_scores);

  free(d.scorer.e_prior.weights);
  for (i=0; i<d.scorer.e_prior.comp_num; i++) {
    free(d.scorer.e_prior.components[i].alpha);
    free(d.scorer.e_prior.components[i].alpha_lookup);
    free(d.scorer.e_prior.components[i].sum_lookup);
  }
  free(d.scorer.e_prior.components);
  free(d.scorer.e_prior.log_weights);
  free(d.scorer.e_prior.scratch);
  free(d.scorer.e_prior.counts);
  free(d.scorer.e_prior.offsets);

  free(d.scorer.d_prior.alpha_lookup.table);
  free(d.scorer.d_prior.beta_lookup.table);
  free(d.scorer.d_prior.sum_lookup.table);

  free(d.scorer.i_prior.alpha_lookup.table);
  free(d.scorer.i_prior.beta_lookup.table);
  free(d.scorer.i_prior.sum_lookup.table);

  free(d.scorer.bg.counts);
  free(d.scorer.bg.probs);
  free(d.scorer.bg.log_probs);

  free(d.alph.prob);

  free(d.col_sampler.offsets);
  free(d.col_sampler.fits);
  free(d.col_sampler.scores);
  free(d.col_sampler.del_probs);
  free(d.col_sampler.col.positions);
  free(d.col_sampler.col.matches);
  free(d.col_sampler.col.emission_counts);
  
  return (output_error ? 1 : 0);
}
