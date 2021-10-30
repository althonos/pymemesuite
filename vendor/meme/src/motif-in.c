#include <assert.h>
#include <string.h>

#include "motif-in.h"
#include "motif-in-dreme-xml.h"
#include "motif-in-meme-xml.h"
#include "motif-in-meme-html.h"
#include "motif-in-meme-json.h"
#include "motif-in-meme-text.h"
#include "motif-in-streme-xml.h"
#include "motif-spec.h"
#include "utils.h"

typedef struct mformat MFORMAT_T;
/*
 * Motif Reader
 */
struct mread {
  char *filename;             // The optional motif filename (may be "-" for stdin)
  int options;                // The options enabled
  ALPH_T *conv_alph;          // Convert the motifs to this alphabet
  double pseudo_total;        // The pseudocount to be added to the motifs
  ARRAY_T *pseudo_bg;         // The background to use to apply the pseudocount
  ARRAY_T *motif_bg;          // The background read from the motifs
  ARRAY_T *other_bg;          // The background derived from another source
  char *other_bg_src;         // The source the background came from (file or keyword)
  MFORMAT_T *formats;         // The motif readers
  int valid;                  // The count of valid readers left
  int total;                  // The total readers left
  bool success;          // Has a reader returned a motif
  FILE *fp;                   // Optional file to read from
  double trim_bits;           // Trim the motifs from the edges
  bool trim;             // Should the motif be trimmed?
  int count;                  // The current count of processed motifs
  bool symmetrical;	      // Should the pseudocount background be made symmetrical if possible?
};

/*
 * Motif Format Reader
 */
struct mformat {
  char *name;
  void *data;
  bool valid; // is the parser still valid
  void (*destroy)(void *data);
  void (*update)(void *data, const char *buffer, size_t size, short end);
  short (*match_score)(void *data); // larger number means better match
  short (*has_warning)(void *data);
  char* (*next_warning)(void *data);
  short (*has_error)(void *data);
  char* (*next_error)(void *data);
  short (*has_motif)(void *data);
  MOTIF_T* (*next_motif)(void *data);
  ALPH_T* (*get_alphabet)(void *data);
  int (*get_strands)(void *data);
  bool (*get_bg)(void *data, ARRAY_T **bg);
  void* (*motif_optional)(void *data, int option); // optional component at motif scope
  void* (*file_optional)(void *data, int option); // optional component at file scope
};

/*
 * Compare the match scores of motif format datastructures
 */
int mformat_cmp(const void *p1, const void *p2) {
  const MFORMAT_T *mf1, *mf2;
  mf1 = (MFORMAT_T*)p1;
  mf2 = (MFORMAT_T*)p2;
  return mf2->match_score(mf2->data) - mf1->match_score(mf1->data);
}

/*
 * Add a motif format datastructure to the list of formats
 */
void add_format(
  MREAD_T *mread, 
  char *name,
  void* (*init)(const char *filename, int options), // use the filename to get a type hint
  void (*destroy)(void *data),
  void (*update)(void *data, const char *buffer, size_t size, short end),
  short (*match_score)(void *data),
  short (*has_warning)(void *data),
  char* (*next_warning)(void *data),
  short (*has_error)(void *data),
  char* (*next_error)(void *data),
  short (*has_motif)(void *data),
  MOTIF_T* (*next_motif)(void *data),
  ALPH_T* (*get_alphabet)(void *data),
  int (*get_strands)(void *data),
  bool (*get_bg)(void *data, ARRAY_T **bg),
  void* (*motif_optional)(void *data, int option), // optional component at motif scope
  void* (*file_optional)(void *data, int option) // optional component at file scope
) {
  MFORMAT_T *format;
  mread->formats = mm_realloc(mread->formats, (mread->total + 1) * sizeof(MFORMAT_T));
  format = mread->formats+(mread->total);
  mread->total++;
  mread->valid++;

  format->name = name;
  format->data = init(mread->filename, mread->options);
  format->valid = true;
  format->destroy = destroy;
  format->update = update;
  format->match_score = match_score;
  format->has_warning = has_warning;
  format->next_warning = next_warning;
  format->has_error = has_error;
  format->next_error = next_error;
  format->has_motif = has_motif;
  format->next_motif = next_motif;
  format->get_alphabet = get_alphabet;
  format->get_strands = get_strands;
  format->get_bg = get_bg;
  format->motif_optional = motif_optional;
  format->file_optional = file_optional;
}
  
/*
 * Sets the background for the pseudo-counts
 * but requires the alphabet to do it.
 */
static void set_pseudo_bg(
  MREAD_T *mread,
  ALPH_T *alph
) {
  if (!alph) alph = mread->formats->get_alphabet(mread->formats->data);
  assert(alph != NULL);
  if (!mread->other_bg) {
    // no change to other_bg
    if (mread->other_bg_src == NULL) mread->other_bg_src = strdup("--nrdb--");
    if (strcmp(mread->other_bg_src, "--nrdb--") == 0) {
      // Note: get_nrdb_frequencies returns NULL if alph is not standard DNA or protein.
      mread->other_bg = get_nrdb_frequencies(alph, NULL);
    } else if (strcmp(mread->other_bg_src, "--uniform--") == 0) {
      mread->other_bg = get_uniform_frequencies(alph, NULL);
    } else if (strcmp(mread->other_bg_src, "--motif--") == 0 
      || strcmp(mread->other_bg_src, "motif-file") == 0
      || strcmp(mread->other_bg_src, "--query--") == 0) {
      // motif_bg is loaded elsewhere
    } else {
      int order = 0;
      char *bg_syms = NULL;
      mread->other_bg = load_markov_model_without_alph(mread->other_bg_src, &order, &bg_syms); 
      // check that the symbols in the background model match those in the alphabet
      STR_T *alph_letters_string = str_create(100);
      const char *alph_syms = alph_string(alph, alph_letters_string);
      if (strcmp(alph_syms, bg_syms) != 0) {
        die("Error: The symbols in background file (%s) don't match those in the alphabet (%s).\n", bg_syms, alph_syms);
      }
      str_destroy(alph_letters_string, false);
      free(bg_syms);
    }
  }
  if (mread->other_bg) {
    mread->pseudo_bg = mread->other_bg;
  } else {
    // other_bg wasn't set yet above.
    // Can't use NRDB, so use background in motif.
    if (mread->other_bg_src) free(mread->other_bg_src);
    mread->other_bg_src = strdup("--motif--");
    mread->pseudo_bg = mread->motif_bg;
  }
  // TLB 26-Jul-2017; Make symmetrical if scanning both strands
  if (alph_has_complement(alph) && mread->symmetrical) 
    average_freq_with_complement(alph, mread->pseudo_bg);
  // TLB 8-Feb-2017; Make sure the background doesn't have any zeros.
  normalize_frequencies(mread->conv_alph ? mread->conv_alph : alph, mread->pseudo_bg, 0.0000005);
} // set_pseudo_bg

/*
 * When the parser has been selected do some processing
 */
static void parser_selected(
  MREAD_T *mread
) {
  ALPH_T* alph;
  MFORMAT_T* format;
  format = mread->formats;
  // get the alphabet
  alph = format->get_alphabet(mread->formats->data);
  // get the background
  if (format->get_bg(format->data, &(mread->motif_bg))) {
    normalize_subarray(0, alph_size_core(alph), 0.0, mread->motif_bg);
  } else {
    mread->motif_bg = get_uniform_frequencies(alph, mread->motif_bg);
  }
  set_pseudo_bg(mread, alph);
}

/*
 * Update the parser and reorganise the formats list
 * if an error occurs or the parser produces a motif
 * and is selected.
 */
static void update_parser(
  MREAD_T *mread, 
  MFORMAT_T *format,  
  const char *buffer, 
  size_t size, 
  short end
) {
  // update the parser
  format->update(format->data, buffer, size, end);
  // reorganise the parser list
  if (format->has_error(format->data)) {
    // error occurred! parser not valid
    format->valid = false;
    if (format != mread->formats+(mread->valid - 1)) {
      // not the last valid in the list so we have to swap
      MFORMAT_T temp;
      temp = *format;
      *format = mread->formats[mread->valid -1];
      mread->formats[mread->valid -1] = temp;
    }
    // decrement valid count
    mread->valid--;
  } else if (!mread->success && format->has_motif(format->data)) {
    // a parser is chosen!
    if (format != mread->formats) {
      // not the first in the list so swap into first position
      MFORMAT_T temp;
      temp = *format;
      *format = mread->formats[0];
      mread->formats[0] = temp;
    }
    mread->valid = 1;
    mread->success = true;
    parser_selected(mread);
  }
} // update_parser

/*
 * loop through all the warnings in the format
 * and print them to stderr.
 */
static void print_parser_warnings(MFORMAT_T *format) {
  char *warning;
  size_t len;
  // skip any format that does not support warnings
  if (format->has_warning == NULL || format->next_warning == NULL) return;
  while (format->has_warning(format->data)) {
    warning = format->next_warning(format->data);
    len = strlen(warning);
    fputs(warning, stderr);
    if (len > 0 && warning[len - 1] != '\n') fputs("\n", stderr);
    free(warning);
  }
}

/*
 * loop through all the errors in the format
 * and print them to stderr.
 */
static void print_parser_errors(MFORMAT_T *format) {
  char *error;
  size_t len;

  while (format->has_error(format->data)) {
    error = format->next_error(format->data);
    len = strlen(error);
    fputs(error, stderr);
    if (len > 0 && error[len - 1] != '\n') fputs("\n", stderr);
    free(error);
  }
}

/*
 * Output any errors produced by the selected parser
 * or if all content has been parsed then output
 * errors produced by the most likely parser matches.
 */
static void output_errors(MREAD_T *mread, short end) {
  MFORMAT_T* format;
  int i;
  short best;
  if (mread->success) {
    // a single parser has been selected so output
    // any errors that it has produced
    print_parser_warnings(mread->formats);
    print_parser_errors(mread->formats);
  } else if (end) {
    // no parsers were selected so try to print the most relevant errors
    // sort parsers by file match
    qsort(mread->formats, mread->total, sizeof(MFORMAT_T), mformat_cmp);
    // get the best format match score
    best = mread->formats->match_score(mread->formats->data);
    // check that the best format match is plausible
    if (best < 1) {
      fputs("There were no convincing matches to any MEME Suite motif "
          "format.\n", stderr);
      return;
    }
    // print out the warnings and errors for the best matching formats
    for (i = 0; i < mread->total; i++) {
      format = mread->formats+i;
      if (format->match_score(format->data) < best) break;
      if (format->has_warning != NULL && format->has_warning(format->data)) {
        fprintf(stderr, "Warnings from %s parser:\n", format->name);
        print_parser_warnings(format);
      }
      if (format->has_error(format->data)) {
        fprintf(stderr, "Errors from %s parser:\n", format->name);
        print_parser_errors(format);
      }
    }
  }
}

/*
 * Do some post processing of the motif to add pseudocounts, calculate
 * the complexity and trim it
 */
#define PSEUDO 0.01
#define DEFAULT_SITES 100
static MOTIF_T* post_process_motif(MREAD_T *mread, MOTIF_T *motif) {
  ARRAY_T *bg;
  MATRIX_T *temp_freqs;
  int site_count;
  if (mread->options & SKIP_POST_PROCESSING) return motif;
  if (motif == NULL) return NULL;
  assert(motif->alph != NULL);
  if (motif->freqs) normalize_motif(motif, 0.00001);
  site_count = (motif->num_sites > 0 ? motif->num_sites : DEFAULT_SITES);
  if (motif->freqs != NULL && motif->scores != NULL) {
    // validate? May not be possible for protein motifs as MEME tweaks the PSSM
  } else if (motif->scores != NULL) {
    // calculate the freqs
    motif->freqs = convert_scores_into_freqs(motif->alph, motif->scores, 
      mread->motif_bg, site_count, PSEUDO);
  } else if (motif->freqs != NULL) {
    // calculate the scores
    motif->scores = convert_freqs_into_scores(motif->alph, motif->freqs,
      //mread->motif_bg, site_count, PSEUDO);  // this is always from the motif, ignoring bfile
      mread->pseudo_bg, site_count, PSEUDO); 	// TLB changed 1-Feb-2017 so -bfile affects scores for MAST
  } else {
    die("Motif with no PSPM or PSSM should not get here!\n");
  }
  // if requested we convert the alphabet
  if (mread->conv_alph) {
    temp_freqs = motif->freqs;
    motif->freqs = convert_matrix_alphabet(temp_freqs, 0.0, motif->alph, mread->conv_alph);
    free_matrix(temp_freqs);
    alph_release(motif->alph);
    motif->alph = alph_hold(mread->conv_alph);
    free_matrix(motif->scores);
    // can't use motif_bg here because it is for the old alphabet
    motif->scores = convert_freqs_into_scores(motif->alph, motif->freqs,
      mread->pseudo_bg, site_count, PSEUDO);
  }
  apply_pseudocount_to_motif(motif, mread->pseudo_bg, mread->pseudo_total);
  motif->complexity = compute_motif_complexity(motif);
  motif->idx = ++(mread->count);
  if (mread->options & CALC_AMBIGS) calc_motif_ambigs(motif);
  if (mread->trim) trim_motif_by_bit_threshold(motif, mread->trim_bits);

  // Compute a single-letter consensus for the motif.
  STR_T *cons_buf = str_create(50);
  str_clear(cons_buf);
  motif2consensus(motif, cons_buf, true);
  motif->consensus = str_internal(cons_buf);
  free(cons_buf);

  return motif;
} // post_process_motif

/*
 * create the reader passing an optional
 * filename to help determine the most
 * likely file type.
 */
MREAD_T* mread_create(
  const char *filename, 
  int options, 
  bool symmetrical	// make background symmetrical if alphabet is complementable
) {
  MREAD_T * mread = mm_malloc(sizeof(MREAD_T));
  memset(mread, 0, sizeof(MREAD_T));
  if (filename) mread->filename = strdup(filename);
  if (options & OPEN_MFILE) {
    if (strcmp(filename, "-") == 0) {
      mread->fp = stdin;
    } else {
      mread->fp = fopen(filename, "r");
      if (!mread->fp) die("Failed to open motif file: %s\n", filename);
    }
  }
  mread->options = options;
  mread->valid = 0;
  mread->total = 0;
  mread->count = 0;
  mread->conv_alph = NULL;
  mread->pseudo_total = 0;
  mread->other_bg_src = strdup("--motif--");
  mread->symmetrical = symmetrical;
  // create readers of each format
  add_format(mread, "MEME XML", mxml_create, mxml_destroy, mxml_update,
      mxml_has_format_match, mxml_has_warning, mxml_next_warning,
      mxml_has_error, mxml_next_error, 
      mxml_has_motif, mxml_next_motif, mxml_get_alphabet, mxml_get_strands, 
      mxml_get_bg, mxml_motif_optional, mxml_file_optional);
  add_format(mread, "MEME HTML", mhtml_create, mhtml_destroy, mhtml_update,
      mhtml_has_format_match, NULL, NULL, mhtml_has_error, mhtml_next_error, 
      mhtml_has_motif, mhtml_next_motif, mhtml_get_alphabet, mhtml_get_strands,
      mhtml_get_bg, mhtml_motif_optional, mhtml_file_optional);
  add_format(mread, "MEME HTML+JSON", mhtml2_create, mhtml2_destroy, mhtml2_update,
      mhtml2_has_format_match, mhtml2_has_warning, mhtml2_next_warning,
      mhtml2_has_error, mhtml2_next_error, 
      mhtml2_has_motif, mhtml2_next_motif, mhtml2_get_alphabet, mhtml2_get_strands,
      mhtml2_get_bg, mhtml2_motif_optional, mhtml2_file_optional);
  add_format(mread, "MEME text", mtext_create, mtext_destroy, mtext_update,
      mtext_has_format_match, mtext_has_warning, mtext_next_warning, mtext_has_error, mtext_next_error,
      mtext_has_motif, mtext_next_motif, mtext_get_alphabet, mtext_get_strands,
      mtext_get_bg, mtext_motif_optional, mtext_file_optional);
  add_format(mread, "DREME XML", dxml_create, dxml_destroy, dxml_update,
      dxml_has_format_match, dxml_has_warning, dxml_next_warning, dxml_has_error, dxml_next_error, 
      dxml_has_motif, dxml_next_motif, dxml_get_alphabet, dxml_get_strands,
      dxml_get_bg, dxml_motif_optional, dxml_file_optional);
  add_format(mread, "STREME XML", sxml_create, sxml_destroy, sxml_update,
      sxml_has_format_match, sxml_has_warning, sxml_next_warning, sxml_has_error, sxml_next_error, 
      sxml_has_motif, sxml_next_motif, sxml_get_alphabet, sxml_get_strands,
      sxml_get_bg, sxml_motif_optional, sxml_file_optional);
  return mread;
}

/*
 * Convert the motifs to the given alphabet (or die trying).
 *
 * This cancels any previously set background option.
 */
void mread_set_conversion(
  MREAD_T *mread, 
  ALPH_T *alph, 
  const ARRAY_T *bg
) {

  // Sanity check.
  if (bg != NULL) {
    if (get_array_length(bg) < alph_size_core(alph)) {
      die("Background from '%s' too short for alphabet '%s'.\n"
       "       Has %d entries but needs at least %d.", 
        (mread->other_bg_src ? mread->other_bg_src : "null"),
        alph_name(alph), get_array_length(bg), alph_size_core(alph));
    }
  }

  // clean up old bg
  if (mread->other_bg != NULL) free_array(mread->other_bg);
  mread->other_bg = NULL;
  if (mread->other_bg_src) free(mread->other_bg_src);
  mread->other_bg_src = NULL;

  // copy bg
  if (bg != NULL) {
    if (mread->other_bg_src != NULL) free(mread->other_bg_src);
    mread->other_bg = allocate_array(get_array_length(bg));
    copy_array(bg, mread->other_bg);
  } else {
    // no background so generate one based on the alphabet
    mread->other_bg = get_uniform_frequencies(alph, NULL);
  }

  // copy alphabet
  mread->conv_alph = alph_hold(alph);

  // set the pseudocount background
  if (mread->success) set_pseudo_bg(mread, mread->conv_alph);
} // mread_set_conversion

/*
 * Set the background to be used for pseudocounts.
 *
 * This cancels a conversion alphabet.
 */
void mread_set_background(MREAD_T *mread, const ARRAY_T *bg, ALPH_T *alph) {
  // clean up old bg
  if (mread->other_bg != NULL) free_array(mread->other_bg);
  mread->other_bg = NULL;
  if (mread->other_bg_src != NULL) free(mread->other_bg_src);
  mread->other_bg_src = NULL;
  // clean up conversion alphabet
  if (mread->conv_alph != NULL) alph_release(mread->conv_alph);
  mread->conv_alph = NULL;
  // copy the passed background
  if (bg != NULL) {
    mread->other_bg = allocate_array(get_array_length(bg));
    copy_array(bg, mread->other_bg);
  }
  // check if we've chosen a parser already
  if (mread->success) set_pseudo_bg(mread, alph);
}

/*
 * Set the background source to be used for pseudocounts.
 * The background is resolved by the following rules:
 *   If source is:
 *     null or "--nrdb--"			use nrdb frequencies
 *     "--uniform--"				use uniform frequencies
 *     "motif-file", "--motif--" or "--query"	use frequencies in motif file
 *     file name                		read frequencies from bg file
 *
 *   This cancels a conversion alphabet.
 */
void mread_set_bg_source(
  MREAD_T *mread, 
  const char *source,
  ALPH_T *alph
) {
  // clean up old bg
  if (mread->other_bg != NULL) free_array(mread->other_bg);
  mread->other_bg = NULL;
  if (mread->other_bg_src != NULL) free(mread->other_bg_src);
  mread->other_bg_src = NULL;

  // clean up conversion alphabet
  if (mread->conv_alph != NULL) alph_release(mread->conv_alph);
  mread->conv_alph = NULL;

  // copy the passed source
  if (source != NULL) {
    if (mread->other_bg_src != NULL) free(mread->other_bg_src);
    mread->other_bg_src = strdup(source);
  }

  // Set the background for pseudocounts.
  if (mread->success) set_pseudo_bg(mread, alph);
} // mread_set_bg_source

/*
 * Set the pseudo-count to be applied to the motifs.
 * Uses the motif background unless a background has been set.
 */
void mread_set_pseudocount(MREAD_T *mread, double pseudocount) {
  mread->pseudo_total = pseudocount;
}

/*
 * Set the trimming threshold to be applied to the motifs
 */
void mread_set_trim(MREAD_T *mread, double trim_bits) {
  mread->trim_bits = trim_bits;
  mread->trim = true;
}

void mread_destroy(MREAD_T *mread) {
  int i;
  for (i = 0; i < mread->total; ++i) {
    MFORMAT_T* format = mread->formats+i;
    format->destroy(format->data);
  }
  memset(mread->formats, 0, sizeof(MFORMAT_T) * mread->total);
  free(mread->formats);
  if (mread->options & OPEN_MFILE && mread->fp != stdin) {
    fclose(mread->fp);
  }
  if (mread->conv_alph) alph_release(mread->conv_alph);
  if (mread->other_bg) free_array(mread->other_bg);
  if (mread->other_bg_src) free(mread->other_bg_src);
  if (mread->motif_bg) free_array(mread->motif_bg);
  if (mread->filename) free(mread->filename);
  memset(mread, 0, sizeof(MREAD_T));
  free(mread);
}

static bool motif_avaliable(MREAD_T *mread) {
  // check that we know the format
  if (mread->success) {
    MFORMAT_T *format = mread->formats;
    return format->valid && format->has_motif(format->data);
  }
  return false;
}

#define MREAD_BUFFER_SIZE 100
static inline void attempt_read_motif(MREAD_T *mread) {
  int size, terminate;
  char chunk[MREAD_BUFFER_SIZE + 1];
  if (mread->fp) {
    terminate = feof(mread->fp);
    while (!motif_avaliable(mread) && !terminate) {
      size = fread(chunk, sizeof(char), MREAD_BUFFER_SIZE, mread->fp);
      chunk[size] = '\0';
      terminate = feof(mread->fp);
      mread_update(mread, chunk, size, terminate);
    }
  }
}

/*
 * Is there another motif to return
 */
bool mread_has_motif(MREAD_T *mread) {
  attempt_read_motif(mread);
  return motif_avaliable(mread);
}

/*
 * Get the next motif. If there is no motif ready then
 * returns null.
 */
MOTIF_T* mread_next_motif(MREAD_T *mread) {
  if (mread_has_motif(mread)) {
    return post_process_motif(mread, 
        mread->formats->next_motif(mread->formats->data));
  }
  return NULL;
}

/*
 * Update the parsers.
 *
 * Send all parsers the information until they error or one produces a motif.
 * Once a parser has produced a motif then assume it is the correct parser and
 * only send data to it. If no parsers produce motifs by then end of the file
 * then sort the parsers by the most likely file match, eliminate any which
 * have an unlikely match, and print the errors for the most likely file type 
 * matches (in case of a tie multiple are allowed).
 */
void mread_update(MREAD_T *mread, const char *buffer, size_t size, short end) {
  int i;
  for (i = (mread->valid - 1); i >= 0; --i) {
    update_parser(mread, mread->formats+i, buffer, size, end);
    if (mread->success) break; // avoid sending the same chunk twice
  }
  output_errors(mread, end);
}

/*
 * Loads all the currently buffered motifs into a list.
 * If the file is set then this will read all the motifs in the
 * file into the list. If a list is not passed then
 * it will create a new one. 
 * returns the list.
 */
ARRAYLST_T* mread_load(MREAD_T *mread, ARRAYLST_T *motifs) {
  MOTIF_T *motif;
  if (motifs == NULL) motifs = arraylst_create();
  while ((motif = mread_next_motif(mread)) != NULL) {
    arraylst_add(motif, motifs);
  }
  return motifs;
}

/*
 * Get the background source string as a free-able copy.
*/
char* mread_get_other_bg_src(MREAD_T *mread) {
  if (mread->other_bg_src)
    return strdup(mread->other_bg_src);
  else
    return NULL;
}

/*
 * Get the alphabet. If the alphabet is unknown then returns NULL.
 * WARNING: does not automatically prevent the alphabet from being deallocated.
 * If users wish to keep the alphabet they must call alph_hold on it themselves.
 */
ALPH_T* mread_get_alphabet(MREAD_T *mread) {
  attempt_read_motif(mread);
  if (mread->valid == 1) {
    MFORMAT_T *format = mread->formats;
    return format->get_alphabet(format->data);
  }
  return NULL;
}

/*
 * Get the strands. 
 * DNA motifs have 1 or 2 strands,
 * otherwise return 0 (for both protein and unknown)
 */
int mread_get_strands(MREAD_T *mread) {
  attempt_read_motif(mread);
  if (mread->valid == 1) {
    MFORMAT_T *format = mread->formats;
    return format->get_strands(format->data);
  }
  return 0;
}

/*
 * Get the background used to assign pseudo-counts. 
 * If the background has not yet been determined it will return null.
 */
ARRAY_T *mread_get_background(MREAD_T *mread) {
  ARRAY_T *bg;
  ALPH_T *alph;
  alph = mread_get_alphabet(mread);
  if (mread->pseudo_bg) {
    bg = allocate_array(get_array_length(mread->pseudo_bg));
    copy_array(mread->pseudo_bg, bg);
    if (mread->options & CALC_AMBIGS) calc_ambigs(alph, false, bg);
    return bg;
  }
  return NULL;
}

/*
 * Get the motif occurrences. If the motif occurrences are unavailable
 * then returns null. Only returns once.
 */
ARRAYLST_T *mread_get_motif_occurrences(MREAD_T *mread) {
  // check that we know the format
  if (mread->valid == 1) {
    MFORMAT_T *format = mread->formats;
    return (ARRAYLST_T*)format->file_optional(format->data, SCANNED_SITES);
  }
  return NULL;
}

/*
 * Destroy the motif occurrences.
 */
void destroy_motif_occurrences(ARRAYLST_T *motif_occurrences) {
  if (motif_occurrences != NULL) {
    arraylst_destroy(sseq_destroy, motif_occurrences);
  }
}
