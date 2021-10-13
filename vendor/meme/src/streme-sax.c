#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "alphabet.h"
#include "streme-sax.h"
#include "linked-list.h"
#include "red-black-tree.h"
#include "sax-parser-utils.h"
#include "utils.h"

/*****************************************************************************
 * Parser States
 ****************************************************************************/
enum parser_state {
  PS_ERROR,
  PS_START,
  PS_IN_STREME,
  PS_IN_MODEL,
  PS_IN_COMMAND_LINE,
  PS_IN_TRAIN_POSITIVES,
  PS_IN_TRAIN_NEGATIVES,
  PS_IN_TEST_POSITIVES,
  PS_IN_TEST_NEGATIVES,
  PS_IN_ALPHABET,
  PS_IN_ALPHABET_LETTER,
  PS_IN_STRANDS,
  PS_IN_BACKGROUND,
  PS_IN_SEQUENCE_DB,
  PS_IN_BACKGROUND_FREQUENCIES,
  PS_IN_BF_ALPHABET_ARRAY,
  PS_IN_BF_AA_VALUE,
  PS_IN_STOP,
  PS_IN_OBJFUN,
  PS_IN_TEST,
  PS_IN_MINW,
  PS_IN_MAXW,
  PS_IN_KMER,
  PS_IN_HOFRACT,
  PS_IN_NEVAL,
  PS_IN_NREF,
  PS_IN_NITER,
  PS_IN_PATIENCE,
  PS_IN_SEED,
  PS_IN_USEER,
  PS_IN_MINSCORE,
  PS_IN_IGNORE_DEPTH,
  PS_IN_NSUBSETS,
  PS_IN_MIN_PAL_RATIO,
  PS_IN_MAX_PAL_ED,
  PS_IN_CAND,
  PS_IN_EXPERIMENTAL,
  PS_IN_TOTALLENGTH,
  PS_IN_ALIGN,
  PS_IN_HOST,
  PS_IN_DESCRIPTION,
  PS_IN_MOTIFS,
  PS_IN_MOTIF,
  PS_IN_POS,
  PS_IN_REASON_FOR_STOPPING,
  PS_IN_RUN_TIME,
  PS_END
};
typedef enum parser_state PS_EN;

/*****************************************************************************
 * Parser State Names (for error messages)
 ****************************************************************************/
static char const * const state_names[] = {
  "ERROR",
  "START",
  "IN_STREME",
  "IN_MODEL",
  "IN_COMMAND_LINE",
  "IN_TRAIN_POSITIVES",
  "IN_TRAIN_NEGATIVES",
  "IN_TEST_POSITIVES",
  "IN_TEST_NEGATIVES",
  "IN_ALPHABET",
  "IN_ALPHABET_LETTER",
  "IN_STRANDS",
  "IN_SEQUENCE_DB",
  "IN_BACKGROUND",
  "IN_BACKGROUND_FREQUENCIES",
  "IN_STOP",
  "IN_OBJFUN",
  "IN_TEST",
  "IN_MINW",
  "IN_MAXW",
  "IN_KMER",
  "IN_HOFRACT",
  "IN_NEVAL",
  "IN_NREF",
  "IN_NITER",
  "IN_PATIENCE",
  "IN_SEED",
  "IN_USEER",
  "IN_MINSCORE",
  "IN_IGNORE_DEPTH",
  "IN_NSUBSETS",
  "IN_MIN_PAL_RATIO",
  "IN_MAX_PAL_ED",
  "IN_CAND",
  "IN_EXPERIMENTAL",
  "IN_HOST",
  "IN_MOTIFS",
  "IN_MOTIF",
  "IN_POS",
  "IN_REASON_FOR_STOPPING",
  "IN_RUN_TIME",
  "END"
};

/*****************************************************************************
 * Expected State Repeats
 ****************************************************************************/
enum expected_state {
  ES_ZERO_OR_ONE,
  ES_ONCE,
  ES_ONE_OR_MORE,
  ES_ANY
};
typedef enum expected_state ES_EN;

/*****************************************************************************
 * Used to keep track of what elements are expected and enforces limits on the
 * number of times a valid element can be repeated. Also avoids the 
 * possibility that elements could be left out.
 ****************************************************************************/
typedef struct es ES_T;
struct es {
  PS_EN state;    // the element that is expected
  ES_EN expected; // the number of times the element is expected (see ES codes)
  int found;    // the number of times the element has been found already
};

/*****************************************************************************
 * Stores the parser state.
 ****************************************************************************/
typedef struct ps PS_T;
struct ps {
  int state;                // the state code (see PS codes)
  int udepth;               // the depth of unknown tags (disabled)
  ATTRBUF_T attrbuf;
  char *letter_id_buf;      // a buffer to store the current letter id
  int letter_id_buf_len;    // the length of the letter id buffer
  CHARBUF_T characters;     // a buffer to store characters between elements
  LINKLST_T *expected_stack;// a queue of expected states
  STREME_IO_XML_CALLBACKS_T *callbacks; // callbacks to the user
  struct prog_version ver;
  bool seen_alphabet;       // distinguish between old and new
  bool seen_ambig;          // ambiguous symbols must only appear after core
  RBTREE_T *alph_ids;       // list of alphabet ids and their positions
  double *freqs;            // buffer of frequencies
  char *motif_id;           // motif id
  int motif_width;            // motif length
  int last_pos;             // last seen motif pos
  void *user_data;          // the user data
};



/*****************************************************************************
 * Report an error to the user.
 ****************************************************************************/
static void error(PS_T *ps, char *fmt, ...) {
  va_list  argp;
  if (ps->callbacks->error) {
    va_start(argp, fmt);
    ps->callbacks->error(ps->user_data, fmt, argp);
    va_end(argp);
  }
  ps->state = PS_ERROR;
}

/*****************************************************************************
 * Report an attribute parse error to the user.
 ****************************************************************************/
void streme_attr_parse_error(void *state, int errcode, const char *tag, const char *attr, const char *value) {
  PS_T *ps = (PS_T*)state;
  if (errcode == PARSE_ATTR_DUPLICATE) {
    error(ps, "Duplicate attribute %s::%s.\n", tag, attr);
  } else if (errcode == PARSE_ATTR_BAD_VALUE) {
    error(ps, "Bad value \"%s\" for attribute %s::%s.\n", value, tag, attr);
  } else if (errcode == PARSE_ATTR_MISSING) {
    error(ps, "Missing required attribute %s::%s.\n", tag, attr);
  }
}

/*****************************************************************************
 * Push an expected state on the stack
 ****************************************************************************/
void streme_push_es(PS_T *ps, PS_EN expected_state, ES_EN expected_occurances) {
  ES_T *es;
  if (expected_state < PS_ERROR ||  expected_state > PS_END) {
    die("Bad state code!\n");
  }
  es = mm_malloc(sizeof(ES_T));
  es->state = expected_state;
  es->expected = expected_occurances;
  es->found = 0;
  linklst_push(es, ps->expected_stack);
}

/*****************************************************************************
 * At the start of a new element check that all previous expected elements
 * have been found and that the current element hasn't been repeated too 
 * many times.
 ****************************************************************************/
int streme_update_es(PS_T *ps, PS_EN next_state) {
  ES_T *es, old_es;
  if (next_state < PS_ERROR ||  next_state > PS_END) {
    die("Bad state code!\n");
  }
  while ((es = (ES_T*)linklst_pop(ps->expected_stack)) != NULL && es->state != next_state) {
    old_es = *es;
    free(es);
    es = NULL;
    switch (old_es.expected) {
      case ES_ONCE:
      case ES_ONE_OR_MORE:
        if (old_es.found < 1){
          error(ps, "Expected state %s not found!\n", state_names[old_es.state]);
          return false;
        }
      default:
        break;
    }
  }
  if (es != NULL) {
    es->found += 1;
    linklst_push(es, ps->expected_stack);
    switch (es->expected) {
      case ES_ONCE:
      case ES_ZERO_OR_ONE:
        if (es->found > 1) {
          error(ps, "Expected state %s only once!\n", state_names[es->state]);
          return false;
        }
      default:
        break;
    }
  } else {
    error(ps, "The state %s was not expected!\n", state_names[next_state]);
    return false;
  }
  return true;
}

/*****************************************************************************
 * Reads frequency attributes into the pre-allocated freqs array.
 ****************************************************************************/
static void parse_freq_attrs(PS_T *ps, const char* tag, const xmlChar **attrs) {
  int i, ncore, seen, *idx;
  char *end_ptr;
  double value, sum;
  RBNODE_T *node;
  bool seen_bad;
  ncore = rbtree_size(ps->alph_ids);
  // initilize the freqs array
  if (ps->freqs == NULL) ps->freqs = mm_malloc(sizeof(double) * ncore);
  // reset freqs array;
  for (i = 0; i < ncore; i++) ps->freqs[i] = -1;
  seen = 0;
  seen_bad = false;
  sum = 0.0;
  // iterate over attributes
  for (i = 0; attrs[i] != NULL; i += 2) {
    idx = (int*)rbtree_get(ps->alph_ids, attrs[i]);
    if (idx != NULL) {
      assert(*idx < ncore);
      if (ps->freqs[*idx] != -1) {
        streme_attr_parse_error(ps, PARSE_ATTR_DUPLICATE, tag, (const char*)attrs[i], NULL);
        continue;
      }
      seen++;
      errno = 0; // reset because we're about to check it
      value = strtod((const char*)attrs[i+1], &end_ptr);
      // allow out of range values, mainly because freqs can be very close to zero
      if (end_ptr == (const char*)attrs[i+1] || (errno && errno != ERANGE) || value < 0 || value > 1) {
        streme_attr_parse_error(ps, PARSE_ATTR_BAD_VALUE, tag, (const char*)attrs[i], (const char*)attrs[i+1]);
        ps->freqs[*idx] = 0; // mark frequence as seen, even though it's bad
        seen_bad = true;
        continue;
      }
      ps->freqs[*idx] = value;
      sum += value;
    }
  }
  // check we got everthing
  if (seen < ncore) {
    // identify what we're missing
    for (node = rbtree_first(ps->alph_ids); node != NULL; node = rbtree_next(node)) {
      idx = (int*)rbtree_value(node);
      if (ps->freqs[*idx] == -1) {
        streme_attr_parse_error(ps, PARSE_ATTR_MISSING, tag, (char*)rbtree_key(node), NULL);
      }
    }
  } else if (!seen_bad) {
    // check the frequencies sum to 1
    double delta = sum - 1;
    delta = (delta < 0 ? -delta : delta);
    if (delta > (0.001 * ncore)) {
      // STREME writes background probabilities to 3 decimal places so assuming 
      // the error on each is at maximum 0.001 then the total error for the 
      // sum must be less than or equal to 0.004
      error(ps, "Probabilities of %s do not sum to 1, got %g .\n", tag, sum);
    }
  }
}

/*****************************************************************************
 * STREME 
 * contains: model, motifs, reason_for_stopping, run_time
 *
 *  release         the release date.
 *  version         the program version.
 ****************************************************************************/
static void start_ele_streme(PS_T *ps, const xmlChar **attrs) {
  char *release;

  char* names[2] = {"release", "version"};
  int (*parsers[2])(char*, void*) = {ld_str, ld_version};
  void *data[2] = {&release, &(ps->ver)};
  bool required[2] = {true, true};
  bool done[2];

  parse_attributes(streme_attr_parse_error, ps, "streme", attrs, 2, names, parsers, data, required, done);
  
  if (ps->callbacks->start_streme && ps->state != PS_ERROR) {
    ps->callbacks->start_streme(ps->user_data, ps->ver.major, ps->ver.minor, ps->ver.patch, release);
  }
  streme_push_es(ps, PS_IN_RUN_TIME, ES_ONCE);
  streme_push_es(ps, PS_IN_REASON_FOR_STOPPING, ES_ONCE);
  streme_push_es(ps, PS_IN_MOTIFS, ES_ONCE);
  streme_push_es(ps, PS_IN_MODEL, ES_ONCE);
}

/*****************************************************************************
 * /STREME
 ****************************************************************************/
static void st_end_ele_streme(PS_T *ps) {
  if (ps->callbacks->end_streme) {
    ps->callbacks->end_streme(ps->user_data);
  }
}

/*****************************************************************************
 * STREME > model
 * contains: command_line, train_positives, train_negatives, test_positives,
 *          test_negatives, [sequence_db,] alphabet, strands, [background,] 
 *          sequence_db, background_frequencies, 
 *          stop, objfun, test, minw, maxw, kmer, hofract, neval, nref,
 *          niter, patience, seed, useer, minscore, ignore_depth, 
 *	    nsubsets, min_pal_ratio, max_pal_ed, cand, experimental, 
 *          totallength, [align], host, description?
 ****************************************************************************/
static void start_ele_model(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_model && ps->state != PS_ERROR) {
    ps->callbacks->start_model(ps->user_data);
  }
  streme_push_es(ps, PS_IN_DESCRIPTION, ES_ZERO_OR_ONE);
  streme_push_es(ps, PS_IN_HOST, ES_ONCE);
  streme_push_es(ps, PS_IN_ALIGN, ES_ZERO_OR_ONE);
  streme_push_es(ps, PS_IN_TOTALLENGTH, ES_ONCE);
  streme_push_es(ps, PS_IN_EXPERIMENTAL, ES_ONCE);
  streme_push_es(ps, PS_IN_CAND, ES_ONCE);
  streme_push_es(ps, PS_IN_MAX_PAL_ED, ES_ONCE);
  streme_push_es(ps, PS_IN_MIN_PAL_RATIO, ES_ONCE);
  streme_push_es(ps, PS_IN_NSUBSETS, ES_ONCE);
  streme_push_es(ps, PS_IN_IGNORE_DEPTH, ES_ONCE);
  streme_push_es(ps, PS_IN_MINSCORE, ES_ONCE);
  streme_push_es(ps, PS_IN_USEER, ES_ONCE);
  streme_push_es(ps, PS_IN_SEED, ES_ONCE);
  streme_push_es(ps, PS_IN_PATIENCE, ES_ONCE);
  streme_push_es(ps, PS_IN_NITER, ES_ONCE);
  streme_push_es(ps, PS_IN_NREF, ES_ONCE);
  streme_push_es(ps, PS_IN_NEVAL, ES_ONCE);
  streme_push_es(ps, PS_IN_HOFRACT, ES_ONCE);
  streme_push_es(ps, PS_IN_KMER, ES_ONCE);
  streme_push_es(ps, PS_IN_MAXW, ES_ONCE);
  streme_push_es(ps, PS_IN_MINW, ES_ONCE);
  streme_push_es(ps, PS_IN_TEST, ES_ONCE);
  streme_push_es(ps, PS_IN_OBJFUN, ES_ONCE);
  streme_push_es(ps, PS_IN_STOP, ES_ONCE);
  streme_push_es(ps, PS_IN_BACKGROUND_FREQUENCIES, ES_ZERO_OR_ONE);	// added after 5.3.0
  streme_push_es(ps, PS_IN_SEQUENCE_DB, ES_ZERO_OR_ONE);		// added after 5.3.0
  streme_push_es(ps, PS_IN_BACKGROUND, ES_ZERO_OR_ONE);			// removed after 5.3.0
  streme_push_es(ps, PS_IN_STRANDS, ES_ONCE);
  streme_push_es(ps, PS_IN_ALPHABET, ES_ONCE);
  streme_push_es(ps, PS_IN_SEQUENCE_DB, ES_ZERO_OR_ONE);		// removed after 5.3.0
  streme_push_es(ps, PS_IN_TEST_NEGATIVES, ES_ZERO_OR_ONE);
  streme_push_es(ps, PS_IN_TEST_POSITIVES, ES_ZERO_OR_ONE);
  streme_push_es(ps, PS_IN_TRAIN_NEGATIVES, ES_ZERO_OR_ONE);
  streme_push_es(ps, PS_IN_TRAIN_POSITIVES, ES_ONCE);
  streme_push_es(ps, PS_IN_COMMAND_LINE, ES_ONCE);
}

/*****************************************************************************
 * STREME > /model
 ****************************************************************************/
static void st_end_ele_model(PS_T *ps) {
  if (ps->callbacks->end_model && ps->state != PS_ERROR) {
    ps->callbacks->end_model(ps->user_data);
  }
}

/*****************************************************************************
 * STREME > model > /command_line
 * the command-line used to run the program.
 ****************************************************************************/
static void st_end_ele_command_line(PS_T *ps) {
  if (ps->callbacks->handle_command_line && ps->state != PS_ERROR) {
    ps->callbacks->handle_command_line(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * STREME > model > train_positives
 *
 *  count           the number of sequences in the positive training set
 *  positions 	    the number of characters in the positive training set
 *  file            the file containing the positive training set
 ****************************************************************************/
static void start_ele_train_positives(PS_T *ps, const xmlChar **attrs) {
  long count, positions;
  char *file;

  // Attributes must be listed in alphabetical order.
  char* names[3] = {"count", "file", "positions"};
  int (*parsers[3])(char*, void*) = {ld_long, ld_str, ld_long};
  void *data[3] = {&count, &file, &positions};
  bool required[3] = {true, true, true};
  bool done[3];

  parse_attributes(streme_attr_parse_error, ps, "train_positives", attrs, 3, names, parsers, data, required, done);

  if (ps->callbacks->handle_train_positives && ps->state != PS_ERROR) {
    ps->callbacks->handle_train_positives(ps->user_data, count, positions, file);
  }
}

/*****************************************************************************
 * STREME > model > train_negatives
 *
 *  count           the number of sequences in the negative dataset
 *  positions 	    the number of characters in the negative training set
 *  file            the file containing the negative training set
 *  from            the source of the negative dataset (eg shuffled positives)
 ****************************************************************************/
static void start_ele_train_negatives(PS_T *ps, const xmlChar **attrs) {
  long count, positions;
  char *file = NULL;
  int from;

  // Attributes must be listed alphabetically.
  char* from_options[3] = {"file", "none", "shuffled"};
  int from_values[3] = {STREME_NEG_FILE, STREME_NEG_SHUFFLED, STREME_NEG_NONE};
  MULTI_T from_multi = {.count = 3, .options = from_options, .outputs = from_values, .target = &(from)};

  // Attributes must be listed alphabetically.
  char* names[4] = {"count", "file", "from", "positions"};
  int (*parsers[4])(char*, void*) = {ld_long, ld_str, ld_multi, ld_long};
  void *data[4] = {&count, &file, &from_multi, &positions};
  bool required[4] = {true, false, true, true};
  bool done[4];

  parse_attributes(streme_attr_parse_error, ps, "train_negatives", attrs, 4, names, parsers, data, required, done);

  if (ps->state != PS_ERROR && from == STREME_NEG_FILE) {
    if (file == NULL) {
      streme_attr_parse_error(ps, PARSE_ATTR_MISSING, "train_negatives", "file", NULL);
    } 
  }

  if (ps->callbacks->handle_train_negatives && ps->state != PS_ERROR) {
    ps->callbacks->handle_train_negatives(ps->user_data, count, positions, (STREME_NEG_EN)from, file);
  }
}

/*****************************************************************************
 * STREME > model > test_positives
 *
 *  count           the number of sequences in the positive training set
 *  positions 	    the number of characters in the positive training set
 ****************************************************************************/
static void start_ele_test_positives(PS_T *ps, const xmlChar **attrs) {
  long count, positions;

  char* names[2] = {"count", "positions"};
  int (*parsers[2])(char*, void*) = {ld_long, ld_long};
  void *data[2] = {&count, &positions};
  bool required[2] = {true, true};
  bool done[2];

  parse_attributes(streme_attr_parse_error, ps, "test_positives", attrs, 2, names, parsers, data, required, done);

  if (ps->callbacks->handle_test_positives && ps->state != PS_ERROR) {
    ps->callbacks->handle_test_positives(ps->user_data, count, positions);
  }
}

/*****************************************************************************
 * STREME > model > test_negatives
 *
 *  count           the number of sequences in the negative training set
 *  positions 	    the number of characters in the negative training set
 ****************************************************************************/
static void start_ele_test_negatives(PS_T *ps, const xmlChar **attrs) {
  long count, positions;

  char* names[2] = {"count", "positions"};
  int (*parsers[2])(char*, void*) = {ld_long, ld_long};
  void *data[2] = {&count, &positions};
  bool required[2] = {true, true};
  bool done[2];

  parse_attributes(streme_attr_parse_error, ps, "test_negatives", attrs, 2, names, parsers, data, required, done);

  if (ps->callbacks->handle_test_negatives && ps->state != PS_ERROR) {
    ps->callbacks->handle_test_negatives(ps->user_data, count, positions);
  }
}

/*****************************************************************************
 * STREME > model > alphabet
 ****************************************************************************/
static void start_ele_alphabet(PS_T *ps, const xmlChar **attrs) {
  char *name;
  int extends;

  char* extends_options[3] = {"dna", "protein", "rna"};
  int extends_values[3] = {ALPH_FLAG_EXTENDS_DNA, ALPH_FLAG_EXTENDS_PROTEIN, ALPH_FLAG_EXTENDS_RNA};
  MULTI_T extends_multi = {.count = 3, .options = extends_options, 
    .outputs = extends_values, .target = &(extends)};

  char* names[2] = {"like", "name"};
  int (*parsers[2])(char*, void*) = {ld_multi, ld_str};
  void *data[2] = {&extends_multi, &name};
  bool required[2] = {false, false};
  bool done[2];
  // just so we know later on when reading the background which used to set the alphabet
  ps->seen_alphabet = true;

  // defaults
  name = NULL;
  extends = 0;
  parse_attributes(streme_attr_parse_error, ps, "alphabet", attrs, 2, names, parsers, data, required, done);

  if (ps->callbacks->start_alphabet && ps->state != PS_ERROR) {
    ps->callbacks->start_alphabet(ps->user_data, name, extends);
  }
  streme_push_es(ps, PS_IN_ALPHABET_LETTER, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * STREME > model > /alphabet
 ****************************************************************************/
static void st_end_ele_alphabet(PS_T *ps) {
  if (ps->callbacks->end_alphabet && ps->state != PS_ERROR) {
    ps->callbacks->end_alphabet(ps->user_data);
  }
}

/*****************************************************************************
 * STREME > model > alphabet > letter
 ****************************************************************************/
static void start_ele_alphabet_letter(PS_T *ps, const xmlChar **attrs) {
  char *aliases, *id, *name, *equals, symbol, complement;
  int colour, idx;

  char* names[7] = {"aliases", "colour", "complement", "equals", "id", "name", "symbol"};
  int (*parsers[7])(char*, void*) = {ld_str, ld_hex, ld_char, ld_str, ld_str, ld_str, ld_char};
  void *data[7] = {&aliases, &colour, &complement, &equals, &id, &name, &symbol};
  bool required[7] = {false, false, false, false, true, false, true};
  bool done[7];

  aliases = NULL;
  name = NULL;
  equals = NULL;
  complement = '\0';
  colour = -1;
  parse_attributes(streme_attr_parse_error, ps, "letter", attrs, 7, names, parsers, data, required, done);

  if (ps->seen_ambig) {
    if (equals == NULL) {
      error(ps, "All core symbols must appear before any ambigous symbols.\n");
    }
  } else if (equals == NULL) {
    idx = rbtree_size(ps->alph_ids);
    rbtree_make(ps->alph_ids, id, &idx);
  } else {
    ps->seen_ambig = true;
  }

  if (ps->callbacks->handle_alphabet_letter && ps->state != PS_ERROR) {
    ps->callbacks->handle_alphabet_letter(ps->user_data, id, symbol, aliases, complement, equals, name, colour);
  }
}

/*****************************************************************************
 * STREME > model > strands\
 ****************************************************************************/
void st_end_ele_strands(PS_T *ps) {
  STREME_STRANDS_EN strands = STREME_STRANDS_FORWARD; // stop compiler complaining
  if (strcmp("both", ps->characters.buffer) == 0) {
    strands = STREME_STRANDS_BOTH;
  } else if (strcmp("forward", ps->characters.buffer) == 0) {
    strands = STREME_STRANDS_FORWARD;
  } else if (strcmp("none", ps->characters.buffer) == 0) {
    strands = STREME_STRANDS_NONE;
  } else {
    error(ps, "Strands value must be both, forward or none.\n");
  }
  if (ps->callbacks->handle_strands && ps->state != PS_ERROR) {
    ps->callbacks->handle_strands(ps->user_data, strands);
  }
}

/*****************************************************************************
 * STREME > model > background	(pre 5.3.0)
 *
 *  type            is the alphabet DNA or RNA (optional when custom alphabet specified)
 *  <symbol>        frequency of <symbol> from core alphabet
 *  from            from the negative dataset or a background file
 *  file            the background file (optional)
 *  last_mod_date   the last modified date of the background file (optional)
 ****************************************************************************/
static void start_ele_background(PS_T *ps, const xmlChar **attrs) {
  int from;

  char* from_options[1] = {"dataset"};
  int from_values[1] = {STREME_BG_FROM_DATASET};
  MULTI_T from_multi = {.count = 1, .options = from_options,
    .outputs = from_values, .target = &(from)};

  char* names[1] = {"from"};
  int (*parsers[1])(char*, void*) = {ld_multi};
  void *data[1] = {&from_multi};
  bool required[1] = {true};
  bool done[1];

  parse_attributes(streme_attr_parse_error, ps, "background", attrs, 1, names, parsers, data, required, done);
  parse_freq_attrs(ps, "background", attrs);

  if (ps->callbacks->handle_background && ps->state != PS_ERROR) {
    ps->callbacks->handle_background(ps->user_data, rbtree_size(ps->alph_ids), ps->freqs, from);
  }
}

/*****************************************************************************
 * STREME > model > background_frequencies
 *
 *  type            is the alphabet DNA or RNA (optional when custom alphabet specified)
 *  <symbol>        frequency of <symbol> from core alphabet
 *  from            from the negative dataset or a background file
 *  file            the background file (optional)
 *  last_mod_date   the last modified date of the background file (optional)
 ****************************************************************************/
static void start_ele_background_frequencies(PS_T *ps, const xmlChar **attrs) {
  char *bgsource;

  char* names[1] = {"source"};
  int (*parsers[1])(char*, void*) = {ld_str};
  void *data[1] = {&bgsource};
  bool required[1] = {true};
  bool done[1];

  parse_attributes(streme_attr_parse_error, ps, "background_frequencies", attrs, 1, names, parsers, data, required, done);

  if (ps->callbacks->start_background_frequencies && ps->state != PS_ERROR) {
    ps->callbacks->start_background_frequencies(ps->user_data, bgsource);
  }
  streme_push_es(ps, PS_IN_BF_ALPHABET_ARRAY, ES_ONCE);
}

/*****************************************************************************
 * STREME > model > background_frequencies
 ****************************************************************************/
static void end_ele_background_frequencies(PS_T *ps) {
  if (ps->callbacks->end_background_frequencies && ps->state != PS_ERROR) {
    ps->callbacks->end_background_frequencies(ps->user_data);
  }
}

/*****************************************************************************
 * STREME > model > background_frequencies > alphabet_array
 ****************************************************************************/
static void start_ele_bf_alphabet_array(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_bf_alphabet_array && ps->state != PS_ERROR) {
    ps->callbacks->start_bf_alphabet_array(ps->user_data);
  }
  streme_push_es(ps, PS_IN_BF_AA_VALUE, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * STREME > model > background_frequencies > /alphabet_array
 ****************************************************************************/
static void end_ele_bf_alphabet_array(PS_T *ps) {
  if (ps->callbacks->end_bf_alphabet_array && ps->state != PS_ERROR) {
    ps->callbacks->end_bf_alphabet_array(ps->user_data);
  }
}

/*****************************************************************************
 * STREME > model > background_frequencies > alphabet_array > value
 ****************************************************************************/
static void start_ele_value(PS_T *ps, const xmlChar **attrs) {
  char *letter_id;
  int letter_id_len;

  char* names[1] = {"letter_id"};
  int (*parsers[1])(char*, void*) = {ld_str};
  void *data[1] = {&letter_id};
  bool required[1] = {true};
  bool done[1];

  parse_attributes(streme_attr_parse_error, ps, "value", attrs, 1, names, parsers, data, required, done);

  letter_id_len = strlen(letter_id);
  if (ps->letter_id_buf_len <= letter_id_len) {
    ps->letter_id_buf = mm_realloc(ps->letter_id_buf, sizeof(char) * (letter_id_len + 1));
    ps->letter_id_buf_len = letter_id_len + 1;
  }
  strncpy(ps->letter_id_buf, letter_id, ps->letter_id_buf_len);
}

/*****************************************************************************
 * STREME > model > background_frequencies > alphabet_array > /value
 ****************************************************************************/
static void end_ele_bf_aa_value(PS_T *ps) {
  double value;
  if (ld_double(ps->characters.buffer, &value)) {
    error(ps, "Couldn't parse value from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_bf_aa_value && ps->state != PS_ERROR) {
    ps->callbacks->handle_bf_aa_value(ps->user_data, ps->letter_id_buf, value);
  }
}

/*****************************************************************************
 * STREME > model > stop
 *
 *  pvt		the maximum p-value threshold for motifs
 *  nmotifs	the maximum number of motifs to find
 *  time 	the maximum time to run (seconds)
 ****************************************************************************/
static void start_ele_stop(PS_T *ps, const xmlChar **attrs) {
  double pvt, time;
  long nmotifs;

  char* names[3] = {"pvt", "nmotifs", "time"};
  int (*parsers[3])(char*, void*) = {ld_double, ld_long, ld_double};
  void *data[3] = {&pvt, &nmotifs, &time};
  bool required[3] = {false, false, false};
  bool done[3];

  parse_attributes(streme_attr_parse_error, ps, "stop", attrs, 3, names, parsers, data, required, done);

  if (ps->callbacks->handle_stop && ps->state != PS_ERROR) {
    ps->callbacks->handle_stop(ps->user_data, pvt, nmotifs, time);
  }
}

/*****************************************************************************
 * STREME > model > /objfun
 ****************************************************************************/
void st_end_ele_objfun(PS_T *ps) {
  STREME_OBJFUN_EN objfun = STREME_DIFFERENTIAL_ENRICHMENT; // stop compiler complaining
  if (strcmp("Differential Enrichment", ps->characters.buffer) == 0) {
    objfun = STREME_DIFFERENTIAL_ENRICHMENT;
  } else if (strcmp("Central Distance", ps->characters.buffer) == 0) {
    objfun = STREME_CENTRAL_DISTANCE;
  } else {
    error(ps, "Objfun value must be Differential Enrichment or Central Distance.\n");
  }
  if (ps->callbacks->handle_objfun && ps->state != PS_ERROR) {
    ps->callbacks->handle_objfun(ps->user_data, objfun);
  }
}

/*****************************************************************************
 * STREME > model > /test
 ****************************************************************************/
void st_end_ele_test(PS_T *ps) {
  STREME_TEST_EN test = STREME_FISHER_EXACT_TEST; // stop compiler complaining
  if (strcmp("Fisher Exact Test", ps->characters.buffer) == 0) {
    test = STREME_FISHER_EXACT_TEST;
  } else if (strcmp("Binomial Test", ps->characters.buffer) == 0) {
    test = STREME_BINOMIAL_TEST;
  } else if (strcmp("Cumulative Bates Distribution", ps->characters.buffer) == 0) {
    test = STREME_CUMULATIVE_BATES_DISTRIBUTION;
  } else {
    error(ps, "Test value must be Fisher Exact Test, Binomial Test or Cumulative Bates Distribution.\n");
  }
  if (ps->callbacks->handle_test && ps->state != PS_ERROR) {
    ps->callbacks->handle_test(ps->user_data, test);
  }
}

/*****************************************************************************
 * STREME > model > /minw
 ****************************************************************************/
static void st_end_ele_minw(PS_T *ps) {
  int minw;

  if (ld_int(ps->characters.buffer, &minw)) {
    error(ps, "Bad value \"%s\" for minw.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_minw && ps->state != PS_ERROR) {
    ps->callbacks->handle_minw(ps->user_data, minw);
  }
}

/*****************************************************************************
 * STREME > model > /maxw
 ****************************************************************************/
static void st_end_ele_maxw(PS_T *ps) {
  int maxw;

  if (ld_int(ps->characters.buffer, &maxw)) {
    error(ps, "Bad value \"%s\" for maxw.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_maxw && ps->state != PS_ERROR) {
    ps->callbacks->handle_maxw(ps->user_data, maxw);
  }
}

/*****************************************************************************
 * STREME > model > /kmer
 ****************************************************************************/
static void st_end_ele_kmer(PS_T *ps) {
  int kmer;

  if (ld_int(ps->characters.buffer, &kmer)) {
    error(ps, "Bad value \"%s\" for kmer.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_kmer && ps->state != PS_ERROR) {
    ps->callbacks->handle_kmer(ps->user_data, kmer);
  }
}

/*****************************************************************************
 * STREME > model > /hofract
 ****************************************************************************/
static void st_end_ele_hofract(PS_T *ps) {
  double hofract;

  if (ld_double(ps->characters.buffer, &hofract)) {
    error(ps, "Bad value \"%s\" for hofract.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_hofract && ps->state != PS_ERROR) {
    ps->callbacks->handle_hofract(ps->user_data, hofract);
  }
}

/*****************************************************************************
 * STREME > model > /neval
 ****************************************************************************/
static void st_end_ele_neval(PS_T *ps) {
  int neval;

  if (ld_int(ps->characters.buffer, &neval)) {
    error(ps, "Bad value \"%s\" for neval.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_neval && ps->state != PS_ERROR) {
    ps->callbacks->handle_neval(ps->user_data, neval);
  }
}

/*****************************************************************************
 * STREME > model > /nref
 ****************************************************************************/
static void st_end_ele_nref(PS_T *ps) {
  int nref;

  if (ld_int(ps->characters.buffer, &nref)) {
    error(ps, "Bad value \"%s\" for nref.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_nref && ps->state != PS_ERROR) {
    ps->callbacks->handle_nref(ps->user_data, nref);
  }
}

/*****************************************************************************
 * STREME > model > /niter
 ****************************************************************************/
static void st_end_ele_niter(PS_T *ps) {
  int niter;

  if (ld_int(ps->characters.buffer, &niter)) {
    error(ps, "Bad value \"%s\" for niter.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_niter && ps->state != PS_ERROR) {
    ps->callbacks->handle_niter(ps->user_data, niter);
  }
}

/*****************************************************************************
 * STREME > model > patience
 ****************************************************************************/
static void st_end_ele_patience(PS_T *ps) {
  long patience;

  if (ld_long(ps->characters.buffer, &patience)) {
    error(ps, "Bad value \"%s\" for patience.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_patience && ps->state != PS_ERROR) {
    ps->callbacks->handle_patience(ps->user_data, patience);
  }
}

/*****************************************************************************
 * STREME > model > /seed
 ****************************************************************************/
static void st_end_ele_seed(PS_T *ps) {
  long seed;

  if (ld_long(ps->characters.buffer, &seed)) {
    error(ps, "Bad value \"%s\" for seed.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_seed && ps->state != PS_ERROR) {
    ps->callbacks->handle_seed(ps->user_data, seed);
  }
}

/*****************************************************************************
 * STREME > model > /useer
 ****************************************************************************/
static void st_end_ele_useer(PS_T *ps) {
  bool useer = false;

  if (strcmp("no", ps->characters.buffer) == 0) {
    useer = false;
  } else if (strcmp("yes", ps->characters.buffer) == 0) {
    useer = true;
  } else {
    error(ps, "Bad value \"%s\" for useer.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_useer && ps->state != PS_ERROR) {
    ps->callbacks->handle_useer(ps->user_data, useer);
  }
}

/*****************************************************************************
 * STREME > model > /minscore
 ****************************************************************************/
static void st_end_ele_minscore(PS_T *ps) {
  double minscore;

  if (ld_double(ps->characters.buffer, &minscore)) {
    error(ps, "Bad value \"%s\" for minscre.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_minscore && ps->state != PS_ERROR) {
    ps->callbacks->handle_minscore(ps->user_data, minscore);
  }
}

/*****************************************************************************
 * STREME > model > /ignore_depth
 ****************************************************************************/
static void st_end_ele_ignore_depth(PS_T *ps) {
  long ignore_depth;

  if (ld_long(ps->characters.buffer, &ignore_depth)) {
    error(ps, "Bad value \"%s\" for ignore_depth.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_ignore_depth && ps->state != PS_ERROR) {
    ps->callbacks->handle_ignore_depth(ps->user_data, ignore_depth);
  }
}

/*****************************************************************************
 * STREME > model > /nsubsets
 ****************************************************************************/
static void st_end_ele_nsubsets(PS_T *ps) {
  long nsubsets;

  if (ld_long(ps->characters.buffer, &nsubsets)) {
    error(ps, "Bad value \"%s\" for nsubsets.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_nsubsets && ps->state != PS_ERROR) {
    ps->callbacks->handle_nsubsets(ps->user_data, nsubsets);
  }
}

/*****************************************************************************
 * STREME > model > /min_pal_ratio
 ****************************************************************************/
static void st_end_ele_min_pal_ratio(PS_T *ps) {
  double min_pal_ratio;

  if (ld_double(ps->characters.buffer, &min_pal_ratio)) {
    error(ps, "Bad value \"%s\" for min_pal_ratio.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_min_pal_ratio && ps->state != PS_ERROR) {
    ps->callbacks->handle_min_pal_ratio(ps->user_data, min_pal_ratio);
  }
}

/*****************************************************************************
 * STREME > model > /max_pal_ed
 ****************************************************************************/
static void st_end_ele_max_pal_ed(PS_T *ps) {
  double max_pal_ed;

  if (ld_double(ps->characters.buffer, &max_pal_ed)) {
    error(ps, "Bad value \"%s\" for max_pal_ed.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_max_pal_ed && ps->state != PS_ERROR) {
    ps->callbacks->handle_max_pal_ed(ps->user_data, max_pal_ed);
  }
}

/*****************************************************************************
 * STREME > model > /cand
 ****************************************************************************/
static void st_end_ele_cand(PS_T *ps) {
  bool cand = false;

  if (strcmp("no", ps->characters.buffer) == 0) {
    cand = false;
  } else if (strcmp("yes", ps->characters.buffer) == 0) {
    cand = true;
  } else {
    error(ps, "Bad value \"%s\" for cand.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_cand && ps->state != PS_ERROR) {
    ps->callbacks->handle_cand(ps->user_data, cand);
  }
}

/*****************************************************************************
 * STREME > model > /experimental
 ****************************************************************************/
static void st_end_ele_experimental(PS_T *ps) {
  bool experimental = false;

  if (strcmp("no", ps->characters.buffer) == 0) {
    experimental = false;
  } else if (strcmp("yes", ps->characters.buffer) == 0) {
    experimental = true;
  } else {
    error(ps, "Bad value \"%s\" for experimental.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_experimental && ps->state != PS_ERROR) {
    ps->callbacks->handle_experimental(ps->user_data, experimental);
  }
}

/*****************************************************************************
 * STREME > model > /totallength
 ****************************************************************************/
static void st_end_ele_totallength(PS_T *ps) {
  long totallength;

  if (ld_long(ps->characters.buffer, &totallength)) {
    error(ps, "Bad value \"%s\" for totallength.\n", ps->characters.buffer);
  }

  if (ps->callbacks->handle_totallength && ps->state != PS_ERROR) {
    ps->callbacks->handle_totallength(ps->user_data, totallength);
  }
}

/*****************************************************************************
 * STREME > model > /align
 ****************************************************************************/
static void st_end_ele_align(PS_T *ps) {
  if (ps->callbacks->handle_host && ps->state != PS_ERROR) {
    ps->callbacks->handle_host(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * STREME > model > /host
 * the name of the computer which ran STREME.
 ****************************************************************************/
static void st_end_ele_host(PS_T *ps) {
  if (ps->callbacks->handle_host && ps->state != PS_ERROR) {
    ps->callbacks->handle_host(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * STREME > model > /description
 * an optional description of the experiment.
 ****************************************************************************/
static void end_ele_description(PS_T *ps) {
  if (ps->callbacks->handle_description && ps->state != PS_ERROR) {
    ps->callbacks->handle_description(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * STREME > motifs
 * contains: motif*
 ****************************************************************************/
static void start_ele_motifs(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_motifs && ps->state != PS_ERROR) {
    ps->callbacks->start_motifs(ps->user_data);
  }
  streme_push_es(ps, PS_IN_MOTIF, ES_ANY);
}

/*****************************************************************************
 * STREME > /motifs
 ****************************************************************************/
static void st_end_ele_motifs(PS_T *ps) {
  if (ps->callbacks->end_motifs && ps->state != PS_ERROR) {
    ps->callbacks->end_motifs(ps->user_data);
  }
}

/*****************************************************************************
 * STREME > motifs > motif
 *
 *  id                the identifier used by STREME
 *  alt		      alternate identifier used by STREME
 *  width	      the width of the motif
 *  initial_width     the width of the prefix seed for the motif
 *  seed 	      the seed word for the motif
 *  score_threshold   the log-odds score threshold for calling sites
 *  train_pos_count   positive training sites
 *  train_neg_count   negative training sites
 *  train_log_pvalue  log10 p-value on training set
 *  train_pvalue      p-value on training set
 *  train_dtc         CD only: average site distance to center in training set
 *  train_bernoulli   Binomial Test only: Bernoulli prob. in training set
 *  test_pos_count    positive test sites
 *  test_neg_count    negative test sites
 *  test_log_pvalue   log10 p-value on test set
 *  test_pvalue       p-value on test set
 *  test_dtc          CD only: average site distance to center in test set
 *  test_bernoulli    Binomial Test only: Bernoulli prob. in test set
 *  is_palindromic    is the motif a perfect palindrome?
 *  elapsed time      time elapsed until motif found (seconds)
 ****************************************************************************/
static void start_ele_motif(PS_T *ps, const xmlChar **attrs) {
  char *id, *alt, *seed;
  int width, initial_width;
  double score_threshold, elapsed_time;
  long train_pos_count, train_neg_count;
  double train_log_pvalue, train_dtc, train_bernoulli;
  char *train_pvalue;
  long test_pos_count, test_neg_count;
  double test_log_pvalue, test_dtc, test_bernoulli;
  char *test_pvalue;
  char *is_palindromic;

  // must be listed alphabetically
  char* names[20] = {
    "alt", "elapsed_time", "id", "initial_width", "is_palindromic", "score_threshold", "seed", 
    "test_bernoulli", "test_dtc", "test_log_pvalue", "test_neg_count", "test_pos_count", "test_pvalue", 
    "train_bernoulli", "train_dtc", "train_log_pvalue", "train_neg_count", "train_pos_count", "train_pvalue", 
    "width"
  };
  int (*parsers[20])(char*, void*) = {
    ld_str, ld_double, ld_str, ld_int, ld_str, ld_double, ld_str,
    ld_double, ld_double, ld_double, ld_long, ld_long, ld_str,
    ld_double, ld_double, ld_double, ld_long, ld_long, ld_str,
    ld_int
  };
  void *data[20] = {
    &alt, &elapsed_time, &id, &initial_width, &is_palindromic, &score_threshold, &seed, 
    &test_bernoulli, &test_dtc, &test_log_pvalue, &test_neg_count, &test_pos_count, &test_pvalue, 
    &train_bernoulli, &train_dtc, &train_log_pvalue, &train_neg_count, &train_pos_count, &train_pvalue, 
    &width
  };
  bool required[20] = {
    true, true, true, true, true, true, true, 
    true, true, true, true, true, true, 
    true, true, true, true, true, true, 
    true
  };
  bool done[20];

  parse_attributes(streme_attr_parse_error, ps, "motif", attrs, 20, names, 
      parsers, data, required, done);

  // copy the motif id so we can use it in any error messages
  if (ps->state != PS_ERROR) {
    int len = strlen(id);
    ps->motif_id = mm_malloc(sizeof(char) * (len + 1));
    strcpy(ps->motif_id, id);
    ps->last_pos = 0;
    ps->motif_width = width;
  }

  if (ps->callbacks->start_motif && ps->state != PS_ERROR) {
    ps->callbacks->start_motif(ps->user_data, 
      id, alt, width, initial_width, seed, score_threshold, 
      train_pos_count, train_neg_count, train_log_pvalue, train_pvalue, train_dtc, train_bernoulli,
      test_pos_count, test_neg_count, test_log_pvalue, test_pvalue, test_dtc, test_bernoulli,
      is_palindromic, elapsed_time);
  }
  streme_push_es(ps, PS_IN_POS, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * STREME > motifs > /motif
 * contains: pos+
 ****************************************************************************/
static void st_end_ele_motif(PS_T *ps) {
  if (ps->state != PS_ERROR && ps->last_pos != ps->motif_width) {
    error(ps, "Motif %s has length %d which does not match the "
        "number of pos elements %d.", ps->motif_width, ps->last_pos);
  }

  if (ps->callbacks->end_motif && ps->state != PS_ERROR) {
    ps->callbacks->end_motif(ps->user_data);
  }
  // free stored motif id
  if (ps->motif_id) free(ps->motif_id);
  ps->motif_id = NULL;
}

/*****************************************************************************
 * STREME > motifs > motif > pos
 *
 *  i                 index of the motif position (optional)
 *  <symbol>          frequency of <symbol>
 ****************************************************************************/
static void start_ele_pos(PS_T *ps, const xmlChar **attrs) {
  int pos = ps->last_pos + 1;
  if (!ps->seen_alphabet) {
    // attribute "i" only exists in the older specification without custom alphabets
    char* names[1] = {"i"};
    int (*parsers[1])(char*, void*) = {ld_int};
    void *data[1] = {&pos};
    bool required[1] = {true};
    bool done[1];
    parse_attributes(streme_attr_parse_error, ps, "pos", attrs, 1, names, parsers, data, required, done);
    if (pos != (ps->last_pos + 1)) {
      error(ps, "Motif %s did not have pos %d but instead "
          "has pos %d in its place.\n", ps->motif_id, ps->last_pos + 1, pos);
    }
  }
  ps->last_pos = pos;
  parse_freq_attrs(ps, "pos", attrs);
  if (ps->callbacks->handle_pos && ps->state != PS_ERROR) {
    ps->callbacks->handle_pos(ps->user_data, pos, rbtree_size(ps->alph_ids), ps->freqs);
  }
}

/*****************************************************************************
 * Handle the document start
 ****************************************************************************/
void handle_streme_start_doc(void *ctx) {
  PS_T *ps = (PS_T*)ctx;
  ps->state = PS_START;
  streme_push_es(ps, PS_END, ES_ONCE);
  streme_push_es(ps, PS_IN_STREME, ES_ONCE);
}

/*****************************************************************************
 * Handle the document end
 ****************************************************************************/
void handle_streme_end_doc(void *ctx) {
  PS_T *ps = (PS_T*)ctx;
  ES_T *es;
  if (ps->state != PS_ERROR) streme_update_es(ps, PS_END);
  while ((es = (ES_T*)linklst_pop(ps->expected_stack)) != NULL) {
    free(es);
  }
}

/*****************************************************************************
 * Handle the characters within elements. Only accepts characters
 * that pass the character filter. If a filter is not set then it doesn't
 * store anything.
 ****************************************************************************/
void handle_streme_characters(void *ctx, const xmlChar *ch, int len) {
  PS_T *ps = (PS_T*)ctx;
  if (ps->state == PS_ERROR) return;
  if (ps->udepth) {//we don't know where we are!
  } else {
    store_xml_characters(&(ps->characters), (const char*)ch, len);
  }
}

/*****************************************************************************
 * Macros to greatly reduce the boilerplate code needed to keep track of the
 * parser state.
 *
 * DO_START_ELE compares the tag name to a possible option and if it matches
 * then validates the state stack to ensure that it really is a possible
 * option. If everything passes then it changes the state to the specified
 * transition state and calls the callback. It also changes the character
 * accept function.
 * CHECK_START_ELE does the same except it does not have a callback function.
 ****************************************************************************/
#define DO_START_ELE(_expected_,_call_,_transition_,_char_accept_) \
  if (strcasecmp((char*)name, #_expected_ ) == 0) { \
    if (streme_update_es(ps, _transition_ )) { \
      ps->state = _transition_; \
      _call_ (ps, translate_attributes(&(ps->attrbuf), attrs)); \
      ps->characters.accept = _char_accept_; \
    } \
    break; \
  }

#define CHECK_START_ELE(_expected_,_transition_,_char_accept_) \
  if (strcasecmp((char*)name, #_expected_ ) == 0) { \
    if (streme_update_es(ps, _transition_ )) { \
      ps->state = _transition_; \
      ps->characters.accept = _char_accept_; \
    } \
    break; \
  }

/*****************************************************************************
 * Handle a starting element.
 ****************************************************************************/
void handle_streme_start_ele(void *ctx, const xmlChar *name, const xmlChar **attrs) {
  PS_T *ps = (PS_T*)ctx;
  int known, allow_unknown;
  if (ps->state == PS_ERROR) return;
  if (ps->udepth) {//we don't know where we are!
    ps->udepth += 1;
  } else {
    //reset the character buffer to the begining
    ps->characters.pos = 0;
    ps->characters.buffer[0] = '\0';
    known = 1; //assume we can find it
    allow_unknown = 1; //are unknowns allowed in this state?
    switch (ps->state) {
      case PS_START:
        DO_START_ELE(streme, start_ele_streme, PS_IN_STREME, IGNORE);
        known = 0;
        break;
      case PS_IN_STREME:
        DO_START_ELE(model, start_ele_model, PS_IN_MODEL, IGNORE);
        DO_START_ELE(motifs, start_ele_motifs, PS_IN_MOTIFS, IGNORE);
        CHECK_START_ELE(reason_for_stopping, PS_IN_REASON_FOR_STOPPING, IGNORE);
        CHECK_START_ELE(run_time, PS_IN_RUN_TIME, IGNORE);
        known = 0;
        break;
      case PS_IN_MODEL:
        CHECK_START_ELE(command_line, PS_IN_COMMAND_LINE, ALL_CHARS);
        DO_START_ELE(train_positives, start_ele_train_positives, PS_IN_TRAIN_POSITIVES, IGNORE);
        DO_START_ELE(train_negatives, start_ele_train_negatives, PS_IN_TRAIN_NEGATIVES, IGNORE);
        DO_START_ELE(test_positives, start_ele_test_positives, PS_IN_TEST_POSITIVES, IGNORE);
        DO_START_ELE(test_negatives, start_ele_test_negatives, PS_IN_TEST_NEGATIVES, IGNORE);
        DO_START_ELE(alphabet, start_ele_alphabet, PS_IN_ALPHABET, IGNORE);
        CHECK_START_ELE(strands, PS_IN_STRANDS, ALL_BUT_SPACE);
        DO_START_ELE(background, start_ele_background, PS_IN_BACKGROUND, IGNORE);
        CHECK_START_ELE(sequence_db, PS_IN_SEQUENCE_DB, IGNORE);
        DO_START_ELE(background_frequencies, start_ele_background_frequencies, PS_IN_BACKGROUND_FREQUENCIES, IGNORE);
        DO_START_ELE(stop, start_ele_stop, PS_IN_STOP, IGNORE);
        CHECK_START_ELE(objfun, PS_IN_OBJFUN, ALL_CHARS);
        CHECK_START_ELE(test, PS_IN_TEST, ALL_CHARS);
        CHECK_START_ELE(minw, PS_IN_MINW, ALL_BUT_SPACE);
        CHECK_START_ELE(maxw, PS_IN_MAXW, ALL_BUT_SPACE);
        CHECK_START_ELE(kmer, PS_IN_KMER, ALL_BUT_SPACE);
        CHECK_START_ELE(hofract, PS_IN_HOFRACT, ALL_BUT_SPACE);
        CHECK_START_ELE(neval, PS_IN_NEVAL, ALL_BUT_SPACE);
        CHECK_START_ELE(nref, PS_IN_NREF, ALL_BUT_SPACE);
        CHECK_START_ELE(niter, PS_IN_NITER, ALL_BUT_SPACE);
        CHECK_START_ELE(patience, PS_IN_PATIENCE, ALL_BUT_SPACE);
        CHECK_START_ELE(seed, PS_IN_SEED, ALL_BUT_SPACE);
        CHECK_START_ELE(useer, PS_IN_USEER, ALL_BUT_SPACE);
        CHECK_START_ELE(minscore, PS_IN_MINSCORE, ALL_BUT_SPACE);
        CHECK_START_ELE(ignore_depth, PS_IN_IGNORE_DEPTH, ALL_BUT_SPACE);
        CHECK_START_ELE(nsubsets, PS_IN_NSUBSETS, ALL_BUT_SPACE);
        CHECK_START_ELE(min_pal_ratio, PS_IN_MIN_PAL_RATIO, ALL_BUT_SPACE);
        CHECK_START_ELE(max_pal_ed, PS_IN_MAX_PAL_ED, ALL_BUT_SPACE);
        CHECK_START_ELE(cand, PS_IN_CAND, ALL_BUT_SPACE);
        CHECK_START_ELE(experimental, PS_IN_EXPERIMENTAL, ALL_BUT_SPACE);
        CHECK_START_ELE(totallength, PS_IN_TOTALLENGTH, ALL_BUT_SPACE);
        CHECK_START_ELE(align, PS_IN_ALIGN, ALL_BUT_SPACE);
        CHECK_START_ELE(host, PS_IN_HOST, ALL_BUT_SPACE);
        CHECK_START_ELE(description, PS_IN_DESCRIPTION, ALL_CHARS);
        break;
      case PS_IN_ALPHABET:
        DO_START_ELE(letter, start_ele_alphabet_letter, PS_IN_ALPHABET_LETTER, IGNORE);
        break;
      case PS_IN_MOTIFS:
        DO_START_ELE(motif, start_ele_motif, PS_IN_MOTIF, IGNORE);
        known = 0;
        break;
      case PS_IN_MOTIF:
        DO_START_ELE(pos, start_ele_pos, PS_IN_POS, IGNORE);
        known = 0;
        break;
      case PS_IN_BACKGROUND_FREQUENCIES:
        DO_START_ELE(alphabet_array, start_ele_bf_alphabet_array, PS_IN_BF_ALPHABET_ARRAY, IGNORE);
        known = 0;
        break;
      case PS_IN_BF_ALPHABET_ARRAY:
        DO_START_ELE(value, start_ele_value, PS_IN_BF_AA_VALUE, ALL_BUT_SPACE);
        known = 0;
        break;
      default:
        known = 0;
        break;
    }
    if (!known) {
      if (allow_unknown) {
        ps->udepth = 1;
      } else {
        error(ps, "Encountered illegal tag %s while in state %s\n",(char*)name, state_names[ps->state]);
      }
    }
  }
}

/*****************************************************************************
 * Macros to greatly reduce the boilerplate code needed to keep track of the
 * parser state.
 *
 * DO_END_ELE compares the state and if it matches then tries to match the
 * element name. If that passes then it calls the callback and if that 
 * succeeds it transitions to the parent state.
 * CHECK_END_ELE does the same except it does not have a callback function.
 ****************************************************************************/
#define DO_END_ELE(_expected_,_call_,_state_,_transition_) \
  case _state_: \
    if (strcasecmp((char*)name, #_expected_) == 0) { \
      known = 1; \
      _call_ (ps); \
      if (ps->state == _state_) ps->state = _transition_; \
    } \
    break
#define CHECK_END_ELE(_expected_,_state_,_transition_) \
  case _state_: \
    if (strcasecmp((char*)name, #_expected_) == 0) { \
      known = 1; \
      ps->state = _transition_; \
    } \
    break

/*****************************************************************************
 * Handle a closing element
 ****************************************************************************/
void handle_streme_end_ele(void *ctx, const xmlChar *name) {
  PS_T *ps = (PS_T*)ctx;
  int known = 0;
  if (ps->state == PS_ERROR) return;
  if (ps->udepth) {
    ps->udepth -= 1; 
  } else {
    switch (ps->state) {
      DO_END_ELE(streme, st_end_ele_streme, PS_IN_STREME, PS_END);
      DO_END_ELE(model, st_end_ele_model, PS_IN_MODEL, PS_IN_STREME);
      DO_END_ELE(command_line, st_end_ele_command_line, PS_IN_COMMAND_LINE, PS_IN_MODEL);
      CHECK_END_ELE(train_positives, PS_IN_TRAIN_POSITIVES, PS_IN_MODEL);
      CHECK_END_ELE(train_negatives, PS_IN_TRAIN_NEGATIVES, PS_IN_MODEL);
      CHECK_END_ELE(test_positives, PS_IN_TEST_POSITIVES, PS_IN_MODEL);
      CHECK_END_ELE(test_negatives, PS_IN_TEST_NEGATIVES, PS_IN_MODEL);
      DO_END_ELE(alphabet, st_end_ele_alphabet, PS_IN_ALPHABET, PS_IN_MODEL);
      CHECK_END_ELE(letter, PS_IN_ALPHABET_LETTER, PS_IN_ALPHABET);
      CHECK_END_ELE(strands, PS_IN_STRANDS, PS_IN_MODEL);
      CHECK_END_ELE(background, PS_IN_BACKGROUND, PS_IN_MODEL);
      CHECK_END_ELE(sequence_db, PS_IN_SEQUENCE_DB, PS_IN_MODEL);
      DO_END_ELE(background_frequencies, end_ele_background_frequencies, PS_IN_BACKGROUND_FREQUENCIES, PS_IN_MODEL);
      DO_END_ELE(alphabet_array, end_ele_bf_alphabet_array, PS_IN_BF_ALPHABET_ARRAY, PS_IN_BACKGROUND_FREQUENCIES);
      DO_END_ELE(value, end_ele_bf_aa_value, PS_IN_BF_AA_VALUE, PS_IN_BF_ALPHABET_ARRAY);
      CHECK_END_ELE(stop, PS_IN_STOP, PS_IN_MODEL);
      DO_END_ELE(objfun, st_end_ele_objfun, PS_IN_OBJFUN, PS_IN_MODEL);
      DO_END_ELE(test, st_end_ele_test, PS_IN_TEST, PS_IN_MODEL);
      DO_END_ELE(minw, st_end_ele_minw, PS_IN_MINW, PS_IN_MODEL);
      DO_END_ELE(maxw, st_end_ele_maxw, PS_IN_MAXW, PS_IN_MODEL);
      DO_END_ELE(kmer, st_end_ele_kmer, PS_IN_KMER, PS_IN_MODEL);
      DO_END_ELE(hofract, st_end_ele_hofract, PS_IN_HOFRACT, PS_IN_MODEL);
      DO_END_ELE(neval, st_end_ele_neval, PS_IN_NEVAL, PS_IN_MODEL);
      DO_END_ELE(nref, st_end_ele_nref, PS_IN_NREF, PS_IN_MODEL);
      DO_END_ELE(niter, st_end_ele_niter, PS_IN_NITER, PS_IN_MODEL);
      DO_END_ELE(patience, st_end_ele_patience, PS_IN_PATIENCE, PS_IN_MODEL);
      DO_END_ELE(seed, st_end_ele_seed, PS_IN_SEED, PS_IN_MODEL);
      DO_END_ELE(useer, st_end_ele_useer, PS_IN_USEER, PS_IN_MODEL);
      DO_END_ELE(minscore, st_end_ele_minscore, PS_IN_MINSCORE, PS_IN_MODEL);
      DO_END_ELE(ignore_depth, st_end_ele_ignore_depth, PS_IN_IGNORE_DEPTH, PS_IN_MODEL);
      DO_END_ELE(nsubsets, st_end_ele_nsubsets, PS_IN_NSUBSETS, PS_IN_MODEL);
      DO_END_ELE(min_pal_ratio, st_end_ele_min_pal_ratio, PS_IN_MIN_PAL_RATIO, PS_IN_MODEL);
      DO_END_ELE(max_pal_ed, st_end_ele_max_pal_ed, PS_IN_MAX_PAL_ED, PS_IN_MODEL);
      DO_END_ELE(cand, st_end_ele_cand, PS_IN_CAND, PS_IN_MODEL);
      DO_END_ELE(experimental, st_end_ele_experimental, PS_IN_EXPERIMENTAL, PS_IN_MODEL);
      DO_END_ELE(totallength, st_end_ele_totallength, PS_IN_TOTALLENGTH, PS_IN_MODEL);
      DO_END_ELE(align, st_end_ele_align, PS_IN_ALIGN, PS_IN_MODEL);
      DO_END_ELE(host, st_end_ele_host, PS_IN_HOST, PS_IN_MODEL);
      DO_END_ELE(description, end_ele_description, PS_IN_DESCRIPTION, PS_IN_MODEL);
      DO_END_ELE(motifs, st_end_ele_motifs, PS_IN_MOTIFS, PS_IN_STREME);
      DO_END_ELE(motif, st_end_ele_motif, PS_IN_MOTIF, PS_IN_MOTIFS);
      CHECK_END_ELE(pos, PS_IN_POS, PS_IN_MOTIF);
      CHECK_END_ELE(reason_for_stopping, PS_IN_REASON_FOR_STOPPING, PS_IN_STREME);
      CHECK_END_ELE(run_time, PS_IN_RUN_TIME, PS_IN_STREME);
    }
    if (!known) {
      error(ps, "Unexpected end tag %s in state %s\n", (char*)name, state_names[ps->state]);
    }
  }
}

/*****************************************************************************
 * Register handlers on the xmlSAXHandler structure
 ****************************************************************************/
void register_streme_io_xml_sax_handlers(xmlSAXHandler *handler) {
  memset(handler, 0, sizeof(xmlSAXHandler));
  handler->startDocument = handle_streme_start_doc;
  handler->endDocument = handle_streme_end_doc;
  handler->characters = handle_streme_characters;
  handler->startElement = handle_streme_start_ele;
  handler->endElement = handle_streme_end_ele;
}

/*****************************************************************************
 * Creates the data to be passed to the SAX parser
 ****************************************************************************/
void* create_streme_io_xml_sax_context(void *user_data, STREME_IO_XML_CALLBACKS_T *callbacks) {
  PS_T *ps;
  CHARBUF_T *buf;
  ps = (PS_T*)mm_malloc(sizeof(PS_T));
  memset(ps, 0, sizeof(PS_T));
  ps->state = PS_START;
  ps->udepth = 0; 
  ps->callbacks = callbacks;
  ps->user_data = user_data;
  ps->seen_alphabet = false;
  ps->seen_ambig = false;
  ps->alph_ids = rbtree_create(rbtree_strcmp, rbtree_strcpy, free, rbtree_intcpy, free);
  ps->freqs = NULL;
  ps->motif_id = NULL;
  ps->last_pos = 0;
  // setup attribute buffer
  attrbuf_init(&(ps->attrbuf));
  //set up character buffer
  buf = &(ps->characters);
  buf->buffer = mm_malloc(sizeof(char)*10);
  buf->buffer[0] = '\0';
  buf->size = 10;
  buf->pos = 0;
  //set up attribute buffer (needed for an alphabet array value)
  ps->letter_id_buf = mm_malloc(sizeof(char)*10);
  ps->letter_id_buf[0] = '\0';
  ps->letter_id_buf_len = 10;
  //set up expected queue
  ps->expected_stack = linklst_create();
  return ps;
}

/*****************************************************************************
 * Destroys the data used by the SAX parser
 ****************************************************************************/
void destroy_streme_io_xml_sax_context(void *ctx) {
  PS_T *ps = (PS_T*)ctx;
  if (ps->motif_id) free(ps->motif_id);
  if (ps->freqs) free(ps->freqs);
  rbtree_destroy(ps->alph_ids);
  free(ps->characters.buffer);
  attrbuf_free(&(ps->attrbuf));
  linklst_destroy_all(ps->expected_stack, free);
  free(ps);
}
