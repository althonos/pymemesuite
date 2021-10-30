#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "alphabet.h"
#include "linked-list.h"
#include "meme-sax.h"
#include "mtype.h"
#include "sax-parser-utils.h"
#include "utils.h"

/*****************************************************************************
 * Used to keep track of what elements are expected and enforces limits on the
 * number of times a valid element can be repeated. Also avoids the 
 * possibility that elements could be left out.
 ****************************************************************************/
#define ES_ZERO_OR_ONE 1
#define ES_ONCE 2
#define ES_ONE_OR_MORE 3
#define ES_ANY 4

/*****************************************************************************
 * Parser State Constants 
 ****************************************************************************/
enum parser_state {
  PS_ERROR,
  PS_START,
  PS_IN_MEME,
  PS_IN_TRAINING_SET,
  PS_IN_MODEL,
  PS_IN_MOTIFS,
  PS_IN_SCANNED_SITES_SUMMARY,
  PS_IN_ALPHABET,
  PS_IN_AMBIGS,
  PS_IN_SEQUENCE,
  PS_IN_LETTER_FREQUENCIES,
  PS_IN_COMMAND_LINE,
  PS_IN_HOST,
  PS_IN_TYPE,
  PS_IN_NMOTIFS,
  PS_IN_EVALUE_THRESHOLD,
  PS_IN_OBJECT_FUNCTION,
  PS_IN_SPFUN,
  PS_IN_USE_LLR,
  PS_IN_MIN_WIDTH,
  PS_IN_MAX_WIDTH,
  PS_IN_MINIC,
  PS_IN_WG,
  PS_IN_WS,
  PS_IN_ENDGAPS,
  PS_IN_SUBSTRING,
  PS_IN_MINSITES,
  PS_IN_MAXSITES,
  PS_IN_WNSITES,
  PS_IN_PROB,
  PS_IN_SPMAP,
  PS_IN_SPFUZZ,
  PS_IN_PRIOR,
  PS_IN_BETA,
  PS_IN_MAXITER,
  PS_IN_DISTANCE,
  PS_IN_NUM_SEQUENCES,
  PS_IN_NUM_POSITIONS,
  PS_IN_SEED,
  PS_IN_SEQFRAC,
  PS_IN_HSFRAC,
  PS_IN_SEARCHSIZE,
  PS_IN_MAXWORDS,
  PS_IN_MAXSIZE,
  PS_IN_NORAND,
  PS_IN_CSITES,
  PS_IN_STRANDS,
  PS_IN_BRIEF,
  PS_IN_PSP_FILE,
  PS_IN_PRIORS_FILE,
  PS_IN_REASON_FOR_STOPPING,
  PS_IN_BACK_ORDER,
  PS_IN_BACKGROUND_FREQUENCIES,
  PS_IN_MOTIF,
  PS_IN_SCANNED_SITES,
  PS_IN_ALPHABET_LETTER,
  PS_IN_AMBIGS_LETTER,
  PS_IN_LF_ALPHABET_ARRAY,
  PS_IN_BF_ALPHABET_ARRAY,
  PS_IN_SCORES,
  PS_IN_PROBABILITIES,
  PS_IN_REGULAR_EXPRESSION,
  PS_IN_CONTRIBUTING_SITES,
  PS_IN_SCANNED_SITE,
  PS_IN_LF_AA_VALUE,
  PS_IN_BF_AA_VALUE,
  PS_IN_SC_ALPHABET_MATRIX,
  PS_IN_PR_ALPHABET_MATRIX,
  PS_IN_CONTRIBUTING_SITE,
  PS_IN_SC_AM_ALPHABET_ARRAY,
  PS_IN_PR_AM_ALPHABET_ARRAY,
  PS_IN_LEFT_FLANK,
  PS_IN_SITE,
  PS_IN_RIGHT_FLANK,
  PS_IN_SC_AM_AA_VALUE,
  PS_IN_PR_AM_AA_VALUE,
  PS_IN_LETTER_REF,
  PS_END
};
typedef enum parser_state PS_EN;

/*****************************************************************************
 * Parser State Names (for error messages)
 ****************************************************************************/
static char const * const state_names[] = {
  "ERROR",
  "START",
  "IN_MEME",
  "IN_TRAINING_SET",
  "IN_MODEL",
  "IN_MOTIFS",
  "IN_SCANNED_SITES_SUMMARY",
  "IN_ALPHABET",
  "IN_AMBIGS",
  "IN_SEQUENCE",
  "IN_LETTER_FREQUENCIES",
  "IN_COMMAND_LINE",
  "IN_HOST",
  "IN_TYPE",
  "IN_NMOTIFS",
  "IN_EVALUE_THRESHOLD",
  "IN_OBJECT_FUNCTION",
  "IN_SPFUN",
  "IN_USE_LLR",
  "IN_MIN_WIDTH",
  "IN_MAX_WIDTH",
  "IN_MINIC",
  "IN_WG",
  "IN_WS",
  "IN_ENDGAPS",
  "IN_SUBSTRING",
  "IN_MINSITES",
  "IN_MAXSITES",
  "IN_WNSITES",
  "IN_PROB",
  "IN_SPMAP",
  "IN_SPFUZZ",
  "IN_PRIOR",
  "IN_BETA",
  "IN_MAXITER",
  "IN_DISTANCE",
  "IN_NUM_SEQUENCES",
  "IN_NUM_POSITIONS",
  "IN_SEED",
  "IN_SEQFRAC",
  "IN_HSFRAC",
  "IN_SEARCHSIZE",
  "IN_MAXWORDS",
  "IN_MAXSIZE",
  "IN_NORAND",
  "IN_CSITES",
  "IN_STRANDS",
  "IN_BRIEF", 
  "IN_PSP_FILE",
  "IN_PRIORS_FILE",
  "IN_REASON_FOR_STOPPING",
  "IN_BACK_ORDER",
  "IN_BACKGROUND_FREQUENCIES",
  "IN_MOTIF",
  "IN_SCANNED_SITES",
  "IN_ALPHABET_LETTER",
  "IN_AMBIGS_LETTER",
  "IN_LF_ALPHABET_ARRAY",
  "IN_BF_ALPHABET_ARRAY",
  "IN_SCORES",
  "IN_PROBABILITIES",
  "IN_REGULAR_EXPRESSION",
  "IN_CONTRIBUTING_SITES",
  "IN_SCANNED_SITE",
  "IN_LF_AA_VALUE",
  "IN_BF_AA_VALUE",
  "IN_SC_ALPHABET_MATRIX",
  "IN_PR_ALPHABET_MATRIX",
  "IN_CONTRIBUTING_SITE",
  "IN_SC_AM_ALPHABET_ARRAY",
  "IN_PR_AM_ALPHABET_ARRAY",
  "IN_LEFT_FLANK",
  "IN_SITE",
  "IN_RIGHT_FLANK",
  "IN_SC_AM_AA_VALUE",
  "IN_PR_AM_AA_VALUE",
  "IN_LETTER_REF",
  "END"
};

/*****************************************************************************
 * Stores the expected state
 ****************************************************************************/
typedef struct es ES_T;
struct es {
  PS_EN state;    // the element that is expected
  int expected; // the number of times the element is expected (see ES codes)
  int found;    // the number of times the element has been found already
};

/*****************************************************************************
 * Stores the parser state.
 ****************************************************************************/
typedef struct ps PS_T;
struct ps {
  PS_EN state;                // the state code (see PS codes)
  int udepth;               // the depth of unknown tags (disabled)
  CHARBUF_T characters;     // a buffer to store characters between elements
  ATTRBUF_T attrbuf;
  char *letter_id_buf;      // a buffer to store the current letter id 
  int letter_id_buf_len;    // the length of the letter id buffer
  int reported_alen;        // reported length of the alphabet, to confirm
  int seen_alen;            // number of alphabet entries seen
  LINKLST_T *expected_stack;// a queue of expected states
  MEME_IO_XML_CALLBACKS_T *callbacks; // callbacks to the user
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
void meme_attr_parse_error(void *state, int errcode, const char *tag, const char *attr, const char *value) {
  PS_T *ps = (PS_T*)state;
  if (errcode == PARSE_ATTR_DUPLICATE) {
    error(ps, "MEME IO XML parser error: Duplicate attribute %s::%s.\n", tag, attr);
  } else if (errcode == PARSE_ATTR_BAD_VALUE) {
    error(ps, "MEME IO XML parser error: Bad value \"%s\" for attribute %s::%s.\n", value, tag, attr);
  } else if (errcode == PARSE_ATTR_MISSING) {
    error(ps, "MEME IO XML parser error: Missing required attribute %s::%s.\n", tag, attr);
  }
}

/*****************************************************************************
 * Push an expected state on the stack
 ****************************************************************************/
void meme_push_es(PS_T *ps, PS_EN expected_state, int expected_occurances) {
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
int meme_update_es(PS_T *ps, PS_EN next_state) {
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
        if (old_es.found < 1) {
          error(ps, "Expected state %s not found!\n", state_names[old_es.state]);
          return false;
        }
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
 * MEME
 ****************************************************************************/
static void start_ele_meme(PS_T *ps, const xmlChar **attrs) {
  struct prog_version ver;
  char *release;

  char* names[2] = {"release", "version"};
  int (*parsers[2])(char*, void*) = {ld_str, ld_version};
  void *data[2] = {&release, &ver};
  bool required[2] = {true, true};
  bool done[2];

  parse_attributes(meme_attr_parse_error, ps, "meme", attrs, 2, names, parsers, data, required, done);
  
  if (ps->callbacks->start_meme && ps->state != PS_ERROR) {
    ps->callbacks->start_meme(ps->user_data, ver.major, ver.minor, ver.patch, release);
  }
  meme_push_es(ps, PS_IN_SCANNED_SITES_SUMMARY, ES_ZERO_OR_ONE);
  meme_push_es(ps, PS_IN_MOTIFS, ES_ONCE);
  meme_push_es(ps, PS_IN_MODEL, ES_ONCE);
  meme_push_es(ps, PS_IN_TRAINING_SET, ES_ONCE);
}

/*****************************************************************************
 * /MEME
 ****************************************************************************/
static void end_ele_meme(PS_T *ps) {
  if (ps->callbacks->end_meme) {
    ps->callbacks->end_meme(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > training_set
 ****************************************************************************/
static void start_ele_training_set(PS_T *ps, const xmlChar **attrs) {
  char *primary_sequences;
  int primary_count;
  int primary_positions;
  char *control_sequences;
  int control_count;

  //
  // Handle both the old and new MEME formats.
  //
  // Attributes must be listed alphabetically.
  char* names[7] = {
    "control_count",		// new MEME format 4-Aug-2017
    "control_sequences",	// new MEME format 4-Aug-2017
    "datafile",			// old name for primary_sequences
    "length",			// old name for primary_count
    "primary_count",
    "primary_positions",
    "primary_sequences"
  };
  int (*parsers[7])(char*, void*) = {ld_int, ld_str, ld_str, ld_int, ld_int, ld_int, ld_str};
  void *data[7] = {&control_count, &control_sequences, &primary_sequences, &primary_count, &primary_count, &primary_positions, &primary_sequences};
  bool required[7] = {false, false, false, false, false, false, false};
  bool done[7];

  parse_attributes(meme_attr_parse_error, ps, "training_set", attrs, 7, names, parsers, data, required, done);

  // Handle the missing fields in the old MEME format.
  if (! done[0]) control_count = 0;
  if (! done[1]) control_sequences = "--none--";
  if (! done[5]) primary_positions = 0;

  // Handle the different names for primary_sequences and primary_count and primary_positions in the old format.
  //   primary_sequences missing?
  if (! (done[2] || done[6])) meme_attr_parse_error(ps, PARSE_ATTR_MISSING, "training_set", names[6], NULL);
  //   primary_count missing?
  if (! (done[3] || done[4])) meme_attr_parse_error(ps, PARSE_ATTR_MISSING, "training_set", names[4], NULL);

  if (ps->callbacks->start_training_set && ps->state != PS_ERROR) {
    ps->callbacks->start_training_set(ps->user_data, primary_sequences, primary_count, primary_positions);
  }
  meme_push_es(ps, PS_IN_LETTER_FREQUENCIES, ES_ONCE);
  meme_push_es(ps, PS_IN_SEQUENCE, ES_ANY);		// Allow to be missing if -brief
  meme_push_es(ps, PS_IN_AMBIGS, ES_ZERO_OR_ONE);
  meme_push_es(ps, PS_IN_ALPHABET, ES_ONCE);
}

/*****************************************************************************
 * MEME > /training_set
 ****************************************************************************/
static void end_ele_training_set(PS_T *ps) {
  if (ps->callbacks->end_training_set) {
    ps->callbacks->end_training_set(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > model
 ****************************************************************************/
static void start_ele_model(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_model && ps->state != PS_ERROR) {
    ps->callbacks->start_model(ps->user_data);
  }
  meme_push_es(ps, PS_IN_BACKGROUND_FREQUENCIES, ES_ONCE);
  meme_push_es(ps, PS_IN_BACK_ORDER, ES_ZERO_OR_ONE);	// not present in new format 4-Aug-2017
  meme_push_es(ps, PS_IN_REASON_FOR_STOPPING, ES_ONCE);
  meme_push_es(ps, PS_IN_PRIORS_FILE, ES_ONCE);
  meme_push_es(ps, PS_IN_PSP_FILE, ES_ZERO_OR_ONE);	// not present in old format 4-Aug-2017
  meme_push_es(ps, PS_IN_BRIEF, ES_ZERO_OR_ONE);	// not present in old format 4-Aug-2017
  meme_push_es(ps, PS_IN_STRANDS, ES_ONCE);
  meme_push_es(ps, PS_IN_MAXWORDS, ES_ZERO_OR_ONE);	// removed in Version 5.0
  meme_push_es(ps, PS_IN_CSITES, ES_ZERO_OR_ONE);
  meme_push_es(ps, PS_IN_NORAND, ES_ZERO_OR_ONE);
  meme_push_es(ps, PS_IN_MAXSIZE, ES_ZERO_OR_ONE);
  meme_push_es(ps, PS_IN_SEARCHSIZE, ES_ZERO_OR_ONE);
  meme_push_es(ps, PS_IN_HSFRAC, ES_ZERO_OR_ONE);
  meme_push_es(ps, PS_IN_SEQFRAC, ES_ZERO_OR_ONE);	// not present in new format 4-Aug-2017
  meme_push_es(ps, PS_IN_SEED, ES_ONCE);
  meme_push_es(ps, PS_IN_NUM_POSITIONS, ES_ONCE);
  meme_push_es(ps, PS_IN_NUM_SEQUENCES, ES_ZERO_OR_ONE); // not present in new format 4-Aug-2017
  meme_push_es(ps, PS_IN_DISTANCE, ES_ONCE);
  meme_push_es(ps, PS_IN_MAXITER, ES_ONCE);
  meme_push_es(ps, PS_IN_BETA, ES_ONCE);
  meme_push_es(ps, PS_IN_PRIOR, ES_ONCE);
  meme_push_es(ps, PS_IN_SPFUZZ, ES_ONCE);
  meme_push_es(ps, PS_IN_SPMAP, ES_ONCE);
  meme_push_es(ps, PS_IN_PROB, ES_ZERO_OR_ONE);		// not present in new format
  meme_push_es(ps, PS_IN_WNSITES, ES_ONCE);
  meme_push_es(ps, PS_IN_MAXSITES, ES_ONCE);
  meme_push_es(ps, PS_IN_MINSITES, ES_ONCE);
  meme_push_es(ps, PS_IN_SUBSTRING, ES_ZERO_OR_ONE);
  meme_push_es(ps, PS_IN_ENDGAPS, ES_ZERO_OR_ONE);	// optional now 23-Oct-2017
  meme_push_es(ps, PS_IN_WS, ES_ZERO_OR_ONE);		// optional now 23-Oct-2017
  meme_push_es(ps, PS_IN_WG, ES_ZERO_OR_ONE);		// optional now 23-Oct-2017
  meme_push_es(ps, PS_IN_MINIC, ES_ZERO_OR_ONE);	// not present in new format 4-Aug-2017
  meme_push_es(ps, PS_IN_MAX_WIDTH, ES_ONCE);
  meme_push_es(ps, PS_IN_MIN_WIDTH, ES_ONCE);
  meme_push_es(ps, PS_IN_USE_LLR, ES_ZERO_OR_ONE);	// not present in new format 4-Aug-2017
  meme_push_es(ps, PS_IN_SPFUN, ES_ZERO_OR_ONE);	// not present in old format 4-Aug-2017 
  meme_push_es(ps, PS_IN_OBJECT_FUNCTION, ES_ONCE);
  meme_push_es(ps, PS_IN_EVALUE_THRESHOLD, ES_ONCE);
  meme_push_es(ps, PS_IN_NMOTIFS, ES_ONCE);
  meme_push_es(ps, PS_IN_TYPE, ES_ONCE);
  meme_push_es(ps, PS_IN_HOST, ES_ONCE);
  meme_push_es(ps, PS_IN_COMMAND_LINE, ES_ONCE);
}

/*****************************************************************************
 * MEME > /model
 ****************************************************************************/
static void end_ele_model(PS_T *ps) {
  if (ps->callbacks->end_model && ps->state != PS_ERROR) {
    ps->callbacks->end_model(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > motifs
 ****************************************************************************/
static void start_ele_motifs(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_motifs && ps->state != PS_ERROR) {
    ps->callbacks->start_motifs(ps->user_data);
  }
  meme_push_es(ps, PS_IN_MOTIF, ES_ANY);
}

/*****************************************************************************
 * MEME > /motifs
 ****************************************************************************/
static void end_ele_motifs(PS_T *ps) {
  if (ps->callbacks->end_motifs && ps->state != PS_ERROR) {
    ps->callbacks->end_motifs(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > scanned_sites_summary
 ****************************************************************************/
static void start_ele_scanned_sites_summary(PS_T *ps, const xmlChar **attrs) {
  double p_thresh;

  char* names[1] = {"p_thresh"};
  int (*parsers[1])(char*, void*) = {ld_double};
  void *data[1] = {&p_thresh};
  bool required[1] = {true};
  bool done[1];

  parse_attributes(meme_attr_parse_error, ps, "scanned_sites_summary", attrs, 1, names, parsers, data, required, done);

  if (ps->callbacks->start_scanned_sites_summary && ps->state != PS_ERROR) {
    ps->callbacks->start_scanned_sites_summary(ps->user_data, p_thresh);
  }
  meme_push_es(ps, PS_IN_SCANNED_SITES, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * MEME > /scanned_sites_summary
 ****************************************************************************/
static void end_ele_scanned_sites_summary(PS_T *ps) {
  if (ps->callbacks->end_scanned_sites_summary && ps->state != PS_ERROR) {
    ps->callbacks->end_scanned_sites_summary(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > training_set > alphabet
 ****************************************************************************/
static void start_ele_alphabet(PS_T *ps, const xmlChar **attrs) {
  int length, type, extends;
  char *name;

  char* extends_options[3] = {"dna", "protein", "rna"};
  int extends_values[3] = {ALPH_FLAG_EXTENDS_DNA, ALPH_FLAG_EXTENDS_PROTEIN, ALPH_FLAG_EXTENDS_RNA};
  MULTI_T extends_multi = {.count = 3, .options = extends_options, .outputs = extends_values, .target = &(extends)};

  char* type_options[2] = {"amino-acid", "nucleotide"};
  int type_values[2] = {MEME_IO_XML_ALPH_PROTEIN, MEME_IO_XML_ALPH_DNA};
  MULTI_T type_multi = {.count = 2, .options = type_options, .outputs = type_values, .target = &(type)};

  char* names[4] = {"id", "length", "like", "name"};
  int (*parsers[4])(char*, void*) = {ld_multi, ld_int, ld_multi, ld_str};
  void *data[4] = {&type_multi, &length, &extends_multi, &name};
  bool required[4] = {false, false, false, false};
  bool done[4];

  // defaults
  type = MEME_IO_XML_ALPH_CUSTOM;
  extends = 0;
  name = NULL;
  parse_attributes(meme_attr_parse_error, ps, "alphabet", attrs, 4, names, parsers, data, required, done);
  ps->reported_alen = (done[1] ? length : -1);

  if (ps->state != PS_ERROR) {
    if (done[0] && done[1]) {
      switch ((MEME_IO_XML_ALPH_T)type) {
        case MEME_IO_XML_ALPH_DNA:
          if (length != 4) 
            error(ps, "DNA or RNA has an alphabet length of 4 but "
                "the file reported a length of %d.\n", length);
          break;
        case MEME_IO_XML_ALPH_PROTEIN:
          if (length != 20)
            error(ps, "Amino acid has an alphabet length of 20 but "
                "the file reported a  length of %d.\n", length);
          break;
        case MEME_IO_XML_ALPH_CUSTOM:
          if (length < 2)
            error(ps, "The alphabet length of %d does not make sense.\n", length);
          break; 
      }
    }
  }

  if (ps->callbacks->start_alphabet && ps->state != PS_ERROR) {
    ps->callbacks->start_alphabet(ps->user_data, (MEME_IO_XML_ALPH_T)type, name, extends);
  }
  meme_push_es(ps, PS_IN_ALPHABET_LETTER, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * MEME > training_set > /alphabet
 ****************************************************************************/
static void end_ele_alphabet(PS_T *ps) {
  if (ps->reported_alen != -1 && ps->reported_alen != ps->seen_alen) {
    error(ps, "Alphabet did not contain the number of letters reported. "
        "Reported %d but saw %d.\n", ps->reported_alen, ps->seen_alen);
  }
  if (ps->callbacks->end_alphabet && ps->state != PS_ERROR) {
    ps->callbacks->end_alphabet(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > training_set > ambigs
 ****************************************************************************/
static void start_ele_ambigs(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_ambigs && ps->state != PS_ERROR) {
    ps->callbacks->start_ambigs(ps->user_data);
  }
  meme_push_es(ps, PS_IN_AMBIGS_LETTER, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * MEME > training_set > /ambigs
 ****************************************************************************/
static void end_ele_ambigs(PS_T *ps) {
  if (ps->callbacks->end_ambigs && ps->state != PS_ERROR) {
    ps->callbacks->end_ambigs(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > training_set > sequence
 ****************************************************************************/
static void start_ele_sequence(PS_T *ps, const xmlChar **attrs) {
  char *id, *name;
  int length;
  double weight;

  char* names[4] = {"id", "length", "name", "weight"};
  int (*parsers[4])(char*, void*) = {ld_str, ld_int, ld_str, ld_double};
  void *data[4] = {&id, &length, &name, &weight};
  bool required[4] = {true, true, true, true};
  bool done[4];

  parse_attributes(meme_attr_parse_error, ps, "sequence", attrs, 4, names, parsers, data, required, done);

  if (ps->callbacks->handle_sequence && ps->state != PS_ERROR) {
    ps->callbacks->handle_sequence(ps->user_data, id, name, length, weight);
  }
}

/*****************************************************************************
 * MEME > training_set > letter_frequencies
 ****************************************************************************/
static void start_ele_letter_frequencies(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_letter_frequencies && ps->state != PS_ERROR) {
    ps->callbacks->start_letter_frequencies(ps->user_data);
  }
  meme_push_es(ps, PS_IN_LF_ALPHABET_ARRAY, ES_ONCE);
}

/*****************************************************************************
 * MEME > training_set > /letter_frequencies
 ****************************************************************************/
static void end_ele_letter_frequencies(PS_T *ps) {
  if (ps->callbacks->end_letter_frequencies && ps->state != PS_ERROR) {
    ps->callbacks->end_letter_frequencies(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > model > /command_line
 ****************************************************************************/
static void end_ele_command_line(PS_T *ps) {
  if (ps->callbacks->handle_command_line && ps->state != PS_ERROR) {
    ps->callbacks->handle_command_line(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * MEME > model > /host
 ****************************************************************************/
static void end_ele_host(PS_T *ps) {
  if (ps->callbacks->handle_host && ps->state != PS_ERROR) {
    ps->callbacks->handle_host(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * MEME > model > /type
 ****************************************************************************/
static void end_ele_type(PS_T *ps) {
  int type;
  char* endgaps_options[4] = {"anr", "oops", "tcm", "zoops"};
  int endgaps_values[4] = {Tcm, Oops, Tcm, Zoops}; 
  MULTI_T endgaps_multi = {.count = 4, .options = endgaps_options, .outputs = endgaps_values, .target = &(type)};
  if (ld_multi(ps->characters.buffer, &endgaps_multi)) {
    error(ps, "Couldn't parse type from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_type && ps->state != PS_ERROR) {
    ps->callbacks->handle_type(ps->user_data, type); 
  }
}

/*****************************************************************************
 * MEME > model > /nmotifs
 ****************************************************************************/
static void end_ele_nmotifs(PS_T *ps) {
  int nmotifs;
  if (ld_int(ps->characters.buffer, &nmotifs)) {
    error(ps, "Couldn't parse number of motifs from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_nmotifs && ps->state != PS_ERROR) {
    ps->callbacks->handle_nmotifs(ps->user_data, nmotifs);
  }
}

/*****************************************************************************
 * MEME > model > /evalue_threshold
 ****************************************************************************/
static void end_ele_evalue_threshold(PS_T *ps) {
  double log10_evalue_threshold;
  if (ld_log10_ev(ps->characters.buffer, &log10_evalue_threshold)) {
    error(ps, "Couldn't parse evalue threshold from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_evalue_threshold && ps->state != PS_ERROR) {
    ps->callbacks->handle_evalue_threshold(ps->user_data, log10_evalue_threshold);
  }
}

/*****************************************************************************
 * MEME > model > /object_function
 ****************************************************************************/
static void end_ele_object_function(PS_T *ps) {
  if (ps->callbacks->handle_object_function && ps->state != PS_ERROR) {
    ps->callbacks->handle_object_function(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * MEME > model > /use_spfun
 ****************************************************************************/
static void end_ele_spfun(PS_T *ps) {
  if (ps->callbacks->handle_spfun && ps->state != PS_ERROR) {
    ps->callbacks->handle_spfun(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * MEME > model > /use_llr
 ****************************************************************************/
static void end_ele_use_llr(PS_T *ps) {
  int use_llr;
  if (ld_int(ps->characters.buffer, &use_llr)) {
    error(ps, "Couldn't parse use_llr from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_use_llr && ps->state != PS_ERROR) {
    ps->callbacks->handle_use_llr(ps->user_data, use_llr);
  }
}

/*****************************************************************************
 * MEME > model > /min_width
 ****************************************************************************/
static void end_ele_min_width(PS_T *ps) {
  int min_width;
  if (ld_int(ps->characters.buffer, &min_width)) {
    error(ps, "Couldn't parse motif minimum width from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_min_width && ps->state != PS_ERROR) {
    ps->callbacks->handle_min_width(ps->user_data, min_width);
  }
}

/*****************************************************************************
 * MEME > model > /max_width
 ****************************************************************************/
static void end_ele_max_width(PS_T *ps) {
  int max_width;
  if (ld_int(ps->characters.buffer, &max_width)) {
    error(ps, "Couldn't parse motif maximum width from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_max_width && ps->state != PS_ERROR) {
    ps->callbacks->handle_max_width(ps->user_data, max_width);
  }
}

/*****************************************************************************
 * MEME > model > /minic
 ****************************************************************************/
static void end_ele_minic(PS_T *ps) {
  double minic;
  if (ld_double(ps->characters.buffer, &minic)) {
    error(ps, "Couldn't parse motif minimum information content from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_minic && ps->state != PS_ERROR) {
    ps->callbacks->handle_minic(ps->user_data, minic);
  }
}

/*****************************************************************************
 * MEME > model > /wg
 ****************************************************************************/
static void end_ele_wg(PS_T *ps) {
  double wg;
  if (ld_double(ps->characters.buffer, &wg)) {
    error(ps, "Couldn't parse wg from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_wg && ps->state != PS_ERROR) {
    ps->callbacks->handle_wg(ps->user_data, wg);
  }
}

/*****************************************************************************
 * MEME > model > /ws
 ****************************************************************************/
static void end_ele_ws(PS_T *ps) {
  double ws;
  if (ld_double(ps->characters.buffer, &ws)) {
    error(ps, "Couldn't parse ws from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_ws && ps->state != PS_ERROR) {
    ps->callbacks->handle_ws(ps->user_data, ws);
  }
}

/*****************************************************************************
 * MEME > model > /endgaps
 ****************************************************************************/
static void end_ele_endgaps(PS_T *ps) {
  int use_endgaps;
  char* endgaps_options[2] = {"no", "yes"};
  int endgaps_values[2] = {false, true}; 
  MULTI_T endgaps_multi = {.count = 2, .options = endgaps_options, .outputs = endgaps_values, .target = &(use_endgaps)};
  if (ld_multi(ps->characters.buffer, &endgaps_multi)) {
    error(ps, "Couldn't parse endgaps from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_endgaps && ps->state != PS_ERROR) {
    ps->callbacks->handle_endgaps(ps->user_data, use_endgaps);
  }
}

/*****************************************************************************
 * MEME > model > /substring
 ****************************************************************************/
static void end_ele_substring(PS_T *ps) {
  int use_substring;
  char* substring_options[2] = {"no", "yes"};
  int substring_values[2] = {false, true}; 
  MULTI_T substring_multi = {.count = 2, .options = substring_options, .outputs = substring_values, .target = &(use_substring)};
  if (ld_multi(ps->characters.buffer, &substring_multi)) {
    error(ps, "Couldn't parse substring from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_substring && ps->state != PS_ERROR) {
    ps->callbacks->handle_substring(ps->user_data, use_substring);
  }
}

/*****************************************************************************
 * MEME > model > /minsites
 ****************************************************************************/
static void end_ele_minsites(PS_T *ps) {
  int minsites;
  if (ld_int(ps->characters.buffer, &minsites)) {
    error(ps, "Couldn't parse minsites from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_minsites && ps->state != PS_ERROR) {
    ps->callbacks->handle_minsites(ps->user_data, minsites);
  }
}

/*****************************************************************************
 * MEME > model > /maxsites
 ****************************************************************************/
static void end_ele_maxsites(PS_T *ps) {
  int maxsites;
  if (ld_int(ps->characters.buffer, &maxsites)) {
    error(ps, "Couldn't parser maxsites from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_maxsites && ps->state != PS_ERROR) {
    ps->callbacks->handle_maxsites(ps->user_data, maxsites);
  }
}

/*****************************************************************************
 * MEME > model > /wnsites
 ****************************************************************************/
static void end_ele_wnsites(PS_T *ps) {
  double wnsites;
  if (ld_double(ps->characters.buffer, &wnsites)) {
    error(ps, "Couldn't parse wnsites from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_wnsites && ps->state != PS_ERROR) {
    ps->callbacks->handle_wnsites(ps->user_data, wnsites);
  }
}

/*****************************************************************************
 * MEME > model > /prob
 ****************************************************************************/
static void end_ele_prob(PS_T *ps) {
  double prob;
  if (ld_double(ps->characters.buffer, &prob)) {
    error(ps, "Couldn't parse prob from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_prob && ps->state != PS_ERROR) {
    ps->callbacks->handle_prob(ps->user_data, prob);
  }
}

/*****************************************************************************
 * MEME > model > /spmap
 ****************************************************************************/
static void end_ele_spmap(PS_T *ps) {
  int map_method;
  char* map_options[2] = {"pam", "uni"};
  int map_values[2] = {1, 0}; 
  MULTI_T map_multi = {.count = 2, .options = map_options, .outputs = map_values, .target = &(map_method)};
  if (ld_multi(ps->characters.buffer, &map_multi)) {
    error(ps, "Couldn't parse spmap from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_spmap && ps->state != PS_ERROR) {
    ps->callbacks->handle_spmap(ps->user_data, map_method);
  }
}

/*****************************************************************************
 * MEME > model > /spfuzz
 ****************************************************************************/
static void end_ele_spfuzz(PS_T *ps) {
  double spfuzz;
  if (ld_double(ps->characters.buffer, &spfuzz)) {
    error(ps, "Couldn't parse spfuzz from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_spfuzz && ps->state != PS_ERROR) {
    ps->callbacks->handle_spfuzz(ps->user_data, spfuzz);
  }
}

/*****************************************************************************
 * MEME > model > /prior
 ****************************************************************************/
static void end_ele_prior(PS_T *ps) {
  int prior_method;
  char* prior_options[5] = {"addone", "dirichlet", "dmix", "megadmix", "megap"};
  //I believe I've given these the same values as the meme PTYPE enum but I should really expose this somehow
  //Or maybe I should import meme.h but that just seems messy.
  int prior_values[5] = {4, 3, 2, 0, 1}; 
  MULTI_T prior_multi = {.count = 5, .options = prior_options, .outputs = prior_values, .target = &(prior_method)};
  if (ld_multi(ps->characters.buffer, &prior_multi)) {
    error(ps, "Couldn't parse prior from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_prior && ps->state != PS_ERROR) {
    ps->callbacks->handle_prior(ps->user_data, prior_method);
  }
}

/*****************************************************************************
 * MEME > model > /beta
 ****************************************************************************/
static void end_ele_beta(PS_T *ps) {
  double beta;
  if (ld_double(ps->characters.buffer, &beta)) {
    error(ps, "Couldn't parse beta from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_beta && ps->state != PS_ERROR) {
    ps->callbacks->handle_beta(ps->user_data, beta);
  }
}

/*****************************************************************************
 * MEME > model > /maxiter
 ****************************************************************************/
static void end_ele_maxiter(PS_T *ps) {
  int maxiter;
  if (ld_int(ps->characters.buffer, &maxiter)) {
    error(ps, "Couldn't parse maxiter from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_maxiter && ps->state != PS_ERROR) {
    ps->callbacks->handle_maxiter(ps->user_data, maxiter);
  }
}

/*****************************************************************************
 * MEME > model > /distance
 ****************************************************************************/
static void end_ele_distance(PS_T *ps) {
  double distance;
  if (ld_double(ps->characters.buffer, &distance)) {
    error(ps, "Couldn't parse distance from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_distance && ps->state != PS_ERROR) {
    ps->callbacks->handle_distance(ps->user_data, distance);
  }
}

/*****************************************************************************
 * MEME > model > /num_sequences
 ****************************************************************************/
static void end_ele_num_sequences(PS_T *ps) {
  int num_sequences;
  if (ld_int(ps->characters.buffer, &num_sequences)) {
    error(ps, "Couldn't parse num_sequences from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_num_sequences && ps->state != PS_ERROR) {
    ps->callbacks->handle_num_sequences(ps->user_data, num_sequences);
  }
}

/*****************************************************************************
 * MEME > model > /num_positions
 ****************************************************************************/
static void end_ele_num_positions(PS_T *ps) {
  int num_positions;
  if (ld_int(ps->characters.buffer, &num_positions)) {
    error(ps, "Couldn't parse num_positions from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_num_positions && ps->state != PS_ERROR) {
    ps->callbacks->handle_num_positions(ps->user_data, num_positions);
  }
}

/*****************************************************************************
 * MEME > model > /seed
 ****************************************************************************/
static void end_ele_seed(PS_T *ps) {
  long seed;
  if (ld_long(ps->characters.buffer, &seed)) {
    error(ps, "Couldn't parse seed from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_seed && ps->state != PS_ERROR) {
    ps->callbacks->handle_seed(ps->user_data, seed);
  }
}

/*****************************************************************************
 * MEME > model > /seqfrac
 ****************************************************************************/
static void end_ele_seqfrac(PS_T *ps) {
  double seqfrac;
  if (ld_double(ps->characters.buffer, &seqfrac)) {
    error(ps, "Couldn't parse seed from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_seqfrac && ps->state != PS_ERROR) {
    ps->callbacks->handle_seqfrac(ps->user_data, seqfrac);
  }
}

/*****************************************************************************
 * MEME > model > /hsfrac
 ****************************************************************************/
static void end_ele_hsfrac(PS_T *ps) {
  double hsfrac;
  if (ld_double(ps->characters.buffer, &hsfrac)) {
    error(ps, "Couldn't parse hsfrac from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_hsfrac && ps->state != PS_ERROR) {
    ps->callbacks->handle_hsfrac(ps->user_data, hsfrac);
  }
}

/*****************************************************************************
 * MEME > model > /norand
 ****************************************************************************/
static void end_ele_norand(PS_T *ps) {
  if (ps->callbacks->handle_norand && ps->state != PS_ERROR) {
    ps->callbacks->handle_norand(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * MEME > model > /maxwords
 * Ignored-- Allowed for backwards compatibility.
 ****************************************************************************/
static void end_ele_maxwords(PS_T *ps) {
  double maxwords;
  if (ld_double(ps->characters.buffer, &maxwords)) {
    error(ps, "Couldn't parse maxwords from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_maxwords && ps->state != PS_ERROR) {
    ps->callbacks->handle_maxwords(ps->user_data, maxwords);
  }
}

/*****************************************************************************
 * MEME > model > /csites
 ****************************************************************************/
static void end_ele_csites(PS_T *ps) {
  double csites;
  if (ld_double(ps->characters.buffer, &csites)) {
    error(ps, "Couldn't parse csites from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_csites&& ps->state != PS_ERROR) {
    ps->callbacks->handle_csites(ps->user_data, csites);
  }
}

/*****************************************************************************
 * MEME > model > /searchsize
 ****************************************************************************/
static void end_ele_searchsize(PS_T *ps) {
  double searchsize;
  if (ld_double(ps->characters.buffer, &searchsize)) {
    error(ps, "Couldn't parse searchsize from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_searchsize && ps->state != PS_ERROR) {
    ps->callbacks->handle_searchsize(ps->user_data, searchsize);
  }
}

/*****************************************************************************
 * MEME > model > /maxsize
 ****************************************************************************/
static void end_ele_maxsize(PS_T *ps) {
  double maxsize;
  if (ld_double(ps->characters.buffer, &maxsize)) {
    error(ps, "Couldn't parse maxsize from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_maxsize && ps->state != PS_ERROR) {
    ps->callbacks->handle_maxsize(ps->user_data, maxsize);
  }
}

/*****************************************************************************
 * MEME > model > /strands
 ****************************************************************************/
static void end_ele_strands(PS_T *ps) {
  int strands;
  char* strands_options[3] = {"both", "forward", "none"};
  int strands_values[3] = {2, 1, 0}; 
  MULTI_T strands_multi = {.count = 3, .options = strands_options, .outputs = strands_values, .target = &(strands)};
  if (ld_multi(ps->characters.buffer, &strands_multi)) {
    error(ps, "Couldn't parse strands from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_strands && ps->state != PS_ERROR) {
    ps->callbacks->handle_strands(ps->user_data, strands);
  }
}

/*****************************************************************************
 * MEME > model > /brief
 ****************************************************************************/
static void end_ele_brief(PS_T *ps) {
  int brief;
  if (ld_int(ps->characters.buffer, &brief)) {
    error(ps, "Couldn't parse brief from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_brief && ps->state != PS_ERROR) {
    ps->callbacks->handle_brief(ps->user_data, brief);
  }
}

/*****************************************************************************
 * MEME > model > /psp_file
 ****************************************************************************/
static void end_ele_psp_file(PS_T *ps) {
  if (ps->callbacks->handle_psp_file && ps->state != PS_ERROR) {
    ps->callbacks->handle_psp_file(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * MEME > model > /priors_file
 ****************************************************************************/
static void end_ele_priors_file(PS_T *ps) {
  if (ps->callbacks->handle_priors_file && ps->state != PS_ERROR) {
    ps->callbacks->handle_priors_file(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * MEME > model > /reason_for_stopping
 ****************************************************************************/
static void end_ele_reason_for_stopping(PS_T *ps) {
  if (ps->callbacks->handle_reason_for_stopping && ps->state != PS_ERROR) {
    ps->callbacks->handle_reason_for_stopping(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * MEME > model > /back_order
 ****************************************************************************/
static void end_ele_back_order(PS_T *ps) {
  int back_order;
  if (ld_int(ps->characters.buffer, &back_order)) {
    error(ps, "Couldn't parse back_order from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_back_order && ps->state != PS_ERROR) {
    ps->callbacks->handle_back_order(ps->user_data, back_order);
  }
}

/*****************************************************************************
 * MEME > model > background_frequencies
 ****************************************************************************/
static void start_ele_background_frequencies(PS_T *ps, const xmlChar **attrs) {
  char *bgsource;

  char* names[1] = {"source"};
  int (*parsers[1])(char*, void*) = {ld_str};
  void *data[1] = {&bgsource};
  bool required[1] = {true};
  bool done[1];

  parse_attributes(meme_attr_parse_error, ps, "background_frequencies", attrs, 1, names, parsers, data, required, done);

  if (ps->callbacks->start_background_frequencies && ps->state != PS_ERROR) {
    ps->callbacks->start_background_frequencies(ps->user_data, bgsource);
  }
  meme_push_es(ps, PS_IN_BF_ALPHABET_ARRAY, ES_ONCE);
}

/*****************************************************************************
 * MEME > model > background_frequencies
 ****************************************************************************/
static void end_ele_background_frequencies(PS_T *ps) {
  if (ps->callbacks->end_background_frequencies && ps->state != PS_ERROR) {
    ps->callbacks->end_background_frequencies(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > motifs > motif
 ****************************************************************************/
static void start_ele_motif(PS_T *ps, const xmlChar **attrs) {
  char *id, *name, *url, *alt;
  int width, sites;
  double ic, re, llr, log10evalue, bayes_threshold, elapsed_time;

  // !!! Note: The name list must be alphabetically ordered!!!
  char* names[12] = {"alt", "bayes_threshold", "e_value", "elapsed_time", 
    "ic", "id", "llr", "name", "re", "sites", "url", "width"};
  int (*parsers[12])(char*, void*) = {ld_str, ld_double, ld_log10_ev, ld_double, 
    ld_double, ld_str, ld_double, ld_str, ld_double, ld_int, ld_str, ld_int};
  void *data[12] = {&alt, &bayes_threshold, &log10evalue, &elapsed_time, 
    &ic, &id, &llr, &name, &re, &sites, &url, &width};
  bool required[12] = {false, true, true, true, 
    true, true, true, true, true, true, false, true};
  bool done[12];

  // set url & alt defaults
  url = "";
  alt = "";

  parse_attributes(meme_attr_parse_error, ps, "motif", attrs, 12, names, 
      parsers, data, required, done);

  if (ps->callbacks->start_motif && ps->state != PS_ERROR) {
    ps->callbacks->start_motif(ps->user_data, id, name, alt, width, sites, llr, ic, 
        re, bayes_threshold, log10evalue, elapsed_time, url);
  }
  meme_push_es(ps, PS_IN_CONTRIBUTING_SITES, ES_ONCE);
  meme_push_es(ps, PS_IN_REGULAR_EXPRESSION, ES_ZERO_OR_ONE);
  meme_push_es(ps, PS_IN_PROBABILITIES, ES_ONCE);
  meme_push_es(ps, PS_IN_SCORES, ES_ONCE);
}

/*****************************************************************************
 * MEME > motifs > /motif
 ****************************************************************************/
static void end_ele_motif(PS_T *ps) {
  if (ps->callbacks->end_motif && ps->state != PS_ERROR) {
    ps->callbacks->end_motif(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > scanned_sites_summary > scanned_sites
 ****************************************************************************/
static void start_ele_scanned_sites(PS_T *ps, const xmlChar **attrs) {
  char *sequence_id;
  double log10pvalue;
  int num_sites;

  char* names[3] = {"num_sites", "pvalue", "sequence_id"};
  int (*parsers[3])(char*, void*) = {ld_int, ld_log10_pv, ld_str};
  void *data[3] = {&num_sites, &log10pvalue, &sequence_id};
  bool required[3] = {true, true, true};
  bool done[3];

  parse_attributes(meme_attr_parse_error, ps, "scanned_sites", attrs, 3, names, parsers, data, required, done);

  if (ps->callbacks->start_scanned_sites && ps->state != PS_ERROR) {
    ps->callbacks->start_scanned_sites(ps->user_data, sequence_id, log10pvalue, num_sites);
  }
  meme_push_es(ps, PS_IN_SCANNED_SITE, ES_ANY);
}

/*****************************************************************************
 * MEME > scanned_sites_summary > /scanned_sites
 ****************************************************************************/
static void end_ele_scanned_sites(PS_T *ps) {
  if (ps->callbacks->end_scanned_sites && ps->state != PS_ERROR) {
    ps->callbacks->end_scanned_sites(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > training_set > alphabet > letter
 ****************************************************************************/
static void start_ele_alphabet_letter(PS_T *ps, const xmlChar **attrs) {
  char *aliases, *id, *name, *equals, symbol, complement;
  int colour;

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
  parse_attributes(meme_attr_parse_error, ps, "letter", attrs, 7, names, parsers, data, required, done);

  if (equals == NULL) ps->seen_alen++;

  if (ps->callbacks->handle_alphabet_letter && ps->state != PS_ERROR) {
    ps->callbacks->handle_alphabet_letter(ps->user_data, id, symbol, aliases, complement, equals, name, colour);
  }
}

/*****************************************************************************
 * MEME > training_set > ambigs > letter
 ****************************************************************************/
static void start_ele_ambigs_letter(PS_T *ps, const xmlChar **attrs) {
  char *id, symbol;

  char* names[2] = {"id", "symbol"};
  int (*parsers[2])(char*, void*) = {ld_str, ld_char};
  void *data[2] = {&id, &symbol};
  bool required[2] = {true, true};
  bool done[2];

  parse_attributes(meme_attr_parse_error, ps, "letter", attrs, 2, names, parsers, data, required, done);

  if (ps->callbacks->handle_ambigs_letter && ps->state != PS_ERROR) {
    ps->callbacks->handle_ambigs_letter(ps->user_data, id, symbol);
  }
}

/*****************************************************************************
 * MEME > training_set > letter_frequencies > alphabet_array
 ****************************************************************************/
static void start_ele_lf_alphabet_array(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_lf_alphabet_array && ps->state != PS_ERROR) {
    ps->callbacks->start_lf_alphabet_array(ps->user_data);
  }
  meme_push_es(ps, PS_IN_LF_AA_VALUE, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * MEME > training_set > letter_frequencies > /alphabet_array
 ****************************************************************************/
static void end_ele_lf_alphabet_array(PS_T *ps) {
  if (ps->callbacks->end_lf_alphabet_array && ps->state != PS_ERROR) {
    ps->callbacks->end_lf_alphabet_array(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > model > background_frequencies > alphabet_array
 ****************************************************************************/
static void start_ele_bf_alphabet_array(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_bf_alphabet_array && ps->state != PS_ERROR) {
    ps->callbacks->start_bf_alphabet_array(ps->user_data);
  }
  meme_push_es(ps, PS_IN_BF_AA_VALUE, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * MEME > model > background_frequencies > /alphabet_array
 ****************************************************************************/
static void end_ele_bf_alphabet_array(PS_T *ps) {
  if (ps->callbacks->end_bf_alphabet_array && ps->state != PS_ERROR) {
    ps->callbacks->end_bf_alphabet_array(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > scores
 ****************************************************************************/
static void start_ele_scores(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_scores && ps->state != PS_ERROR) {
    ps->callbacks->start_scores(ps->user_data);
  }
  meme_push_es(ps, PS_IN_SC_ALPHABET_MATRIX, ES_ONCE);
}

/*****************************************************************************
 * MEME > motifs > motif > /scores
 ****************************************************************************/
static void end_ele_scores(PS_T *ps) {
  if (ps->callbacks->end_scores && ps->state != PS_ERROR) {
    ps->callbacks->end_scores(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > probabilities
 ****************************************************************************/
static void start_ele_probabilities(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_probabilities && ps->state != PS_ERROR) {
    ps->callbacks->start_probabilities(ps->user_data);
  }
  meme_push_es(ps, PS_IN_PR_ALPHABET_MATRIX, ES_ONCE);
}

/*****************************************************************************
 * MEME > motifs > motif > /probabilities
 ****************************************************************************/
static void end_ele_probabilities(PS_T *ps) {
  if (ps->callbacks->end_probabilities && ps->state != PS_ERROR) {
    ps->callbacks->end_probabilities(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > regular_expression
 ****************************************************************************/
static void end_ele_regular_expression(PS_T *ps) {
  if (ps->callbacks->handle_regular_expression && ps->state != PS_ERROR) {
    ps->callbacks->handle_regular_expression(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > contributing_sites
 ****************************************************************************/
static void start_ele_contributing_sites(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_contributing_sites && ps->state != PS_ERROR) {
    ps->callbacks->start_contributing_sites(ps->user_data);
  }
  meme_push_es(ps, PS_IN_CONTRIBUTING_SITE, ES_ANY);	// allow to be missing if -brief
}

/*****************************************************************************
 * MEME > motifs > motif > /contributing_sites
 ****************************************************************************/
static void end_ele_contributing_sites(PS_T *ps) {
  if (ps->callbacks->end_contributing_sites && ps->state != PS_ERROR) {
    ps->callbacks->end_contributing_sites(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > scanned_sites_summary > scanned_sites > scanned_site
 ****************************************************************************/
static void start_ele_scanned_site(PS_T *ps, const xmlChar **attrs) {
  char *motif_id;
  int strand, position;
  double log10pvalue;

  strand = 0; //default to none
  char* strand_options[3] = {"minus", "none", "plus"};
  int strand_values[3] = {'-', 0, '+'}; 
  MULTI_T strand_multi = {.count = 3, .options = strand_options, .outputs = strand_values, .target = &(strand)};
  
  char* names[4] = {"motif_id", "position", "pvalue", "strand"};
  int (*parsers[4])(char*, void*) = {ld_str, ld_int, ld_log10_pv, ld_multi};
  void *data[4] = {&motif_id, &position, &log10pvalue, &strand_multi};
  bool required[4] = {true, true, true, false};
  bool done[4];

  parse_attributes(meme_attr_parse_error, ps, "scanned_site", attrs, 4, names, parsers, data, required, done);

  if (ps->callbacks->handle_scanned_site && ps->state != PS_ERROR) {
    ps->callbacks->handle_scanned_site(ps->user_data, motif_id, (char)strand, position, log10pvalue);
  }
}

/*****************************************************************************
 * MEME > training_set > letter_frequencies > alphabet_array > value
 * MEME > model  >   background_frequencies > alphabet_array > value
 * MEME > motifs > motif >
 * (scores|probabilities) > alphabet_matrix > alphabet_array > value
 ****************************************************************************/
static void start_ele_value(PS_T *ps, const xmlChar **attrs) {
  char *letter_id;
  int letter_id_len;

  char* names[1] = {"letter_id"};
  int (*parsers[1])(char*, void*) = {ld_str};
  void *data[1] = {&letter_id};
  bool required[1] = {true};
  bool done[1];

  parse_attributes(meme_attr_parse_error, ps, "value", attrs, 1, names, parsers, data, required, done);

  letter_id_len = strlen(letter_id);
  if (ps->letter_id_buf_len <= letter_id_len) {
    ps->letter_id_buf = mm_realloc(ps->letter_id_buf, sizeof(char) * (letter_id_len + 1));
    ps->letter_id_buf_len = letter_id_len + 1;
  }
  strncpy(ps->letter_id_buf, letter_id, ps->letter_id_buf_len);
}

/*****************************************************************************
 * MEME > training_set > letter_frequencies > alphabet_array > /value
 ****************************************************************************/
static void end_ele_lf_aa_value(PS_T *ps) {
  double value;
  if (ld_double(ps->characters.buffer, &value)) {
    error(ps, "Couldn't parse value from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_lf_aa_value && ps->state != PS_ERROR) {
    ps->callbacks->handle_lf_aa_value(ps->user_data, ps->letter_id_buf, value);
  }
}

/*****************************************************************************
 * MEME > model > background_frequencies > alphabet_array > /value
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
 * MEME > motifs > motif > scores > alphabet_matrix
 ****************************************************************************/
static void start_ele_sc_alphabet_matrix(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_sc_alphabet_matrix && ps->state != PS_ERROR) {
    ps->callbacks->start_sc_alphabet_matrix(ps->user_data);
  }
  meme_push_es(ps, PS_IN_SC_AM_ALPHABET_ARRAY, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * MEME > motifs > motif > scores > /alphabet_matrix
 ****************************************************************************/
static void end_ele_sc_alphabet_matrix(PS_T *ps) {
  if (ps->callbacks->end_sc_alphabet_matrix && ps->state != PS_ERROR) {
    ps->callbacks->end_sc_alphabet_matrix(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > probabilities > alphabet_matrix
 ****************************************************************************/
static void start_ele_pr_alphabet_matrix(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_pr_alphabet_matrix && ps->state != PS_ERROR) {
    ps->callbacks->start_pr_alphabet_matrix(ps->user_data);
  }
  meme_push_es(ps, PS_IN_PR_AM_ALPHABET_ARRAY, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * MEME > motifs > motif > probabilities > /alphabet_matrix
 ****************************************************************************/
static void end_ele_pr_alphabet_matrix(PS_T *ps) {
  if (ps->callbacks->end_pr_alphabet_matrix && ps->state != PS_ERROR) {
    ps->callbacks->end_pr_alphabet_matrix(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > contributing_sites > contributing_site
 ****************************************************************************/
static void start_ele_contributing_site(PS_T *ps, const xmlChar **attrs) {
  char *sequence_id;
  int position, strand;
  double log10pvalue;


  strand = '\0'; //default to none
  char* strand_options[3] = {"minus", "none", "plus"};
  int strand_values[3] = {'-', '\0', '+'}; 
  MULTI_T strand_multi = {.count = 3, .options = strand_options, .outputs = strand_values, .target = &(strand)};
  
  char* names[4] = {"position", "pvalue", "sequence_id", "strand"};
  int (*parsers[4])(char*, void*) = {ld_int, ld_log10_pv, ld_str, ld_multi};
  void *data[4] = {&position, &log10pvalue, &sequence_id, &strand_multi};
  bool required[4] = {true, true, true, false};
  bool done[4];

  parse_attributes(meme_attr_parse_error, ps, "contributing_site", attrs, 4, names, parsers, data, required, done);

  if (ps->callbacks->start_contributing_site && ps->state != PS_ERROR) {
    ps->callbacks->start_contributing_site(ps->user_data, sequence_id, position, (char)strand, log10pvalue);
  }
  meme_push_es(ps, PS_IN_RIGHT_FLANK, ES_ONCE);
  meme_push_es(ps, PS_IN_SITE, ES_ONCE);
  meme_push_es(ps, PS_IN_LEFT_FLANK, ES_ONCE);
}

/*****************************************************************************
 * MEME > motifs > motif > contributing_sites > /contributing_site
 ****************************************************************************/
static void end_ele_contributing_site(PS_T *ps) {
  if (ps->callbacks->end_contributing_site && ps->state != PS_ERROR) {
    ps->callbacks->end_contributing_site(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > scores > alphabet_matrix > alphabet_array
 ****************************************************************************/
static void start_ele_sc_am_alphabet_array(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_sc_am_alphabet_array && ps->state != PS_ERROR) {
    ps->callbacks->start_sc_am_alphabet_array(ps->user_data);
  }
  meme_push_es(ps, PS_IN_SC_AM_AA_VALUE, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * MEME > motifs > motif > scores > alphabet_matrix > /alphabet_array
 ****************************************************************************/
static void end_ele_sc_am_alphabet_array(PS_T *ps) {
  if (ps->callbacks->end_sc_am_alphabet_array && ps->state != PS_ERROR) {
    ps->callbacks->end_sc_am_alphabet_array(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > probabilities > alphabet_matrix > alphabet_array
 ****************************************************************************/
static void start_ele_pr_am_alphabet_array(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_pr_am_alphabet_array && ps->state != PS_ERROR) {
    ps->callbacks->start_pr_am_alphabet_array(ps->user_data);
  }
  meme_push_es(ps, PS_IN_PR_AM_AA_VALUE, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * MEME > motifs > motif > probabilities > alphabet_matrix > /alphabet_array
 ****************************************************************************/
static void end_ele_pr_am_alphabet_array(PS_T *ps) {
  if (ps->callbacks->end_pr_am_alphabet_array && ps->state != PS_ERROR) {
    ps->callbacks->end_pr_am_alphabet_array(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > contributing_sites > contributing_site > /left_flank
 ****************************************************************************/
static void end_ele_left_flank(PS_T *ps) {
  if (ps->callbacks->handle_left_flank && ps->state != PS_ERROR) {
    ps->callbacks->handle_left_flank(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > contributing_sites > contributing_site > site
 ****************************************************************************/
static void start_ele_site(PS_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_site && ps->state != PS_ERROR) {
    ps->callbacks->start_site(ps->user_data);
  }
  meme_push_es(ps, PS_IN_LETTER_REF, ES_ONE_OR_MORE);
}

/*****************************************************************************
 * MEME > motifs > motif > contributing_sites > contributing_site > /site
 ****************************************************************************/
static void end_ele_site(PS_T *ps) {
  if (ps->callbacks->end_site && ps->state != PS_ERROR) {
    ps->callbacks->end_site(ps->user_data);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > contributing_sites > contributing_site > /right_flank
 ****************************************************************************/
static void end_ele_right_flank(PS_T *ps) {
  if (ps->callbacks->handle_right_flank && ps->state != PS_ERROR) {
    ps->callbacks->handle_right_flank(ps->user_data, ps->characters.buffer);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > scores > alphabet_matrix > alphabet_array > /value
 ****************************************************************************/
static void end_ele_sc_am_aa_value(PS_T *ps) {
  double value;
  if (ld_double(ps->characters.buffer, &value)) {
    error(ps, "Couldn't parse value from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_sc_am_aa_value && ps->state != PS_ERROR) {
    ps->callbacks->handle_sc_am_aa_value(ps->user_data, ps->letter_id_buf, value);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > probabilities > alphabet_matrix > 
 *  alphabet_array > /value
 ****************************************************************************/
static void end_ele_pr_am_aa_value(PS_T *ps) {
  double value;
  if (ld_double(ps->characters.buffer, &value)) {
    error(ps, "Couldn't parse value from \"%s\"\n", ps->characters.buffer);
  }
  if (ps->callbacks->handle_pr_am_aa_value && ps->state != PS_ERROR) {
    ps->callbacks->handle_pr_am_aa_value(ps->user_data, ps->letter_id_buf, value);
  }
}

/*****************************************************************************
 * MEME > motifs > motif > contributing_sites > contributing_site > site >
 *  letter_ref
 * Handle letter_ref starting element (EMPTY)
 ****************************************************************************/
static void start_ele_letter_ref(PS_T *ps, const xmlChar **attrs) {
  char *letter_id;

  char* names[1] = {"letter_id"};
  int (*parsers[1])(char*, void*) = {ld_str};
  void *data[1] = {&letter_id};
  bool required[1] = {true};
  bool done[1];

  parse_attributes(meme_attr_parse_error, ps, "letter_ref", attrs, 1, names, parsers, data, required, done);

  if (ps->callbacks->handle_letter_ref && ps->state != PS_ERROR) {
    ps->callbacks->handle_letter_ref(ps->user_data, letter_id);
  }
}

/*****************************************************************************
 * Handle the document start
 ****************************************************************************/
void handle_start_doc(void *ctx) {
  PS_T *ps = (PS_T*)ctx;
  meme_push_es(ps, PS_END, ES_ONCE);
  meme_push_es(ps, PS_IN_MEME, ES_ONCE);
}

/*****************************************************************************
 * Handle the document end
 ****************************************************************************/
void handle_end_doc(void *ctx) {
  PS_T *ps = (PS_T*)ctx;
  ES_T *es;
  if (ps->state != PS_ERROR) meme_update_es(ps, PS_END);
  while ((es = (ES_T*)linklst_pop(ps->expected_stack)) != NULL) {
    free(es);
  }
}

/*****************************************************************************
 * Handle the characters within elements. Only accepts characters
 * that pass the character filter. If a filter is not set then it doesn't
 * store anything.
 ****************************************************************************/
void handle_characters(void *ctx, const xmlChar *ch, int len) {
  PS_T *ps = (PS_T*)ctx;
  if (ps->state == PS_ERROR) return;
  if (ps->udepth) {//we don't know where we are!
    /*
    if (ps->callbacks->characters_unknown) {
      ps->callbacks->characters_unknown(ps->invoker_state, ch, len);
    }
    */
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
    if (meme_update_es(ps, _transition_ )) { \
      ps->state = _transition_; \
      _call_ (ps, translate_attributes(&(ps->attrbuf), attrs)); \
      ps->characters.accept = _char_accept_; \
    } \
    break; \
  }

#define CHECK_START_ELE(_expected_,_transition_,_char_accept_) \
  if (strcasecmp((char*)name, #_expected_ ) == 0) { \
    if (meme_update_es(ps, _transition_ )) { \
      ps->state = _transition_; \
      ps->characters.accept = _char_accept_; \
    } \
    break; \
  }

#define CHECK_START_ELE_ALIAS(_expected_,_alias_,_transition_,_char_accept_) \
  if (strcasecmp((char*)name, #_expected_ ) == 0 || strcasecmp((char*)name, #_alias_) == 0 ) { \
    if (meme_update_es(ps, _transition_ )) { \
      ps->state = _transition_; \
      ps->characters.accept = _char_accept_; \
    } \
    break; \
  }

/*****************************************************************************
 * Handle a starting element.
 ****************************************************************************/
void handle_start_ele(void *ctx, const xmlChar *name, const xmlChar **attrs) {
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
    allow_unknown = 0; //are unknowns allowed in this state?
    switch (ps->state) {
      case PS_START:
        DO_START_ELE(MEME, start_ele_meme, PS_IN_MEME, IGNORE);
        known = 0;
        break;
      case PS_IN_MEME:
        DO_START_ELE(training_set, start_ele_training_set, PS_IN_TRAINING_SET, IGNORE);
        DO_START_ELE(model, start_ele_model, PS_IN_MODEL, IGNORE);
        DO_START_ELE(motifs, start_ele_motifs, PS_IN_MOTIFS, IGNORE);
        DO_START_ELE(scanned_sites_summary, start_ele_scanned_sites_summary, PS_IN_SCANNED_SITES_SUMMARY, IGNORE);
        known = 0;
        break;
      case PS_IN_TRAINING_SET:
        DO_START_ELE(alphabet, start_ele_alphabet, PS_IN_ALPHABET, IGNORE);
        DO_START_ELE(ambigs, start_ele_ambigs, PS_IN_AMBIGS, IGNORE);
        DO_START_ELE(sequence, start_ele_sequence, PS_IN_SEQUENCE, IGNORE);
        DO_START_ELE(letter_frequencies, start_ele_letter_frequencies, PS_IN_LETTER_FREQUENCIES, IGNORE);
        known = 0;
        break;
      case PS_IN_MODEL:
        CHECK_START_ELE(command_line, PS_IN_COMMAND_LINE, ALL_CHARS);
        CHECK_START_ELE(host, PS_IN_HOST, ALL_BUT_SPACE);
        CHECK_START_ELE(type, PS_IN_TYPE, ALL_BUT_SPACE);
        CHECK_START_ELE(nmotifs, PS_IN_NMOTIFS, ALL_BUT_SPACE);
        CHECK_START_ELE(evalue_threshold, PS_IN_EVALUE_THRESHOLD, ALL_BUT_SPACE);
        CHECK_START_ELE(object_function, PS_IN_OBJECT_FUNCTION, ALL_CHARS);
        CHECK_START_ELE(spfun, PS_IN_SPFUN, ALL_CHARS);
        CHECK_START_ELE(use_llr, PS_IN_USE_LLR, ALL_BUT_SPACE);
        CHECK_START_ELE(min_width, PS_IN_MIN_WIDTH, ALL_BUT_SPACE);
        CHECK_START_ELE(max_width, PS_IN_MAX_WIDTH, ALL_BUT_SPACE);
	CHECK_START_ELE(minic, PS_IN_MINIC, ALL_BUT_SPACE);	// for meme_4.10.0
        CHECK_START_ELE(wg, PS_IN_WG, ALL_BUT_SPACE);
        CHECK_START_ELE(ws, PS_IN_WS, ALL_BUT_SPACE);
        CHECK_START_ELE(endgaps, PS_IN_ENDGAPS, ALL_BUT_SPACE);
        CHECK_START_ELE(substring, PS_IN_SUBSTRING, ALL_BUT_SPACE);
        CHECK_START_ELE(minsites, PS_IN_MINSITES, ALL_BUT_SPACE);
        CHECK_START_ELE(maxsites, PS_IN_MAXSITES, ALL_BUT_SPACE);
        CHECK_START_ELE(wnsites, PS_IN_WNSITES, ALL_BUT_SPACE);
        CHECK_START_ELE(prob, PS_IN_PROB, ALL_BUT_SPACE);
        CHECK_START_ELE(spmap, PS_IN_SPMAP, ALL_BUT_SPACE);
        CHECK_START_ELE(spfuzz, PS_IN_SPFUZZ, ALL_BUT_SPACE);
        CHECK_START_ELE(prior, PS_IN_PRIOR, ALL_BUT_SPACE);
        CHECK_START_ELE(beta, PS_IN_BETA, ALL_BUT_SPACE);
        CHECK_START_ELE(maxiter, PS_IN_MAXITER, ALL_BUT_SPACE);
        CHECK_START_ELE(distance, PS_IN_DISTANCE, ALL_BUT_SPACE);
        CHECK_START_ELE(num_sequences, PS_IN_NUM_SEQUENCES, ALL_BUT_SPACE);
        CHECK_START_ELE(num_positions, PS_IN_NUM_POSITIONS, ALL_BUT_SPACE);
        CHECK_START_ELE(seed, PS_IN_SEED, ALL_BUT_SPACE);
        CHECK_START_ELE(seqfrac, PS_IN_SEQFRAC, ALL_BUT_SPACE);
        CHECK_START_ELE_ALIAS(hsfrac, ctfrac, PS_IN_HSFRAC, ALL_BUT_SPACE);
        CHECK_START_ELE(searchsize, PS_IN_SEARCHSIZE, ALL_BUT_SPACE);
        CHECK_START_ELE(maxsize, PS_IN_MAXSIZE, ALL_BUT_SPACE);
        CHECK_START_ELE(norand, PS_IN_NORAND, ALL_CHARS);
        CHECK_START_ELE(csites, PS_IN_CSITES, ALL_BUT_SPACE);
        CHECK_START_ELE(maxwords, PS_IN_MAXWORDS, ALL_BUT_SPACE);
        CHECK_START_ELE(strands, PS_IN_STRANDS, ALL_BUT_SPACE);
        CHECK_START_ELE(brief, PS_IN_BRIEF, ALL_BUT_SPACE);
        CHECK_START_ELE(psp_file, PS_IN_PSP_FILE, ALL_CHARS);
        CHECK_START_ELE(priors_file, PS_IN_PRIORS_FILE, ALL_CHARS);
        CHECK_START_ELE(reason_for_stopping, PS_IN_REASON_FOR_STOPPING, ALL_CHARS);
        CHECK_START_ELE(back_order, PS_IN_BACK_ORDER, ALL_BUT_SPACE);
        DO_START_ELE(background_frequencies, start_ele_background_frequencies, PS_IN_BACKGROUND_FREQUENCIES, IGNORE);
        known = 0;
        break;
      case PS_IN_MOTIFS:
        DO_START_ELE(motif, start_ele_motif, PS_IN_MOTIF, IGNORE);
        known = 0;
        break;
      case PS_IN_SCANNED_SITES_SUMMARY:
        DO_START_ELE(scanned_sites, start_ele_scanned_sites, PS_IN_SCANNED_SITES, IGNORE);
        known = 0;
      case PS_IN_ALPHABET:
        DO_START_ELE(letter, start_ele_alphabet_letter, PS_IN_ALPHABET_LETTER, IGNORE);
        known = 0;
        break;
      case PS_IN_AMBIGS:
        DO_START_ELE(letter, start_ele_ambigs_letter, PS_IN_AMBIGS_LETTER, IGNORE);
        known = 0;
        break;
      case PS_IN_LETTER_FREQUENCIES:
        DO_START_ELE(alphabet_array, start_ele_lf_alphabet_array, PS_IN_LF_ALPHABET_ARRAY, IGNORE);
        known = 0;
        break;
      case PS_IN_BACKGROUND_FREQUENCIES:
        DO_START_ELE(alphabet_array, start_ele_bf_alphabet_array, PS_IN_BF_ALPHABET_ARRAY, IGNORE);
        known = 0;
        break;
      case PS_IN_MOTIF:
        DO_START_ELE(scores, start_ele_scores, PS_IN_SCORES, IGNORE);
        DO_START_ELE(probabilities, start_ele_probabilities, PS_IN_PROBABILITIES, IGNORE);
        CHECK_START_ELE(regular_expression, PS_IN_REGULAR_EXPRESSION, ALL_BUT_SPACE);
        DO_START_ELE(contributing_sites, start_ele_contributing_sites, PS_IN_CONTRIBUTING_SITES, IGNORE);
        known = 0;
        break;
      case PS_IN_SCANNED_SITES:
        DO_START_ELE(scanned_site, start_ele_scanned_site, PS_IN_SCANNED_SITE, IGNORE);
        known = 0;
        break;
      case PS_IN_LF_ALPHABET_ARRAY:
        DO_START_ELE(value, start_ele_value, PS_IN_LF_AA_VALUE, ALL_BUT_SPACE);
        known = 0;
        break;
      case PS_IN_BF_ALPHABET_ARRAY:
        DO_START_ELE(value, start_ele_value, PS_IN_BF_AA_VALUE, ALL_BUT_SPACE);
        known = 0;
        break;
      case PS_IN_SCORES:
        DO_START_ELE(alphabet_matrix, start_ele_sc_alphabet_matrix, PS_IN_SC_ALPHABET_MATRIX, IGNORE);
        known = 0;
        break;
      case PS_IN_PROBABILITIES:
        DO_START_ELE(alphabet_matrix, start_ele_pr_alphabet_matrix, PS_IN_PR_ALPHABET_MATRIX, IGNORE);
        known = 0;
        break;
      case PS_IN_CONTRIBUTING_SITES:
        DO_START_ELE(contributing_site, start_ele_contributing_site, PS_IN_CONTRIBUTING_SITE, IGNORE);
        known = 0;
        break;
      case PS_IN_SC_ALPHABET_MATRIX:
        DO_START_ELE(alphabet_array, start_ele_sc_am_alphabet_array, PS_IN_SC_AM_ALPHABET_ARRAY, IGNORE);
        known = 0;
        break;
      case PS_IN_PR_ALPHABET_MATRIX:
        DO_START_ELE(alphabet_array, start_ele_pr_am_alphabet_array, PS_IN_PR_AM_ALPHABET_ARRAY, IGNORE);
        known = 0;
        break;
      case PS_IN_CONTRIBUTING_SITE:
        CHECK_START_ELE(left_flank, PS_IN_LEFT_FLANK, ALL_BUT_SPACE);
        DO_START_ELE(site, start_ele_site, PS_IN_SITE, IGNORE);
        CHECK_START_ELE(right_flank, PS_IN_RIGHT_FLANK, ALL_BUT_SPACE);
        known = 0;
        break;
      case PS_IN_SC_AM_ALPHABET_ARRAY:
        DO_START_ELE(value, start_ele_value, PS_IN_SC_AM_AA_VALUE, ALL_BUT_SPACE);
        known = 0;
        break;
      case PS_IN_PR_AM_ALPHABET_ARRAY:
        DO_START_ELE(value, start_ele_value, PS_IN_PR_AM_AA_VALUE, ALL_BUT_SPACE);
        known = 0;
        break;
      case PS_IN_SITE:
        DO_START_ELE(letter_ref, start_ele_letter_ref, PS_IN_LETTER_REF, IGNORE);
        known = 0;
        break;
      default:
        known = 0;
    }
    if (!known) {
      if (allow_unknown) {
        ps->udepth = 1;
        /*
        if (ps->callbacks->start_unknown) {
          ps->callbacks->start_unknown(ps->user_data, name, attrs);
        }
        */
      } else {
        error(ps, "MEME IO XML parser encountered illegal tag %s while in state %s\n",(char*)name, state_names[ps->state]);
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
#define DO_END_ELE_ALIAS(_expected_,_alias_,_call_,_state_,_transition_) \
  case _state_: \
    if (strcasecmp((char*)name, #_expected_) == 0 || strcasecmp((char*)name, #_alias_) == 0) { \
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
void handle_end_ele(void *ctx, const xmlChar *name) {
  PS_T *ps = (PS_T*)ctx;
  int known = 0;
  if (ps->state == PS_ERROR) return;
  if (ps->udepth) {
    ps->udepth -= 1; 
    /*
    if (ps->callbacks->end_unknown) {
      ps->callbacks->end_unknown(ps->user_data, name);
    }
    */
  } else {
    switch (ps->state) {
      DO_END_ELE(MEME, end_ele_meme, PS_IN_MEME, PS_END);
      DO_END_ELE(training_set, end_ele_training_set, PS_IN_TRAINING_SET, PS_IN_MEME);
      DO_END_ELE(model, end_ele_model, PS_IN_MODEL, PS_IN_MEME);
      DO_END_ELE(motifs, end_ele_motifs, PS_IN_MOTIFS, PS_IN_MEME);
      DO_END_ELE(scanned_sites_summary, end_ele_scanned_sites_summary, PS_IN_SCANNED_SITES_SUMMARY, PS_IN_MEME);
      DO_END_ELE(alphabet, end_ele_alphabet, PS_IN_ALPHABET, PS_IN_TRAINING_SET);
      DO_END_ELE(ambigs, end_ele_ambigs, PS_IN_AMBIGS, PS_IN_TRAINING_SET);
      CHECK_END_ELE(sequence, PS_IN_SEQUENCE, PS_IN_TRAINING_SET);
      DO_END_ELE(letter_frequencies, end_ele_letter_frequencies, PS_IN_LETTER_FREQUENCIES, PS_IN_TRAINING_SET);
      DO_END_ELE(command_line, end_ele_command_line, PS_IN_COMMAND_LINE, PS_IN_MODEL);
      DO_END_ELE(host, end_ele_host, PS_IN_HOST, PS_IN_MODEL);
      DO_END_ELE(type, end_ele_type, PS_IN_TYPE, PS_IN_MODEL);
      DO_END_ELE(nmotifs, end_ele_nmotifs, PS_IN_NMOTIFS, PS_IN_MODEL);
      DO_END_ELE(evalue_threshold, end_ele_evalue_threshold, PS_IN_EVALUE_THRESHOLD, PS_IN_MODEL);
      DO_END_ELE(object_function, end_ele_object_function, PS_IN_OBJECT_FUNCTION, PS_IN_MODEL);
      DO_END_ELE(spfun, end_ele_spfun, PS_IN_SPFUN, PS_IN_MODEL);
      DO_END_ELE(use_llr, end_ele_use_llr, PS_IN_USE_LLR, PS_IN_MODEL);
      DO_END_ELE(min_width, end_ele_min_width, PS_IN_MIN_WIDTH, PS_IN_MODEL);
      DO_END_ELE(max_width, end_ele_max_width, PS_IN_MAX_WIDTH, PS_IN_MODEL);
      DO_END_ELE(minic, end_ele_minic, PS_IN_MINIC, PS_IN_MODEL);
      DO_END_ELE(wg, end_ele_wg, PS_IN_WG, PS_IN_MODEL);
      DO_END_ELE(ws, end_ele_ws, PS_IN_WS, PS_IN_MODEL);
      DO_END_ELE(endgaps, end_ele_endgaps, PS_IN_ENDGAPS, PS_IN_MODEL);
      DO_END_ELE(substring, end_ele_substring, PS_IN_SUBSTRING, PS_IN_MODEL);
      DO_END_ELE(minsites, end_ele_minsites, PS_IN_MINSITES, PS_IN_MODEL);
      DO_END_ELE(maxsites, end_ele_maxsites, PS_IN_MAXSITES, PS_IN_MODEL);
      DO_END_ELE(wnsites, end_ele_wnsites, PS_IN_WNSITES, PS_IN_MODEL);
      DO_END_ELE(prob, end_ele_prob, PS_IN_PROB, PS_IN_MODEL);
      DO_END_ELE(spmap, end_ele_spmap, PS_IN_SPMAP, PS_IN_MODEL);
      DO_END_ELE(spfuzz, end_ele_spfuzz, PS_IN_SPFUZZ, PS_IN_MODEL);
      DO_END_ELE(prior, end_ele_prior, PS_IN_PRIOR, PS_IN_MODEL);
      DO_END_ELE(beta, end_ele_beta, PS_IN_BETA, PS_IN_MODEL);
      DO_END_ELE(maxiter, end_ele_maxiter, PS_IN_MAXITER, PS_IN_MODEL);
      DO_END_ELE(distance, end_ele_distance, PS_IN_DISTANCE, PS_IN_MODEL);
      DO_END_ELE(num_sequences, end_ele_num_sequences, PS_IN_NUM_SEQUENCES, PS_IN_MODEL);
      DO_END_ELE(num_positions, end_ele_num_positions, PS_IN_NUM_POSITIONS, PS_IN_MODEL);
      DO_END_ELE(seed, end_ele_seed, PS_IN_SEED, PS_IN_MODEL);
      DO_END_ELE(seqfrac, end_ele_seqfrac, PS_IN_SEQFRAC, PS_IN_MODEL);
      DO_END_ELE_ALIAS(hsfrac, ctfrac, end_ele_hsfrac, PS_IN_HSFRAC, PS_IN_MODEL);
      DO_END_ELE(searchsize, end_ele_searchsize, PS_IN_SEARCHSIZE, PS_IN_MODEL);
      DO_END_ELE(maxwords, end_ele_maxwords, PS_IN_MAXWORDS, PS_IN_MODEL);
      DO_END_ELE(maxsize, end_ele_maxsize, PS_IN_MAXSIZE, PS_IN_MODEL);
      DO_END_ELE(norand, end_ele_norand, PS_IN_NORAND, PS_IN_MODEL);
      DO_END_ELE(csites, end_ele_csites, PS_IN_CSITES, PS_IN_MODEL);
      DO_END_ELE(strands, end_ele_strands, PS_IN_STRANDS, PS_IN_MODEL);
      DO_END_ELE(brief, end_ele_brief, PS_IN_BRIEF, PS_IN_MODEL);
      DO_END_ELE(psp_file, end_ele_psp_file, PS_IN_PSP_FILE, PS_IN_MODEL);
      DO_END_ELE(priors_file, end_ele_priors_file, PS_IN_PRIORS_FILE, PS_IN_MODEL);
      DO_END_ELE(reason_for_stopping, end_ele_reason_for_stopping, PS_IN_REASON_FOR_STOPPING, PS_IN_MODEL);
      DO_END_ELE(back_order, end_ele_back_order, PS_IN_BACK_ORDER, PS_IN_MODEL);
      DO_END_ELE(background_frequencies, end_ele_background_frequencies, PS_IN_BACKGROUND_FREQUENCIES, PS_IN_MODEL);
      DO_END_ELE(motif, end_ele_motif, PS_IN_MOTIF, PS_IN_MOTIFS);
      DO_END_ELE(scanned_sites, end_ele_scanned_sites, PS_IN_SCANNED_SITES, PS_IN_SCANNED_SITES_SUMMARY);
      CHECK_END_ELE(letter, PS_IN_ALPHABET_LETTER, PS_IN_ALPHABET);
      CHECK_END_ELE(letter, PS_IN_AMBIGS_LETTER, PS_IN_AMBIGS);
      DO_END_ELE(alphabet_array, end_ele_lf_alphabet_array, PS_IN_LF_ALPHABET_ARRAY, PS_IN_LETTER_FREQUENCIES);
      DO_END_ELE(alphabet_array, end_ele_bf_alphabet_array, PS_IN_BF_ALPHABET_ARRAY, PS_IN_BACKGROUND_FREQUENCIES);
      DO_END_ELE(scores, end_ele_scores, PS_IN_SCORES, PS_IN_MOTIF);
      DO_END_ELE(probabilities, end_ele_probabilities, PS_IN_PROBABILITIES, PS_IN_MOTIF);
      DO_END_ELE(regular_expression, end_ele_regular_expression, PS_IN_REGULAR_EXPRESSION, PS_IN_MOTIF);
      DO_END_ELE(contributing_sites, end_ele_contributing_sites, PS_IN_CONTRIBUTING_SITES, PS_IN_MOTIF);
      CHECK_END_ELE(scanned_site, PS_IN_SCANNED_SITE, PS_IN_SCANNED_SITES);
      DO_END_ELE(value, end_ele_lf_aa_value, PS_IN_LF_AA_VALUE, PS_IN_LF_ALPHABET_ARRAY);
      DO_END_ELE(value, end_ele_bf_aa_value, PS_IN_BF_AA_VALUE, PS_IN_BF_ALPHABET_ARRAY);
      DO_END_ELE(alphabet_matrix, end_ele_sc_alphabet_matrix, PS_IN_SC_ALPHABET_MATRIX, PS_IN_SCORES);
      DO_END_ELE(alphabet_matrix, end_ele_pr_alphabet_matrix, PS_IN_PR_ALPHABET_MATRIX, PS_IN_PROBABILITIES);
      DO_END_ELE(contributing_site, end_ele_contributing_site, PS_IN_CONTRIBUTING_SITE, PS_IN_CONTRIBUTING_SITES);
      DO_END_ELE(alphabet_array, end_ele_sc_am_alphabet_array, PS_IN_SC_AM_ALPHABET_ARRAY, PS_IN_SC_ALPHABET_MATRIX);
      DO_END_ELE(alphabet_array, end_ele_pr_am_alphabet_array, PS_IN_PR_AM_ALPHABET_ARRAY, PS_IN_PR_ALPHABET_MATRIX);
      DO_END_ELE(left_flank, end_ele_left_flank, PS_IN_LEFT_FLANK, PS_IN_CONTRIBUTING_SITE);
      DO_END_ELE(site, end_ele_site, PS_IN_SITE, PS_IN_CONTRIBUTING_SITE);
      DO_END_ELE(right_flank, end_ele_right_flank, PS_IN_RIGHT_FLANK, PS_IN_CONTRIBUTING_SITE);
      DO_END_ELE(value, end_ele_sc_am_aa_value, PS_IN_SC_AM_AA_VALUE, PS_IN_SC_AM_ALPHABET_ARRAY);
      DO_END_ELE(value, end_ele_pr_am_aa_value, PS_IN_PR_AM_AA_VALUE, PS_IN_PR_AM_ALPHABET_ARRAY);
      CHECK_END_ELE(letter_ref, PS_IN_LETTER_REF, PS_IN_SITE);
      case PS_START: // satisify compiler
      case PS_END:
      case PS_ERROR:
      break;
    }
    if (!known) {
      error(ps, "MEME IO XML parser error: unexpected end tag %s in state %s\n", (char*)name, state_names[ps->state]);
    }
  }
}


/*****************************************************************************
 * Register handlers on the xmlSAXHandler structure
 ****************************************************************************/
void register_meme_io_xml_sax_handlers(xmlSAXHandler *handler) {
  memset(handler, 0, sizeof(xmlSAXHandler));
  handler->startDocument = handle_start_doc;
  handler->endDocument = handle_end_doc;
  handler->characters = handle_characters;
  handler->startElement = handle_start_ele;
  handler->endElement = handle_end_ele;
}


/*****************************************************************************
 * Creates the data to be passed to the SAX parser
 ****************************************************************************/
void* create_meme_io_xml_sax_context(void *user_data, MEME_IO_XML_CALLBACKS_T *callbacks) {
  PS_T *ps;
  CHARBUF_T *buf;
  ps = (PS_T*)mm_malloc(sizeof(PS_T));
  memset(ps, 0, sizeof(PS_T));
  ps->udepth = 0;
  ps->state = PS_START;
  ps->callbacks = callbacks;
  ps->user_data = user_data;
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
void destroy_meme_io_xml_sax_context(void *ctx) {
  PS_T *ps = (PS_T*)ctx;
  attrbuf_free(&(ps->attrbuf));
  free(ps->characters.buffer);
  free(ps->letter_id_buf);
  linklst_destroy_all(ps->expected_stack, free);
  free(ps);
}

