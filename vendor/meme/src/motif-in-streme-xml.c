#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include "streme-sax.h"
#include "linked-list.h"
#include "motif-in-streme-xml.h"
#include "motif-spec.h"
#include "string-match.h"
#include "utils.h"
#include "alph-in.h"
#include "alphabet.h"


struct fscope {
  int options_found;
  int options_returned;
  ALPH_T* alphabet;
  STREME_STRANDS_EN strands;
  ARRAY_T *background;
};

typedef struct streme_ctx CTX_T;
struct streme_ctx {
  int options;
  ALPH_READER_T *alph_rdr;
  struct fscope fscope;
  MOTIF_T *motif;
  LINKLST_T *motif_queue;
  LINKLST_T *errors;
  LINKLST_T *warnings;
  short file_type_match;
};

typedef struct sxml DXML_T;
struct sxml {
  CTX_T *data;
  STREME_IO_XML_CALLBACKS_T *callbacks;
  void *sax_context;
  xmlSAXHandler *handler;
  xmlParserCtxt *ctxt;
};

static CTX_T* create_parser_data(int options, const char *optional_filename) {
  CTX_T *data;
  data = mm_malloc(sizeof(CTX_T));
  memset(data, 0, sizeof(CTX_T));
  data->options = options;
  data->motif_queue = linklst_create();
  data->errors = linklst_create();
  data->warnings = linklst_create();
  data->file_type_match = file_name_match("streme", "xml", optional_filename);
  return data;
}

static void destroy_parser_data(CTX_T *data) {
  linklst_destroy_all(data->motif_queue, destroy_motif);
  linklst_destroy_all(data->errors, free);
  linklst_destroy_all(data->warnings, free);
  if (data->alph_rdr) {
    alph_reader_destroy(data->alph_rdr);
  }
  if (data->fscope.alphabet) {
    alph_release(data->fscope.alphabet);
  }
  if (data->fscope.background) {
    free_array(data->fscope.background);
  }
  if (data->motif) destroy_motif(data->motif);
  memset(data, 0, sizeof(CTX_T));
  free(data);
}
/*
 * Callbacks
 */
void sxml_error(void *ctx, const char *fmt, va_list ap) {
  CTX_T *data = (CTX_T*)ctx;
  int len;
  char *msg;
  va_list ap_copy;
  va_copy(ap_copy, ap);
  len = vsnprintf(NULL, 0, fmt, ap_copy);
  va_end(ap_copy);
  msg = mm_malloc(sizeof(char) * (len + 1));
  vsnprintf(msg, len + 1, fmt, ap);
  linklst_add(msg, data->errors);
}

/*****************************************************************************
 * As the other error function takes va_list this is a simply to
 * enable calling the real error function
 ****************************************************************************/
static void local_error(CTX_T *data, char *format, ...) {
  va_list  argp;
  va_start(argp, format);
  sxml_error(data, format, argp);
  va_end(argp);
}

/*****************************************************************************
 * Stores a parser warning message for return
 ****************************************************************************/
static void local_warning(CTX_T *data, char *format, ...) {
  va_list  argp;
  int len;
  char *msg;
  va_start(argp, format);
  len = vsnprintf(NULL, 0, format, argp);
  va_end(argp);
  msg = mm_malloc(sizeof(char) * (len + 1));
  va_start(argp, format);
  vsnprintf(msg, len + 1, format, argp);
  va_end(argp);
  linklst_add(msg, data->warnings);
}

void sxml_start_streme(void *ctx, int major, int minor, int patch, char *release_date) {
  CTX_T *data;
  data = (CTX_T*)ctx;
  data->file_type_match = 3; // probably a streme xml file
}

/*****************************************************************************
 * STREME > model > alphabet
 * Read in the alphabet name.
 ****************************************************************************/
void sxml_start_alphabet(void *ctx, char *name, int extends_flag) {
  CTX_T *data = (CTX_T*)ctx;
  data->alph_rdr = alph_reader_create();
  alph_reader_header(data->alph_rdr, 1, name, extends_flag);
}

/*****************************************************************************
 * STREME > model > alphabet > letter
 * Read in a identifier and symbol pair for a letter.
 ****************************************************************************/
void sxml_alphabet_letter(void *ctx, char *id, char symbol, char *aliases, char complement, char *equals, char *name, int colour) {
  CTX_T *data;
  char sym[2];
  data = (CTX_T*)ctx;
  if (equals == NULL) {
    alph_reader_core(data->alph_rdr, symbol, aliases, name, colour, complement);
  } else {
    alph_reader_ambig(data->alph_rdr, symbol, aliases, name, colour, equals);
  }
}

/*****************************************************************************
 * STREME > model > strands
 * Read in the strand handling method
 ****************************************************************************/
void sxml_handle_strands(void *ctx, STREME_STRANDS_EN strands) {
  CTX_T *data = (CTX_T*)ctx;
  data->fscope.strands = strands;
}

/*****************************************************************************
 * STREME > model > /alphabet
 * Read in the number of symbols in the alphabet and if it is nucleotide or 
 * amino-acid (RNA is apparently classed as nucleotide).
 ****************************************************************************/
void sxml_end_alphabet(void *ctx) {
  PARMSG_T *message;
  CTX_T *data;
  RBNODE_T *node;
  char *id, symbol;
  bool *exists;
  int i;

  data = (CTX_T*)ctx;
  alph_reader_done(data->alph_rdr);
  // report any errors that the alphabet reader found
  while (alph_reader_has_message(data->alph_rdr)) {
    message = alph_reader_next_message(data->alph_rdr);
    if (message->severity == SEVERITY_ERROR) {
      local_error(data, "Alphabet error: %s.\n", message->message);
    } else {
      local_warning(data, "Alphabet warning: %s.\n", message->message);
    }
    parmsg_destroy(message);
  }
  // try to get an alphabet
  data->fscope.alphabet = alph_reader_alphabet(data->alph_rdr);
  alph_reader_destroy(data->alph_rdr);
  data->alph_rdr = NULL;
}

void sxml_handle_background(void *ctx, int nfreqs, double *freqs, STREME_BG_EN source) {
  CTX_T *data;
  MOTIF_T *motif;
  ARRAY_T *bg;
  int i;

  data = (CTX_T*)ctx;
  data->file_type_match = 4; // it must be a streme xml file!
  data->fscope.background = allocate_array(nfreqs);
  bg = data->fscope.background;
  for (i = 0; i < nfreqs; i++) set_array_item(i, freqs[i], bg);
}

void sxml_start_motif(void *ctx, 
  char *id, char *alt, int width, int initial_width, char *seed, double score_threshold,
  long train_pos_count, long train_neg_count, double train_log_pvalue, char *train_pvalue, double train_dtc, double train_bernoulli,
  long test_pos_count, long test_neg_count, double test_log_pvalue, char *test_pvalue, double test_dtc, double test_bernoulli,
  char *is_palindromic, double elapsed_time) 
{
  CTX_T *data;
  MOTIF_T *motif;

  data = (CTX_T*)ctx;
  data->motif = (MOTIF_T*)mm_malloc(sizeof(MOTIF_T));
  motif = data->motif;
  memset(motif, 0, sizeof(MOTIF_T));
  set_motif_id(id, strlen(id), motif);
  set_motif_id2(alt, strlen(alt), motif);
  set_motif_strand('+', motif);
  motif->length = width;
  motif->num_sites = train_pos_count + test_pos_count;
  // Use the test p-value unless there was no test set.
  motif->log_evalue = (test_pos_count+test_neg_count > 0) ? test_log_pvalue : train_log_pvalue;
  motif->evalue = (test_pos_count+test_neg_count > 0) ? pow(10, test_log_pvalue) : pow(10, train_log_pvalue);
  // calculate alphabet size
  motif->alph = alph_hold(data->fscope.alphabet);
  motif->flags = (data->fscope.strands == 2 ? MOTIF_BOTH_STRANDS : 0);
  // allocate the matrix
  motif->freqs = allocate_matrix(motif->length, alph_size_core(motif->alph));
  motif->scores = NULL; // no scores in STREME xml
  // no url in STREME
  motif->url = strdup("");
  // set by postprocessing
  motif->complexity = -1;
  motif->trim_left = 0;
  motif->trim_right = 0;
}

void sxml_end_motif(void *ctx) {
  CTX_T *data;
  MOTIF_T *motif;

  data = (CTX_T*)ctx;
  motif = data->motif;

  linklst_add(motif, data->motif_queue);
  data->motif = NULL;
}

void sxml_handle_pos(void *ctx, int pos, int nfreqs, double *freqs) {
  CTX_T *data;
  MOTIF_T *motif;
  ARRAY_T *row;
  int i;

  data = (CTX_T*)ctx;
  motif = data->motif;
  row = get_matrix_row(pos - 1, motif->freqs);
  for (i = 0; i < nfreqs; i++) set_array_item(i, freqs[i], row);
}

void* sxml_create(const char *optional_filename, int options) {
  DXML_T *parser;
  // allocate the structure to hold the parser information
  parser = (DXML_T*)mm_malloc(sizeof(DXML_T));
  // setup data store
  parser->data = create_parser_data(options, optional_filename);
  // setup meme xml callbacks
  parser->callbacks = (STREME_IO_XML_CALLBACKS_T*)mm_malloc(sizeof(STREME_IO_XML_CALLBACKS_T));
  memset(parser->callbacks, 0, sizeof(STREME_IO_XML_CALLBACKS_T));
  //start callbacks
  parser->callbacks->error = sxml_error;
  parser->callbacks->start_streme = sxml_start_streme;
  // start streme
  // start model
  parser->callbacks->start_alphabet = sxml_start_alphabet;
  parser->callbacks->end_alphabet = sxml_end_alphabet;
  parser->callbacks->handle_alphabet_letter = sxml_alphabet_letter;
  parser->callbacks->handle_strands = sxml_handle_strands;
  parser->callbacks->handle_background = sxml_handle_background;
  // end model
  // start motifs
  parser->callbacks->start_motif = sxml_start_motif;
  parser->callbacks->end_motif = sxml_end_motif;
  // start motif
  parser->callbacks->handle_pos = sxml_handle_pos;
  // end motif
  // end motifs
  // end streme
  // end callbacks
  // setup context for SAX parser
  parser->sax_context = create_streme_io_xml_sax_context(parser->data, parser->callbacks);
  // setup SAX handler
  parser->handler = (xmlSAXHandler*)mm_malloc(sizeof(xmlSAXHandler));
  register_streme_io_xml_sax_handlers(parser->handler);
  // create the push parser context
  parser->ctxt = xmlCreatePushParserCtxt(parser->handler, parser->sax_context, NULL, 0, optional_filename); 
  return parser;
}

void sxml_destroy(void *data) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  // destroy the push parser context
  xmlFreeParserCtxt(parser->ctxt);
  // destroy the push parser callbacks
  free(parser->handler);
  // destroy the sax parser context
  destroy_streme_io_xml_sax_context(parser->sax_context);
  // destroy the sax parser callbacks
  free(parser->callbacks);
  // destroy the data store
  destroy_parser_data(parser->data);
  // destroy the parser information
  free(parser);
}

void sxml_update(void *data, const char *chunk, size_t size, short end) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  xmlParseChunk(parser->ctxt, chunk, size, end);
}

short sxml_has_motif(void *data) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  return linklst_size(parser->data->motif_queue) > 0;
}

short sxml_has_format_match(void *data) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  return parser->data->file_type_match;
}

/*****************************************************************************
 * Check if the parser has found a problem with the file that has not already
 * been returned.
 ****************************************************************************/
short sxml_has_warning(void *data) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  return linklst_size(parser->data->warnings) > 0;
}

/*****************************************************************************
 * Get the next human readable warning message detailing the 
 * problem with the file. The caller is responsible for freeing the memory.
 ****************************************************************************/
char* sxml_next_warning(void *data) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  return (char*)linklst_pop(parser->data->warnings);
}

short sxml_has_error(void *data) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  return linklst_size(parser->data->errors) > 0;
}

char* sxml_next_error(void *data) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  return linklst_pop(parser->data->errors);
}

MOTIF_T* sxml_next_motif(void *data) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  return linklst_pop(parser->data->motif_queue);
}

/*
 * WARNING, does not "hold" alphabet. If users wish to keep the alphabet they
 * must call alph_hold on it themselves.
 */
ALPH_T* sxml_get_alphabet(void *data) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  return parser->data->fscope.alphabet;
}

int sxml_get_strands(void *data) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  return (parser->data->fscope.strands == STREME_STRANDS_BOTH ? 2 : 1);
}

bool sxml_get_bg(void *data, ARRAY_T **bg) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  if (parser->data->fscope.background == NULL) return false;
  *bg = resize_array(*bg, get_array_length(parser->data->fscope.background));
  copy_array(parser->data->fscope.background, *bg);
  return true;
}

void* sxml_motif_optional(void *data, int option) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  if (!(parser->data->options & option)) {
    die("Requested value of optional component which was not requested "
        "during construction.\n");
    return NULL;
  }
  // STREME xml does not currently support optional components
  return NULL;
}

void* sxml_file_optional(void *data, int option) {
  DXML_T *parser;
  parser = (DXML_T*)data;
  if (!(parser->data->options & option)) {
    die("Requested value of optional component which was not requested "
        "during construction.\n");
    return NULL;
  }
  if (parser->data->fscope.options_found & option) {
    if (parser->data->fscope.options_returned & option) {
      die("Sorry, optional values are only returned once. "
          "This is because we cannot guarantee that the "
          "previous caller did not deallocate the memory. "
          "Hence this is a feature to avoid memory bugs.\n");
      return NULL;
    }
    parser->data->fscope.options_returned |= option;
  } else {
    // Not yet found or unsupported
    return NULL;
  }
  switch (option) {
    default:
      die("Option code %d does not match any single option. "
          "This means that it must contain multiple options "
          "and only one is allowed at a time\n.", option);
  }
  return NULL; //unreachable! make compiler happy
}
