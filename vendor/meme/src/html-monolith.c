/*
 * FILE: html-monolith.c
 * AUTHOR: James Johnson 
 * CREATE DATE: 30-November-2011
 * PROJECT: shared
 * COPYRIGHT: UQ, 2009
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: It's easier for the user to handle one html results file than 
 * a html file, a couple of javascript files and a couple of css files. However 
 * when developing it is easier to keep the files seperate. This provides the
 * tools to merge all the disparate resources into one file as well as easily
 * generate JSON data to subsitiute for an included script file at run time.
 */
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>

#include "buffer.h"
#include "html-monolith.h"
#include "red-black-tree.h"
#include "regex-utils.h"
#include "utils.h"

struct monolith_writer {
  char* search_dir;
  FILE* src;
  FILE* dest;
  bool autoclose;
  bool minimize;
  RBTREE_T* replace;
  const char* variable;
  JSONWR_T* json_writer;
  BUF_T* buffer;
  BMSTR_T* landmarks[3];
  regex_t source_re;
};

#define CHUNK 100

const char *SCRIPT_START = "<script src=\"";
const char *SCRIPT_END = "\"></script>";
const char *LINK_START = "<link rel=\"stylesheet\" type=\"text/css\" href=\"";
const char *LINK_END = "\">";
const char *IMG_START = "<img src=\"";
const char *IMG_END = "\"";

const char *JSON_VAR = "//@JSON_VAR";
const char *SCRIPT_OPEN = "<script type=\"text/javascript\">";
const char *SCRIPT_CLOSE = "</script>";
const char *STYLE_OPEN = "<style type=\"text/css\">";
const char *STYLE_CLOSE = "</style>";

const char *SOURCE_PAT = "[a-zA-Z_-]+\\.(js|css|png)$";

static char b64_enc[] = {
  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
  'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
  'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
  'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
  'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
  'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
  'w', 'x', 'y', 'z', '0', '1', '2', '3',
  '4', '5', '6', '7', '8', '9', '+', '/'
};

/*
 * Return true if the letter is a double quote
 */
static int is_quote(void *ignored, int letter) {
  return letter == '"';
}

/*
 * Echo the file to the dest
 */
static void echo_file(HTMLWR_T* htmlwr, bool is_style, char *secondary_source) {
  char buffer[CHUNK];
  size_t size_read;
  char *src_path;
  FILE *src_file;
  src_path = make_path_to_file(htmlwr->search_dir, secondary_source);
  src_file = fopen(src_path, "r");
  if (src_file) {
    fprintf(htmlwr->dest, "%s\n", (is_style ? STYLE_OPEN : SCRIPT_OPEN));
    while ((size_read = fread(buffer, sizeof(char), CHUNK, src_file)) != 0) {
      fwrite(buffer, sizeof(char), size_read, htmlwr->dest);
    }
    fclose(src_file);
    fprintf(htmlwr->dest, "\n    %s", (is_style ? STYLE_CLOSE : SCRIPT_CLOSE));
  } else {
    // can't find it, issue a warning
    fprintf(stderr, "Can not open included file \"%s\".\n", src_path);
  }
  free(src_path);
}

/*
 * Echo the file to the dest in base64 with the specified mime type.
 */
static void base64enc_file(HTMLWR_T* htmlwr, char *mimetype, char *binary_source) {
  char *src_path;
  FILE *src_file;
  int32_t bits24;
  int i, c, pad;
  src_path = make_path_to_file(htmlwr->search_dir, binary_source);
  src_file = fopen(src_path, "r");
  if (src_file) {
    fprintf(htmlwr->dest, "data:%s;base64,", mimetype);
    while (!feof(src_file) && !ferror(src_file)) {
      // create a number containing 3 bytes
      bits24 = 0;
      pad = 0;
      for (i = 0; i < 3; i++) {
        c = fgetc(src_file);
        if (c == EOF) {
          pad = 3 - i;
          bits24 = bits24 << (8 * pad);
          break;
        }
        bits24 = (bits24 << 8) | c;
      }
      switch (pad) {
        case 0:
          fputc(b64_enc[(bits24 >> 18) & 0x3F], htmlwr->dest);
          fputc(b64_enc[(bits24 >> 12) & 0x3F], htmlwr->dest);
          fputc(b64_enc[(bits24 >> 6) & 0x3F], htmlwr->dest);
          fputc(b64_enc[bits24 & 0x3F], htmlwr->dest);
          break;
        case 1:
          fputc(b64_enc[(bits24 >> 18) & 0x3F], htmlwr->dest);
          fputc(b64_enc[(bits24 >> 12) & 0x3F], htmlwr->dest);
          fputc(b64_enc[(bits24 >> 6) & 0x3F], htmlwr->dest);
          fputc('=', htmlwr->dest);
          break;
        case 2:
          fputc(b64_enc[(bits24 >> 18) & 0x3F], htmlwr->dest);
          fputc(b64_enc[(bits24 >> 12) & 0x3F], htmlwr->dest);
          fputc('=', htmlwr->dest);
          fputc('=', htmlwr->dest);
          break;
      }
    }
    if (ferror(src_file)) {
      fprintf(stderr, "Error reading included file \"%s\".\n", src_path);
    }
    fclose(src_file);
  } else {
    // can't find it, issue a warning
    fprintf(stderr, "Can not open included file \"%s\".\n", src_path);
  }
  free(src_path);
}

/*
 * htmlwr_create
 * Creates a html monolith writer.
 * Takes the directory to search for sources, the name of the primary source.
 * Returns null if the primary source does not exist.
 */
HTMLWR_T* htmlwr_create(const char* search_dir, const char* primary_source, bool minimize) {
  HTMLWR_T *htmlwr;
  char *src_path;
  FILE *src_file;
  src_path = make_path_to_file(search_dir, primary_source);
  src_file = fopen(src_path, "r");
  free(src_path);
  if (!src_file) return NULL;
  htmlwr = mm_malloc(sizeof(HTMLWR_T));
  memset(htmlwr, 0, sizeof(HTMLWR_T));
  htmlwr->minimize = minimize;
  htmlwr->search_dir = strdup(search_dir);
  htmlwr->src = src_file;
  htmlwr->replace = rbtree_create(rbtree_strcmp, rbtree_strcpy, free, rbtree_strcpy, free);
  htmlwr->buffer = buf_create(100);
  buf_flip(htmlwr->buffer); // prepare for read
  htmlwr->landmarks[0] = bmstr_create(SCRIPT_START);
  htmlwr->landmarks[1] = bmstr_create(LINK_START);
  htmlwr->landmarks[2] = bmstr_create(IMG_START);
  regcomp_or_die("source name", &(htmlwr->source_re), SOURCE_PAT, REG_EXTENDED | REG_NOSUB);
  return htmlwr;
}

/*
 * htmlwr_destroy
 * Destroys a html monolith writer.
 */
void htmlwr_destroy(HTMLWR_T* htmlwr) {
  free(htmlwr->search_dir);
  fclose(htmlwr->src);
  rbtree_destroy(htmlwr->replace);
  buf_destroy(htmlwr->buffer);
  bmstr_destroy(htmlwr->landmarks[0]);
  bmstr_destroy(htmlwr->landmarks[1]);
  bmstr_destroy(htmlwr->landmarks[2]);
  regfree(&(htmlwr->source_re));
  if (htmlwr->json_writer) jsonwr_close(htmlwr->json_writer);
  if (htmlwr->autoclose && htmlwr->dest) fclose(htmlwr->dest);
  memset(htmlwr, 0, sizeof(HTMLWR_T));
  free(htmlwr);
}

/*
 * htmlwr_set_dest_name
 * Sets the output directory and file name.
 */
void htmlwr_set_dest_name(HTMLWR_T* htmlwr, const char* directory, const char* name) {
  char *dest_path;
  dest_path = make_path_to_file(directory, name);
  htmlwr_set_dest_path(htmlwr, dest_path);
  free(dest_path);
}

/*
 * htmlwr_set_dest_path
 * Sets the full path to the output file.
 */
void htmlwr_set_dest_path(HTMLWR_T* htmlwr, const char* path) {
  FILE *file;
  file = fopen(path, "w");
  if (!file) die("Failed to open destination file \"%s\" %s", path, strerror(errno));
  htmlwr_set_dest_file(htmlwr, file, true);
}

/*
 * htmlwr_set_dest_file
 * Sets the output file and if it should be closed automatically.
 */
void htmlwr_set_dest_file(HTMLWR_T* htmlwr, FILE* file, bool close) {
  if (htmlwr->dest) die("Destination file already set.\n");
  htmlwr->dest = file;
  htmlwr->autoclose = close;
}

/*
 * htmlwr_replace
 * Tell the html monolith writer to stop when it encounters the replaced script
 * and let the user write a Javascript Object Notation variable in its place.
 */
void htmlwr_replace(HTMLWR_T* htmlwr, const char* replaced_script, const char* variable_name) {
  rbtree_put(htmlwr->replace, (void*)replaced_script, (void*)variable_name);
}

/*
 * htmlwr_replacing
 * Returns the variable name of the json object being written.
 */
const char* htmlwr_replacing(HTMLWR_T* htmlwr) {
  return htmlwr->variable;
}

/*
 * htmlwr_output
 * Replaces each <script>file.js</script> with the contents of file.js.
 * Writes until encountering the next script to be replaced with JSON or the end
 * of the primary source. Returns the JSON writer for the next replacement or
 * NULL on EOF.
 */
JSONWR_T* htmlwr_output(HTMLWR_T* htmlwr) {
  int which;
  bool bad_name;
  char *name;
  const char *before, *after;

  if (htmlwr->src == NULL) die("Source file already closed.\n");
  if (htmlwr->dest == NULL) die("Destination file not set.\n");
  
  // this may be a secondary invocation, so check if we were replacing with JSON
  if (htmlwr->json_writer) {
    jsonwr_close(htmlwr->json_writer);
    htmlwr->json_writer = NULL;
    htmlwr->variable = NULL;
    fprintf(htmlwr->dest, ";\n    %s", SCRIPT_CLOSE);
  }
  
  while ((which = buf_fread_until(htmlwr->buffer, htmlwr->src, htmlwr->dest, 3, htmlwr->landmarks))) {
    if (which == -1) die("Error reading source\n");
    if (which == 1) {
      before = SCRIPT_START;
      after = SCRIPT_END;
    } else if (which == 2) {
      before = LINK_START;
      after = LINK_END;
    } else if (which == 3) {
      before = IMG_START;
      after = IMG_END;
    } else {
      die("Impossible state");
      return NULL; // stop compiler complaining
    }
    assert(!buf_unexpected(htmlwr->buffer, before, false));
    name = buf_fread_token(htmlwr->buffer, htmlwr->src, is_quote, NULL, false, NULL, 500, NULL);
    bad_name = !regexec_or_die("source name", &(htmlwr->source_re), name, 0, NULL, 0);
    // ensure we have read enough to check the after chars
    buf_compact(htmlwr->buffer);
    while (buf_position(htmlwr->buffer) < strlen(after)) {
      if (buf_fread(htmlwr->buffer, htmlwr->src) == 0) break;
    }
    buf_flip(htmlwr->buffer);
    // check everything matches
    if (bad_name || buf_unexpected(htmlwr->buffer, after, false)) {
      if (bad_name) {
        fprintf(stderr, "Bad file name.\n");
      } else {
        fprintf(stderr, "Could not find trailing letters.\n");
      }
      fputs(before, htmlwr->dest);
      fputs(name, htmlwr->dest);
      free(name);
      continue;
    }
    if (which == 1 && (htmlwr->variable = rbtree_get(htmlwr->replace, name))) {
      fprintf(htmlwr->dest, "%s\n      %s %s\n      var %s = ", 
          SCRIPT_OPEN, JSON_VAR, htmlwr->variable, htmlwr->variable);
      htmlwr->json_writer = jsonwr_open(htmlwr->dest, htmlwr->minimize, 6, 2, 80);
      free(name);
      return htmlwr->json_writer;
    } else if (which == 3) {
      fputs(IMG_START, htmlwr->dest);
      base64enc_file(htmlwr, "image/png", name);
      fputs(IMG_END, htmlwr->dest);
      free(name);
    } else {
      echo_file(htmlwr, which-1, name);
      free(name);
    }
  }
  return NULL;
}
