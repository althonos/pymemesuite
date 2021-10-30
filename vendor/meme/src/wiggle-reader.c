/********************************************************************
 * FILE: wiggle-reader.c
 * AUTHOR: Charles Grant, Bill Noble, Tim Bailey, James Johnson
 * CREATION DATE: 2012-05-09
 * COPYRIGHT: 2012 UW
 *
 * This file contains the implmentation for the data structures and
 * functions used to read wiggle files.
 ********************************************************************/

#define _GNU_SOURCE
#include <errno.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"
#include "wiggle-reader.h"

/******************************************************************************
 * This structure is the UDT for reading data from wiggle files.
 *****************************************************************************/
struct wiggle_reader {
  INPUT_FMT_T format;
  FILE *wiggle_file;
  long prev_line_position;
  char* filename;
  char *chrom;
  size_t start;
  size_t step;
  size_t span;
};

#define ERROR_MESSAGE_SIZE 500
static char error_message[ERROR_MESSAGE_SIZE];

#define BUFF_SIZE 100
static char buffer[BUFF_SIZE];

static regex_t fixed_format_regex; // Regex describing fixedStep declaration line
static regex_t fixed_value_regex;  // Regex describing fixedStep data line
static regex_t var_format_regex;   // Regex describing variableStep declaration line
static regex_t var_value_regex;    // Regex describing variableStep data line

#define NUM_MATCHES 6
regmatch_t matches[NUM_MATCHES];

/********************************************************************
 * This function tests whether the matched string is to long to 
 * fit in the shared, fixed length, buffer. If the match would
 * overrun the buffer, we die with an error message.
 ********************************************************************/
static void catch_match_too_long(int len, const char *line, const char *filename) {
  if (len >= BUFF_SIZE) {
    die("Wiggle format line in %s longer than expected:\n%s", filename, line);
  }
}

/********************************************************************
 * This function tests whether the passed in line is a wiggle
 * track line, i.e. it starts with the string 'track'
 ********************************************************************/
static bool is_track_line(const char *line) {

  const char *track_marker = "track";
  const size_t track_marker_size = strlen(track_marker);

  if (strncmp(line, track_marker, track_marker_size) == 0) {
    return true;
  }
  else {
    return false;
  }

}

/********************************************************************
 * This function tests whether the passed in line is a wiggle
 * declaration line, i.e. it starts with the string 'fixedStep'
 * or the string 'variableStep'
 ********************************************************************/
static bool is_format_line(const char *line) {

  const char *fixed_marker = "fixedStep";
  const size_t fixed_marker_size = strlen(fixed_marker);
  const char *var_marker = "variableStep";
  const size_t var_marker_size = strlen(var_marker);

  if (strncmp(line, fixed_marker, fixed_marker_size) == 0) {
    return true;
  }
  else if (strncmp(line, var_marker, var_marker_size) == 0) {
    return true;
  }
  else {
    return false;
  }

}

/********************************************************************
 * This function reads the next line from a WIGGLE_READER_T.
 * If a line was read it is passed out in the line parameter.
 * Any IO problems other than EOF will be treated as a fatal error.
 *
 * Returns true if the line was successfully read, false if EOF
 * is reached.
 ********************************************************************/
static bool read_wiggle_line(WIGGLE_READER_T *reader, char **line) {

  errno = 0;
  *line = getline2(reader->wiggle_file);
  if (*line == NULL) {
    if (feof(reader->wiggle_file)) {
      return false;
    }
    else {
      die(
        "Unable to read wiggle file: %s.\nError message: %s.\n", 
        reader->filename, 
        strerror(errno)
      );
    }
  }

  return true;

}

/********************************************************************
 * This function initializes the regular expressions used to parse
 * the elements of a wiggle file.
 ********************************************************************/
static void setup_wiggle_regexp(WIGGLE_READER_T *reader) {

  int status;

  // Initialize regular express for parsing fixed step format record
  // Expected format is:
  //   fixedStep chrom=<str> start=<integer number> [span=<integer number>]
  status = regcomp(
    &fixed_format_regex, 
    "fixedStep[[:space:]]+"
    "chrom=([^[:space:]:]+)[[:space:]]+"
    "start=([[:digit:]]+)[[:space:]]+"
    "step=([[:digit:]]+)"
    "([[:space:]]+span=)?([[:digit:]]+)?",
    REG_EXTENDED | REG_ICASE | REG_NEWLINE
  );
  if (status != 0) {
    regerror(status, &fixed_format_regex, error_message, ERROR_MESSAGE_SIZE);
    die(
      "Error while intitializing regular expression\n"
      "for parsing fixed step wiggle file: %s\n", 
      error_message
    );
  }

  // Initialize regular express for parsing variable step format record
  // Expected format is:
  //   variableStep chrom=<str> [span=<integer number>]
  status = regcomp(
    &var_format_regex, 
    "variableStep[[:space:]]+"
    "chrom=([^[:space:]:]+)"
    "([[:space:]]+span=)?([[:digit:]]+)?",
    REG_EXTENDED | REG_ICASE | REG_NEWLINE
  );
  if (status != 0) {
    regerror(status, &var_format_regex, error_message, ERROR_MESSAGE_SIZE);
    die(
      "Error while intitializing regular expression\n"
      "for parsing variable step wiggle file: %s\n", 
      error_message
    );
  }

  // Initialize regular express for parsing variable step value record
  // Expected format is:
  //   <number> <floating point number>
  status = regcomp(
    &var_value_regex, 
    "([^[:space:]:]+)[[:space:]]+([^[:space:]:]+)", 
    REG_EXTENDED | REG_ICASE | REG_NEWLINE
  );
  if (status != 0) {
    regerror(status, &var_value_regex, error_message, ERROR_MESSAGE_SIZE);
    die(
      "Error while intitializing regular expression\n"
      "for parsing variable step data in wiggle file: %s\n", 
      error_message
    );
  }

  // Initialize regular express for parsing fixed step value record
  // Expected format is:
  //   <floating point number>
  status = regcomp(
    &fixed_value_regex, 
    "([^[:space:]:]+)", 
    REG_EXTENDED | REG_ICASE | REG_NEWLINE
  );
  if (status != 0) {
    regerror(status, &fixed_value_regex, error_message, ERROR_MESSAGE_SIZE);
    die(
      "Error while intitializing regular expression\n"
      "for parsing fixed step data in wiggle file: %s\n", 
      error_message
    );
  }
}

static void parse_wiggle_format_line(WIGGLE_READER_T *reader, const char *line) {

  reader->format = INVALID_FMT;
  int status = regexec(
    &fixed_format_regex,
    line,
    NUM_MATCHES,
    matches,
    0
  );
  if (!status) {
    reader->format = FIXED_WIGGLE_FMT;
  }
  else if (status != REG_NOMATCH ) {
    regerror(status, &fixed_format_regex, error_message, ERROR_MESSAGE_SIZE);
    die(
      "Error checking for fixed step wiggle file: %s\n"
      "error message is %s\n",
      reader->filename, 
      error_message
    );
  }
  else {
    status = regexec(&var_format_regex, line, NUM_MATCHES, matches, 0);
    if (!status) {
      reader->format = VAR_WIGGLE_FMT;
    }
    else if (status != REG_NOMATCH ) {
      regerror(status, &var_format_regex, error_message, ERROR_MESSAGE_SIZE);
      die(
        "Error checking for variable step wiggle file: %s\n"
        "error message is %s\n",
        reader->filename, 
        error_message
      );
    }
  }

  if (reader->format == FIXED_WIGGLE_FMT) {
    const int CHROM_GROUP = 1;
    int len = matches[CHROM_GROUP].rm_eo - matches[CHROM_GROUP].rm_so;
    catch_match_too_long(len, line, reader->filename);
    myfree(reader->chrom);
    reader->chrom = mm_malloc(len + 1);
    char *last 
      = strncpy(reader->chrom, line + matches[CHROM_GROUP].rm_so, len) + len;
    *last = 0; // Terminate string
    const int START_GROUP = 2;
    len = matches[START_GROUP].rm_eo - matches[START_GROUP].rm_so;
    catch_match_too_long(len, line, reader->filename);
    last = strncpy(buffer, line + matches[START_GROUP].rm_so, len) + len;
    *last = 0; // Terminate string
    reader->start = strtod(buffer, NULL);
    const int STEP_GROUP = 3;
    len = matches[STEP_GROUP].rm_eo - matches[STEP_GROUP].rm_so;
    catch_match_too_long(len, line, reader->filename);
    last = strncpy(buffer, line + matches[STEP_GROUP].rm_so, len) + len;
    *last = 0; // Terminate string
    reader->step = strtod(buffer, NULL);
    // Optional match, may not be present
    const int SPAN_GROUP = 5;
    len = matches[SPAN_GROUP].rm_eo - matches[SPAN_GROUP].rm_so;
    if (len > 0) {
      catch_match_too_long(len, line, reader->filename);
      last = strncpy(buffer, line + matches[SPAN_GROUP].rm_so, len) + len;
      *last = 0; // Terminate string
      reader->span = strtod(buffer, NULL);
    }
  }
  else if (reader->format == VAR_WIGGLE_FMT) {
    const int CHROM_GROUP = 1;
    int len = matches[CHROM_GROUP].rm_eo - matches[CHROM_GROUP].rm_so;
    catch_match_too_long(len, line, reader->filename);
    // Check if chromosome has changed since last call.
    if (reader->chrom == NULL 
        || strncmp(reader->chrom, line + matches[CHROM_GROUP].rm_so, len) != 0) {
      myfree(reader->chrom);
      reader->chrom = mm_malloc(len + 1);
    }
    char *last 
      = strncpy(reader->chrom, line + matches[CHROM_GROUP].rm_so, len) + len;
    *last = 0; // Terminate string
    // Optional match, may not be present
    const int SPAN_GROUP = 3;
    len = matches[SPAN_GROUP].rm_eo - matches[SPAN_GROUP].rm_so;
    if (len > 0) {
      catch_match_too_long(len, line, reader->filename);
      last = strncpy(buffer, line + matches[SPAN_GROUP].rm_so, len) + len;
      *last = 0; // Terminate string
      reader->span = strtod(buffer, NULL);
    }
  }
  else {
    die(
      "Unable to determine type of wiggle format in %s.\n",
      reader->filename
    );
  }
}

static void parse_wiggle_value_line(
  WIGGLE_READER_T *reader,
  char *line,
  size_t *start,
  double *value
) {

  char *end_str = NULL;

  if (reader->format == FIXED_WIGGLE_FMT) {
    // Read single floating point number from line
    errno = 0;
    *value = strtod(line, &end_str);
    if (errno == ERANGE || end_str == line) {
      die("Unable to parse line as floating point number:\n%s", line);
    } 
    *start = reader->start;
    reader->start = reader->start + reader->step;
  }
  else if (reader->format == VAR_WIGGLE_FMT) {
    // Read two floating point numbers from line
    char *value_str = line;
    while (isspace(*value_str)) {
      // Move past leading white space
      ++value_str;
    }
    // Read first number as long
    errno = 0;
    *start = strtol(value_str, &end_str, 10);
    if (errno == ERANGE || errno == EINVAL) {
      die("Unable to parse line as two numbers:\n%s", line);
    }
    while (!isspace(*value_str)) {
      // Move past first number
      ++value_str;
    }
    while (isspace(*value_str)) {
      // Move past white space
      ++value_str;
    }
    // Read second number as double
    errno = 0;
    *value = strtod(value_str, NULL);
    if (errno == ERANGE || end_str == value_str) {
      die("Unable to parse line as two numbers:\n%s", line);
    }
    reader->start = *start;
  }
}

WIGGLE_READER_T *new_wiggle_reader(const char *filename) {

    // Allocate reader
    WIGGLE_READER_T *reader = mm_malloc(sizeof(WIGGLE_READER_T));

    // Open wiggle file
    errno = 0;
    reader->wiggle_file = fopen(filename, "r");
    if (reader->wiggle_file == NULL) {
      die(
        "Unable to open wiggle file: %s.\nError message: %s.\n", 
        filename, 
        strerror(errno)
      );
    }

    // Initialize members
    reader->filename = strdup(filename);
    reader->prev_line_position = 0;
    reader->chrom = NULL;
    reader->start = 0;
    reader->step = 0;
    reader->span = WIGGLE_DEFAULT_SPAN;

    setup_wiggle_regexp(reader);
    
    return reader;
}

void free_wiggle_reader(WIGGLE_READER_T *reader) {
  myfree(reader->filename);
  myfree(reader->chrom);
  myfree(reader);
}

void reset_wiggle_reader(WIGGLE_READER_T *reader) {
  errno = 0;
  rewind(reader->wiggle_file);
  if (errno) {
    die(
      "Error while trying to move to start %s:\n%s\n",
      strerror(errno)
    );
  }
  myfree(reader->chrom);
  reader->chrom = NULL;
  reader->start = 0;
  reader->step = 0;
  reader->span = WIGGLE_DEFAULT_SPAN;
}

/******************************************************************************
 * Read from the current position in the file to the first format line 
 * containing a new sequence.  Update the value of the current sequence.
 *
 * Returns true if it was able to advance to the next sequence, false if 
 * EOF reached before the next sequence was found. Dies if other errors
 * encountered.
 *****************************************************************************/
bool go_to_next_sequence_in_wiggle(
  WIGGLE_READER_T *reader
) {

  char *line = NULL;
  char *prev_chrom = NULL;
  bool result = false;

  while((result = read_wiggle_line(reader, &line)) == true) {
    if (is_format_line(line) == true) {
      // Parse format lines and move to next line
      if (reader->chrom) {
        prev_chrom = strdup(reader->chrom);
      }
      parse_wiggle_format_line(reader, line);
      int chrom_match = 1;
      if (prev_chrom) {
        chrom_match = strcmp(prev_chrom, reader->chrom);
        myfree(prev_chrom);
      }
      if (chrom_match != 0) {
        // We've reached a new sequence
        result = true;
        break;
      }
    }
    myfree(line);
  }

  myfree(line);
  return result;
}

/******************************************************************************
 * Reads the next data line from the wiggle file
 * 
 * Returns true if it was able to read the line, false if it reaches
 * a new sequence sequence or EOF.  Dies if other errors encountered.
 *****************************************************************************/
bool get_next_data_line_from_wiggle(
  WIGGLE_READER_T *reader,
  char **chrom,
  size_t *start,
  size_t *step,
  size_t *span,
  double *value,
  bool *found_format_line
) {

  char *line = NULL;
  bool done_looking = false;
  bool found_data = false;
  *found_format_line = false;

  // Save current reader state
  INPUT_FMT_T prev_format = reader->format;
  long prev_prev_line_position = reader->prev_line_position;
  size_t prev_start = reader->start;
  size_t prev_step = reader->step;
  size_t prev_span = reader->span;
  char *prev_chrom = NULL;
  if (reader->chrom) {
    prev_chrom = strdup(reader->chrom);
  }
  long save_pos = ftell(reader->wiggle_file);

  // Keep trying until we've read a data line or reached a new sequence
  while (done_looking == false) {

    bool read_line = read_wiggle_line(reader, &line);

    if (read_line != true) {
      // We've reached EOF
      found_data = false;
      done_looking = true;
      break;
    }

    // A wiggle file line is either a track line, 
    // a format line, or a data line.

    if (is_track_line(line) == true) {
      // Skip track lines
      continue;
    }
    else if (is_format_line(line) == true) {

      // A format line may indicate a new sequence
      // or simply be setting the start, span or step.
      // If it's a new sequence we should just back up
      // and return false since no more data blocs are
      // available for this sequence.
      *found_format_line = true;
      parse_wiggle_format_line(reader, line);

      // Check to see if we've reached a new sequence 
      int chrom_match = -1;
      if (prev_chrom) {
        chrom_match = strcmp(prev_chrom, reader->chrom);
      }

      if (chrom_match != 0) {
        // We've reached a new sequence. Back up one and return false
        errno = 0;
        int e = fseek(reader->wiggle_file, save_pos, SEEK_SET);
        if (e) {
          die(
            "Error when trying to get data line in wiggle file: %s.\nError message: %s.\n", 
            reader->filename, 
            strerror(errno)
          );
        }
        // Restore reader values
        reader->format = prev_format;
        reader->prev_line_position = prev_prev_line_position;
        reader->start = prev_start;
        reader->step = prev_step;
        reader->span = prev_span;
        myfree(reader->chrom);
        reader->chrom = strdup(prev_chrom);
        found_data = false;
        done_looking = true;
      }
      else {
        // Must be updating the step or span
        if (*chrom) {
          myfree(*chrom);
        }
        *chrom = strdup(reader->chrom);
        *step = reader->step;
        *span = reader->span;
        reader->prev_line_position = save_pos;
      }
    }
    else {
      // We've reached a data line
      parse_wiggle_value_line(reader, line, start, value);
      if (*chrom) {
        myfree(*chrom);
      }
      *chrom = strdup(reader->chrom);
      *step = reader->step;
      *span = reader->span;
      reader->prev_line_position = save_pos;
      found_data = true;
      done_looking = true;
    }

    myfree(line);

  }

  myfree(prev_chrom);
  return found_data;
}

/******************************************************************************
 * Sets the state of the reader to the positon before the last data line was
 * read. Only the immediately previous get_next_data_line() can be undone.
 * 
 * Returns true if it was able to undo the get_next_line(), false otherwise.
 * Dies if other errors encountered.
 *****************************************************************************/
bool unget_data_line_from_wiggle(WIGGLE_READER_T *reader) {

  bool result = false;

  if (reader->prev_line_position) {
    errno = 0;
    int e = fseek(reader->wiggle_file, reader->prev_line_position, SEEK_SET);
    if (e) {
      die(
        "Error when trying to unget data line in prio wiggle file: "
        "%s.\nError message: %s.\n", 
        reader->filename, 
        strerror(errno)
      );
    }
    else {
      reader->prev_line_position = 0;
      result = true;
    }
  }

  return result;
}

int get_wiggle_eof(WIGGLE_READER_T *reader) {
  return feof(reader->wiggle_file);
}

const char* get_wiggle_filename(WIGGLE_READER_T *reader) {
  return reader->filename;
}

INPUT_FMT_T get_wiggle_format(WIGGLE_READER_T *reader) {
  return reader->format;
}

size_t get_wiggle_span(WIGGLE_READER_T *reader) {
  return reader->span;
}

size_t get_wiggle_step(WIGGLE_READER_T *reader) {
  return reader->step;
}

char *get_wiggle_seq_name(WIGGLE_READER_T *reader) {
  return strdup(reader->chrom);
}

