/********************************************************************
 * FILE: utils.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 9-8-97
 * PROJECT: shared
 * COPYRIGHT: 1997-2013 WSN, TLB
 * DESCRIPTION: Various useful generic utilities.
 ********************************************************************/
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <libgen.h> // for basename
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <ctype.h>
#include "dir.h"
#include "mtwist.h" 	// must come before utils.h
#include "utils.h"

/**********************************************************************
 * See .h file for description.
 **********************************************************************/
#ifdef NOCLOCK
double myclock() {return(0);}
double mysysclock() {return(0);}
double mytime() {return(0);}
#else

#ifdef crayc90
/* No problem on the CRAY. */
#include <time.h>
double myclock() {return((double)clock());}
double mysysclock() {return((double)clock());}
double mytime() {return((double)time(NULL));}

#else
int getrusage(int who, struct rusage *rusage);

double mysysclock()
{
  static bool first_time = true;
  static double    start_time;
  double           elapsed;
  struct rusage    ru;

  if (first_time) {
    getrusage(RUSAGE_SELF, &ru);
    start_time = (ru.ru_stime.tv_sec * 1.0E6) + ru.ru_stime.tv_usec;
    first_time = false;
    return 0;

  } else {
    getrusage(RUSAGE_SELF, &ru);
    elapsed = (ru.ru_stime.tv_sec * 1.0E6) + ru.ru_stime.tv_usec -
      start_time;
    return elapsed;
  }
}

double myclock()
{
  static bool first_time = true;
  static double    start_time;
  double           elapsed;
  struct rusage    ru;

  if (first_time) {
    getrusage(RUSAGE_SELF, &ru);
    start_time = (ru.ru_utime.tv_sec * 1.0E6) + ru.ru_utime.tv_usec;
    first_time = false;
    return 0;

  } else {
    getrusage(RUSAGE_SELF, &ru);
    elapsed = (ru.ru_utime.tv_sec * 1.0E6) + ru.ru_utime.tv_usec -
      start_time;
    return elapsed;
  }
}

//
// Return the elapsed Wall Time in microseconds since
// the first call to mytime() or the previous call.
//
double mytime(
  bool previous 	// return time since previous call
)
{
  static bool first_time = true;
  static struct timeval start, end, last;
  double elapsed = 0;

  if (first_time) {
    first_time = false;
    gettimeofday(&start, NULL);
    end = start;
  } else if (previous) { 
    last = end;
    gettimeofday(&end, NULL);
    elapsed = (end.tv_sec - last.tv_sec) * 1.0E6 + (end.tv_usec - last.tv_usec);
  } else {
    gettimeofday(&end, NULL);
    elapsed = (end.tv_sec - start.tv_sec) * 1.0E6 + (end.tv_usec - start.tv_usec);
  }
  return(elapsed);
}

#endif /* crayc90 */
#endif /* NOCLOCK */

/************************************************************************
 * char* make_path_to_file
 *
 * Concatenates the file name to the directory name
 * Caller is responisble for freeing concatenated string.
 * This should be suitable for an UNIX like OS, including cygwin
 * but would be inappropriate for native Windows.
 *
 * RETURN: a pointer to the concatenated string.
 ************************************************************************/
char* make_path_to_file(const char *directory, const char* filename) {
    size_t directory_length = strlen(directory);
    size_t filename_length = strlen(filename);
    size_t max_path_length = directory_length + filename_length  + 2;
    char *path = mm_malloc(max_path_length);
    strncpy(path, directory, max_path_length);
    // Check if the terminating '/' is present in the directory
    // if not add it to the path.
    if (path[directory_length - 1] != '/') {
      strncat(path, "/", max_path_length - directory_length);
    }
    strncat(path, filename, max_path_length - directory_length - 1);
    return path;
}

/************************************************************************
 * char* concat
 *
 * Returns the concatenated string of two given strings
 *
 * RETURN: a pointer to the concatenated string.
 ************************************************************************/
char* concat_string(const char *string1, const char* string2) {
    size_t string1_len = strlen(string1);
    size_t string2_len = strlen(string2);
    size_t max_length = string1_len + string2_len  + 1;
    char *c = mm_malloc(max_length);
    strncpy(c, string1, max_length);
    strncat(c, string2, max_length - string1_len);
    return c;
}

/************************************************************************
 * See .h file for description.
 ************************************************************************/
bool open_file
  (const char *    filename,            /* Name of the file to be opened. */
   const char *    file_mode,           /* Mode to be passed to fopen. */
   bool allow_stdin,         /* If true, filename "-" is stdin. */
   const char *    file_description,
   const char *    content_description,
   FILE **         afile)               /* Pointer to the open file. */
{
  if (filename == NULL) {
    fprintf(stderr, "Error: No %s filename specified.\n", file_description);
    return(false);
  } else if ((allow_stdin) && (strcmp(filename, "-") == 0)) {
    if (strchr(file_mode, 'r') != NULL) {
      fprintf(stderr, "Reading %s from stdin.\n", content_description);
      *afile = stdin;
    } else if (strchr(file_mode, 'w') != NULL) {
      fprintf(stderr, "Writing %s to stdout.\n", content_description);
      *afile = stdout;
    } else {
      fprintf(stderr, "Sorry, I can't figure out whether to use stdin ");
      fprintf(stderr, "or stdout for %s.\n", content_description);
      return(false);
    }
  } else if ((*afile = fopen(filename, file_mode)) == NULL) {
    fprintf(stderr, "Error opening file %s.\n", filename);
    return(false);
  }
  return(true);
}

/*************************************************************************
 * Add one string to the end of another using the index at to indicate where
 * the end of the string to be added to is.
 *************************************************************************/
static inline void add(char *to, char *from, int max, int *at) {
  assert(to != NULL);
  assert(from != NULL);
  assert(max > 0);
  assert(at != NULL);
  assert(*at >= 0);
  to += (*at);
  max -= 1;
  while (*from != '\0' && *at < max) {
    *to = *from;
    *at += 1;
    ++to;
    ++from;
  }
  *to = '\0';
}

/*************************************************************************
 * Run a program from a given directory with given arguments.  Return
 * the resulting pipe.
 *************************************************************************/
static FILE* run_program
  (char*      program,     // The program to run in the pipe.
   char*      directory,   // Directory where program resides.
   char*      arguments,   // Program arguments.
   char*      type)        // Read ("r") or write ("w").
{
  char* command;
  FILE* return_value;
  int size, pos;

  size = (strlen(directory) + strlen(program) + strlen(arguments) + 3);

  // Allocate space for the command.
  command = (char*)mm_malloc(sizeof(char) * size);
  //Create the path
  pos = 0;
  add(command, directory, size, &pos);
  if (pos > 0 && command[pos-1] != '/') {
    add(command, "/", size, &pos);
  }
  add(command, program, size, &pos);
  if (access(command, F_OK | X_OK) == 0) {
    // Formulate the command
    add(command, " ", size, &pos);
    add(command, arguments, size, &pos);
    // Run the program.
    return_value = popen(command, type);
  } else {
    return_value = NULL;
  }
  myfree(command);
  return(return_value);
}

/*************************************************************************
 * Attempt to run a given program in a given directory with the given
 * arguments, and check that it gives the expected one-line reply.
 *************************************************************************/
static bool try_to_run
  (char*      program,          // The program to run in the pipe.
   char*      directory,        // Directory to look in.
   char*      test_arguments,   // Arguments used when searching for program.
   char*      expected_reply)   // Expected reply from search.
{
  char* reply;
  FILE* pipe;
  bool return_value;


  // Allocate space for the reply.
  reply = (char*)mm_malloc(sizeof(char) * (strlen(expected_reply) + 1));

  // Run the command.
  pipe = run_program(program, directory, test_arguments, "r");

  // Check the pipe.
  if (pipe == NULL) {
    return_value = false;
  } else {

    // Read from the pipe.
    if (fgets(reply, strlen(expected_reply) + 1, pipe) == NULL) {
      return_value = false;
    } else {
      return_value = (strcmp(reply, expected_reply) == 0);
    }

    // Close the pipe.
    if (pclose(pipe) == -1) {
      return_value = false;
    }
  }


  myfree(reply);
  return(return_value);
}

/*************************************************************************
 * Open a read-only pipe using a given command line.
 *************************************************************************/
FILE* open_command_pipe
  (char*     program,          // The program to run in the pipe.
   char*     directory,        // Directory to look in.
   char*     test_arguments,   // Arguments used when searching for program.
   char*     expected_reply,   // Expected reply from search.
   char*     real_arguments,   // Arguments used when running the program.
   bool stdout_on_error,  // If command fails, return STDOUT?
   char*     error_message)    // Error or warning if command fails.
{
  FILE* return_value;

  // Try to run the command with no directory specified.
  if (try_to_run(program, "", test_arguments, expected_reply)) {
    return_value = run_program(program, "", real_arguments, "w");
  }

  // Try to run the program in the specified directory.
  else if (try_to_run(program, directory, test_arguments, expected_reply)) {
    return_value = run_program(program, directory, real_arguments, "w");

  } else {

    // If we failed, print the error message.
    fprintf(stderr, "%s", error_message);
    if (stdout_on_error) {
      return_value = stdout;
    } else {
      exit(1);
    }
  }

  return(return_value);
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
void die
  (char *format,
   ...)
{
  va_list  argp;

  fprintf(stderr, "FATAL: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);

#ifdef DEBUG
/*
  Cause a crash
*/
  char *crash = NULL;
  *crash = 'x';
  abort();
#else
  exit(1);
#endif
}

/**************************************************************************
 * See .h file for description.
 **************************************************************************/
void myassert
  (bool die_on_error,
   bool test,
   char * const    format,
   ...)
{
  va_list  argp;

  if (!test) {

    if (die_on_error) {
      fprintf(stderr, "FATAL: ");
    } else {
      fprintf(stderr, "WARNING: ");
    }

    /* Issue the error message. */
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fflush(stderr);

    if (die_on_error) {
#ifdef DEBUG
      abort();
#else
      exit(1);
#endif
    }
  }
}

/********************************************************************
 * void mm_malloc, mm_calloc, mm_realloc
 *
 * See .h file for descriptions.
 ********************************************************************/
void *mm_malloc
  (size_t size)
{
  void * temp_ptr;

  if (size == 0)
    size++;

  temp_ptr = malloc(size);
  //memset(temp_ptr, '\0', size);		// Set this to see RAM requirements.

  if (temp_ptr == NULL)
    die("Memory exhausted.  Cannot allocate %d bytes.", (int)size);

  return(temp_ptr);
}

void *mm_calloc
  (size_t nelem,
   size_t size)
{
  void * temp_ptr;

  /* Make sure we allocate something. */
  if (size == 0) {
    size = 1;
  }
  if (nelem == 0) {
    nelem = 1;
  }

  temp_ptr = calloc(nelem, size);

  if (temp_ptr == NULL)
    die("Memory exhausted.  Cannot allocate %d bytes.", (int)size);

  return(temp_ptr);
}

void * mm_realloc
  (void * ptr,
   size_t  size)
{
  void * temp_ptr;

  /* Make sure we allocate something. */
  if (size == 0)
    size = 1;
  assert(size > 0);

  /* Some non-ANSI systems complain about reallocating NULL pointers. */
  if (ptr == NULL) {
    temp_ptr = malloc(size);
  } else {
    temp_ptr = realloc(ptr, size);
  }

  if (temp_ptr == NULL)
    die("Memory exhausted.  Cannot allocate %d bytes.", (int)size);

  return(temp_ptr);
}

/* The lookup table. */
#define LOG_PRECISION 1.0e5
static PROB_T log_table[(int) LOG_PRECISION + 2];

/**********************************************************************
 * Set up lookup table for log(x), 0 <= x <= 1.
 **********************************************************************/
void init_log_prob
  (void)
{
  int    i_table;
  PROB_T table_value;

  log_table[0] = LOG_ZERO;
  for (i_table = 1; i_table <= LOG_PRECISION; i_table++) {
    table_value = (double)(i_table / LOG_PRECISION);
    log_table[i_table] = log(table_value);
  }
  log_table[i_table] = 0;  /* For use in iterpolation when x=1 */
}

/**********************************************************************
 * Efficiently find log(x), when 0 < x <= 1.  Doesn't check bounds.
 **********************************************************************/
PROB_T log_prob
  (PROB_T value)
{
  const PROB_T scaled_value = value * LOG_PRECISION;
  const int    log_index = (int)scaled_value;
  const PROB_T decimal_part = scaled_value - log_index;
  const PROB_T lower_value = log_table[log_index];
  const PROB_T upper_value = log_table[log_index+1];
  const PROB_T interpolation = decimal_part * (lower_value - upper_value);

  if (value == 0.0) {
    return(LOG_ZERO);
  }
  return(lower_value + interpolation);
}

/**************************************************************************
 * Return the nearest double smaller than the given double
 **************************************************************************/
double get_next_smaller_double(double x) {
  // IEEE doubles run in lexigraphical order
  // so if we want the next smaller double
  // we just need to cast to integer type and
  // decrement.
  *(long long *) &x = *(long long *) &x - 1;
  return x;
}

/**************************************************************************
 * Return the nearest double larger than the given double
 **************************************************************************/
double get_next_larger_double(double x) {
  // IEEE doubles run in lexigraphical order
  // so if we want the next larger double
  // we just need to cast to integer type and
  // increment.
  *(long long *) &x = *(long long *) &x + 1;
  return x;
}

/**************************************************************************
 * See .h file for description.
 **************************************************************************/
bool is_zero
  (double    value,
   bool log_form)
{
  if ((log_form) && (value < LOG_SMALL)) {
    return(true);
  } else if ((!log_form) && (value == 0.0)) {
    return(true);
  } else {
    return(false);
  }
}

/**************************************************************************
 * See .h file for description.
 **************************************************************************/
bool almost_equal
  (double value1,
   double value2,
   double slop)
{
  if ((value1 - slop > value2) || (value1 + slop < value2)) {
    return(false);
  } else {
    return(true);
  }
}

/*************************************************************************
 * Convert a boolean to and from a "true" or "false" string.
 *************************************************************************/
const char* boolean_to_string
 (bool the_boolean) {
  return the_boolean == true ? "true" : "false";
}

bool boolean_from_string
  (char* true_or_false)
{
  if (strcmp(true_or_false, "true") == 0) {
    return(true);
  } else if (strcmp(true_or_false, "false") == 0) {
    return(false);
  } else {
    die("Invalid input to boolean_from_string (%s)\n", true_or_false);
  }
  return(false); /* Unreachable. */
}

/*************************************************************************
 * Writes a unicode codepoint in UTF-8 encoding to the start of the buffer
 * and return it, optionally record the length of the written code unit.
 * The buffer must be at least 6 bytes long to hold any valid codepoint
 * though it will use the minimum possible. The output is NOT null terminated.
 *************************************************************************/
char* unicode_to_string(uint32_t code, char *buffer, int *code_unit_length) {
  int bytes, i;
  if (code <= 0x7F) { // 1 byte, max 7 bits
    buffer[0] = code;
    if (code_unit_length != NULL) *code_unit_length = 1;
    return buffer;
  } else if (code <= 0x7FF) { // 2 bytes, max 11 bits
    bytes = 2;
  } else if (code <= 0xFFFF) { // 3 bytes, max 16 bits
    bytes = 3;
  } else if (code <= 0x1FFFFF) { // 4 bytes, max 21 bits
    bytes = 4;
  } else if (code <= 0x3FFFFFF) { // 5 bytes, max 26 bits
    bytes = 5;
  } else if (code <= 0x7FFFFFFF) { // 6 bytes, max 31 bits
    bytes = 6;
  } else {
    die("a unicode codepoint can be at maximum 31 bits.");
    return NULL;
  }
  for (i = bytes-1; i > 0; i--) {
    buffer[i] = 0x80 | (code & 0x3F);
    code = code >> 6;
  }
  buffer[0] = ((0xFF << (8 - bytes)) & 0xFF) | code;
  if (code_unit_length != NULL) *code_unit_length = code;
  return buffer;
}

/*************************************************************************
 * Returns the Unicode codepoint at the start of the string assuming UTF-8.
 * This function handles NUL bytes!
 *
 * If something is wrong it will return one of these error codes:
 * -1   If the string begins with a byte 10xxxxxx indicating a middle byte,
 * -2   If the codepoint would go past the end of the string,
 * -3   If the start byte is the illegal value 0xFE, or 0xFF
 * -4   If there are too few middle bytes following the start byte,
 * -5   If the codepoint uses more bytes than needed.
 *************************************************************************/
int32_t unicode_from_string(const char *str, size_t len, int *code_unit_length) {
  int bytes, bytes_after, i;
  int32_t codepoint, min;
  if (code_unit_length != NULL) *code_unit_length = 1;
  if ((str[0] & 0x80) == 0x00) {
    codepoint = str[0];
    return codepoint;
  } else if ((str[0] & 0xC0) == 0x80) {
    return -1;
  } else if ((str[0] & 0xE0) == 0xC0) {
    bytes = 2;
    codepoint = (str[0] & 0x1F) << 6;
  } else if ((str[0] & 0xF0) == 0xE0) {
    bytes = 3;
    codepoint = (str[0] & 0x0F) << 12;
  } else if ((str[0] & 0xF8) == 0xF0) {
    bytes = 4;
    codepoint = (str[0] & 0x07) << 18;
  } else if ((str[0] & 0xFC) == 0xF8) {
    bytes = 5;
    codepoint = (str[0] & 0x03) << 24;
  } else if ((str[0] & 0xFE) == 0xFC) {
    bytes = 6;
    codepoint = (str[0] & 0x01) << 30;
  } else if ((str[0] & 0xFF) == 0xFE || (str[0] & 0xFF) == 0xFF) {
    return -3;
  } else {
    die("Impossible state!");
    return -6; 
  }
  if (code_unit_length != NULL) *code_unit_length = bytes;
  if (bytes > len) return -2;
  for (i = 1, bytes_after = bytes - 2; bytes_after >= 0; i++, bytes_after--) {
    if ((str[i] & 0xC0) != 0x80) return -4;
    codepoint |= ((str[i] & 0x3F) << (6 * bytes_after));
  }
  if (bytes > 2) {
    min = 1 << ((6 * (bytes - 2)) + (8 - bytes));
  } else {
    min = 0x80;
  }
  if (codepoint < min) return -5;
  return codepoint;
}

/*************************************************************************
 * MEME and DREME motifs can have very small E-values which are impossible
 * to represent using a double. By converting to log values we lose
 * precision but the range possible becomes much greater.
 *************************************************************************/
double log10_evalue_from_string(const char *str) {
  const char * EVALUE_RE = "^[+]?([0-9]*\\.?[0-9]+)([eE]([-+]?[0-9]+))?$";
  const char * INF_RE = "^[+]?inf(inity)?$";
  regex_t re_evalue, re_inf;
  regmatch_t matches[4];
  char *buffer;
  int len_m, len_e, i, j, myerrno;
  double m, e, log_ev;

  myerrno = 0;
  regcomp(&re_evalue, EVALUE_RE, REG_EXTENDED);
  regcomp(&re_inf, INF_RE, REG_EXTENDED | REG_ICASE);
  if (regexec(&re_evalue, str, 4, matches, 0) == 0) {
    len_m = matches[1].rm_eo - matches[1].rm_so;
    len_e = matches[3].rm_eo - matches[3].rm_so;
    buffer = mm_malloc(sizeof(char) * (MAX(len_m, len_e) + 1));
    for (i = 0, j = matches[1].rm_so; i < len_m; ++i, ++j) buffer[i] = str[j];
    buffer[i] = '\0';
    errno = 0; m = strtod(buffer, NULL); if (errno) myerrno = errno;
    e = 0;
    if (len_e) {
      for (i = 0, j = matches[3].rm_so; i < len_e; ++i, ++j) buffer[i] = str[j];
      buffer[i] = '\0';
      errno = 0; e = strtod(buffer, NULL); if (errno) myerrno = errno;
    }
    free(buffer);
    log_ev = log10(m) + e;
  } else if (regexec(&re_inf, str, 0, matches, 0) == 0) {
    log_ev = HUGE_VAL;
  } else {
    log_ev = 0; // seems safest
    myerrno = EINVAL;
  }
  regfree(&re_evalue);
  regfree(&re_inf);
  errno = myerrno;
  return log_ev;
}

/*************************************************************************
 * MEME and DREME motifs can have very small E-values which are impossible
 * to represent using a double. Internally we represent such small evalues
 * as log evalues. When outputting we have to convert those log values
 * back into a string representation taking into account the possiblity
 * that the user specified an e-value of zero (which we actually recommended
 * in documentation for quite a while).
 *
 * This function writes the E-value to a string buffer in scientific
 * notation.
 *
 * This returns the number of letters that would be written (not including
 * the nul byte) if the buffer was large enough, thus if it returns a number
 * larger than or equal to size it has been truncated and is not nul
 * terminated.
 *************************************************************************/
int log10_evalue_to_string(double log10_ev, int prec, char *buffer, int size) {
  int size_req;
  double m, e;
  if (log10_ev > -HUGE_VAL && log10_ev < HUGE_VAL) { // normal value
    e = floor(log10_ev);
    m = pow(10.0, log10_ev - e);
    // check that rounding up won't cause a 9.9999 to go to a 10
    if (m + (.5 * pow(10,-prec)) >= 10) {
      m = 1;
      e += 1;
    }
    size_req = snprintf(buffer, size, "%.*fe%+04.0f", prec, m, e);
  } else if (log10_ev >= HUGE_VAL) { // infinity
    strncpy(buffer, "inf", size);
    size_req = 3;
  } else { // negative infinity
    size_req = snprintf(buffer, size, "%.*fe+000", prec, 0.0);
  }
  return size_req;
}

/*************************************************************************
 * MEME and DREME motifs can have very small E-values which are impossible
 * to represent using a double. Internally we represent such small evalues
 * as log evalues. When outputting we have to convert those log values
 * back into a string representation taking into account the possiblity
 * that the user specified an e-value of zero (which we actually recommended
 * in documentation for quite a while).
 *
 * This function writes the E-value to an expandable string buffer in scientific
 * notation. The buffer, expanded to fit the number if required, is returned.
 *************************************************************************/
char* log10_evalue_to_string2(double log10_ev, int prec,
    char **expandable_buffer, int *buffer_size, int offset) {
  char *buffer, *expanded_buffer;
  int space, written, needed;
  assert(*buffer_size >= 0);
  assert(offset >= 0);
  // calculate the size of the buffer
  buffer = (*expandable_buffer)+offset;
  space = *buffer_size - offset;
  // attempt to write the E-value
  written = log10_evalue_to_string(log10_ev, prec, buffer, space);
  // check for success
  if (written >= space) {
    needed = *buffer_size + (written - space) + 1;
    expanded_buffer = realloc(*expandable_buffer, sizeof(char) * needed);
    if (expanded_buffer) {
      *expandable_buffer = expanded_buffer;
      *buffer_size = needed;
    } else {
      //memory allocation failure
      fprintf(stderr, "FATAL: log10_evalue_to_string2 - realloc failed to "
          "expand buffer by %d bytes.\n", (written - space) + 1);
      exit(1);
    }
    // calculate the new size of the buffer
    buffer = (*expandable_buffer)+offset;
    space = *buffer_size - offset;
    written = log10_evalue_to_string(log10_ev, prec, buffer, space);
  }
  return buffer;
}

/**************************************************************************
 * Does a given character appear in a given string?
 **************************************************************************/
bool char_in_string
  (const char* a_string,
   char        a_char)
{
  int  i_string;    /* Index into the string. */
  char string_char; /* Character appearing at that index. */

  i_string = 0;
  string_char = a_string[i_string];
  while (string_char != '\0') {
    if (string_char == a_char) {
      return(true);
    }
    i_string++;
    string_char = a_string[i_string];
  }
  return(false);
}

/**************************************************************************
 * Generic functions for converting between integer and string
 * representations of an enumerated type.
 *
 * Assumes that the longest string representation of the enumerated
 * type does not exceed 100 characters.
 *
 * Assumes that the zeroth enumerated type element is invalid.
 **************************************************************************/
char * convert_enum_type
  (int     enum_type, /* The enumerated type object to be converted. */
   char *  enum_strs[], /* String values associated with this type. */
   int     num_enums) /* Number of values of the type. */
{
  if ((enum_type <= 0) || (enum_type >= num_enums)) {
    die("Illegal enumerated type value (%d).", enum_type);
  }

  return(enum_strs[enum_type]);
}

int convert_enum_type_str
  (char *  enum_type_str, /* String to be converted. */
   int     default_value, /* Value to return if first arg is null. */
   char ** enum_strs,     /* String values associated with this type. */
   int     num_enums)     /* Number of values of the type. */
{
  int i_enum;

  /* If no string was given, return the default. */
  if (enum_type_str == NULL) {
    return(default_value);
  }

  /* Search for the value corresponding to the given string. */
  for (i_enum = 0; i_enum < num_enums; i_enum++) {
    if (strcmp(enum_type_str, enum_strs[i_enum]) == 0) {
      return(i_enum);
    }
  }
  die("Illegal value (%s).", enum_type_str);
  return(0); /* Unreachable. */
}

/****************************************************************************
 * Get the name of the CPU.
 ****************************************************************************/
#define HOST_LENGTH 100
const char* hostname
  ()
{
  FILE *           hostname_stream;
  static char      the_hostname[HOST_LENGTH];
  static bool	   first_time = true;
  int              num_scanned;

  if (first_time) {
    hostname_stream = (FILE *)popen("hostname", "r"); /* SGI needs cast. */
    num_scanned = fscanf(hostname_stream, "%s", the_hostname);
    assert(num_scanned == 1);
    num_scanned = pclose(hostname_stream);
    assert(num_scanned == 0);
  }
  return(the_hostname);
}

/****************************************************************************
 * Get the current date and time.
 ****************************************************************************/
const char* date_and_time
  ()
{
  FILE *           date_stream;
  static char      the_date[HOST_LENGTH];
  static bool first_time = true;

  if (first_time) {
    char *date_str;
    date_stream = (FILE *)popen("date", "r"); /* SGI needs cast. */
    if (date_stream == NULL) die("Running date failed\n");
    date_str = fgets(the_date, HOST_LENGTH, date_stream);
    if (date_str == NULL) die("Read from date via pipe failed\n");
    pclose(date_stream);
    /* Remove the EOL. */
    assert(date_str[strlen(date_str)-1] == '\n');
    date_str[strlen(date_str)-1] = '\0';
    /* Use cached value next time */
    first_time = false;
  }

  return(the_date);
}

/****************************************************************************
 * Get the last modified time for a file.
 * Date buffer must be at least 26 char long.
 ****************************************************************************/
char *get_last_modified_time(char *filename, char *date_buffer) {

  char *date_string = NULL;
  struct stat stbuf; // buffer for stat call
  stat(filename, &stbuf);
  date_string = ctime_r(&stbuf.st_mtime, date_buffer);

  if (date_string) {
    // Remove newline
    int len = strlen(date_string);
    assert(date_string[len - 1] == '\n');
    date_string[len - 1] = '\0';
  }

  return date_string;
}

/****************************************************************************
 * Copy a string, with allocation.
 ****************************************************************************/
void copy_string
 (char** target,
  const char*  source)
{
  if (source == NULL) {
    *target = NULL;
  } else {
    *target = (char *)mm_calloc(strlen(source) + 1, sizeof(char));
    strcpy(*target, source);
  }
}

/************************************************************************
 * Test whether a string is empty or contains only white space
 * Assumes a null terminated string.
 ************************************************************************/
bool is_empty_string
  (char *s)
{
  bool result = true;

  assert(s != NULL);

  while(*s != '\0') {
    if (isspace(*s)) {
      s++;
    } else {
      // Found a non-white space character
      result = false;
      break;
    }
  };

  return result;
}

/************************************************************************
 * Copy an array of integers.
 ************************************************************************/
void copy_int_array
 (int  nelems,
  int* source,
  int* target)
{
  int i;

  for (i = 0; i < nelems; i++)
    target[i] = source[i];
}

/************************************************************************
 * Read characters from a file until a whitespace character is encounterd
 ************************************************************************/
void get_non_blank
(FILE * infile,
 char * a_char)
{
  do {
    *a_char = getc(infile);
    if (*a_char == EOF) {
      die("Premature end of file.\n");
    }
  } while ((*a_char == ' ') || (*a_char == '\n') || (*a_char == '\t'));
}

/************************************************************************
 *  Is the next character in the file the eoln character
 ************************************************************************/
bool is_eoln
  (FILE * infile)
{
  register long ch;

  ch = getc(infile);
  if (ch == EOF) {
    return(true);
  }
  ungetc(ch, infile);
  return(ch == '\n');
}

/************************************************************************
 *  Does a file exist?
 ************************************************************************/
bool file_exists(char* file_path) {

  struct stat stat_buffer;
  bool result = false;

  // Does the file exist?
  if (stat(file_path, &stat_buffer)) {
    if (errno == ENOENT) {
      // stat failed because the path doesn't exist.
      result = false;
    }
    else {
      // stat failed for some other reason
      die(
        "Unable to check for status of file '%s'.\n"
        "Error: %s.\n",
        file_path,
        strerror(errno)
      );
    }
  }
  else {
    result = true;
  }

  return result;
}

/************************************************************************
 *  Is a file executable?
 ************************************************************************/
bool file_executable(char* file_path) {

  struct stat file_status;

  // Is there a path at all?
  if (file_path == NULL) {
    return false;
  }
  
  // Does the file exist?
  if (stat(file_path, &file_status)) {
    if (errno == ENOENT) {
      // stat failed because the path doesn't exist.
      return false;
    } else {
      // stat failed for some other reason
      die(
        "Unable to check for status of file '%s'.\n"
        "Error: %s.\n",
        file_path,
        strerror(errno)
      );
    }
  }
  // stat worked, now check if it's a regular executable
  if (!(S_ISREG(file_status.st_mode) && access(file_path, X_OK) == 0)) return false;
  return true;
} // file_executable

/************************************************************************
 *  Is a file a directory?
 ************************************************************************/
bool file_directory(char* file_path) {

  struct stat file_status;

  // Is there a path at all?
  if (file_path == NULL) {
    return false;
  }
  
  // Does the file exist?
  if (stat(file_path, &file_status)) {
    if (errno == ENOENT) {
      // stat failed because the path doesn't exist.
      return false;
    } else {
      // stat failed for some other reason
      die(
        "Unable to check for status of file '%s'.\n"
        "Error: %s.\n",
        file_path,
        strerror(errno)
      );
    }
  }
  // stat worked, now check if it's a directory
  return (S_ISDIR(file_status.st_mode));
} // file_directory

/******************************************************************
 * This function copies the contents of the file named by the first
 * argument into the file named by the second argument.
 *
 * Returns true if successful, false otherwise.
******************************************************************/
bool copy_file(char *from_filename, char *to_filename) {

  FILE *from = NULL;
  FILE *to = NULL;
  int c = 0;
  int result = false;

  // Open the files
  from = fopen(from_filename, "r");
  if (!from) {
    // Couldn't open the input file.
    fprintf(stderr, "Couldn't open %s for input.\n", from_filename);
    return false;
  }
  to = fopen(to_filename, "w");
  if (!to) {
    // Couldn't open the output file.
    fprintf(stderr, "Couldn't open %s for output.\n", to_filename);
    fclose(from);
    return false;
  }

  // Copy the bytes
  while ((c = fgetc(from)) != EOF) {
    c = fputc(c, to);
    if (c == EOF) {
      break; // Error writing byte
    }
  }

  // Should be at EOF for input with no errors on output.
  if (feof(from) && !ferror(to)) {
    // Success
    result = true;
  }
  else {
    //  Failure
    if (ferror(from)) {
      fprintf(stderr, "Error reading from %s.\n", from_filename);
    }
    if (ferror(to)) {
      fprintf(stderr, "Error writing to %s.\n", to_filename);
    }
    result = false;
  }

  fclose(from);
  fclose(to);
  return result;

}

/******************************************************************
 * This function checks the environment for a variable pointing
 * to where the MEME etc files are installed. If the environment
 * variable is not found, it returns the definitions from the configuration.
******************************************************************/
const char* get_meme_data_dir() {
  const char *data_dir;
  data_dir = getenv("MEME_DATA_DIR");
  return data_dir != NULL ? data_dir : DATA_DIR;
}

/******************************************************************
 * This function checks the environment for a variable pointing
 * to where the MEME executable files are installed. If the environment
 * variable is not found it returns the definition from the configuration.
******************************************************************/
const char* get_meme_bin_dir() {
  const char *bin_dir;
  bin_dir = getenv("MEME_BIN_DIR");
  return bin_dir != NULL ? bin_dir : BIN_DIR;
}

/******************************************************************
 * This function checks the environment for a variable pointing
 * to where the MEME temporary files should be placed. If the environment
 * variable is not found it returns the definition from the configuration
 * and if that is not given then it searches the environment variables
 * TMPDIR, TMP, TEMP and TEMPDIR. If none of those exist it defaults
 * to "/tmp"
******************************************************************/
const char* get_meme_temp_dir() {
  const char *temp_dir;
  // first check the MEME_TEMP_DIR environment variable
  temp_dir = getenv("MEME_TEMP_DIR");
  if (temp_dir != NULL) return temp_dir;
  // second check if a temporary directory was configured
  if (TEMP_DIR[0] != '\0') return TEMP_DIR;
  // finally search other environment variable options, 
  // apparently this is the order checked by the C++ boost::filesystem library
  if ((temp_dir = getenv("TMPDIR")) != NULL) return temp_dir;
  if ((temp_dir = getenv("TMP")) != NULL) return temp_dir;
  if ((temp_dir = getenv("TEMP")) != NULL) return temp_dir;
  if ((temp_dir = getenv("TEMPDIR")) != NULL) return temp_dir;
  // when all else fails default to "/tmp"
  return "/tmp";
}

/******************************************************************
 * This function checks the colon separated dirs for a file with
 * the given name.
 * If the file is found then an allocated path is returned to it,
 * alternatively NULL is returned if no such file exists.
 * The caller is responsible for freeing memory.
******************************************************************/
char* get_meme_dirs_file(const char* dirs, const char* file_name) {
  int start, end, dir_len, sep_len, file_name_len;
  struct stat stat_buffer;
  char *path;
  if (dirs == NULL || file_name == NULL) return NULL;
  file_name_len = strlen(file_name);
  end = 0;
  while (dirs[end] != '\0') {
    start = end;
    for (; dirs[end] != '\0'; end++) {
      if (dirs[end] == ':') break;
    }
    dir_len = end - start;
    sep_len = (dir_len > 0 && dirs[end-1] != '/' ? 1 : 0);
    path = mm_malloc(sizeof(char) * (dir_len + sep_len + file_name_len + 1));
    if (end > start) strncpy(path, dirs+start, dir_len);
    if (sep_len) path[dir_len] = '/';
    strcpy(path+(dir_len + sep_len), file_name);
    path[dir_len + sep_len + file_name_len] = '\0';
    if (stat(path, &stat_buffer) == 0) return path;
    free(path);
    if (dirs[end] == ':') end++;
  }
  return NULL;
}

/******************************************************************
 * This function checks the environment variables MEME_BIN_DIRS and
 * MEME_BIN_DIR for a file with the given name. If neither of those
 * is set then it will look for the file in the compiled BIN_DIR.
 * If the file is found then an allocated path is returned to it,
 * alternatively NULL is returned if no such file exists.
 * The caller is responsible for freeing memory.
******************************************************************/
char* get_meme_bin_file(const char* file_name) {
  struct stat stat_buffer;
  const char *dirs;
  char *path;
  dirs = getenv("MEME_BIN_DIRS");
  if (dirs == NULL) dirs = getenv("MEME_BIN_DIR");
  if (dirs != NULL) {
    return get_meme_dirs_file(dirs, file_name);
  }
  path = make_path_to_file(BIN_DIR, file_name);
  if (stat(path, &stat_buffer) == 0) return path;
  free(path);
  return NULL;
}

/******************************************************************
 * This function checks the environment variables MEME_BIN_DIRS and
 * MEME_LIBEXEC_DIR for a file with the given name. If neither of those
 * is set then it will look for the file in the compiled LIBEXEC_DIR.
 * If the file is found then an allocated path is returned to it,
 * alternatively NULL is returned if no such file exists.
 * The caller is responsible for freeing memory.
******************************************************************/
char* get_meme_libexec_file(const char* file_name) {
  int i;
  struct stat stat_buffer;
  const char *dirs;
  char *path;
  char *dlist[2] = {"MEME_BIN_DIRS", "MEME_LIBEXEC_DIR"};

  for (i=0; i<2; i++) {
    dirs = getenv(dlist[i]);
    if (dirs != NULL) {
      path = get_meme_dirs_file(dirs, file_name);
      if (stat(path, &stat_buffer) == 0) {
        return path;
      } else {
        free(path);
      }
    }
  }

  path = make_path_to_file(LIBEXEC_DIR, file_name);
  if (stat(path, &stat_buffer) == 0) {
    return path;
  } else {
    free(path);
  }

  return NULL;
}

/******************************************************************
 *
 * This function creates a string from the command-line arguments.
 * Caller is responsible for freeing the string returned.
 *
******************************************************************/
char *get_command_line(int argc, char* argv[]) {

  int buffer_size = 200; // Allocate enough for a typical command line
  char *command_line = mm_malloc(sizeof(char) * buffer_size);
  command_line[0] = 0;
  int total_arg_length = 0;
  int i = 0;
  // Work through each of the arguments
  for (i = 0; i < argc; i++) {
    char *nextword = (i == 0) ? basename(argv[0]) : argv[i];
    int arg_length = strlen(nextword) + 2; // +1 for leading space,
                                           // +1 for trailing null
    // Put argument in quotes if it contains white space.
    bool has_whitespace = strchr(nextword, ' ') || strchr(nextword, '\t');
    if (has_whitespace) arg_length += 2;
    total_arg_length += arg_length;
    // Do we have space left in the buffer?
    if (total_arg_length > buffer_size) {
      buffer_size = 2 * total_arg_length;
      command_line =
        mm_realloc(command_line, buffer_size * sizeof(char));
    }
    if (i > 0) {
      // Add leading space if not in the first arugment.
      strcat(command_line, " ");
    }
    if (has_whitespace) strcat(command_line, "'");
    strcat(command_line, nextword);
    if (has_whitespace) strcat(command_line, "'");
  } // i

  return command_line;
}

/******************************************************************
 * Get a description from a file or a string.
 *
 * If there is a file containing the description then read up to
 * DESC_LIMIT characters from the file.
 * If there is a string, then copy it.
 * If neither a file or string is avaliable then return null.
 * Convert windows new lines (CRLF) into into line feed. Convert old mac
 * new lines (CR) into line feed, note that LFCR is not supported. 
 * Merge 3 or more contiguous line feeds into two line feeds.
 * Remove leading and trailing space.
 *
 * The caller is responsible for freeing the returned string.
******************************************************************/
static const size_t DESC_LIMIT = 500;
char *get_job_description(const char* desc_file, const char* desc_message) {
  char *desc;
  int i, shift, found;
  if (desc_file) {
    FILE *fp;
    size_t loaded;
    desc = mm_malloc(sizeof(char) * (DESC_LIMIT + 1));
    fp = fopen(desc_file, "r");
    if (fp == NULL) {
      fprintf(stderr, "Warning: Could not read job description file (%s).\n",
          desc_file);
      return NULL;
    }
    loaded = fread(desc, sizeof(char), DESC_LIMIT, fp);
    fclose(fp);
    desc[loaded] = '\0'; // null terminate
  } else if (desc_message) {
    // duplicate the command-line option so we can edit it
    desc = strdup(desc_message);
  } else {
    return NULL;
  }
  // convert different new-line styles into line-feed
  shift = 0;
  for (i = 0; desc[i] != '\0'; i++) {
    if (desc[i] == '\r') {
      desc[i - shift] = '\n';
      if (desc[i + 1] == '\n') {
        i++;
        shift++;
      }
    } else if (shift > 0) {
      desc[i - shift] = desc[i];
    }
  }
  desc[i - shift] = '\0';
  // remove leading whitespace
  shift = 0;
  for (i = 0; desc[i] != '\0' && isspace(desc[i]); i++) {
    shift++;
  }
  // merge 3 or more new lines into two
  found = 0;
  for (; desc[i] != '\0'; i++) {
    if (desc[i] != '\n') found = 0;
    else found++;
    if (found > 2) {
      shift++;
      continue;
    }
    desc[i - shift] = desc[i];
  }
  desc[i - shift] = '\0';
  // remove trailing whitespace
  for (i = i - shift - 1; i >= 0 && isspace(desc[i]); i--) desc[i] = '\0';
  // resize memory allocation
  desc = mm_realloc(desc, sizeof(char) * (i + 1));
  // return description
  return desc;
}

/* Initialize Mersenne Twist random number generator. */
void srand_mt(
  uint32_t seed
) {
  mt_seed32new(seed);
}

/* Random double in [0.0, 1.0); Mersenne Twist */
/* Uses program-wide state. */
double drand_mt() {
  return mt_ldrand();
}

/* Random integer between 0 (inclusive) and 2^31-1; Mersenne Twist */
/* Uses program-wide state. */
long random_mt() {
  long value;
  // get unsigned value and convert to positive signed by zeroing the top bit
  value = mt_lrand() & 0x7FFFFFFF;
  return value;
}

/* Random integer between 0 (inclusive) and n (exclusive) */
/* n must be > 0 and <= RAND_MAX+1; Mersenne Twist */
/* Uses program-wide state. */
int rand_int(const unsigned n) {  /* from C FAQ */
  assert(n > 0);
  assert(n <= RAND_MAX + 1u);
  unsigned r = mt_ldrand() * n;
  return r;
}

/* Random double between 0 (inclusive) and n (exclusive) */
/* n must be > 0; Mersenne Twist */
/* Uses program-wide state. */
double rand_dbl(const double n) {
  assert(n == n);  /* check for NaN */
  assert(n > 0);
  assert(n <= DBL_MAX);  /* check for inf */
  double r = mt_ldrand() * n;
  return r;
}

/********************************************************************************
Swap to memory regions.
********************************************************************************/
void memswap(void *s, void *t, size_t n) {
  char *a = s;
  char *b = t;

  while (n--) {  /* if n > 0, do-while might be faster? */
    const char tmp = *a;
    *a++ = *b;
    *b++ = tmp;
  }
}

/********************************************************************************
 * Arrange the N elements of ARRAY in random order.
********************************************************************************/
void shuffle(void *base, size_t n, size_t size) {
  for (; n > 1; --n) {
    const int r = rand_int(n);
    memswap((char *)base + (n-1) * size, (char *)base + r * size, size);
  }
}

/*****************************************************************************
 * Outputs the file name when given a path to that file.
 *
 * Searches the path string for '/' chars.  Removes the extension
 * if requested.  Replaces underscores with spaces if requested.
 *
 * This is not intended to be perfect as I'm sure there's some valid way to have
 * a '/' in a file name but I expect it will work most of the time.
 *****************************************************************************/
char *file_name_from_path(
  char *path, 
  bool remove_ext, 
  bool remove_meme_ext, 
  bool replace_underscores
) {

  // scan forwards and find the last '/'
  char *start = path;
  char *where = path;
  for (; *where != '\0'; where++) {
    if (*where == '/') start = where+1;
  }
  char *end = where;		// at null at end of string

  // Remove extension if requested.
  // Work backwards and find the last '.' after the last '/'
  if (remove_ext || remove_meme_ext) {
    for (; where > start; where--) {
      if (*where == '.') {
        if (!remove_meme_ext || strcmp(where, ".meme")==0) {
	  end = where;
        }
	break;
      }
    }
  }
  int len = end - start;

  // Create name.
  char *name = mm_malloc(sizeof(char) * (len + 1));
  char *c;
  int i;
  for (i = 0, c = start; c < end; c++, i++) {
    name[i] = (replace_underscores && *c == '_') ? ' ' : *c;
  }
  name[len] = '\0';

  return name;
} // file_name_from_path

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
