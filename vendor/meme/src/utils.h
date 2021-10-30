/********************************************************************
 * FILE: utils.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 9-8-97
 * PROJECT: shared
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Various useful generic utilities.
 ********************************************************************/
#ifndef UTILS_H
#define UTILS_H
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "macros.h"

typedef int VERBOSE_T;
#define INVALID_VERBOSE 0
#define QUIET_VERBOSE 1
#define NORMAL_VERBOSE 2
#define HIGH_VERBOSE 3
#define HIGHER_VERBOSE 4
#define DUMP_VERBOSE 5

#define BIG HUGE_VAL
#define LITTLE -BIG
#define SMALL_POS 1e-300

#define SHUFFLE(base, n) shuffle(base, n, sizeof *(base))

extern VERBOSE_T verbosity;

/*
 * Support LLVM as a compiler.
 */
#ifdef __llvm__
#define inline
#endif

/********************************************************************
 * double myclock
 *
 * Return number of user microseconds since first call to myclock().
 * This corrects the bug in the system version of clock that causes it
 * to loop after about 36 minutes.
 *
 * (Taken from Tim Bailey's MEME package.)
 ********************************************************************/
double myclock(void);
/********************************************************************
 * double mysysclock
 *
 * Return number of system microseconds since first call to my_sys_clock().
 * This corrects the bug in the system version of clock that causes it
 * to loop after about 36 minutes.
 *
 * (Taken from Tim Bailey's MEME package.)
 ********************************************************************/
double mysysclock(void);

//
// Return the elapsed Wall Time in microseconds since
// the first call to mytime() (or the previous call).
//
double mytime(
  bool previous		// if True, time since previous call, not first call
);

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
char* make_path_to_file(const char *directory, const char* filename);

/************************************************************************
 * char* concat
 *
 * Returns the concatenated string of two given strings
 *
 * RETURN: a pointer to the concatenated string.
 ************************************************************************/
char* concat_string( const char *string1, const char* string2);

/************************************************************************
 * int open_file
 *
 * Open a file gracefully.
 *
 * RETURN: Was the open successful?
 ************************************************************************/
bool open_file
  (const char*     filename,            /* Name of the file to be opened. */
   const char*     file_mode,           /* Mode to be passed to fopen. */
   bool allow_stdin,         /* If true, filename "-" is stdin. */
   const char*     file_description,
   const char*     content_description,
   FILE **         afile);              /* Pointer to the open file. */

/*************************************************************************
 * Open a write-only pipe using a given command line.
 *
 * The program argument is the name of a target program to be
 * executed.  This function first searches for the target program in
 * the current working directory, then in the specified directory.  If
 * the target program is still not found, this function tries to
 * locate it via the operating system's PATH variable, by calling the
 * target program with the user-provided test arguments.  The results
 * of this call are compared to the user-provided "expected reply,"
 * and if they match, then a pipe is opened with the real arguments.
 *
 * If the program is not found, the function prints a message to
 * stderr and then either aborts or returns the stdout stream,
 * depending upon the value of stdout_on_error.
 *
 * The "expected reply" is assumed to be no more than one line long.
 *************************************************************************/
FILE* open_command_pipe
  (char*     program,          // The program to run in the pipe.
   char*     directory,        // Directory to look in.
   char*     test_arguments,   // Arguments used when searching for program.
   char*     expected_reply,   // Expected reply from search.
   char*     real_arguments,   // Arguments used when running the program.
   bool stdout_on_error,  // If command fails, return STDOUT?
   char*     error_message);   // Error or warning if command fails.


/********************************************************************
 * DEBUG_CODE (macro)
 *
 * Allow debugging code to be included or excluded from a compiled
 * program.
 ********************************************************************/
#ifdef DEBUG
#define DEBUG_CODE( debug_value, code_fragment ) \
   { if (debug_value) { code_fragment } }
#else
#define DEBUG_CODE( debug_value, code_fragment )
#endif

/********************************************************************
 * void die()
 *
 * Print an error message and die. The arguments are formatted exactly
 * like arguments to printf().
 *
 * (Taken from Sean Eddy's HMMER package.)
 ********************************************************************/
void die
  (char* format,
   ...);

/**************************************************************************
 * Make an assertion, and print the given message if the assertion fails.
 *
 * If the first parameter is set to true, then die if the assertion
 * doesn't go through.  Otherwise, just issue the warning.
 *
 * On exit, dump core if DEBUG is defined.
 **************************************************************************/
void myassert
  (bool die_on_error,
   bool test,
   char*  const    format,
   ...);

/********************************************************************
 * Allocate dynamic memory. Die gracefully if memory is exhausted.
 ********************************************************************/
void *mm_malloc
  (size_t size);
void *mm_calloc
  (size_t nelem,
   size_t size);
void * mm_realloc
  (void * ptr,
   size_t size);

/***************************************************************************
 * Dynamically create or grow an array;
 * P = pointer, N = new size, T = type
 **************************************************************************/
#define mm_resize(P,N,T) { \
  void *new_P; \
  new_P = (P) ? realloc((malloc_t)(P), (N)*sizeof(T)) : malloc((N)*sizeof(T)); \
  if (!new_P) { \
    fprintf(stderr, "mm_resize(" #P "," #N "," #T ") failed!\n"); \
    exit(1); \
  } \
  (P) = (T *) new_P; \
}

/********************************************************************
 * Math macros.
 ********************************************************************/
/* Note that the following type must be the  same as the MTYPE and ATYPE
   defined in 'matrix.h' and 'array.h'. */
typedef double PROB_T;       // Type definition for probability/frequency.
#define PROB_SCAN " %lf"     // Scanf string for PROB_T.

#define LOG_ZERO  (-1.0E10)  // Zero on the log scale.
#define LOG_SMALL (-0.5E10)  // Threshold below which everything is zero.
#define MM_BITS      (33.2)     // = LOG2(-LOG_ZERO)

#ifndef MIN
#define MIN(a,b)         (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b)         (((a)>(b))?(a):(b))
#endif

/***************************************************************************
 * Find the nearest integer.
 ***************************************************************************/
#define nint(x) ((int)((x) >= 0 ? ((x) + 0.5) : ((x) - 0.5)))

/**************************************************************************
 * Compute the logarithm of x, when 0 <= x <= 1.
 **************************************************************************/
void init_log_prob
  (void);

PROB_T log_prob
  (PROB_T value);

#define log_prob2(x)   (log_prob(x) * 1.44269504)

#define LOG_VALUE(logx) \
( ( (logx) < LOG_SMALL ) ? \
    LOG_ZERO : \
    (logx) \
)

#define EXP2(x)                          \
( ( (x) < LOG_SMALL) ?                   \
  0.0 :                                  \
  (exp((x) * 0.69314718 ))               \
)

/**************************************************************************
 * Compute the logarithm of x.  Returns SMALL_POS if x==0.
 **************************************************************************/
static inline PROB_T my_log_helper(PROB_T x, bool base2) {
  double lx;
  if (x > 0.0) {
    lx = log(x);
    if (base2) {
      return(LOG_VALUE(lx) * 1.44269504);
    } else {
      return(LOG_VALUE(lx));
    }
  } else if (x < 0.0) {
    die("Tried to take the log of a negative value (%g).", x);
  } // else if (x == 0.0)
  return(SMALL_POS);
}

#define my_log(x)   my_log_helper(x, false)
#define my_log2(x)  my_log_helper(x, true)


/**************************************************************************
 * Given the logs (in base 2) of two numbers, return the log of their
 * sum.
 *
 * This function is optimized based upon the following formula:
 *
 *      log(x+y) = log(x) + log(1 + exp(log(y) - log(x)))
 *
 **************************************************************************/

#define LOG_SUM1(logx, logy) \
( \
  ( ( (logx) - (logy) ) > MM_BITS ) ? \
    LOG_VALUE(logx) : \
    (logx) + my_log2( 1 + EXP2((logy) - (logx) ) ) \
)

#define LOG_SUM(logx, logy) \
( \
  ( (logx) > (logy) ) ? \
    LOG_SUM1( (logx), (logy) ) : \
    LOG_SUM1( (logy), (logx) ) \
)

/**************************************************************************
 * Return the nearest double smaller than the given double
 **************************************************************************/
double get_next_smaller_double(double x);

/**************************************************************************
 * Return the nearest double larger than the given double
 **************************************************************************/
double get_next_larger_double(double x);

/**************************************************************************
 * Test for zero on a value that may be either a log or a raw float.
 **************************************************************************/
bool is_zero
  (double    value,
   bool log_form);

/**************************************************************************
 * Test to see if two values are approximately equal.
 **************************************************************************/
bool almost_equal
  (double value1,
   double value2,
   double slop);

/*************************************************************************
 * Convert a boolean to and from a "true" or "false" string.
 *************************************************************************/
const char*  boolean_to_string(bool the_boolean);

bool boolean_from_string(char* true_or_false);

/*************************************************************************
 * Writes a unicode codepoint in UTF-8 encoding to the start of the buffer
 * and return it, optionally record the length of the written code unit.
 * The buffer must be at least 6 bytes long to hold any valid codepoint
 * though it will use the minimum possible. The output is NOT null terminated.
 *************************************************************************/
char* unicode_to_string(uint32_t code, char *buffer, int *code_unit_length);

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
int32_t unicode_from_string(const char *str, size_t len, int *code_unit_length);

/*************************************************************************
 * MEME and DREME motifs can have very small E-values which are impossible
 * to represent using a double. By converting to log values we lose 
 * precision but the range possible becomes much greater.
 *************************************************************************/
double log10_evalue_from_string(const char *str);

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
 * This returns the number of characters printed (excluding the nul byte).
 * If the returned value is larger than or equal to size then it has
 * been truncated and not nul terminated.
 *************************************************************************/
int log10_evalue_to_string(double log10_ev, int prec, char *buffer, int size);

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
    char **expandable_buffer, int *buffer_size, int offset);

/**************************************************************************
 * Does a given character appear in a given string?
 **************************************************************************/
bool char_in_string
  (const char* a_string,
   char        a_char);

/**************************************************************************
 * Generic functions for converting between integer and string
 * representations of an enumerated type.
 *
 * Assumes that the longest string representation of the enumerated
 * type does not exceed 100 characters.
 *
 * Assumes that the zeroth enumerated type element is invalid.
 **************************************************************************/
char*  convert_enum_type
  (int     enum_type,  /* The enumerated type object to be converted. */
   char*  enum_strs[],  /* String values associated with this type. */
   int     num_enums); /* Number of values of the type. */

int convert_enum_type_str
  (char*   enum_type_str, /* String to be converted. */
   int     default_value, /* Value to return if first arg is null. */
   char**  enum_strs,     /* String values associated with this type. */
   int     num_enums);    /* Number of values of the type. */

/**************************************************************************
 * Get the name of the CPU.
 **************************************************************************/
const char* hostname
  ();

/**************************************************************************
 * Get the current date and time.
 **************************************************************************/
const char* date_and_time
  ();

/****************************************************************************
 * Get the last modified time for a file.
 * Date buffer must be at least 26 char long.
 ****************************************************************************/
char *get_last_modified_time(char *filename, char *date_buffer);

/**************************************************************************
 * Copy a string, with allocation.
 **************************************************************************/
void copy_string
 (char**  target,
  const char*   source);

/************************************************************************
 * Test whether a string is empty or contains only white space
 * Assumes a null terminated string.
 ************************************************************************/
bool is_empty_string
  (char *s);

/************************************************************************
 * Copy an array of integers.
 ************************************************************************/
void copy_int_array
 (int  nelems,
  int* source,
  int* target);


/************************************************************************
 * Read characters from a file until a whitespace character is encounterd
 ************************************************************************/
void get_non_blank
(FILE * infile,
 char * a_char);

/************************************************************************
 *  Is the next character in the file the eoln character
 ************************************************************************/
bool is_eoln
  (FILE * infile);

/************************************************************************
 *  Does a file exist?
 ************************************************************************/
bool file_exists(char* file_path);

/************************************************************************
 *  Is a file executable?
 ************************************************************************/
bool file_executable(char* file_path);

/************************************************************************
 *  Is a file a directory?
 ************************************************************************/
bool file_directory(char* file_path);

/******************************************************************
 * This function copies the contents of the file named by the first
 * argument into the file named by the second argument.
 *
 * Returns true if successful, false otherwise.
******************************************************************/
bool copy_file(char *from_filename, char *to_filename);

/******************************************************************
 * This function checks the environment for a variable pointing
 * to where the MEME data files are installed. If the environment
 * variable is not found, it returns the definitions from the configuration.
******************************************************************/
const char* get_meme_data_dir();

/******************************************************************
 * This function checks the environment for a variable pointing
 * to where the MEME executable files are installed. If the environment
 * variable is not found it returns the definition from the configuration.
******************************************************************/
const char* get_meme_bin_dir();

/******************************************************************
 * This function checks the environment for a variable pointing
 * to where the MEME temporary files should be placed. If the environment
 * variable is not found it returns the definition from the configuration
 * and if that is not given then it searches the environment variables
 * TMPDIR, TMP, TEMP and TEMPDIR. If none of those exist it defaults
 * to "/tmp"
******************************************************************/
const char* get_meme_temp_dir();

/******************************************************************
 * This function checks the colon separated dirs for a file with
 * the given name.
 * If the file is found then an allocated path is returned to it,
 * alternatively NULL is returned if no such file exists.
 * The caller is responsible for freeing memory.
******************************************************************/
char* get_meme_dirs_file(const char* dirs, const char* file_name);

/******************************************************************
 * This function checks the environment variables MEME_BIN_DIRS and
 * MEME_BIN_DIR for a file with the given name. If neither of those
 * is set then it will look for the file in the compiled BIN_DIR.
 * If the file is found then an allocated path is returned to it,
 * alternatively NULL is returned if no such file exists.
 * The caller is responsible for freeing memory.
******************************************************************/
char* get_meme_bin_file(const char* file_name);

/******************************************************************
 * This function checks the environment variables MEME_BIN_DIRS and
 * MEME_LIBEXEC_DIR for a file with the given name. If neither of those
 * is set then it will look for the file in the compiled LIBEXEC_DIR.
 * If the file is found then an allocated path is returned to it,
 * alternatively NULL is returned if no such file exists.
 * The caller is responsible for freeing memory.
******************************************************************/
char* get_meme_libexec_file();

/******************************************************************
 *
 * This function creates a string from the command-line arguments.
 * Caller is responsible for freeing the string returned.
 *
******************************************************************/
char *get_command_line(int argc, char* argv[]);

/******************************************************************
 * Get a description from a file or a string.
 *
 * If there is a file containing the description then read up to
 * 500 characters from the file.
 * If there is a string, then copy it.
 * If neither a file or string is avaliable then return null.
 * Convert windows new lines (CRLF) into into line feed. Convert old mac
 * new lines (CR) into line feed, note that LFCR is not supported. 
 * Merge 3 or more contiguous line feeds into two line feeds.
 * Remove leading and trailing space.
 *
 * The caller is responsible for freeing the returned string.
******************************************************************/
char *get_job_description(const char* desc_file, const char* desc_message);

/******************************************************************
 * Randomly shuffle the elements of an array.
******************************************************************/
void shuffle(void *base, size_t n, size_t size);

/* Initialize Mersenne Twist random number generator. */
void srand_mt(
  uint32_t seed
);

/* Random double in [0.0, 1.0); Mersenne Twist */
/* Uses program-wide state. */
double drand_mt();

/* Random integer between 0 (inclusive) and 2^31-1; Mersenne Twist */
long random_mt();

/* Random integer between 0 (inclusive) and n (exclusive) */
/* n must be > 0 and <= RAND_MAX+1; Mersenne Twist */
int rand_int(const unsigned n);

/* Random double between 0 (inclusive) and n (exclusive) */
/* n must be > 0; Mersenne Twist */
double rand_dbl(const double n);

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
);
#endif
