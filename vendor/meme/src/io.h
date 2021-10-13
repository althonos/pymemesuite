#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdbool.h>

extern char *getline2(
  FILE *stream 						/* input stream */
);

int create_output_directory(
  char *output_dirname,	/* Name of the output directory to create */
  bool clobber,	/* Whether or not to overwrite an existing dir */
  bool warn	/* Print warning/informative messages to stderr? */
);

/*************************************************************************
 * Output a log value in scientific notation to the given
 * file pointer, or return a string (that must be freed by caller).
 *************************************************************************/
char *print_log_value(FILE *file, double loge_val, int prec);

#endif
