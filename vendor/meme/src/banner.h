#ifndef BANNER_H
#define BANNER_H

#include "macros.h"

extern void banner(
  char *program, /* name of program */
  char *info,	 /* information on program */
  char *cite,     /* how to cite this program */
  FILE *outfile  /* destination for output */
);


#endif
