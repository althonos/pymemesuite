/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1994, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/
/*
        Print name, version, date, reference.
*/

#include "meme.h"
#include "projrel.h"

void banner(
  char *program,  /* name of program */
  char *info, 	  /* information on program */
  char *cite,	  /* how to cite this program */
  FILE *outfile   /* destination for output */
) 
{

  /* announce the program */
  PSTARS(outfile); 
  fprintf(outfile, "%s\n", program); 
  PSTARS(outfile);
  fprintf(outfile, "MEME version " VERSION " (Release date: " ARCHIVE_DATE ")\n\n%s\n", info);
  PSTARS(outfile);
  fprintf(outfile, "\n\n");

  /* print reference citation */
  PSTARS(outfile); 
  fprintf(outfile, "REFERENCE\n"); 
  PSTARS(outfile);
  fprintf(outfile, "%s", cite);
  PSTARS(outfile);
  fprintf(outfile, "\n\n");
}
