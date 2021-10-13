/********************************************************************
 * FILE: metameme.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 05-19-98
 * PROJECT: MHMM
 * COPYRIGHT: 2001-2008, WSN
 * DESCRIPTION: Global functions for the Meta-MEME toolkit.
 ********************************************************************/
#include <stdio.h>
#include <time.h>
#include "utils.h"
#include "metameme.h"
#include "citation.js.h"

/*************************************************************************
 * Print the program header, including help info, if necessary.
 *************************************************************************/
void write_header
  (char*    program,     // Which program called this function? */
   char*    contents,    // What is in this file?
   char*    description, // Text description.
   char*    meme_file,   // Motif filename.
   char*    model_file,  // HMM filename.
   char*    seq_file,    // Name of the input sequence file.
   FILE*    outfile)     // Stream to write to.
{
  time_t current_time;
  fprintf(outfile, "# The MEME Suite of motif-based sequence analysis tools.\n");
  fprintf(outfile, "# \n");
  fputs("# version " VERSION " (Release date: " ARCHIVE_DATE ")\n", outfile);
  fprintf(outfile, "# \n");
  fprintf(outfile, "# Program: %s\n", program);
  fprintf(outfile, "# File contents: %s\n", contents);
  if (description != NULL) {
    fprintf(outfile, "# Model description: %s\n", description);
  }
  if (meme_file != NULL) {
    fprintf(outfile, "# Motif file: %s\n", meme_file);
  }
  if (model_file != NULL) {
    fprintf(outfile, "# HMM file: %s\n", model_file);
  }
  if (seq_file != NULL) {
    fprintf(outfile, "# Sequence file: %s\n", seq_file);
  }
  current_time = time(0);
  fprintf(outfile, "# Create date: %s", ctime(&current_time));
  fprintf(outfile, "# \n");
  fprintf(outfile, "# For further information on how to interpret these ");
  fprintf(outfile, "results please access %s.\n", SITE_URL);
  fprintf(outfile, "# To get a copy of the MEME Suite software, ");
  fprintf(outfile, "please access %s.\n", SOURCE_URL);
  fprintf(outfile, "# \n");
  fprintf(outfile, MCAST_CITE);
#ifdef BAD
  fprintf(outfile, "# If you use MCAST in your research, please cite\n");
  fprintf(outfile, "# \n");
  fprintf(
    outfile, 
    "#  Timothy Bailey and William Stafford Noble,\n"
    "#  \"Searching for statistically significant regulatory modules\",\n"
    "#  Bioinformatics (Proceedings of the European Conference on "
    "Computational Biology),\n"
    "#  19(Suppl. 2):ii16-ii25, 2003.\n"
  );
#endif
  fprintf(outfile, "# \n");
  fprintf(
    outfile, 
    "# ******************************************"
    "***********************************\n"
  );
}
