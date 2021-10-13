/****************************************************************************
 * FILE: gtf-io.h
 * AUTHOR: Timothy L. Bailey
 * CREATE DATE: 3/28/19
 * COPYRIGHT: 2019 UNR
 ****************************************************************************/
#ifndef GTF_IO_H
#define GTF_IO_H
#include "gtf.h"

typedef enum {CAGE, LONGPAP, GENCODE} GTF_TYPE;
#define GTF_TYPE_STRINGS "Cage LongPap GenCode"

//
// Read a GTF file into memory.
// Coordinates are assumed to be 1-based, closed.
// Returns the number of entries or -1 if there was a problem.
//
int read_gtf_file(
  FILE *gtf_file,
  bool convert_coords,          // convert to 0-based, half-open format
  GTF_TYPE type,                // type of GTF file
  GTF_ENTRY_T*** gtf_entries
);
#endif
