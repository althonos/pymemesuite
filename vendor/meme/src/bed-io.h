/****************************************************************************
 * FILE: bed-io.h
 * AUTHOR: Timothy L. Bailey
 * CREATE DATE: 3/28/19
 * COPYRIGHT: 2019 UNR
 ****************************************************************************/
#ifndef BED_IO_H
#define BED_IO_H
#include "bed.h"
//
// Read a BED file into memory.
// Coordinates are assumed to be 0-based, half-open.
// Returns the number of entries or -1 if there was a problem.
//
int read_bed_file(
  FILE *bed_file,
  bool convert_coords,		// convert to 1-based, closed
  BED_TYPE type,                // type of BED file
  BED_ENTRY_T*** bed_entries
);
#endif
