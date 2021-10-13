/****************************************************************************
 * FILE: bed.c
 * AUTHOR: Timothy L. Bailey
 * CREATE DATE: 3/28/19
 * COPYRIGHT: 2019 UNR
 ****************************************************************************/
#include "bed.h"

// Free one BED entry.
void free_bed_entry(
  BED_ENTRY_T *entry
) {
  if (entry) {
    myfree(entry->chrom);
    myfree(entry);
  }
} // free_bed_entry

// Compare the content of two BED entry pointers by genomic coordinates.
int compare_bed_coords(
 const void *p1,
 const void *p2
) {
  const BED_ENTRY_T *item1 = *(BED_ENTRY_T **) p1;
  const BED_ENTRY_T *item2 = *(BED_ENTRY_T **) p2;

  // Compare chromosome names.
  int cmp = strcmp(item1->chrom, item2->chrom);
  if (cmp != 0) return(cmp);

  // Compare starts.
  cmp = item1->chrom_start - item2->chrom_start;
  if (cmp != 0) return(cmp);

  // Compare ends.
  cmp = item1->chrom_end - item2->chrom_end;
  if (cmp != 0) return(cmp);

  // Compare input pointers.
  cmp = p1 - p2;
  if (cmp != 0) return(cmp);

  die("Result of compare_bed_coords() was ambiguous. p1 = %X p2 = %X\n", (long)p1, (long)p2);
  return(cmp);
} // compare_bed_coords
