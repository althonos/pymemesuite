/****************************************************************************
 * FILE: gtf.c
 * AUTHOR: Timothy L. Bailey
 * CREATE DATE: 3/28/19
 * COPYRIGHT: 2019 UNR
 ****************************************************************************/
#include "gtf.h"

// Free one GTF entry.
void free_gtf_entry(
  GTF_ENTRY_T *entry
) {
  if (entry) {
    myfree(entry->seqname);
    myfree(entry->gene_id);
    myfree(entry->gene_name);
    myfree(entry->transcript_id);
    myfree(entry->transcript_type);
    myfree(entry);
  }
} // free_gtf_entry

// Compare the content of two GTF entry pointers by genomic coordinates.
int compare_gtf_coords(
 const void *p1,
 const void *p2
) {
  const GTF_ENTRY_T *item1 = *(GTF_ENTRY_T **) p1;
  const GTF_ENTRY_T *item2 = *(GTF_ENTRY_T **) p2;

  // Compare chromosome names.
  int cmp = strcmp(item1->seqname, item2->seqname);
  if (cmp != 0) return(cmp);

  // Compare starts.
  cmp = item1->start - item2->start;
  if (cmp != 0) return(cmp);

  // Compare ends.
  cmp = item1->end - item2->end;
  if (cmp != 0) return(cmp);

  // Compare transcript IDs
  cmp = strcmp(item1->transcript_id, item2->transcript_id);
  if (cmp != 0) return(cmp);

  die("Result of compare_gtf_coords() was ambiguous.\n");
  return(cmp);
} // compare_gtf_coords
