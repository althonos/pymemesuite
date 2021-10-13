#ifndef GTF_H
#define GTF_H
#include "macros.h"
#include "utils.h"

typedef struct Gtf_entry_t {
  char *seqname;
  long start;
  long end;
  char strand;
  char *gene_id;
  char *gene_name;
  char *transcript_id;
  char *transcript_type;
  double expr;
} GTF_ENTRY_T;

// Free one GTF entry.
void free_gtf_entry(
  GTF_ENTRY_T *entry
);

// Compare the content of two GTF entry pointers by genomic coordinates.
int compare_gtf_coords(
 const void *p1,
 const void *p2
);

#endif
