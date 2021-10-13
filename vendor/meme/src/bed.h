#ifndef BED_H
#define BED_H
#include "macros.h"
#include "utils.h"

typedef enum {BED3, BROADPEAK, NARROWPEAK} BED_TYPE;

#define get_bed_chrom(ptr) (ptr)->chrom
#define get_bed_start(ptr) (ptr)->chrom_start
#define get_bed_end(ptr) (ptr)->chrom_end
#define get_bed_signal(ptr) (ptr)->signal_value
typedef struct Bed_entry_t {
  char *chrom;
  long chrom_start;
  long chrom_end;
  double signal_value;
} BED_ENTRY_T;

void free_bed_entry(
  BED_ENTRY_T *entry
);

// Compare the content of two BED entry pointers by genomic coordinates.
int compare_bed_coords(
 const void *p1,
 const void *p2
);

#endif
