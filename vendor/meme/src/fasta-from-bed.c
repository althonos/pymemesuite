#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <stdbool.h>
#include "array-list.h"
#include "hash_table.h"
#include "io.h"
#include "string-list.h"
#include "utils.h"

typedef struct {
  char *name;
  long start_offset;
  long length;
  long line_length;
  long line_length_bytes;
} INDEX_ENTRY;

char* program_name = "fasta-from-bed";
VERBOSE_T verbosity = NORMAL_VERBOSE;

INDEX_ENTRY *find_index_entry(HASH_TABLE entry_hash_table, char *name) {

  HASH_TABLE_ENTRY *hte = hash_lookup_str(name, entry_hash_table);
  INDEX_ENTRY *entry = (INDEX_ENTRY *) hash_get_entry_value(hte);

  return entry;
}

void print_sequence_data(
  FILE *output, 
  char *data, 
  INDEX_ENTRY *entry,
  long start, 
  long end
) {
  long num_lines = start / entry->line_length;
  long position = entry->start_offset 
    + (num_lines * entry->line_length_bytes)
    + (start % entry->line_length);

  // First line may be a partial line
  fwrite(data+position, 1, entry->line_length - (start % entry->line_length), output);
  position += (entry->line_length - (start % entry->line_length)) 
    + (entry->line_length_bytes - entry->line_length);

  // Middle lines
  num_lines = end / entry->line_length;
  long final_position = entry->start_offset + num_lines * entry->line_length_bytes;
  while (position < final_position) {
    fwrite(data+position, 1, entry->line_length, output);
    position += entry->line_length_bytes;
  }

  // Possible parial last line
  fwrite(data+position, 1, end % entry->line_length, output);
  fputs("\n", output);
  
}

int main(int argc, char *argv[]) {

  const char *bed_file_name = argv[1];
  const char *genome_file_name = argv[2];
  const char *index_file_name = argv[3];
  const char *fasta_file_name = argv[4];

  const char* usage = "Usage: fasta_from_bed <BED file> <Genome filename> <index filename> <FASTA output filename>\n";

  if (argc != 5) {
    fprintf(stderr, "%s", usage);
    exit(EXIT_SUCCESS);
  }

  // Open genome file as memory mapped 
  errno = 0; 
  FILE *genome_file  = fopen(genome_file_name, "r"); 
  if (genome_file == NULL) { 
    perror("Error:");
    die("Unable to open the genome FASTA file %s\n", genome_file_name);
  }
  // Get length of file by fseeking to the EOF
  errno = 0;
  int status = fseek(genome_file,  0L, SEEK_END);
  if (status < 0) {
    perror("Error:");
    die("Unable to determine size of the genome FASTA file %s\n", genome_file_name);
  }
  long file_size = ftell(genome_file);
  int genome_fd = fileno(genome_file);
  char *genome_data = mmap((caddr_t) 0, file_size, PROT_READ, MAP_PRIVATE, genome_fd, 0);

  // Open bed file for reading

  errno = 0;
  FILE *bed_file = fopen(bed_file_name, "r");
  if (bed_file == NULL) {
    perror("Error:");
    die("Unable to open the BED file %s\n", bed_file_name);
  }

  // Open index file for reading
  errno = 0;
  FILE *index_file = fopen(index_file_name, "r");
  if (index_file == NULL) {
    perror("Error:");
    die("Unable to open the index file %s\n", index_file_name);
  }

  // Open FASTA file for writing
  errno = 0;
  FILE *fasta_file = fopen(fasta_file_name, "w");
  if (fasta_file == NULL) {
    perror("Error:");
    die("Unable to open the FASTA file %s\n", fasta_file_name);
  }

  // Read the index into memory
  char* line = NULL;
  HASH_TABLE entry_hash_table = hash_create(10000, free);
  while ((line = getline2(index_file)) != NULL) {
    STRING_LIST_T* entry_values  = new_string_list_char_split('\t', line);
    INDEX_ENTRY *entry = calloc(1, sizeof(INDEX_ENTRY));
    entry->name = strdup(get_nth_string(0, entry_values));
    entry->length = strtol(get_nth_string(1, entry_values), NULL, 10);
    entry->start_offset = strtol(get_nth_string(2, entry_values), NULL, 10);
    entry->line_length = strtol(get_nth_string(3, entry_values), NULL, 10);
    entry->line_length_bytes = strtol(get_nth_string(4, entry_values), NULL, 10);
    hash_insert_str_value(entry->name, (void *) entry, entry_hash_table);
    myfree(line);
    myfree(entry_values);
  }

  // Iterate through the BED file
  while ((line = getline2(bed_file)) != NULL) {
    STRING_LIST_T* bed_values  = new_string_list_char_split('\t', line);
    char *chrom = get_nth_string(0, bed_values);
    long chrom_start = strtol(get_nth_string(1, bed_values), NULL, 10);
    long chrom_end = strtol(get_nth_string(2, bed_values), NULL, 10);
    INDEX_ENTRY *entry = find_index_entry(entry_hash_table, chrom);
    if (entry) {
      if (chrom_end < entry->length) {
        fprintf(fasta_file, ">%s:%ld-%ld\n", chrom, chrom_start, chrom_end);
        print_sequence_data(
            fasta_file, 
            genome_data, 
            entry,
            chrom_start, 
            chrom_end
        );
      }
      else {
        fprintf(
          stderr, 
          "Feature (%s:%ld-%ld) beyond length of %s size (%ld bp). Skipping.\n", 
          chrom, chrom_start, chrom_end, 
          entry->name, entry->length
        );
      }
    }
    else {
        fprintf(
          stderr, 
          "Feature (%s:%ld-%ld) not found in genome file %s. Skipping.\n", 
          chrom, chrom_start, chrom_end, genome_file_name
        );
    }
    myfree(line);
    myfree(bed_values);
  }

  hash_destroy(entry_hash_table);
  fclose(fasta_file);
  fclose(genome_file);
  fclose(index_file);

  return 0;
}
