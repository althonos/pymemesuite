#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <stdbool.h>
#include "utils.h"

typedef enum { LF = 0, CR, CRLF, NONE} eol_type;
char eol_char[] = { '\n', '\r', '\n', '\0' };

eol_type find_eol_type(char* fasta_data, long file_size) {
  long position = 0;
  if (fasta_data == NULL || file_size < 0) {
    // Invalid inputs
    return NONE;
  }
  while(position < file_size) {
    switch(fasta_data[position]) {
      case '\n': 
        return LF;
      case '\r': 
        if (position < (file_size - 1) && fasta_data[position + 1] == '\n') {
          return CRLF;
        }
        else {
          return CR;
        }
    }
    position++;
  }
  return NONE;
}

long find_next_eol(char *fasta_data, long current_pos, long file_size, char eol) {
  long position = current_pos;
  if (fasta_data == NULL || current_pos < 0 || file_size <= 0|| current_pos >= file_size) {
    // Invalid inputs
    return -1;
  }
  while (position < file_size) {
    if (fasta_data[position] == eol) {
      break;
    }
    position++;
  }
  return position;
}

bool is_valid_sequence_start(char* fasta_data, long position, long file_size, char eol) {

  if (fasta_data == NULL || position < 0 || file_size <= 0|| position >= file_size) {
    // Invalid inputs
    return false;
  }

  // A valid start of sequence is a '>' at the start of a file, or immediately following an end of line
  bool result;
  if (fasta_data[position] == '>' && (position == 0 || fasta_data[position - 1] == eol)) {
    result = true;
  }
  else {
    result = false;
  }

  return result;
}

long get_first_sequence(char *fasta_data, long file_size, char eol) {

  if (fasta_data == NULL || file_size <= 0) {
    // Invalid inputs
    return -1;
  }
  long position = 0;
  while (position < file_size) {
    if (is_valid_sequence_start(fasta_data, position, file_size, eol)) {
      break;
    }
    else if (isspace(fasta_data[position])) {
      // Skip over white space
      position++;
      continue;
    }
    else {
      // Non-whitespace, non-> character before first sequence
      // Not a valid FASTA file.
      position = -1;
      break;
    }
  }

  return position;
}

long get_next_sequence(char *fasta_data, long position, long file_size, char eol) {

  if (fasta_data == NULL || position < 0 || file_size <= 0|| position >= file_size) {
    // Invalid inputs
    return -1;
  }

  while (position < file_size) {
    if (is_valid_sequence_start(fasta_data, position, file_size, eol)) {
      break;
    }
    position++;
  }

  return position;

}

char *get_sequence_name(char *fasta_data, long seq_start, long file_size, char eol) {

  if (fasta_data == NULL || seq_start < 0 || file_size <= 0|| seq_start >= file_size) {
    // Invalid inputs
    return NULL;
  }

  long position = seq_start + 1;
  char* name = NULL;
  while (position < file_size) {
    // The sequence name is from '>' to first white space
    if (isspace(fasta_data[position])) {
       long name_len = position - (seq_start + 1);
       name = calloc(name_len + 1, 1);
       if (name == NULL) {
         die("Failed to allocate memory for storing sequence name.");
       }
       strncpy(name, fasta_data + seq_start + 1, name_len);
       break;
    }
    position++;
  }
  return name;
}

long get_sequence_first_base(char *fasta_data, long seq_start, long file_size, char eol) {

  if (fasta_data == NULL || seq_start < 0 || file_size <= 0|| seq_start >= file_size) {
    // Invalid inputs
    return -1;
  }

  long position = seq_start;

  // Find end of line
  while (position < file_size) {
    if (fasta_data[position] == eol) {
      break;
    }
    position++;
  }
  if (position < file_size) {
    // Move past the end of the line
    position++;
  }
  else {
    // No end of line found
    position = -1;
  }

  return position;
}

long get_sequence_length(char *fasta_data, long seq_start, long file_size, char eol) {

  if (fasta_data == NULL || seq_start < 0 || file_size <= 0|| seq_start >= file_size) {
    // Invalid inputs
    return -1;
  }

  // Move to first end of line in sequence
  long position = seq_start;
  while (position < file_size) {
    if (fasta_data[position] == eol) {
      break;
    }
   position++;
  }
  if (position >= file_size) {
    // Didn't find first end of line
    return -1;
  }
  // Move past end of first line
  position++;

  // Scan until end of sequence indicated by '>' following and end of line or reaching end of data array
  long length = 0;
  while (position < file_size) {
    if (fasta_data[position] == '>' && fasta_data[position - 1] == eol) {
      break;
    }
    if (!isspace(fasta_data[position])) {
      length++;
    }
    position++;
  }
  
  // We're at the start of the next sequence or just past the end of the data array, so back up 1
  return length;
}

long get_sequence_line_length(char *fasta_data, long seq_start, long file_size, char eol) {

  if (fasta_data == NULL || seq_start < 0 || file_size <= 0|| seq_start >= file_size) {
    // Invalid inputs
    return -1;
  }

  // Move to first end of line in sequence
  long position = seq_start;
  while (position < file_size) {
    if (fasta_data[position] == eol) {
      break;
    }
   position++;
  }
  if (position == file_size) {
    // Didn't find first end of line
    return -1;
  }
  // Move past first end of line
  position++;

  // Scan until second end of line or reaching end of data array
  long first_line_length = 0;
  long line_length = 0;
  while (position < file_size) {
    if (fasta_data[position] == '>') {
      break;
    }
    else if (fasta_data[position] == eol) {
      if (first_line_length == 0) {
        // Set length of first line
        first_line_length = line_length;
      }
      else if (line_length != first_line_length) {
        return -1;
      }
      line_length = 0;
      break;
    }
    else  {
      line_length++;
    }
    position++;
  }
  
  return first_line_length;
}

void create_fasta_index_file(const char* fasta_file_name, const char* index_file_name) { 
  // Open fasta file as memory mapped 
  errno = 0; 
  FILE *fasta_file  = fopen(fasta_file_name, "r"); 
  if (fasta_file == NULL) { 
    perror("Error:");
    die("Unable to open the FASTA file %s\n", fasta_file_name);
  }
  // Get length of file by fseeking to the EOF
  errno = 0;
  int status = fseek(fasta_file,  0L, SEEK_END);
  if (status < 0) {
    perror("Error:");
    die("Unable to determine size of the FASTA file %s\n", fasta_file_name);
  }
  long file_size = ftell(fasta_file);
  int fasta_fd = fileno(fasta_file);
  char *fasta_data = mmap((caddr_t) 0, file_size, PROT_READ, MAP_PRIVATE, fasta_fd, 0);

  // Open index file for writing

  errno = 0;
  FILE *index_file = fopen(index_file_name, "w");
  if (index_file == NULL) {
    perror("Error:");
    die("Unable to open the index file %s\n", index_file_name);
  }

  eol_type eol = find_eol_type(fasta_data, file_size);
  if (eol == NONE) {
    die("%s doesn't seem to be a valid FASTA file. No end of line character found.\n", fasta_file_name);
  }

  long sequence_start = get_first_sequence(fasta_data, file_size, eol_char[eol]);
  char *name = get_sequence_name(fasta_data, sequence_start, file_size, eol_char[eol]);
  long offset = get_sequence_first_base(fasta_data, sequence_start, file_size, eol_char[eol]);
  long length = get_sequence_length(fasta_data, sequence_start, file_size, eol_char[eol]);
  long line_length = get_sequence_line_length(fasta_data, sequence_start, file_size, eol_char[eol]);
  long line_length_bytes = (line_length + 2) ? (line_length + 1) : eol == CRLF;
  fprintf(index_file, "%s\t%ld\t%ld\t%ld\t%ld\n", name, length, offset, line_length, line_length_bytes);
  free(name);

  while (true) {
    sequence_start = get_next_sequence(fasta_data, sequence_start + 1, file_size, eol_char[eol]);
    if (sequence_start < 0) {
      die("Invalid inputs when searching for next sequence in %s", fasta_file_name);
    }
    if (sequence_start == file_size) {
      break;
    }
    name = get_sequence_name(fasta_data, sequence_start, file_size, eol_char[eol]);
    offset = get_sequence_first_base(fasta_data, sequence_start, file_size, eol_char[eol]);
    length = get_sequence_length(fasta_data, sequence_start, file_size, eol_char[eol]);
    line_length = get_sequence_line_length(fasta_data, sequence_start, file_size, eol_char[eol]);
    line_length_bytes = (line_length + 2) ? (line_length + 1) : eol == CRLF;
    fprintf(index_file, "%s\t%ld\t%ld\t%ld\t%ld\n", name, length, offset, line_length, line_length_bytes);
    free(name);
  }

  fclose(fasta_file);
  fclose(index_file);
}

int main(int argc, char *argv[]) {

  const char *fasta_filename = argv[1];
  const char *index_filename = argv[2];

  const char* usage = "Usage: index_fasta_file <FASTA filename> <index filename>\n";

  if (argc != 3) {
    fprintf(stderr, "%s", usage);
    exit(EXIT_SUCCESS);
  }
  create_fasta_index_file(fasta_filename, index_filename);
  return 0;
}
