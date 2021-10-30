/****************************************************************************
 * FILE: fasta-io.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 9/14/98
 * PROJECT: MHMM
 * DESCRIPTION: Read sequences from a FASTA file into memory.
 * COPYRIGHT: 1998-2008, UCSD, UCSC, UW
 ****************************************************************************/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include "utils.h"
#include "alphabet.h"
#include "fasta-io.h"
#include "io.h"
#include "seq-reader-from-fasta.h"
#include "prior-reader-from-psp.h"
#include "seq.h"

/****************************************************************************
 * Find the beginning of the next sequence, and read the sequence name
 * and description. The name is everything after the '>' to the
 * first white space. The description is everything after all white
 * space following the name to the end of the line.
 * After calling *name will point to a buffer containing the name, 
 * and *description will point to a buffer containing the description. 
 * The caller is * responsible for freeing the memory in these buffers.
 *
 * The return value is true if a title line was read and false otherwise.
 ****************************************************************************/
static bool read_title_line
  (FILE* fasta_file,
   char** name,
   char** description)
{
  char *id_line = NULL;  // Line containing the ID and comment.

  // Read until the first occurrence of ">".
  do {
    if (id_line)free(id_line);
    id_line = getline2(fasta_file);
  } while (id_line != NULL && id_line[0] != '>');

  if (id_line == NULL) {
    // We never found a title line
    return false;
  }

  // id_line will contain trailing '\n'
  size_t line_length = strlen(id_line) - 1;
  id_line[line_length] = 0;
  if (line_length == 1) {
    // The sequence must have a name.
    // but the line contained only a '>'
    die("Error reading sequence ID.\n");
  }

  // Separate sequence name from description
  size_t name_start = 1;
  size_t description_start = line_length;
  size_t name_length = 0;
  size_t description_length = 0;
  size_t line_index = 0;

  // Find end of sequence name (ends at first white space)
  for (line_index = 1; line_index < line_length; ++line_index) {
    if (isspace(id_line[line_index])) {
      description_start = line_index;
      break;
    }
    name_length = line_index;
  }

  // Skip over white space at the start of description
  for(; description_start < line_length; ++description_start) {
    if (!isspace(id_line[description_start])) {
      break;
    }
  }

  if (description_start < line_length) {
    description_length = line_length - description_start;
  }

  *name = mm_malloc(sizeof(char) * (name_length + 1));
  strncpy(*name, id_line + 1, name_length);
  (*name)[name_length] = 0;
  if (description_length > 0) {
    *description = mm_malloc(sizeof(char) * (description_length + 1));
    strncpy(*description, id_line + description_start, description_length);
    (*description)[description_length] = 0;
  }
  else {
    *description = mm_malloc(sizeof(char));
    (*description)[0] = '\0';
  }

  myfree(id_line);
  return(true);
}

/****************************************************************************
 * Read raw sequence until a '>' is encountered or too many letters
 * are read.  The new sequence is appended to the end of the given
 * sequence.
 *
 * Return: Was the sequence read completely?
 ****************************************************************************/
static bool read_raw_sequence
  (ALPH_T* alph,
   FILE* fasta_file,   // Input Fasta file.
   char* name,         // Sequence ID (used in error messages).
   int   max_chars,    // Maximum number of allowed characters.
   char* raw_sequence) // Pre-allocated sequence.
{
  // tlb; change a_char to integer so it will compile on SGI
  int a_char;
  int i_seq;
  bool return_value = true;

  // Start at the end of the given sequence.
  i_seq = strlen(raw_sequence);
  assert((int)strlen(raw_sequence) <= max_chars);

  // Read character by character.
  while ((a_char = getc(fasta_file)) != EOF) {

    // Check for the beginning of the next sequence.
    if (a_char == '>') {
      // Put the ">" back onto the stream for the next call to find.
      ungetc(a_char, fasta_file);
      break;
    }

    // Skip non-alphabetic characters.
    if (!isalnum(a_char) && a_char != '-' && a_char != '*' && a_char != '.') {
      if ((a_char != ' ') && (a_char != '\t') && (a_char != '\n') && (a_char != '\r')) {
        fprintf(stderr, "Warning: Skipping character %c in sequence %s.\n",
                a_char, name);
      }
    } else {
      // skip check if unknown alph
      if (alph != NULL && !alph_is_known(alph, a_char)) {
        fprintf(stderr, "Warning: Converting illegal character %c to %c ",
                a_char, alph_wildcard(alph));
        fprintf(stderr, "in sequence %s.\n", name);
        a_char = alph_wildcard(alph);
      }
      raw_sequence[i_seq] = (char)a_char;
      i_seq++;
    }
    if (i_seq >= max_chars) {
      return_value = false;
      break;
    }
  }
  raw_sequence[i_seq] = '\0';

  return(return_value);
}

/****************************************************************************
 * Read one sequence from a file in Fasta format.
 *
 * Return: Was a sequence successfully read?
 ****************************************************************************/
bool read_one_fasta
  (ALPH_T*   alph,
   FILE*     fasta_file,
   int       max_seq_length,
   SEQ_T**   sequence)
{
  char *name = NULL;     // Just the sequence ID.
  char *desc = NULL;     // Just the comment field.
  char* buffer = NULL;         // The sequence, as it's read in.

  // Read the title line.
  if (!read_title_line(fasta_file, &name, &desc)) {
    return(false);
  }
  
  // Create local static storage space.
  if (buffer == NULL) {
    buffer = mm_malloc(sizeof(char) * (max_seq_length+1));
  }

  // Read the sequence.
  buffer[0] = '\0';
  if (!read_raw_sequence(alph, fasta_file, name, max_seq_length, buffer)) {
    //die("Sequence %s is too long (maximum=%d).\n", name, max_seq_length);
    fprintf(stderr, "Skipping sequence %s because it is too long (maximum=%d).\n", name, max_seq_length);
    myfree(buffer);
    myfree(name);
    myfree(desc);
    return(false);
  }
    
  // Allocate the new sequence object.
  *sequence = allocate_seq(name, desc, 0, buffer);
  myfree(name);
  myfree(desc);
  myfree(buffer);

  // Tell the user what's up.
  if (verbosity >= DUMP_VERBOSE) {
    fprintf(stderr, 
	"Read sequence %s of length %d.\n",
	get_seq_name(*sequence), get_seq_length(*sequence));
  }
  return(true);
}

/****************************************************************************
 * Read up to N letters of one sequence from a file in Fasta format.
 * If the given sequence is NULL, start from the beginning of the next
 * sequence.  Otherwise, append to the given sequence.
 *
 * Return: Was a sequence segment successfully read?
 ****************************************************************************/
bool read_one_fasta_segment
  (ALPH_T*   alph,
   FILE*     fasta_file,
   int       max_chars,
   SEQ_T**   sequence)
{
  char *name;     // Just the sequence ID.
  char *desc;     // Just the comment field.
  int offset = 0; // Offset of this sequence.
  char *new_sequence = (char*)mm_calloc(max_chars + 1, sizeof(char));
  bool is_complete;  // Is this a complete sequence?

  // Extract information from the current sequence.
  if (*sequence != NULL) {
    // FIXME CEG This is inefficient. We should really add a 
    // method to allow setting the sequence in SEQ_T objects
    char *old_name = get_seq_name(*sequence);
    int len = strlen(old_name) + 1;
    name = mm_malloc(sizeof(char) * len);
    strcpy(name, old_name);
    char *old_desc = get_seq_description(*sequence);
    len = strlen(old_desc) + 1;
    desc = mm_malloc(sizeof(char) * len);
    strcpy(desc, old_desc);
    assert(get_seq_length(*sequence) <= max_chars);
    strcpy(new_sequence, get_raw_sequence(*sequence));
    offset = get_seq_offset(*sequence);
    free_seq(*sequence);
    *sequence = NULL;
  }
  
  // Read information from the file.
  else {
    if (!read_title_line(fasta_file, &name, &desc)) {
      myfree(new_sequence);
      return(false);
    }
  }

  // Read until we run out of space.
  is_complete = read_raw_sequence(alph, fasta_file, name, max_chars,
				  new_sequence);

  // Store the new sequence.
  *sequence = allocate_seq(name, desc, offset, new_sequence);
  set_complete(is_complete, *sequence);
  myfree(name);
  myfree(desc);
  myfree(new_sequence);

  return(true);
}

/****************************************************************************
 * Read all the sequences from a FASTA file at once.
   Multiple files can be appended by calling this more than once.
 ****************************************************************************/
#define NUM_ALLOC 1000 /* Allocate how many sequences at once? */
void read_many_fastas
  (ALPH_T*    alph,
   FILE*      fasta_file,
   int        max_seq, // Maximum sequence length.
   int*       num_seqs,	// Must be set to zero unless appending to sequences.
   SEQ_T***   sequences)
{
  int i_seq;         /* Index of the current sequence. */
  int num_allocated; /* Number of pointers currently allocated. */

  /* Allocate initial memory. */
  num_allocated = *num_seqs + NUM_ALLOC;
  if (*num_seqs > 0) {
    *sequences = (SEQ_T**)mm_realloc(*sequences, sizeof(SEQ_T*) * num_allocated);
  } else {
    *sequences = (SEQ_T**)mm_malloc(sizeof(SEQ_T*) * num_allocated);
  }

  /* Read the sequences one by one. */
  i_seq = *num_seqs;
  while (read_one_fasta(alph, fasta_file, max_seq, &((*sequences)[i_seq]))) {
    i_seq++;

    /* Allocate more space, if need be. */
    if (i_seq >= num_allocated) {
      num_allocated += NUM_ALLOC;
      *sequences = (SEQ_T**)mm_realloc(*sequences, sizeof(SEQ_T*) * num_allocated);
    }
  }

  /* Record the total number of sequences. */
  *num_seqs = i_seq;

  /* Complain if nothing was read. */
  if (i_seq == 0) {
    die("Failed to read a single sequence from the given FASTA file.\n");
  }

}

/****************************************************************************
 * Read all the sequences from a FASTA file at once, parsing the
 * chromosome and position information from the sequence name, and
 * storing the starting position in the SEQ_T "offset" element.
 ****************************************************************************/
void new_read_many_fastas
  (ALPH_T*    alph,
   char*      fasta_filename,
   int        max_seq, 		// Maximum sequence length.
   int*       num_seqs,
   SEQ_T***   sequences)
{
  int i_seq;         /* Index of the current sequence. */
  int num_allocated; /* Number of pointers currently allocated. */
  char *desc = "";	// doesn't read FASTA descriptions (FIXME)

  /* Allocate initial memory. */
  num_allocated = NUM_ALLOC;
  *sequences = (SEQ_T**)mm_malloc(sizeof(SEQ_T*) * num_allocated);

  // Allocate the DATA_BLOCK_READER
  DATA_BLOCK_READER_T *fasta_reader
    = new_seq_reader_from_fasta(true /*parse genomic coord.*/, alph, fasta_filename);

  /* Read the sequences one by one. */
  i_seq = 0;
  while (
    fasta_reader->go_to_next_sequence(fasta_reader) != false
  ) {
    char *fasta_seq_name = NULL;
    bool fasta_result
      = fasta_reader->get_seq_name(fasta_reader, &fasta_seq_name);
    DATA_BLOCK_T *seq_block = new_sequence_block(max_seq);
    bool more_to_read = 
      fasta_reader->get_next_block(fasta_reader, seq_block);
    int start = get_start_pos_for_data_block(seq_block);
    if (more_to_read) {
      if (start != 1) {
        die("Sequence %s:%d-? is too long (maximum=%d).\n", 
          fasta_seq_name, start, max_seq);
      } else {
        die("Sequence %s is too long (maximum=%d).\n", 
          fasta_seq_name, max_seq);
      }
    }
    char *raw_seq = get_sequence_from_data_block(seq_block);
    // save the start of the sequence (from the chr:start-stop) in ->offset
    (*sequences)[i_seq] = allocate_seq(fasta_seq_name, desc, start, raw_seq);
    free(fasta_seq_name);
    free_data_block(seq_block);

    /* Allocate more space, if need be. */
    i_seq++;
    if (i_seq >= num_allocated) {
      num_allocated += NUM_ALLOC;
      *sequences = (SEQ_T**)mm_realloc(*sequences, 
				      sizeof(SEQ_T*) * num_allocated);
    }
  }
  fasta_reader->close(fasta_reader);
  free_data_block_reader(fasta_reader);

  /* Record the total number of sequences. */
  *num_seqs = i_seq;

  /* Complain if nothing was read. */
  if (i_seq == 0) {
    die("Failed to read a single sequence from the given FASTA file.\n");
  }
}

/****************************************************************************
 * Print a single sequence in FASTA format.
 ****************************************************************************/
#define FASTA_LINE 50
void print_fasta
  (SEQ_T*    sequence,
   FILE*     outfile)  
{
  int   i_seq;
  int   seq_length;
  char* raw_sequence;

  fprintf(outfile, ">%s %s\n", get_seq_name(sequence),
	  get_seq_description(sequence));

  raw_sequence = get_raw_sequence(sequence);
  i_seq = 0;
  seq_length = strlen(raw_sequence);
  while (seq_length - i_seq > FASTA_LINE) {
    fprintf(outfile, "%.*s\n", FASTA_LINE, &(raw_sequence[i_seq]));
    i_seq += FASTA_LINE;
  }
  fprintf(outfile, "%s\n\n", &(raw_sequence[i_seq]));
}

/****************************************************************************
 *  Create a sequence object from the first sequence in a FASTA file.
 ****************************************************************************/
SEQ_T* read_sequence_from_file(ALPH_T* alph, char* filename) {
  
  bool result = false;
  FILE* seq_file = NULL;
  SEQ_T* seq = NULL;

  result = open_file(
    filename, 
    "r", 
    1, 
    "sequence", 
    "sequence", 
    &seq_file
  );
  if (result != true) {
    die("Couldn't open the file %s.\n", filename);
  }
  if (read_one_fasta(alph, seq_file, MAX_SEQ, &seq) == false) {
    die("Error while reading %s as a FASTA file.\n", filename);
  } 
  if (seq_file != NULL) {
    fclose(seq_file);
  }

  return seq;
}

/****************************************************************************
 *  Create a zero order background from the sequences in a FASTA file.
 *  Note that the alphabet must be set in alphabet.h.
 ****************************************************************************/
ARRAY_T* calc_bg_from_file(
  ALPH_T* alph, 
  char* filename, 
  bool symmetrical,		// make background symmetrical
  bool allow_stdin
) {
    SEQ_T *seq;
    BGCALC_T *calc;
    FILE *fasta_fp;
    int status, order;
    bool join;
    double eplison;

    order = 0;
    eplison = 1.0;
    status = open_file(filename, "r", allow_stdin, "FASTA sequences", 
        "sequences", &fasta_fp);
    if (!status) exit(EXIT_FAILURE);

    seq = NULL;
    calc = NULL;
    join = true;
    while (read_one_fasta_segment(alph, fasta_fp, MAX_SEQ, &seq)) {
      calculate_markov_model(alph, order, eplison, join, get_raw_sequence(seq), &calc);
      if (!is_complete(seq)) {
        shift_sequence(get_seq_length(seq), seq);
        join = true;
      } else {
        free_seq(seq);
        seq = NULL;
        join = false;
      }
    }
    fclose(fasta_fp);
    ARRAY_T *bg = calculate_markov_model(alph, order, eplison, join, NULL, &calc);
    if (symmetrical) average_freq_with_complement(alph, bg);
    return(bg);
}

#ifdef MAIN
#include "simple-getopt.h"

VERBOSE_T verbosity = NORMAL_VERBOSE;

/*************************************************************************
 * int main
 *************************************************************************/
int main(int argc, char *argv[])
{
  char*     fasta_filename = NULL;
  FILE*     fasta_file;
  SEQ_T*    sequence = NULL;
  SEQ_T**   sequences;
  int       num_seqs;
  int       i_seq;
  bool seq_reader = false;
  bool many = false;
  bool dna = false;
  int       max_chars = -1;
  ALPH_T*   alph;

  // Define command line options.
  cmdoption const options[] = {
    {"seq-reader", NO_VALUE},
    {"many", NO_VALUE},
    {"dna", NO_VALUE},
    {"blocksize", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE}
  };
  int option_count = 5;
  int option_index = 0;

  // Define the usage message.
  char      usage[1100] = "";
  strcat(usage, "Usage: fasta-io [options] <sequences>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     -seq-reader\n");
  strcat(usage, "     -many\n");
  strcat(usage, "     -dna\n");
  strcat(usage, "     -blocksize <n>\n");
  strcat(usage, "     -verbosity 1|2|3|4|5 (default=2)\n");
  strcat(usage, "\n");

	simple_setopt(argc, argv, option_count, options);

  // Parse the command line.
  while (1) { 
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
		const char *  message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
    	simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "seq-reader") == 0) {
      seq_reader = true;
    } else if (strcmp(option_name, "many") == 0) {
      many = true;
    } else if (strcmp(option_name, "dna") == 0) {
      dna = true;
    } else if (strcmp(option_name, "blocksize") == 0) {
      max_chars = atoi(option_value);
    } else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = (VERBOSE_T)atoi(option_value);
    }
  }

  // Read the single required argument.
  if (option_index + 1 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  fasta_filename = argv[option_index];

  if (dna) {
    alph = alph_dna();
  } else {
    alph = alph_protein();
  }

  // Open the file for reading.
  if (open_file(fasta_filename, "r", true, "FASTA", "sequences", &fasta_file)
      == 0) {
    exit(1);
  }

  if (seq_reader) {

    // Read all the sequences into memory.
    new_read_many_fastas(alph, fasta_filename, 1000000, &num_seqs, &sequences);

    // Print each one and then free it.
    for (i_seq = 0; i_seq < num_seqs; i_seq++) {
      print_fasta(sequences[i_seq], stdout);
      free_seq(sequences[i_seq]);
    }
    myfree(sequences);

  } else if (many) {

    // Read all the sequences into memory.
    num_seqs = 0;
    read_many_fastas(alph, fasta_file, MAX_SEQ, &num_seqs, &sequences);

    // Print each one and then free it.
    for (i_seq = 0; i_seq < num_seqs; i_seq++) {
      print_fasta(sequences[i_seq], stdout);
      free_seq(sequences[i_seq]);
    }
    myfree(sequences);


  } else if (max_chars == -1) {

    // Read each sequence one at a time, print and free it.
    while (read_one_fasta(alph, fasta_file, MAX_SEQ, &sequence)) {
      print_fasta(sequence, stdout);
      free_seq(sequence);
    }
  } else {

    // Read each sequence one block at a time, print and free it.
    while (read_one_fasta_segment(alph, fasta_file, max_chars, &sequence)) {
      print_fasta(sequence, stdout);
      if (is_complete(sequence)) {
	free_seq(sequence);
	sequence = NULL;
      } else {
	shift_sequence(get_seq_length(sequence), sequence);
      }
    }
  }
    
  fclose(fasta_file);
  alph_release(alph);
  return(0);
}
#endif
