#include "read_sequence.h"

/* local functions */
static long read_sequence_data(
  ALPH_T *alph,
  FILE *data_file,                /* data file of sequences */
  char **sequence,                 /* sequence data */
  char *name                        /* name of sequence */
);
static int read_sequence_de(
  FILE *data_file,                /* data file of sequences */
  char **descriptor                /* sequence descriptor */
);

/**********************************************************************/
/*
        read_sequence

        Read a single sequence from the data file.
        Returns false on EOF or bad sequence data.

        Supports FASTA and SWISSPROT formats.

        Note:
        setup_hash_alphabet must have been called prior to first call to this.
*/
/**********************************************************************/
bool read_sequence(
  ALPH_T *alph,
  FILE *data_file,                /* file containing sequences */
  char **sample_name,                /* unique identifier of sequence */
  char **sample_de,                /* descriptor of sequence */
  char **sequence,                /* protein or DNA letters */
  long *length                        /* length of sequence */
)
{
  int i, c;
  FORMAT_TYPE format = FASTA;                /* assume FASTA format */
  char *name = NULL;

  /* skip anything until first sample name */
  c = ' '; 
  while(c != EOF) { 
    if((c=fgetc(data_file)) == '>') {        /* FASTA format */
      break;
    } else if (c == 'I') {                  /* swiss-prot format? */
      if ((c = fgetc(data_file)) == 'D') {
        format = SWISSPROT;
        c = fgetc(data_file);
        Skip_whi(c, data_file);                /* skip whitespace to name */
        break;
      }
    } 
    Skip_eol(c, data_file);                /* go to end of line */
  }
  if (c==EOF) return false;                /* no more sequences */

  /* get the sample name */
  /* read to first blank/tab/ or end of line/file */
  for (i=0; (c=fgetc(data_file))!=EOF; ) {
    if ((c==' ') && i==0) {        /* skip blanks until name starts */
      continue;
    } else if (c==' ' || c=='\t' || c=='\n' || c=='\r') {
      break;                                /* blank or nl ends name */
    } else {
      if ((i % RCHUNK) == 0) {
        Resize(name, i+RCHUNK, char);
      }
      // change any non-ASCII codes to '_'
      if (c > 127) {
        fprintf(stderr, "FASTA sequence ID contains the non-ASCII character code %d; converted to '_'\n", c);
        c = '_';
      }
      name[i++] = c;                        /* non-blank: add to name */
    }
  }
  Resize(name, i+1, char);
  name[i] = '\0';

  /* read in description */
  *sample_de = NULL;                        /* in case no description found */
  if (format == FASTA) {
    if (c != '\n' && c != '\r' && c != EOF) { 
       Skip_whi(c, data_file);                /* skip whitespace to name */
       (void) read_sequence_de(data_file, sample_de);
    }
  } else if (format == SWISSPROT) {
    Skip_eol(c, data_file);                                /* go to end of line */
    while (c != EOF) {                                        /* read all DE lines */
      if ((c=fgetc(data_file)) == 'D') {                 /* start of DE line? */
        if ((c=fgetc(data_file)) == 'E') {
          c=fgetc(data_file);
          Skip_whi(c, data_file);                        /* skip white space */
          (void) read_sequence_de(data_file, sample_de);/* read description */
          c = '\n';                                        /* at end of line */
        }
      } else if (c == 'S') {                                 /* start of SQ line? */
        if ((c=fgetc(data_file)) == 'Q') {
          Skip_eol(c, data_file);                        /* go to end of line */
          break;
        }
      } else {
        Skip_eol(c, data_file);                                /* go to end of line */
      }
    } /* read all DE lines */
  } /* format */

  /* read in the actual sequence data */
  *length = read_sequence_data(alph, data_file, sequence, name);

  /* sequence had bad data */
  if (*length < 0) {
    myfree(name);
    myfree(*sample_de);
    return false;
  }

  *sample_name = name;
  /* insure that the description string exists */
  if (*sample_de == NULL) {
    Resize(*sample_de, 1, char);
    (*sample_de)[0] = '\0';
  }

  return true;
}

/**********************************************************************/
/*
        read_sequence_data

        Read the sequence data into a dynamic array.
        Converts to upper case.
        Checks for illegal sequence characters.

        Returns the length of the sequence or -1 on error.
*/
/**********************************************************************/
static long read_sequence_data(
  ALPH_T *alph,
  FILE *data_file,                /* data file of sequences */
  char **sequence,                 /* sequence data */
  char *name                        /* name of sequence */
)
{
  long length; 
  int c;
  char *seq = NULL;

  /* 
    read sample sequence 
  */
  for(length=0; (c=fgetc(data_file))!=EOF; ) {
    if (c == '>') {                         /* end of FASTA sequence */
      ungetc(c,data_file); 
      break; 
    } else if (c == '/') {                /* end of SWISSPROT sequence */
      c = fgetc(data_file);
      if (c != '/') {
        fprintf(stderr, "\nError reading SWISSPROT database.\n");
        exit(1);
      }
      while ((c=fgetc(data_file))!=EOF && c != '\n') ;        /* find end of line */
      break;
    } else {
      if (isspace(c)) continue;                        /* remove whitespace from seq */
      if (!alph_is_known(alph, c)) {                        /* illegal character? */
        fprintf(stderr, "\nIllegal character `%c' in sequence %s.\n", c, name);
        fprintf(stderr, "Change alphabet or fix data file.\n");
        if (seq != NULL) myfree(seq);
        return -1;
      }
      /* alocate space for ASCII version of sequence */
      if ((length % RCHUNK) == 0) {
        Resize(seq, length+RCHUNK, char);
      }
      seq[length++] = c;
    } /* read sample sequence */
  }
  /* put NULL at the end of the sequence */
  Resize(seq, length+1, char);
  seq[length] = '\0';

  *sequence = seq;

  return length;
}

/**********************************************************************/
/*
        read_sequence_de

        Read the sequence descriptor into a dynamic array.
*/
/**********************************************************************/
static int read_sequence_de(
  FILE *data_file,                /* data file of sequences */
  char **descriptor                /* sequence descriptor; if not NULL,
                                   string will be appended */
)
{
  static int seen_truncate_warning = 0;
  int length, real_length, c;
  char *de = *descriptor;

  /* find current length of descriptor */
  length = 0;
  real_length = 0;
  if (de != NULL) {                        /* continuation of prev. descriptor */
    length = strlen(de);                /* old length */
    /* alocate space for descriptor string */
    Resize(de, (1+length/RCHUNK)*RCHUNK, char);
    de[length++] = ' ';                        /* put a blank at the end */
  }

  /* read in sample descriptor */
  while ( (c=fgetc(data_file)) != EOF) {
    if (c == '\n') { 
      break; 
    } else if (length <= MAXDELEN) {
      /* alocate space for descriptor string */
      if ((length % RCHUNK) == 0) Resize(de, length+RCHUNK, char);
      // change any non-ASCII codes to '_'
      if (c > 127) {
        fprintf(stderr, "FASTA sequence descriptor contains the non-ASCII character code %d; converted to '_'\n", c);
        c = '_';
      }
      de[length++] = c;                        /* append character to de */
    }
    real_length++;
  }
  /* put a NULL at the end of the de */
  Resize(de, length+1, char);
  de[length] = '\0';

  *descriptor = de;
  if (length != real_length && !seen_truncate_warning) {
    fprintf(stderr, "\nWarning: The sequence description was could not be "
        "completely stored. MEME stored %d characters but the description was "
        "%d characters long. To fix this problem increase the constant "
        "MAXDELEN in the file read_sequence.h and recompile. This warning will "
        "not be repeated for any further description truncations.\n", length, 
        real_length);
    seen_truncate_warning = 1;
  }

  return length;
}
