#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#include "st_multidef.h"
#include "streme-utils.h"

// Record for use in sorting sequences in alphabetical order.
typedef struct {
  Uchar *id;		// pointer to sequence ID
  Uint idlen; 		// length of sequence ID
  Uchar *seq;		// pointer to sequence data
  Uint seqlen;		// length of sequence data
} Seqptr;
DECLAREARRAYSTRUCT(Seqptr);

//
// Compare two Seqptr objects in increasing alphabetical order.
// Returns:
//   +1 if s1->seq > s2->seq 
//   -1 if s1->seq < s2->seq 
// Breaks ties by seqlen, then sequence name, then address.
//
int compare_seqptrs(
  void *v1,
  void *v2
) {
  const Seqptr *s1 = (const Seqptr *) v1;
  const Seqptr *s2 = (const Seqptr *) v2;
  int len = MIN(s1->seqlen, s2->seqlen);
  int cmp = strncmp((char *) s1->seq, (char *) s2->seq, len);
  if (cmp > 0) {
    return(+1);
  } else if (cmp < 0) {
    return(-1);
  } else if (s1->seqlen > s2->seqlen) {
    return(+1);
  } else if (s1->seqlen < s2->seqlen) {
    return(-1);
  } else if (s1->idlen > s2->idlen) {
    return(+1);
  } else if (s1->idlen < s2->idlen) {
    return(-1);
  } else {
    cmp = strncmp((char *) s1->id, (char *) s2->id, s1->idlen);
    if (cmp > 0) {
      return(+1);
    } else if (cmp < 0) {
      return(-1);
    } else {
      return(s1-s2);
    }
  }
} // compare_seqptrs

void initmultiseq(
  Multiseq *multiseq
) {
  multiseq->startdesc = NULL;
  INITARRAY(&multiseq->markpos, Uint);
  INITARRAY(&multiseq->descspace, Uchar);
  multiseq->sequence = NULL;
  multiseq->rcsequence = NULL;
  multiseq->logcumback = NULL;
  multiseq->numofsequences = 0;
  multiseq->totallength = 0;
} // initmultiseq

void freemultiseq(
  Multiseq *multiseq
) {
  if (DELETEMEMORYMAP(multiseq->startdesc) != 0) FREESPACE_TLB(multiseq->startdesc);
  if (DELETEMEMORYMAP(multiseq->markpos.spaceUint) != 0) FREEARRAY_TLB(&multiseq->markpos,Uint);
  if (DELETEMEMORYMAP(multiseq->descspace.spaceUchar) != 0) FREEARRAY_TLB(&multiseq->descspace,Uchar);
  if (multiseq->originalsequence != NULL &&
     multiseq->originalsequence != multiseq->sequence
  ) {
    if (DELETEMEMORYMAP(multiseq->originalsequence) != 0) FREESPACE_TLB(multiseq->originalsequence);
  }
  free(multiseq->sequence);
  if (multiseq->logcumback) free(multiseq->logcumback);
  FREESPACE_TLB(multiseq->rcsequence);
  if (multiseq->freqs) free(multiseq->freqs);
  if (multiseq->background) free(multiseq->background);
  if (multiseq->lcbp) free(multiseq->lcbp);
  if (multiseq->seqstarts) free(multiseq->seqstarts);
  if (multiseq->seqlengths) free(multiseq->seqlengths);
} // freemultiseq

Sint getseqnum(
  Multiseq *multiseq,
  Uint position
) {
  Uint *recordseps = multiseq->markpos.spaceUint;
  Uint numofrecords = multiseq->numofsequences;
  Uint totalwidth = multiseq->totallength;
  Uint *leftptr, *midptr = NULL, *rightptr, len;
  if (numofrecords == UintConst(1) || position < recordseps[0]) return 0;
  if (position > recordseps[numofrecords-2]) {
    if (position < totalwidth) return numofrecords - 1;
    fprintf(stderr, "ERROR: cannot find position %lu\n",(Showuint) position);
    return -1;
  }
  leftptr = recordseps;
  rightptr = recordseps + numofrecords - 2;
  while (leftptr<=rightptr) {
    len = (Uint) (rightptr-leftptr);
    midptr = leftptr + DIV2(len);
    // TLB; There was a bug here if the position is actually between 2 records.
    if (position == *midptr) return -1;
    if (*midptr < position) {
      if (position < *(midptr+1)) {
        return (Sint) (midptr - recordseps + 1);
      }
      leftptr = midptr + 1;
    } else {
      if (*(midptr-1) < position) {
        return (Sint) (midptr - recordseps);
      }
      rightptr = midptr-1;
    }
  }
  fprintf(stderr, "ERROR: cannot find position %lu\n",(Showuint) position);
  return -1;
} // getseqnum

Sint pos2pospair(
  Multiseq *multiseq,
  PairUint *pos,Uint position
) {
  Sint retcode = getseqnum(multiseq,position);
  if (retcode < 0) return -1;
  pos->uint0 = (Uint) retcode;
  if (pos->uint0 == 0) {
    pos->uint1 = position;
  } else {
    pos->uint1 = position - multiseq->markpos.spaceUint[pos->uint0-1] - 1;
  }
  return 0;
} // pos2pospair

//
// Convert "T" to "U" and "t" to "u" if motifs are RNA
// but sequences are DNA.
//
void inline convert_dna_to_rna(
  Seqptr sptr 			// record describing memory mapped sequence
) {
  Uint i;
  Uchar *seq = sptr.seq;
  Uint len = sptr.seqlen;

  for (i=0; i<len; i++) {
    Uchar a = seq[i];
    if (a == 't') {
      seq[i] = 'u';
    } else if (a == 'T') {
      seq[i] = 'U';
    }
  }
} // convert_dna_to_rna

//
// Convert a sequence in place to its reverse-complement if
// that is alphabetically smaller.
//
void inline convert_to_rc_if_smaller(
  Seqptr sptr,			// record describing memory mapped sequence
  int alen,			// length of alphabet
  Uchar *compalph		// lookup table of complementary letters
) {
  Uint i, j;

  // Determine if RC sequence is smaller than given sequence.
  BOOL rc_less = False;
  Uchar *seq = sptr.seq;
  Uint len = sptr.seqlen;
  for (i=0, j=len-1; i<j; i++, j--) {
    if (seq[i] < compalph[(Uint) seq[j]]) {
      break;
    } else if (seq[i] > compalph[(Uint) seq[j]]) {
      rc_less = True;
      break;
    } 
  }

  // Convert sequence to RC if less. Conversion is done in-place.
  if (rc_less) {
    for (i=0, j=len-1; i<=j; i++, j--) {
      Uchar a = seq[i];
      seq[i] = compalph[(Uint) seq[j]];
      if (i != j) seq[j] = compalph[(Uint) a];
    }
  }

} // convert_to_rc_if_smaller(

//
// Save a legal sequence.
//
inline void save_legal_sequence(
  ArraySeqptr *seqptrs,		// stack of pointers to saved sequences
  Seqptr sptr,			// sequence pointer record to save
  Uint alen,			// alphabet length
  Uchar *compalph,		// lookup table of letter complements
  BOOL use_rc,			// ensure results are independent of strand
  BOOL is_rna,			// motifs are RNA, sequences may be DNA so convert
  Uint minlength		// minimum allowed sequence length
) {
  if (sptr.seqlen >= minlength) {
    if (is_rna) convert_dna_to_rna(sptr);
    if (use_rc) convert_to_rc_if_smaller(sptr, alen, compalph);
    PUSHARRAY(seqptrs, Seqptr, 4096, sptr); 
  } else {
    DEBUG_FMT(QUIET_VERBOSE, "# Warning: Ignoring sequence '>%*.*s' because it is too short (%d < %d).\n",
      sptr.idlen, sptr.idlen, sptr.id, sptr.seqlen, minlength);
  }
} // save_legal_sequence

//
//  Read a FASTA file into one or two multiseq objects
//  Results are independent of the order of the sequences and
//  which strand is provided (for complementable alphabets).
//  If maxtotallength > 0, the total length of all sequences
//  in the multiseq objects will be truncated at <= maxtotallength.
//
Sint read_fasta_to_multiseqs(
  Multiseq **train_multiseq_ptr,// the training sequences	// OUT
  Multiseq **test_multiseq_ptr,	// the test sequences		// OUT
  double hofract,		// put this fraction of sequences in test_multiseq
  int min_ho_size,		// minimum allowed number of test sequences
  BOOL use_rc,			// ensure results are independent of strand
  BOOL is_rna,			// motifs are RNA, sequences may be DNA so convert
  BOOL allow_ambigs,		// don't convert ambiguous characters to SEPARATOR
  char *filename,		// the FASTA file name to read from
  ALPH_T *alph,			// MEME-style alphabet
  Uchar *compalph,		// lookup table of letter complements
  Uint minlength,		// minimum allowed sequence length
  Uint maxtotallength		// truncate total length of sequences to this length
) {
  Uint alen = alph_size_core(alph);     // length of alphabet
  Uchar *input;			// ptr to memory mapped file
  Uint inputlen;		// size of memory mapped file
  Uchar *inputptr,		// points to a suffix of the input
        *newptr,		// points to the transformed
        tmpchar;		// temporary character
  BOOL indesc = False,          // inside description part of sequence
       inid = False; 		// currently in the ID
  BOOL case_insensitive = alph_is_case_insensitive(alph);
  ArraySeqptr seqptrs;		// stack of pointers to sequences in the FASTA file
  INITARRAY(&seqptrs, Seqptr);	// initialize the sequence stack
  Seqptr sptr;			// current sequence pointer record
  Multiseq *train_multiseq=NULL, *test_multiseq=NULL;
  *train_multiseq_ptr = *test_multiseq_ptr = NULL;

  // Memory map the FASTA file.
  input = CREATEMEMORYMAP (filename, True, &inputlen);
  if (input == NULL || inputlen == 0) {
    fprintf(stderr, "ERROR: cannot open file \"%s\" or file \"%s\" is empty\n", filename, filename);
    return -1;
  }

  //
  // Create an index of the sequences. 
  // In-place remove spaces, 
  // convert to uppercase if needed,
  // and convert non-core letters to the SEPARATOR.
  //
  int num_sequences_read = 0;
  for (inputptr=newptr=input; inputptr<input+inputlen; inputptr++) {
    if (indesc) {
      if (inid) {
        if (isspace((Ctypeargumenttype) *inputptr)) {
          inid = False;		// Done with the ID
        } else {
          sptr.idlen++;		// Still reading the ID
        }
      }
      if (*inputptr == '\n') {
        indesc = False;		// Done with description completely.
        newptr = inputptr + 1;
      }
    } else {
      if (*inputptr == FASTASEPARATOR) {
	// Save the previous sequence record on the stack and start the new one.
	if (num_sequences_read > 0) {
          save_legal_sequence(&seqptrs, sptr, alen, compalph, use_rc, is_rna, minlength);
        }
        num_sequences_read++;
	sptr.id = inputptr+1;		// Don't include the ">" in the ID.
	sptr.idlen = 0;
	sptr.seq = NULL;
	sptr.seqlen = 0;
        indesc = True;
        inid = True;
      } else {
        tmpchar = *inputptr;
        if (!isspace ((Ctypeargumenttype) tmpchar)) {
          // Check that letter is in the alphabet.
	  if (!alph_is_known(alph, tmpchar)) {
	    fprintf(stderr, "ERROR: unknown letter '%c' in sequence %d of FASTA file.\n",
              tmpchar, num_sequences_read);
	    return -2;
          }
	  // Convert to uppercase if the alphabet is not case-sensitive.
	  if (case_insensitive) tmpchar = (Uchar) toupper((Ctypeargumenttype) tmpchar);
	  // Convert non-core letters to the SEPARATOR (unless told not to).
	  if (! allow_ambigs && ! alph_is_core(alph, (Ctypeargumenttype) tmpchar)) tmpchar = SEPARATOR;
	  *newptr = tmpchar;
          if (sptr.seqlen == 0) sptr.seq = newptr;
          newptr++;
          sptr.seqlen++;
        } // ! isspace
      } // ! FASTASEPARATOR
    } // ! indesc
  } // inputptr

  // Push the last sequence onto the sequence pointer array.
  if (num_sequences_read > 0) {
    save_legal_sequence(&seqptrs, sptr, alen, compalph, use_rc, is_rna, minlength);
  }

  // Check that there were sequences read.
  Uint nseqs = TOPINDEXARRAY(&seqptrs, Seqptr) + 1;
  if (nseqs == 0) {
    fprintf(stderr, "ERROR: No (valid) sequences in multiple FASTA file.\n");
    return -2;
  }

  // Sort the index alphabetically by content.
  SORTARRAY(&seqptrs, Seqptr, compare_seqptrs);

  // Create a permuted list of integers from 0 to nseqs-1.
  int *permuted_indices = (int *) malloc(nseqs * sizeof(int));
  int i, j;
  for (i=0; i<nseqs; i++) permuted_indices[i] = i;
  for (i=nseqs-1; i>0; i--) {
    j = drand_mt() * (i + 1);
    int tmp = permuted_indices[i]; permuted_indices[i] = permuted_indices[j]; permuted_indices[j] = tmp; // swap
  }

  // Determine the length of the training and test multiseq arrays.
  Uint train_len = 0;	// total length of training sequences
  Uint test_len = 0;	// total length of length sequences
  int ntest = (hofract * nseqs);
  if (hofract > 0 && ntest < min_ho_size) ntest = 0;	// prevent too small test set
  int ntrain = nseqs - ntest;

  // Get the total amount of storage required for the sequences in the multiseq object(s).
  for (i=0; i<nseqs; i++) {
    j = permuted_indices[i];
    sptr = *(PEEKARRAY(&seqptrs, j, Seqptr));
    if (i < ntrain) {
      train_len += sptr.seqlen;
    } else {
      test_len += sptr.seqlen;
    }
  }
  if (train_len == 0) {
    fprintf(stderr, "ERROR: all sequences in the test set are length 0.\n");
    return -3;
  }
  if (ntest > 0 && test_len == 0) {
    fprintf(stderr, "ERROR: all sequences in the test set are length 0.\n");
    return -3;
  }
  train_len += ntrain - 1;
  test_len += ntest ? ntest - 1 : 0;

  // Create and initialize the multiseq object(s).
  train_multiseq = (Multiseq *) calloc(1, sizeof(Multiseq));
  initmultiseq(train_multiseq);
  train_multiseq->originalsequence = NULL;
  train_multiseq->alph = alph;	// MEME-style alphabet
  train_multiseq->freqs = (double *) calloc(alen, sizeof(double));
  train_multiseq->sequence = (Uchar *) malloc(train_len * sizeof(Uchar));
  train_multiseq->numofsequences = 0;
  train_multiseq->startdesc = ALLOCSPACE_TLB(train_multiseq->startdesc, Uint, ntrain+1);
  if (ntest > 0) {
    test_multiseq = (Multiseq *) calloc(1, sizeof(Multiseq));
    initmultiseq(test_multiseq);
    test_multiseq->originalsequence = NULL;
    test_multiseq->alph = alph;	// MEME-style alphabet
    test_multiseq->freqs = (double *) calloc(alen, sizeof(double));
    test_multiseq->sequence = (Uchar *) malloc(test_len * sizeof(Uchar));
    test_multiseq->numofsequences = 0;
    test_multiseq->startdesc = ALLOCSPACE_TLB(test_multiseq->startdesc, Uint, ntest+1);
  } else {
    test_multiseq = NULL;
  }

  //
  // Concatenate the sequences into the multiseq object(s) in shuffled order,
  // putting the first ntrain into the training set.
  // Truncate the total length of the sequences if maxtotallength > 0.
  //
  Uint train_maxlength = (1-hofract) * maxtotallength;
  Uint test_maxlength = hofract * maxtotallength;
  int i_seq;			// sequence index in current dataset
  Uint totallength = 0;		// total length of sequences in multiseq object
  double *freqs = (double *) calloc(alen, sizeof(double));
  Uint nletters = 0;
  for (i=i_seq=0; i<nseqs; i++, i_seq++) {
    // The first ntrain (shuffled order) sequences are training; rest are test.
    Multiseq *multiseq = (i < ntrain) ? train_multiseq : test_multiseq;
    Uint maxlength = (i < ntrain) ? train_maxlength : test_maxlength;
    if (i==0 || i==ntrain) { 	// start train or test sequences?
      newptr = multiseq->sequence;
      i_seq = 0;
      totallength = 0;
    }
    sptr = *(PEEKARRAY(&seqptrs, permuted_indices[i], Seqptr));
    // Skip sequence if it would cause total length to be exceeded.
    if (maxlength > 0 && totallength + sptr.seqlen > maxlength) continue;
    // Update the number of sequences stored and their total length.
    multiseq->numofsequences++;
    totallength += sptr.seqlen;
    // Store the sequence ID.
    multiseq->startdesc[i_seq] = multiseq->descspace.nextfreeUchar;
    for (j=0; j<sptr.idlen; j++) {
      STOREINARRAY(&multiseq->descspace, Uchar, 4096, sptr.id[j]);
    }
    STOREINARRAY(&multiseq->descspace, Uchar, 4096, (Uchar) '\0');
    // Place the separator and store the markpos (except for first sequence).
    if (i != 0 && i != ntrain) {
      STOREINARRAY(&multiseq->markpos, Uint, 128, (Uint) (newptr - multiseq->sequence));
      *newptr++ = SEPARATOR;
    }
    // Copy the sequence into the multiseq object and update the letter counts,
    // ignoring wildcards.
    for (j=0; j<sptr.seqlen; j++) {
      Uchar c = sptr.seq[j];
      *newptr++ = c;
      Uint ci = A2I(c);
      if (ci < alen) {
        freqs[ci]++;
        nletters++;
      }
    }
    // Record the total length of the multiseq object.
    multiseq->totallength = (Uint) (newptr - multiseq->sequence);
    // Mark the end of the last training descriptor.
    if (i==ntrain-1) {
      multiseq->startdesc[ntrain] = multiseq->descspace.nextfreeUchar;
    }
    // Mark the end of the last test descriptor.
    if (ntest && i==nseqs-1) {
      multiseq->startdesc[ntest] = multiseq->descspace.nextfreeUchar;
    }
  } // nseqs

  // Compute the letter frequencies in the sequences and put same values in train/test.
  // Average complementary letter frequencies if required.
  for (i=0; i<alen; i++) {
    train_multiseq->freqs[i] = use_rc ? (freqs[i]+freqs[I2CI(i)]) / (2.0*nletters) : freqs[i] / nletters;
    if (test_multiseq) test_multiseq->freqs[i] = train_multiseq->freqs[i];
  }

  // Free test_multiseq if it has too few sequences (due to totallength constraint).
  if (test_multiseq && test_multiseq->numofsequences < min_ho_size) {
    freemultiseq(test_multiseq);
    test_multiseq = NULL;
  }

  // Free space.
  free(permuted_indices);
  free(freqs);
  FREEARRAY_TLB(&seqptrs, Seqptr);
  DELETEMEMORYMAP(input);

  // Return multiseq objects.
  *train_multiseq_ptr = train_multiseq;
  *test_multiseq_ptr = test_multiseq;
  return 0;
} // read_fasta_to_multiseqs
