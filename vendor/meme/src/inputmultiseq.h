#ifndef INPUTMULTISEQ_H
#define INPUTMULTISEQ_H
#include "st_multidef.h"

//
//  Read a FASTA file into one or two multiseq objects
//
Sint read_fasta_to_multiseqs(
  Multiseq **train_multiseq_ptr,// the training sequences       // OUT
  Multiseq **test_multiseq_ptr, // the test sequences           // OUT
  double hofract,               // put this fraction of sequences in test_multiseq
  int min_ho_size,              // minimum allowed number of test sequences
  BOOL use_rc,			// ensure results are independent of strand
  BOOL is_rna,			// motifs are RNA, sequences may be DNA so convert
  BOOL allow_ambigs,		// don't convert ambiguous characters to SEPARATOR
  char *filename,               // the FASTA file name to read from
  ALPH_T *alph,                 // MEME-style alphabet
  Uchar *compalph,		// lookup table of letter complements
  Uint minlength,		// minimum allowed sequence length
  Uint maxtotallength           // truncate total length of sequences to this length
);

Sint pos2pospair(
  Multiseq *multiseq,
  PairUint *pos,Uint position
);

void freemultiseq(
  Multiseq *multiseq
);

#endif
