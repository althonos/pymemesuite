/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef MULTIDEF_H
#define MULTIDEF_H
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "st_arraydef.h"
// TLB added 31-May-2020.
#include "alphabet.h"
#include "array.h"

// TLB moved here.
#define SEPARATOR UCHAR_MAX

//}

/*
  This file defines the datatype \texttt{Multiseq} which stores information 
  about \(k\)-sequences \(T_{0}\), \(\ldots\), \(T_{k-1}\):
  \begin{enumerate}
  \item
  For each \(i\in[0,k-1]\), \texttt{startdesc[i]} stores the index in 
  \texttt{descspace.spaceUchar}, where a textual description for sequence 
  \(T_{i}\) starts. A description for sequence \(T_{i}\)
  ends with a newline character at index \texttt{startdesc[i+1]-1}.
  The description can e.g.\ be the text following the symbol 
  \texttt{>} in a fasta formatted file. 
  \item
  For each \(i\in[0,k-2]\), \texttt{markpos[i]} is the position of a 
  \emph{separator character} between sequence \(T_{i}\) and \(T_{i+1}\).
  \item
  Let \(i\in[0,k-1]\).
  If \(i=0\), then \(T_{i}\) is stored in the component \texttt{sequence}
  from index \(0\) to index \(\Size{T_{i}}-1\). 
  If \(i>0\), then \(T_{i}\) is stored in the component \texttt{sequence}
  from index \(\texttt{markpos[i-1]+1}\) to index 
  \(\texttt{markpos[i-1]}+1+\Size{T_{i}}\).
  \item
  \texttt{numofsequences} is the number \(k\) of sequences stored.
  \item
  \texttt{totallength} is the total length of the stored sequences 
  including the \(k-1\) separator characters.
  \end{enumerate}
*/

/*
  The following defines the separator symbol for fasta files.
*/

#define FASTASEPARATOR '>'

/*
  For a given multiseq and sequence number, the following macros deliver
  a pointer to the first character of the description, and the 
  length of the description.
*/

#define DESCRIPTIONSTARTDESC(MS,SN)\
        ((MS)->startdesc[SN])

#define DESCRIPTIONPTR(MS,SN)\
        ((MS)->descspace.spaceUchar + DESCRIPTIONSTARTDESC(MS,SN))

//TLB removed
/* #define DESCRIPTIONLENGTH(MS,SN)\
        (DESCRIPTIONSTARTDESC(MS,(SN)+1) - DESCRIPTIONSTARTDESC(MS,SN))
*/

/*
  The following macros specifies a default initialization of a structure
  of type \texttt{Showdescinfo}.
*/

#define ASSIGNDEFAULTSHOWDESC(DESC)\
            (DESC)->defined = True;\
            (DESC)->skipprefix = 0;\
            (DESC)->maxlength = 0;\
            (DESC)->replaceblanks = False;\
            (DESC)->untilfirstblank = False

/*
  The following defines the undefined file separator position
*/

#define UNDEFFILESEP 0

typedef struct 
{
  ArrayPosition markpos;	       // positions of SEPARATORS
  Uint *startdesc,                     // of length numofsequences + 1
       numofsequences,                 // the number of sequences
       totallength;                    // the total length of all sequences
  ArrayCharacters descspace;           // the space for the descriptions
  Uchar *sequence,                     // the concatenated sequences
        *rcsequence,                   // NULL or points to 
                                       // reverse complemented sequences
        *originalsequence;             // NULL or points to orig. sequence
  // TLB added:
  ALPH_T *alph;			       // MEME-style alphabet
  Uint npos;			       // Number of positive sequences
  Uint nneg;			       // Number of negative sequences
  Uint pos_length;		       // Total length of positive sequences.
  Uint neg_length;		       // Total length of negative sequences.
  BOOL do_rc;			       // Reverse complements have been appended 
  double *freqs;		       // letter frequencies in sequences
  double *background;		       // 0-order background model of sequences (may be from --bfile)
  double *lcbp;			       // (kmer-1)-order log_2 background conditional probabilites Pr(a | w)
  double *logcumback;		       // log cumulative background probabilies of sequences (and rc)
  int bg_order;			       // Markov order of background model
  Uint *seqstarts;		       // start position of each sequence
  Uint *seqlengths;		       // length of each sequence
  Uint min_poslen;		       // minimum positive sequence length
  Uint min_neglen;		       // minimum negative sequence length
  Uint max_poslen;		       // maximum positive sequence length
  Uint max_neglen;		       // maximum negative sequence length
  double avg_poslen;		       // average length of positive sequences
  double avg_neglen;		       // average length of negative sequences
  BOOL use_binomial;		       // use the Cumulative binomial distribution
  double bernoulli;		       // probability of a hit in a positive sequence
} Multiseq;                  // \Typedef{Multiseq}

/*
  The following type describes how to format a sequence description.
*/

typedef struct
{
  BOOL defined,          // show a description
       replaceblanks,    // replaceblanks by underscore
       untilfirstblank;  // only show sequence until first blank
  Uint skipprefix,       // always skip this number of prefixes
       maxlength;        // maximal number of chars of description to be shown
} Showdescinfo;

/*
  The following type is used to store some basic information about
  a sequence stored in a \texttt{Multiseq}-record.
*/

typedef struct
{
  Uint seqnum,       // the sequence number in multiseq
       seqstartpos,  // the position of the first character in multiseq.sequence
       seqlength,    // the length of the sequence
       relposition;  // the relative position of the sequence
} Seqinfo;           // \Typedef{Seqinfo}

//\Ignore{

#endif

//}
