/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef ARRAYDEF_H
#define ARRAYDEF_H
#include "st_types.h"
#include "st_errordef.h"
#include "st_spacedef.h"

//}

/*
  This file defines macros to conveniently declare and 
  manipulate dynamic arrays whose size grow on demand. Each dynamic 
  array over some type \texttt{T}
  is implemented by a structure consisting of three components:
  \begin{enumerate}
  \item
  \texttt{space\#\#T} is a pointer to the space block of type \texttt{T}
  allocated for the array.
  \item
  \texttt{allocated\#\#T} is an \texttt{Uint} storing the number 
  of entries in the array currently allocated.
  \item
  \texttt{nextfree\#\#T} holds the smallest index of the array where no 
  value is stored.
  \end{enumerate}
  Here \texttt{\#\#} is the concatenation operator of the C-preprocessor.
  The following macro expands to a corresponding type definition over 
  some given \texttt{TYPE}.
*/

#define DECLAREARRAYSTRUCT(TYPE)\
        typedef struct\
        {\
          TYPE *space##TYPE;\
          Uint allocated##TYPE, nextfree##TYPE;\
        } Array##TYPE

/*
  \texttt{INITARRAY} initializes an empty array.
*/

#define INITARRAY(A,TYPE)\
        { \
        (A)->space##TYPE = NULL;\
        (A)->allocated##TYPE = (A)->nextfree##TYPE = 0;\
        }

/*
  \texttt{CHECKARRAYSPACE} checks if the integer \texttt{nextfree\#\#T}
  points to an index for which the space is not allocated yet. If this is 
  the case, the number of cells allocated is incremented by \texttt{L}. The 
  contents of the previously filled array elements is of course maintained.
*/

// TLB added: use malloc
#define CHECKARRAYSPACE_TLB(A,TYPE,L) \
        if((A)->nextfree##TYPE + (L) >= (A)->allocated##TYPE) \
        { \
          (A)->allocated##TYPE += L; \
          Uint size = (Uint) sizeof(TYPE) * (Uint) (A)->allocated##TYPE; \
          (A)->space##TYPE = (TYPE *) realloc((A)->space##TYPE, sizeof(TYPE) * (Uint) (A)->allocated##TYPE); \
        }\
        NOTSUPPOSEDTOBENULL((A)->space##TYPE) \


#define CHECKARRAYSPACE(A,TYPE,L)\
        if((A)->nextfree##TYPE >= (A)->allocated##TYPE)\
        {\
          (A)->allocated##TYPE += L;\
          (A)->space##TYPE\
             = (TYPE *) allocandusespaceviaptr(__FILE__,(Uint) __LINE__,\
                                               (A)->space##TYPE,\
                                               (Uint) sizeof(TYPE),\
                                               (A)->allocated##TYPE);\
        }\
        NOTSUPPOSEDTOBENULL((A)->space##TYPE)

/*
  The next macro is a variation of \texttt{CHECKARRAYSPACE}, 
  which checks if the next
  \texttt{L} cells have been allocated. If not, then this is done.
*/

#define CHECKARRAYSPACEMULTI(A,TYPE,L)\
        if((A)->nextfree##TYPE + (L) >= (A)->allocated##TYPE)\
        {\
          (A)->allocated##TYPE += L;\
          (A)->space##TYPE\
             = (TYPE *) allocandusespaceviaptr(__FILE__,(Uint) __LINE__,\
                                               (A)->space##TYPE,\
                                               (Uint) sizeof(TYPE),\
                                               (A)->allocated##TYPE);\
        }\
        NOTSUPPOSEDTOBENULL((A)->space##TYPE)

/*
  This macro checks the space and delivers a pointer \texttt{P} 
  to the next free element in the array.
*/

#define GETNEXTFREEINARRAY(P,A,TYPE,L)\
        CHECKARRAYSPACE_TLB(A,TYPE,L);\
        P = (A)->space##TYPE + (A)->nextfree##TYPE++;

/*
  This macro checks the space and stores \texttt{V} in the 
  \texttt{nextfree}-component of the array. \texttt{nextfree}
  is incremented.
*/

#define STOREINARRAY(A,TYPE,L,VAL)\
        CHECKARRAYSPACE_TLB(A,TYPE,L);\
        (A)->space##TYPE[(A)->nextfree##TYPE++] = VAL; \

// TLB added: push 
#define PUSHARRAY(A, TYPE, L, VAL)\
	STOREINARRAY(A,TYPE,L,VAL)
// TLB added: pop top item and return its pointer
#define POPARRAY(A, TYPE)\
        ((A)->nextfree##TYPE==0 ? (TYPE *) NULL : ((A)->space##TYPE + --(A)->nextfree##TYPE))
// TLB added: peak at given item with given index in stack
#define PEEKARRAY(A, INDEX, TYPE)\
        ((INDEX) >= (A)->nextfree##TYPE ? (TYPE *) NULL : ((A)->space##TYPE + (INDEX)))
// TLB added: peak at top item (top of stack pointer)
#define TOPARRAY(A, TYPE)\
        ((A)->nextfree##TYPE==0 ? (TYPE *) NULL : ((A)->space##TYPE + (A)->nextfree##TYPE - 1))
// TLB added: get top of stack index
#define TOPINDEXARRAY(A, TYPE)\
        ((A)->nextfree##TYPE - 1)
// TLB added: check if stack is empty
#define EMPTYARRAY(A, TYPE)\
        ((A)->nextfree##TYPE == 0)
// TLB added: sort an array
#define SORTARRAY(A, TYPE, COMPARE)\
	qsort((A)->space##TYPE, (A)->nextfree##TYPE, sizeof(TYPE), (Qsortcomparefunction)COMPARE)
// TLB added: sort an array
#define SORTARRAY_R(A, TYPE, COMPARE, THUNK)\
	qsort_r((A)->space##TYPE, (A)->nextfree##TYPE, sizeof(TYPE), &(THUNK), (Qsortrcomparefunction)COMPARE)

/*
  This macro frees the space for an array if it is not \texttt{NULL}.
*/

// TLB added: use malloc
#define FREEARRAY_TLB(A,TYPE)\
        if((A)->space##TYPE != NULL)\
        {\
          free((A)->space##TYPE);\
        }

#define FREEARRAY(A,TYPE)\
        if((A)->space##TYPE != NULL)\
        {\
          FREESPACE((A)->space##TYPE);\
        }

/* 
  Some declarations for the most common array types.
*/

DECLAREARRAYSTRUCT(Uchar);
DECLAREARRAYSTRUCT(Ushort);
DECLAREARRAYSTRUCT(char);
DECLAREARRAYSTRUCT(Uint);
DECLAREARRAYSTRUCT(Sint);
DECLAREARRAYSTRUCT(PairUint);
DECLAREARRAYSTRUCT(ThreeUint);

/*
  And some type synonyms.
*/

typedef ArrayUint  ArrayPosition;       // \Typedef{ArrayPosition}
typedef ArrayUchar ArrayCharacters;     // \Typedef{ArrayCharacters}

/*
  The following array type has some extra components. However, it can be 
  manipulated by the macros above since the record-components
  \texttt{spaceStrings}, \texttt{nextfreeStrings}, and 
  \texttt{allocatedStrings} is declared appropriately.
*/

typedef struct
{
  Stringtype *spaceStrings;
  Uchar *stringbuffer;
  Uint stringbufferlength, nextfreeStrings, allocatedStrings;
} ArrayStrings;   // \Typedef{ArrayStrings}

//\Ignore{

#endif

//}
