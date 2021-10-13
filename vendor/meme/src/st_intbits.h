/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef INTBITS_H
#define INTBITS_H
#include <limits.h>
#include "st_types.h"
#include "st_errordef.h"
#include "st_spacedef.h"

//}

/*
  This file contains some definitions manipulating bitvectors represented
  by a \texttt{Uint}. In the comment lines we use $w$ for the word size
  and \texttt{\symbol{94}} for exponentiation of the previous character.
*/

#define INTWORDSIZE\
        (UintConst(1) << LOGWORDSIZE)     // # of bits in Uint = w
#define FIRSTBIT\
        (UintConst(1) << (INTWORDSIZE-1)) // \(10^{w-1}\)
#define ISBITSET(S,I)\
        (((S) << (I)) & FIRSTBIT)         // is \(i\)th bit set?
#define ITHBIT(I)\
        (FIRSTBIT >> (I))                 // \(0^{i}10^{w-i-1}\) 
#define SECONDBIT\
        (FIRSTBIT >> 1)                   // \(010^{w-2}\)
#define THIRDBIT\
        (FIRSTBIT >> 2)                   // \(0010^{w-3}\)
#define FIRSTTWOBITS\
        (UintConst(3) << (INTWORDSIZE-2)) // \(11^{w-2}\)
#define EXCEPTFIRSTBIT\
        (~FIRSTBIT)                       // \(01^{w-1}\)
#define EXCEPTFIRSTTWOBITS\
        (EXCEPTFIRSTBIT >> 1)             // \(001^{w-2}\)
#define EXCEPTFIRSTTHREEBITS\
        (EXCEPTFIRSTBIT >> 2)             // \(0001^{w-3}\)
#define DIVWORDSIZE(I)\
        ((I) >> LOGWORDSIZE)              // \((I) div w\)
#define MODWORDSIZE(I)\
        ((I) & (INTWORDSIZE-1))           // \((I) mod w\)
#define MULWORDSIZE(I)\
        ((I) << LOGWORDSIZE)              // \((I) * w\)

/*
  The following macro allocates a bitarray of \texttt{N} bits. All bits
  are off.
*/

#define INITBITTAB(TAB,N)\
        {\
          Uint *tabptr, tabsize = 1 + DIVWORDSIZE(N);\
          TAB = ALLOCSPACE(NULL,Uint,tabsize);\
          for(tabptr = TAB; tabptr < (TAB) + tabsize; tabptr++)\
          {\
            *tabptr = 0;\
          }\
        }

// TLB: added for clarity
#define FREEBITTAB(TAB)\
	{\
	if ((TAB) != NULL) FREESPACE_TLB(TAB) \
        (TAB) = NULL; \
	}
	
// TLB: redefined for speed
#define INITBITTAB_TLB(TAB,N)\
        {\
          Uint *tabptr, tabsize = 1 + DIVWORDSIZE(N);\
          TAB = (Uint *)calloc(tabsize, sizeof(Uint));\
        }

/*
  The following macro inititalizes a bitarray such tha all bits
  are off.
*/

#define CLEARBITTAB(TAB,N)\
        {\
          Uint *tabptr, tabsize = 1 + DIVWORDSIZE(N);\
          for(tabptr = TAB; tabptr < TAB + tabsize; tabptr++)\
          {\
            *tabptr = 0;\
          }\
        }

// TLB: ADDED
// Copy the right operand into the left operand.
#define COPYBITTAB(TAB1, TAB2, N)\
        {\
          Uint *tabptr1, *tabptr2, tabsize = 1 + DIVWORDSIZE(N);\
          for(tabptr1=TAB1,tabptr2=TAB2; tabptr1 < TAB1 + tabsize; tabptr1++, tabptr2++)\
          {\
            *tabptr1 = *tabptr2;\
          }\
        }

// TLB added:
// Bitwise-OR two bit tables of the same size putting the result in the left operand.
#define ORBITTABS(TAB1,TAB2,N)\
        {\
          Uint *tabptr1, *tabptr2, tabsize = 1 + DIVWORDSIZE(N);\
          for(tabptr1=TAB1,tabptr2=TAB2; tabptr1 < TAB1 + tabsize; tabptr1++, tabptr2++)\
          {\
            *tabptr1 |= *tabptr2;\
          }\
        }

// TLB added:
// Bitwise-AND two bit tables of the same size putting the result in the left operand.
#define ANDBITTABS(TAB1,TAB2,N)\
        {\
          Uint *tabptr1, *tabptr2, tabsize = 1 + DIVWORDSIZE(N);\
          for(tabptr1=TAB1,tabptr2=TAB2; tabptr1 < TAB1 + tabsize; tabptr1++, tabptr2++)\
          {\
            *tabptr1 &= *tabptr2;\
          }\
        }

// TLB added:
// Bitwise-AND TAB1 with Bitwise-NOT of TAB2 and put result in TAB1, and vice-versa.
// TAB3 is used for scratch.
#define MUTUALNANDBITTABS(TAB1,TAB2,TAB3,N)\
        {\
          Uint *tabptr1, *tabptr2, *tabptr3, tabsize = 1 + DIVWORDSIZE(N);\
          for(tabptr1=TAB1,tabptr2=TAB2,tabptr3=TAB3; tabptr1 < TAB1 + tabsize; tabptr1++, tabptr2++, tabptr3++)\
          {\
            *tabptr3 = *tabptr1 & (~*tabptr2);\
            *tabptr2 = *tabptr2 & (~*tabptr1);\
            *tabptr1 = *tabptr3; \
          }\
        }

// TLB added:
// Count the number of set bits in a bit table using Brian Kernighan's algorithm.
#define COUNTBITTAB(CNT,TAB,N)\
        {\
	  CNT = 0;\
          Uint *tabptr, tabsize = 1 + DIVWORDSIZE(N);\
          for(tabptr = TAB; tabptr < TAB + tabsize; tabptr++)\
          {\
	    Uint n = *tabptr;\
	    while (n)\
	      {\
                n &= (n-1);\
		CNT++;\
	      }\
          }\
        }

// TLB added:
// Return a list of bit numbers of the set bits in a bit table.
#define LISTBITTAB(LIST,TAB,N)\
        {\
          Uint *tabptr, tabsize = 1 + DIVWORDSIZE(N);\
	  Uint offset = 0;\
	  Uint index = 0;\
          for(tabptr = TAB; tabptr < TAB + tabsize; tabptr++, offset+=INTWORDSIZE)\
          {\
	    Uint n = *tabptr, tmp;\
	    while (n)\
	    {\
	       /* Get a mask with just the lowest "1" bit from n set. */ \
	       tmp = (n ^ (n-1)) + 1; \
	       tmp = (tmp==0) ? ITHBIT(0) : tmp >> 1; \
	       /* Get the position of the "1" bit. */ \
	       int lo = 0, mid = INTWORDSIZE/2, hi = INTWORDSIZE-1; \
	       while (tmp != power2_table[mid]) { \
		 if (tmp < power2_table[mid]) { \
		   hi = mid; \
		   mid = lo + (hi-lo)/2; \
		 } else { \
		   lo = mid; \
		   mid = hi - (hi-lo)/2; \
		 } \
	       } \
	       int bitno = INTWORDSIZE-1-mid; \
	       /* Mask out the lowest "1" bit in n. */ \
	       n &= ~tmp; \
	       LIST[index++] = bitno + offset; \
	    }\
          }\
        }

// TLB added:
// Initialize the power table needed by binary search in LISTBITTAB
#include "macros.h"
EXTERN Uint power2_table[INTWORDSIZE];
#define INITPOWER2TABLE()\
	{\
	  int i; \
	  power2_table[0] = 1; \
	  for (i=1; i<INTWORDSIZE; i++) power2_table[i] = 2 * power2_table[i-1]; \
	}

// TLB added:
// Print the set bit indices.
#define PRINTBITTAB(TAB,N)\
        {\
          Uint *tabptr, tabsize = 1 + DIVWORDSIZE(N), bitno, w = INTWORDSIZE;\
          for(tabptr = TAB, bitno = 0; tabptr < TAB + tabsize; tabptr++, bitno+=w)\
          {\
	    Uint n = *tabptr, offset = 0;\
            for (i=0; n && i<w; i++) {\
	      {\
		if (n & FIRSTBIT) {\
	          fprintf(stderr, " %u", bitno+i); \
		} \
		n = n << 1;\
	      }\
            }\
          }\
        }

/*
  \texttt{SETIBIT(TAB,I)} sets the \texttt{I}-th bit in bitarray 
  \texttt{TAB} to 1.
*/

#define SETIBIT(TAB,I)    (TAB)[DIVWORDSIZE(I)] |= ITHBIT(MODWORDSIZE(I))

/*
  \texttt{UNSSETIBIT(TAB,I)} sets the \texttt{I}-th bit in bitarray 
  \texttt{TAB} to 0.
*/

#define UNSETIBIT(TAB,I)  (TAB)[DIVWORDSIZE(I)] &= ~(ITHBIT(MODWORDSIZE(I)))

/*
  \texttt{ISIBITSET(TAB,I)} checks if the \texttt{I}-th bit in bitarray 
  \texttt{TAB} is 1.
*/

#define ISIBITSET(TAB,I)  ((TAB)[DIVWORDSIZE(I)] & ITHBIT(MODWORDSIZE(I)))

//\Ignore{

#ifdef __cplusplus
  extern "C" {
#endif
  
char *intbits2string(Uint bs);

#ifdef __cplusplus
}
#endif

#endif

//}
