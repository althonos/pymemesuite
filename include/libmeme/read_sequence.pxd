from libc.stdio cimport FILE

from libmeme.alphabet cimport ALPH_T


cdef extern from "read_sequence.h" nogil:

    const size_t RCHUNK   = 100
    const size_t MAXDELEN = 10000

    ctypedef enum FORMAT_TYPE:
        FASTA
        SWISSPROT

    cdef bint read_sequence(
      ALPH_T *alph,
      FILE *data_file,
      char **sample_name,
      char **sample_de,
      char **sequence,
      long *length
    )
