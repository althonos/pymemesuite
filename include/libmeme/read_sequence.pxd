from libc.stdio cimport FILE

from libmeme.alphabet cimport ALPH_T


cdef extern from "read_sequence.h" nogil:

    const size_t RCHUNK
    const size_t MAXDELEN

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
