from libc.math cimport HUGE_VAL, INFINITY
from libc.stdint cimport uint32_t

cdef extern from "utils.h" nogil:
    const double BIG
    const double LITTLE
    const double SMALL_POS

    ctypedef int VERBOSE_T

    cdef enum:
        INVALID_VERBOSE = 0
        QUIET_VERBOSE   = 1
        NORMAL_VERBOSE  = 2
        HIGH_VERBOSE    = 3
        HIGHER_VERBOSE  = 4
        DUMP_VERBOSE    = 5

    cdef extern VERBOSE_T verbosity

    cdef srand_mt(uint32_t seed)
