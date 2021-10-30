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

    cdef void srand_mt(uint32_t seed)
    cdef void shuffle(void *base, size_t n, size_t size)
    cdef double drand_mt()
    cdef long random_mt()
    cdef int rand_int(const unsigned n)
    cdef double rand_dbl(const double n)
