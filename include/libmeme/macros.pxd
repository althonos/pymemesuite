cdef extern from "macros.h" nogil:

    const double BITS
    const double LOGZERO

    cdef double LOGEV(double logn, double logp)
    cdef double exp10_logx(double logx, double m, double e, double prec)

    cdef void free_2array(void** array, int rows)
