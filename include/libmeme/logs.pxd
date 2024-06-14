from libmeme.macros cimport BITS


cdef extern from "logs.h" nogil:

    const double log_precision
    const double log_table[2*(<const int> log_precision) + 2]

    cdef double LOGL_Y(double x)
    cdef double LOGL_I(double x)
    cdef double LOGL_LOW(double x)
    cdef double LOGL_HI(double x)
    cdef double LOGL(double x)

    const double exp_precision
    const double exp_table[(<const int> BITS) * (<const int> exp_precision) + 2]

    cdef double EXPL_Y(double x)
    cdef double EXPL_I(double x)
    cdef double EXPL_LOW(double x)
    cdef double EXPL_HI(double x)
    cdef double EXPL(double x)

    cdef double LOGL_SUM1(double logx, double logy)
    cdef double LOG_SUM(double logx, double logy)

    cdef void init_log()
    cdef void init_exp()
