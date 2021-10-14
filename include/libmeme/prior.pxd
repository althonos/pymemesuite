from libmeme.alphabet cimport ALPH_T

# "prior.h" does not have include guards, so importing it
# here will cause errors about redefinitions - we use `extern from *`
# since "prior.h" will already be included in "meme.h"
cdef extern from * nogil:

    ctypedef struct PriorLib:
        ALPH_T*  alph
        int      L
        double*  Mix
        double*  B
        double** Distr
        int*     FullUpdate
        int*     QUpdate
        char**   StructID
        char**   Comment

    ctypedef struct PriorLib:
        pass

    cdef PriorLib* alloc_PriorLib(int L, ALPH_T* alph)
    cdef void free_PriorLib(PriorLib* lib)
    cdef PriorLib* read_PriorLib(char* plib_name, double desired_beta, ALPH_T* custom_alph)

    cdef void mixture_regularizer(double* freq, PriorLib* Lib, double* reg)

    double logpajgy(double* y, PriorLib* Lib, int j, int Calc)
    double logpygaj(double* y, double* a, int AlphLength)
