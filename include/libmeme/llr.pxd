from libmeme.array cimport ARRAY_T


cdef extern from "llr.h" nogil:
    cdef void init_llr_pv_tables(int min, int max, int alength, ARRAY_T* back, bint pal)
    cdef double get_llr_pv(double llr, double n, int w, int range, double frac, int alength, ARRAY_T* back)
    cdef double get_llr_mean(double n)
