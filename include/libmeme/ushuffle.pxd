cdef extern from "ushuffle.h" nogil:
    cdef void ushuffle(const char* s, char* t, long l, long k)
    cdef void ushuffle1(const char* s, long l, long k)
    cdef void ushuffle2(char* t)

    ctypedef long(*randfunc_t)()
    cdef void set_randfunc(randfunc_t randfunc)
    cdef void permutec(char* t, long l)
