from libmeme.meme cimport SAMPLE

cdef extern from "message.h" nogil:
    cdef void balance_loop(SAMPLE** samples, int n_samples)
    cdef void balance_loop1(SAMPLE** samples, int n_samples)
