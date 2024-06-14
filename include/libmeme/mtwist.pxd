from libc.stdint cimport uint32_t
from libc.stdio cimport FILE


cdef extern from * nogil:
# cdef extern from "mtwist.h" nogil:

    const size_t MT_STATE_SIZE

    ctypedef struct mt_state:
        uint32_t statevec[MT_STATE_SIZE]
        int      stateptr
        int      initialize

    cdef void mts_mark_initialized(mt_state* state)
    cdef void mts_seed32(mt_state* state, uint32_t seed)
    cdef void mts_seed32new(mt_state* state, uint32_t seed)
    cdef void mts_seedfull(mt_state* state, uint32_t seeds[MT_STATE_SIZE])
    cdef uint32_t mts_seed(mt_state* state)
    cdef uint32_t mts_goodseed(mt_state* state)
    cdef void mts_bestseed(mt_state* state)
    cdef void mts_refresh(mt_state* state)
    cdef int mts_savestate(FILE* statefile, mt_state* state)
    cdef int mts_loadstate(FILE* statefile, mt_state* state)

    cdef void mt_seed32(uint32_t seed)
    cdef void mt_seed32new(uint32_t seed)
    cdef void mt_seedfull(uint32_t seeds[MT_STATE_SIZE])
    cdef uint32_t mt_seed()
    cdef uint32_t mt_goodseed()
    cdef void mt_bestseed()
    cdef mt_state* mt_getstate()
    cdef int mt_savestate(FILE* statefile)
    cdef int mt_loadstate(FILE* statefile);

    cdef extern double mt_drand()
    cdef extern double mt_ldrand()
