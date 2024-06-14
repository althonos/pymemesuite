from libcpp cimport bool

from libmeme.alphabet cimport ALPH_T
from libmeme.array cimport ARRAY_T
from libmeme.array_list cimport ARRAYLST_T
from libmeme.motif cimport MOTIF_T


cdef extern from "motif-in-flags.h" nogil:

    const int OPEN_MFILE
    const int CALC_AMBIGS
    const int SCANNED_SITES
    const int SKIP_POST_PROCESSING


cdef extern from "motif-in.h" nogil:

    cdef struct mread:
        pass
    ctypedef mread MREAD_T

    cdef MREAD_T* mread_create(const char* filename, int options, bool symmetrical)
    void mread_destroy(MREAD_T* mread)

    void mread_update(MREAD_T* mread, const char* buffer, size_t size, bint end)

    ARRAYLST_T* mread_load(MREAD_T* mread, ARRAYLST_T* motifs)

    bool mread_has_motif(MREAD_T* mread)
    MOTIF_T* mread_next_motif(MREAD_T* mread)

    void mread_set_conversion(MREAD_T* mread, ALPH_T* alph, const ARRAY_T* bg)
    void mread_set_background(MREAD_T* mread, const ARRAY_T* bg, ALPH_T* alph)
    void mread_set_bg_source(MREAD_T* mread, const char* source, ALPH_T* alph)
    void mread_set_pseudocount(MREAD_T* mread, double pseudocount)
    void mread_set_trim(MREAD_T* mread, double trim_bits)

    char* mread_get_other_bg_src(MREAD_T* mread)
    ALPH_T* mread_get_alphabet(MREAD_T* mread)
    int mread_get_strands(MREAD_T* mread)
    ARRAY_T* mread_get_background(MREAD_T* mread)
    ARRAYLST_T* mread_get_motif_occurrences(MREAD_T* mread)
    void destroy_motif_occurrences(ARRAYLST_T* motif_occurrences)
