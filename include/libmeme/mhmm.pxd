from libcpp cimport bool

from libmeme.alphabet cimport ALPH_T
from libmeme.array cimport ARRAY_T
from libmeme.matrix cimport MATRIX_T
from libmeme.motif cimport MAX_MOTIF_ID_LENGTH


cdef extern from "pssm.h" nogil:
    cdef struct pssm:
        pass
    ctypedef pssm PSSM_T


cdef extern from "mhmm-state.h" nogil:

    ctypedef enum STATE_T:
        INVALID_STATE
        START_STATE
        START_MOTIF_STATE
        MID_MOTIF_STATE
        END_MOTIF_STATE
        SPACER_STATE
        END_STATE
    cdef extern char** STATE_STRS
    cdef extern int NUM_STATE_T

    cdef struct mhmm_state:
        STATE_T type
        int ntrans_out
        int* itrans_out
        ARRAY_T* trans_out
        int ntrans_in
        int* itrans_in
        ARRAY_T* trans_in
        ARRAY_T* emit
        double num_sites
        int i_motif
        int w_motif
        char motif_id[MAX_MOTIF_ID_LENGTH + 2]
        int i_position
        char id_char
        ARRAY_T* emit_odds
        double alpha
        PSSM_T* pssm
        PSSM_T* npssm
        double min_sig_score
        double min_pvalue
    ctypedef mhmm_state MHMM_STATE_T

    ctypedef enum HMM_T:
        INVALID_HMM
        LINEAR_HMM
        COMPLETE_HMM
        STAR_HMM
    cdef extern char** HMM_STRS
    cdef extern int NUM_HMM_T

    cdef struct mhmm:
        HMM_T type
        bool log_odds
        int num_motifs
        int num_states
        int num_spacers
        int spacer_states
        ALPH_T* alph
        ARRAY_T* background
        char* description
        char* motif_file
        char* sequence_file
        MHMM_STATE_T* states
        MATRIX_T* trans
        int* hot_states
        int num_hot_states
    ctypedef mhmm MHMM_T

    cdef MHMM_T* allocate_mhmm(ALPH_T* alph, int num_states)
    cdef char* get_state_motif_id(bool stranded, MHMM_STATE_T* this_state)
    cdef void copy_mhmm(MHMM_T* src_mhmm, MHMM_T* dst_mhmm)
    cdef void mix_mhmm(const float mixing, MHMM_T* const new_hmm, MHMM_T* old_mhmm)
    cdef void free_mhmm(MHMM_T* an_mhmm)

    cdef void compute_ins_and_outs(MHMM_T* the_hmm, const bool log_form)
    cdef char choose_consensus(ALPH_T* alph, bool log_space, MHMM_STATE_T* a_state)
    cdef bool has_fims(MHMM_T* the_hmm)
    cdef int get_num_motifs(MHMM_T* the_hmm)
