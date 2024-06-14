from libcpp cimport bool

from libmeme.alphabet cimport ALPH_T
from libmeme.array cimport ARRAY_T, ATYPE
from libmeme.matrix cimport MATRIX_T
from libmeme.mhmm cimport MHMM_T
from libmeme.motif cimport MOTIF_T
from libmeme.prior_dist cimport PRIOR_DIST_T

cdef extern from "pssm.h" nogil:

    const size_t PSSM_RANGE

    cdef struct pssm:
        MATRIX_T* matrix
        ALPH_T* alph
        int alphsize
        int w
        bool matrix_is_log
        bool matrix_is_scaled
        double scale
        double offset
        int range
        ARRAY_T* pv
        int num_gc_bins
        ARRAY_T** gc_pv
        int min_score
        int max_score
        double nolog_max_score
        MOTIF_T* motif
    ctypedef pssm PSSM_T

    cdef struct pssm_pair:
        PSSM_T* pos_pssm
        PSSM_T* neg_pssm
        int num_gc_bins
        MATRIX_T** gc_n_pv_lookup
        ARRAY_T* scaled_to_ama
        MOTIF_T* motif
    ctypedef pssm_pair PSSM_PAIR_T

    cdef int get_pssm_w(PSSM_T* pssm)
    cdef ALPH_T* get_pssm_alph(PSSM_T* pssm)
    cdef int get_pssm_alphsize(PSSM_T* pssm)
    cdef double get_pssm_scale(PSSM_T* pssm)
    cdef double get_pssm_offset(PSSM_T* pssm)
    cdef int get_pssm_pv_length(PSSM_T* pssm)
    cdef ATYPE get_pssm_pv(int score, PSSM_T* pssm)
    cdef ATYPE get_pssm_score(int row, int col, PSSM_T* pssm)
    cdef ATYPE get_pssm_min_pv(PSSM_T* pssm)
    cdef double get_pssm_max_raw_score(PSSM_T* pssm)

    cdef void set_up_pssms_and_pvalues (
        bool motif_scoring,
        double p_threshold,
        bool use_both_strands,
        bool allow_weak_motifs,
        MHMM_T *the_hmm,
        PRIOR_DIST_T *prior_dist,
        double alpha
    )

    cdef void compute_motif_score_matrix(
        bool use_pvalues,
        double p_threshold,
        int* int_sequence,
        size_t seq_length,
        double *priors,
        size_t num_priors,
        double alpha,
        MHMM_T* the_hmm,
        MATRIX_T** motif_score_matrix
    )

    cdef void scale_score_matrix(
        MATRIX_T* matrix,
        int rows,
        int cols,
        PRIOR_DIST_T* prior_dist,
        double alpha,
        int range,
        double* offset_p,
        double* scale_p
    )

    cdef void scale_pssm(
        PSSM_T* pssm,
        PRIOR_DIST_T* prior_dist,
        double alpha,
        int range
    )

    cdef ARRAY_T *scale_prior_dist(
        ARRAY_T* priors,
        int range,
        double scale,
        double offset
    )

    cdef void get_pv_lookup_pos_dep(
        PSSM_T* pssm,
        MATRIX_T* background_matrix,
        ARRAY_T* scaled_prior_dist
    )

    cdef void get_pv_lookup(
        PSSM_T* pssm,
        ARRAY_T* background,
        ARRAY_T* scaled_prior_dist
    )

    cdef double scaled_to_raw(double x, double w, double scale, double offset)
    cdef double pssm_scaled_to_raw(double x, PSSM_T* pssm)
    cdef double raw_to_scaled(double x, double w, double scale, double offset)

    cdef double get_unscaled_pssm_score(double score, PSSM_T* pssm)
    cdef double get_scaled_pssm_score(double score, PSSM_T* pssm)

    cdef PSSM_T* build_motif_pssm(
        MOTIF_T* motif,
        ARRAY_T* bg_freqs,
        ARRAY_T* pv_bg_freqs,
        PRIOR_DIST_T* prior_dist,
        double alpha,
        int range,
        int num_gc_bins,
        bool no_log
    )
    cdef PSSM_T* build_matrix_pssm(
        ALPH_T* alph,
        MATRIX_T* matrix,
        ARRAY_T* bg_freqs,
        PRIOR_DIST_T *prior_dist,
        double alpha,
        int range
    )

    cdef double get_ama_pv(
        double ama_score,
        int seqlen,
        double seq_gc,
        PSSM_PAIR_T* pssm_pair
    )

    PSSM_PAIR_T* create_pssm_pair(PSSM_T* pos_pssm, PSSM_T* neg_pssm)
    cdef void free_pssm_pair(PSSM_PAIR_T *pssm_pair)

    cdef PSSM_T* allocate_pssm(
        ALPH_T* alph,
        int w,
        int alphsize,
        int num_gc_bins
    );

    cdef void free_pssm(PSSM_T* pssm)

    cdef double pssm_best_match_score(PSSM_T* pssm)
    cdef void print_pssm(PSSM_T* pssm)

    cdef double pssm_pvalue_to_score(PSSM_T *pssm, double pvalue)
