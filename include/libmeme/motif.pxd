from libcpp cimport bool
from libc.stdio cimport FILE

from libmeme.alphabet cimport ALPH_T
from libmeme.array cimport ARRAY_T
from libmeme.array_list cimport ARRAYLST_T
from libmeme.matrix cimport MATRIX_T, MTYPE
from libmeme.mtwist cimport mt_state
from libmeme.red_black_tree cimport RBTREE_T
from libmeme.string_builder cimport STR_T

cdef extern from "motif.h" nogil:

    const char NON_MOTIF_ID_CHAR
    const size_t MAX_MOTIF_ID_LENGTH
    const size_t MAX_MOTIF_URL_LENGTH

    cdef struct motif_t:
        pass
    ctypedef motif_t MOTIF_T

    cdef int get_motif_idx(MOTIF_T* motif)
    cdef char* get_motif_id(MOTIF_T* motif)
    cdef char* get_motif_st_id(MOTIF_T* motif)
    cdef char* get_motif_id2(MOTIF_T* motif)
    cdef char* get_motif_consensus(MOTIF_T* motif)
    cdef char get_motif_strand(MOTIF_T* motif)
    cdef MATRIX_T* get_motif_freqs(MOTIF_T* motif)
    cdef MATRIX_T* get_motif_scores(MOTIF_T* motif)
    cdef int get_motif_length(MOTIF_T* motif)
    cdef int get_motif_trimmed_length(MOTIF_T* motif)
    cdef ALPH_T* get_motif_alph(MOTIF_T* motif)
    cdef bool get_motif_strands(MOTIF_T* motif)
    cdef void set_motif_strands(MOTIF_T* motif)
    cdef int get_motif_alph_size(MOTIF_T* motif)
    cdef int get_motif_ambiguous_size(MOTIF_T* motif)
    cdef double get_motif_evalue(MOTIF_T* motif)
    cdef double get_motif_log_evalue(MOTIF_T* motif)
    cdef double get_motif_complexity(MOTIF_T* motif)
    cdef double get_motif_nsites(MOTIF_T* motif)
    cdef MTYPE get_motif_counts(int position, MOTIF_T* motif)
    cdef char* get_motif_url(MOTIF_T* motif)
    cdef bool has_motif_url(MOTIF_T* motif)
    cdef int get_motif_trim_left(MOTIF_T* motif)
    cdef int get_motif_trim_right(MOTIF_T* motif)

    cdef void clear_motif_trim(MOTIF_T* motif)
    cdef bool has_motif_zeros(MOTIF_T* motif)

    cdef bool have_motif(char* motif_id, int num_motifs, MOTIF_T* motifs)

    cdef void copy_motif(MOTIF_T* source, MOTIF_T* dest)

    cdef MATRIX_T* convert_matrix_alphabet(MATRIX_T* int, MTYPE value, ALPH_T* source_alp, ALPH_T* target_alph)

    cdef void shuffle_motif(MOTIF_T*, mt_state* prng)

    cdef MATRIX_T* convert_freqs_into_scores(ALPH_T* alph, MATRIX_T* freqs, ARRAY_T* bg, int site_count, double pseudo_count)
    cdef MATRIX_T* convert_scores_into_freqs(ALPH_T* alph, MATRIX_T* scores, ARRAY_T* bg, int site_count, double pseudo_count)

    cdef void reverse_complement_motif(MOTIF_T* motif)

    cdef void apply_pseudocount_to_motif(MOTIF_T* motif, ARRAY_T* background, double pseudo_count)

    void calc_motif_ambigs(MOTIF_T* motif)
    void normalize_motif(MOTIF_T* motif, double tolerance)

    cdef void trim_motif_by_bit_threshold(MOTIF_T* motif, double threshold_bits)

    cdef double compute_motif_complexity(MOTIF_T* motif)

    cdef int get_info_content_position(bool from_start, float threshold, ARRAY_T* background, MOTIF_T* motif)

    cdef char* get_best_possible_match(MOTIF_T* motif)

    cdef MOTIF_T* dup_rc_motif(MOTIF_T* motif)

    cdef MOTIF_T* duplicate_motif(MOTIF_T* motif)

    cdef void free_motif(MOTIF_T* a_motif)
    cdef void destroy_motif(void* a_motif)

    cdef void motif_list_to_array(ARRAYLST_T* motif_list, MOTIF_T** motif_array, int* num)
    cdef void motif_tree_to_array(RBTREE_T* motif_tree, MOTIF_T** motif_array, int* num);

    cdef MOTIF_T* motif_at(MOTIF_T* array_of_motifs, int index)

    cdef void free_motif_array(MOTIF_T* motif_array, int num)

    cdef void free_motifs(ARRAYLST_T* motif_list)

    cdef void add_reverse_complements(ARRAYLST_T* motifs)

    cdef void dump_motif_freqs(FILE* out, MOTIF_T* m)

    cdef void read_motif_alphabets(ARRAYLST_T* motif_sources, bool xalph, ALPH_T** alph)

    cdef int ap_cmp(const void* p1, const void* p2)

    cdef int idx_cmp(const void* p1, const void* p2)

    cdef void motif2consensus(MOTIF_T* motif, STR_T* consensus, bool single_letter)


cdef extern from "motif-spec.h" nogil:

    cdef MOTIF_T* allocate_motif(char* id, char* id2, ALPH_T* alph, MATRIX_T* freqs, MATRIX_T* scores)
    cdef void set_motif_id(const char* id, int len, MOTIF_T* motif)
    cdef void set_motif_id2(const char* id, int len, MOTIF_T* motif)
    cdef void set_motif_strand(char strand, MOTIF_T *motif)
    cdef void set_motif_nsites(MOTIF_T* motif, int num_sites)
    cdef void set_motif_url(char* url, MOTIF_T *motif)
