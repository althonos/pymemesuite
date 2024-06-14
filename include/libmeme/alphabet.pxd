from libc.stdint cimport uint8_t, uint32_t, uint64_t, int64_t
from libc.stdio cimport FILE
from libc.limits cimport UCHAR_MAX

from libmeme.array cimport ARRAY_T
from libmeme.data_types cimport LCB_T
from libmeme.json_writer cimport JSONWR_T
from libmeme.string_builder cimport STR_T


cdef extern from "alphabet.h" nogil:

    ctypedef enum ALPHABET_T:
        Dna
        Rna
        Protein
        Custom

    ctypedef struct ALPH_T:
        uint64_t  references
        int       flags
        char*     alph_name
        int       ncore
        int       nfull
        char*     symbols
        char**    aliases
        char**    names
        uint32_t* colours
        uint8_t*  ncomprise
        uint8_t** comprise
        uint8_t*  complement
        uint8_t   encode[UCHAR_MAX+1]
        uint8_t   encode2core[UCHAR_MAX+1]
        uint8_t   encodesafe[UCHAR_MAX+1]
        uint8_t   encodesafe2core[UCHAR_MAX]

    ctypedef struct XLATE_T:
        ALPH_T*   src_alph
        ALPH_T*   dest_alph
        uint8_t   src_nsym
        uint8_t   dest_nsym
        uint32_t* xlate

    cdef struct bgcalc:
        pass
    ctypedef bgcalc BGCALC_T


    const int ALPH_FLAG_BUILTIN
    const int ALPH_FLAG_EXTENDS
    const int ALPH_FLAG_EXTENDS_DNA
    const int ALPH_FLAG_EXTENDS_RNA
    const int ALPH_FLAG_EXTENDS_PROTEIN
    const int ALPH_CASE_INSENSITIVE

    cdef int alph_sym_cmp(const void* v1, const void* v2)
    cdef int alph_str_cmp(const void* v1, const void* v2)
    cdef ALPH_T* alph_dna()
    cdef ALPH_T* alph_rna()
    cdef ALPH_T* alph_protein()
    cdef ALPH_T* alph_generic()
    cdef ALPH_T* alph_load()
    cdef int alph_pick(int nalphs, ALPH_T** alphs, char* symbols, int64_t* counts)
    cdef ALPH_T* alph_guess(char* symbols, int64_t* counts)
    cdef ALPH_T* alph_hold(ALPH_T* alphabet)
    cdef void alph_release(ALPH_T* alphabet)
    cdef void alph_print_head(ALPH_T* alphabet, FILE* out)
    cdef void alph_print(ALPH_T* alphabet, bint header, FILE* out)
    cdef void alph_print_xml(ALPH_T* alphabet, char* tag, char* pad, char* ident, FILE* out)
    cdef void alph_print_json(ALPH_T* alph, JSONWR_T* jsonwr)
    cdef bint alph_equal(const ALPH_T* a1, const ALPH_T* a2)
    cdef int alph_core_subset(const ALPH_T* sub_alph, const ALPH_T* super_alph)

    cdef bint alph_extends(const ALPH_T* alph)
    cdef bint alph_extends_rna(const ALPH_T* alph)
    cdef bint alph_extends_dna(const ALPH_T* alph)
    cdef bint alph_extends_protein(const ALPH_T* alph)

    cdef bint alph_is_builtin(const ALPH_T* alph)
    cdef bint alph_is_builtin_rna(const ALPH_T* alph)
    cdef bint alph_is_builtin_dna(const ALPH_T* alph)
    cdef bint alph_is_builtin_protein(const ALPH_T* alph)

    cdef const char* alph_name(const ALPH_T* alph)

    cdef int alph_size_pairs(const ALPH_T* a)
    cdef int alph_size_core(const ALPH_T* a)
    cdef int alph_size_wild(const ALPH_T* a)
    cdef int alph_size_all(const ALPH_T* a)
    cdef int alph_size_ambig(const ALPH_T* a)
    cdef int alph_size_full(const ALPH_T* a)

    cdef bint alph_has_complement(const ALPH_T* alph)
    cdef bint alph_is_case_insensitive(const ALPH_T* alph)
    cdef uint8_t alph_index(const ALPH_T* alph, uint8_t letter)
    cdef uint8_t alph_indexc(const ALPH_T* alph, uint8_t letter)
    cdef uint8_t alph_encode(const ALPH_T* alph, uint8_t letter)
    cdef uint8_t alph_encodec(const ALPH_T* alph, uint8_t letter)

    cdef int alph_wild(const ALPH_T* alph)
    cdef uint8_t alph_complement(const ALPH_T* alph, int index)
    cdef char alph_char(const ALPH_T* alph, int index)

    cdef bint alph_has_sym_name(ALPH_T* alph, int index)
    cdef char* alph_sym_name(ALPH_T* alph)

    cdef char* alph_xml_id(ALPH_T* alph, int index, STR_T* buffer)

    cdef uint8_t alph_ncomprise(const ALPH_T* alph, int index)
    cdef uint8_t* alph_comprise(ALPH_T* alph, int index, int cindex)

    cdef uint32_t alph_colour(const ALPH_T* alph, int index)
    cdef uint8_t alph_colour_r(const ALPH_T* alph, index)
    cdef uint8_t alph_colour_g(const ALPH_T* alph, index)
    cdef uint8_t alph_colour_b(const ALPH_T* alph, index)

    cdef char* alph_aliases(ALPH_T* alph, int index)
    cdef const char* alph_string(ALPH_T* alph, STR_T* buffer)
    cdef char alph_wildcard(const ALPH_T* alph)
    cdef bint alph_is_known(const ALPH_T* alph, uint8_t letter)
    cdef bint alph_is_concrete(const ALPH_T* alph, uint8_t letter)
    cdef bint alph_is_core(const ALPH_T* alph, uint8_t letter)
    cdef bint alph_is_ambiguous(const ALPH_T* alph, uint8_t letter)
    cdef bint alph_is_wildcard(const ALPH_T* alph, uint8_t letter)
    cdef bint alph_is_prime(ALPH_T* alph, char letter)

    cdef int get_num_ambiguous_letters(ALPH_T* alph, char* string, int length)

    cdef bint alph_test(ALPH_T** alpha, int index, char letter)
    cdef ALPH_T* alph_type(const char* alphabet, int max)

    cdef void calc_ambigs(ALPH_T* alph, bint log_space, ARRAY_T* array)
    cdef void average_freq_with_complement(ALPH_T* alph, ARRAY_T* freqs)
    cdef ARRAY_T* get_mast_frequencies(ALPH_T* alph, bint has_ambigs, bint translate)
    cdef ARRAY_T* get_nrdb_frequencies(ALPH_T* alph, ARRAY_T* freqs)
    cdef ARRAY_T* get_uniform_frequencies(ALPH_T* alph, ARRAY_T* freqs)
    cdef void normalize_frequencies(ALPH_T* alph, ARRAY_T* freqs, double pseudo)
    cdef ARRAY_T* get_file_frequencies(ALPH_T* alph, char* bg_filename)

    cdef ARRAY_T* load_markov_model(ALPH_T* alph, int* order, const char* bg_filename) except NULL
    cdef ARRAY_T* load_markov_model_without_alph(const char* bg_filename, int* order, char** syms)
    cdef ARRAY_T* calculate_markov_model(ALPH_T* alph, int order, double epsilon, bint join_seq, const char* seq, BGCALC_T** save)
    cdef void average_rc_markov_model(ALPH_T* alph, int order, ARRAY_T* bg)
    cdef void resize_markov_model(int asize0, int asize1, ARRAY_T* tuples, int* order_p)

    cdef enum ambig_calc:
        SUM_FREQS
        SUM_LOGS
        AVG_FREQS
    ctypedef ambig_calc AMBIG_CALC_EN

    cdef void extend_markov_model(ALPH_T* alph, bint wildcard_only, AMBIG_CALC_EN method, ARRAY_T* tuples)
    cdef void extrapolate_markov_model(int asize0, int asize1, double ambig_fraction, ARRAY_T* tuples)

    cdef double calculate_log_cumulative_background(ALPH_T* alph, bint wildcard_only, int order, ARRAY_T* a_cp, const char* seq, LCB_T* logcumback)

    cdef ARRAY_T* get_background(ALPH_T* alph, char* bg_filename)

    cdef void dist_ambigs(ALPH_T* alph, ARRAY_T* freqs)

    cdef void complement_swap_freqs(ALPH_T* alph, ARRAY_T* a1, ARRAY_T* a2)

    const int ALPH_NO_ALIASES
    const int ALPH_NO_AMBIGS
    const int ALPH_NO_UNKNOWN

    cdef void translate_seq(ALPH_T* alph, char* sequence, int flags)
    cdef void invcomp_seq(ALPH_T* alph, char* sequence, long length)
    cdef bint is_palindrome(ALPH_T* alph, char* word)
    cdef void seq2r(ALPH_T* alph, uint8_t* res, char* seq, int len)
    cdef char comp_sym(const ALPH_T* alph, uint8_t sym)

    cdef int* dhash_seq(ALPH_T* alph, XLATE_T* trans, bint full_alph, char* sequence, long length)
    cdef XLATE_T* xlate_dna2protein()
    cdef void xlate_destroy(XLATE_T* translator)
    cdef void xlate_print(XLATE_T* translator, FILE* out)
    cdef ALPH_T* xlate_src_alph(XLATE_T* translator)
    cdef uint8_t xlate_src_nsyms(XLATE_T* translator)
    cdef ALPH_T* xlate_dest_alph(XLATE_T* translator)
    cdef int xlate_pos(XLATE_T* xlate, bint invcomp, const char* sequence)

    cdef char* fasta_get_markov(int argc, char** argv, bint tmp_file)
    cdef ARRAY_T* get_markov_from_sequences(char* seqfile, int* order, double pseudo, ALPH_T* alph, char* alph_file, ALPHABET_T alphabet_type, bint rc)

    cdef uint32_t xlate_index(XLATE_T* translator, bint invcomp, const char* sequence)
    cdef uint32_t xlate_index2(XLATE_T* translator, uint32_t position)

    cdef bint alph_check(ALPH_T* alph, char* syms)
