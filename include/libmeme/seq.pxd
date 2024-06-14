from libcpp cimport bool

from libc.stdint cimport int8_t

from libmeme.alphabet cimport ALPH_T
from libmeme.array cimport ARRAY_T


cdef extern from "seq.h" nogil:

    cdef struct seq:
        pass
    ctypedef seq SEQ_T

    const int SEQ_KEEP
    const int SEQ_NOAMBIG

    cdef SEQ_T* allocate_seq(char* name, char* description, unsigned int offset, char* sequence)

    cdef char* get_seq_name(SEQ_T* sequence)
    cdef char* get_seq_description(SEQ_T* sequence)
    cdef unsigned int get_seq_length(SEQ_T* sequence)
    cdef unsigned int get_seq_offset(SEQ_T* sequence)
    cdef void set_seq_offset(unsigned int offset, SEQ_T* sequence)
    cdef unsigned int get_seq_starting_coord(SEQ_T* sequence)
    cdef void set_seq_starting_coord(unsigned int start, SEQ_T* sequence)
    cdef float get_seq_weight(SEQ_T* sequence)
    cdef void set_seq_weight(float weight, SEQ_T* sequence)
    cdef unsigned char get_seq_char(int index, SEQ_T* sequence)
    cdef void set_seq_char(int index, char character, SEQ_T* sequences)
    cdef int get_seq_int(int index, SEQ_T* sequence)
    cdef bool has_seq_gc(SEQ_T* sequence)
    cdef int get_seq_gc(int index, SEQ_T* sequence)

    cdef char* get_raw_sequence(SEQ_T* sequence)
    cdef void set_raw_sequence(char* raw_sequence, bool is_complete, SEQ_T* sequence)
    cdef int8_t* get_isequence(SEQ_T* sequence)
    cdef int* get_int_sequence(SEQ_T* sequence)
    cdef int* get_gc_sequence(SEQ_T* sequence)
    cdef double get_total_gc_sequence(SEQ_T* sequence)
    cdef void set_total_gc_sequence(SEQ_T* sequence, double gc)

    cdef unsigned int get_seq_num_priors(SEQ_T* sequence)
    cdef double* get_seq_priors(SEQ_T* sequence)
    cdef void set_seq_priors(double* priors, SEQ_T* sequence)

    cdef SEQ_T* copy_sequence(SEQ_T* sequence)

    cdef void index_sequence(SEQ_T* sequence, ALPH_T* alpha, int options)

    cdef void set_complete(bool is_complete, SEQ_T* sequence)
    cdef bool is_complete(SEQ_T* sequence)

    cdef void remove_flanking_xs(SEQ_T* sequence)

    cdef void prepare_sequence(SEQ_T* sequence, ALPH_T* alpha, bool hard_mask)

    cdef void shift_sequence(int offset, SEQ_T* sequence)

    cdef int get_max_seq_length(int num_seqs, SEQ_T** sequences)
    cdef int get_max_seq_name(int num_seqs, SEQ_T** sequences)

    cdef void set_sequence_weights(char* weight_filename, int num_seqs, SEQ_T** sequences)

    cdef ARRAY_T* get_sequence_freqs(SEQ_T* seqs, ALPH_T* alph)

    cdef char* get_raw_subsequence(int start, int stop, SEQ_T* sequence)

    cdef SEQ_T* shuffle_seq(SEQ_T* sequence, int kmer, int i_copy)

    cdef void free_seq(SEQ_T* sequence)
