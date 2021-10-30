from libc.stdio cimport FILE

from libmeme.alphabet cimport ALPH_T
from libmeme.meme cimport DATASET, SAMPLE
from libmeme.array_list cimport ARRAYLST_T

cdef extern from "read_sequence.h" nogil:
    cdef DATASET* read_seq_file(char* file_name, ALPH_T* alph, bint use_comp, bint ignore_dup)
    cdef DATASET* create_meme_dataset_from_momo(ARRAYLST_T* seq_array, ALPH_T* alph, int width, bint eliminate_repeats)
    cdef SAMPLE* get_sample_by_name(char* sample_name)
    cdef void add_control_samples(DATASET* dataset, DATASET* control, int c1, int c2, bint use_comp)
    cdef void shuffle_dataset_order(DATASET* dataset)
    cdef void shuffle_dataset_letters(DATASET* dataset, int kmer, int revcomp)
    cdef void free_sample(SAMPLE *sample);
