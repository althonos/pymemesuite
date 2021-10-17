from libc.stdio cimport FILE

from libmeme.alphabet cimport ALPH_T
from libmeme.meme cimport DATASET, MODEL, CANDIDATE, MOTIF_SUMMARY, THETA


cdef extern from "display.h" nogil:
    cdef void print_results(DATASET* dataset, MODEL* model, CANDIDATE* candidates, FILE* outfile)
    cdef void record_results(DATASET* dataset, MODEL* model, MOTIF_SUMMARY* motif_summaries)

    cdef void print_meme_file_xml(MODEL* model, DATASET* dataset, DATASET* neg_dataset, int nmotifs, MOTIF_SUMMARY* motif_summaries, char* stopping_reason, char* xml_filename)

    cdef void print_meme_file_html(char* stylesheet_file_path, char* input_file_path, char* output_file_path)

    cdef void print_theta(int imotif, char* consensus, int format, int nsites, THETA theta, int w, double log_ev, char* str_space, DATASET* dataset, FILE* outfile)

    cdef void print_zij(DATASET* dataset, MODEL* model)
    cdef void print_wij(DATASET* dataset)

    cdef char* get_consensus(THETA theta, int w, DATASET* dataset, int N, double min_prob)

    cdef void print_command_summary(MODEL* model, DATASET* dataset, FILE* outfile)

    cdef void print_dataset_summary (DATASET* dataset, FILE* outfile)

    cdef void print_summary(MODEL* model, DATASET* dataset, int nmotifs, FILE* outfile)

    cdef void print_meme_doc(FILE* outfile)

    cdef char *get_single_letter_consensus(MODEL* model, ALPH_T* alphabet)
