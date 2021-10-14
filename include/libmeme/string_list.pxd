from libc.stdio cimport FILE

from libmeme.array cimport ARRAY_T


cdef extern from "string-list.h" nogil:

    cdef struct string_list_t:
        pass
    ctypedef string_list_t STRING_LIST_T

    cdef STRING_LIST_T* new_string_list()
    cdef int max_string_length(STRING_LIST_T* a_list)
    cdef int get_num_strings(STRING_LIST_T* a_list)
    cdef char* get_nth_string(int n, STRING_LIST_T* a_list)
    cdef void set_nth_string(char* new_string, int n, STRING_LIST_T* a_list)
    cdef double get_nth_score(int n, STRING_LIST_T* a_list)
    cdef void set_nth_score(double new_score, int n, STRING_LIST_T* a_list)
    cdef ARRAY_T* get_string_list_scores(STRING_LIST_T* a_list)
    cdef void append_to_nth_string(char* new_string, int n, STRING_LIST_T* a_list)
    cdef void prepend_to_strings(char* new_string, STRING_LIST_T* a_list)
    cdef bint have_space(STRING_LIST_T* a_list)
    cdef void add_string(char* a_string, STRING_LIST_T* a_list)
    cdef void add_string_with_score(char* a_string, STRING_LIST_T* a_list, double score)
    cdef void add_strings(STRING_LIST_T* source_list, STRING_LIST_T* target_list)
    cdef void remove_string(char* a_string, STRING_LIST_T* a_list)
    cdef void remove_strings(STRING_LIST_T* source_list, STRING_LIST_T* target_list)
    cdef bint have_string(char* a_string, STRING_LIST_T* a_list)
    cdef bint have_substring(char* a_string, STRING_LIST_T* a_list)
    cdef bint have_regex(char* a_regex, STRING_LIST_T* a_list)
    cdef STRING_LIST_T* get_matches_in_string_list(char* a_regex, STRING_LIST_T* a_list)
    cdef int get_index_in_string_list(char* a_string, STRING_LIST_T* a_list)
    cdef STRING_LIST_T* copy_string_list(STRING_LIST_T* string_list)
    cdef void sort_string_list(STRING_LIST_T* a_list)
    cdef void sort_string_list_by_score(STRING_LIST_T* a_list, bint reverse)
    cdef void overlap_string_list(STRING_LIST_T* list_a, STRING_LIST_T* list_b, STRING_LIST_T** intersection, STRING_LIST_T** a_minus_b, STRING_LIST_T** b_minus_a)
    cdef char* combine_string_list(STRING_LIST_T* a_list, char* separator)
    cdef STRING_LIST_T* read_string_list(FILE* infile)
    cdef STRING_LIST_T* read_string_list_from_file(char* filename)
    cdef STRING_LIST_T* new_string_list_char_split(char separator, char* string)
    cdef void write_string_list(char* separator, STRING_LIST_T* a_list, FILE* outfile)
    cdef void free_string_list(STRING_LIST_T* a_list)
    cdef bint equal_string_lists(STRING_LIST_T* a_list, STRING_LIST_T* b_list)
    cdef bint has_duplicates(char* message, STRING_LIST_T* my_list)
