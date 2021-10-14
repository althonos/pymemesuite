from libc.stdio cimport FILE


cdef extern from "array.h" nogil:

    ctypedef double ATYPE

    cdef struct array_t:
        int    num_items
        ATYPE  key
        ATYPE* items
    ctypedef array_t ARRAY_T

    cdef ARRAY_T* allocate_array(const int length)
    cdef ARRAY_T* resize_array(ARRAY_T* array, const int num_items)
    cdef ARRAY_T* resize_init_array(ARRAY_T* array, const int num_items, const ATYPE default_value)

    cdef int get_array_length(const ARRAY_T* array)
    cdef ATYPE get_array_item(const int index, ARRAY_T* array)
    cdef void  set_array_item(const int index, ATYPE value, ARRAY_T* array)
    cdef void  incr_array_item(const int index, ATYPE value, ARRAY_T* array)
    cdef ATYPE get_array_minimum(ARRAY_T* array)
    cdef ATYPE get_array_maximum(ARRAY_T* array)

    cdef ATYPE* raw_array(ARRAY_T* array)
    cdef void init_array(ATYPE value, ARRAY_T* array)
    cdef void remove_array_item(int item_index, ARRAY_T* array)
    cdef void fill_array(ATYPE* raw_array, ARRAY_T* array)
    cdef void copy_array(const ARRAY_T* source_array, ARRAY_T* target_array)
    cdef bint equal_arrays(ATYPE close_enough, ARRAY_T* array1, ARRAY_T* array2)
    cdef void sum_array(ARRAY_T* array1, ARRAY_T* array2)
    cdef void element_product(ARRAY_T* array1, ARRAY_T* array2)
    cdef double dot_product(ARRAY_T* array1, ARRAY_T* array2)
    cdef double euclidean_distance(ARRAY_T* array1, ARRAY_T* array2)
    cdef void scalar_add(ATYPE value, ARRAY_T* array)
    cdef void scalar_mult(ATYPE value, ARRAY_T* array)
    cdef ATYPE array_total(ARRAY_T* array)
    cdef ATYPE total_subarray(int start_index, int length, ARRAY_T* array)
    cdef bint is_sorted(bint good_score_is_low, ARRAY_T* my_array)
    cdef void sort_array(bint reverse_sort, ARRAY_T* array)
    cdef void set_array_key(ATYPE key, ARRAY_T* array)
    cdef void get_array_key(ARRAY_T* array)
    cdef ATYPE sum_of_squares(ARRAY_T* array)
    cdef ATYPE sum_of_square_diffs(ARRAY_T* array1, ARRAY_T* array2)
    cdef ATYPE sum_of_square_diffs(ARRAY_T* array1, ARRAY_T* array2)
    cdef ATYPE compute_median(ARRAY_T* array)

    cdef void read_array(FILE* infile, ARRAY_T* array)
    cdef ARRAY_T* read_array_from_file(const char* filename)
    cdef void print_array(ARRAY_T* array, int width, int precision, bint eol, FILE* outfile)
    cdef void print_sub_array(int start, int end, ARRAY_T* array, int width, int precision, bint eol, FILE* outfile)

    cdef void free_array(ARRAY_T* array)

    cdef void dot_divide(ARRAY_T* array1, ARRAY_T* array2)
    cdef ARRAY_T* bootstrap_array(ARRAY_T* source_array, int maxsize)
    cdef ATYPE ave_array(ARRAY_T* array)
    cdef ATYPE array_variance(ARRAY_T* array)
    cdef void sum_to_zero(ARRAY_T* array)
    cdef void variance_one_array(ARRAY_T* array)
    cdef void normalize(ATYPE close_enough, ARRAY_T* array)
    cdef ATYPE correlation(ARRAY_T* array1, ARRAY_T* array2)
    cdef void log_array(ARRAY_T* array)
    cdef void unlog_array(ARRAY_T* array)
    cdef ATYPE log_array_total(ARRAY_T* array)
    cdef void log_normalize(ATYPE close_enough, ARRAY_T* array)
    cdef void convert_to_from_log_array(bint to_lot, ARRAY_T* source_array, ARRAY_T* target_array)
    cdef void mix_log_arrays(float mixing, ARRAY_T* array1, ARRAY_T* array2)
    cdef void randomize_array(ATYPE magnitude, ARRAY_T* array)
    cdef void add_noise(ATYPE magnitude, ARRAY_T* array)
    cdef void all_positive(ARRAY_T* array)
    cdef void normalize_subarray(int start_index, int length, double tolerance, ARRAY_T* array)
    cdef ARRAY_T* extract_subarray(ARRAY_T* array, const int start, const int end)
