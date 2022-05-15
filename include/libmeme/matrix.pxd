from libcpp cimport bool
from libc.stdio cimport FILE

from libmeme.array cimport ARRAY_T, ATYPE


cdef extern from "matrix.h" nogil:

    ctypedef ATYPE MTYPE

    cdef struct matrix_t:
        int num_rows
        int num_cols
        ARRAY_T** rows
    ctypedef matrix_t MATRIX_T

    cdef MATRIX_T* allocate_matrix(int num_rows, int num_cols)
    cdef void free_matrix(MATRIX_T* matrix)

    cdef MATRIX_T* duplicate_matrix(MATRIX_T* matrix)
    cdef MATRIX_T* convert_matrix(double** theta, int w, int alength)
    cdef void grow_matrix(ARRAY_T* one_row, MATRIX_T* matrix)
    cdef void remove_matrix_row(int row_index, MATRIX_T* matrix)
    cdef void remove_matrix_col(int col_index, MATRIX_T* matrix)

    cdef void resize_matrix(int row_count, int col_count, MTYPE value, MATRIX_T* matrix)

    cdef int get_num_rows(MATRIX_T* matrix)
    cdef int get_num_cols(MATRIX_T* matrix)
    cdef ARRAY_T* get_matrix_row(int row, MATRIX_T* matrix)
    cdef void set_matrix_row(int row, ARRAY_T* one_row, MATRIX_T* matrix)

    cdef MTYPE get_matrix_cell(int row, int col, MATRIX_T* matrix)
    cdef void set_matrix_cell(int row, int col, MTYPE value, MATRIX_T* matrix)
    cdef void incr_matrix_cell(int row, int col, MTYPE value, MATRIX_T* matrix)

    cdef ARRAY_T* get_matrix_column(int i_col, MATRIX_T* matrix)
    cdef ARRAY_T* get_matrix_column2(int i_col, ARRAY_T* buffer, MATRIX_T* matrix)
    cdef void set_matrix_column(ARRAY_T* column, int i_col, MATRIX_T* matrix)
    cdef MATRIX_T* array_to_matrix(bool one_row, ARRAY_T* array)

    cdef void copy_matrix(MATRIX_T* source_matrix, MATRIX_T* target_matrix)

    cdef void init_matrix(MTYPE value, MATRIX_T* matrix)

    cdef void fill_matrix(MTYPE* raw_matrix, MATRIX_T* matrix)

    cdef void sum_matrices(MATRIX_T* matrix1, MATRIX_T* matrix2)

    cdef ARRAY_T* extract_diagonal(MATRIX_T* matrix)

    cdef bool is_symmetric(bool verbose, MTYPE slop, MATRIX_T* matrix)

    cdef MATRIX_T* average_across_diagonal(MATRIX_T* matrix1, MATRIX_T* matrix2)

    cdef void add_to_diagonal(MTYPE value, MATRIX_T* matrix)

    cdef bool equal_matrices(ATYPE close_enough, MATRIX_T* matrix1, MATRIX_T* matrix2)

    cdef void scalar_mult_matrix(MTYPE value, MATRIX_T* matrix)
    cdef void scalar_add_matrix(MTYPE value, MATRIX_T* matrix)

    cdef void mult_matrix(MATRIX_T* matrix1, MATRIX_T* matrix2)

    cdef void convert_to_from_log_matrix(bool to_log, MATRIX_T* the_matrix, MATRIX_T* the_log_matrix)
    cdef void mix_log_matrices(float mixing, MATRIX_T* matrix1, MATRIX_T* matrix2)

    cdef MATRIX_T* read_matrix(FILE* infile)
    cdef MATRIX_T* read_known_matrix(int num_rows, int num_cols, FILE* infile)
    cdef void print_matrix(MATRIX_T* matrix, int width, int precision, bool print_titles, FILE* outfile)

    cdef void randomize_matrix(MTYPE max_value, MATRIX_T* matrix)
    cdef MTYPE sum_of_matrix(MATRIX_T* matrix)
    cdef MTYPE sum_of_squares_matrix(MATRIX_T* matrix)
    cdef MTYPE sum_of_square_diff_matrices(MATRIX_T* matrix1, MATRIX_T* matrix2)

    cdef void zero_mean_matrix_rows(MATRIX_T* matrix)
    cdef void zero_mean_matrix_cols(MATRIX_T* matrix)

    cdef void variance_one_matrix_rows(MATRIX_T* matrix)

    cdef MATRIX_T* matrix_multiply(MATRIX_T* matrix1, MATRIX_T* matrix2)

    cdef void normalize_rows(float tolerance, MATRIX_T* matrix)
    cdef void normalize_matrix(float tolerance, MATRIX_T* matrix)

    cdef ARRAY_T* get_matrix_row_sums(MATRIX_T* matrix)
    cdef ARRAY_T* get_matrix_col_sums(MATRIX_T* matrix)

    cdef void sort_matrix_rows(bool reverse_sort, ARRAY_T* keys, MATRIX_T* matrix)

    cdef void permute_matrix(MATRIX_T* matrix, bool cols, int* permutation, int count)
    cdef void shuffle_matrix_cols(MATRIX_T* matrix, bool cols)
    cdef void shuffle_matrix(MATRIX_T* matrix)

    cdef void get_matrix_rows(int, int, MATRIX_T*, MATRIX_T*)
