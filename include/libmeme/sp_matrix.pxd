from libc.stdio cimport FILE

from libmeme.alphabet cimport ALPH_T
from libmeme.hash_table cimport HASH_TABLE
from libmeme.heap cimport HEAP
from libmeme.meme cimport DATASET, S_POINT, MAXSITE, MODEL, BRANCH_PARAMS
from libmeme.partition cimport PARTITION


cdef extern from "sp_matrix.h" nogil:

    cdef struct sp_matrix:
        S_POINT**   matrix
        int*        central_ws
        int         n_ws
        int*        central_ns
        int         n_ns
        PARTITION** partitions
        int         n_parts
    ctypedef sp_matrix SP_MATRIX

    cdef SP_MATRIX* create_spoint_matrix(int* central_ws, int n_ws, int* central_ns, int n_ns, DATASET* dataset)
    cdef S_POINT* create_spoint_row(int width, int* central_ws, int n_ws, int* central_ns, int n_ns, DATASET* dataset)
    cdef int get_hs(int curr_width, int min_width, int max_width, int w_dist, int main_w_hs, DATASET* dataset)
    cdef void transfer_final_scores(ALPH_T* alph, SP_MATRIX* sp_matrix)
    cdef void branching_search(BRANCH_PARAMS* branch_params, MODEL* model, DATASET* dataset, SP_MATRIX* sp_mat, HASH_TABLE evaluated_seed_ht)

    cdef int sp_get_num_rows(SP_MATRIX* sp_mat)
    cdef int sp_get_num_cols(SP_MATRIX* sp_mat)
    cdef S_POINT** get_sp_matrix(SP_MATRIX* sp_mat)
    cdef int get_num_central_widths(SP_MATRIX* sp_mat)
    cdef int get_num_central_nsites(SP_MATRIX* sp_mat)
    cdef int* get_central_ws(SP_MATRIX* sp_mat)
    cdef int* get_central_ns(SP_MATRIX* sp_mat)
    cdef int get_min_width(SP_MATRIX* sp_mat)
    cdef int get_max_width(SP_MATRIX* sp_mat)
    cdef int get_min_nsites(SP_MATRIX* sp_mat)
    cdef int get_max_nsites(SP_MATRIX* sp_mat)
    cdef S_POINT* get_sp_arr(SP_MATRIX* sp_mat, int row_idx)
    cdef S_POINT* get_spoint(SP_MATRIX* sp_mat, int row_idx, int col_idx)
    cdef S_POINT* get_s_point(SP_MATRIX* sp_mat, int width, int nsites)

    cdef void update_s_point_heaps(S_POINT* s_points, char* str_seed, int n_nsites0)
    cdef S_POINT* choose_em_starts(SP_MATRIX* sp_matrix, DATASET* data, MODEL* model, int* n_starts)

    cdef void free_s_point(S_POINT* sp)
    cdef void copy_s_point(S_POINT* s_point, S_POINT* sp_copy)
    cdef void print_s_point(S_POINT* s_point, FILE* out)
    cdef void free_sp_matrix(SP_MATRIX* sp_mat)
    cdef void print_sp_matrix(SP_MATRIX* sp_mat, FILE* out)

    const size_t SEED_PACKET_LENGTH
    cdef struct seed_packet:
        double score
        int width
        int iseq
        int ipos
        int num_seed_packets
        int nsites0
        char seed[MAXSITE]
    ctypedef seed_packet SEED_PACKET

    cdef void reduce_across_heaps(S_POINT* s_points, int n_nsites0)
    cdef HEAP* create_heap_from_sp_matrix(SP_MATRIX* sp_mat)
