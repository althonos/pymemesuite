from libc.stdio cimport FILE


cdef extern from "partition.h" nogil:

    cdef struct partition:
        int min_w
        int max_w
        int central_w
        int min_n
        int max_n
        int central_n
    ctypedef partition PARTITION

    cdef PARTITION* new_partition(int part_min_w, int part_max_w, int central_w, int part_min_n, int part_max_n, int central_n)
    cdef void print_partition(PARTITION* part, FILE* out)
