from libc.stdio cimport FILE

from libmeme.hash_table cimport HASH_TABLE


cdef extern from "heap.h" nogil:

    cdef struct heap:
        int max_size
        int cur_size
        int height
        int next_node
        void** node_list
        int (*compare)(void *p1, void* p2)
        void *(*copy)(void* p)
        void (*destroy)(void* p)
        char* (*get_key)(void* p)
        # void (*print)(FILE* f, void* p)
        HASH_TABLE ht
    ctypedef heap HEAP

    cdef void* add_node_heap(HEAP* heap, void* node)
    cdef void *pop_heap_root(HEAP* h)
    cdef HEAP* copy_heap(HEAP* h)
    cdef void  destroy_heap(HEAP* h)
    cdef int   get_num_nodes(HEAP* h)
    cdef void* get_node(HEAP* h, int i)
    cdef int   get_max_size(HEAP* h)
    cdef int   get_best_node(HEAP* h)
    cdef void  print_heap(FILE* outfile, HEAP* heap)
