from libc.stdio cimport FILE


cdef extern from "array-list.h" nogil:

    cdef struct arraylst_t:
        pass
    ctypedef arraylst_t ARRAYLST_T

    cdef ARRAYLST_T *arraylst_create_sized(int size);
    cdef ARRAYLST_T *arraylst_create();
    cdef void arraylst_destroy(void (*optional_item_destructor)(void*), ARRAYLST_T *arraylst);
    cdef void arraylst_free(ARRAYLST_T *arraylst);
    cdef void arraylst_preallocate(int size, ARRAYLST_T *arraylst);
    cdef void arraylst_fit(ARRAYLST_T *arraylst);
    cdef bint arraylst_is_empty(ARRAYLST_T *arraylst);
    cdef int arraylst_size(ARRAYLST_T *arraylst);
    cdef void arraylst_put_n(int times, int index, void *item, ARRAYLST_T *arraylst);
    cdef void arraylst_add_n(int times, void *item, ARRAYLST_T *arraylst);
    cdef void arraylst_put(int index, void *item, ARRAYLST_T *arraylst);
    cdef void arraylst_add(void *item, ARRAYLST_T *arraylst);
    cdef void arraylst_set(int index, void *item, ARRAYLST_T *arraylst);
    cdef void arraylst_swap(int index1, int index2, ARRAYLST_T *arraylst);
    cdef void *arraylst_get(int index, ARRAYLST_T *arraylst);
    cdef void *arraylst_peek(ARRAYLST_T *arraylst);
    cdef void *arraylst_remove_range(int index1, int index2, void (*optional_item_destructor)(void*), ARRAYLST_T *arraylst);
    cdef void arraylst_clear(void (*optional_item_destructor)(void*), ARRAYLST_T *arraylst);
    cdef void *arraylst_remove(int index, ARRAYLST_T *arraylst);
    cdef void *arraylst_take(ARRAYLST_T *arraylst);
    cdef void arraylst_map_range(void* (*map_fun)(void*), int index, int count, ARRAYLST_T *arraylst);
    cdef void arraylst_map(void* (*map_fun)(void*), ARRAYLST_T *arraylst);
    cdef void arraylst_apply_range(void (*fun)(void*), int index, int count, ARRAYLST_T *arraylst);
    cdef void arraylst_apply(void (*fun)(void*), ARRAYLST_T *arraylst);
    cdef void* arraylst_accumulate_range(void (*accumulator_fun)(void*, void*), void *initval, int index, int count, ARRAYLST_T *arraylst);
    cdef void* arraylst_accumulate(void (*accumulator_fun)(void*, void*), void *initval, ARRAYLST_T *arraylst);
    cdef void arraylst_qsort(int(*compar)(const void *, const void *), ARRAYLST_T *arraylst);
    cdef int arraylst_bsearch(const void *key, int (*compar)(const void *, const void *), ARRAYLST_T *arraylst);
    cdef int arraylst_compar_txt(const void *v1, const void *v2);
