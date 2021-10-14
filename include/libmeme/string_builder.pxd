from libc.stdint cimport int32_t


cdef extern from "string-builder.h" nogil:

    cdef struct str:
        pass
    ctypedef str STR_T

    cdef STR_T* str_create(int capacity)
    cdef char*  str_destroy(STR_T* strb, int keep_string)
    cdef void   str_append(STR_T* strb, const char* str, int len)
    cdef void   str_append2(STR_T* strb, const char* str)
    # cdef void   str_vappendf(STR_T *strb, const char *fmt, va_list ap)
    # cdef void str_appendf(STR_T *strb, const char *fmt, ...)
    cdef void str_append_code(STR_T *strb, int32_t code_point)
    cdef void str_insert(STR_T* strb, int offset, const char* str, int len)
    cdef void str_insert2(STR_T* strb, int offset, const char* str)
    # cdef void str_vinsertf(STR_T *strb, int offset, const char *fmt, va_list ap)
    # cdef void str_insertf(STR_T *strb, int offset, const char *fmt, ...)
    cdef void str_insert_code(STR_T *strb, int offset, int32_t code_point)
    cdef void str_replace(STR_T *strb, int start, int end, const char *str, int len)
    cdef void str_replace2(STR_T *strb, int start, int end, const char *str)
    # cdef void str_vreplacef(STR_T *strb, int start, int end, const char *fmt, va_list ap)
    # cdef void str_replacef(STR_T *strb, int start, int end, const char *fmt, ...)
    cdef void str_set(STR_T *strb, const char *str, int len)
    # cdef void str_vsetf(STR_T *strb, const char *fmt, va_list ap)
    # cdef void str_setf(STR_T *strb, const char *fmt, ...)
    cdef void str_delete(STR_T *strb, int start, int end)
    cdef void str_truncate(STR_T *strb, int length)
    cdef void str_clear(STR_T *strb)
    cdef size_t str_len(STR_T *strb)
    cdef char* str_internal(STR_T *strb)
    cdef char str_char(STR_T *strb, int pos)
    cdef char* str_copy(STR_T *strb)
    cdef char* str_subcopy(STR_T *strb, int start, int end)
    cdef int str_cmp(const STR_T *strb, const char *s2)
    cdef int str_ncmp(const STR_T *strb, const char *s2, size_t n)
    cdef int str_casecmp(const STR_T *strb, const char *s2)
    cdef void str_fit(STR_T *strb)
    cdef void str_ensure(STR_T *strb, int new_capacity)
    # cdef void str_append_path(STR_T *strb, int segments, ...)
    cdef void str_append_evalue(STR_T *strb, double log10_ev, int prec)
    cdef char* str_evalue(STR_T *strb, double log10_ev, int prec)
