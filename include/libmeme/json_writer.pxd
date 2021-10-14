from libc.stdio cimport FILE
from posix.time cimport time_t


cdef extern from "json-writer.h" nogil:

    cdef struct json_writer:
        pass
    ctypedef json_writer JSONWR_T

    cdef JSONWR_T* jsonwr_open(FILE *dest, bint minimize, int min_cols, int tab_cols, int line_cols)
    cdef void jsonwr_close(JSONWR_T* jsonwr)
    cdef void jsonwr_property(JSONWR_T* jsonwr, const char* property)
    cdef void jsonwr_start_object_value(JSONWR_T* jsonwr)
    cdef void jsonwr_end_object_value(JSONWR_T* jsonwr)
    cdef void jsonwr_start_array_value(JSONWR_T* jsonwr)
    cdef void jsonwr_end_array_value(JSONWR_T* jsonwr)
    cdef void jsonwr_null_value(JSONWR_T* jsonwr)
    cdef void jsonwr_str_value(JSONWR_T* jsonwr, const char* value)
    cdef void jsonwr_strn_value(JSONWR_T* jsonwr, const char* value, size_t len)
    cdef void jsonwr_nstr_value(JSONWR_T* jsonwr, const char* value)
    cdef void jsonwr_lng_value(JSONWR_T* jsonwr, long value)
    cdef void jsonwr_dbl_value(JSONWR_T* jsonwr, double value)
    cdef void jsonwr_log10num_value(JSONWR_T* jsonwr, double value, int prec)
    cdef void jsonwr_bool_value(JSONWR_T* jsonwr, int value)

    cdef void jsonwr_str_array_value(JSONWR_T* jsonwr, char** values, int count)
    cdef void jsonwr_nstr_array_value(JSONWR_T* jsonwr, char** values, int count)
    cdef void jsonwr_lng_array_value(JSONWR_T* jsonwr, long* values, int count)
    cdef void jsonwr_dbl_array_value(JSONWR_T* jsonwr, double* values, int count)
    cdef void jsonwr_log10num_array_value(JSONWR_T* jsonwr, double* values, int prec, int count);
    cdef void jsonwr_bool_array_value(JSONWR_T* jsonwr, int* values, int count)
    cdef void jsonwr_null_prop(JSONWR_T* jsonwr, const char* property)
    cdef void jsonwr_str_prop(JSONWR_T* jsonwr, const char* property, const char* value)
    cdef void jsonwr_strn_prop(JSONWR_T* jsonwr, const char* property, const char* value, size_t len)
    cdef void jsonwr_nstr_prop(JSONWR_T* jsonwr, const char* property, const char* value)
    cdef void jsonwr_lng_prop(JSONWR_T* jsonwr, const char* property, long value)
    cdef void jsonwr_dbl_prop(JSONWR_T* jsonwr, const char* property, double value)
    cdef void jsonwr_log10num_prop(JSONWR_T* jsonwr, const char* property, double value, int prec)
    cdef void jsonwr_bool_prop(JSONWR_T* jsonwr, const char* property, int value)
    cdef void jsonwr_str_array_prop(JSONWR_T* jsonwr, const char* property, char** values, int count)
    cdef void jsonwr_nstr_array_prop(JSONWR_T* jsonwr, const char* property, char** values, int count)
    cdef void jsonwr_lng_array_prop(JSONWR_T* jsonwr, const char* property, long* values, int count)
    cdef void jsonwr_dbl_array_prop(JSONWR_T* jsonwr, const char* property, double* values, int count)
    cdef void jsonwr_log10num_array_prop(JSONWR_T* jsonwr, const char* property, double* values, int prec, int count);
    cdef void jsonwr_bool_array_prop(JSONWR_T* jsonwr, const char* property, int* values, int count)
    cdef void jsonwr_args_prop(JSONWR_T* jsonwr, const char* property, int argc, char** argv)
    cdef void jsonwr_desc_prop(JSONWR_T* jsonwr, const char* property, const char* desc_file, const char* desc_message);
    cdef void jsonwr_time_prop(JSONWR_T* jsonwr, const char* property, time_t* when)
