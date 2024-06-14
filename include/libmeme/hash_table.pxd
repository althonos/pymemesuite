from libmeme.string_list cimport STRING_LIST_T


cdef extern from "hash_table.h" nogil:

    const size_t MAX_BINS

    cdef struct hash_table_entry:
        pass
    ctypedef hash_table_entry HASH_TABLE_ENTRY

    cdef struct hash_table:
        pass
    ctypedef hash_table* HASH_TABLE

    cdef int get_num_entries(HASH_TABLE ht)
    cdef HASH_TABLE hash_create(int n, void (*free_value)(void*))
    cdef void hash_destroy(void* ht)
    cdef void hash_entry_destroy(HASH_TABLE_ENTRY* hte, void (*free_value)(void*))
    cdef bint hash_insert_str(char* key, HASH_TABLE ht)
    cdef void* hash_get_entry_value(HASH_TABLE_ENTRY* hte)
    cdef char* hash_get_entry_key(HASH_TABLE_ENTRY* hte)
    cdef void hash_set_entry_value(void* value, HASH_TABLE_ENTRY* hte)
    cdef bint hash_insert(char* key1, int key2, HASH_TABLE ht)
    cdef bint hash_insert_str_value(char* key1, void* value, HASH_TABLE ht)
    cdef bint hash_insert_value(char* key1, int key2, void* value, HASH_TABLE ht)
    cdef bint hash_remove_str(char* key1, HASH_TABLE ht)
    cdef bint hash_remove(char* key1, int key2, HASH_TABLE ht)
    cdef HASH_TABLE_ENTRY* hash_lookup(char* key1, int key2, HASH_TABLE ht)
    cdef HASH_TABLE_ENTRY* hash_lookup_str(char* key1, HASH_TABLE ht)
    cdef STRING_LIST_T* hash_get_keys(HASH_TABLE ht)
