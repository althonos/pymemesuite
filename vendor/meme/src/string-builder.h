/**************************************************************************
 * FILE: string-builder.h
 * AUTHOR: James Johnson 
 * CREATE DATE: 02 June 2011
 * PROJECT: shared
 * COPYRIGHT: UQ, 2011
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: Enables easier assembling of strings from multiple parts.
 **************************************************************************/

#ifndef STRING_BUILDER_H
#define STRING_BUILDER_H

#include <stdarg.h>
#include <stdint.h>

typedef struct str STR_T;

/*
 * str_create
 * Creates a string builder with the intial 
 * capacity as specified.
 */
STR_T* str_create(int capacity);

/*
 * str_destroy
 * Destroys the string builder. If keep_string
 * is specified it returns the string otherwise
 * it destroys the string and returns null.
 */
char* str_destroy(STR_T *strb, int keep_string);

/*
 * str_append
 * Appends a string to the string builder.
 */
void str_append(STR_T *strb, const char *str, int len);

/*
 * str_append
 * Appends a string to the string builder.
 */
void str_append2(STR_T *strb, const char *str);

/*
 * str_vappendf
 * Appends a formatted string to the string builder.
 */
void str_vappendf(STR_T *strb, const char *fmt, va_list ap);

/*
 * str_appendf
 * Appends a formatted string to the string builder.
 */
void str_appendf(STR_T *strb, const char *fmt, ...);

/*
 * str_append_code
 * Appends a Unicode code point in UTF-8 encoding
 */
void str_append_code(STR_T *strb, int32_t code_point);

/*
 * str_insert
 * Inserts a string at the offset to the string builder.
 */
void str_insert(STR_T *strb, int offset, const char *str, int len);

/*
 * str_insert2
 * Inserts a string at the offset to the string builder.
 */
void str_insert2(STR_T *strb, int offset, const char *str);

/*
 * str_vinsertf
 * Inserts a formatted string at the offset to the string builder.
 */
void str_vinsertf(STR_T *strb, int offset, const char *fmt, va_list ap); 

/*
 * str_insertf
 * Inserts a formatted string at the offset to the string builder.
 */
void str_insertf(STR_T *strb, int offset, const char *fmt, ...);

/*
 * str_insertf
 * Inserts a UTF-8 code unit at the offset to the string builder.
 */
void str_insert_code(STR_T *strb, int offset, int32_t code_point);

/*
 * str_replace
 * Replaces a string in the range specified with a passed string in 
 * the string builder.
 */
void str_replace(STR_T *strb, int start, int end, const char *str, int len);

/*
 * str_replace
 * Replaces a string in the range specified with a passed string in 
 * the string builder.
 */
void str_replace2(STR_T *strb, int start, int end, const char *str);

/*
 * str_vreplacef
 * Replaces a string in the range specified with a formatted string in 
 * the string builder.
 */
void str_vreplacef(STR_T *strb, int start, int end, const char *fmt, va_list ap);

/*
 * str_replacef
 * Replaces a string in the range specified with a formatted string in 
 * the string builder.
 */
void str_replacef(STR_T *strb, int start, int end, const char *fmt, ...);

/*
 * str_set
 * Set the stored string to the passed string in
 * the string builder.
 */
void str_set(STR_T *strb, const char *str, int len);

/*
 * str_vsetf
 * Set the stored string to the passed string in
 * the string builder.
 */
void str_vsetf(STR_T *strb, const char *fmt, va_list ap);

/*
 * str_setf
 * Set the stored string to the passed string in
 * the string builder.
 */
void str_setf(STR_T *strb, const char *fmt, ...);

/*
 * str_delete
 * Removes a string in the range specified from the string builder.
 */
void str_delete(STR_T *strb, int start, int end);

/*
 * str_truncate
 * Truncates the string in the string builder to the
 * specified length. If the value is negative
 * then the string will be shortened by that amount.
 */
void str_truncate(STR_T *strb, int length);

/*
 * str_clear
 * Truncates the string to length 0.
 */
void str_clear(STR_T *strb);

/*
 * str_len
 * Returns the length of the string under construction.
 */
size_t str_len(STR_T *strb);

/*
 * str_internal
 * Returns the internal string.
 *
 * !!WARNING!!
 * You should only use this in cases where the string is not going to be 
 * modified. Also be aware that when the string builder is destroyed the string 
 * will also be freed unless the destructor is specifically instructed not to.
 */
char* str_internal(STR_T *strb);

/*
 * str_char
 * Returns the character at the given position. Negative positions
 * are counted from the end so a pos of -1 returns the character 
 * before the null byte.
 */
char str_char(STR_T *strb, int pos);

/*
 * str_copy
 * Returns a copy of the internal string. The caller is responsible for
 * freeing the memory.
 */
char* str_copy(STR_T *strb);

/*
 * str_subcopy
 * Returns the specified substring copy of the internal string. The caller is
 * responsible for freeing the memory.
 */
char* str_subcopy(STR_T *strb, int start, int end);

/*
 * str_cmp
 * Applies strcmp to the strings.
 */
int str_cmp(const STR_T *strb, const char *s2);

/*
 * str_ncmp
 * Applies strncmp to the strings.
 */
int str_ncmp(const STR_T *strb, const char *s2, size_t n);

/*
 * str_casecmp
 * Applies strcasecmp to the strings.
 */
int str_casecmp(const STR_T *strb, const char *s2);

/*
 * str_fit
 * Shrinks the allocated memory to fit the string
 * exactly. Removes any minimum set by ensure.
 */
void str_fit(STR_T *strb);

/*
 * str_ensure
 * Ensures the requested capacity.
 */
void str_ensure(STR_T *strb, int new_capacity);

/*
 * str_append_path
 * Appends multiple path segments together separating them with / when the
 * segments do not already end with a /. If the last segment ends with a /
 * then it is removed, unless of course if the path is just /.
 *
 * Zero length segments are skipped.
 */
void str_append_path(STR_T *strb, int segments, ...);

/*
 * str_append_evalue
 * Appends a log10 evalue.
 */
void str_append_evalue(STR_T *strb, double log10_ev, int prec);

/*
 * str_evalue
 * sets the string to the log10 evalue and returns the string.
 */
char* str_evalue(STR_T *strb, double log10_ev, int prec);

#endif
