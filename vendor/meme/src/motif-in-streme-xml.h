#ifndef MOTIF_IN_STREME_XML_H
#define MOTIF_IN_STREME_XML_H

#include "motif.h"
#include "motif-in-common.h"

void* sxml_create(const char *optional_filename, int options);
void sxml_destroy(void *data);
void sxml_update(void *data, const char *chunk, size_t size, short end);
short sxml_has_motif(void *data);
short sxml_has_format_match(void *data);
short sxml_has_warning(void *data);
char* sxml_next_warning(void *data);
short sxml_has_error(void *data);
char* sxml_next_error(void *data);
MOTIF_T* sxml_next_motif(void *data);
ALPH_T* sxml_get_alphabet(void *data);
int sxml_get_strands(void *data);
bool sxml_get_bg(void *data, ARRAY_T **bg);
void* sxml_motif_optional(void *data, int option);
void* sxml_file_optional(void *data, int option);

#endif


