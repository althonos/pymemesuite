/**
 * @file hash_table.h
 *
 * This module implements the hash table datatype.
 *
 * $Id$
*/

// hash_table.h
//
//
// HASH TABLE object definitions
//
#ifndef HASH_H
#define HASH_H
#include "macros.h"
#include "string-list.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif

// Data structures

#define MAX_BINS 10000000 // 10 million

/* hash table entry */
typedef struct hash_table_entry HASH_TABLE_ENTRY;

/* hash table */
typedef struct hash_table *HASH_TABLE;

// External Functions
extern int get_num_entries(
  HASH_TABLE ht
);

HASH_TABLE hash_create(
  int n,
  void (*free_value)(void *)		// function for freeing value;
					// caller must free them if NULL
);

void hash_destroy(
  void *ht			/* hash table to destroy */
);

static void hash_entry_destroy(
  HASH_TABLE_ENTRY *hte,	/* hash table entry to destroy */
  void (*free_value)(void *)		// function for freeing value;
					// caller must free them if NULL
);

bool hash_insert_str(
  char *key1,			/* character key */
  HASH_TABLE ht			/* the hash table */
);

void* hash_get_entry_value(
  HASH_TABLE_ENTRY *hte         /* hash table entry of interest */
);

char* hash_get_entry_key(
  HASH_TABLE_ENTRY *hte         /* hash table entry of interest */
);

void hash_set_entry_value(
  void* value,			/* the value */
  HASH_TABLE_ENTRY *hte         /* hash table entry to equip with an value */
);

bool hash_insert(
  char *key1,    		/* character key */
  int key2,			/* integer key */
  HASH_TABLE ht         	/* the hash table */
);

bool hash_insert_str_value(
  char *key1,			/* character key */
  void *value,			/* the value */
  HASH_TABLE ht			/* the hash table */
);

bool hash_insert_value(
  char *key1,    		/* character key */
  int key2,			/* integer key */
  void *value,			/* the value */
  HASH_TABLE ht         	/* the hash table */
);

bool hash_remove_str(
  char *key1,			/* character key */
  HASH_TABLE ht			/* the hash table */
);

bool hash_remove(
  char *key1,                   /* character key */
  int key2,                     /* integer key */
  HASH_TABLE ht                 /* the hash table */
);

HASH_TABLE_ENTRY *hash_lookup(
  char *key1,			/* character key */
  int key2, 			/* integer key */
  HASH_TABLE ht         	/* the hash table */
);

HASH_TABLE_ENTRY *hash_lookup_str(
  char *key1,			/* character key */
  HASH_TABLE ht			/* the hash table */
);

STRING_LIST_T *hash_get_keys(
  HASH_TABLE ht                 // the hash table
); 

#endif
