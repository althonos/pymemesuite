/**
 * @file hash_table.c
 *
 * This module implements the hash-table datatype.
 *
 * $Id: hash_table.c 1589 2007-02-06 04:17:07Z tom $
 * $Log$
 * Revision 1.1.2.5  2006/01/13 02:41:39  twhitington
 * Added CVS Log: line to header.
 *
*/

/* hash_table.c */
//
// HASH TABLE function definitions
//
#include "hash_table.h"
#include <assert.h>

#define two24 	8388610			/* 2**24 */
#define DUMMY_INT 0                     /* The dummy value */

// Structures
struct hash_table_entry {
  char *key1;			// character key
  int key2;			// integer key
  int hash_value;		// hash value of <keq1, key2>
  void* value;			// the value assigned to a key (optional)
  HASH_TABLE_ENTRY *next;	// link to collision list
  HASH_TABLE_ENTRY *prev;	// backward link in collision list
};

struct hash_table {
  int n;			/* number of bins in hash table */
  int n_entries;                /* total number of entries in the hash table */
  HASH_TABLE_ENTRY **table;	/* array of entry pointers */
  void (*free_value)(void *);	// function for freeing entry value;
				// caller must free them if NULL
};

// Functions

// Prototype for internal hash lookup function
HASH_TABLE_ENTRY *hash_lookup_internal(char *key1, int key2, HASH_TABLE ht, int *hash_value);

/**********************************************************************/
/*
	hash

	Hashing function for <string, int> items.
*/
/**********************************************************************/
static int hash_keys(
  char *key1, 				/* character key */
  int key2,				/* integer key */
  int n					/* modulo */
)
{
  int i, p, d;

  if (!key1) return key2 % n;		// no key1

  for (i=0, d=1, p=key2; key1[i] != 0; i++, d *= 256) {
    if (d >= two24) d = 1;		/* d<2**23 insures d*key1<2**31 */
    p += (d * key1[i]) % n;
  }
  return p % n;
}


/**********************************************************************/
/*
	get_num_entries

	Get the number of entries in the hash table.
*/
/**********************************************************************/
int get_num_entries(
  HASH_TABLE ht
)
{
  return ht->n_entries;
}


/**********************************************************************/
/*
	hash_create

	Create a hash table.
*/
/**********************************************************************/
HASH_TABLE hash_create(
  int n,
  void (*free_value)(void *)		// function for freeing value;
					// caller must free them if NULL
)
{
  // PRECONDITION: Caller must specify a number of bins >= 1:
  assert(n>0);

  int i;
  HASH_TABLE ht = (HASH_TABLE) mm_malloc(sizeof(struct hash_table));

  /* initialize hash table */
  ht->n = n;
  ht->table = (HASH_TABLE_ENTRY **) mm_malloc(n * sizeof(HASH_TABLE_ENTRY *));
  ht->n_entries = 0;
  for (i=0; i<n; i++) ht->table[i] = NULL;
  ht->free_value = free_value;

  return ht;
}

/**********************************************************************/
/*
	hash_destroy

	Destroy a hash table.
*/
/**********************************************************************/
void hash_destroy(
  void *htv 				/* hash table to destroy */
)
{
  int i;
  HASH_TABLE ht = htv;
  for (i=0; i < ht->n; i++) {
    HASH_TABLE_ENTRY *hte = ht->table[i];
    while (hte != NULL) {
      HASH_TABLE_ENTRY *next = hte->next;
      hash_entry_destroy(hte, ht->free_value);
      hte = next;
    }
  }
  myfree(ht->table);
  myfree(ht);
} // hash_destroy

/**********************************************************************/
/*
	hash_entry_destroy

	Destroy a hash table entry.

*/
/**********************************************************************/

static void hash_entry_destroy(
  HASH_TABLE_ENTRY *hte, 	/* hash table entry to destroy */
  void (*free_value)(void *) 	// function for freeing entry value;
				// caller must free them if NULL
)
{
  myfree(hte->key1);
  if (free_value && hte->value) free_value(hte->value);
  myfree(hte);
} // hash_entry_destroy

/**********************************************************************/
/*
	hash_insert_str

	Insert a <string, int> item into a hash table, where the int
        is equal to DUMMY_INT. This method allows the user to just specify
        a string, and yet to be sure that future hash table queries with
        that string will correctly recognise the prior presence of the
        string.

	Returns true if successful.
	Returns false if the item already is present.
*/
/**********************************************************************/
bool hash_insert_str(
  char *key1,			/* character key */
  HASH_TABLE ht			/* the hash table */
) {
  return hash_insert_value(key1, DUMMY_INT, NULL, ht);
}

/**********************************************************************/
/*
	hash_set_entry_value

	Adds a value (void*) to the specified entry.

	Attention: former values are not freed by this method but must be
	freed by the user
*/
/**********************************************************************/
void hash_set_entry_value(
  void* value,			/* the value */
  HASH_TABLE_ENTRY *hte         /* hash table entry to equip with an value */
){
  assert(hte != NULL);
  hte->value = value;
}

/**********************************************************************/
/*
	hash_get_entry_value

	Returns the value attached to a given hash entry
*/
/**********************************************************************/
void* hash_get_entry_value(
  HASH_TABLE_ENTRY *hte         /* hash table entry of interest */
){
  assert(hte != NULL);
  return (hte->value);
}

/**********************************************************************/
/*
	hash_get_entry_key

	Returns the key attached to a given hash entry
*/
/**********************************************************************/
char* hash_get_entry_key(
  HASH_TABLE_ENTRY *hte         /* hash table entry of interest */
){
  assert(hte != NULL);
  return (hte->key1);
}


/**********************************************************************/
/*
	hash_insert_str_value

	Insert a <string, int> item into a hash table, where the int
        is equal to DUMMY_INT. This method allows the user to just specify
        a string, and yet to be sure that future hash table queries with
        that string will correctly recognise the prior presence of the
        string.

	Returns true if successful.
	Returns false if the item already is present.
*/
/**********************************************************************/
bool hash_insert_str_value(
  char *key1,			/* character key */
  void *value,			/* the value */
  HASH_TABLE ht			/* the hash table */
) {
  return hash_insert_value(key1, DUMMY_INT, value, ht);
}

/**********************************************************************/
/*
	hash_insert

	Insert a <string, int> item into a hash table.

	Returns true if successful.
	Returns false if the item already is present.
*/
/**********************************************************************/
bool hash_insert(
  char *key1,			/* character key */
  int key2,			/* integer key */
  HASH_TABLE ht			/* the hash table */
)
{
	return hash_insert_value(key1, key2, NULL ,ht);
} // hash_insert


/**********************************************************************/
/*
	hash_insert_value

	Insert a <string, int> item with a value into a hash table.

	Returns true if successful.
	Returns false if the item already is present.
*/
/**********************************************************************/
bool hash_insert_value(
  char *key1,			/* character key */
  int key2,			/* integer key */
  void *value,			/* the value */
  HASH_TABLE ht			/* the hash table */
)
{
  HASH_TABLE_ENTRY *hte;
  int hash_value;

  if ((hte = hash_lookup_internal(key1, key2, ht, &hash_value)) == NULL) {
    /* make a hash-table entry */
    hte = (HASH_TABLE_ENTRY *) mm_malloc(sizeof(HASH_TABLE_ENTRY));
    if (key1) {
      hte->key1 = (char *) mm_malloc (strlen(key1) + 1);
      strcpy(hte->key1, key1);
    } else {
      hte->key1 = NULL;
    }
    hte->key2 = key2;
    hte->hash_value = hash_value;
    hte->value = value;

    // Keep count of the number of entries in the hash table:
    ht->n_entries++;

    /* insert the entry into the table at the head of the collision list */
    hte->next = ht->table[hash_value];	// insert at head; next may be NULL
    hte->prev = NULL;			// new head of list; no previous node
    // If there was a collision, the head of the collision list points back
    // to the new entry.
    if (ht->table[hash_value] != NULL) ht->table[hash_value]->prev = hte;
    ht->table[hash_value] = hte;	// now put at head of list
    return(true);
  } else {
    return(false);
  }
} // hash_insert

/**********************************************************************/
/*
	hash_remove

	Remove a <string, int> item from a hash table.

	ATTENTION: If a value was passed to the hash_table the user has to
	ensure that the memory for this value is/will be freed by himself!

	Returns true if successful.
	Returns false if the item was not present.
*/
/**********************************************************************/
bool hash_remove(
  char *key1,			/* character key */
  int key2,			/* integer key */
  HASH_TABLE ht			/* the hash table */
)
{
  HASH_TABLE_ENTRY *hte;
  int hash_value;

  if ((hte = hash_lookup_internal(key1, key2, ht, &hash_value)) != NULL) {

    // delete the entry from the table
    if (hte->prev == NULL) {		// at head of list
      ht->table[hash_value] = hte->next;
      if (hte->next != NULL) (hte->next)->prev = NULL;	// not also at tail
    } else {				// not at head of list
      (hte->prev)->next = hte->next;
      if (hte->next != NULL) {		// in middle of list
        (hte->next)->prev = hte->prev;
      }
    }
    hash_entry_destroy(hte, ht->free_value);	// destroy the entry

    // Keep count of the number of entries in the hash table:
    ht->n_entries--;

    return(true);
  } else {
    return(false);
  }
} // hash_insert

/**********************************************************************/
/*
	hash_remove_str

	Remove a <string, int> item from a hash table, where the int
        is equal to DUMMY_INT. This method allows the user to just specify
        a string, and yet to be sure that future hash table queries with
        that string will correctly recognise the abscence of the
        string.

	ATTENTION: If a value was passed to the hash_table the user has to
	ensure that the memory for this value is/will be freed by himself!

	Returns true if successful.
	Returns false if the item was not present.
*/
/**********************************************************************/
bool hash_remove_str(
  char *key1,			/* character key */
  HASH_TABLE ht			/* the hash table */
) {
  return hash_remove(key1, DUMMY_INT, ht);
}

/**********************************************************************/
/*
	hash_lookup_internal

	Look up a <string, int> item in a hash table.

	Returns hash value of item in hash_value.

	Returns ptr to entry if the item is present, NULL otherwise.
*/
/**********************************************************************/
HASH_TABLE_ENTRY *hash_lookup_internal(
  char *key1,			/* character key */
  int key2,			/* integer key */
  HASH_TABLE ht,		/* the hash table */
  int *hash_value		// hash value of key
)
{
  /* compute the position in hash table of entry */
  *hash_value = hash_keys(key1, key2, ht->n);
  /* get collision list for this hash index */
  HASH_TABLE_ENTRY *hte = ht->table[*hash_value];
  while (hte != NULL) {
    //fprintf(stderr, "looking up %s\n", key1);
    if ((hte->key2 == key2) && 		// integer keys match AND
      ((!key1 && !hte->key1) || 	// string keys are NULL OR
      (key1 && hte->key1 && !strcmp(hte->key1, key1)))) // string keys match
        return hte;
    hte = hte->next;
  }
  return NULL;
} // hash_lookup_internal

/**********************************************************************/
/*
	hash_lookup

	Look up a <string, int> item in a hash table.

	Returns ptr to entry if the item is present, NULL otherwise.
*/
/**********************************************************************/
HASH_TABLE_ENTRY *hash_lookup(
  char *key1,			/* character key */
  int key2,			/* integer key */
  HASH_TABLE ht			/* the hash table */
)
{
  int hash;
  return hash_lookup_internal(key1, key2, ht, &hash);
} // hash_lookup

/**********************************************************************/
/*
	hash_lookup_str

	Look up a <string, int> item in a hash table, where the int
        is equal to DUMMY_INT. The user can be sure that the presence
        of a <string> item will be detected, as long as they have only
        used "hash_insert_str" to insert strings.

	Returns ptr to entry if the item is present, NULL otherwise.
*/
/**********************************************************************/
HASH_TABLE_ENTRY *hash_lookup_str(
  char *key1,			/* character key */
  HASH_TABLE ht			/* the hash table */
) {
  int hash;
  return hash_lookup_internal(key1, DUMMY_INT, ht, &hash);
}

/**********************************************************************/
/*
	hash_get_keys

	Return the (string) keys of a hash table.
*/
/**********************************************************************/
STRING_LIST_T *hash_get_keys(
  HASH_TABLE ht			// the hash table
) {

  int i;
  STRING_LIST_T *keys = new_string_list();

  for (i=0; i < ht->n; i++) {
    HASH_TABLE_ENTRY *hte = ht->table[i];
    while (hte != NULL) {
      add_string(hte->key1, keys);
      HASH_TABLE_ENTRY *next = hte->next;
      hte = next;
    }
  }

  return(keys);

} // hash_get_keys
