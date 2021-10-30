/**************************************************************************
 * FILE: linked-list.h
 * AUTHOR: James Johnson 
 * CREATE DATE: 17-August-2009 
 * PROJECT: shared
 * COPYRIGHT: TBA 
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: Data structure for manipulating a list
 **************************************************************************/

#ifndef LINKED_LIST_H
#define LINKED_LIST_H

/*
 * Types
 */
typedef struct linklst_t LINKLST_T;
typedef struct link_t LL_LINK_T;


/*
 * linklst_create
 * creates and returns an empty linked list
 */
LINKLST_T *linklst_create();

/*
 * linklst_destroy_all
 * destroys a linked list using the passed function 
 * to deallocate all contained items
 */
void linklst_destroy_all(LINKLST_T *linklst, void(*free_item)(void *));

/*
 * linklst_destroy
 * destroys a linked list
 * caller must deallocate items themselves
 */
void linklst_destroy(LINKLST_T *linklst);

/*
 * Return a copy of a linked list.
 */
LINKLST_T *linklst_copy(LINKLST_T *linklst);

/*
 * linklst_size
 * returns the number of items in the list
 */
int linklst_size(LINKLST_T *linklst);

/*
 * linklst_first
 * gets the first link
 */
LL_LINK_T *linklst_first(LINKLST_T *linklst);

/*
 * linklst_last
 * gets the last link
 */
LL_LINK_T *linklst_last(LINKLST_T *linklst);

/*
 * linklst_add_after
 * adds an item after the specified link
 */
LL_LINK_T *linklst_add_after(void *item, LL_LINK_T *before, LINKLST_T *linklst);

/*
 * linklst_add_before
 * adds an item before the specified link
 */
LL_LINK_T *linklst_add_before(void *item, LL_LINK_T *after, LINKLST_T *linklst);

/*
 * linklst_add
 * adds an item at the end of the list
 */
LL_LINK_T *linklst_add(void *item, LINKLST_T *linklst);

/*
 * Add the elements in the second list to the first list
 * and return first list.
 */
LINKLST_T *linklst_plus_equals(LINKLST_T *linklst, LINKLST_T *linklst2);

/*
 * Return a new list consisting of the elements
 * of the first list followed by those of the second list.
 */
LINKLST_T *linklst_plus(LINKLST_T *linklst, LINKLST_T *linklst2);

/*
 * linklst_take
 * removes an item from the end of the list
 * and returns a pointer to the item
 */
void *linklst_take(LINKLST_T *linklst);

/*
 * linklst_push
 * adds an item to the front of the list
 */
LL_LINK_T *linklst_push(void *item, LINKLST_T *linklst);

/*
 * linklst_pop
 * removes an item from the front of the list
 * and returns a pointer to the item
 */
void *linklst_pop(LINKLST_T *linklst);

/*
 * Returns the item at the front of the linked list
 */
void *linklst_peek(LINKLST_T *linklst);

/*
 * linklst_remove
 * removes the passed link from the list
 * and returns a pointer to the item.
 */
void *linklst_remove(LL_LINK_T *link, LINKLST_T *linklst);

/*
 * linklst_sort
 * sorts the list using the passed
 * comparator to compare items.
 */
void linklst_sort(int (*comparator)(void*,void*), LINKLST_T *linklst);

/*
 * linklst_next
 * gets the next link
 */
LL_LINK_T *linklst_next(LL_LINK_T *link);

/*
 * linklst_prev
 * gets the previous link
 */
LL_LINK_T *linklst_prev(LL_LINK_T *link);

/*
 * linklst_get
 * gets the item contained in the link
 */
void *linklst_get(LL_LINK_T *link);

/*
 * linklst_set
 * sets the item contained in the link
 */
void linklst_set(void *item, LL_LINK_T *link);

#endif

