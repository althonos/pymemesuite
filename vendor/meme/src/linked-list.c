/**************************************************************************
 * FILE: linked-list.c
 * AUTHOR: James Johnson 
 * CREATE DATE: 17-August-2009 
 * PROJECT: shared
 * COPYRIGHT: TBA 
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: Data structure for manipulating a list
 **************************************************************************/
#include <stdlib.h>
#include <strings.h>
#include "linked-list.h"
#include "utils.h"


/*
 * This is a link in a linked list.
 * The value 'before' will always point to
 * the previous link in the list or NULL if
 * it is the first link. The value 'after'
 * will always point to the next link in 
 * the list or NULL if it is the last link.
 */
struct link_t {
  LL_LINK_T *before;
  LL_LINK_T *after;
  void *item;
};

/*
 * This is the controling structure for the
 * linked list. The value 'size' holds the
 * number of links in the list. If 'size' is
 * zero then the 'first' and 'last' values will
 * be NULL, otherwise they will point to the
 * first and last links in the list. Size will 
 * never be negative.
 */
struct linklst_t {
  int size;
  LL_LINK_T *first;
  LL_LINK_T *last;
};

/*
 * Creates a link
 */
LL_LINK_T *link_create(LL_LINK_T *before, LL_LINK_T *after, void *item) {
  LL_LINK_T *link = (LL_LINK_T*)mm_malloc(sizeof(LL_LINK_T));
  link->before = before;
  link->after = after;
  link->item = item;
  if (before != NULL) before->after = link;
  if (after != NULL) after->before = link;
  return link;
}

/*
 * Destroys a link
 */
void link_destroy(LL_LINK_T *link) {
    memset(link, 0, sizeof(LL_LINK_T));
    free(link);
}

/*
 * Creates a linked list
 */
LINKLST_T *linklst_create() {
  LINKLST_T *linklst = (LINKLST_T*)mm_malloc(sizeof(LINKLST_T));
  linklst->size = 0;
  linklst->first = NULL;
  linklst->last = NULL;
  return linklst;
}

/*
 * Destroys a linked list, using the passed function to free the items
 */
void linklst_destroy_all(LINKLST_T *linklst, void(*free_item)(void *)) {
  LL_LINK_T *current = linklst->first; 
  LL_LINK_T *next;
  while (current != NULL) {
    next = current->after; 
    free_item(current->item);
    link_destroy(current);
    current = next;
  }
  memset(linklst, 0, sizeof(LINKLST_T));
  free(linklst);
}

/*
 * Literally do nothing. Used by the linklst_destroy function as an argument
 * to linklst_destroy_all
 */
static void do_nothing(void *ptr) {}

/*
 * Destroys a linked list but does not free the items
 */
void linklst_destroy(LINKLST_T * linklst) {
  linklst_destroy_all(linklst, do_nothing);
}

/*
 * Gets the size of the list
 */
int linklst_size(LINKLST_T *linklst) {
  return linklst->size;
}

/*
 * Gets the first link in the list
 */
LL_LINK_T *linklst_first(LINKLST_T *linklst) {
  return linklst->first;
}

/*
 * Gets the last link in the list
 */
LL_LINK_T *linklst_last(LINKLST_T *linklst) {
  return linklst->last;
}

/*
 * Add the item to the list after the passed link
 * and return the new link.
 * Note that the link must have come from this linked
 * list or bad things will happen!
 */
LL_LINK_T *linklst_add_after(void *item, LL_LINK_T *before, LINKLST_T *linklst) {
  LL_LINK_T *after = before->after;
  LL_LINK_T *link =  link_create(before, after, item); 
  if (after == NULL) linklst->last = link;
  linklst->size += 1;
  return link;
}

/*
 * Add the item to the list before the passed link
 * and return the new link.
 * Note that the link must have come from this linked
 * list or bad things will happen!
 */
LL_LINK_T *linklst_add_before(void *item, LL_LINK_T *after, LINKLST_T *linklst) {
  LL_LINK_T *before = after->before;
  LL_LINK_T *link = link_create(before, after, item);
  if (before == NULL) linklst->first = link;
  linklst->size += 1;
  return link;
}

/*
 * Add the item to the end of the linked list
 * and return the new link.
 */
LL_LINK_T *linklst_add(void *item, LINKLST_T *linklst) {
  LL_LINK_T *before = linklst->last;
  LL_LINK_T *after = NULL;
  LL_LINK_T *link = link_create(before, after, item);
  if (before == NULL) linklst->first = link;
  linklst->last = link;
  linklst->size += 1;
  return link;
}

/*
 * Return a copy of a linked list.
 */
LINKLST_T *linklst_copy(LINKLST_T *linklst) {
  LINKLST_T *new_list = linklst_create();
  LL_LINK_T *link = linklst->first;
  while (link) {
    linklst_add(link->item, new_list);
    link = link->after;
  }
  return new_list;
}

/*
 * Add the elements in the second list to the end of the
 * first list. Returns the first list.
 */
LINKLST_T *linklst_plus_equals(LINKLST_T *linklst, LINKLST_T *linklst2) {
  LL_LINK_T *link, *before, *after, *add_link;
  add_link = linklst2->first;
  while (add_link) {
    before = linklst->last;
    after = NULL;
    link = link_create(before, after, add_link->item);
    if (before == NULL) linklst->first = link;
    linklst->last = link;
    linklst->size += 1;
    add_link = add_link->after;
  }
  return linklst;
}

/*
 * Return a new list consisting of the elements 
 * of the first list followed by those of the second list.
 */
LINKLST_T *linklst_plus(LINKLST_T *linklst, LINKLST_T *linklst2) {
  return linklst_plus_equals(linklst_copy(linklst), linklst2);
}

/*
 * Takes a link from the end of the linked list
 * and returns the item
 */
void *linklst_take(LINKLST_T *linklst) {
  if (linklst->size == 0) return NULL;
  LL_LINK_T *link = linklst->last;
  LL_LINK_T *new_last = link->before;
  if (new_last != NULL) new_last->after = NULL;
  linklst->last = new_last;
  linklst->size -= 1;
  void *item = link->item;
  link_destroy(link);
  return item;
}

/*
 * Pushes the item on the front of the linked list
 * and returns the new link.
 */
LL_LINK_T *linklst_push(void *item, LINKLST_T *linklst) {
  LL_LINK_T *after = linklst->first;
  LL_LINK_T *before = NULL;
  LL_LINK_T *link = link_create(before, after, item);
  linklst->first = link;
  if (after == NULL) linklst->last = link;
  linklst->size += 1;
  return link;
}

/*
 * Pops an item off the front of the linked list
 * and returns the item
 */
void *linklst_pop(LINKLST_T *linklst) {
  if (linklst->size == 0) return NULL;
  LL_LINK_T *link = linklst->first;
  LL_LINK_T *new_first = link->after;
  if (new_first != NULL) new_first->before = NULL;
  else linklst->last = NULL;
  linklst->first = new_first;
  linklst->size -= 1;
  void *item = link->item;
  link_destroy(link);
  return item;
}

/*
 * Returns the item at the front of the linked list
 */
void *linklst_peek(LINKLST_T *linklst) {
  if (linklst->size == 0) return NULL;
  return linklst->first->item;
}

/*
 * Removes the passed link from the list and 
 * returns the item.
 * Note that the link must be a part of the
 * list or bad things will happen!
 */
void *linklst_remove(LL_LINK_T *link, LINKLST_T *linklst) {
  LL_LINK_T *before = link->before;
  LL_LINK_T *after = link->after;
  if (before != NULL) before->after = after;
  if (after != NULL) after->before = before;
  if (link == linklst->first) linklst->first = after;
  if (link == linklst->last) linklst->last = before;
  linklst->size -= 1;
  void *item = link->item;
  link_destroy(link);
  return item;
}

/*
 * Does a merge sort on the list
 *
 */
void linklst_sort(int (*comparator)(void*, void*), LINKLST_T *linklst) {
  LL_LINK_T *p, *q, *temp;
  int k, i, psize, qsize;
  if (linklst->size <= 1) return;
  k = 1;
  p = linklst->first;
  while (true) {
    q = p;
    for (i = 0; q != NULL && i < k; ++i) q = q->after; 
    psize = i;
    qsize = k;
    //merge list p of size psize and list q of maximum size qsize
    while (psize > 0 && qsize > 0 && q != NULL) {
      if (comparator(q->item, p->item) < 0) {
        //move q to be infront of p
        //save the item after q so we can move to it
        temp = q->after;
        //there is guaranteed to be something after p and something before q
        q->before->after = q->after;
        if (q->after) q->after->before = q->before;
        else linklst->last = q->before;
        q->after = p;
        q->before = p->before;
        if (p->before) p->before->after = q;
        else linklst->first = q;
        p->before = q;
        // move the pointer along the q list
        q = temp;
        qsize--;
      } else {
        //leave the item at p
        p = p->after;
        psize--;
      }
    }
    //skip over any already sorted stuff on the end of the q list
    for (; qsize > 0 && q != NULL; --qsize) q = q->after;
    //check if we've done a full pass
    if (q == NULL) {
      if (k >= linklst->size) return;
      k *= 2;
      p = linklst->first;
    } else {
      p = q;
    }
  }
}

/*
 * Gets the next link
 */
LL_LINK_T *linklst_next(LL_LINK_T *link) {
  return link->after;
}

/*
 * Gets the previous link
 */
LL_LINK_T *linklst_prev(LL_LINK_T *link) {
  return link->before;
}

/*
 * Gets the item
 */
void *linklst_get(LL_LINK_T *link) {
  return link->item;
}

/*
 * Sets the item
 */
void linklst_set(void *item, LL_LINK_T *link) {
  link->item = item;
}

