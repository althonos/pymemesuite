/* Copyright (c) 2007
 *   Minghui Jiang, James Anderson, Joel Gillespie, and Martin Mayne.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *      this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 * 3. The names of its contributors may not be used to endorse or promote
 *      products derived from this software without specific prior written
 *      permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 *      ushuffle.c - uShuffle library code
 *      Mon Apr 23 14:35:21 MDT 2007
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ushuffle.h"
#include "utils.h"

// set random function 

static randfunc_t randfunc = random;

void set_randfunc(randfunc_t func) {
  randfunc = func;
}

// global variables for the Euler algorithm 

static const char *s_ = NULL;
static long l_ = 0;
static long k_ = 0;

typedef struct vertex {
  long *indices;
  long n_indices;
  long i_indices;
  long intree, next;
  long i_sequence;
} vertex;

static vertex *vertices = NULL;
static long n_vertices;
static long *indices = NULL;
static long root;

// hashtable utility 

typedef struct hentry {
  struct hentry *next;
  long i_sequence;
  long i_vertices;
} hentry;

static hentry *entries = NULL;
static hentry **htable = NULL;
static long htablesize = 0;
static double hmagic;

static long hcode(long i_sequence) {
  double f = 0.0;
  long i;

  for (i = 0; i < k_ - 1; i++) {
    f += s_[i_sequence + i];
    f *= hmagic;
  }
  if (f < 0.0)
    f = -f;
  return (long) (htablesize * f) % htablesize;
}

static void hinit(long size) {
  entries = mm_calloc(size, sizeof(hentry));
  htable = mm_calloc(size, sizeof(hentry *));
  htablesize = size;
  hmagic = (sqrt(5.0) - 1.0) / 2.0;
}

static void hcleanup() {
  free(entries);
  entries = NULL;
  free(htable);
  htable = NULL;
  htablesize = 0;
}

static void hinsert(long i_sequence) {
  long code = hcode(i_sequence);
  hentry *e, *e2 = &entries[i_sequence];

  for (e = htable[code]; e; e = e->next)
    if (strncmp(&s_[e->i_sequence], &s_[i_sequence], k_ - 1) == 0) {
      e2->i_sequence = e->i_sequence;
      e2->i_vertices = e->i_vertices;
      return;
    }
  e2->i_sequence = i_sequence;
  e2->i_vertices = n_vertices++;
  e2->next = htable[code];
  htable[code] = e2;
}

// the Euler algorithm 

void ushuffle1(const char *s, long l, long k) {
  long i, j, n_lets;

  s_ = s;
  l_ = l;
  k_ = k;
  if (k_ >= l_ || k_ <= 1) // two special cases 
    return;

  // use hashtable to find distinct vertices 
  n_lets = l_ - k_ + 2; // number of (k-1)-lets 
  n_vertices = 0;
  hinit(n_lets);
  for (i = 0; i < n_lets; i++)
    hinsert(i);
  root = entries[n_lets - 1].i_vertices; // the last let 
  if (vertices)
    free(vertices);
  vertices = mm_calloc(n_vertices, sizeof(vertex));

  // set i_sequence and n_indices for each vertex 
  for (i = 0; i < n_lets; i++) { // for each let 
    hentry *ev = &entries[i];
    vertex *v = &vertices[ev->i_vertices];

    v->i_sequence = ev->i_sequence;
    if (i < n_lets - 1) // not the last let 
      v->n_indices++;
  }

  // distribute indices for each vertex 
  if (indices)
    free(indices);
  indices = mm_calloc((n_lets - 1), sizeof(long));
  j = 0;
  for (i = 0; i < n_vertices; i++) { // for each vertex 
    vertex *v = &vertices[i];

    v->indices = indices + j;
    j += v->n_indices;
  }

  // populate indices for each vertex 
  for (i = 0; i < n_lets - 1; i++) { // for each edge 
    hentry *eu = &entries[i];
    hentry *ev = &entries[i + 1];
    vertex *u = &vertices[eu->i_vertices];

    u->indices[u->i_indices++] = ev->i_vertices;
  }
  hcleanup();
}

void permutec(char *t, long l) {
  long i, j;
  char tmp;

  for (i = l - 1; i > 0; i--) {
    j = (*randfunc)() % (i + 1);
    tmp = t[i]; t[i] = t[j]; t[j] = tmp; // swap 
  }
}

static void permutei(long *t, long l) {
  long i, j;
  long tmp;

  for (i = l - 1; i > 0; i--) {
    j = (*randfunc)() % (i + 1);
    tmp = t[i]; t[i] = t[j]; t[j] = tmp; // swap 
  }
}

void ushuffle2(char *t) {
  vertex *u, *v;
  long i, j;

  // exact copy case 
  if (k_ >= l_) {
    strncpy(t, s_, l_);
    return;
  }

  // simple permutation case 
  if (k_ <= 1) {
    strncpy(t, s_, l_);
    permutec(t, l_);
    return;
  }

  // the Wilson algorithm for random arborescence 
  for (i = 0; i < n_vertices; i++)
    vertices[i].intree = 0;
  vertices[root].intree = 1;
  for (i = 0; i < n_vertices; i++) {
    u = &vertices[i];
    while (!u->intree) {
      u->next = (*randfunc)() % u->n_indices;
      u = &vertices[u->indices[u->next]];
    }
    u = &vertices[i];
    while (!u->intree) {
      u->intree = 1;
      u = &vertices[u->indices[u->next]];
    }
  }

  // ushuffle indices to prepare for walk 
  for (i = 0; i < n_vertices; i++) {
    u = &vertices[i];
    if (i != root) {
      j = u->indices[u->n_indices - 1]; // swap the last one 
      u->indices[u->n_indices - 1] = u->indices[u->next];
      u->indices[u->next] = j;
      permutei(u->indices, u->n_indices - 1); // permute the rest 
    } else
      permutei(u->indices, u->n_indices);
    u->i_indices = 0; // reset to zero before walk 
  }

  // walk the graph 
  strncpy(t, s_, k_ - 1); // the first let remains the same 
  u = &vertices[0];
  i = k_ - 1;
  while (u->i_indices < u->n_indices) {
    v = &vertices[u->indices[u->i_indices]];
    j = v->i_sequence + k_ - 2;
    t[i++] = s_[j];
    u->i_indices++;
    u = v;
  }
}

void ushuffle(const char *s, char *t, long l, long k) {
  ushuffle1(s, l, k);
  ushuffle2(t);
}
