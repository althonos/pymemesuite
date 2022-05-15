# coding: utf-8
# cython: language_level=3, linetrace=True
"""Internal API common to all MEME tools.
"""

from libmeme.alphabet cimport ALPH_T
from libmeme.array cimport ARRAY_T
from libmeme.hash_table cimport HASH_TABLE
from libmeme.data_types cimport WEIGHTS_T, Z_T, LCB_T
from libmeme.motif cimport MOTIF_T
from libmeme.meme cimport CANDIDATE, DATASET, MODEL, SAMPLE, THETA, MOTIF_SUMMARY
from libmeme.prior_dist cimport PRIOR_DIST_T
from libmeme.pssm cimport PSSM_T
from libmeme.seq cimport SEQ_T
from libmeme.user cimport MINSITES

# --- Python imports ---------------------------------------------------------

import os
import warnings

from .errors import AllocationError


# --- Utils ------------------------------------------------------------------

cdef void* matrix_create(size_t nrows, size_t ncols, size_t itemsize) except NULL
cdef void matrix_free(void** matrix)


# --- Alphabet ---------------------------------------------------------------

cdef class Alphabet:
    cdef ALPH_T* _alph


# --- Array ------------------------------------------------------------------

cdef class Array:
    cdef          ARRAY_T*   _array


# --- Candidate model --------------------------------------------------------

cdef class Candidate:
    cdef CANDIDATE* _cd


# --- Dataset ----------------------------------------------------------------

cdef class Dataset:
    cdef          DATASET* _ds
    cdef readonly Alphabet alphabet
    cdef          list     samples

    cdef void _allocate(self) except *
    cdef void _deallocate(self) except *
    cdef void _append(self, Sample sample) except *
    cdef void _set_database_freq_and_entropy(self)

    cpdef void shuffle(self) except *


# --- Model ------------------------------------------------------------------

cdef class Model:
    cdef          MODEL*   _mm
    cdef readonly object   _owner
    cdef readonly Alphabet alphabet

    # cdef void _allocate(self) except *
    # cdef void _deallocate(self) except *
    cpdef Model copy(self)


# --- Motif ------------------------------------------------------------------

cdef class Motif:
    cdef          MOTIF_T*   _motif
    cdef readonly Alphabet alphabet


# --- Motif summary ----------------------------------------------------------

cdef class MotifSummaries:
    cdef readonly Alphabet       alphabet
    cdef          MOTIF_SUMMARY* _motif_summaries
    cdef          size_t         _nmotifs

cdef class MotifSummary:
    cdef          MOTIF_SUMMARY* _ms
    cdef readonly object         _owner
    cdef readonly Alphabet       alphabet


# --- PriorDistribution ------------------------------------------------------

cdef class PriorDistribution:
    cdef          PRIOR_DIST_T* _pd


# --- PSSM -------------------------------------------------------------------

cdef class PSSM:
    cdef          PSSM_T* _pssm


# --- Sample -----------------------------------------------------------------

cdef class Sample:
    cdef          SAMPLE*  _sm
    cdef readonly object   _owner
    cdef readonly Alphabet alphabet
    cdef readonly bint     use_complement

    cdef void _allocate(self) except *
    cdef void _deallocate(self) except *
    cdef void _encode_sequence(self) except *
    cdef void _compute_sequence_weights(self) except *
    cdef void _count_residues(self) except *

    cpdef Sample copy(self)


# --- Sequence ---------------------------------------------------------------

cdef class Sequence:
  cdef          SEQ_T*  _seq
