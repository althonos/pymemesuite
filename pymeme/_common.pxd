# coding: utf-8
# cython: language_level=3, linetrace=True
"""Internal API common to all MEME tools.
"""

from libmeme.alphabet cimport ALPH_T
from libmeme.array cimport ARRAY_T
from libmeme.hash_table cimport HASH_TABLE
from libmeme.data_types cimport WEIGHTS_T, Z_T, LCB_T
from libmeme.matrix cimport MATRIX_T
from libmeme.motif cimport MOTIF_T
from libmeme.motif_in cimport MREAD_T
from libmeme.meme cimport CANDIDATE, DATASET, MODEL, SAMPLE, THETA, MOTIF_SUMMARY
from libmeme.prior_dist cimport PRIOR_DIST_T
from libmeme.pssm cimport PSSM_T
from libmeme.reservoir cimport RESERVOIR_SAMPLER_T
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
    cdef          ARRAY_T* _array
    cdef readonly object   _owner

    cpdef Array copy(self)


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


# --- Matrix -----------------------------------------------------------------

cdef class Matrix:
    cdef          MATRIX_T* _mx
    cdef readonly object    _owner


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

    cpdef PSSM build_pssm(
        self,
        Array bg_freqs,
        Array pv_freqs,
        PriorDist prior_dist = ?,
        double alpha = *,
        int range = *,
        int num_gc_bins = *,
        bint no_log = *,
    )


# --- MotifFile --------------------------------------------------------------

cdef class MotifFile:
    cdef          MREAD_T*  _reader
    cdef          bint      _close
    cdef readonly object    handle
    cdef readonly bytearray buffer

    cpdef void close(self)
    cpdef Motif read(self)


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

cdef class PriorDist:
    cdef          PRIOR_DIST_T* _pd


# --- PSSM -------------------------------------------------------------------

cdef class PSSM:
    cdef          PSSM_T*  _pssm
    cdef readonly Motif    motif


# --- ReservoirSampler -------------------------------------------------------

cdef class ReservoirSampler:
    cdef          RESERVOIR_SAMPLER_T* _reservoir


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
