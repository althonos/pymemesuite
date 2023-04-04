# coding: utf-8
# cython: language_level=3, linetrace=True
"""Internal API common to all MEME tools.
"""

from libmeme.alphabet cimport ALPH_T
from libmeme.array cimport ATYPE, ARRAY_T
from libmeme.hash_table cimport HASH_TABLE
from libmeme.data_types cimport WEIGHTS_T, Z_T, LCB_T
from libmeme.matrix cimport MATRIX_T, MTYPE
from libmeme.motif cimport MOTIF_T
from libmeme.motif_in cimport MREAD_T

from libmeme.prior_dist cimport PRIOR_DIST_T
from libmeme.pssm cimport PSSM_T
from libmeme.reservoir cimport RESERVOIR_SAMPLER_T
from libmeme.seq cimport SEQ_T
from libmeme.user cimport MINSITES

# --- Python imports ---------------------------------------------------------

import os
import warnings

from .errors import AllocationError


# --- Alphabet ---------------------------------------------------------------

cdef class Alphabet:
    cdef ALPH_T* _alph


# --- Array ------------------------------------------------------------------

cdef class Array:
    cdef          ARRAY_T*   _array
    cdef          Py_ssize_t _len
    cdef readonly object     _owner

    cpdef Array copy(self)
    cpdef ATYPE sum(self)


# --- Background -------------------------------------------------------------

cdef class Background:
    cdef readonly Alphabet alphabet
    cdef readonly Array    frequencies

    cpdef Background copy(self)


# --- Matrix -----------------------------------------------------------------

cdef class Matrix:
    cdef          MATRIX_T* _mx
    cdef readonly object    _owner

    cdef Array _get_row(self, int row)
    cdef MTYPE _get_element(self, int row, int col)
    cdef int _set_element(self, int row, int col, MTYPE value) except -1


# --- Motif ------------------------------------------------------------------

cdef class Motif:
    cdef          MOTIF_T*   _motif
    cdef readonly Alphabet alphabet

    cpdef PSSM build_pssm(
        self,
        Background background,
        Background background_p = ?,
        PriorDist prior_dist = ?,
        double alpha = *,
        int range = *,
        int num_gc_bins = *,
        bint no_log = *,
    )
    cpdef Motif reverse_complement(self)


# --- MotifFile --------------------------------------------------------------

cdef class MotifFile:
    cdef          MREAD_T*  _reader
    cdef          bint      _close
    cdef readonly object    handle
    cdef readonly bytearray buffer

    cpdef void close(self)
    cpdef Motif read(self)


# --- PriorDistribution ------------------------------------------------------

cdef class PriorDist:
    cdef          PRIOR_DIST_T* _pd


# --- PSSM -------------------------------------------------------------------

cdef class PSSM:
    cdef          PSSM_T*  _pssm
    cdef readonly Alphabet alphabet
    cdef readonly Motif    motif

    cpdef PSSM copy(self)
    cpdef PSSM reverse_complement(self)


# --- ReservoirSampler -------------------------------------------------------

cdef class ReservoirSampler:
    cdef          RESERVOIR_SAMPLER_T* _reservoir


# --- Sequence ---------------------------------------------------------------

cdef class Sequence:
  cdef          SEQ_T*  _seq
