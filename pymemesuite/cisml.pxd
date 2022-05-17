# coding: utf-8
# cython: language_level=3, linetrace=True

from libmeme.cisml cimport (
    CISML_T,
    MATCHED_ELEMENT_T,
    MULTI_PATTERN_T,
    PATTERN_T,
    SCANNED_SEQUENCE_T,
)

# --- CisML ------------------------------------------------------------------

cdef class CisML:
    cdef CISML_T* _cisml


# --- MatchedElement ---------------------------------------------------------

cdef class MatchedElements:
    cdef readonly object              _owner
    cdef readonly int                 _length
    cdef          MATCHED_ELEMENT_T** _elements

cdef class MatchedElement:
    cdef readonly object             _owner
    cdef          MATCHED_ELEMENT_T* _me


# --- MultiPattern -----------------------------------------------------------

cdef class MultiPattern:
    cdef readonly object           _owner
    cdef          MULTI_PATTERN_T* _mp


# --- Pattern ----------------------------------------------------------------

cdef class Pattern:
    cdef readonly object          _owner
    cdef          PATTERN_T*      _pattern


# --- ScannedSequence --------------------------------------------------------

cdef class ScannedSequences:
    cdef readonly object               _owner
    cdef readonly int                  _length
    cdef          SCANNED_SEQUENCE_T** _sequences

cdef class ScannedSequence:
    cdef readonly object              _owner
    cdef          SCANNED_SEQUENCE_T* _sseq
