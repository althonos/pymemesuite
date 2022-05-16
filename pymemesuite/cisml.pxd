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

cdef class MatchedElement:
    cdef readonly object             owner
    cdef          MATCHED_ELEMENT_T* _me


# --- MultiPattern -----------------------------------------------------------

cdef class MultiPattern:
    cdef readonly object           owner
    cdef          MULTI_PATTERN_T* _mp


# --- Pattern ----------------------------------------------------------------

cdef class MatchedElements:
    cdef readonly Pattern             owner
    cdef          MATCHED_ELEMENT_T** _elements

cdef class Pattern:
    cdef readonly object          owner
    cdef          PATTERN_T*      _pattern


# --- ScannedSequence --------------------------------------------------------

cdef class ScannedSequence:
    cdef readonly Pattern             owner
    cdef          SCANNED_SEQUENCE_T* _sseq
