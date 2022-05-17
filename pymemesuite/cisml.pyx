# coding: utf-8
# cython: language_level=3, linetrace=True
"""Types for the CisML data model.

CisML is an XML-based format proposed by Haverty & Weng for reporting
detected sequence motifs. It separates the results about a motif scan
in three main components: *pattern*, *scanned sequence*, and *matched
element*.

See Also:
    - Haverty, Peter M., and Zhiping Weng. *CisML: An XML-Based Format for
    Sequence Motif Detection Software*.
    Bioinformatics 20, no. 11 (22 July 2004): 1815–17.
    :doi:`10.1093/bioinformatics/bth162`.

"""

cimport cython

cimport libmeme.cisml
from libmeme.cisml cimport CISML_T

# --- Python imports ---------------------------------------------------------

from .errors import AllocationError


# --- CisML ------------------------------------------------------------------

cdef class CisML:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._cisml = NULL

    def __init__(self):
        self._cisml = libmeme.cisml.allocate_cisml(b"", b"", b"", b"")

    def __dealloc__(self):
        if self._cisml is not NULL:
            libmeme.cisml.free_cisml(self._cisml)

    def __len__(self):
        assert self._cisml is not NULL
        return libmeme.cisml.get_cisml_num_patterns(self._cisml)

    def __getitem__(self, int index):
        assert self._cisml is not NULL

        cdef Pattern pattern
        cdef int     i       = index
        cdef int     length  = libmeme.cisml.get_cisml_num_patterns(self._cisml)

        if i < 0:
            i += length
        if i < 0 or i >= length:
            raise IndexError(index)

        pattern = Pattern.__new__(Pattern)
        pattern._owner = self
        pattern._pattern = libmeme.cisml.get_cisml_patterns(self._cisml)[i]
        return pattern


# --- MatchedElement ---------------------------------------------------------

@cython.freelist(8)
@cython.no_gc_clear
cdef class MatchedElements:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._owner = None
        self._elements = NULL

    def __init__(self, object owner):
        if isinstance(owner, Pattern):
            self._owner = owner
            self._elements = libmeme.cisml.get_pattern_stored_matches((<Pattern> owner)._pattern)
            self._length = libmeme.cisml.get_pattern_num_stored_matches((<Pattern> owner)._pattern)
        else:
            ty = type(owner).__name__
            raise TypeError(f"Expected `Pattern`, found {ty}")

    def __len__(self):
        assert self._owner is not None
        return self._length

    def __getitem__(self, int index):
        assert self._elements is not NULL
        assert self._owner is not None

        cdef MatchedElement me
        cdef int            i      = index
        cdef int            length = len(self)

        if i < 0:
            i += length
        if i < 0 or i >= length:
            raise IndexError(index)

        me = MatchedElement.__new__(MatchedElement)
        me._owner = self
        me._me = self._elements[i]
        return me


@cython.freelist(8)
@cython.no_gc_clear
cdef class MatchedElement:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._me = NULL
        self._owner = None

    def __dealloc__(self):
        if self._owner is None and self._me is not NULL:
            libmeme.cisml.free_matched_element(self._me)

    def __init__(
        self,
        int start,
        int end,
        ScannedSequence scanned_sequence not None,
    ):
        self._owner = scanned_sequence
        self._me = libmeme.cisml.allocate_matched_element(start, end, scanned_sequence._sseq)

    # --- Properties ---------------------------------------------------------

    @property
    def start(self):
        assert self._me is not NULL
        return libmeme.cisml.get_matched_element_start(self._me)

    @property
    def stop(self):
        assert self._me is not NULL
        return libmeme.cisml.get_matched_element_stop(self._me)

    @property
    def score(self):
        assert self._me is not NULL
        if not libmeme.cisml.has_matched_element_score(self._me):
            return None
        return libmeme.cisml.get_matched_element_score(self._me)

    @property
    def pvalue(self):
        assert self._me is not NULL
        if not libmeme.cisml.has_matched_element_pvalue(self._me):
            return None
        return libmeme.cisml.get_matched_element_pvalue(self._me)

    @property
    def qvalue(self):
        assert self._me is not NULL
        if not libmeme.cisml.has_matched_element_qvalue(self._me):
            return None
        return libmeme.cisml.get_matched_element_qvalue(self._me)

    @property
    def strand(self):
        assert self._me is not NULL
        return chr(libmeme.cisml.get_matched_element_strand(self._me))

    @property
    def sequence(self):
        assert self._me is not NULL
        cdef const char* seq = libmeme.cisml.get_matched_element_sequence(self._me)
        return None if seq is NULL else seq.decode('ascii')

    @property
    def source(self):
        assert self._me is not NULL
        cdef ScannedSequence sseq = ScannedSequence.__new__(ScannedSequence)
        sseq._owner = self
        sseq._sseq = libmeme.cisml.get_matched_element_scanned_seq(self._me)
        return sseq


# --- MultiPattern -----------------------------------------------------------

@cython.no_gc_clear
cdef class MultiPattern:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._mp = NULL
        self._owner = None

    def __dealloc__(self):
        if self._owner is None and self._mp is not NULL:
            libmeme.cisml.free_multi_pattern(self._mp)

    # --- Properties ---------------------------------------------------------


# --- Pattern ----------------------------------------------------------------

@cython.no_gc_clear
cdef class Pattern:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._pattern = NULL
        self._owner = None

    def __dealloc__(self):
        if self._owner is None and self._pattern is not NULL:
            libmeme.cisml.free_pattern(self._pattern)

    def __init__(self, bytes accession not None, bytes name not None):
        self._owner = None
        self._pattern = libmeme.cisml.allocate_pattern(accession, name)
        if self._pattern is NULL:
            raise AllocationError("PATTERN_T", sizeof(PATTERN_T*))

    # --- Properties ---------------------------------------------------------

    @property
    def name(self):
        assert self._pattern is not NULL
        return libmeme.cisml.get_pattern_name(self._pattern)

    @property
    def accession(self):
        assert self._pattern is not NULL
        return libmeme.cisml.get_pattern_accession(self._pattern)

    @property
    def pvalue(self):
        assert self._pattern is not NULL
        if not libmeme.cisml.has_pattern_pvalue(self._pattern):
            return None
        return libmeme.cisml.get_pattern_pvalue(self._pattern)

    @pvalue.setter
    def pvalue(self, object value):
        assert self._pattern is not NULL
        if value is None:
            libmeme.cisml.clear_pattern_pvalue(self._pattern)
        libmeme.cisml.set_pattern_pvalue(self._pattern, value)

    @pvalue.deleter
    def pvalue(self):
        assert self._pattern is not NULL
        libmeme.cisml.clear_pattern_pvalue(self._pattern)

    @property
    def score(self):
        assert self._pattern is not NULL
        if not libmeme.cisml.has_pattern_score(self._pattern):
            return None
        return libmeme.cisml.get_pattern_score(self._pattern)

    @score.setter
    def score(self, object value):
        assert self._pattern is not NULL
        if value is None:
            libmeme.cisml.clear_pattern_score(self._pattern)
        libmeme.cisml.set_pattern_score(self._pattern, value)

    @score.deleter
    def score(self):
        assert self._pattern is not NULL
        libmeme.cisml.clear_pattern_score(self._pattern)

    @property
    def matched_elements(self):
        assert self._pattern is not NULL
        if not libmeme.cisml.get_pattern_is_complete(self._pattern):
            raise RuntimeError("Cannot access matched elements until the pattern is complete")
        return MatchedElements(self)

    @property
    def scanned_sequences(self):
        assert self._pattern is not NULL
        if not libmeme.cisml.get_pattern_is_complete(self._pattern):
            raise RuntimeError("Cannot access scanned sequences until the pattern is complete")
        return ScannedSequences(self)


# --- ScannedSequence --------------------------------------------------------

@cython.freelist(8)
@cython.no_gc_clear
cdef class ScannedSequences:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._owner = None
        self._sequences = NULL

    def __init__(self, object owner):
        if isinstance(owner, Pattern):
            self._owner = owner
            self._sequences = libmeme.cisml.get_pattern_scanned_sequences((<Pattern> owner)._pattern)
            self._length = libmeme.cisml.get_pattern_num_scanned_sequences((<Pattern> owner)._pattern)
        else:
            ty = type(owner).__name__
            raise TypeError(f"Expected `Pattern`, found {ty}")

    def __len__(self):
        assert self._owner is not None
        return self._length

    def __getitem__(self, int index):
        assert self._sequences is not NULL
        assert self._owner is not None

        cdef ScannedSequence ss
        cdef int             i      = index
        cdef int             length = len(self)

        if i < 0:
            i += length
        if i < 0 or i >= length:
            raise IndexError(index)

        ss = ScannedSequence.__new__(ScannedSequence)
        ss._owner = self
        ss._sseq = self._sequences[i]
        return ss


@cython.freelist(8)
@cython.no_gc_clear
cdef class ScannedSequence:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._sseq = NULL
        self._owner = None

    def __dealloc__(self):
        if self._owner is None and self._sseq is not NULL:
            libmeme.cisml.free_scanned_sequence(self._sseq)

    def __init__(
        self,
        bytes accession not None,
        bytes name not None,
        Pattern pattern not None,
    ):
        self._owner = pattern
        self._sseq = libmeme.cisml.allocate_scanned_sequence(accession, name, pattern._pattern)

    # --- Properties ---------------------------------------------------------

    @property
    def accession(self):
        assert self._sseq is not NULL
        return <bytes> libmeme.cisml.get_scanned_sequence_accession(self._sseq)

    @property
    def name(self):
        assert self._sseq is not NULL
        return <bytes> libmeme.cisml.get_scanned_sequence_name(self._sseq)

    @property
    def score(self):
        assert self._sseq is not NULL
        if not libmeme.cisml.has_scanned_sequence_score(self._sseq):
            return None
        return libmeme.cisml.get_scanned_sequence_score(self._sseq)

    @property
    def pvalue(self):
        assert self._sseq is not NULL
        if not libmeme.cisml.has_scanned_sequence_pvalue(self._sseq):
            return None
        return libmeme.cisml.get_scanned_sequence_pvalue(self._sseq)
