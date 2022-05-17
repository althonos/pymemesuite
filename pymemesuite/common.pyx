# coding: utf-8
# cython: language_level=3, linetrace=True
"""Internal API common to all MEME tools.
"""

cimport cython
from cpython.buffer cimport PyBUF_FORMAT, PyBUF_READ, PyBUF_WRITE, PyBuffer_FillInfo
from cpython.bytes cimport PyBytes_FromString, PyBytes_FromStringAndSize, PyBytes_AsString
from cpython.mem cimport PyMem_Free, PyMem_Malloc, PyMem_Realloc
from cpython.memoryview cimport PyMemoryView_FromMemory, PyMemoryView_GET_BUFFER
from libc.errno cimport errno
from libc.limits cimport LONG_MAX
from libc.math cimport INFINITY
from libc.stdint cimport uint8_t
from libc.stdio cimport FILE, fopen, fclose
from libc.string cimport strdup, strcmp, strcpy, strncpy, memcpy, memset
from libc.stdlib cimport calloc, free, malloc, realloc

cimport libmeme.alphabet
cimport libmeme.array
cimport libmeme.macros
cimport libmeme.matrix
cimport libmeme.meme
cimport libmeme.motif
cimport libmeme.motif_in
cimport libmeme.pssm
cimport libmeme.read_sequence
cimport libmeme.read_seq_file
cimport libmeme.reservoir
cimport libmeme.seq
from libmeme.alphabet cimport ALPH_T
from libmeme.array cimport ATYPE
from libmeme.hash_table cimport HASH_TABLE
from libmeme.data_types cimport WEIGHTS_T, Z_T, LCB_T
from libmeme.meme cimport CANDIDATE, DATASET, MODEL, OBJTYPE, SAMPLE, THETA, P_POINT, BRANCH_PARAMS, POINT_BRANCHES
from libmeme.mtype cimport MOTYPE
from libmeme.user cimport MINSITES, BFACTOR, HSIZE, HS_DECREASE

# --- Python imports ---------------------------------------------------------

import array
import os
import warnings

from .errors import AllocationError


# --- Constants --------------------------------------------------------------

cdef dict MEME_OBJECTIVE_FUNCTIONS = {
    "classic": OBJTYPE.Classic,
    "nc": OBJTYPE.NC,
    "se": OBJTYPE.SE,
    "smhg": OBJTYPE.SE,
    "de": OBJTYPE.DE,
    "hs": OBJTYPE.DE,
    "cv": OBJTYPE.DE,
    "nz": OBJTYPE.NZ,
    "ce": OBJTYPE.CE,
    "cd": OBJTYPE.CD,
}

cdef dict MEME_MODEL_TYPES = {
    "anr": MOTYPE.Tcm,
    "tcm": MOTYPE.Tcm,
    "oops": MOTYPE.Oops,
    "zoops": MOTYPE.Zoops,
}

# --- Utilities --------------------------------------------------------------

cdef void* allocate(size_t size, str typename) except NULL:
    cdef void* tmp = malloc(size)
    if tmp == NULL:
        raise AllocationError(typename, size)
    return tmp

cdef void* matrix_create(size_t nrows, size_t ncols, size_t itemsize) except NULL:
    cdef size_t  r
    cdef void** matrix

    # allocate array of pointers
    matrix = <void**> malloc(sizeof(void*) * nrows)
    if matrix == NULL:
        raise AllocationError("double*", sizeof(double), nrows)

    # allocate contiguous memory
    matrix[0] = <void*> malloc(itemsize * nrows * ncols)
    if matrix[0] == NULL:
        free(matrix)
        raise AllocationError("double", sizeof(double), nrows * ncols)

    # update pointer offsets
    for r in range(nrows):
        matrix[r] = matrix[0] + r*ncols*itemsize
    return matrix

cdef void matrix_free(void** matrix):
    if matrix != NULL and matrix[0] != NULL:
        free(matrix[0])
    if matrix != NULL:
        free(matrix)


# --- Alphabet ---------------------------------------------------------------

cdef class Alphabet:

    # --- Class methods ------------------------------------------------------

    @classmethod
    def protein(cls):
        """protein(cls)\n--

        Create a default protein alphabet.

        """
        cdef Alphabet alphabet = Alphabet.__new__(Alphabet)
        alphabet._alph = libmeme.alphabet.alph_protein()
        if alphabet._alph == NULL:
            raise RuntimeError("Failed to create protein alphabet")
        return alphabet

    @classmethod
    def amino(cls):
        """amino(cls)\n--

        Create a default protein alphabet (alias for `Alphabet.protein`).

        """
        return cls.protein()

    @classmethod
    def dna(cls):
        """dna(cls)\n--

        Create a default DNA alphabet.

        """
        cdef Alphabet alphabet = Alphabet.__new__(Alphabet)
        alphabet._alph = libmeme.alphabet.alph_dna()
        if alphabet._alph == NULL:
            raise RuntimeError("Failed to create DNA alphabet")
        return alphabet

    @classmethod
    def rna(cls):
        """rna(cls)\n--

        Create a default RNA alphabet.

        """
        cdef Alphabet alphabet = Alphabet.__new__(Alphabet)
        alphabet._alph = libmeme.alphabet.alph_rna()
        if alphabet._alph == NULL:
            raise RuntimeError("Failed to create DNA alphabet")
        return alphabet

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._alph = NULL

    def __dealloc__(self):
        libmeme.alphabet.alph_release(self._alph)

    def __repr__(self):
        assert self._alph is not NULL
        if libmeme.alphabet.alph_is_builtin_rna(self._alph):
            return "Alphabet.rna()"
        elif libmeme.alphabet.alph_is_builtin_dna(self._alph):
            return "Alphabet.dna()"
        elif libmeme.alphabet.alph_is_builtin_protein(self._alph):
            return "Alphabet.protein()"
        else:
            return "Alphabet()"

    # --- Properties ---------------------------------------------------------

    @property
    def size(self):
        """`int`: The number of core letters in the alphabet.
        """
        assert self._alph is not NULL
        return libmeme.alphabet.alph_size_core(self._alph)

    @property
    def size_wild(self):
        """`int`: The number of core and wildcard letters in the alphabet.
        """
        assert self._alph is not NULL
        return libmeme.alphabet.alph_size_wild(self._alph)

    @property
    def size_full(self):
        """`int`: The number of core and wildcard letters in the alphabet.
        """
        assert self._alph is not NULL
        return libmeme.alphabet.alph_size_full(self._alph)

    @property
    def wildcard(self):
        """`str`: The wildcard letter of the alphabet.
        """
        assert self._alph is not NULL
        return chr(libmeme.alphabet.alph_wildcard(self._alph))

    @property
    def symbols(self):
        """`str`: The symbols in the alphabet.
        """
        assert self._alph is not NULL
        return self._alph.symbols.decode('ascii')

# --- Array ------------------------------------------------------------------

@cython.freelist(8)
cdef class Array:
    """A 1D vector of fixed size with double-precision elements.
    """

    # --- Class methods ------------------------------------------------------

    @classmethod
    def zeros(cls, int n):
        """zeros(cls, n)\n--

        Create a new array of size ``n`` filled with zeros.

        """
        if n < 0:
            raise ValueError("Cannot create a vector with negative size")

        cdef Array array = Array.__new__(Array)
        array._owner = None
        array._array = libmeme.array.allocate_array(n)
        if array._array is NULL:
            raise AllocationError("ARRAY_T", sizeof(ARRAY_T))
        with nogil:
            libmeme.array.init_array(0, array._array)

        return array

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._array = NULL
        self._owner = None

    def __init__(self, object iterable = ()):
        """__init__(self, iterable=())\n--

        Create a new array from the given iterable of values.

        """
        cdef int        i
        cdef ATYPE      item
        cdef ATYPE[::1] view
        cdef int        length = len(iterable)

        self._owner = None
        self._array = libmeme.array.allocate_array(length)
        if self._array is NULL:
            raise AllocationError("ARRAY_T", sizeof(ARRAY_T))

        try:
            # attempt to use a memory view to copy the data in one go
            view = iterable
            with nogil:
                libmeme.array.fill_array(&view[0], self._array)
        except ValueError:
            # if the iterable is not in contiguous memory, revert back
            # to manually setting array elements one-by-one
            for i, item in enumerate(iterable):
                libmeme.array.set_array_item(i, item, self._array)

    def __dealloc__(self):
        if self._owner is None:
            libmeme.array.free_array(self._array)

    def __bool__(self):
        self._array
        return libmeme.array.get_array_length(self._array) != 0

    def __len__(self):
        assert self._array is not NULL
        return libmeme.array.get_array_length(self._array)

    def __copy__(self):
        return self.copy()

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        assert self._array is not NULL

        if flags & PyBUF_FORMAT:
            buffer.format = b"d"
        else:
            buffer.format = NULL

        buffer.buf = libmeme.array.raw_array(self._array)
        buffer.internal = NULL
        buffer.itemsize = sizeof(double)
        buffer.len = libmeme.array.get_array_length(self._array)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = NULL
        buffer.suboffsets = NULL
        buffer.strides = NULL

    def __getitem__(self, int index):
        assert self._array is not NULL

        cdef int length = libmeme.array.get_array_length(self._array)
        cdef int i      = index

        if i < 0:
            i += index
        if i < 0 or i >= length:
            raise IndexError(index)

        return libmeme.array.get_array_item(i, self._array)

    def __repr__(self):
        cdef type ty   = type(self)
        cdef str  name = ty.__name__
        cdef str  mod  = ty.__module__
        return f"{mod}.{name}({list(self)!r})"

    def __sizeof__(self):
        assert self._array is not NULL
        return (
            sizeof(self)
          + sizeof(ARRAY_T)
          + libmeme.array.get_array_length(self._array) * sizeof(ATYPE)
        )

    def __reduce__(self):
        assert self._array is not NULL
        cdef object buffer = array.array("d")
        buffer.frombytes(memoryview(self).cast("b"))
        return Array, (buffer,)

    # --- Properties ---------------------------------------------------------

    @property
    def itemsize(self):
        return sizeof(ATYPE)

    @property
    def format(self):
        return "d"

    # --- Methods ------------------------------------------------------------

    cpdef Array copy(self):
        """copy(self)\n--

        Create a copy of the array, allocating a new buffer.

        """
        assert self._array is not NULL

        cdef int   length = libmeme.array.get_array_length(self._array)
        cdef Array copy   = Array.__new__(Array)

        copy._array = libmeme.array.allocate_array(length)
        if copy._array is NULL:
            raise AllocationError("ARRAY_T", sizeof(ARRAY_T))
        with nogil:
            libmeme.array.copy_array(self._array, copy._array)

        return copy


# --- Matrix -----------------------------------------------------------------

cdef class Matrix:

    # --- Class methods ------------------------------------------------------

    @classmethod
    def zeros(cls, int m, int n):
        if m < 0 or n < 0:
            raise ValueError("Cannot create a matrix with negative dimension")

        cdef Matrix matrix = Matrix.__new__(Matrix)
        matrix._owner = None
        matrix._mx = libmeme.matrix.allocate_matrix(m, n)
        if matrix._mx is NULL:
            raise AllocationError("MATRIX_T", sizeof(MATRIX_T))
        with nogil:
            libmeme.matrix.init_matrix(0, matrix._mx)

        return matrix

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._owner = None
        self._mx = NULL

    def __dealloc__(self):
        if self._owner is None:
            libmeme.matrix.free_matrix(self._mx)


# --- Motif ------------------------------------------------------------------

cdef class Motif:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._motif = NULL
        self.alphabet = None

    def __dealloc__(self):
        libmeme.motif.free_motif(self._motif)

    def __init__(
        self,
        Alphabet alphabet not None,
        *,
        Matrix frequencies = None,
        Matrix scores = None,
        bytes name = None,
        bytes accession = None,
    ):
        if frequencies is None and scores is None:
            raise ValueError("Either `frequencies` or `scores` are required to create a `Motif`")

        self.alphabet = alphabet
        self._motif = libmeme.motif.allocate_motif(
            b"",
            b"",
            alphabet._alph,
            NULL if frequencies is None else frequencies._mx,
            NULL if scores is None else scores._mx,
        )

    # --- Properties ---------------------------------------------------------

    @property
    def accession(self):
        assert self._motif is not NULL
        return libmeme.motif.get_motif_id(self._motif)

    @property
    def name(self):
        assert self._motif is not NULL
        return libmeme.motif.get_motif_id2(self._motif)

    @property
    def strand(self):
        assert self._motif is not NULL
        return chr(libmeme.motif.get_motif_strand(self._motif))

    @property
    def width(self):
        assert self._motif is not NULL
        return libmeme.motif.get_motif_length(self._motif)

    @property
    def frequencies(self):
        assert self._motif is not NULL
        cdef Matrix matrix = Matrix.__new__(Matrix)
        matrix._owner = self
        matrix._mx = libmeme.motif.get_motif_freqs(self._motif)
        return matrix

    @property
    def scores(self):
        assert self._motif is not NULL
        cdef Matrix matrix = Matrix.__new__(Matrix)
        matrix._owner = self
        matrix._mx = libmeme.motif.get_motif_scores(self._motif)
        return matrix

    @property
    def evalue(self):
        assert self._motif is not NULL
        return libmeme.motif.get_motif_evalue(self._motif)

    @property
    def log_evalue(self):
        assert self._motif is not NULL
        return libmeme.motif.get_motif_log_evalue(self._motif)

    @property
    def consensus(self):
        assert self._motif is not NULL
        return libmeme.motif.get_motif_consensus(self._motif).decode('ascii')

    @property
    def url(self):
        assert self._motif is not NULL
        if not libmeme.motif.has_motif_url(self._motif):
            return None
        return libmeme.motif.get_motif_url(self._motif).decode('ascii')

    @property
    def complexity(self):
        assert self._motif is not NULL
        return libmeme.motif.get_motif_complexity(self._motif)

    # --- Methods ------------------------------------------------------------

    cpdef PSSM build_pssm(
        self,
        Array bg_freqs,
        Array pv_freqs,
        PriorDist prior_dist = None,
        double alpha = 1.0,
        int range = 1,
        int num_gc_bins = 0,
        bint no_log = False,
    ):
        """build_pssm(self, bg_freqs, pv_freqs, prior_dist=None, alpha=1.0, range=1, num_gc_bins=0, no_log=False)

        Build a `PSSM` from this motif.

        Arguments:
            bg_freqs (`Array`): The background frequencies.
            pv_freqs (`Array`): The background frequencies for the p-values.
            prior_dist (`PriorDist`, optional): The distribution of priors,
                or `None`.
            alpha (`float`): The scale factor for non-specific priors. Only
                used when ``prior_dist`` is not `None`.
            range (`int`): The range of scaled scores.
            num_gc_bins (`int`): A number of GC bins to use to create the
                p-value tables instead of using the ``pv_freqs`` array.
            no_log (`bool`): Set to `True` to build a likelihood ratio PSSM.

        """
        assert self._motif is not NULL

        if range <= 0:
            raise ValueError("``range`` must be strictly positive")
        if len(bg_freqs) != self.alphabet.size:
            raise ValueError("``bg_freqs`` length is inconsistent with motif alphabet")
        if len(bg_freqs) != self.alphabet.size:
            raise ValueError("``pv_freqs`` length is inconsistent with motif alphabet")

        cdef PSSM pssm = PSSM.__new__(PSSM)
        with nogil:
            pssm._pssm = libmeme.pssm.build_motif_pssm(
                self._motif,
                bg_freqs._array,
                pv_freqs._array,
                NULL if prior_dist is None else prior_dist._pd,
                alpha,
                range,
                num_gc_bins,
                no_log
            )
            pssm._pssm.motif = self._motif
        pssm.motif = self
        pssm.alphabet = Alphabet.__new__(Alphabet)
        pssm.alphabet._alph = libmeme.alphabet.alph_hold(pssm._pssm.alph)

        return pssm

    cpdef Motif reverse_complement(self):
        """reverse_complement(self)\n--

        Create a new motif with the reverse-complement of this motif.

        """
        assert self._motif is not NULL
        assert self.alphabet is not None

        if not libmeme.alphabet.alph_has_complement(libmeme.motif.get_motif_alph(self._motif)):
            raise ValueError("Cannot reverse-complement a motif in a non-complementable alphabet")

        cdef Motif rc = Motif.__new__(Motif)
        with nogil:
            rc._motif = libmeme.motif.dup_rc_motif(self._motif)
        if rc._motif is NULL:
            raise AllocationError("MOTIF_T", sizeof(MOTIF_T*))
        rc.alphabet = Alphabet.__new__(Alphabet)
        rc.alphabet._alph = libmeme.alphabet.alph_hold(libmeme.motif.get_motif_alph(rc._motif))
        return rc


# --- MotifFile --------------------------------------------------------------

cdef class MotifFile:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._reader = NULL
        self._close = False
        self.handle = None
        self.buffer = None

    def __init__(
        self,
        object handle,
        *,
        bint symmetrical=False,
        double pseudocount=0.1,
    ):
        """__init__(self, handle, *, symmetrical=False, pseudocount=0.1)\n--

        Create a new motif reader from a path or file-like object.

        Arguments:
            handle (`str`, `os.PathLike` or file-like object): The handle to
                the file to read the motif data from.

        Keyword Arguments:
            symmetrical (`bool`): Set to `True` to make the background
                symmetrical, provided the motif alphabet is complementable.
            pseudocount (`float`): The pseudo-count to be applied to the
                motif being read.

        """
        self.buffer = bytearray(2048)

        try:
            self.handle = open(handle, "rb")
            path = handle
            self._close = True
        except TypeError:
            self.handle = handle
            path = getattr(handle, "name", None)
            self._close = False

        if path is not None:
            name = os.fsencode(path)
            self._reader = libmeme.motif_in.mread_create(name, 0, symmetrical)
        else:
            self._reader = libmeme.motif_in.mread_create(NULL, 0, symmetrical)
        if self._reader is NULL:
            raise AllocationError("MREAD_T*", sizeof(MREAD_T*))

        libmeme.motif_in.mread_set_bg_source(self._reader, NULL, NULL)
        libmeme.motif_in.mread_set_pseudocount(self._reader, pseudocount)

    def __dealloc__(self):
        if self._close and not self.handle.closed:
            warnings.warn("unclosed motif file", ResourceWarning)
            self.close()
        libmeme.motif_in.mread_destroy(self._reader)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        cdef Motif motif = self.read()
        if motif is None:
            raise StopIteration
        return motif

    # --- Properties ---------------------------------------------------------

    @property
    def alphabet(self):
        assert self._reader is not NULL

        cdef ALPH_T*  alph
        cdef Alphabet alphabet

        alph = libmeme.motif_in.mread_get_alphabet(self._reader)
        if alph is NULL:
            return None

        alphabet = Alphabet.__new__(Alphabet)
        alphabet._alph = libmeme.alphabet.alph_hold(alph)
        return alphabet

    @property
    def background(self):
        assert self._reader is not NULL
        cdef Array array = Array.__new__(Array)
        array._owner = None
        array._array = libmeme.motif_in.mread_get_background(self._reader)
        return None if array._array is NULL else array

    # --- Methods ------------------------------------------------------------

    cpdef void close(self):
        """close(self)\n--

        Close the file and free the resources used by the reader.

        """
        if self._close:
            self.handle.close()

    cpdef Motif read(self):
        assert self._reader is not NULL

        cdef int       bytes_read = 1
        cdef char[::1] view       = self.buffer
        cdef Motif     motif      = Motif.__new__(Motif)

        while bytes_read > 0 and not libmeme.motif_in.mread_has_motif(self._reader):
            bytes_read = self.handle.readinto(self.buffer)
            libmeme.motif_in.mread_update(self._reader, &view[0], bytes_read, bytes_read == 0)

        with nogil:
            motif._motif = libmeme.motif_in.mread_next_motif(self._reader)
        if motif._motif is NULL:
            return None

        motif.alphabet = Alphabet.__new__(Alphabet)
        motif.alphabet._alph = libmeme.alphabet.alph_hold(libmeme.motif.get_motif_alph(motif._motif))
        return motif


# --- PriorDistribution ------------------------------------------------------

cdef class PriorDist:

    # --- Magic Methods ------------------------------------------------------

    def __cinit__(self):
        self._pd = NULL


# --- PSSM -------------------------------------------------------------------

cdef class PSSM:

    # --- Magic Methods ------------------------------------------------------

    def __cinit__(self):
        self.motif = None
        self._pssm = NULL

    def __dealloc__(self):
        libmeme.pssm.free_pssm(self._pssm)

    def __copy__(self):
        assert self._pssm is not NULL
        return self.copy()

    # --- Properties ---------------------------------------------------------

    @property
    def width(self):
        assert self._pssm is not NULL
        return libmeme.pssm.get_pssm_w(self._pssm)

    @property
    def matrix(self):
        assert self._pssm is not NULL
        cdef Matrix matrix = Matrix.__new__(Matrix)
        matrix._mx = self._pssm.matrix
        matrix._owner = self
        return matrix

    # --- Methods ------------------------------------------------------------

    cpdef PSSM copy(self):
        assert self._pssm is not NULL

        cdef int  i
        cdef int  length
        cdef PSSM copy   = PSSM.__new__(PSSM)

        with nogil:
            copy._pssm = libmeme.pssm.allocate_pssm(
                libmeme.pssm.get_pssm_alph(self._pssm),
                libmeme.pssm.get_pssm_w(self._pssm),
                libmeme.pssm.get_pssm_alphsize(self._pssm),
                self._pssm.num_gc_bins,
            )
            libmeme.matrix.copy_matrix(self._pssm.matrix, copy._pssm.matrix)

            copy._pssm.matrix_is_log = self._pssm.matrix_is_log
            copy._pssm.matrix_is_scaled = self._pssm.matrix_is_scaled
            copy._pssm.scale = self._pssm.scale
            copy._pssm.offset = self._pssm.offset
            copy._pssm.range = self._pssm.range

            if self._pssm.pv is not NULL:
                copy._pssm.pv = libmeme.array.allocate_array(libmeme.array.get_array_length(self._pssm.pv))
                libmeme.array.copy_array(self._pssm.pv, copy._pssm.pv)
            else:
                copy._pssm.pv = NULL

            if self._pssm.gc_pv is not NULL:
                copy._pssm.gc_pv = <ARRAY_T**> calloc(self._pssm.num_gc_bins, sizeof(ARRAY_T*))
                if copy._pssm.gc_pv is NULL:
                    raise AllocationError("ARRAY_T*", sizeof(ARRAY_T*), self._pssm.num_gc_bins)
                for i in range(self._pssm.num_gc_bins):
                    length = libmeme.array.get_array_length(self._pssm.gc_pv[i])
                    copy._pssm.gc_pv[i] = libmeme.array.allocate_array(length)
                    libmeme.array.copy_array(self._pssm.gc_pv[i], copy._pssm.gc_pv[i])
            else:
                copy._pssm.gc_pv = NULL

            copy._pssm.num_gc_bins = self._pssm.num_gc_bins
            copy._pssm.min_score = self._pssm.min_score
            copy._pssm.max_score = self._pssm.max_score
            copy._pssm.nolog_max_score = self._pssm.nolog_max_score

        copy.alphabet = Alphabet.__new__(Alphabet)
        copy.alphabet._alph = libmeme.alphabet.alph_hold(copy._pssm.alph)

        if self._pssm.motif is not NULL:
            copy.motif = Motif.__new__(Motif)
            copy.motif._motif = libmeme.motif.duplicate_motif(self._pssm.motif)
            copy._pssm.motif = copy.motif._motif
            copy.motif.alphabet = copy.alphabet
        else:
            copy.motif = None
            copy._pssm.motif = NULL

        return copy

    cpdef PSSM reverse_complement(self):
        assert self._pssm is not NULL

        cdef int      i
        cdef ARRAY_T* left_scores
        cdef ARRAY_T* right_scores
        cdef ALPH_T*  alph         = libmeme.pssm.get_pssm_alph(self._pssm)
        cdef PSSM     rc           = self.copy()
        cdef int      length       = libmeme.matrix.get_num_rows(rc._pssm.matrix)

        with nogil:
            if rc._pssm.motif is not NULL:
                libmeme.motif.reverse_complement_motif(rc._pssm.motif)
            for i in range((length + 1) // 2):
                left_scores = libmeme.matrix.get_matrix_row(i, rc._pssm.matrix)
                right_scores = libmeme.matrix.get_matrix_row(length - i - 1, rc._pssm.matrix)
                libmeme.alphabet.complement_swap_freqs(alph, left_scores, right_scores)

        return rc


# --- ReservoirSampler -------------------------------------------------------

cdef class ReservoirSampler:

    # --- Magic Methods ------------------------------------------------------

    def __cinit__(self):
        self._reservoir = NULL

    def __init__(self, size_t size):
        self._reservoir = libmeme.reservoir.new_reservoir_sampler(size, NULL)
        if self._reservoir is NULL:
            raise AllocationError("RESERVOIR_SAMPLER_T", sizeof(RESERVOIR_SAMPLER_T*))

    def __dealloc__(self):
        libmeme.reservoir.free_reservoir(self._reservoir)

    # --- Properties ---------------------------------------------------------

    @property
    def samples_seen(self):
        """`int`: The number of samples seen by the reservoir sampler.
        """
        assert self._reservoir is not NULL
        return libmeme.reservoir.get_reservoir_num_samples_seen(self._reservoir)

    @property
    def samples_retained(self):
        """`int`: The number of samples retained in the reservoir sampler.
        """
        assert self._reservoir is not NULL
        return libmeme.reservoir.get_reservoir_num_samples_retained(self._reservoir)

    @property
    def values(self):
        assert self._reservoir is not NULL
        cdef size_t  length  = libmeme.reservoir.get_reservoir_num_samples_retained(self._reservoir)
        cdef double* samples = libmeme.reservoir.get_reservoir_samples(self._reservoir)
        cdef object mem = PyMemoryView_FromMemory(
            <char*> samples,
            length * sizeof(double),
            PyBUF_READ | PyBUF_WRITE
        )

        # DANGER (@althonos): We need to make sure that the object is not
        #                     deallocated before the memoryview, so we
        #                     register the `ReservoirSampler` object (self) as
        #                     the memoryview exporter. This needs a bit of
        #                     tweaking with the internal Py_buffer, because
        #                     `PyMemoryView_FromMemory` doesn't allow setting
        #                     up an exporter, and directly setting the `obj`
        #                     attribute is unsafe.
        PyBuffer_FillInfo(
            PyMemoryView_GET_BUFFER(mem),
            self,
            <char*> samples,
            length * sizeof(double),
            True,
            PyBUF_READ,
        )

        return mem.cast('d')


# --- Sequence ---------------------------------------------------------------

cdef class Sequence:

    def __cinit__(self):
        self._seq = NULL

    def __init__(
        self,
        str sequence not None,
        bytes name = None,
        bytes description = None,
        unsigned int offset = 0,
    ):
        self._seq = libmeme.seq.allocate_seq(
            NULL if name is None else <char*> name,
            NULL if description is None else <char*> description,
            offset,
            sequence.encode('ascii')
        )
        if self._seq is NULL:
            raise AllocationError("SEQ_T", sizeof(SEQ_T*))

    def __len__(self):
        assert self._seq is not NULL
        return libmeme.seq.get_seq_length(self._seq)

    def __getitem__(self, int index):
        assert self._seq is not NULL

        cdef int length = libmeme.seq.get_seq_length(self._seq)
        cdef int i      = index

        if i < 0:
            i += index
        if i < 0 or i >= length:
            raise IndexError(index)

        return chr(libmeme.seq.get_seq_char(i, self._seq))

    @property
    def name(self):
        assert self._seq is not NULL
        return <bytes> libmeme.seq.get_seq_name(self._seq)

    @property
    def description(self):
        assert self._seq is not NULL
        return <bytes> libmeme.seq.get_seq_description(self._seq)

    @property
    def offset(self):
        assert self._seq is not NULL
        return libmeme.seq.get_seq_offset(self._seq)

    @offset.setter
    def offset(self, unsigned int offset):
        assert self._seq is not NULL
        libmeme.seq.set_seq_offset(offset, self._seq)

    @property
    def sequence(self):
        assert self._seq is not NULL

        cdef object mem = PyMemoryView_FromMemory(
            libmeme.seq.get_raw_sequence(self._seq),
            libmeme.seq.get_seq_length(self._seq) * sizeof(char),
            PyBUF_READ | PyBUF_WRITE
        )

        # DANGER (@althonos): We need to make sure that the object is not
        #                     deallocated before the memoryview, so we
        #                     register the `Sequence` object (self) as the
        #                     memoryview exporter. This needs a bit of
        #                     tweaking with the internal Py_buffer, because
        #                     `PyMemoryView_FromMemory` doesn't allow setting
        #                     up an exporter, and directly setting the `obj`
        #                     attribute is unsafe.
        PyBuffer_FillInfo(
            PyMemoryView_GET_BUFFER(mem),
            self,
            libmeme.seq.get_raw_sequence(self._seq),
            libmeme.seq.get_seq_length(self._seq) * sizeof(char),
            False,
            PyBUF_READ | PyBUF_WRITE,
        )

        return mem.cast('B')
