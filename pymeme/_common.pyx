# coding: utf-8
# cython: language_level=3, linetrace=True
"""Internal API common to all MEME tools.
"""

from cpython.buffer cimport PyBUF_FORMAT
from cpython.bytes cimport PyBytes_FromString, PyBytes_FromStringAndSize, PyBytes_AsString
from cpython.mem cimport PyMem_Free, PyMem_Malloc, PyMem_Realloc
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
cimport libmeme.meme
cimport libmeme.pssm
cimport libmeme.read_sequence
cimport libmeme.read_seq_file
from libmeme.alphabet cimport ALPH_T
from libmeme.hash_table cimport HASH_TABLE
from libmeme.data_types cimport WEIGHTS_T, Z_T, LCB_T
from libmeme.meme cimport CANDIDATE, DATASET, MODEL, OBJTYPE, SAMPLE, THETA, P_POINT, BRANCH_PARAMS, POINT_BRANCHES
from libmeme.mtype cimport MOTYPE
from libmeme.user cimport MINSITES, BFACTOR, HSIZE, HS_DECREASE

# --- Python imports ---------------------------------------------------------

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

    def __cinit__(self):
        self._alph = NULL

    def __dealloc__(self):
        libmeme.alphabet.alph_release(self._alph)

    def __repr__(self):
        assert self._alph != NULL
        if libmeme.alphabet.alph_is_builtin_rna(self._alph):
            return "Alphabet.rna()"
        elif libmeme.alphabet.alph_is_builtin_dna(self._alph):
            return "Alphabet.dna()"
        elif libmeme.alphabet.alph_is_builtin_protein(self._alph):
            return "Alphabet.protein()"
        else:
            return "Alphabet()"

    @property
    def size(self):
        """`int`: The number of core letters in the alphabet.
        """
        return libmeme.alphabet.alph_size_core(self._alph)

    @property
    def size_wild(self):
        """`int`: The number of core and wildcard letters in the alphabet.
        """
        return libmeme.alphabet.alph_size_wild(self._alph)


# --- Array ------------------------------------------------------------------

cdef class Array:
    """A 1D vector of fixed size with double-precision elements.
    """

    @classmethod
    def zeros(cls, int n):
        """zeros(cls, n)\n--

        Create a new array of size ``n`` filled with zeros.

        """
        if n < 0:
            raise ValueError("Cannot create a vector with negative size")

        cdef Array array = Array.__new__(Array)
        array._array = libmeme.array.allocate_array(n)
        if array._array is NULL:
            raise AllocationError("ARRAY_T", sizeof(ARRAY_T))
        with nogil:
            libmeme.array.init_array(0, array._array)

        return array

    def __cinit__(self):
        self._array = NULL

    def __dealloc__(self):
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


# --- Candidate model --------------------------------------------------------

cdef class Candidate:
    pass


# --- Dataset ----------------------------------------------------------------

cdef class Dataset:

    @classmethod
    def load(
        self,
        str filename,
        Alphabet alphabet,
        *,
        bint ignore_duplicates=True,
        bint use_complement=True,
        bint palindromes=False,
        int heapsize=HSIZE,
        int heapsize_decrease=HS_DECREASE
    ):
        cdef char* sample_name
        cdef char* sample_desc
        cdef char* sequence
        cdef long  length
        cdef bytes sample_name_py
        cdef bytes sample_desc_py
        cdef bytes sequence_py

        cdef FILE*   file      = fopen(os.fsencode(filename), "r")
        cdef Dataset dataset   = Dataset(
            alphabet,
            use_complement=use_complement,
            heapsize=heapsize,
            heapsize_decrease=heapsize_decrease
        )
        if file == NULL:
            raise OSError(errno, f"Failed to open file: {filename!r}")

        # create a hash table of sequence names (use a Python dict here)
        cdef set seq_names = set()

        try:
            # read all sequences in the file
            while libmeme.read_sequence.read_sequence(
                alphabet._alph,
                file,
                &sample_name,
                &sample_desc,
                &sequence,
                &length,
            ):
                # skip sequence if an error occured (FIXME?)
                if length < 0:
                    warnings.warn("Error encountered.")
                    continue

                # copy strings to Python heap
                sequence_py    = PyBytes_FromStringAndSize(sequence, length)
                sample_desc_py = PyBytes_FromString(sample_desc)
                sample_name_py = PyBytes_FromString(sample_name)
                free(sequence)
                free(sample_name)
                free(sample_desc)
                # parse weights if given (TODO)
                if sample_name_py == "WEIGHTS":
                    continue
                # ignore duplicates
                if sample_name_py in seq_names:
                    warnings.warn("Duplicate sequence name in input file.")
                    free(sample_name)
                    free(sample_desc)
                    free(sequence)
                    continue
                else:
                    seq_names.add(sample_name_py)
                # create the sample and add it to the dataset
                sample = Sample(alphabet, sequence_py, sample_name_py, sample_desc_py)
                dataset._append(sample)
        finally:
            fclose(file)

        # resize the array of samples
        if dataset._ds.n_samples > 0:
            dataset._ds.samples = <SAMPLE**> PyMem_Realloc(
                dataset._ds.samples,
                dataset._ds.n_samples * sizeof(SAMPLE*)
            )
            if dataset._ds.samples == NULL:
                raise AllocationError("SAMPLE*", sizeof(SAMPLE*), dataset._ds.n_samples)
        # initialize group sizes
        dataset._ds.n_region[0] = dataset._ds.max_slength
        dataset._ds.n_group = (dataset._ds.n_samples, 0, 0)
        # set the left flank region to be the entire sequence
        dataset._ds.region_last_pos[0] = dataset._ds.max_slength - 1
        dataset._ds.region_last_pos[1] = dataset._ds.max_slength - 1
        dataset._ds.region_last_pos[2] = dataset._ds.max_slength - 1
        # record the input order
        dataset._ds.input_order = dataset._ds.samples
        # record the filename
        dataset._ds.datafile = strdup(os.fsencode(filename))
        # set the database residue frequencies and entropy
        dataset._set_database_freq_and_entropy()
        # return the finalized dataset
        return dataset

    def __cinit__(self):
        self._ds = NULL
        self.alphabet = None
        self.samples = []

    def __init__(
        self,
        Alphabet alphabet,
        *,
        bint use_complement=True,
        bint palindromes=False,
        int heapsize=HSIZE,
        int heapsize_decrease=HS_DECREASE
    ):
        # make sure calling __init__ more than once doesn't cause a memory leak
        if self._ds == NULL:
            # allocate memory and check allocation was succesful
            self._ds = <DATASET*> allocate(sizeof(DATASET), "DATASET")
        else:
            # free any previous resource in use by the internal structure
            self._deallocate()
        # record alphabet and configuration from kwargs
        self.alphabet = alphabet
        self._ds.alph = libmeme.alphabet.alph_hold(alphabet._alph)
        # allocate internal buffers
        self._allocate()
        # reset counters
        self._ds.total_res = 0
        self._ds.wgt_total_res = 0.0
        self._ds.samples = NULL
        self._ds.input_order = NULL
        self._ds.n_samples = 0
        self._ds.n_group = (0, 0, 0)
        self._ds.seq_weights = NULL
        self._ds.group_last_idx = (0, 0, 0)
        self._ds.n_wgts = 0
        self._ds.max_slength = 0
        self._ds.min_slength = LONG_MAX
        self._ds.ce_frac = 0.0
        self._ds.n_region = (0, 0, 0)
        self._ds.region_last_pos = (0, 0, 0)
        self._ds.ce_max_dist = 0
        self._ds.psp_w = 0
        self._ds.log_psp_w = 0
        self._ds.mpi = False
        self._ds.invcomp = use_complement
        self._ds.pal = palindromes
        self._ds.bfile = NULL
        self._ds.print_pred = False
        self._ds.print_heaps = False
        self._ds.seed = 0
        self._ds.pspfile = NULL

    def __dealloc__(self):
        assert self._ds != NULL
        self._deallocate()
        free(self._ds)

    def __iter__(self):
        assert self._ds != NULL
        return iter(self.samples)

    def __len__(self):
        assert self._ds != NULL
        return self._ds.n_samples

    cdef void _allocate(self) except *:
        """allocate(self)\n--

        Allocate internal memory buffers for this dataset.

        """
        assert self._ds != NULL
        # clear memory
        # memset(self._ds, 0, sizeof(DATASET))
        # reset pointers
        self._ds.samples = NULL
        self._ds.input_order = NULL
        self._ds.seq_weights = NULL
        self._ds.res_freq = NULL
        self._ds.map = NULL
        self._ds.lomap = NULL

    cdef void _deallocate(self) except *:
        assert self._ds != NULL
        free(self._ds.samples)
        self._ds.samples = NULL
        # free(self._ds.input_order)
        self._ds.input_order = NULL
        free(self._ds.seq_weights)
        self._ds.seq_weights = NULL
        free(self._ds.res_freq)
        self._ds.res_freq = NULL
        matrix_free(<void**> self._ds.map)
        self._ds.map = NULL
        matrix_free(<void**> self._ds.lomap)
        self._ds.lomap = NULL

    cdef void _append(self, Sample sample) except *:
        assert self._ds != NULL

        # create a copy of the sample that we are going to edit and mark
        # its internal data as owned by the Dataset
        cdef Sample copy = sample.copy()

        # record the sample original index
        copy._sm.orig_index = self._ds.n_samples
        copy._sm.group = 0

        # record the sequence length
        if copy._sm.length > self._ds.max_slength:
            self._ds.max_slength = copy._sm.length
        if copy._sm.length < self._ds.min_slength:
            self._ds.min_slength = copy._sm.length

        # reallocate the sample array if too small
        if self._ds.n_samples % libmeme.read_sequence.RCHUNK == 0:
            self._ds.samples = <SAMPLE**> PyMem_Realloc(
                self._ds.samples,
                (self._ds.n_samples + libmeme.read_sequence.RCHUNK) * sizeof(SAMPLE*)
            )
            if self._ds.samples == NULL:
                raise AllocationError("SAMPLE*", sizeof(SAMPLE*), (self._ds.n_samples + libmeme.read_sequence.RCHUNK) * sizeof(SAMPLE*))

        # record the sample
        self._ds.samples[self._ds.n_samples] = copy._sm
        self._ds.n_samples += 1
        self.samples.append(copy)

    cdef void _set_database_freq_and_entropy(self):
        cdef int     i
        cdef int     j
        cdef int     k
        cdef long    slen
        cdef double  sw
        cdef double  core_count
        cdef double* seq_weights = self._ds.seq_weights
        cdef int     n_wgts      = self._ds.n_wgts
        cdef ALPH_T* alph        = self._ds.alph

        # reset the frequency array
        free(self._ds.res_freq)
        # count residues
        self._ds.res_freq = <double*> allocate(self.alphabet.size * sizeof(double), "double*")
        for i in range(self.alphabet.size):
            self._ds.res_freq[i] = 0
        self._ds.total_res = 0
        self._ds.wgt_total_res = 0
        for i in range(self._ds.n_samples):
            slen = self._ds.samples[i].length
            if n_wgts > i:
                sw = self._ds.samples[i].sw = seq_weights[i]
            else:
                sw = self._ds.samples[i].sw = 1
            self._ds.total_res += slen
            core_count = 0
            for j in range(self.alphabet.size):
                core_count += self._ds.samples[i].counts[j]
                if self._ds.invcomp:
                    k = libmeme.alphabet.alph_complement(self.alphabet._alph, j)
                    self._ds.res_freq[j] += sw * self._ds.samples[i].counts[j]/2.0
                    self._ds.res_freq[k] += sw * self._ds.samples[i].counts[j]/2.0
                else:
                    self._ds.res_freq[j] += sw * self._ds.samples[i].counts[j]

            self._ds.wgt_total_res += sw * core_count

        # convert counts to frequencies */
        for i in range(self.alphabet.size):
            self._ds.res_freq[i] /= self._ds.wgt_total_res

    cpdef void shuffle(self) except *:
        assert self._ds != NULL
        libmeme.read_seq_file.shuffle_dataset_order(self._ds)


# --- Model ------------------------------------------------------------------

cdef class Model:

    def __cinit__(self):
        self._mm = NULL
        self._owner = None
        self.alphabet = None

    def __init__(
        self,
        Alphabet alphabet not None,
        unsigned int maximum_width,
        *,
        str model_type = "zoops",
        str objective_function = "classic",
        bint use_complement=True,
    ):
        assert alphabet._alph != NULL

        # make sure the objective function is valid
        if objective_function not in MEME_OBJECTIVE_FUNCTIONS:
            raise ValueError(f"Unknown or invalid objective function: {objective_function!r}")
        if model_type not in MEME_MODEL_TYPES:
            raise ValueError(f"Unknown or invalid model type: {model_type!r}")
        # make sure alphabet can be complemented
        if not libmeme.alphabet.alph_has_complement(alphabet._alph):
            if use_complement:
                raise ValueError("A complementable alphabet must be provided if `use_complement=True`")

        # allow calling __init__ more than once without leaking
        if self._mm != NULL:
            libmeme.meme.free_model(self._mm)
        # allocate a new model
        self._mm = libmeme.meme.create_model(
            MEME_MODEL_TYPES[model_type],
            use_complement,
            maximum_width,
            alphabet._alph,
            MEME_OBJECTIVE_FUNCTIONS[objective_function],
        )
        if self._mm == NULL:
            raise AllocationError("MODEL", sizeof(MODEL))

        # keep a reference to the alphabet
        self.alphabet = alphabet

    def __dealloc__(self):
        assert self._mm != NULL
        # if self._owner is None:
        #     libmeme.meme.free_model(self._mm)

    # cdef void _allocate(self) except *:
    #     assert self._mm != NULL
    #     # clear old pointers
    #     self._mm.theta = NULL
    #     self._mm.logtheta = NULL
    #     self._mm.logtheta_rc = NULL
    #     self._mm.obs = NULL
    #     # allocate new buffers
    #     self._mm.theta = <THETA> matrix_create(self._mm.max_w + 1, self.alphabet.size_wild, sizeof(double))
    #     self._mm.logtheta = <THETA> matrix_create(self._mm.max_w + 1, self.alphabet.size_wild, sizeof(double))
    #     self._mm.logtheta_rc = <THETA> matrix_create(self._mm.max_w + 1, self.alphabet.size_wild, sizeof(double))
    #     self._mm.obs = <THETA> matrix_create(self._mm.max_w + 1, self.alphabet.size, sizeof(double))
    #
    # cdef void _deallocate(self) except *:
    #     assert self._mm != NULL
    #     matrix_free(<void**> self._mm.theta)
    #     self._mm.theta = NULL
    #     matrix_free(<void**> self._mm.logtheta)
    #     self._mm.logtheta = NULL
    #     matrix_free(<void**> self._mm.logtheta_rc)
    #     self._mm.logtheta_rc = NULL
    #     matrix_free(<void**> self._mm.obs)
    #     self._mm.obs = NULL

    cpdef Model copy(self):
        """copy(self)\n--

        Create a copy of the model with the same parameters.

        """
        assert self._mm != NULL
        # create new model and allocate internal memory
        cdef Model new = Model.__new__(Model)
        new._mm = libmeme.meme.create_model(
            self._mm.mtype,
            self._mm.invcomp,
            self._mm.max_w,
            self.alphabet._alph,
            self._mm.objfun,
        )
        # store reference to alphabet
        new.alphabet = self.alphabet
        # copy attributes
        libmeme.meme.copy_model(self._mm, new._mm, self.alphabet._alph)
        # return the copy
        return new


# --- Motif ------------------------------------------------------------------

cdef class Motif:

    def __cinit__(self):
        self._motif = NULL
        self.alphabet = None

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

        cdef PSSM pssm = PSSM.__new__(PSSM)
        pssm.motif = self
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

        return pssm


# --- Motif summary ----------------------------------------------------------

cdef class MotifSummaries:

    def __cinit__(self):
        self.alphabet = None
        self._motif_summaries = NULL
        self._nmotifs = 0

    def __dealloc__(self):
        cdef size_t i
        for i in range(self._nmotifs):
            # FIXME
            # libmeme.macros.free_2array(<void**> self._motif_summaries[i].pssm, self._motif_summaries[i].width + 1)
            # libmeme.macros.free_2array(<void**> self._motif_summaries[i].psfm, self._motif_summaries[i].width + 1)
            # free(self._motif_summaries[i].sites)
            # free(self._motif_summaries[i].regexp)
            # free(self._motif_summaries[i].consensus)
            pass
        free(self._motif_summaries)

    def __len__(self):
        return self._nmotifs

    def __getitem__(self, ssize_t index):
        cdef MotifSummary ms

        if index < 0:
            index += self._nmotifs
        if index >= self._nmotifs or index < 0:
            raise IndexError("list index out of range")

        ms = MotifSummary.__new__(MotifSummary)
        ms._owner = self
        ms._ms = &self._motif_summaries[index]
        ms.alphabet = self.alphabet
        return ms


cdef class MotifSummary:

    def __cinit__(self):
        self._ms = NULL
        self._owner = None
        self.alphabet = None

    @property
    def width(self):
        """`int`: The width of the motif.
        """
        assert self._ms != NULL
        return self._ms.width

    @property
    def num_sites(self):
        """`int`: The number of sites of the motif.
        """
        assert self._ms != NULL
        return self._ms.num_sites

    @property
    def information_content(self):
        """`float`: The information content of the motif.
        """
        assert self._ms != NULL
        return self._ms.ic

    @property
    def relative_entropy(self):
        """`float`: The relative entropy of the motif.
        """
        assert self._ms != NULL
        return self._ms.re

    @property
    def e_value(self):
        """`float`: The E-value of the motif.
        """
        assert self._ms != NULL
        return self._ms.e_value_mant * 10 ** self._ms.e_value_exp

    @property
    def consensus(self):
        """`str`: The motif as a IUPAC consensus sequence.
        """
        assert self._ms != NULL
        return self._ms.consensus.decode("ascii")

    @property
    def regexp(self):
        """`str`: The motif as a regular expression
        """
        assert self._ms != NULL
        return self._ms.regexp.decode("ascii")


# --- PriorDistribution ------------------------------------------------------

cdef class PriorDist:

    def __cinit__(self):
        self._pd = NULL


# --- PSSM -------------------------------------------------------------------

cdef class PSSM:

    def __cinit__(self):
        self.motif = None
        self._pssm = NULL

    def __dealloc__(self):
        libmeme.pssm.free_pssm(self._pssm)

    @property
    def width(self):
        assert self._pssm is not NULL
        return libmeme.pssm.get_pssm_w(self._pssm)


# --- Sample -----------------------------------------------------------------

cdef class Sample:

    def __cinit__(self):
        self._sm = NULL
        self._owner = None
        self.alphabet = None
        self.use_complement = True

    def __init__(
        self,
        Alphabet alphabet,
        bytes sequence,
        bytes name,
        bytes description,
        bint use_complement = True,
    ):
        cdef int    e
        cdef size_t i
        cdef size_t length = len(sequence)

        # make sure the sequence is not empty
        if len(sequence) <= 0:
            raise ValueError("Cannot create a zero-length sample sequence.")

        # make sure calling __init__ more than once doesn't cause a memory leak
        if self._sm == NULL:
            # allocate memory and make sure pointers are set to NULL to avoid a double free
            self._sm = <SAMPLE*> allocate(sizeof(SAMPLE), "SAMPLE")
        else:
            self._deallocate()

        # mark the sample as standalone (memory is not owned by some other object)
        self._owner = None
        # record parameters that will influence allocation
        self._sm.length = length
        self.alphabet = alphabet
        # allocate internal numeric buffers
        self._allocate()
        # copy sequence
        memcpy(self._sm.seq, PyBytes_AsString(sequence), length)
        self._sm.seq[length] = b'\0'
        # fill internal numeric buffers from the sequence contents
        self._encode_sequence()
        self._compute_sequence_weights()
        self._count_residues()
        # reset housekeeping
        self._sm.max_log_psp = -INFINITY
        self._sm.nsites = 0
        # store metadata
        self.name = name
        self.description = description
        self.use_complement = use_complement

    def __dealloc__(self):
        assert self._sm != NULL
        if self._owner is None:
            self._deallocate()
            free(self._sm)

    @property
    def name(self):
        assert self._sm != NULL
        return <bytes> self._sm.sample_name

    @name.setter
    def name(self, bytes name):
        assert self._sm != NULL
        # free previous resource if there is any
        if self._sm.sample_name != NULL:
            free(self._sm.sample_name)
            self._sm.sample_name = NULL
        # copy new string
        self._sm.sample_name = <char*> allocate((len(name) + 1) * sizeof(char), "char*")
        strcpy(self._sm.sample_name, PyBytes_AsString(name))

    @property
    def description(self):
        assert self._sm != NULL
        return <bytes> self._sm.descript

    @description.setter
    def description(self, bytes desc):
        assert self._sm != NULL
        # free previous resource if there is any
        if self._sm.descript != NULL:
            free(self._sm.descript)
            self._sm.descript = NULL
        # copy new string
        self._sm.descript = <char*> allocate((len(desc) + 1) * sizeof(char), "char*")
        strcpy(self._sm.descript, PyBytes_AsString(desc))

    @property
    def sequence(self):
        assert self._sm != NULL
        return PyBytes_FromStringAndSize(self._sm.seq, self._sm.length)

    cdef void _allocate(self) except *:
        """allocate(self)\n--

        Allocate internal memory buffers for this sample.

        """
        assert self._sm != NULL

        cdef size_t length    = <size_t> self._sm.length
        cdef size_t core_size = <size_t> self.alphabet.size

        # set all pointers to NULL in case we encounter an exception while
        # not all arrays have been allocated
        self._sm.sample_name = NULL
        self._sm.descript = NULL
        self._sm.seq = NULL
        self._sm.res = NULL
        self._sm.resic = NULL
        self._sm.weights = NULL
        self._sm.not_o = NULL
        self._sm.log_not_o = NULL
        self._sm.pY = NULL
        self._sm.pYic = NULL
        self._sm.z_buf = self._sm.z = NULL
        self._sm.counts = NULL
        self._sm.logcumback = NULL
        self._sm.log_psp_buf = self._sm.log_psp = NULL
        self._sm.psp_original_buf = self._sm.psp_original = NULL
        self._sm.sites = NULL

        # allocate sequence
        self._sm.seq = <char*> allocate((length + 1) * sizeof(char), "char*")
        # allocate residues
        self._sm.res = <uint8_t*> allocate(length * sizeof(uint8_t), "uint8_t*")
        if self.use_complement:
            self._sm.resic = <uint8_t*> allocate(length * sizeof(uint8_t), "uint8_t*")
        # allocate weights
        self._sm.weights = <WEIGHTS_T*> allocate(length * sizeof(WEIGHTS_T), "WEIGHTS_T*")
        # allocate erasing arrays
        self._sm.not_o = <double*> allocate(length * sizeof(double), "double*")
        self._sm.log_not_o = <int*> allocate(length * sizeof(int), "int*")
        # allocate set up arrays to hold posterior probabilites
        self._sm.pY = <int**> matrix_create(3, length, sizeof(int))
        self._sm.pYic = <char*> allocate(length * sizeof(char), "char*")
        # allocate set up Z array and log_psp arrays
        if self.use_complement:
            self._sm.z_buf = <Z_T*> allocate(sizeof(Z_T) * (2 * length + 1), "Z_T*")
            self._sm.z = self._sm.z_buf + length
        else:
            self._sm.z_buf = <Z_T*> allocate(sizeof(Z_T) * (length + 1), "Z_T*")
            self._sm.z = self._sm.z_buf
        # allocate array to hold character counts and get character counts
        self._sm.counts = <double*> allocate(core_size * sizeof(double), "double*")
        # allocate up remaining arrays
        self._sm.logcumback = <LCB_T*> allocate((length+1) * sizeof(LCB_T), "LCB_T*")

    cdef void _encode_sequence(self) except *:
        """_encode_sequence(self)\n--

        Encode the text sequence to digits using ``self.alphabet``.

        """
        assert self._sm != NULL
        assert self._sm.seq != NULL
        assert self._sm.res != NULL
        assert self.alphabet._alph != NULL

        cdef size_t i
        cdef int    e
        cdef size_t length = <size_t> self._sm.length

        for i in range(length):
            e = libmeme.alphabet.alph_encodec(self.alphabet._alph, self._sm.seq[i])
            self._sm.res[i] = e
            if self.use_complement:
                self._sm.resic[length - i - 1] = libmeme.alphabet.alph_complement(self.alphabet._alph, e)

    cdef void _compute_sequence_weights(self) except *:
        """_compute_sequence_weights(self)\n--

        Compute the sequence weights.

        """
        assert self._sm != NULL
        assert self._sm.seq != NULL
        assert self._sm.res != NULL
        assert self.alphabet._alph != NULL

        cdef long i

        for i in range(self._sm.length):
            self._sm.weights[i] = 1.0

    cdef void _count_residues(self) except *:
        """_count_residues(self)\n--

        Count the residues in the input sequence.

        """
        assert self._sm != NULL
        assert self._sm.res != NULL
        assert self._sm.counts != NULL
        assert self.alphabet._alph != NULL

        cdef long i
        cdef int  e

        # reset counts to zero
        for i in range(self.alphabet.size):
            self._sm.counts[i] = 0.0

        # count residues
        for i in range(self._sm.length):
            e = self._sm.res[i]
            if e != libmeme.alphabet.alph_wild(self.alphabet._alph):
                self._sm.counts[e] += 1

    cdef void _deallocate(self) except *:
        """_deallocate(self)\n--

        Free memory allocate for the internal buffers of this sample.

        """
        assert self._sm != NULL
        free(self._sm.sample_name)
        self._sm.sample_name = NULL
        free(self._sm.descript)
        self._sm.descript = NULL
        free(self._sm.seq)
        self._sm.seq = NULL
        free(self._sm.res)
        self._sm.res = NULL
        free(self._sm.resic)
        self._sm.resic = NULL
        free(self._sm.weights)
        self._sm.weights = NULL
        free(self._sm.not_o)
        self._sm.not_o = NULL
        free(self._sm.log_not_o)
        self._sm.log_not_o = NULL
        matrix_free(<void**> self._sm.pY)
        self._sm.pY = NULL
        free(self._sm.pYic)
        self._sm.pYic = NULL
        free(self._sm.z_buf)
        self._sm.z_buf = self._sm.z = NULL
        free(self._sm.psp_original_buf)
        self._sm.psp_original = self._sm.psp_original_buf = NULL
        free(self._sm.log_psp_buf)
        self._sm.log_psp_buf = self._sm.log_psp = NULL
        free(self._sm.counts)
        self._sm.counts = NULL
        free(self._sm.logcumback)
        self._sm.logcumback = NULL
        free(self._sm.sites)
        self._sm.sites = NULL

    # cdef void _copy_to(self, SAMPLE* dest) except *:
    #     """_copy_to(self, SAMPLE* dest)\n--
    #
    #     Copy this sample to a
    #
    #     """
    #     assert self._sm != NULL
    #     # copy name
    #     dest.sample_name = strdup(self._sm.sample_name)
    #     if dest.sample_name == NULL:
    #         raise AllocationError("char*")
    #     # copy description
    #     dest.descript = strdup(self._sm.descript)
    #     if dest.descript == NULL:
    #         raise AllocationError("char*")
    #     # copy sequence
    #     dest.length = self._sm.length
    #     dest.seq = strdup(self._sm.seq)
    #     if dest.seq == NULL:
    #         raise AllocationError("char*")
    #     # copy residues
    #     dest.res = <uint8_t*> allocate(self._sm.length * sizeof(uint8_t), "uint8_t*")
    #     memcpy(dest.res, self.res, length * sizeof(uint8_t))
    #     # copy reverse complement
    #     dest.resic = <uint8_t*> allocate(self._sm.length * sizeof(uint8_t), "uint8_t*")
    #     memcpy(dest.res, self.res, length * sizeof(uint8_t))
    #     # copy sequence weight
    #     dest.sw = self._sm.sw
    #     dest.weights = <WEIGHTS_T*> allocate(self._sm.length * sizeof(WEIGHTS_T), "WEIGHTS_T*")
    #     memcpy(dest.weights, self._sm.weights, self._sm.length * sizeof(WEIGHTS_T))
    #     # copy erasing arrays
    #     dest.not_o = <double*> allocate(self._sm.length * sizeof(double), "double*")
    #     memcpy(self._sm )

    cpdef Sample copy(self):
        """copy(self)\n--

        Return an exact duplicate of the sample.

        """
        assert self._sm != NULL
        return Sample(
            self.alphabet,
            self.sequence,
            self.name,
            self.description,
            self._sm.resic != NULL,
        )
