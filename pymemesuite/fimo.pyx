# coding: utf-8
# cython: language_level=3, linetrace=True
"""Implementation of the Find Individual Motif Occurences (FIMO) method.

FIMO searches a set of sequences for occurrences of known motifs, treating
each motif independently. It uses a dynamic programming algorithm to convert
log-odds scores into *p-values*, assuming a zero-order background model, and
then into *q-values* following the method of Benjamini & Hochberg.

See Also:
    `Grant, Charles E., Timothy L. Bailey, and William Stafford Noble.
    ‘FIMO: Scanning for Occurrences of a given Motif’. Bioinformatics 27, no. 7
    (1 April 2011): 1017–18. <https://doi.org/10.1093/bioinformatics/btr064>`_.

"""

from cpython.pythread cimport (
    PyThread_type_lock,
    PyThread_allocate_lock,
    PyThread_free_lock,
    PyThread_acquire_lock,
    PyThread_release_lock,
    WAIT_LOCK,
    PY_LOCK_FAILURE,
)

from libc.math cimport NAN, isnan, log2
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.stdint cimport int8_t
from libc.stdio cimport printf

cimport libmeme.alphabet
cimport libmeme.array
cimport libmeme.cisml
cimport libmeme.matrix
cimport libmeme.motif
cimport libmeme.pssm
cimport libmeme.reservoir
cimport libmeme.seq
from libmeme.alphabet cimport ALPH_T
from libmeme.array cimport ARRAY_T
from libmeme.matrix cimport MATRIX_T
from libmeme.motif cimport MOTIF_T
from libmeme.pssm cimport PSSM_RANGE, PSSM_T
from libmeme.reservoir cimport RESERVOIR_SAMPLER_T
from libmeme.seq cimport SEQ_T
from libmeme.cisml cimport PATTERN_T, SCANNED_SEQUENCE_T, MATCHED_ELEMENT_T

from .common cimport Array, Background, Motif, PSSM, Sequence, ReservoirSampler, PriorDist
from .cisml cimport Pattern, MatchedElement, ScannedSequence

# --- Python imports ---------------------------------------------------------

import warnings

from .errors import AllocationError


cdef extern from "alloca.h" nogil:
    cdef void* alloca(size_t size)


cdef struct site_score_t:
    double score
    double pvalue
    bint   scoreable


cdef class FIMO:

    cdef double       alpha
    cdef double       threshold
    cdef bint         both_strands
    cdef unsigned int max_stored_scores

    def __init__(
        self,
        *,
        double alpha = 1.0,
        bint both_strands = True,
        double threshold = 1e-4,
        int max_stored_scores = 100000,
    ):
        """__init__(self, *, alpha=1.0, both_strands=True, threshold=1e-4, max_stored_scores=100000)\n--

        Create a new FIMO runner.

        Keyword Arguments:
            alpha (`float`): The scale factor for non-specific priors,
            both_strands (`bool`): Whether or not to run the motif scan on
                both segments.
            threshold (`float`): The p-value above which matches are
                marked as non-significant and discarded.
            max_scored_scores (`int`): The maximum number of scored sites
                to keep in the result `Pattern`. Using a smaller number
                will reduce memory consumption at the cost of statistical
                accuracy, and significant results may be lost.

        """
        self.alpha = alpha
        self.both_strands = both_strands
        self.threshold = threshold
        self.max_stored_scores = max_stored_scores

    cdef bint _record_score(
        self,
        PATTERN_T* pattern,
        SCANNED_SEQUENCE_T* scanned_seq,
        RESERVOIR_SAMPLER_T* reservoir,
        MATCHED_ELEMENT_T* match,
    ) noexcept nogil:
        cdef double pvalue   = libmeme.cisml.get_matched_element_pvalue(match)
        cdef bint   recorded = False

        libmeme.cisml.add_scanned_sequence_scanned_position(scanned_seq)
        libmeme.reservoir.reservoir_sample(reservoir, pvalue)

        if pvalue < self.threshold:
            recorded = libmeme.cisml.add_pattern_matched_element(pattern, match)

        return recorded

    cdef site_score_t _score_site(
        self,
        const int8_t* sequence,
        PSSM_T* pssm,
        double prior,
    ) noexcept nogil:
        cdef site_score_t site_score
        cdef int          motif_position
        cdef double       prior_log_odds
        cdef double       score
        cdef double       pvalue
        cdef int8_t       aindex
        cdef int          w               = pssm.w
        cdef ARRAY_T*     pv_lookup       = pssm.pv
        cdef MATRIX_T*    pssm_matrix     = pssm.matrix
        cdef ALPH_T*      alphabet        = pssm.alph
        cdef double       scaled_log_odds = 0.0

        site_score.scoreable = True
        for motif_position in range(w):
            aindex = sequence[motif_position]
            if aindex == -1:
                site_score.scoreable = False
                break
            scaled_log_odds += libmeme.matrix.get_matrix_cell(motif_position, aindex, pssm_matrix)

        if site_score.scoreable:

            score = libmeme.pssm.get_unscaled_pssm_score(scaled_log_odds, pssm)
            if not isnan(prior):
                w += 1
                prior_log_odds = self.alpha * prior
                prior_log_odds = log2(prior_log_odds / (1.0 - prior_log_odds))
                score += prior_log_odds
                scaled_log_odds = libmeme.pssm.raw_to_scaled(score, w, pssm.scale, pssm.offset)

            max_log_odds = libmeme.array.get_array_length(pv_lookup) - 1
            if scaled_log_odds < 0.0:
                # warnings.warn(f"Scaled log-odds out of range ({scaled_log_odds}), using {0}")
                scaled_log_odds = 0.0
            elif scaled_log_odds > max_log_odds:
                # warnings.warn(f"Scaled log-odds out of range ({scaled_log_odds}), using {max_log_odds}")
                scaled_log_odds = max_log_odds

            site_score.score = libmeme.pssm.scaled_to_raw(scaled_log_odds, w, pssm.scale, pssm.offset)
            site_score.pvalue = libmeme.array.get_array_item(<int> scaled_log_odds, pv_lookup)

        return site_score

    cdef int _score_sequence(
        self,
        RESERVOIR_SAMPLER_T* reservoir,
        PATTERN_T* pattern,
        SEQ_T* seq,
        PSSM_T* pssm,
        PSSM_T* pssm_rev,
        int8_t** buffer,
        size_t*  buflen,
    ) except 1 nogil:
        cdef int                 offset
        cdef int                 start_bwd
        cdef int                 start_fwd
        cdef int                 stop_bwd
        cdef int                 stop_fwd
        cdef site_score_t        scores_fwd
        cdef site_score_t        scores_bwd
        cdef SCANNED_SEQUENCE_T* scanned_seq = NULL
        cdef MATCHED_ELEMENT_T*  match_fwd   = NULL
        cdef MATCHED_ELEMENT_T*  match_bwd   = NULL
        cdef ALPH_T*             alphabet    = libmeme.pssm.get_pssm_alph(pssm)
        cdef unsigned char*      seqdata     = <unsigned char*> libmeme.seq.get_raw_sequence(seq)
        cdef int                 length      = libmeme.seq.get_seq_length(seq)
        cdef int                 width       = pssm.w

        # create a new CisML scanned sequence to store results
        scanned_seq = libmeme.cisml.allocate_scanned_sequence(
            libmeme.seq.get_seq_name(seq),
            libmeme.seq.get_seq_name(seq),
            pattern
        )

        # use the temporary buffer to index the sequence
        if buflen[0] < length:
            buffer[0] = <int8_t*> realloc(buffer[0], length * sizeof(int8_t))
            if buffer[0] == NULL:
                raise AllocationError("uint8_t", sizeof(int8_t), length)
            buflen[0] = length
        for offset in range(length):
            buffer[0][offset] = libmeme.alphabet.alph_indexc(alphabet, seqdata[offset])

        # iterate the sliding window
        for offset in range(length - width):
            # compute block coordinates
            stop_bwd = start_fwd = offset + 1
            stop_fwd = start_bwd = start_fwd + width - 1
            # recycle match memory if possible
            if match_fwd is NULL:
                match_fwd = libmeme.cisml.allocate_matched_element(start_fwd, stop_fwd, scanned_seq)
            else:
                libmeme.cisml.set_matched_element_start(match_fwd, start_fwd)
                libmeme.cisml.set_matched_element_stop(match_fwd, stop_fwd)
            # match motif on forward strand
            scores_fwd = self._score_site(&buffer[0][offset], pssm, NAN)
            if scores_fwd.scoreable:
                libmeme.cisml.set_matched_element_pvalue(match_fwd, scores_fwd.pvalue)
                if self._record_score(pattern, scanned_seq, reservoir, match_fwd):
                    libmeme.cisml.set_matched_element_sequence(match_fwd, <char*> &seqdata[offset])
                    libmeme.cisml.set_matched_element_strand(match_fwd, b'+')
                    libmeme.cisml.set_matched_element_score(match_fwd, scores_fwd.score)
                    match_fwd = NULL
            # if reverse strand is supported, match motif on reverse strand
            if pssm_rev is not NULL:
                # recycle match memory if possible
                if match_bwd is NULL:
                    match_bwd = libmeme.cisml.allocate_matched_element(start_bwd, stop_bwd, scanned_seq)
                else:
                    libmeme.cisml.set_matched_element_start(match_bwd, start_bwd)
                    libmeme.cisml.set_matched_element_stop(match_bwd, stop_bwd)
                # match motif on reverse strand
                scores_bwd = self._score_site(&buffer[0][offset], pssm_rev, NAN)
                if scores_bwd.scoreable:
                    libmeme.cisml.set_matched_element_pvalue(match_bwd, scores_bwd.pvalue)
                    if self._record_score(pattern, scanned_seq, reservoir, match_bwd):
                        libmeme.cisml.set_matched_element_sequence(match_bwd, <char*> &seqdata[offset])
                        libmeme.cisml.set_matched_element_strand(match_bwd, b'-')
                        libmeme.cisml.set_matched_element_score(match_bwd, scores_bwd.score)
                        # DANGER(@althonos): This bit is completely unsafe, as we exploit
                        #                    the fact that `set_matched_element_sequence`
                        #                    already allocated a buffer where we can
                        #                    reverse-complement inplace
                        libmeme.alphabet.invcomp_seq(
                            pssm.alph,
                            <char*> libmeme.cisml.get_matched_element_sequence(match_bwd),
                            width
                        )
                        match_bwd = NULL

        # free matched elements that were not recorded
        if match_fwd is not NULL:
            libmeme.cisml.free_matched_element(match_fwd)
        if match_bwd is not NULL:
            libmeme.cisml.free_matched_element(match_bwd)

        return 0

    cpdef Pattern score_pssm(
        self,
        PSSM pssm,
        list sequences,
    ):
        """score_pssm(self, pssm, sequences)\n--

        Score sequences with a position-specific scoring matrix.

        Arguments:
            pssm (`~pymemesuite.common.PSSM`): The position-specific scoring
                matrix to score the sequences with.
            sequences (`list` of `~pymemesuite.common.Sequence`): A list
                containing the target sequences.

        Returns:
            `~pymemesuite.cisml.Pattern`: A pattern object storing
            significant matches of the PSSM to the sequences.

        """
        cdef int              i
        cdef Sequence         sequence
        cdef ARRAY_T          values
        cdef PSSM             pssm_rev     = None
        cdef ReservoirSampler reservoir    = ReservoirSampler(10000)
        cdef Pattern          pattern      = Pattern(pssm.motif.accession, pssm.motif.name)
        cdef SEQ_T**          seqptr       = NULL
        cdef int              length       = len(sequences)
        cdef PSSM_T*          pssm_raw_fwd = pssm._pssm
        cdef PSSM_T*          pssm_raw_bwd = NULL
        cdef int8_t*          seqbuffer    = NULL
        cdef size_t           seqbuflen    = 0

        # store sequences in contiguous array to process them without the GIL
        seqptr = <SEQ_T**> calloc(length, sizeof(SEQ_T*))
        if seqptr is NULL:
            raise AllocationError("SEQ_T*", sizeof(SEQ_T*), length)
        for i, sequence in enumerate(sequences):
            seqptr[i] = sequence._seq

        try:
            # record options in the pattern object
            libmeme.cisml.set_pattern_max_stored_matches(pattern._pattern, self.max_stored_scores)
            libmeme.cisml.set_pattern_max_pvalue_retained(pattern._pattern, self.threshold)
            # compute reverse complement
            if self.both_strands:
                pssm_rev = pssm.reverse_complement()
                pssm_raw_bwd = pssm_rev._pssm
            with nogil:
                # score all sequences with the given motif
                for i in range(length):
                    self._score_sequence(
                        reservoir._reservoir,
                        pattern._pattern,
                        seqptr[i],
                        pssm_raw_fwd,
                        pssm_raw_bwd,
                        &seqbuffer,
                        &seqbuflen,
                    )
                # complete pattern
                libmeme.cisml.set_pattern_is_complete(pattern._pattern)
                # compute q-values
                values.num_items = libmeme.reservoir.get_reservoir_num_samples_retained(reservoir._reservoir)
                values.items = libmeme.reservoir.get_reservoir_samples(reservoir._reservoir)
                libmeme.cisml.pattern_calculate_qvalues(pattern._pattern, &values)
        finally:
            free(seqptr)
            free(seqbuffer)
    
        return pattern

    cpdef Pattern score_motif(
        self,
        Motif motif,
        list sequences,
        Background background,
        Background background_p = None,
        PriorDist prior_dist = None,
    ):
        """score_motif(self, pssm, sequences)\n--

        Score sequences with a motif.

        The motif is first converted into a PSSM with `Motif.build_pssm`,
        and then sequences are scored with the `~FIMO.score_pssm` method.

        Arguments:
            Motif (`~pymemesuite.common.Motif`): The MEME motif to score
                the sequences with.
            sequences (`list` of `~pymemesuite.common.Sequence`): A list
                containing the target sequences.
            background (`~pymemesuite.common.Background`): The background
                letter frequencies for building the odds ratio for each
                position of the PSSM.
            background_p (`~pymemesuite.common.Background`, optional): The
                background letter frequencies for building the p-value
                lookup table of the PSSM. If `None`, use the `background`
                values.
            prior_dist (`~pymemesuite.common.PriorDist`, optional): A
                distribution of priors for building the PSSM.

        Returns:
            `~pymemesuite.cisml.Pattern`: A pattern object storing
            significant matches of the PSSM to the sequences.

        """
        # build PSSMs from motif
        pssm = motif.build_pssm(
            background,
            background_p,
            prior_dist,
            self.alpha,
            PSSM_RANGE,
            0,
            False
        )
        # score PSSM
        return self.score_pssm(pssm, sequences)
