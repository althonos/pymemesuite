# coding: utf-8
# cython: language_level=3, linetrace=True
"""Implementation of the Find Individual Motif Occurences (FIMO) method.

FIMO searches a set of sequences for occurrences of known motifs, treating
each motif independently. It uses a dynamic programming algorithm to convert
log-odds scores into *p-values*, assuming a zero-order background model, and
then into *q-values* following the method of Benjamini & Hochberg.

"""

from libc.math cimport NAN, isnan, log2

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

from .common cimport Array, Motif, PSSM, Sequence, ReservoirSampler, PriorDist
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

    cdef double alpha
    cdef double threshold
    cdef bint   both_strands

    def __init__(
        self,
        *,
        double alpha = 1.0,
        bint both_strands = True,
        double threshold = 1e-4,
    ):
        self.alpha = alpha
        self.both_strands = both_strands
        self.threshold = threshold

    cdef bint _record_score(
        self,
        PATTERN_T* pattern,
        SCANNED_SEQUENCE_T* scanned_seq,
        RESERVOIR_SAMPLER_T* reservoir,
        MATCHED_ELEMENT_T* match,
    ) nogil:
        cdef double pvalue   = libmeme.cisml.get_matched_element_pvalue(match)
        cdef bint   recorded = False

        libmeme.cisml.add_scanned_sequence_scanned_position(scanned_seq)
        libmeme.reservoir.reservoir_sample(reservoir, pvalue)

        if pvalue < self.threshold:
            recorded = libmeme.cisml.add_pattern_matched_element(pattern, match)

        return recorded

    cdef site_score_t _score_site(
        self,
        const unsigned char* sequence,
        PSSM_T* pssm,
        double prior,
    ) nogil:
        cdef site_score_t site_score
        cdef int          motif_position
        cdef double       prior_log_odds
        cdef double       score
        cdef double       pvalue
        cdef int          w               = pssm.w
        cdef ARRAY_T*     pv_lookup       = pssm.pv
        cdef MATRIX_T*    pssm_matrix     = pssm.matrix
        cdef ALPH_T*      alphabet        = pssm.alph
        cdef double       scaled_log_odds = 0.0

        site_score.scoreable = True
        for motif_position in range(w):
            c = sequence[motif_position]
            aindex = libmeme.alphabet.alph_indexc(alphabet, c)
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

    cdef void _score_sequence(
        self,
        RESERVOIR_SAMPLER_T* reservoir,
        PATTERN_T* pattern,
        SEQ_T* seq,
        PSSM_T* pssm,
        PSSM_T* pssm_rev,
    ) nogil:
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
        cdef unsigned char*      seqdata     = <unsigned char*> libmeme.seq.get_raw_sequence(seq)
        cdef int                 length      = libmeme.seq.get_seq_length(seq)
        cdef int                 width       = pssm.w

        scanned_seq = libmeme.cisml.allocate_scanned_sequence(
            libmeme.seq.get_seq_name(seq),
            libmeme.seq.get_seq_name(seq),
            pattern
        )

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
            scores_fwd = self._score_site(&seqdata[offset], pssm, NAN)
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
                match_bwd = libmeme.cisml.allocate_matched_element(start_fwd, stop_fwd, scanned_seq)
                scores_bwd = self._score_site(&seqdata[offset], pssm_rev, NAN)
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

        if match_fwd is not NULL:
            libmeme.cisml.free_matched_element(match_fwd)
        if match_bwd is not NULL:
            libmeme.cisml.free_matched_element(match_bwd)

    cpdef Pattern score_pssm(
        self,
        PSSM pssm,
        list sequences,
    ):
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

        # store sequences in contiguous array to process them without the GIL
        seqptr = <SEQ_T**> alloca(length * sizeof(SEQ_T*))
        if seqptr is NULL:
            raise AllocationError("SEQ_T*", sizeof(SEQ_T*), length)
        for i, sequence in enumerate(sequences):
            seqptr[i] = sequence._seq
        # record options in the pattern object
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
                )
            # complete pattern
            libmeme.cisml.set_pattern_is_complete(pattern._pattern)
            # compute q-values
            values.num_items = libmeme.reservoir.get_reservoir_num_samples_retained(reservoir._reservoir)
            values.items = libmeme.reservoir.get_reservoir_samples(reservoir._reservoir)
            libmeme.cisml.pattern_calculate_qvalues(pattern._pattern, &values)


        return pattern

    cpdef Pattern score_motif(
        self,
        Motif motif,
        list sequences,
        Array bg_freqs,
        PriorDist prior_dist = None,
    ):
        # build PSSMs from motif
        pssm = motif.build_pssm(
            bg_freqs,
            bg_freqs,
            prior_dist,
            self.alpha,
            PSSM_RANGE,
            0,
            False
        )
        # score PSSM
        return self.score_pssm(pssm, sequences)