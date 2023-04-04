import pickle
import unittest
import sys
import tempfile
import textwrap

from pymemesuite.common import Alphabet, Background, Sequence


class TestBackground(unittest.TestCase):

    def test_from_sequences(self):
        dna = Alphabet.dna()
        seq1 = Sequence("AACGT")
        seq2 = Sequence("ACATA")
        bg = Background.from_sequences(dna, seq1, seq2)
        for freq, exp in zip( bg.frequencies, [0.3495, 0.1505, 0.1505, 0.3495]):
            self.assertAlmostEqual(freq, exp, places=2)

    def test_from_sequences_norc(self):
        dna = Alphabet.dna()
        seq1 = Sequence("AACGT")
        seq2 = Sequence("ACATA")
        bg = Background.from_sequences(dna, seq1, seq2, rc=False)
        for freq, exp in zip( bg.frequencies, [0.4975, 0.2005, 0.1015, 0.2005]):
            self.assertAlmostEqual(freq, exp, places=2)

    def test_from_sequences_empty(self):
        dna = Alphabet.dna()
        bg = Background.from_sequences(dna)
        for freq, exp in zip( bg.frequencies, [0.25, 0.25, 0.25, 0.25]):
            self.assertAlmostEqual(freq, exp, places=2)

    def test_from_file(self):
        with tempfile.NamedTemporaryFile(mode="w") as tmpfile:
            tmpfile.write(textwrap.dedent("""
                # 0-order Markov frequencies from file /tmp/x.fasta
                # seqs: 1    min: 4    max: 4    avg: 4.0    sum: 4    alph: DNA
                # order 0
                A 2.500e-01
                C 2.500e-01
                G 2.500e-01
                T 2.500e-01
            """.strip()))
            tmpfile.flush()
            bg = Background.from_file(Alphabet.dna(), tmpfile.name)
            for freq, exp in zip( bg.frequencies, [0.25, 0.25, 0.25, 0.25]):
                self.assertAlmostEqual(freq, exp, places=2)