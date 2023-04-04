import pickle
import unittest
import sys

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
