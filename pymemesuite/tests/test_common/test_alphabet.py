import pickle
import unittest
import sys

from pymemesuite.common import Alphabet


class TestAlphabet(unittest.TestCase):

    def test_size(self):
        dna = Alphabet.dna()
        self.assertEqual(dna.size, 4)
        rna = Alphabet.rna()
        self.assertEqual(rna.size, 4)
        amino = Alphabet.amino()
        self.assertEqual(amino.size, 20)

    def test_wildcard(self):
        dna = Alphabet.dna()
        self.assertEqual(dna.wildcard, 'N')
