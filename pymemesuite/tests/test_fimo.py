import csv
import unittest

try:
    import importlib.resources as importlib_resources
except ImportError:
    import importlib_resources

import pymemesuite
from pymemesuite.common import MotifFile, Sequence
from pymemesuite.fimo import FIMO

from . import fasta


class TestFIMO(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with importlib_resources.open_text("pymemesuite.tests.data.fimo", "mibig-genes.fna") as f:
            cls.sequences = [
                Sequence(
                    record.seq,
                    name=record.id.encode(),
                    description=record.description.encode()
                )
                for record in fasta.parse(f)
            ]

    def assertMatchEqual(self, match, line):
        self.assertEqual( match.source.name.decode(), line[2] )
        self.assertEqual( min(match.start, match.stop), int(line[3]) )
        self.assertEqual( max(match.start, match.stop), int(line[4]) )
        self.assertEqual( match.strand, line[5] )
        self.assertAlmostEqual( match.score, float(line[6]), places=2)
        self.assertAlmostEqual( match.pvalue, float(line[7]), places=2 )
        self.assertAlmostEqual( match.qvalue, float(line[8]), places=2 )
        self.assertEqual( match.sequence, line[9] )

    def test_both_strands(self):

        with importlib_resources.open_binary("pymemesuite.tests.data.fimo", "prodoric_mx000001_meme.txt") as f:
            motif_file = MotifFile(f, symmetrical=True)
            motif = motif_file.read()

        fimo = FIMO(both_strands=True)
        pattern = fimo.score_motif(motif, self.sequences, motif_file.background)

        with importlib_resources.open_text("pymemesuite.tests.data.fimo", "results.tsv") as f:
            reader = csv.reader(f.readlines()[1:-4], dialect="excel-tab")
            hits = [ line for line in reader ]
            hits.sort(key=lambda line: float(line[7]))

        self.assertEqual(len(pattern.matched_elements), len(hits))
        for match, line in zip(pattern.matched_elements, hits):
            self.assertMatchEqual(match, line)

    def test_single_strand(self):

        with importlib_resources.open_binary("pymemesuite.tests.data.fimo", "prodoric_mx000001_meme.txt") as f:
            motif_file = MotifFile(f, symmetrical=False)
            motif = motif_file.read()

        fimo = FIMO(both_strands=False)
        pattern = fimo.score_motif(motif, self.sequences, motif_file.background)

        with importlib_resources.open_text("pymemesuite.tests.data.fimo", "results-norc.tsv") as f:
            reader = csv.reader(f.readlines()[1:-4], dialect="excel-tab")
            hits = [ line for line in reader ]
            hits.sort(key=lambda line: float(line[7]))

        self.assertEqual(len(pattern.matched_elements), len(hits))
        for match, line in zip(pattern.matched_elements, hits):
            self.assertMatchEqual(match, line)
