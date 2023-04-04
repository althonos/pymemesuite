import os
import io
import pickle
import shutil
import sys
import unittest
import tempfile
import warnings

try:
    import importlib.resources as importlib_resources
except ImportError:
    import importlib_resources

from pymemesuite.common import MotifFile


class TestMotifFile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.motif_data = importlib_resources.read_binary("pymemesuite.tests.data.fimo", "prodoric_mx000001_meme.txt")

    def test_fileobject(self):
        fileobj = io.BytesIO(self.motif_data)
        motif_file = MotifFile(fileobj)
        motif = motif_file.read()
        self.assertEqual(motif.accession, b"MX000001")
        motif = motif_file.read()
        self.assertIs(motif, None)

    def test_filename(self):
        try:
            fd, filename = tempfile.mkstemp(suffix=".txt")
            with os.fdopen(fd, "wb") as f:
                f.write(self.motif_data)
            with MotifFile(filename) as motif_file:
                motif = motif_file.read()
                self.assertEqual(motif.accession, b"MX000001")
                motif = motif_file.read()
                self.assertIs(motif, None)
        finally:
            os.remove(filename)

    def test_background(self):
        # try reading background before reading a motif
        fileobj = io.BytesIO(self.motif_data)
        motif_file = MotifFile(fileobj)
        self.assertIsNot(motif_file.background, None)
        # try reading background after reading a motif
        fileobj = io.BytesIO(self.motif_data)
        motif_file = MotifFile(fileobj)
        motif_file.read()
        self.assertIsNot(motif_file.background, None)
