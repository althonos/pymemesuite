# 🐍Ⓜ️ PyMEMEsuite [![Stars](https://img.shields.io/github/stars/althonos/pymemesuite.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pymemesuite/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to the [MEME suite](https://meme-suite.org), a collection of tools for the analysis of sequence motifs.*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pymemesuite/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pymemesuite/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pymemesuite?logo=codecov&style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/pymemesuite/)
[![PyPI](https://img.shields.io/pypi/v/pymemesuite.svg?logo=pypi&style=flat-square&maxAge=3600)](https://pypi.org/project/pymemesuite)
[![Wheel](https://img.shields.io/pypi/wheel/pymemesuite.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pymemesuite/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pymemesuite.svg?logo=python&style=flat-square&maxAge=3600)](https://pypi.org/project/pymemesuite/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pymemesuite.svg?logo=python&style=flat-square&maxAge=3600&label=impl)](https://pypi.org/project/pymemesuite/#files)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pymemesuite/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pymemesuite.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pymemesuite/issues)
[![Docs](https://img.shields.io/readthedocs/pymemesuite/latest?style=flat-square&maxAge=600)](https://pymemesuite.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pymemesuite/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpymemesuite)](https://pepy.tech/project/pymemesuite)


**🚩 If you are interested in running sequence motif analyses with PSSMs, I recommend using the [`lightmotif`](https://github.com/althonos/lightmotif/tree/main/lightmotif-py) project, which is a re-implementation of the algorithm used in FIMO, which can process sequences at speeds of GB/s.**

## 🗺️ Overview

The [MEME suite](https://meme-suite.org/) is a collection of tools used for
discovery and analysis of biological sequence motifs.

`pymemesuite` is a Python module, implemented using the [Cython](https://cython.org)
language, that provides bindings to the [MEME](https://meme-suite.org/) suite.
It directly interacts with the MEME internals, which has the following
advantages over CLI wrappers:

- **single dependency**: If your software or your analysis pipeline is
  distributed as a Python package, you can add `pymemesuite` as a dependency
  to your project, and stop worrying about the MEME binaries being properly
  setup on the end-user machine.
- **no intermediate files**: Everything happens in memory, in Python objects
  you have control on, making it easier to pass your inputs to MEME without
  needing to write them to a temporary file. Output retrieval is also done
  in memory.

*This library is still a work-in-progress, and in an experimental stage,
but it should already pack enough features to run biological analyses or
workflows involving [FIMO](https://meme-suite.org/meme/doc/fimo.html).*


## 🔧 Installing

`pymemesuite` can be installed from [PyPI](https://pypi.org/project/pymemesuite/),
which hosts some pre-built CPython wheels for x86-64 Linux, as well as the
code required to compile from source with Cython:
```console
$ pip install pymemesuite
```

<!-- ## 📖 Documentation

A complete [API reference](https://pymemesuite.readthedocs.io/en/stable/api/) can
be found in the [online documentation](https://pymemesuite.readthedocs.io/), or
directly from the command line using
[`pydoc`](https://docs.python.org/3/library/pydoc.html):
```console
$ pydoc pymemesuite
``` -->

## 💡 Example

Use `MotifFile` to load a motif from a MEME motif file, and display the
consensus motif sequence followed by the letter frequencies:

```python
from pymemesuite.common import MotifFile

with MotifFile("tests/data/fimo/prodoric_mx000001_meme.txt") as motif_file:
    motif = motif_file.read()

print(motif.name.decode())
print(motif.consensus)

for row in motif.frequencies:
    print(" ".join(f'{freq:.2f}' for freq in row))
```

Then use `FIMO` to find occurences of this particular motif in a collection of
sequences, and show coordinates of matches:

```python
import Bio.SeqIO
from pymemesuite.common import Sequence
from pymemesuite.fimo import FIMO

sequences = [
    Sequence(str(record.seq), name=record.id.encode())
    for record in Bio.SeqIO.parse("tests/data/fimo/mibig-genes.fna", "fasta")
]

fimo = FIMO(both_strands=False)
pattern = fimo.score_motif(motif, sequences, motif_file.background)

for m in pattern.matched_elements:
    print(
        m.source.accession.decode(),
        m.start,
        m.stop,
        m.strand,
        m.score,
        m.pvalue,
        m.qvalue
    )
```

You should then get a single matched element on the forward strand:
```
BGC0002035.1_3425_15590 6700 6714 + 9.328571428571422 1.1024163606971822e-05 0.6174858127445146
```

## 📋 Features

### 🧶 Thread-safety

`FIMO` objects are thread-safe, and the `FIMO.score_motif` and `FIMO.score_pssm`
methods are re-entrant. This means you can search occurences of several
motifs in parallel with a `ThreadPool` and a single `FIMO` instance:
```python
from multiprocessing.pool import ThreadPool
from pymemesuite.fimo import FIMO

fimo = FIMO()
with ThreadPool() as pool:
    patterns = pool.map(
        lambda motif: fimo.score_motif(motif, sequences, background),
        motifs
    )
```

### 📌 Roadmap

- [ ] **error management**: Make sure to catch exceptions raised by the MEME core without exiting forcefully.
- [ ] **transfac**: Support for TRANSFAC motifs in addition to MEME motifs.
- [ ] **meme**: Motif discovery through enrichment analysis between two collections of sequences.


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/pymemesuite/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See [`CONTRIBUTING.md`](https://github.com/althonos/pymemesuite/blob/master/CONTRIBUTING.md) for more details.


## ⚖️ License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/).
The MEME suite code is available under an academic license which allows
distribution and non-commercial usage. See `vendor/meme/COPYING` for more
information.

Test sequence data were obtained from [MIBiG](https://mibig.secondarymetabolites.org/)
and are distributed under the [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
license. Test motifs were obtained from [PRODORIC](https://www.prodoric.de) and are
distributed under the [CC BY-NC 4.0](https://creativecommons.org/licenses/by/4.0/)
license.

*This project is in no way affiliated, sponsored, or otherwise endorsed by
the [original MEME suite authors](https://meme-suite.org/meme/doc/authors.html).
It was developed by [Martin Larralde](https://github.com/althonos/pymemesuite)
during his PhD project at the [European Molecular Biology Laboratory](https://www.embl.de/)
in the [Zeller team](https://github.com/zellerlab).*
