# 🐍Ⓜ️ pyMEMEsuite [![Stars](https://img.shields.io/github/stars/althonos/pymemesuite.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pymemesuite/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to the [MEME suite](https://meme-suite.org), a collection of tools for the analysis of sequence motifs.*

<!-- [![Actions](https://img.shields.io/github/workflow/status/althonos/pymemesuite/Test/master?logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pymemesuite/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pymemesuite?logo=codecov&style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/pymemesuite/)
[![PyPI](https://img.shields.io/pypi/v/pymemesuite.svg?logo=pypi&style=flat-square&maxAge=3600)](https://pypi.org/project/pymemesuite)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pymemesuite?logo=anaconda&style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/pymemesuite)
[![AUR](https://img.shields.io/aur/version/python-pymemesuite?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pymemesuite)
[![Wheel](https://img.shields.io/pypi/wheel/pymemesuite.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pymemesuite/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pymemesuite.svg?logo=python&style=flat-square&maxAge=3600)](https://pypi.org/project/pymemesuite/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pymemesuite.svg?logo=python&style=flat-square&maxAge=3600&label=impl)](https://pypi.org/project/pymemesuite/#files)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pymemesuite/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/pymemesuite/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pymemesuite.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pymemesuite/issues)
[![Docs](https://img.shields.io/readthedocs/pymemesuite/latest?style=flat-square&maxAge=600)](https://pymemesuite.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pymemesuite/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpymemesuite)](https://pepy.tech/project/pymemesuite)
-->


## 🗺️ Overview

The MEME suite is a collection of tools used for discovery and analysis of
biological sequence motifs.

`pymemesuite` is a Python module, implemented using the [Cython](https://cython.org)
language, that provides bindings to the MEME suite. It directly interacts with the
MEME internals, which has the following advantages over CLI wrappers:

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
workflows involving `fimo`.*


<!-- ## 🔧 Installing

`pymemesuite` can be installed from [PyPI](https://pypi.org/project/pymemesuite/),
which hosts some pre-built CPython wheels for x86-64 Linux, as well as the
code required to compile from source with Cython:
```console
$ pip install pymemesuite
```

A [Bioconda](https://bioconda.github.io/) package is also available:
```console
$ conda install -c bioconda pymemesuite
``` -->

<!-- ## 📖 Documentation

A complete [API reference](https://pymemesuite.readthedocs.io/en/stable/api/) can
be found in the [online documentation](https://pymemesuite.readthedocs.io/), or
directly from the command line using
[`pydoc`](https://docs.python.org/3/library/pydoc.html):
```console
$ pydoc pymemesuite
``` -->


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/pymemesuite/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

<!-- ### 🏗️ Contributing

Contributions are more than welcome! See [`CONTRIBUTING.md`](https://github.com/althonos/pymemesuite/blob/master/CONTRIBUTING.md) for more details. -->


## ⚖️ License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/).
The MEME suite code is available under an academic license which allows
distribution and non-commercial usage. See `vendor/meme/COPYING` for more
information.

*This project is in no way affiliated, sponsored, or otherwise endorsed by
the [original MEME suite authors](https://meme-suite.org/meme/doc/authors.html).
It was developed by [Martin Larralde](https://github.com/althonos/pymemesuite)
during his PhD project at the [European Molecular Biology Laboratory](https://www.embl.de/)
in the [Zeller team](https://github.com/zellerlab).*
