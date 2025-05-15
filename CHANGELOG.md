# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pymemesuite/compare/v0.1.0-a4...HEAD


## [v0.1.0-a4] - 2024-06-14
[v0.1.0-a4]: https://github.com/althonos/pymemesuite/compare/v0.1.0-a3...v0.1.0-a4

### Fixed
- Multithreading issue caused by Cython 3.0 `noexcept` code ([#2](https://github.com/althonos/pymemesuite/issues/2)).


## [v0.1.0-a3] - 2024-06-14
[v0.1.0-a3]: https://github.com/althonos/pymemesuite/compare/v0.1.0-a2...v0.1.0-a3

### Added
- `Matrix` constructor from an iterable of rows.
- Wheels for CPython 3.12 and PyPy 3.10.

### Changed
- Bumped Cython dependency to `v3.0`.


## [v0.1.0-a2] - 2023-04-04
[v0.1.0-a2]: https://github.com/althonos/pymemesuite/compare/v0.1.0-a1...v0.1.0-a2

### Added
- `pymemesuite.common.Background` to store the background model frequencies.

### Changed
- `FIMO.score_motif` and `Motif.build_pssm` now expect `Background` objects instead of raw `Array` objects with the background model frequencies. 


## [v0.1.0-a1] - 2022-05-19
[v0.1.0-a1]: https://github.com/althonos/pymemesuite/compare/02ff01e...v0.1.0-a1

Initial alpha release (test deployment to PyPI).
