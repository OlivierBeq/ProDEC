# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-07-13

### Changed

- Modernized packaging: metadata now lives entirely in `pyproject.toml`
  (PEP 621); `setup.cfg` and `setup.py` were removed.
- Raised minimum supported Python version to 3.11; added support for 3.12 and 3.13.
- Replaced Travis CI (which was no longer running) with GitHub Actions.
- Replaced flake8 and its plugin stack with `ruff`; added `mypy` type checking.
- Added type hints throughout `src/prodec/` and a `py.typed` marker.
- Replaced the Sphinx/ReadTheDocs documentation with plain Markdown docs
  under `docs/`.
- Rewrote `README.md` with badges, a quickstart, and links to the full docs.

### Fixed

- Declared `tqdm` as a runtime dependency — it was imported by
  `prodec.utils` but never listed in package metadata.
- Fixed tests that referenced a since-renamed `Descriptor.Scales_names`/
  `Scales_values` API (now `Descriptor.definition`), and a test's hardcoded
  descriptor list that was missing two descriptors (`PhysChem`,
  `Zscale van Westen`) already present in `data.json`. These went unnoticed
  because CI (`travis.yml`, missing its leading dot) never actually ran.

### Added

- `CONTRIBUTING.md`, GitHub issue/PR templates, and this changelog.
- `pytest-cov` based coverage reporting.

## [1.0.2.post5] - 2023-01-16

### Fixed

- Prevent gap removal in PDT calculation.
- Fix ADFQ descriptor.

## [1.0.2.post4] - 2022-11-17

### Fixed

- Fix `Transform.is_compatible`.

## [1.0.2.post3] - 2022-10-28

- Metadata-only re-release of 1.0.2.post2.

## [1.0.2.post2] - 2022-10-28

### Fixed

- Apply FFT per dimension.
- Fix `Transform` instantiation.
- Raise exception in `pandas_get` on invalid input.

## [1.0.2.post1] - 2022-07-15

- Metadata-only re-release of 1.0.2.

## [1.0.2] - 2022-07-15

### Added

- Add Zscale Sandberg descriptor truncated to 3 principal components.

### Fixed

- Fix typos.

## [1.0.1.post3] - 2022-05-23

### Fixed

- Fix `Transform.available_transforms`.

## [1.0.1.post2] - 2022-04-03

- Metadata-only re-release of 1.0.1.

## [1.0.1] - 2022-04-03

Initial release published to PyPI.

### Added

- Core `ProteinDescriptors`, `Descriptor`, `Transform` and `TransformType`
  classes, with the AVG, ACC, PDT and FFT transforms.
- PhysChem descriptor.
- `pandas_get` methods returning results as a `pandas.DataFrame`, with
  support for user-supplied sequence identifiers.
- Packaging for distribution on PyPI.
