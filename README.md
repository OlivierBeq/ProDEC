# ProDEC

[![PyPI version](https://img.shields.io/pypi/v/prodec.svg)](https://pypi.org/project/prodec/)
[![Python versions](https://img.shields.io/pypi/pyversions/prodec.svg)](https://pypi.org/project/prodec/)
[![CI](https://github.com/OlivierBeq/ProDEC/actions/workflows/ci.yml/badge.svg)](https://github.com/OlivierBeq/ProDEC/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A package to easily calculate descriptors of protein sequences and their common transforms
(domain averages, auto-cross covariances, physicochemical distance transformations, and the
fast Fourier transform).

## Table of contents

- [Installation](#installation)
- [Quickstart](#quickstart)
- [Documentation](#documentation)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

## 📦 Installation

```bash
pip install prodec
```

## Quickstart

```python
from prodec import ProteinDescriptors, Transform, TransformType

sequence = 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTL'

# Load the catalog of available descriptors and pick one
pdescs = ProteinDescriptors()
zscales = pdescs.get_descriptor('Zscale Hellberg')

# Get raw per-residue descriptor values
raw_values = zscales.get(sequence)

# Apply a transform, e.g. domain averages over 5 domains
avg_zscale = Transform(TransformType.AVG, zscales)
avg_values = avg_zscale.get(sequence, domains=5)
```

## Documentation

- [`docs/usage.md`](docs/usage.md) — full walkthrough, gap handling, non-standard
  amino acids, transform compatibility, and how to add new descriptors.
- [`docs/api.md`](docs/api.md) — reference for `ProteinDescriptors`, `Descriptor`,
  `Transform`, and `TransformType`.

## 📖 Citation

If you use ProDEC in your work, please cite this repository:

> Béquignon, O. J. M. *ProDEC: a package to calculate protein sequence descriptors and their
> common transforms.* https://github.com/OlivierBeq/ProDEC

The bundled descriptors and transforms themselves originate from the literature (e.g. Zscales
by Hellberg *et al.*, and the Raychaudhury descriptor) — see `Descriptor.summary` for the
citation of a specific descriptor, and refer to the corresponding original publication when
using it in your own research.

## Contributing

Contributions are welcome — see [`CONTRIBUTING.md`](CONTRIBUTING.md) for how to set up a
development environment and run the test/lint/type-check suite locally.

## License

ProDEC is distributed under the [MIT License](LICENSE).
