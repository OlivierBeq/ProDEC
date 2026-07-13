# Contributing to ProDEC

Thanks for considering a contribution! This project is small and welcomes
bug reports, descriptor additions, and code improvements.

## Development setup

```bash
git clone https://github.com/OlivierBeq/ProDEC.git
cd ProDEC
python -m venv .venv
source .venv/bin/activate  # .venv\Scripts\activate on Windows
pip install -e ".[testing,lint]"
```

## Running the checks locally

```bash
# Lint
ruff check .

# Type-check
mypy src/prodec

# Tests, with coverage
pytest --cov=prodec tests/
```

All three run in CI (`.github/workflows/ci.yml`) against Python 3.11, 3.12
and 3.13 on Linux and Windows. Please make sure they pass before opening a
pull request.

## Adding a descriptor

Descriptors are described declaratively in `src/prodec/data.json`. See
[`docs/usage.md`](docs/usage.md#adding-new-descriptors) for the expected
format, and add a test in `tests/test_descriptor.py` covering the new entry.

## Pull requests

- Keep changes focused; unrelated cleanup belongs in a separate PR.
- Add or update tests for any behavioural change.
- Update `CHANGELOG.md` under `[Unreleased]`.
