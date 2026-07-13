# API reference

Hand-written reference for ProDEC's public classes. For a task-oriented
walkthrough, see [`usage.md`](usage.md).

## `ProteinDescriptors`

Loads and caches the catalog of available protein descriptors.

```python
ProteinDescriptors(json: Optional[str] = None)
```

- `json` — optional path to a custom descriptor JSON file. Defaults to the
  bundled `data.json`.

| Member | Description |
| --- | --- |
| `available_descriptors` (property) | Sorted list of loaded descriptor IDs (or names, if a descriptor has no ID). |
| `update_available_descriptors(json)` | Reload descriptors from a path or file handle, replacing the currently loaded set. |
| `get_descriptor(id: str) -> Descriptor` | Instantiate (and cache) a `Descriptor` for the given ID. |

## `Descriptor`

Represents a single descriptor and computes raw per-residue values for a sequence.

```python
Descriptor(desc_data: dict)
```

Normally instantiated via `ProteinDescriptors.get_descriptor`, not directly.

| Member | Description |
| --- | --- |
| `summary` (property) | Dict of descriptor metadata: authors, year, journal, DOI, PMID, patent. |
| `definition` (property) | Tuple of `(scale names, scale values)` as defined for the descriptor. |
| `is_sequence_valid(sequence: str) -> bool` | Whether every residue in `sequence` is supported by this descriptor. |
| `get(sequence, flatten=True, gaps=0, prec=60, power=-4, dtype=np.float16, fast=False, **kwargs)` | Raw descriptor values for `sequence`. See [`usage.md`](usage.md#advanced-descriptor-usage) for the meaning of each argument. |
| `pandas_get(sequences, ids=None, ..., nproc=None, quiet=False, ipynb=False)` | Same as `get`, computed over many sequences in parallel and returned as a `pandas.DataFrame`. |

## `TransformType`

`Enum` identifying the transform to perform: `AVG`, `ACC`, `PDT`, `FFT`.

| Member | Description |
| --- | --- |
| `fullname` (property) | Human-readable name of the transform. |
| `constant_length` (property) | Whether the transform produces output of constant length regardless of sequence length. |
| `binary` (property) | Whether the transform is restricted to binary descriptors. |
| `available()` (classmethod) | List of transform names and members that can be used to instantiate a `Transform`. |

## `Transform`

Applies a `TransformType` on top of a `Descriptor`'s raw values.

```python
Transform(type: Union[str, TransformType], descriptor: Descriptor)
```

| Member | Description |
| --- | --- |
| `summary` (property) | Metadata dict for the current transform type. |
| `available_transforms()` (staticmethod) | Same as `TransformType.available()`. |
| `is_compatible(transform, descriptor) -> bool` (staticmethod) | Whether a transform type and a descriptor can be combined. |
| `get(sequence, flatten=True, lag=1, domains=2, normalize=True, **kwargs)` | Transformed values for `sequence`. `kwargs` are forwarded to `Descriptor.get`. See [`usage.md`](usage.md#advanced-transform-usage). |
| `pandas_get(sequences, ids=None, lag=1, domains=2, nproc=None, quiet=False, ipynb=False)` | Same as `get`, computed over many sequences in parallel and returned as a `pandas.DataFrame`. |
