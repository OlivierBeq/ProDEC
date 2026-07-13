# Usage guide

This guide walks through loading descriptors, calculating raw values for a
protein sequence, and applying the built-in transforms. For a quick class-by-class
reference, see [`api.md`](api.md).

## Table of contents

- [Getting started](#getting-started)
- [Advanced descriptor usage](#advanced-descriptor-usage)
  - [Flattening raw values](#flattening-raw-values)
  - [Dealing with gaps](#dealing-with-gaps)
  - [Non-standard amino acids](#non-standard-amino-acids)
  - [Raychaudhury's descriptor](#raychaudhurys-descriptor)
- [Advanced transform usage](#advanced-transform-usage)
  - [Compatibility](#compatibility)
  - [Transforms and advanced descriptor arguments](#transforms-and-advanced-descriptor-arguments)
- [Adding new descriptors](#adding-new-descriptors)
- [Checking descriptor support for amino acids](#checking-descriptor-support-for-amino-acids)

## Getting started

ProDEC is organised around four classes:

1. **`ProteinDescriptors`** — loads all available descriptors and lets you instantiate them.
2. **`Descriptor`** — instantiated from the former, retrieves raw descriptor values for a sequence.
3. **`Transform`** — computes domain averages (AVG), auto-cross covariances (ACC),
   physicochemical distance transformations (PDT), and the fast Fourier transform (FFT).
4. **`TransformType`** — identifies the transform to be performed.

Let's fetch a protein sequence from UniProt to work with:

```python
import urllib.request

url = 'https://rest.uniprot.org/uniprotkb/A0A5A9P0L4.fasta'
with urllib.request.urlopen(url) as data:
    sequence = ''.join(line.decode('ascii').strip() for line in data if not line.startswith(b'>'))
```

First, load the available descriptors:

```python
from prodec import *

pdescs = ProteinDescriptors()
```

and print out their IDs:

```python
print(pdescs.available_descriptors)
```

Identify the descriptor ID corresponding to Zscales (Hellberg *et al.* 1987):

```python
zscales = pdescs.get_descriptor('Zscale Hellberg')
```

Get information about the descriptor as defined in the original article:

```python
print(zscales.summary)
```

and the values defined for each amino acid:

```python
print(zscales.definition)
```

Now obtain the descriptor values for the protein sequence:

```python
raw_values = zscales.get(sequence)
```

To transform raw values, first check the available transforms (a static method):

```python
print(Transform.available_transforms())
```

Instantiate the desired transform (here, domain averages):

```python
avg_zscale = Transform(TransformType.AVG, zscales)
```

and obtain 50 domain averages (defaults to 2 if not specified):

```python
avg_values = avg_zscale.get(sequence, domains=50)
```

You can get information about the transform itself:

```python
print(avg_zscale.summary)
```

Similarly, ACC, PDT and FFT can be obtained with:

```python
acc_zscale = Transform(TransformType.ACC, zscales)
# or Transform('ACC', zscales)
acc_values = acc_zscale.get(sequence, lag=10)  # default lag=1

pdt_zscale = Transform(TransformType.PDT, zscales)
# or Transform('PDT', zscales)
pdt_values = pdt_zscale.get(sequence, lag=100)  # default lag=1

fft_zscale = Transform(TransformType.FFT, zscales)
# or Transform('FFT', zscales)
fft_values = fft_zscale.get(sequence)
```

## Advanced descriptor usage

### Flattening raw values

When multiple values are defined for one amino acid, the resulting sequence
descriptors are flattened by default: you get a list in which values for each
amino acid are contiguous. This can be turned off, resulting in a list of
lists, one per dimension (e.g. for Zscales Hellberg, a list of 3 sub-lists,
the first sub-list holding the values of the first dimension for the whole
sequence):

```python
zscales.get(sequence, flatten=False)
```

### Dealing with gaps

For aligned sequences, you may want to handle gaps explicitly. By default,
gaps are considered and given a value of `0.0`. Gaps can be omitted:

```python
zscales.get(sequence, gaps='omit')
```

or given an arbitrary value:

```python
zscales.get(sequence, gaps=-1)
```

### Non-standard amino acids

If working with a dictionary other than the 20 standard amino acids, you can
supply the ones you're working with — this only works if the descriptor you
use supports them:

```python
pdescs = ProteinDescriptors()
mydesc = pdescs.get_descriptor('Descriptor supporting Selenocysteine and Pyrrolysine')
mydesc.get(sequence, dictionary=list('ACDEFGHIKLMNOPQRSTUVWY'))
```

### Raychaudhury's descriptor

Raychaudhury *et al.*'s values can be weighted by different powers (default: `-4`):

```python
pdescs = ProteinDescriptors()
raych = pdescs.get_descriptor('Raychaudhury')
raych.get(sequence, power=-3)
```

Calculating Raychaudhury's values is **O(n²)**. A sliding-window optimization
brings this down to **O(n)**. By default, the window width is 120, giving
accuracy to the third decimal place. Change the width via `prec` (half the
window size):

```python
raych.get(sequence, prec=80)  # window size = 160
```

To turn the optimization off and get full precision:

```python
raych.get(sequence, prec=0)
```

## Advanced transform usage

### Compatibility

Some transforms cannot be calculated for binary descriptors, and vice versa.
Check compatibility between a transform and a descriptor:

```python
psm = pdescs.get_descriptor('PSM')
Transform.is_compatible('AVG', psm)
```

### Transforms and advanced descriptor arguments

Any argument accepted by `Descriptor.get` can also be passed to a transform's
`get` method:

```python
pdt_zscale.get(sequence, lag=10, average=False, flatten=False)

raych = pdescs.get_descriptor('Raychaudhury')
acc_raych = Transform('ACC', raych)
acc_raych.get(sequence, power=-3, gaps='omit', prec=100, flatten=False, lag=12)
```

## Adding new descriptors

Supplied descriptors are described in `data.json`, under `src/prodec`. The
list of available descriptors is loaded from that file when
`ProteinDescriptors` is instantiated. Add your own descriptor to the list,
respecting the existing format and giving it a unique ID, to make it
available.

## Checking descriptor support for amino acids

Check whether your engineered descriptor is compatible with a given sequence:

```python
vstv = pdescs.get_descriptor('VSTV')
vstv.is_sequence_valid('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
```
