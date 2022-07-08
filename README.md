# ProDEC

A package to easily calculate descriptors of protein sequences and their common transforms.

## Installation

    pip install prodec


## Getting started

ProDEC is organised in three classes:
 1. **ProteinDescripors** - loads all available descriptors and allows you to instantiate them
 2. **Descriptor** - instantiated from the latter, allows retrieval of raw descriptor values
 3. **Transform** - to calculate domain averages, auto-cross covariances (ACC), physicochemical distance transformations (PDT) and fast Fourier transform (FFT)
 4. **TransformType** - to identify the transform to be performed

Let us get the largest protein sequence from uniprot (*as of May 29th, 2020*).

    import urllib.request
    
    url = 'https://www.uniprot.org/uniprot/A0A5A9P0L4.fasta'
    with urllib.request.urlopen(url) as data:
        sequence = ''.join([line.decode('ascii').strip() for line in data][1:])
First load available descriptors:

    from  prodec import *
    pdescs = ProteinDescriptors()
and print out their ID:

    print(pdescs.available_descriptors)
Identify the descriptor ID corresponding to Zscales (Hellberg *et al.* 1987). 

    zscales = pdescs.get_descriptor('Zscale Hellberg')
Get information about the descriptor as defined in the original article

    print(zscales.summary)
and values defined for each amino acid.

    print(zscales.definition)

Now, obtain such descriptor values for the protein sequence.

    raw_values = zscales.get(sequence)

To transform raw values, first identify available transforms (static method).

    print(Transform.available_transforms())
Let us instantiate the desired transform (here domain averages)

    avg_zscale = Transform(TransformType.AVG, zscales)
and obtain 50 domain averages (defaults to 2 if not specified).

    avg_values = avg_zscale.get(sequence, domains=50)

One can get information about the transform.

    print(avg_zscale.summary)

Similarly, ACC, PDT and FFT can be obtain with

    acc_zscale = Transform(TransformType.ACC, zscales)
    # or Transform('ACC', zscales)
    acc_values = acc_zscale.get(sequence, lag=10) # default lag=1
    pdt_zscale = Transform(TransformType.PDT, zscales)
    # or Transform('PDT', zscales)
    pdt_values = pdt_zscale.get(sequence, lag=100) # default lag=1
    fft_zscale = Transform(TransformType.FFT, zscales)
    # or Transform('FFT', zscales)
    fft_values = pdt_zscale.get(sequence)

## Advanced usage
### Descriptors

 - ***Flattening raw values***

In the case of multiple values being defined for one amino acid, the resulting sequence descriptors are flattened by default. This means that one gets a list in which values for each amino acid are contiguous.
This feature can be turned off, resulting in a list of lists, each dimension being separate from the other (e.g. for Zscales Hellberg, a list containing 3 sub-lists: the first sub-list with values of the first dimension for the whole sequence). 

    zscales.get(sequence, flatten=False)

 - ***Dealing with gaps***

In the case of aligned sequences, one may want to omit gaps. By default, gaps are considered and given a value of  0.0 . Gaps can either be omitted like so:

    zscales.get(sequence, gaps='omit')
or given any arbitrary value

    zscales.get(sequence, gaps=-1)

 - ***Non-standard amino acids***

If working with another dictionary than the 20 standard amino acids, one can provide the ones they are working with. This is only possible if the user defines their own descriptor supporting these aminoacids.

    pdescs = prodec.ProteinDescriptors()
    mydesc = pdescs.get('Descriptor supporting Selenocysteine and Pyrrolysine')
    mydesc.get(sequence, dictionary=list('ACDEFGHIKLMNOPQRSTUVWY'))
 

 - ***Raychaudhury's descriptor***

Rachaudhury  *et al.*'s values can be weighted by different powers (*default: -4*).

    pdescs = prodec.ProteinDescriptors()
    raych = pdescs.get('Raychaudhury')
    raych.get(sequence, power=-3)

Calculation of Raychaudhury's values is ***O(n²)*** . To speed this calculation, a sliding window optimization has been made, resulting in an  ***O(n)*** algorithm. By default the window width is set to 120 giving accuracy to the third decimal place. One may change the width by specifying the precision (half of the window size).

    raych.get(sequence, prec=80) # Window size = 160
To turn the optimization off and get full precision:

    raych.get(sequence, prec=0)

### Transfoms

 - ***Compatibility***

Some transforms cannot be calculated for binary descriptors. Some others can only be calculated with binary descriptors. One can check for compatibility between a transform and a descriptor.

    psm = pdescs.get_descriptor('PSM')
    prodec.Transform.is_compatible('AVG', 'PSM')

 - ***Transforms and advanced descriptor arguments***

All arguments a *Descriptor* accepts can be supplied to a transform's *get* method.

    pdt_zscale.get(sequence, lag=10, average=False, flatten=False)
    raych = pdescs.get('Raychaudhury')
    acc_raych = prodec.Transform('ACC', raych)
    acc_raych.get(sequence, power=-3, gaps='omit', prec=100, flatten=False, lag=12)

### Adding new descriptors
Supplied descriptors are described in the file named *data.json* under the *src* folder.
The list of available descriptors is loaded from the *data.json* file when **ProteinDescriptors** is instantiated.
Add your favorite descriptor to the list, respecting the format of the file and giving it a unique ID, for it to be available.

### Checking descriptor for amino acids support
One can check the compatibility of their engineered descriptor with any sequence.

    vstv= pdescs.get_descriptor('VSTV')
    vstv.is_sequence_valid('ABCDEFGHIJKLMNOPQRSTUVWXYZ')

