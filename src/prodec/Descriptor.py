# -*- coding: utf-8 -*-

"""Class handling amino acid description in a sequence."""

import itertools
import multiprocessing
import warnings
from typing import Callable, List, Optional, Type, Union


import numpy as np
import pandas as pd

from .utils import alpha_carbon_distance_edge_vector_memory, \
                   alpha_carbon_distance_edge_vector_speed, \
                   enough_avail_memory, get_values, \
                   raychaudhury_memory, raychaudhury_speed, \
                   _multiprocess_get


class Descriptor:
    """A class used for calculating protein descriptors."""

    def __init__(self, desc_data: dict):
        """Instantiate a Descriptor.

        :param desc_data: Descriptor data
        """
        self.ID = desc_data['ID']
        self.Type = desc_data['Type']
        self.Name = desc_data['Name']
        self._scales_names = desc_data['Scales']['Names']
        self._scales_values = desc_data['Scales']['Values']
        self.Info = {'Authors': desc_data['Authors'],
                     'Year': desc_data['Year'],
                     'Journal': desc_data['Journal'],
                     'DOI': desc_data['DOI'],
                     'PMID': desc_data['PMID'],
                     'Patent': desc_data['Patent']}
        first_val = list(self._scales_values.values())[0]
        if isinstance(first_val, list):
            self.Size = len(first_val)
        elif isinstance(first_val, dict):
            self.Size = 1
        else:
            self.Size = 1
        self.first_val_type = type(first_val)
        self.Binary = set(get_values(self._scales_values)) == {0, 1}

    @property
    def summary(self):
        """Get summary of the descriptor."""
        return self.Info

    @property
    def definition(self):
        """Get values defining the descriptor."""
        return self._scales_names, self._scales_values

    def is_sequence_valid(self, sequence: str):
        """Check if a sequence can be fully described using the current descriptor.

        :param sequence: Protein sequence
        """
        for letter in sequence:
            if letter not in self._scales_values.keys():
                return False
        return True

    def get(self, sequence: str, flatten: bool=True, gaps: Union[int, str]=0,
            prec: Optional[int]=60, power: int=-4, dtype: Type=np.float16, fast: bool=False,
            **kwargs):
        """Get the raw values of the provided sequence.

        :param sequence : protein sequence
        :param flatten  : not to return dimensions separate
        :param gaps     : how should gaps be considered.
                          Allowed values: 'omit' or 0, ...+inf
        :param dtype    : data type for memory efficiency
        :param prec     : max number of amino acids to cosider 
                          before and after the current Calpha
        :param power    : power the topological distance is raised to
        :param fast     : whether to speed up at the cost of intense memory use
        """
        # Dealing with gaps
        if gaps == 'omit':
            sequence = ''.join(filter(str.isalpha, sequence))
            replacer = 0
        else:
            sequence = ''.join([x if str.isalpha(x)
                                else '-' for x in sequence])
            replacer = gaps
            if self.Type == 'Linear':
                if self.first_val_type in [float, int]:
                    self._scales_values['-'] = replacer
                else:
                    self._scales_values['-'] = [replacer] * self.Size
            elif self.Type == 'Distance':
                self._scales_values['-'] = {}
                for key, value in self._scales_values.items():
                    self._scales_values[key]['-'] = value
                    self._scales_values['-'][key] = value
            else:
                self._scales_values['-'] = replacer
        # Checking sequence
        if not self.is_sequence_valid(sequence):
            raise Exception('Sequence has unsupported amino acid')
        # Create final array
        if self.Size > 1:
            dtype_ = int if self.Binary else float
            values = np.zeros((len(sequence), self.Size), dtype=dtype_)
            values.fill(replacer)
        elif self.Type == 'Distance':
            values = [replacer] * (len(sequence) - 1)
        else:
            values = [replacer] * len(sequence)
        if self.Type == 'Linear':
            for i in range(len(sequence)):
                if self.Size > 1:
                    for j in range(self.Size):
                        values[i, j] = self._scales_values[sequence[i]][j]
                else:
                    values[i] = self._scales_values[sequence[i]]
        elif self.Type == 'Distance':
            for i in range(len(sequence) - 1):
                values[i] = self._scales_values[sequence[i]][sequence[i+1]]
        else:
            if self.ID == 'Raychaudhury':
                def mapping(x):
                    return self._scales_values[x]
                enough_ram = enough_avail_memory(len(sequence), dtype)
                if not prec:
                    if fast and enough_ram:
                        values = raychaudhury_speed(sequence, mapping,
                                                    power, dtype)
                    else:
                        if not enough_ram:
                            warnings.warn('Not enough memory available, reverting to slow calculation.')
                        values = raychaudhury_memory(sequence, mapping, power, dtype)
                else:
                    if fast and enough_ram:
                        values = alpha_carbon_distance_edge_vector_speed(
                                                                 sequence,
                                                                 mapping,
                                                                 1, power,
                                                                 dtype,
                                                                 prec, 1)
                    else:
                        if not enough_ram:
                            warnings.warn('Not enough memory available, reverting to slow calculation.')
                        values = alpha_carbon_distance_edge_vector_memory(sequence,
                                                                    mapping,
                                                                    1, power,
                                                                    dtype,
                                                                    prec, 1)
        if isinstance(values, list):
            return values
        if not flatten:
            return values.tolist()
        else:
            return values.flatten(order='C').tolist()

    def pandas_get(self, sequences: List[str],
                   ids: Optional[List[str]]=None,
                   gaps: Union[int, str]=0,
                   prec: Optional[int]=60,
                   power: int=-4,
                   dtype: Type=np.float16,
                   fast: bool=False,
                   nproc: Optional[int]=None,
                   quiet=False,
                   ipynb=False, **kwargs) -> pd.DataFrame:
        """Get the raw values of the provided sequences in a pandas DataFrame.

        :param sequences: protein sequences
        :param gaps     : how should gaps be considered.
                          Allowed values: 'omit' or 0, ...+inf
        :param prec     : max number of amino acids to cosider 
                          before and after the current Calpha
        :param power    : power the topological distance is raised to
        :param dtype    : data type for memory efficiency
        :param fast     : whether to speed up at the cost of intense memory use
        :param nproc    : number of concurrent processes to run
        :param quiet    : whether to report progress
        :param ipynb    : whether the function is run from a notebook
        """
        values = pd.DataFrame(_multiprocess_get(self, sequences=sequences, ids=ids, nproc=nproc, ipynb=ipynb, quiet=quiet,
                                                gaps=gaps, prec=prec, power=power,
                                                fast=fast, dtype=dtype))
        if ids:
            values.columns = ['ID'] + [f'{self.ID.split()[0]}_{x}' for x in range(1, len(values.columns))]
        else:
            values.columns = [f'{self.ID.split()[0]}_{x}' for x in range(1, len(values.columns) + 1)]
        return values

                