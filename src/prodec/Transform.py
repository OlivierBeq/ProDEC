# -*- coding: utf-8 -*-

"""Class handling transformation of descriptor raw values."""

import math
from typing import List, Optional, Union

import numpy as np
import pandas as pd

from .Descriptor import Descriptor
from .utils import _multiprocess_get

TRANSFORMS = {'AVG': {'Fullname': 'Domain average',
                      'Constant Length': True, 'Binary': False},
              'ACC': {'Fullname': 'Auto-cross covariances',
                      'Constant Length': True, 'Binary': False},
              'PDT': {'Fullname': 'Physicochemical Distance Tranformation',
                      'Constant Length': True, 'Binary': False}}


class Transform(Descriptor):
    """A class used to process raw values of protein descriptors."""

    def __init__(self, type: str, descriptor: Descriptor):
        """Instanciate a Transform.

        :param type: Type of Transform
        :param descriptor: Descriptor to be transformed
        """
        if type not in TRANSFORMS.keys():
            raise Exception(f'Transform type {type} can only be '
                            f'of the following values '
                            f'{list(TRANSFORMS.keys())}')
        self.Type = type
        self.Descriptor = descriptor
        self.Binary = TRANSFORMS[type]['Binary']

    @staticmethod
    def available_transforms():
        """Get possible transforms for all descriptors."""
        return list(TRANSFORMS.keys())

    @property
    def summary(self):
        """Get information on current transform."""
        return TRANSFORMS[self.Type]

    @staticmethod
    def is_compatible(type: str, descriptor: Descriptor) -> bool:
        """Say if types of tranform and descriptor are compatible.

        :param type: Type of transform
        :param descriptor: Descriptor to be transformed
        """
        if type not in TRANSFORMS.keys():
            raise NotImplementedError(f'Transform type {type}'
                                      ' is not implemented')
        return descriptor.Binary == (TRANSFORMS[type]['Binary'] or
                                     TRANSFORMS[type]['Binary'] == 'All')

    def get(self, sequence: str, flatten: bool = True, lag: int = 1,
            domains: int = 2, **kwargs):
        """Transform raw protein descriptor values.

        :param sequence : protein sequence
        :param flatten  : not to return dimensions separate, only for AVG and ACC
        :param lag      : lag value
        :param domains  : number of averaged domains
        """
        if not Transform.is_compatible(self.Type, self.Descriptor):
            raise Exception(f'Type of transform ({self.Type}) is incompatible'
                            f' with descriptor type ({self.Descriptor.Type})')
        if self.Type == 'AVG':
            return self.__average__(sequence, domains, flatten, **kwargs)
        elif self.Type == 'ACC':
            return self.__autocrosscovariance__(sequence, lag, flatten,
                                                **kwargs)
        elif self.Type == 'PDT':
            return self.__physicochem_distance_tansform__(sequence, lag,
                                                          **kwargs)
        else:
            raise NotImplementedError(f'Transform type {self.Type} '
                                      'is not implemented')

    def __average__(self, sequence: str, domains: int,
                    flatten: bool, **kwargs):
        """Calculate domain mean average values.

        :param sequence : protein sequence
        :param domains: number of domains to split the sequence into
        :param flatten: whether to give global average
                        or average per dimension.
        """
        length = len(sequence)
        if domains < 1 or domains > length:
            raise Exception(f'Number of domains ({domains}) '
                            'has to be greater or equal to 1 '
                            ' and lower or equal than the length '
                            f'of the sequence ({length})')
        # Keep track of original shape if not flat result
        if not flatten:
            raw = np.array(self.Descriptor.get(sequence, flatten=False,
                                               **kwargs))
            shape = raw.shape
        raw = self.Descriptor.get(sequence, flatten=True, **kwargs)
        average = np.zeros((domains * self.Descriptor.Size), dtype=np.float64)
        counts = np.zeros((domains * self.Descriptor.Size), dtype=np.int64)
        scaling_factor = ((len(raw) / self.Descriptor.Size) + 1) / domains
        for i in range(len(raw)):
            # Get scaled position
            seq_pos = int(i / self.Descriptor.Size)
            desc_pos = i % self.Descriptor.Size
            scaled_pos = (math.floor((seq_pos+1)/scaling_factor) *
                          self.Descriptor.Size + desc_pos)
            # Sum values and keep counts
            average[scaled_pos] += raw[i]
            counts[scaled_pos] += 1
        # Average values
        for i in range(len(average)):
            average[i] /= counts[i]
        # Reshape if not flatten
        if not flatten and len(shape) > 1:
            for row in range(self.Descriptor.Size):
                # Add mean average of each dimension
                weights = np.identity(self.Descriptor.Size)[row, :].tolist() \
                          * domains
                average = np.append(average, np.average(average[:domains *
                                    self.Descriptor.Size], weights=weights))
            average = np.reshape(np.round(average, 6), (domains+1, shape[1])).tolist()
            # Add global mean average
            average.append(np.round(np.average(average), 6))
        else:
            # Add global mean average
            average = np.append(np.round(average, 6),
                                np.round(np.average(average), 6)).tolist()
        return average

    def __autocrosscovariance__(self, sequence: str, lag: int,
                                flatten: bool, **kwargs):
        """Calculate auto-cross covariances values.

        :param sequence : protein sequence
        :param lag: Lag between amino acids
        :param flatten: Whether to give global average or average per dimension
        """
        length = (len(sequence) - 1
                  if self.Descriptor.Type == 'Distance'
                  else len(sequence))
        if lag < 1 or lag >= length:
            if self.Descriptor.Type != 'Distance':
                raise Exception(f'Lag ({lag}) has to be greater or equal to 1 '
                                'and lower than the length of the sequence '
                                f'({length})')
            else:
                raise Exception(f'Lag ({lag}) has to be greater or equal to 1 '
                                ' and lower than the length of the sequence '
                                f'minus one ({length-1})')
        acc = np.zeros((self.Descriptor.Size, self.Descriptor.Size),
                       dtype=np.float64)
        raw = self.Descriptor.get(sequence, flatten=True, **kwargs)
        for j in range(self.Descriptor.Size):
            for m in range(self.Descriptor.Size):
                for i in range(length - lag):
                    acc[m, j] += ((raw[i*self.Descriptor.Size+j] *
                                   raw[(i+lag)*self.Descriptor.Size+m]) /
                                  (length-lag))
        acc = np.round(acc, 6)
        if flatten:
            return acc.flatten(order='F').tolist()
        return acc.tolist()


    def __physicochem_distance_tansform__(self, sequence: str, lag: int, **kwargs):
        """Calculate physicochemical distance transform values.

        Gaps will automatically be omitted.
        :param sequence : protein sequence
        :param lag: lag between amino acids
        """
        length = (len(sequence) - 1
                  if self.Descriptor.Type == 'Distance'
                  else len(sequence))
        if lag >= length:
            if self.Descriptor.Type != 'Distance':
                raise Exception(f'Lag ({lag}) has to be greater or equal to 1'
                                ' and lower than the length of the sequence '
                                f'({length})')
            else:
                raise Exception(f'Lag ({lag}) has to be greater or equal to 1 '
                                'and lower than the length of the sequence '
                                f'minus one ({length-1})')
        # Normalisation of scores
        std_scores = np.array(self.Descriptor.get('ACDEFGHIKLMNPQRSTVWY',
                                                  flatten=False))
        mean_score = np.average(std_scores, axis=0)
        mean_difference = np.sqrt(np.sum(np.power(std_scores - mean_score, 2),
                                         axis=0))
        # Adjust raw values
        raw = np.array(self.Descriptor.get(sequence, flatten=False, **kwargs))
        raw = (raw - mean_score) / mean_difference
        # Removing gaps
        sequence = ''.join(filter(str.isalpha, sequence))
        # Get PDT
        pdt = np.zeros(self.Descriptor.Size)
        for i in range(length - lag):
            if self.Descriptor.Size == 1:
                pdt[0] += math.pow(raw[i] - raw[i+lag], 2)
            else:
                for m in range(self.Descriptor.Size):
                    pdt[m] += math.pow(raw[i, m] - raw[i+lag, m], 2)
        pdt /= length - lag
        pdt = np.round(pdt, 6)
        return pdt.tolist()


    def pandas_get(self, sequences: List[str],
                   ids: Optional[List[str]]=None,
                   lag: int = 1,
                   domains: int = 2,
                   average: bool = True,
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
        if domains > min(map(len, sequences)) or domains < 1:
            raise ValueError(f'Number of domains ({domains}) has to be greater or equal to 1'
                            ' and lower or equal than the length of the smallest sequence '
                            f'({min(map(len, sequences))})')
        if lag >= min(map(len, sequences)) or lag < 1:
            raise ValueError(f'Lag ({lag}) has to be greater or equal to 1 and '
                              'lower than the length of the smallest sequence '
                              f'({min(map(len, sequences))})')
        values = pd.DataFrame(_multiprocess_get(self, sequences=sequences, ids=ids, nproc=nproc, ipynb=ipynb, quiet=quiet,
                                                lag=lag, domains=domains, average=average))
        info = f'domains{domains}' if self.Type == "AVG" else f'lag{lag}'
        if ids:
            values.columns = ['ID'] + [f'{self.Type}_{info}_{self.Descriptor.ID.split()[0]}_{x}' for x in range(1, len(values.columns))]
        else:
            values.columns = [f'{self.Type}_{info}_{self.Descriptor.ID.split()[0]}_{x}' for x in range(1, len(values.columns) + 1)]
        return values
