# -*- coding: utf-8 -*-

"""Class handling transformation of descriptor raw values."""

import math
from enum import Enum, unique
from itertools import chain
from typing import List, Optional, Union

import numpy as np
import pandas as pd

from .Descriptor import Descriptor
from .utils import _multiprocess_get


@unique
class TransformType(Enum):
    """Types of transforms one can apply to descriptors."""
    AVG = 0
    ACC = 1
    PDT = 2
    FFT = 3

    def __repr__(self):
        """Create the transform type representation."""
        return f'<{self.__class__.__name__}.{self.name}>'

    @property
    def data(self):
        transforms = [{'Fullname': 'Domain average', 'Constant Length': True, 'Binary': False},
                      {'Fullname': 'Auto-cross covariances', 'Constant Length': True, 'Binary': False},
                      {'Fullname': 'Physicochemical Distance Tranformation', 'Constant Length': True, 'Binary': False},
                      {'Fullname': 'Fast Fourier Transform', 'Constant Length': False, 'Binary': False}
                      ]
        return transforms[self.value]

    @property
    def fullname(self):
        return self.data['Fullname']

    @property
    def constant_length(self):
        return self.data['Constant Length']

    @property
    def binary(self):
        return self.data['Binary']

    @classmethod
    def available(cls):
        return list(cls.__members__.keys()) + list(cls.__members__.values())


class Transform:
    """A class used to process raw values of protein descriptors."""

    def __init__(self, type: Union[str, TransformType], descriptor: Descriptor):
        """Instantiate a Transform.

        :param type: Type of Transform
        :param descriptor: Descriptor to be transformed
        """
        if not isinstance(descriptor, Descriptor):
            raise ValueError('descriptor must br a valid ProDEC descriptor instance')
        if not isinstance(type, (str, TransformType)):
            raise ValueError('transform type must be either a str or a TransformType')
        if isinstance(type, str):
            if type not in TransformType.available():
                raise ValueError('transform type is not supported')
            self.Type = TransformType[type]
        else:
            self.Type = type
        self.Descriptor = descriptor
        self.Binary = self.Type.binary

    @staticmethod
    def available_transforms():
        """Get possible transforms for all descriptors."""
        return [x for x in TransformType.available()]

    @property
    def summary(self):
        """Get information on current transform."""
        return self.Type.data

    @staticmethod
    def is_compatible(transform: Union[str, TransformType], descriptor: Descriptor) -> bool:
        """Say if types of transforms and descriptors are compatible.

        :param transform: Type of transform
        :param descriptor: Descriptor to be transformed
        """
        if not isinstance(transform, (str, TransformType)):
            raise ValueError('transform must be either a str or a TransformType')
        if transform not in TransformType.available():
            raise ValueError(f'Transform type {transform} is not implemented')
        if isinstance(transform, str):
            return descriptor.Binary == TransformType[transform].data['Binary'] or TransformType[transform].data[
                'Binary'] == 'All'
        return descriptor.Binary == transform.data['Binary'] or transform.data['Binary'] == 'All'

    def get(self, sequence: str, flatten: bool = True, lag: int = 1,
            domains: int = 2, normalize: bool = True, **kwargs):
        """Transform raw protein descriptor values.

        :param sequence : Protein sequence
        :param flatten  : Should dimensions not be returned separately (only for AVG and ACC)
        :param lag      : Lag between amino acids (only for ACC and PDT)
        :param domains  : Number of domains to split the sequence into (only for AVG)
        :param normalize: Should FFT amplitudes be normalized
        :param kwargs: keyword arguments passed to the get method of the underlying descriptor
        """
        if not Transform.is_compatible(self.Type, self.Descriptor):
            raise Exception(f'Type of transform ({self.Type}) is incompatible'
                            f' with descriptor type ({self.Descriptor.Type})')
        if self.Type is TransformType.AVG:
            return self.__average__(sequence, domains, flatten, **kwargs)
        elif self.Type is TransformType.ACC:
            return self.__autocrosscovariance__(sequence, lag, flatten,
                                                **kwargs)
        elif self.Type is TransformType.PDT:
            return self.__physicochem_distance_transform__(sequence, lag,
                                                           **kwargs)
        elif self.Type is TransformType.FFT:
            try:
                return self.__fast_fourier_transform__(sequence, normalize, flatten, **kwargs)
            except Exception as e:
                raise RuntimeException('something occured') from e
        else:
            raise NotImplementedError(f'Transform type {self.Type} is not implemented')

    def __average__(self, sequence: str, domains: int,
                    flatten: bool, **kwargs):
        """Calculate domain mean average values.

        :param sequence : Protein sequence
        :param domains: Number of domains to split the sequence into
        :param flatten: Should dimensions not be returned separately
        :param kwargs: keyword arguments passed to the get method of the underlying descriptor
        """
        length = (len(sequence) - 1
                  if self.Descriptor.Type == 'Distance'
                  else len(sequence))
        if domains < 1 or domains > length:
            if self.Descriptor.Type != 'Distance':
                raise Exception(f'Number of domains ({domains}) '
                                'has to be greater or equal to 1 '
                                ' and lower or equal than the length '
                                f'of the sequence ({length})')
            else:
                raise Exception(f'Number of domains ({domains}) '
                                'has to be greater or equal to 1 '
                                ' and lower or equal than the length '
                                f'of the sequence minus one ({length})')
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

        :param sequence : Protein sequence
        :param lag: Lag between amino acids
        :param flatten: Should dimensions not be returned separately
        :param kwargs: keyword arguments passed to the get method of the underlying descriptor
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
            return acc.flatten(order='C').tolist()
        return acc.tolist()

    def __physicochem_distance_transform__(self, sequence: str, lag: int, **kwargs):
        """Calculate physicochemical distance transform values.

        Gaps will automatically be omitted.
        :param sequence : Protein sequence
        :param lag: Lag between amino acids
        :param kwargs: keyword arguments passed to the get method of the underlying descriptor
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

    def __fast_fourier_transform__(self, sequence: str, normalize: bool = True, flatten: bool = True, **kwargs):
        """
        Calculate Fast Fourier transform values.

        Code is adapted from Mehdi D. Davari's PyPEF (https://github.com/mdavari/PyPEF)
        PyPEF—An Integrated Framework for Data-Driven Protein Engineering
        Niklas E. Siedhoff, Alexander-Maurice Illig, Ulrich Schwaneberg, and Mehdi D. Davari
        Journal of Chemical Information and Modeling 2021 61 (7), 3463-3476
        DOI: 10.1021/acs.jcim.1c00099

        :param sequence: sequence to get transformed values of
        :param kwargs: keyword arguments passed to the get method of descriptor
        :param normalize: whether amplitudes should be normalized
        :retun: amplitudes of the Fast Fourier Transform of descriptor values of the sequence
        """
        # Get exponent k
        k = np.log2(len(sequence))
        # Get raw scores
        raw = np.array(self.Descriptor.get(sequence, flatten=False, **kwargs))
        # Obtain tranformed values
        result = []
        for dim in range(raw.shape[1]):
            # Calculate mean to substract from raw values
            # to avoid artefacts of FFT
            mean = np.mean(raw[:, dim], axis=0)
            centered = np.subtract(raw[:, dim], mean)
            # Check if len raw = 2^k with k an integer
            if math.isclose(int(k), k, rel_tol=1e-8):
                # Pad with 0s
                reshaped = np.zeros(pow(2, (int(k) + 1)))  # reshape array
                for index, value in enumerate(centered):
                    reshaped[index] = value
                centered = reshaped
            # Calculate FFT
            fourier_transformed = np.fft.fft(centered)  # FFT
            ft_real = np.real(fourier_transformed)
            ft_imag = np.imag(fourier_transformed)
            # Get amplitude of frequencies
            amplitude = ft_real * ft_real + ft_imag * ft_imag
            if normalize and max(amplitude) != 0:
                amplitude = np.true_divide(amplitude, max(amplitude))  # normalization of amplitude
            # Return only the first half of the FFT values:
            # DFT amplitudes are inherently mirrored passed
            # half the maximum frequency
            amplitude = amplitude[:math.ceil(len(amplitude) / 2)]
            result.append(np.round(amplitude, 6).tolist())
        if flatten:
            return list(chain.from_iterable(result))
        return result

    def pandas_get(self, sequences: List[str],
                   ids: Optional[List[str]] = None,
                   lag: int = 1,
                   domains: int = 2,
                   nproc: Optional[int] = None,
                   quiet=False,
                   ipynb=False, **kwargs) -> pd.DataFrame:
        """Get the raw values of the provided sequences in a pandas DataFrame.

        :param sequences: Protein sequences
        :param ids      : Unique identifiers of sequences
        :param lag      : Lag between amino acids
        :param domains  : Number of domains to split the sequence into
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
        values = pd.DataFrame(
            _multiprocess_get(self, sequences=sequences, ids=ids, nproc=nproc, ipynb=ipynb, quiet=quiet,
                              lag=lag, domains=domains))
        if self.Type is TransformType.AVG:
            info = f'domains{domains}'
        elif self.Type is TransformType.ACC or self.Type is TransformType.PDT:
            info = f'lag{lag}'
        elif self.Type is TransformType.FFT:
            info = ''
        if ids:
            values.columns = ['ID'] + [f'{self.Type.name}_{info}_{self.Descriptor.ID.replace(" ", "-")}_{x}' for x in
                                       range(1, len(values.columns))]
        else:
            values.columns = [f'{self.Type.name}_{info}_{self.Descriptor.ID.replace(" ", "-")}_{x}' for x in
                              range(1, len(values.columns) + 1)]
        return values
