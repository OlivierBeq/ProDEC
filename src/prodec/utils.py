# -*- coding: utf-8 -*-

"""Utility functions for ProDEC."""

import functools
import sys
from collections import deque
from itertools import islice
from numbers import Number
from multiprocessing import cpu_count, Manager, Pool
from typing import Callable, List, Type, Union

import numpy as np
from psutil import virtual_memory


std_amino_acids = list('ACDEFGHIKLMNPQRSTVWY')


def alpha_carbon_distance_edge_vector_speed(sequence: str,
                                            mapping: Callable[[str], Number],
                                            distance: Number = 1,
                                            power: Number = -4,
                                            dtype: Type = np.half,
                                            window: int = 60,
                                            offset: int = 1) -> List[Number]:
    """Calculate alphaCarbon distance-normalized edge vector.

    This implementation is fast but requires a lot of memory.

    :param sequence : protein sequence
    :param mapping  : function mapping amino acid to scale value
    :param distance : distance between alpha carbons
    :param power    : power scaling of distance matrix
    :param dtype    : value type for memory efficiency
    :param window   : sliding window in which amino acids are considered
    :param offset   : distance value offset between an alphaC and itself
    """
    # Create distance edge matrix
    len_ = len(sequence)
    de_mat = np.zeros((len_, len_), dtype=dtype)
    # Create distances
    r_ = np.power(np.abs(np.array(range(-len_, len_), dtype=dtype)) *
                  distance + offset, power)
    # Fill in distance edge matrix
    for i in range(len_):
        strt_window = max(i-window, 0)
        end_window = min(i+window, len_)
        strt_roll = -len_+strt_window
        end_roll = -len_+end_window if end_window < len_ else None
        de_mat[i, strt_window:end_window] = np.roll(r_, i)[strt_roll:end_roll]
    # Create weight matrix
    w_mat = np.array([mapping(sequence[i]) for i in range(len_)], dtype=dtype)
    res = np.dot(de_mat, w_mat)
    return res


def alpha_carbon_distance_edge_vector_memory(sequence: str,
                                             mapping: Callable[[str], Number],
                                             distance: Number = 1,
                                             power: Number = -4,
                                             dtype: Type = np.half,
                                             window: int = 60,
                                             offset: int = 1) -> List[Number]:
    """Calculate alphaCarbon distance-normalized edge vector.

    This implementation uses little memory but is slow.

    :param sequence : protein sequence
    :param mapping  : function mapping amino acid to scale value
    :param distance : distance between alpha carbons
    :param power    : power scaling of distance matrix
    :param dtype    : value type for memory efficiency
    :param window   : sliding window in which amino acids are considered
    :param offset   : distance value offset between an alphaC and itself
    """
    # Create distance edge vector
    len_ = len(sequence)
    dists = np.zeros((len_, 1), dtype=dtype)
    for i in range(len_):
        dists[i, 0] = np.power(dtype(i) * distance + offset, power)
    # Calculate indices
    values = np.zeros((len_, 1), dtype=dtype)
    for i in range(len_):
        raw_value = mapping(sequence[i])
        for j in range(max(i-window, 0), min(i+window, len_)):
            values[j] += raw_value * dists[abs(j - i)]
    return values


def raychaudhury_speed(sequence: str,
                       mapping: Callable[[str], Number],
                       power: int,
                       dtype: Type = np.half) -> List[Number]:
    """Calculate Raychaudhury's distance-normalized edge vector.

    This implementation is fast but requires a lot of memory.
    :param sequence : protein sequence
    :param mapping  : function mapping amino acid to scale value
    :param distance : distance between alpha carbons
    :param power    : power scaling of distance matrix
    :param dtype    : value type for memory efficiency
    """
    # Create distance edge matrix
    len_ = len(sequence)
    de_mat = np.zeros((len_, len_), dtype=dtype)
    # Create distances
    r_ = np.power(np.abs(np.array(range(-len_, len_), dtype=dtype)) + 1, power)
    # Fill in distance edge matrix
    for i in range(len_):
        de_mat[i, :] = np.roll(r_, i)[-len_:]
    # Multiply by Raychaudhury's scale (weight matrix)
    w_mat = np.array([mapping(sequence[i]) for i in range(len(sequence))],
                     dtype=dtype)
    res = np.dot(de_mat, w_mat)
    return res


def raychaudhury_memory(sequence: str,
                        mapping: Callable[[str], Number],
                        power: int,
                        dtype: Type = np.half) -> List[Number]:
    """Calculate Raychaudhury's distance-normalized edge vector.

    This implementation uses little memory but is slow.
    :param sequence : protein sequence
    :param mapping  : function mapping amino acid to scale value
    :param distance : distance between alpha carbons
    :param power    : power scaling of distance matrix
    :param dtype    : value type for memory efficiency
    """
    # Create distance edge vector
    len_ = len(sequence)
    dists = np.zeros((len_, 1), dtype=dtype)
    for i in range(len_):
        dists[i, 0] = np.power(dtype(i + 1), power)
    # Calculate indices
    values = np.zeros((len_, 1), dtype=dtype)
    for i in range(len_):
        raw_value = mapping(sequence[i])
        for j in range(len_):
            values[j] += raw_value * dists[abs(j - i)]
    return values


def get_values(scale_values: dict):
    """Get recursively the inner most values from a recursive dict."""
    values = []
    for key in scale_values.keys():
        if isinstance(scale_values[key], dict):
            values.extend(get_values(scale_values[key]))
        elif isinstance(scale_values[key], list):
            values.extend(scale_values[key])
        else:
            values.append(scale_values[key])
    return values


def enough_avail_memory(array_size: int, dtype: Type, margin: float=0.1) -> bool:
    """Check if enough memory is available.

    This is to use fast in memory computation of distance edge vectors.
    :param array_size : size of the sequence
    :param dtype: uderlying data type of the array
    :param margin: makes sure that the RAM available is at least (1 + margin)
                   times the predicted RAM required to hold the array.
    """
    def f(x): return 2 * x ** 2 + 112
    def g(x): return 4 * x ** 2 + 112
    def h(x): return 8 * x ** 2 + 112
    avail_ram = virtual_memory().available
    obj_size = sys.getsizeof(dtype(1))
    if obj_size <= 26:
        return f(array_size) * (1 + margin) < avail_ram
    if obj_size <= 28:
        return g(array_size) * (1 + margin) < avail_ram
    else: # usually with obj_size >= 32:
        return h(array_size) * (1 + margin) < avail_ram
    return False


class DescriptorPool:
    """Multiprofcessing class to calculate descriptors and transforms."""
    def __init__(self, calc,
                 nproc: int, **kwargs):
        """Instantiate a DescriptorPool.
        
        :param calc: a prodec.Descriptor or prodec.Transform object
        :param nproc: number of concurrent processes
        :param kwargs: any prodec.Descriptor or prodec.Transform parameter
        """
        self.pool = Pool(nproc)
        self.mgr = Manager()
        self.calc = functools.partial(calc.get, **kwargs)
        self.nproc = nproc

    def __enter__(self):
        """Start iteration."""
        self.mgr.__enter__()
        return self

    def __exit__(self, *args, **kwargs):
        """Terminate iteration."""
        self.pool.terminate()
        self.mgr.__exit__(*args, **kwargs)

    def map(self, seqs):
        """Distribute calculation for sequences to a PoolIterator.
        
        :param seqs: protein sequences
        """
        return PoolIterator(self, seqs, self.nproc * 2 + 10)

    def submit(self, seq):
        """Submit calculation to the internal Pool.

        :param seq: protein sequence
        """
        return self.pool.apply_async(self.calc, (seq,))


class PoolIterator:
    """Multiprocessing iterator to be used with DescriptorPool"""
    def __init__(self, pool: DescriptorPool, seqs: List[str], buf: int):
        """Instantiate a PoolIterator.
        
        :param pool: the DescriptorPool to submit calculations to
        :param seqs: list of protein sequences
        :param buf: buffer of sequences per process
        """
        self.pool = pool
        self.futures = deque()
        self.seqs = zip(range(len(seqs)), seqs)

        for id, seq in islice(self.seqs, buf):
            self.submit(id, seq)

    def submit(self, id, seq):
        """Add sequence id and future to internal futures.
        
        :param id: ID of the sequence
        :param seq: protein sequence
        """
        self.futures.append((id, self.pool.submit(seq)))

    def __iter__(self):
        """Iterate to create a container."""
        return self

    def __next__(self):
        """Submit job and return."""
        try:
            id, seq = next(self.seqs)
            self.submit(id, seq)
        except StopIteration:
            pass

        try:
            id, fut = self.futures.popleft()
            return id, fut.get()
        except IndexError:
            raise StopIteration

    next = __next__


def _multiprocess_get(descriptor,  # a Descriptor or Transform object
                      sequences: List[str],
                      nproc: int=1,
                      ipynb: bool=False,
                      quiet:bool=False, **kwargs):
    """Calculate protein descriptors or transforms and return a pandas dataframe.
    
    :param descriptor: prodec.Descriptor or prodec.Transform
    :param sequences: protein sequences
    :param nproc: number of concurrent processes
    :param ipynb: whether it is used in a notebook
    :param quiet: whether to show progress or not
    """
    if not isinstance(nproc, int) or nproc < 1:
        nproc = cpu_count()
    if ipynb:
        import tqdm.notebook
        pbar = tqdm.notebook.tqdm(total=len(sequences), disable=quiet)
    else:
        import tqdm
        pbar = tqdm.tqdm(total=len(sequences), disable=quiet)
    with DescriptorPool(descriptor, nproc, **kwargs) as pool:
        for seq, res in pool.map(sequences):
            yield res
            if not quiet:
                pbar.update()
    pbar.close()