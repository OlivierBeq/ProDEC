# -*- coding: utf-8 -*-

"""Utility functions for ProDEC."""

from numbers import Number
from sys import getsizeof
from typing import Callable, List, Type

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
    obj_size = getsizeof(dtype(1))
    if obj_size <= 26:
        return f(array_size) * (1 + margin) < avail_ram
    if obj_size <= 28:
        return g(array_size) * (1 + margin) < avail_ram
    else: # usually with obj_size >= 32:
        return h(array_size) * (1 + margin) < avail_ram
    return False
