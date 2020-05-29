import numpy as np
from sys import getsizeof
from psutil import virtual_memory
from typing import Callable, List, Type, Tuple
from numbers import Number

std_amino_acids = list('ACDEFGHIKLMNPQRSTVWY')

def alphaC_distance_edge_vector_speed(sequence: str, 
                                      mapping: Callable[[str], Number], 
                                      distance: Number=1,
                                      power: Number=-4,
                                      dtype: Type = np.float32,
                                      window: int = 60,
                                      offset: int=1
                                     ) -> List[Number] :
    '''
    Fast method to calculate alphaCarbon distance-normalized
    edge vector from mapping values
    
    :param sequence : protein sequence
    :param mapping  : function mapping amino acid to scale value
    :param distance : distance between alpha carbons
    :param power    : power scaling of distance matrix
    :param dtype    : value type for memory efficiency
    :param window   : sliding window in which amino acids are considered
    :param offset   : distance value offset between an alphaC and itself
    '''
    # Create distance edge matrix
    len_ = len(sequence)
    deM = np.zeros((len_, len_), dtype=dtype)
    # Create distances 
    r_ = np.power(np.abs(np.array(range(-len_,len_), dtype=dtype)) * distance + offset, power)
    # Fill in distance edge matrix
    for i in range(len_):
        strt_window = max(i-window, 0)
        end_window = min(i+window, len_)
        strt_roll = -len_+strt_window
        end_roll = -len_+end_window if end_window < len_ else None
        deM[i, strt_window:end_window] = np.roll(r_, i)[strt_roll:end_roll]
    # Create weight matrix
    wM = np.array([mapping(sequence[i]) for i in range(len_)], dtype=dtype)
    res = np.dot(deM, wM)
    return res

def alphaC_distance_edge_vector_memory(sequence: str, 
                                       mapping: Callable[[str], Number], 
                                       distance: Number=1,
                                       power: Number=-4,
                                       dtype: Type = np.float32,
                                       window: int = 60,
                                       offset: int = 1
                                      ) -> List[Number] :
    '''
    Memory efficient method to calculate
    alphaCarbon distance-normalized
    edge vector from mapping values
    
    :param sequence : protein sequence
    :param mapping  : function mapping amino acid to scale value
    :param distance : distance between alpha carbons
    :param power    : power scaling of distance matrix
    :param dtype    : value type for memory efficiency
    :param window   : sliding window in which amino acids are considered
    :param offset   : distance value offset between an alphaC and itself
    '''
    # Create distance edge vector
    len_ = len(sequence)
    dists = np.zeros((len_, 1), dtype=dtype)
    for i in range(len_):
        dists[i, 0] = np.power(i * distance + offset, power)
    # Calculate indices
    values = np.zeros((len_, 1), dtype=dtype)
    for i in range(len_):
        raw_value = mapping(sequence[i])
        for j in range(max(i-window, 0), min(i+window, len_)):
            values[j] += raw_value * dists[abs(j -i)]
    return values

def raychaudhury_speed(sequence: str,
                       mapping: Callable[[str], Number],
                       power: int,
                       dtype: Type = np.float32
                      ) -> List[Number]:
    '''
    Fast  method to calculate
    Raychaudhury's distance-normalized
    edge vector from mapping values
    
    :param sequence : protein sequence
    :param mapping  : function mapping amino acid to scale value
    :param distance : distance between alpha carbons
    :param power    : power scaling of distance matrix
    :param dtype    : value type for memory efficiency
    '''
    # Create distance edge matrix
    len_ = len(sequence)
    deM = np.zeros((len_, len_), dtype=dtype)
    # Create distances 
    r_ = np.power(np.abs(np.array(range(-len_,len_), dtype=dtype)) + 1, power)
    # Fill in distance edge matrix
    for i in range(len_):
        deM[i, :] = np.roll(r_, i)[-len_:]
    # Multiply by Raychaudhury's scale (weight matrix)
    wM = np.array([mapping(sequence[i]) for i in range(len(sequence))], dtype=dtype)
    res = np.dot(deM, wM)
    return res

def raychaudhury_memory(sequence: str,
                        mapping: Callable[[str], Number],
                        power: int,
                        dtype: Type = np.float32
                       ) -> List[Number]:
    '''
    Memory efficient method to calculate
    Raychaudhury's distance-normalized
    edge vector from mapping values
    
    :param sequence : protein sequence
    :param mapping  : function mapping amino acid to scale value
    :param distance : distance between alpha carbons
    :param power    : power scaling of distance matrix
    :param dtype    : value type for memory efficiency
    '''
    # Create distance edge vector
    len_ = len(sequence)
    dists = np.zeros((len_, 1), dtype=dtype)
    for i in range(1, len_):
        dists[i, 0] = np.power(i, power)
    # Calculate indices
    values = np.zeros((len_, 1), dtype=dtype)
    for i in range(len_):
        raw_value = mapping(sequence[i])
        for j in range(len_):
            values[j] += raw_value * dists[abs(j -i)]
        print(i)
    return values

def get_values(scale_values: dict):
    '''
    Get recursively the inner most values from a
    dictionnary of dictionaries of ... of dictionaries
    '''
    values= []
    keys = scale_values.keys()
    for key in scale_values.keys():
        if isinstance(scale_values[key], dict):
            values.extend(get_values(scale_values[key]))
        elif isinstance(scale_values[key], list):
            values.extend(scale_values[key])
        else:
            values.append(scale_values[key])
    return values

def enough_avail_memory(array_size: int, dtype: Type) -> bool:
    '''
    Whether enough memory is available to  use fast 
    in memory computation of distance edge vectors
    
    :param array_size : size of the sequence
    '''
    e = lambda x : 1 * x ** 2 + 112
    f = lambda x : 2 * x ** 2 + 112
    g = lambda x : 4 * x ** 2 + 112
    h = lambda x : 8 * x ** 2 + 112
    avail_ram = virtual_memory().available
    obj_size = getsizeof(dtype(1))
    if obj_size <= 25:
        return e(array_size) * 1.5 < avail_ram
    if obj_size == 26:
        return f(array_size) * 1.5 < avail_ram
    if obj_size == 28:
        return g(array_size) * 1.5 < avail_ram
    if obj_size >= 32:
        return h(array_size) * 1.5 < avail_ram
    return False
