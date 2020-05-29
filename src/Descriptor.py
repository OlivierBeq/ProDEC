import itertools
import numpy as np
from .utils import *

class Descriptor:
    '''
    A class used for calculating protein descriptors
    '''
    def __init__(self, desc_data):
        self.ID = desc_data['ID']
        self.Type = desc_data['Type']
        self.Name = desc_data['Name']
        self.Scales_names = desc_data['Scales']['Names']
        self.Scales_values = desc_data['Scales']['Values']
        self.Info = {'Authors':desc_data['Authors'], 'Year':desc_data['Year'],
                     'Journal':desc_data['Journal'], 'DOI':desc_data['DOI'],
                     'PMID':desc_data['PMID'], 'Patent':desc_data['Patent']}
        first_val = list(self.Scales_values.values())[0]
        if isinstance(first_val, list):
            self.Size = len(first_val)
        elif isinstance(first_val, dict):
            self.Size = 1
        else:
            self.Size = 1
        self.Binary = set(get_values(self.Scales_values)) == {0, 1}
    
    
    @property
    def summary(self):
        return self.Info
    
    @property
    def definition(self):
        return self.Scales
    
    def is_sequence_valid(self, sequence: str):
        for letter in sequence:
            if letter not in self.Scales_values.keys():
                return False
        return True
    
    def get(self, sequence: str, flatten=True, gaps=0, dictionary=std_amino_acids, dtype=np.float64, prec=60, power=-4, **kwargs):
        '''
        Get the raw values of the provided sequence.
        
        :param sequence  : protein sequence
        :param flatten   : not to return dimensions separate
        :param gaps      : how should gaps be considered.
                           Allowed values: 'omit' or 0, ...+inf
        :param dictionary: list of valid letters
        '''
        # Dealing with gaps
        if gaps == 'omit':
            sequence = ''.join(filter(str.isalpha, sequence))
            replacer = 0
        else:
            sequence = ''.join([x if str.isalpha(x) else '-' for x in sequence])
            replacer = gaps
            if self.Type == 'Linear':
                self.Scales_values['-'] = [replacer] * self.Size
            elif self.Type == 'Distance':
                self.Scales_values['-'] = {}
                for key, value in self.Scales_values.items():
                    self.Scales_values[key]['-'] = value
                    self.Scales_values['-'][key] = value
            else:
                self.Scales_values['-'] = replacer
        # Checking sequence
        if not self.is_sequence_valid(sequence):
            raise Exception('Sequence has unsupported amino acid')
        # Create final array
        if self.Size > 1:
            dtype = int if self.Binary else float
            values = np.zeros((len(sequence), self.Size), dtype=dtype)
            values.fill(replacer)
        elif self.Type == 'Distance':
            values = [replacer] * (len(sequence) -1)
        else:
            values = [replacer] * len(sequence)
        if self.Type == 'Linear':
            for i in range(len(sequence)):
                if self.Size > 1:
                    for j in range(self.Size):
                        values[i, j] = self.Scales_values[sequence[i]][j]
                else:
                    values[i] = self.Scales_values[sequence[i]]
        elif self.Type == 'Distance':
            for i in range(len(sequence) - 1):
                values[i] = self.Scales_values[sequence[i]][sequence[i+1]]
        else:
            if self.ID == 'Raychaudhury':
                mapping = lambda x: self.Scales_values[x]
                if not prec:
                    if enough_avail_memory(len(sequence), dtype):
                        return raychaudhury_speed(sequence, mapping, power, dtype)
                    return raychaudhury_memory(sequence, mapping, power, dtype)
                else:
                    if enough_avail_memory(len(sequence), dtype):
                        return alphaC_distance_edge_vector_speed(sequence, mapping, 1, power, dtype, prec, 1)
                    return alphaC_distance_edge_vector_memory(sequence, mapping, 1, power, dtype, prec, 1)
        if isinstance(values, list):
            return values
        if not flatten:
            return values.tolist()
        else:
            return values.flatten(order='C').tolist()
