import math
import numpy as np
from .Descriptor import Descriptor


TRANSFORMS = {'AVG': {'Fullname': 'Domain average', 'Constant Length': True, 'Binary': False},
              'ACC': {'Fullname': 'Auto-cross covariances', 'Constant Length': True, 'Binary': False},
              'PDT': {'Fullname': 'Physicochemical Distance Tranformation', 'Constant Length': True, 'Binary': False}
             }

class Transform(Descriptor):
    '''
    A class used to process raw values of protein descriptors
    '''
    def __init__(self, type: str, descriptor: Descriptor):
        self.Type = type
        self.Descriptor = descriptor
        self.Binary = TRANSFORMS[type]['Binary']
    
    @staticmethod
    def available_transforms():
        return list(TRANSFORMS.keys())
    
    @property
    def summary(self):
        return TRANSFORMS[self.Type]
    
    @staticmethod
    def is_compatible(transform_name: str, descriptor: Descriptor) -> bool:
        if transform_name not in TRANSFORMS.keys():
            raise NotImplemented(f'Transform type {self.Type} is not implemented')
        return descriptor.Binary == TRANSFORMS[transform_name]['Binary'] or TRANSFORMS[transform_name]['Binary'] == 'All'
    
    def get(self, sequence: str, 
                  flatten: bool =True,
                  lag: int =0, 
                  domains: int =0, 
                  average: bool =True, 
                  **kwargs):
        '''
        Transform raw protein descriptor values.
        
        :param sequence : raw descriptor values
        :param flatten  : not to return dimensions separate, only for AVG and ACC
        :param lag      : lag value
        :param domains  : number of averaged domains
        :param average  : global average of dimensions, only for PDT
        '''
        if not Transform.is_compatible(self.Type, self.Descriptor):
            raise Exception(f'Type of transform ({self.Type}) is not compatible'
                            f' with descriptor type ({self.Descriptor.Type})')
        if self.Type == 'AVG':
            return self.__average__(sequence, domains, flatten, **kwargs)
        elif self.Type == 'ACC':
            return self.__autocrosscovariance__(sequence, lag, flatten, **kwargs)
        elif self.Type == 'PDT':
            return self.__physicochem_distance_tansform__(sequence, lag, average, **kwargs)
        else:
            raise NotImplemented(f'Transform type {self.Type} is not implemented')
    
    def __average__(self, sequence: str, domains: int, flatten: bool, **kwargs):
        '''
        Calculate domain mean average values.
        If flatten, mean averages of each dimension are added,
        otherwise, global average (all dimensions) is added.
        '''
        length = len(sequence)
        if domains < 2 or domains >= length:
            raise Exception(f'Number of domains ({domains}) has to be greater or equal to 1'
                            f' and lower or equal than the length of the sequence ({length})')
        # Keep track of original shape if not flat result
        if not flatten:
            raw = np.array(self.Descriptor.get(sequence, flatten=False, **kwargs))
            shape = raw.shape
        raw = self.Descriptor.get(sequence, flatten=True, **kwargs)
        average = np.zeros((domains * self.Descriptor.Size), dtype=np.float64)
        counts = np.zeros((domains * self.Descriptor.Size), dtype=np.int64)
        scaling_factor = ((len(raw) / self.Descriptor.Size) + 1) / domains
        for i in range(len(raw)):
            # Get scaled position
            seq_pos = int(i / self.Descriptor.Size)
            desc_pos = i % self.Descriptor.Size
            scaled_pos = math.floor((seq_pos+1)/scaling_factor) * self.Descriptor.Size + desc_pos
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
                weights = np.identity(self.Descriptor.Size)[row,:].tolist() * domains
                average = np.append(average, np.average(average[:domains*self.Descriptor.Size], weights=weights))
            average = np.reshape(average, (domains+1, shape[1])).tolist()
        else:
            # Add global mean average
            average = np.append(np.round(average, 6), np.average(np.round(average, 6))).tolist()
        return average
    
    def __autocrosscovariance__(self, sequence: str, lag: int, flatten: bool, **kwargs):
        length = len(sequence) - 1 if self.Descriptor.Type == 'Distance' else len(sequence)
        if lag >= length:
            if self.Descriptor.Type != 'Distance':
                raise Exception(f'Lag ({lag}) has to be greater or equal to 1 and lower'
                                f' than the length of the sequence ({length})')
            else:
                raise Exception(f'Lag ({lag}) has to be greater or equal to 1 and lower'
                                f' than the length of the sequence minus one ({length-1})')
        acc = np.zeros((self.Descriptor.Size, self.Descriptor.Size), dtype=np.float64)
        raw = self.Descriptor.get(sequence, flatten=True, **kwargs)
        for j in range(self.Descriptor.Size):
            for m in range(self.Descriptor.Size):
                for i in range(length - lag):
                    acc[m, j] += (raw[i*self.Descriptor.Size+j] * raw[(i+lag)*self.Descriptor.Size+m]) / (length-lag)
        if flatten:
            return acc.flatten(order='F')
        return acc
    
    def __physicochem_distance_tansform__(self, sequence: str, lag: int, average: bool, **kwargs):
        length = len(sequence) - 1 if self.Descriptor.Type == 'Distance' else len(sequence)
        if lag >= length:
            if self.Descriptor.Type != 'Distance':
                raise Exception(f'Lag ({lag}) has to be greater or equal to 1 and lower'
                                f' than the length of the sequence ({length})')
            else:
                raise Exception(f'Lag ({lag}) has to be greater or equal to 1 and lower'
                                f' than the length of the sequence minus one ({length-1})')
        # Normalisation of scores
        stdScores = np.array(self.Descriptor.get('ACDEFGHIKLMNPQRSTVWY', flatten=False))
        meanScore = np.average(stdScores, axis=0)
        meanDifference = np.sqrt(np.sum(np.power(stdScores - meanScore, 2), axis=0))
        # Adjust raw values
        raw = np.array(self.Descriptor.get(sequence, flatten=False, **kwargs))
        raw = (raw - meanScore) / meanDifference
        # Taking gaps into account
        if ('gaps', 'omit') in kwargs.items():
            sequence = ''.join(filter(str.isalpha, sequence))
        # Get PDT
        pdt = np.zeros(self.Descriptor.Size)
        for i in range(length - lag):
            if self.Descriptor.Size == 1:
                pdt[0] += math.pow(raw[i] - raw[i+lag], 2)
            else:
                for m in range(self.Descriptor.Size):
                    pdt[m] += math.pow(raw[i,m] - raw[i+lag,m], 2)
        pdt /= length - lag
        if average:
            return np.average(pdt)
        return pdt
