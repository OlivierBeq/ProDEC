import orjson
import collections
from os import path
from .Descriptor import Descriptor

class ProteinDescriptors:
    '''
    A class used for caching available protein descripors 
    '''
    def __init__(self):
        self.update_cache(path.join(path.dirname(__file__), 'data.json'))
        self.cache = {}

    @property
    def available_descriptors(self):
        return sorted([descriptor['ID'] if 'ID' in descriptor.keys() 
                                        else descriptor['Name'] 
                       for descriptor in self.descriptors], key=str.casefold)
    
    def update_cache(self, cache):
        with open(cache, 'rb') as input:
            descs = input.read().decode('utf-8')
        self.descriptors = orjson.loads(descs)
        self.__check_uniqueness__()
    
    def __check_uniqueness__(self):
        ids = (descriptor['ID'] for descriptor in self.descriptors)
        if len(list(filter(lambda x: x[1] > 1, collections.Counter(ids).items()))) != 0:
            raise Exception('Non unique descriptor ID')
    
    def get_descriptor(self, id):
        if id not in self.cache.keys():
            for descriptor in self.descriptors:
                if descriptor['ID'] == id:
                    self.cache[id] = Descriptor(descriptor)
        return self.cache[id]
