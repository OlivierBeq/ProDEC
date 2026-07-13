# -*- coding: utf-8 -*-

"""Default protein descriptor loader."""

import collections
from os import path
from typing import Any, BinaryIO, Dict, List, Optional, Union

import orjson

from .Descriptor import Descriptor


class ProteinDescriptors:
    """A class used for caching available protein descripors."""

    def __init__(self, json: Optional[str] = None) -> None:
        """Instantiate a new ProteinDesciptors.

        :param json: Path to json file or file handle describing descriptors
        """
        self.descriptors: List[Dict[str, Any]] = []
        if json is None:
            self.update_available_descriptors(path.join(path.dirname(__file__),
                                                        'data.json'))
        else:
            if not path.isfile(json):
                raise Exception(f'File {json} does not exist')
            else:
                self.update_available_descriptors(json)
        self.cache: Dict[str, Descriptor] = {}

    @property
    def available_descriptors(self) -> List[str]:
        """Give descriptors loaded for calculation."""
        return sorted((descriptor['ID'] if 'ID' in descriptor.keys()
                       else descriptor['Name']
                       for descriptor in self.descriptors), key=str.casefold)

    def update_available_descriptors(self, json: Union[str, BinaryIO]) -> None:
        """Read descriptor file and update available descriptors.

        :param json: Path to json file or file handle describing descriptors
        """
        if isinstance(json, str):
            with open(json, 'rb') as input:
                descs = input.read().decode('utf-8')
        else:
            descs = json.read().decode('utf-8')
        self.descriptors = orjson.loads(descs)
        self.__check_uniqueness__()

    def __check_uniqueness__(self) -> None:
        """Check IDs of descriptors are unique."""
        ids = (descriptor['ID'] for descriptor in self.descriptors)
        non_unique_ids = filter(lambda x: x[1] > 1,
                                collections.Counter(ids).items())
        if len(list(non_unique_ids)) != 0:
            raise Exception('Non unique descriptor ID')

    def get_descriptor(self, id: str) -> Descriptor:
        """Get Descriptor instance from ID.

        :param id: ID or the descriptor
        """
        if id not in self.cache.keys():
            for descriptor in self.descriptors:
                if descriptor['ID'] == id:
                    self.cache[id] = Descriptor(descriptor)
        return self.cache[id]
