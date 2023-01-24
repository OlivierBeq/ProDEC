# -*- coding: utf-8 -*-


"""Tests for ProteinDescriptors."""

import unittest
from typing import *
from numbers import Number

from os import path

import prodec
from tests.constants import *


class TestFileExists(unittest.TestCase):
    """Tests for the presence of minimal data."""

    def test_data(self):
        """Test default ProteinDescriptor data is present."""
        self.assertTrue(path.isfile(DFLT_DATA))


class TestProteinDescriptors(unittest.TestCase):
    """Tests for ProteinDescriptors."""

    def setUp(self) -> None:
        """Load default ProteinDescriptor data"""
        self.pdescs = prodec.ProteinDescriptors()
        self.default_descs = ["Sneath", "HPI", "Kidera", 
                              "Zscale Sjöström", "Zscale Hellberg", 
                              "Zscale Jonsson", "Independent descriptors", 
                              "Combined descriptors", "GRID tscore", 
                              "ISA-ECI", "Contact energies", "Zscale Sandberg", 
                              "Raychaudhury", "MS-WHIM", "E-scale", "VSTV", 
                              "c-scales", "VHSE", "PSM", "SSIA AM1", 
                              "GH-scale", "SSIA PM3", "SSIA HF", "SSIA DFT", 
                              "FASGAI", "Tscale", "VSGETAWAY", "SVRDF", "VTSA", 
                              "SZOTT", "V-scale", "HSEHPCSV", "VSW", "VHSEH", 
                              "BLOSUM", "VARIMAX", "HESH", "DPPS", "STscale", 
                              "CBFQ", "CDFQ", "CUFQ", "ADFQ", "SVEEVA", "SVRG", 
                              "G-scales", "ProtFP hash", "P-scale", "SVWG", "VSTPV", 
                              "QCP", "ProtFP PCA", "Zscale binary", "SVMW", "SVHEHS", "SVGER"]

    def test_proteindescriptors_loading_type(self):
        """Test ProteinDescriptors loads default data properly"""
        self.assertIsInstance(self.pdescs.available_descriptors, List)
    
    def test_proteindescriptors_loading_size(self):
        self.assertGreater(len(self.pdescs.available_descriptors), 0)
    
    def test_proteindescriptors_loading_value(self):
        self.assertListEqual(sorted(self.pdescs.available_descriptors), sorted(self.default_descs))
    
    def test_proteindescriptors_descriptor(self):
        for desc in self.default_descs:
            self.assertEqual(desc, self.pdescs.get_descriptor(desc).ID)

    def test_proteindescriptors_scales(self):
        """Test descriptor sizes are loaded properly""" 
        for desc in self.pdescs.available_descriptors:
            desc = self.pdescs.get_descriptor(desc)
            dtype = desc.Type
            values = desc.Scales_values
            self.assertIsInstance(values, dict)
            self.assertEqual(len(values.keys()), 20)
            if dtype == 'Linear':
                for key, value in values.items():
                    self.assertIn(key, DFLT_AA)
                    self.assertTrue(isinstance(value, list) or isinstance(value, Number))
            elif dtype == 'Distance':
                for key1, value1 in values.items():
                    self.assertIn(key1, DFLT_AA)
                    self.assertIsInstance(value1, dict)
                    for key2, value2 in value1.items():
                        self.assertIn(key2, DFLT_AA)
                        self.assertIsInstance(value2, Number)
            else:
                self.assertEqual(dtype, 'Other')
                self.assertEqual(desc.ID, 'Raychaudhury')
                for key, value in values.items():
                    self.assertIn(key, DFLT_AA)
                    self.assertIsInstance(value, Number)
    
    def test_proteindescriptors_info(self):
        """Test if info is correctly given for descriptor"""
        for desc in self.pdescs.available_descriptors:
            desc = self.pdescs.get_descriptor(desc)
            self.assertGreater(len(desc.Info['Authors']), 0)
            self.assertIsNotNone(desc.Info['Year'])
            self.assertEqual(len(str(desc.Info['Year'])), 4)
            self.assertIsInstance(desc.Info['Year'], Number)
            # Either Journal or Patent and if journal, eihter DOI or PMID
            self.assertTrue((desc.Info['Journal'] is not None and
                             (desc.Info['DOI'] is not None or
                              desc.Info['PMID'] is not None)
                             ) or desc.Info['Patent'] is not None)

