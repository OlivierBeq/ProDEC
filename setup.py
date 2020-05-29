# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:00:00 2020

@author: OliBeq
"""

import setuptools

setuptools.setup(
    name="proDEC", 
    version="0.1",
    author="Olivier J. M. Béquignon",
    author_email="olivier.bequignon.maintainer@gmail.com",
    description="A package to calculate protein sequence descriptors",
    url="https://gitlab.services.universiteitleiden.nl/cdd/prodec",
    packages = ['proDEC'],
    package_dir={'proDEC':'src'},
    package_data={'proDEC':['*.json']},
    install_requires=[
        'orjson',
        'numpy',
        'psutil'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
