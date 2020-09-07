#!/bin/env python3

"""Setup.py for the HistAstro Python package."""


# Package version:
version="0.0.3"

# Get long description from README.md:
with open("README.md", "r") as fh:
    long_description = fh.read()

# Set package properties:
from setuptools import setup
setup(
    name='histastro',
    description='A Python package for historical-astronomy calculations of Sun, Moon and planets',
    author='Marc van der Sluys',
    url='http://astro.ru.nl/~sluys/HistAstro',
    
    packages=['histastro'],
    install_requires=['numpy','fortranformat'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    
    version=version,
    license='GPLv3+',
    keywords=['astronomy','history','ephemeris','sun','moon','planets'],
    
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Sociology :: History",
    ]
)

