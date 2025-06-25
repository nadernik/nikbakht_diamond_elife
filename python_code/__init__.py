"""
Rat Vision Analysis Package

A Python implementation of the rat vision analysis from Nikbakht & Diamond (2021).
"""

__version__ = "1.0.0"
__author__ = "Nader Nikbakht"
__email__ = "nader.nikbakht@sissa.it"

from .beautify import beautify
from .utils import load_matlab_data, cumulative_gaussian_lapse, fit_psychometric_function

__all__ = [
    'beautify',
    'load_matlab_data', 
    'cumulative_gaussian_lapse',
    'fit_psychometric_function'
] 