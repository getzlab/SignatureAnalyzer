# Modules
from . import plotting as pl
from . import pathways as pw

# Independent imports
from .bnmf import ardnmf
from .utils import postprocess_msigs
from .consensus import consensus_cluster

from .signatureanalyzer import run_maf
from .signatureanalyzer import run_spectra
from .signatureanalyzer import run_matrix

__version__ = '0.0.7'
