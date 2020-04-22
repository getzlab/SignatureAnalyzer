# Modules
from . import plotting as pl
from . import pathways as pw
from . import spectra as spectra

# Independent imports
from .bnmf import ardnmf
from .semi_supervised_bnmf import ss_ardnmf
from .utils import postprocess_msigs
from .consensus import consensus_cluster

from .signatureanalyzer import run_maf
from .signatureanalyzer import run_spectra
from .signatureanalyzer import run_matrix

__version__ = '0.0.7'
