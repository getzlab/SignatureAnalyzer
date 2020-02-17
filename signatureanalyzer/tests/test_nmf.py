import unittest
import pandas as pd
import numpy as np
import os

from signatureanalyzer.signatureanalyzer import run_spectra
from signatureanalyzer.bnmf import ardnmf
from signatureanalyzer.utils import file_loader

SPECTRA_ARROW = "../../examples/example_luad_spectra_1.tsv"
SPECTRA_WORD = "../../examples/example_luad_spectra_2.tsv"

LAM_TRUE_1 = [0.7218399,0.41164818,1.0229105,0.5949916,0.587698, 1.0327462,
    0.48455864,0.8585669,0.5379705,1.4745888]

class TestNmf(unittest.TestCase):
    """
    Test ARD-NMF.
    """
    def test_ardnmf_poisson(self):
        spectra = file_loader(SPECTRA_ARROW)

        np.random.seed(0)
        res = ardnmf(spectra, K0=10, max_iter=50, verbose=False)

        for i,_ in enumerate(LAM_TRUE_1):
            self.assertAlmostEqual(list(res['lam'])[i], LAM_TRUE_1[i], places=7)

if __name__ == '__main__':
    unittest.main()
