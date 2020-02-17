import unittest
import pandas as pd
import numpy as np
import os
import tempfile
import shutil

from signatureanalyzer.signatureanalyzer import run_spectra
from signatureanalyzer.bnmf import ardnmf
from signatureanalyzer.utils import file_loader

SPECTRA_ARROW = "../../examples/example_luad_spectra_1.tsv"
SPECTRA_WORD = "../../examples/example_luad_spectra_2.tsv"

class TestMapping(unittest.TestCase):
    """
    Test Mapping
    """
    def test_sbs_cosmic2(self):
        """Test SBS Cosmic2"""
        dirpath = tempfile.mkdtemp()

        np.random.seed(0)
        spectra = file_loader(SPECTRA_ARROW)
        run_spectra(spectra, outdir=dirpath, cosmic='cosmic2', nruns=1, K0=10, max_iter=100, plot_results=False)
        cosine_df_arrow = pd.read_hdf(os.path.join(dirpath,'nmf_output.h5'),"cosine")
        shutil.rmtree(dirpath)

        ref_cosine = np.load("refs/test_mapping_cosmic2.npy")
        self.assertEqual(np.linalg.norm(ref_cosine - cosine_df_arrow.sum(1).values),0)

        np.random.seed(0)
        spectra = file_loader(SPECTRA_WORD)
        run_spectra(spectra, outdir=dirpath, cosmic='cosmic2', nruns=1, K0=10, max_iter=100, plot_results=False)
        cosine_df_word = pd.read_hdf(os.path.join(dirpath,'nmf_output.h5'),"cosine")
        shutil.rmtree(dirpath)

        self.assertEqual(np.linalg.norm(cosine_df_arrow.values - cosine_df_word.values),0)

    def test_sbs_cosmic3(self):
        """Test SBS Cosmic3"""
        dirpath = tempfile.mkdtemp()

        np.random.seed(0)
        spectra = file_loader(SPECTRA_ARROW)
        run_spectra(spectra, outdir=dirpath, cosmic='cosmic3', nruns=1, K0=10, max_iter=100, plot_results=False)
        cosine_df_arrow = pd.read_hdf(os.path.join(dirpath,'nmf_output.h5'),"cosine")
        shutil.rmtree(dirpath)

        ref_cosine = np.load("refs/test_mapping_cosmic3.npy")
        self.assertEqual(np.linalg.norm(ref_cosine - cosine_df_arrow.sum(1).values),0)

        np.random.seed(0)
        spectra = file_loader(SPECTRA_WORD)
        run_spectra(spectra, outdir=dirpath, cosmic='cosmic3', nruns=1, K0=10, max_iter=100, plot_results=False)
        cosine_df_word = pd.read_hdf(os.path.join(dirpath,'nmf_output.h5'),"cosine")
        shutil.rmtree(dirpath)

        self.assertEqual(np.linalg.norm(cosine_df_arrow.values - cosine_df_word.values),0)

    def test_sbs_cosmic3exome(self):
        """Test SBS Cosmic3 Exome"""
        dirpath = tempfile.mkdtemp()

        np.random.seed(0)
        spectra = file_loader(SPECTRA_ARROW)
        run_spectra(spectra, outdir=dirpath, cosmic='cosmic3_exome', nruns=1, K0=10, max_iter=100, plot_results=False)
        cosine_df_arrow = pd.read_hdf(os.path.join(dirpath,'nmf_output.h5'),"cosine")
        shutil.rmtree(dirpath)

        ref_cosine = np.load("refs/test_mapping_cosmic3_exome.npy")
        self.assertEqual(np.linalg.norm(ref_cosine - cosine_df_arrow.sum(1).values),0)

        np.random.seed(0)
        spectra = file_loader(SPECTRA_WORD)
        run_spectra(spectra, outdir=dirpath, cosmic='cosmic3_exome', nruns=1, K0=10, max_iter=100, plot_results=False)
        cosine_df_word = pd.read_hdf(os.path.join(dirpath,'nmf_output.h5'),"cosine")
        shutil.rmtree(dirpath)

        self.assertEqual(np.linalg.norm(cosine_df_arrow.values - cosine_df_word.values),0)

if __name__ == '__main__':
    unittest.main()
