import unittest
import pandas as pd
import numpy as np

from signatureanalyzer.spectra import get_spectra_from_maf
from signatureanalyzer.utils import file_loader

MAF_TEST_FILE = "../../examples/example_luad_maf.tsv"
HG_FILE = "../../examples/hg19.2bit"

class TestSpectra(unittest.TestCase):
    """
    Test Spectra Creation.
    """
    def test_cosmic_sbs(self):
        """
        Test single-base substituion mappers.
        """
        maf_df = file_loader(MAF_TEST_FILE)
        spectra_ref = pd.read_parquet("refs/test_cosmic2_spectra.parquet")

        _,spectra = get_spectra_from_maf(maf_df, cosmic='cosmic2', hgfile=HG_FILE)
        self.assertEqual(np.linalg.norm(spectra.values - spectra_ref.values),0)

    def test_cosmic_dbs(self):
        """
        Test doublet-base substitution mappers.
        """
        maf_df = file_loader(MAF_TEST_FILE)
        spectra_ref = pd.read_parquet("refs/test_cosmic3dbs_spectra.parquet")

        _,spectra = get_spectra_from_maf(maf_df, cosmic='cosmic3_DBS')
        self.assertEqual(np.linalg.norm(spectra.values - spectra_ref.values),0)

    def test_cosmic_id(self):
        """
        Test insertion-deletion mappers.
        """
        maf_df = file_loader(MAF_TEST_FILE)
        spectra_ref = pd.read_parquet("refs/test_cosmic3id_spectra.parquet")

        _,spectra = get_spectra_from_maf(maf_df, cosmic='cosmic3_ID', hgfile=HG_FILE)
        self.assertEqual(np.linalg.norm(spectra.values - spectra_ref.values),0)


if __name__ == '__main__':
    unittest.main()
