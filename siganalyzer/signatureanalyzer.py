import sys
import argparse
import os
import pkg_resources
import pandas as pd

from typing import Union
from .spectra import get_spectra_from_maf
from .bnmf import ardnmf

def run_maf(
    maf: Union[str, pd.DataFrame],
    outdir: str = '.',
    cosmic: str = 'cosmic2',
    hg_build: Union[str, None] = 'hg19',
    nruns: int = 250,
    verbose: bool = False,
    **nmf_kwargs
    ):
    """
    Run maf.
    """
    if outdir is not ".": os.makedirs(outdir, exist_ok=True)

    # Human Genome Build
    if hg_build is not None:
        hg_build = pkg_resources.resource_filename('siganalyzer', 'ref/{}.2bit'.format(hg_build))

    # Cosmic Signatures
    if cosmic == 'cosmic2':
        cosmic = pd.read_csv(pkg_resources.resource_filename('siganalyzer', 'ref/cosmic_v2/cosmic_v2.txt'), sep='\t').dropna(1)
        cosmic_index = "Somatic Mutation Type"
    else:
        raise Exception("Not yet implemented for {}".format(cosmic))

    # Generate Spectra from Maf
    maf, spectra = get_spectra_from_maf(
        pd.read_csv(maf, sep='\t'),
        hgfile=hg_build
    )

    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'w')

    for n_iter in range(nruns):
        store['X'] = spectra

        res = ardnmf(
            spectra,
            **nmf_kwargs,
            tag="{}: ".format(n_iter)
        )

        postprocess_msigs(res, cosmic, cosmic_index)
        res["W"]["lambda"] = res["lambda"]
        store["run{}/H".format(n_iter)] = res["H"]
        store["run{}/W".format(n_iter)] = res["W"]
        store["run{}/Hraw".format(n_iter)] = res["Hraw"]
        store["run{}/Wraw".format(n_iter)] = res["Wraw"]
        store["run{}/markers".format(n_iter)] = res["markers"]
        store["run{}/signatures".format(n_iter)] = res["signatures"]
        store["run{}/log".format(n_iter)] = res["log"]
        store["run{}/cosine".format(n_iter)] = res["cosine"]

    store.close()

    df = get_nlogs_from_output(os.path.join(outdir,'nmf_output.h5'))




def run_spectra(
    spectra: Union[str, pd.DataFrame],
    outdir: str = '.',
    cosmic: str = 'cosmic2',
    nruns: int = 250,
    verbose: bool = False,
    **nmfkwargs
    ):
    """
    Run spectra.
    """
    if outdir is not ".": os.makedirs(outdir, exist_ok=True)
    pass

def run_rna(
    rna: Union[str, pd.DataFrame],
    outdir: str = '.',
    nruns: int = 20,
    verbose: bool = False,
    **nmfkwargs
    ):
    """
    Run rna.
    """
    if outdir is not ".": os.makedirs(outdir, exist_ok=True)
    pass
