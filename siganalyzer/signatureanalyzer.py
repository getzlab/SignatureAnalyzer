import sys
import argparse
import os
import pkg_resources
import pandas as pd
from typing import Union

from .utils import postprocess_msigs
from .utils import get_nlogs_from_output
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
    [nmf_kwargs.pop(key) for key in ['input', 'type']]

    if outdir is not ".":
        print("   * Creating output dir at {}".format(outdir))
        os.makedirs(outdir, exist_ok=True)

    # Human Genome Build
    if hg_build is not None:
        print("   * Using {} build".format(hg_build))
        hg_build = pkg_resources.resource_filename('siganalyzer', 'ref/{}.2bit'.format(hg_build))

    # Cosmic Signatures
    if cosmic == 'cosmic2':
        print("   * Using {} signatures".format(cosmic))
        cosmic = pd.read_csv(pkg_resources.resource_filename('siganalyzer', 'ref/cosmic_v2/cosmic_v2.txt'), sep='\t').dropna(1)
        cosmic_index = "Somatic Mutation Type"
    else:
        raise Exception("Not yet implemented for {}".format(cosmic))

    # Generate Spectra from Maf
    print("   * Loading spectra from {}".format(maf))
    maf, spectra = get_spectra_from_maf(
        pd.read_csv(maf, sep='\t'),
        hgfile=hg_build
    )

    print("   * Saving ARD-NMF outputs to {}".format(os.path.join(outdir,'nmf_output.h5')))
    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'w')

    print("   * Running ARD-NMF...")
    for n_iter in range(nruns):
        store['X'] = spectra

        res = ardnmf(
            spectra,
            tag="\t{}: ".format(n_iter),
            verbose=verbose,
            **nmf_kwargs
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
    best_run = int(df.obj.idxmin())

    print("   * Run {} had the best objective function with K = {}.".format(best_run, df.loc[best_run]['K']))

    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'a')
    store["H"] = store["run{}/H".format(best_run)]
    store["W"] = store["run{}/W".format(best_run)]
    store["Hraw"] = store["run{}/Hraw".format(best_run)]
    store["Wraw"] = store["run{}/Wraw".format(best_run)]
    store["markers"] = store["run{}/markers".format(best_run)]
    store["signatures"] = store["run{}/signatures".format(best_run)]
    store["log"] = store["run{}/log".format(best_run)]
    store["cosine"] = store["run{}/cosine".format(best_run)]
    store["aggr"] = df
    store.close()

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
    [nmf_kwargs.pop(key) for key in ['input', 'type']]

    # Load spectra
    spectra = pd.read_csv(spectra, sep="\t", index_col=0)

    if outdir is not ".":
        print("   * Creating output dir at {}".format(outdir))
        os.makedirs(outdir, exist_ok=True)

    # Cosmic Signatures
    if cosmic == 'cosmic2':
        print("   * Using {} signatures".format(cosmic))
        cosmic = pd.read_csv(pkg_resources.resource_filename('siganalyzer', 'ref/cosmic_v2/cosmic_v2.txt'), sep='\t').dropna(1)
        cosmic_index = "Somatic Mutation Type"
    else:
        raise Exception("Not yet implemented for {}".format(cosmic))

    print("   * Saving ARD-NMF outputs to {}".format(os.path.join(outdir,'nmf_output.h5')))
    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'w')

    print("   * Running ARD-NMF...")
    for n_iter in range(nruns):
        store['X'] = spectra

        res = ardnmf(
            spectra,
            tag="\t{}: ".format(n_iter),
            verbose=verbose,
            **nmf_kwargs
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
    best_run = int(df.obj.idxmin())

    print("   * Run {} had the best objective function with K = {}.".format(best_run, df.loc[best_run]['K']))

    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'a')
    store["H"] = store["run{}/H".format(best_run)]
    store["W"] = store["run{}/W".format(best_run)]
    store["Hraw"] = store["run{}/Hraw".format(best_run)]
    store["Wraw"] = store["run{}/Wraw".format(best_run)]
    store["markers"] = store["run{}/markers".format(best_run)]
    store["signatures"] = store["run{}/signatures".format(best_run)]
    store["log"] = store["run{}/log".format(best_run)]
    store["cosine"] = store["run{}/cosine".format(best_run)]
    store["aggr"] = df
    store.close()

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
