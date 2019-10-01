import sys
import argparse
import os
import pkg_resources
import pandas as pd
from typing import Union
import numpy as np
import matplotlib.pyplot as plt

from .utils import postprocess_msigs, get_nlogs_from_output, file_loader
from .consensus import consensus_cluster

from .plotting import k_dist, consensus_matrix
from .plotting import signature_barplot, stacked_bar
from .plotting import marker_heatmap

from .spectra import get_spectra_from_maf
from .bnmf import ardnmf

def run_maf(
    maf: Union[str, pd.DataFrame],
    outdir: str = '.',
    cosmic: str = 'cosmic2',
    hg_build: Union[str, None] = 'hg19',
    nruns: int = 10,
    verbose: bool = False,
    **nmf_kwargs
    ):
    """
    Args:
        * maf: input .maf file format
        * outdir: output directory to save files
        * cosmic: cosmic signature set to use
        * hg_build: human genome build for generating reference context
        * nruns: number of iterations for ARD-NMF
        * verbose: bool

    NMF_kwargs:
        * K0: starting number of latent components
        * objective: objective function for optimizaiton
        * max_iter: maximum number of iterations for algorithm
        * del_: n/a
        * tolerance: stop point for optimization
        * phi: dispersion parameter
        * a: shape parameter
        * b: shape parameter
        * prior_on_W: L1 or L2
        * prior_on_H: L1 or L2
        * report_freq: how often to print stats
        * active_thresh: threshold for a latent component's impact on
            signature if the latent factor is less than this, it does not contribute
        * cut_norm: min normalized value for mean signature
            (used in post-processing)
        * cut_diff: difference between mean signature and rest of signatures
            for marker selction
            (used in post-processing)
        * cuda_int: GPU to use. Defaults to 0. If "None" or if no GPU available,
            will perform decomposition using CPU.
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
            tag="\t{}/{}: ".format(n_iter,nruns-1),
            verbose=verbose,
            **nmf_kwargs
        )

        postprocess_msigs(res, cosmic, cosmic_index)
        lam = pd.DataFrame(data=res["lam"], columns=["lam"])
        lam.index.name = "K0"

        store["run{}/H".format(n_iter)] = res["H"]
        store["run{}/W".format(n_iter)] = res["W"]
        store["run{}/lam".format(n_iter)] = lam
        store["run{}/Hraw".format(n_iter)] = res["Hraw"]
        store["run{}/Wraw".format(n_iter)] = res["Wraw"]
        store["run{}/markers".format(n_iter)] = res["markers"]
        store["run{}/signatures".format(n_iter)] = res["signatures"]
        store["run{}/log".format(n_iter)] = res["log"]
        store["run{}/cosine".format(n_iter)] = res["cosine"]

    store.close()

    # Select Best Result
    aggr = get_nlogs_from_output(os.path.join(outdir,'nmf_output.h5'))
    max_k = aggr.groupby("K").size().idxmax()
    max_k_iter = aggr[aggr['K']==max_k].shape[0]
    best_run = int(aggr[aggr['K']==max_k].obj.idxmin())
    print("   * Run {} had lowest objective with mode (n={:g}) K = {:g}.".format(best_run, max_k_iter, aggr.loc[best_run]['K']))

    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'a')
    store["H"] = store["run{}/H".format(best_run)]
    store["W"] = store["run{}/W".format(best_run)]
    store["lam"] = store["run{}/lam".format(best_run)]
    store["Hraw"] = store["run{}/Hraw".format(best_run)]
    store["Wraw"] = store["run{}/Wraw".format(best_run)]
    store["markers"] = store["run{}/markers".format(best_run)]
    store["signatures"] = store["run{}/signatures".format(best_run)]
    store["log"] = store["run{}/log".format(best_run)]
    store["cosine"] = store["run{}/cosine".format(best_run)]
    store["aggr"] = aggr
    store.close()

    # Plots
    print("   * Saving report plots to {}".format(outdir))
    H = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "H")
    W = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "W")

    _ = signature_barplot(W, contributions=np.sum(H))
    plt.savefig(os.path.join(outdir, "signature_contributions.pdf"), dpi=100, bbox_inches='tight')
    _ = stacked_bar(H)
    plt.savefig(os.path.join(outdir, "signature_stacked_barplot.pdf"), dpi=100, bbox_inches='tight')
    _ = k_dist(np.array(aggr.K, dtype=int))
    plt.savefig(os.path.join(outdir, "k_dist.pdf"), dpi=100, bbox_inches='tight')

def run_spectra(
    spectra: Union[str, pd.DataFrame],
    outdir: str = '.',
    cosmic: str = 'cosmic2',
    nruns: int = 10,
    verbose: bool = False,
    **nmf_kwargs
    ):
    """
    Args:
        * spectra: filepath or pd.DataFrame of input spectra file (context x samples)
            NOTE: index should be context in the following format (1234): 3[1>2]4
        * outdir: output directory to save files
        * cosmic: cosmic signature set to use
        * nruns: number of iterations for ARD-NMF
        * verbose: bool

    NMF_kwargs:
        * K0: starting number of latent components
        * objective: objective function for optimizaiton
        * max_iter: maximum number of iterations for algorithm
        * del_: n/a
        * tolerance: stop point for optimization
        * phi: dispersion parameter
        * a: shape parameter
        * b: shape parameter
        * prior_on_W: L1 or L2
        * prior_on_H: L1 or L2
        * report_freq: how often to print stats
        * active_thresh: threshold for a latent component's impact on
            signature if the latent factor is less than this, it does not contribute
        * cut_norm: min normalized value for mean signature
            (used in post-processing)
        * cut_diff: difference between mean signature and rest of signatures
            for marker selction
            (used in post-processing)
        * cuda_int: GPU to use. Defaults to 0. If "None" or if no GPU available,
            will perform decomposition using CPU.
    """
    [nmf_kwargs.pop(key) for key in ['input', 'type', 'hg_build']]

    # Load spectra
    if isinstance(spectra, str):
        spectra = file_loader(spectra)

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
            tag="\t{}/{}: ".format(n_iter,nruns-1),
            verbose=verbose,
            **nmf_kwargs
        )

        postprocess_msigs(res, cosmic, cosmic_index)
        lam = pd.DataFrame(data=res["lam"], columns=["lam"])
        lam.index.name = "K0"

        store["run{}/H".format(n_iter)] = res["H"]
        store["run{}/W".format(n_iter)] = res["W"]
        store["run{}/lam".format(n_iter)] = lam
        store["run{}/Hraw".format(n_iter)] = res["Hraw"]
        store["run{}/Wraw".format(n_iter)] = res["Wraw"]
        store["run{}/markers".format(n_iter)] = res["markers"]
        store["run{}/signatures".format(n_iter)] = res["signatures"]
        store["run{}/log".format(n_iter)] = res["log"]
        store["run{}/cosine".format(n_iter)] = res["cosine"]

    store.close()

    # Select Best Result
    aggr = get_nlogs_from_output(os.path.join(outdir,'nmf_output.h5'))
    max_k = aggr.groupby("K").size().idxmax()
    max_k_iter = aggr[aggr['K']==max_k].shape[0]
    best_run = int(aggr[aggr['K']==max_k].obj.idxmin())
    print("   * Run {} had lowest objective with mode (n={:g}) K = {:g}.".format(best_run, max_k_iter, aggr.loc[best_run]['K']))

    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'a')
    store["H"] = store["run{}/H".format(best_run)]
    store["W"] = store["run{}/W".format(best_run)]
    store["lam"] = store["run{}/lam".format(best_run)]
    store["Hraw"] = store["run{}/Hraw".format(best_run)]
    store["Wraw"] = store["run{}/Wraw".format(best_run)]
    store["markers"] = store["run{}/markers".format(best_run)]
    store["signatures"] = store["run{}/signatures".format(best_run)]
    store["log"] = store["run{}/log".format(best_run)]
    store["cosine"] = store["run{}/cosine".format(best_run)]
    store["aggr"] = aggr
    store.close()

    # Plots
    print("   * Saving report plots to {}".format(outdir))
    H = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "H")
    W = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "W")

    _ = signature_barplot(W, contributions=np.sum(H))
    plt.savefig(os.path.join(outdir, "signature_contributions.pdf"), dpi=100, bbox_inches='tight')
    _ = stacked_bar(H)
    plt.savefig(os.path.join(outdir, "signature_stacked_barplot.pdf"), dpi=100, bbox_inches='tight')
    _ = k_dist(np.array(aggr.K, dtype=int))
    plt.savefig(os.path.join(outdir, "k_dist.pdf"), dpi=100, bbox_inches='tight')

def run_matrix(
    matrix: Union[str, pd.DataFrame],
    outdir: str = '.',
    nruns: int = 20,
    verbose: bool = False,
    **nmf_kwargs
    ):
    """
    Args:
        * matrix: expression matrix; this should be normalized to accomodate
            Gaussian noise assumption (log2-norm) (n_features x n_samples)

            NOTE: recommended to filter out lowly expressed genes for RNA:
            *************** example filtering ***************
            tpm = tpm[
                (np.sum(tpm >= 0.1, 1) > tpm.shape[1]*0.2) &
                (np.sum(counts.iloc[:,1:] >= 6, 1) > tpm.shape[1]*0.2)
            ]
            *************************************************

            NOTE: reccomended to select a set of highly variable genes following
                this (~ 2000 - 7500 genes)

        * outdir: output directory to save files
        * cosmic: cosmic signature set to use
        * nruns: number of iterations for ARD-NMF
        * verbose: bool

    NMF_kwargs:
        * K0: starting number of latent components
        * objective: objective function for optimizaiton
        * max_iter: maximum number of iterations for algorithm
        * del_: n/a
        * tolerance: stop point for optimization
        * phi: dispersion parameter
        * a: shape parameter
        * b: shape parameter
        * prior_on_W: L1 or L2
        * prior_on_H: L1 or L2
        * report_freq: how often to print stats
        * active_thresh: threshold for a latent component's impact on
            signature if the latent factor is less than this, it does not contribute
        * cut_norm: min normalized value for mean signature
            (used in post-processing)
        * cut_diff: difference between mean signature and rest of signatures
            for marker selction
            (used in post-processing)
        * cuda_int: GPU to use. Defaults to 0. If "None" or if no GPU available,
            will perform decomposition using CPU.
    """
    [nmf_kwargs.pop(key) for key in ['input', 'type', 'hg_build', 'cosmic']]

    # Load matrix
    if isinstance(matrix, str):
        matrix = file_loader(matrix)

    if outdir is not ".":
        print("   * Creating output dir at {}".format(outdir))
        os.makedirs(outdir, exist_ok=True)

    print("   * Saving ARD-NMF outputs to {}".format(os.path.join(outdir,'nmf_output.h5')))
    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'w')

    print("   * Running ARD-NMF...")
    for n_iter in range(nruns):
        store['X'] = matrix

        res = ardnmf(
            matrix,
            tag="\t{}/{}: ".format(n_iter,nruns-1),
            verbose=verbose,
            **nmf_kwargs
        )

        lam = pd.DataFrame(data=res["lam"], columns=["lam"])
        lam.index.name = "K0"

        store["run{}/H".format(n_iter)] = res["H"]
        store["run{}/W".format(n_iter)] = res["W"]
        store["run{}/lam".format(n_iter)] = lam
        store["run{}/Hraw".format(n_iter)] = res["Hraw"]
        store["run{}/Wraw".format(n_iter)] = res["Wraw"]
        store["run{}/markers".format(n_iter)] = res["markers"]
        store["run{}/signatures".format(n_iter)] = res["signatures"]
        store["run{}/log".format(n_iter)] = res["log"]

    store.close()

    # Select Best Result
    aggr = get_nlogs_from_output(os.path.join(outdir,'nmf_output.h5'))
    max_k = aggr.groupby("K").size().idxmax()
    max_k_iter = aggr[aggr['K']==max_k].shape[0]
    best_run = int(aggr[aggr['K']==max_k].obj.idxmin())
    print("   * Run {} had lowest objective with mode (n={:g}) K = {:g}.".format(best_run, max_k_iter, aggr.loc[best_run]['K']))

    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'a')
    store["H"] = store["run{}/H".format(best_run)]
    store["W"] = store["run{}/W".format(best_run)]
    store["lam"] = store["run{}/lam".format(best_run)]
    store["Hraw"] = store["run{}/Hraw".format(best_run)]
    store["Wraw"] = store["run{}/Wraw".format(best_run)]
    store["markers"] = store["run{}/markers".format(best_run)]
    store["signatures"] = store["run{}/signatures".format(best_run)]
    store["log"] = store["run{}/log".format(best_run)]
    store["aggr"] = aggr
    store.close()

    # Plots
    print("   * Saving report plots to {}".format(outdir))
    markers = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "markers")
    H = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "H")

    _ = k_dist(np.array(aggr.K, dtype=int))
    plt.savefig(os.path.join(outdir, "k_dist.pdf"), dpi=100, bbox_inches='tight')
    _ = marker_heatmap(markers, H)
    plt.savefig(os.path.join(outdir, "marker_heatmap.pdf"), dpi=100, bbox_inches='tight')

    print("   * Computing consensus matrix")
    _ = consensus_cluster(os.path.join(outdir, 'nmf_output.h5'))
    f,d = consensus_matrix(cmatrix, n_clusters=max_k_iter)
    plt.savefig(os.path.join(outdir, 'consensus_matrix.pdf'), dpi=100, bbox_inches='tight')

    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'a')
    store['consensus'] = d
    store.close()
