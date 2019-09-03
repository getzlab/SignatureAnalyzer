from anndata import AnnData
import scanpy as sc
import numpy as np
import pandas as pd
import sys
import os
import torch

sys.path.append(os.path.join(os.path.dirname(__file__), '.', 'SignatureAnalyzer-GPU'))
from ARD_NMF import ARD_NMF, run_method_engine

# Relative Imports
from .utils import compute_phi, transfer_weights, select_signatures, select_markers

# ---------------------------------
# NMF Wrapper
# ---------------------------------
def ARD_NMF(X, K0=None, objective='poisson', max_iter=10000, del_=1, \
        tolerance=1e-6, phi=1, a=10.0, b=None, prior_on_W='L1', prior_on_H='L1', \
        report_frequency=100, active_thresh=1e-5, cut_norm=0.5, cut_diff=1.0):
        """
        AnnData wrapper for ARD-NMF.
        ------------------------
        Inputs
            * X: dataframe
                (genes x samples)

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
            * report_frequency: how often to print stats
            * parameters: parameters file
            * cut_norm
            * cut_diff: difference between mean signature and rest of signatures
                for marker selction
            * active_thresh: threshold for a latent component's impact on signature
                if the latent factor is less than this, it does not contribute
            * inplace: whether or not to edit the AnnData object directly
        Outputs
            It is highly reccomended to use only highly variable genes for factorization.
            Default parameters are to use highly variable genes and filter out mitochondrial
            or ribosomal genes.
        """
        if objective == 'poisson':
            Beta = 1
        elif objective == 'gaussian':
            Beta = 2
        else:
            ValueError("Objective should be either 'gaussian' or 'poisson'")

        if phi is None:
            phi = compute_phi(np.mean(X.values), np.var(X.values), Beta)

        data = ARD_NMF(X, objective)
        channel_names = data.channel_names
        sample_names = data.sample_names

        # ---------------------------------
        # Run NMF
        # ---------------------------------
        W, H, cost = run_method_engine(
            data, \
            a, \
            phi, \
            b, \
            Beta, \
            prior_on_W, \
            prior_on_H, \
            K0, \
            tolerance, \
            max_iter \
        )

        W,H,nsig = transfer_weights(W, H, active_thresh=active_thresh)
        sig_names = [str(i) for i in range(1,nsig+1)]
        W = pd.DataFrame(data=W, index=channel_names, columns=sig_names)
        H = pd.DataFrame(data=H, index=sig_names, columns=sample_names)

        W,H = select_signatures(W,H)
        markers, gene_signatures = select_markers(X,W,H,cut_norm=cut_norm,cut_diff=cut_diff)

        return H,W,markers,gene_signatures
