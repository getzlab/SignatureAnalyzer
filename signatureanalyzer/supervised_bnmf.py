from sys import stdout
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from typing import Union
SEloss = nn.MSELoss(reduction = 'sum')

from .utils import select_signatures, select_markers, transfer_weights

def supervised_ardnmf(
    X_i: pd.DataFrame,
    W_ref_i: pd.DataFrame,
    lam_ref: Union[None,np.ndarray] = None,
    update_lam: bool = False,
    objective: str = 'poisson',
    max_iter: int = 10000,
    del_: int = 1,
    tolerance: float = 1e-6,
    phi: float = 1.0,
    a: float = 10.0,
    b: Union[float, None] = None ,
    prior_on_W: str = 'L1',
    prior_on_H: str = 'L1',
    report_freq: int = 100,
    active_thresh: float = -1.0,
    cut_norm: float = 0.5,
    cut_diff: float = 1.0,
    cuda_int: Union[int, None] = 0,
    verbose: bool = True,
    tag: str = ""
    ) -> dict:
    """
    Supervised ARD-NMF
    ------------------------------------------------------------------------
    Args:
        * X: input matrix (features x samples)
        * W_ref: reference dataframe
        * lam_ref: reference Lambda
        * update_lam: update Lambda during optimization
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
        * verbose: verbose reporting of algorithm convergence
    Returns:
        * results - dict with following keys:
            {'H', 'W', 'Wraw', 'Hraw', 'markers', 'signatures', 'objective', 'log', 'lam'}
    """
    assert objective in ('poisson','gaussian'), \
        "Unable to use {}; specify either poisson or gaussian objective.".format(objective)

    if objective == 'poisson': Beta = 1
    if objective == 'gaussian': Beta = 2

    assert prior_on_W in ('L1','L2'), \
        "Unable to use {}; use either L1 or L2 prior on W.".format(prior_on_W)
    assert prior_on_H in ('L1','L2'), \
        "Unable to use {}; use either L1 or L2 prior on H.".format(prior_on_H)

    # ---------------------------------
    # Intersecting features
    # ---------------------------------
    intersecting_genes = np.intersect1d(X_i.index,W_ref_i.index)
    print("{} intersecting features found.".format(intersecting_genes.shape[0]))

    X = X_i.loc[intersecting_genes]
    W_ref = W_ref_i.loc[intersecting_genes]

    # Fix naming
    sig_map = {"S"+str(i):x for i,x in enumerate(W_ref.columns)}
    sig_num_map = {i:x for i,x in enumerate(W_ref.columns)}
    inv_sig_map = {x:"S"+str(i) for i,x in enumerate(W_ref.columns)}
    W_ref.columns = W_ref.columns.map(inv_sig_map)

    # ---------------------------------
    # Load data into tensors
    # ---------------------------------
    results = ARD_NMF(X, objective, verbose=verbose)
    channel_names = results.channel_names
    sample_names = results.sample_names
    sig_names = W_ref.columns
    sig_names = [x.split("S")[-1] for x in sig_names]

    # initalize the NMF run
    results.initalize_data(a, phi, b, prior_on_W, prior_on_H, Beta, W_ref.shape[1])
    W_ref = torch.tensor(W_ref.loc[results.channel_names].values, dtype=results.dtype, requires_grad=False)
    results.W = W_ref

    # specify GPU
    cuda_string = 'cuda:'+str(cuda_int)

    # copy data to GPU
    if torch.cuda.device_count() > 0 and cuda_int is not None:
        print("   * Using GPU: {}".format(cuda_string))
        W, H, V, Lambda, C, b0, eps_, phi = results.W.cuda(cuda_string),results.H.cuda(cuda_string),results.V.cuda(cuda_string),results.Lambda.cuda(cuda_string),results.C.cuda(cuda_string),results.b.cuda(cuda_string),results.eps_.cuda(cuda_string),results.phi.cuda(cuda_string)
    else:
        W, H, V, Lambda, C, b0, eps_, phi = results.W,results.H,results.V,results.Lambda,results.C,results.b,results.eps_,results.phi
        print("   * Using CPU")

    # tracking variables
    deltrack = 1000
    times = list()
    report = dict()
    iter = 0

    if lam_ref is not None:
        Lambda = torch.squeeze(torch.tensor(lam_ref))
        if torch.cuda.device_count() > 0 and cuda_int is not None:
            Lambda = Lambda.cuda(cuda_string)

    H_previous = H

    # set method
    method = SS_NMF_algorithim(Beta, prior_on_H, prior_on_W, update_lam=update_lam)

    while deltrack >= tolerance and iter < max_iter:
        # compute updates
        H,W,Lambda = method.forward(W,H,V,Lambda,C,b0,eps_,phi)

        # compute objective and cost
        l_ = beta_div(Beta,V,W,H,eps_)
        cost_ = calculate_objective_function(Beta,V,W,H,Lambda,C,eps_,phi,results.K0)

        # update tracking
        deltrack = torch.max(torch.div(torch.abs(H-H_previous), H_previous+1e-30))
        H_previous = H

        # ---------------------------- Reporting ---------------------------- #
        if iter % report_freq == 0:
            report[iter] = {
                'K': torch.sum((torch.sum(H,1) * torch.sum(W,0))>active_thresh).cpu().numpy(),
                'obj': cost_.cpu().numpy(),
                'b_div': l_.cpu().numpy(),
                'lam': torch.sum(Lambda).cpu().numpy(),
                'del': deltrack.cpu().numpy(),
                'W_sum': torch.sum(W).cpu().numpy(),
                'H_sum': torch.sum(H).cpu().numpy()
            }
            print_report(iter,report,verbose,tag)
        # ------------------------------------------------------------------- #
        iter+=1

    # --------------------------- Final Report --------------------------- #
    report[iter] = {
        'K': torch.sum((torch.sum(H,1) * torch.sum(W,0))>active_thresh).cpu().numpy(),
        'obj': cost_.cpu().numpy(),
        'b_div': l_.cpu().numpy(),
        'lam': torch.sum(Lambda).cpu().numpy(),
        'del': deltrack.cpu().numpy(),
        'W_sum': torch.sum(W).cpu().numpy(),
        'H_sum': torch.sum(H).cpu().numpy()
    }
    print_report(iter, report, verbose, tag)

    if not verbose:
        stdout.write("\n")

    final_report = pd.DataFrame.from_dict(report).T
    final_report.index.name = 'iter'

    cost = cost_.cpu().numpy()
    Lambda = Lambda.cpu().numpy()
    W = W.cpu().numpy()
    H = H.cpu().numpy()

    Hraw = pd.DataFrame(data = H, index=sig_names, columns=sample_names)
    Wraw = pd.DataFrame(data = W, index=channel_names, columns=sig_names)

    # Transfer weights
    W, H, nsig, nonzero_idx = transfer_weights(W, H, active_thresh=active_thresh)
    H = pd.DataFrame(data = H, index=sig_names, columns=sample_names)
    W = pd.DataFrame(data = W, index=channel_names, columns=sig_names)

    # Fix log typing
    final_report['K'] = final_report['K'].astype(int)
    final_report['obj'] = final_report['obj'].astype('float')
    final_report['b_div'] = final_report['b_div'].astype('float')
    final_report['lam'] = final_report['lam'].astype('float')
    final_report['del'] = final_report['del'].astype('float')
    final_report['W_sum'] = final_report['W_sum'].astype('float')
    final_report['H_sum'] = final_report['H_sum'].astype('float')

    W,H = select_signatures(W,H)
    markers, signatures = select_markers(X, W, H, cut_norm=cut_norm, cut_diff=cut_diff, verbose=verbose)

    # Assign Signatures to Original Input Names
    W = W.rename(columns=sig_map)
    W['max_id'] = W['max_id'].apply(lambda x: sig_num_map[x])
    H = H.rename(columns=sig_map)
    H['max_id'] = H['max_id'].apply(lambda x: sig_num_map[x])

    Hraw.index = Hraw.index.astype(int)
    Hraw = Hraw.rename(index=sig_num_map)
    Wraw.columns = Wraw.columns.astype(int)
    Wraw = Wraw.rename(columns=sig_num_map)

    return {
        'H': H,
        'W': W,
        'Wraw':Wraw,
        'Hraw':Hraw,
        'markers': markers,
        'signatures': signatures,
        'objective': cost,
        'log': final_report,
        'lam': Lambda
    }

from .signatureanalyzer_gpu.ARD_NMF import ARD_NMF
from .signatureanalyzer_gpu.ARD_NMF import print_report
from .signatureanalyzer_gpu.NMF_functions import beta_div
from .signatureanalyzer_gpu.NMF_functions import calculate_objective_function
from .signatureanalyzer_gpu.NMF_functions import update_H_poisson_L1
from .signatureanalyzer_gpu.NMF_functions import update_H_poisson_L2
from .signatureanalyzer_gpu.NMF_functions import update_H_gaussian_L1
from .signatureanalyzer_gpu.NMF_functions import update_H_gaussian_L2
from .signatureanalyzer_gpu.NMF_functions import update_W_poisson_L1
from .signatureanalyzer_gpu.NMF_functions import update_W_poisson_L2
from .signatureanalyzer_gpu.NMF_functions import update_W_gaussian_L1
from .signatureanalyzer_gpu.NMF_functions import update_W_gaussian_L2
from .signatureanalyzer_gpu.NMF_functions import update_del
from .signatureanalyzer_gpu.NMF_functions import update_lambda_L1
from .signatureanalyzer_gpu.NMF_functions import update_lambda_L2
from .signatureanalyzer_gpu.NMF_functions import update_lambda_L1_L2
from .signatureanalyzer_gpu.NMF_functions import update_lambda_L2_L1

class SS_NMF_algorithim(nn.Module):
    ''' implements ARD NMF from https://arxiv.org/pdf/1111.6085.pdf '''
    def __init__(self, Beta, H_prior, W_prior, update_lam=False):
        super(SS_NMF_algorithim, self).__init__()
        # Beta paramaterizes the objective function
        # Beta = 1 induces a poisson objective
        # Beta = 2 induces a gaussian objective
        # Priors on the component matrices are Exponential (L1) and half-normal (L2)

        self.update_lam = update_lam

        if Beta == 1 and H_prior == 'L1' and W_prior == 'L1' :
            self.update_W = update_W_poisson_L1
            self.update_H = update_H_poisson_L1
            self.lambda_update = update_lambda_L1

        elif Beta == 1 and H_prior == 'L1' and W_prior == 'L2':
            self.update_W = update_W_poisson_L2
            self.update_H = update_H_poisson_L1
            self.lambda_update = update_lambda_L2_L1

        elif Beta == 1 and H_prior == 'L2' and W_prior == 'L1':
            self.update_W = update_W_poisson_L1
            self.update_H = update_H_poisson_L2
            self.lambda_update = update_lambda_L1_L2

        elif Beta == 1 and H_prior == 'L2' and W_prior == 'L2':
            self.update_W = update_W_poisson_L2
            self.update_H = update_H_poisson_L2
            self.lambda_update = update_lambda_L2

        if Beta == 2 and H_prior == 'L1' and W_prior == 'L1':
            self.update_W = update_W_gaussian_L1
            self.update_H = update_H_gaussian_L1
            self.lambda_update = update_lambda_L1

        elif Beta == 2 and H_prior == 'L1' and W_prior == 'L2':
            self.update_W = update_W_gaussian_L2
            self.update_H = update_H_gaussian_L1
            self.lambda_update = update_lambda_L2_L1

        elif Beta == 2 and H_prior == 'L2' and W_prior == 'L1':
            self.update_W = update_W_gaussian_L1
            self.update_H = update_H_gaussian_L2
            self.lambda_update = update_lambda_L1_L2

        elif Beta == 2 and H_prior == 'L2' and W_prior == 'L2':
            self.update_W = update_W_gaussian_L2
            self.update_H = update_H_gaussian_L2
            self.lambda_update = update_lambda_L2

    def forward(self, W, H, V, lambda_, C, b0, eps_, phi):
        h_ = self.update_H(H, W, lambda_, phi, V, eps_)

        if self.update_lam:
            lambda_ = self.lambda_update(W, h_, b0, C, eps_)

        return h_, W, lambda_
