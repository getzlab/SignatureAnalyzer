import numpy as np
import pandas as pd
import sys
from tqdm import tqdm

def compute_phi(mu, var, beta):
    """
    Compute Phi
    ------------------------
    Compute the Dispersion parameter.
    """
    return var / (mu ** (2-beta))

def transfer_weights(W, H, active_thresh=1e-5):
    """
    Transfers weights from output of NMF.
    ------------------------
    Args:
        * W: input W matrix (K x n_features)
        * H: input H matrix (n_samples x K)
        * active_thresh: active threshold to consider a factor loading significant

    Returns:
        * W_final: normalized W matrix (K x n_features)
        * H_final: normalized H matrix (n_samples x K)
        * nsig: number of signatures found
    """
    nonzero_idx = (np.sum(H, axis=1) * np.sum(W, axis=0)) > active_thresh
    W_active = W[:, nonzero_idx]
    H_active = H[nonzero_idx, :]
    nsig = np.sum(nonzero_idx)

    # Normalize W and transfer weight to H matrix
    W_weight = np.sum(W_active, axis=0)
    W_final = W_active / W_weight
    H_final = W_weight[:, np.newaxis] * H_active

    return W_final, H_final, nsig

def select_signatures(W, H):
    """
    Scales NMF output by sample and feature totals to select Signatures.
    ------------------------
    Args:
        * W: input W matrix (K x n_features)
        * H: input H matrix (n_samples x K)

    Returns:
        * W: output W matrix with max_id, max, and max_norm columns
        * H: output H matrix with max_id, max, and max_norm columns
    """
    Wnorm = W.copy()
    Hnorm = H.copy()

    # Scale Matrix
    for j in range(W.shape[1]):
        Wnorm.iloc[:,j] *= H.sum(1).values[j]
        Hnorm.iloc[j,:] *= W.sum(0).values[j]

    # Normalize
    Wnorm = Wnorm.div(Wnorm.sum(1),axis=0)
    Hnorm = Hnorm.div(Hnorm.sum(0),axis=1)

    H = H.T
    Hnorm = Hnorm.T

    # Get Max Values
    H_max_id = H.idxmax(axis=1, skipna=True).astype('int')
    H['max'] = H.max(axis=1, skipna=True)
    H['max_id'] = H_max_id
    Hnorm['max_norm']=Hnorm.max(axis=1, skipna=True)

    W_max_id = W.idxmax(axis=1, skipna=True).astype('int')
    W['max'] = W.max(axis=1, skipna=True)
    W['max_id'] = W_max_id
    Wnorm['max_norm']=Wnorm.max(axis=1, skipna=True)

    H['max_norm'] = Hnorm['max_norm']
    W['max_norm'] = Wnorm['max_norm']

    H = H.rename(columns={x:'S'+x for x in list(H)[:-3]})
    W = W.rename(columns={x:'S'+x for x in list(H)[:-3]})

    return W,H

def select_markers(X, W, H, cut_norm=0.5, cut_diff=1.0, verbose=False):
    """
    Marker selection from NMF.
    ------------------------
    Args:
        * X: Input X matrix (n_samples x n_features)
        * W: input W matrix (K x n_features) with max_id, max, and max_norm columns
        * H: input H matrix (n_samples x K) with max_id, max, and max_norm columns
        * cut_norm: minimum normalized signature strength
        * cut_diff: minimum difference between selected signature and other signatures

    Returns:
        * Pandas Dataframe of NMF markers
        * Pandas Dataframe of full W matrix
    """
    markers = list()
    full = list()

    pd.options.mode.chained_assignment = None

    for n in tqdm(np.unique(W['max_id']), desc='Clusters: ', disable=not verbose):
        if H[H['max_id']==n].shape[0] > 0:
            tmp = W[W['max_id']==n]
            tmp.loc[:,'mean_on'] = X.loc[np.array(tmp.index), H[H['max_id']==n].index].mean(axis=1)
            tmp.loc[:,'mean_off'] = X.loc[np.array(tmp.index), H[H['max_id']!=n].index].mean(axis=1)
            tmp.loc[:,'diff'] = tmp.loc[:,'mean_on'] - tmp.loc[:,'mean_off']

            tmp.sort_values('diff', ascending=False, inplace=True)
            full.append(tmp)
            markers.append(tmp[(tmp['diff'] > cut_diff) & (tmp['max_norm'] >= cut_norm)])

    nmf_markers = X.loc[pd.concat(markers).index,H.max_id.sort_values().index]
    nmf_markers.index.name = 'feat'

    return nmf_markers, pd.concat(full)
