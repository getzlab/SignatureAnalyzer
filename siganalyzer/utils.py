import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity

COMPL = {"A":"T","T":"A","G":"C","C":"G"}

# ---------------------------------
# NMF Utils
# ---------------------------------
def compute_phi(mu: float, var: float, beta: float):
    """
    Compute Phi
    ------------------------
    Compute the Dispersion parameter.
    """
    return var / (mu ** (2-beta))

def transfer_weights(W: pd.DataFrame, H: pd.DataFrame, active_thresh:float = 1e-5):
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

def select_signatures(W: pd.DataFrame, H: pd.DataFrame):
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

    _rename = {x:'S'+x for x in list(H)[:-3]}
    H = H.rename(columns=_rename)
    W = W.rename(columns=_rename)

    return W,H

def select_markers(
    X: pd.DataFrame, \
    W: pd.DataFrame, \
    H: pd.DataFrame, \
    cut_norm: float = 0.5, \
    cut_diff: float = 1.0, \
    verbose: bool = False \
    ):
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

# ---------------------------------
# Mutational Signature Utils
# ---------------------------------
def compl(seq: str):
    return ''.join([COMPL[x] if x in COMPL.keys() else x for x in seq])

def _map_sigs(W: pd.DataFrame, cosmic: pd.DataFrame, sub_index:str = 'Substitution Type'):
    """
    Map signatures.
    """
    def _check_to_flip(x, ref):
        if x[2:-2] in ref:
            return x
        else:
            return compl(x)

    context_s = W['context96.word'].apply(lambda x: x[2]+'['+x[0]+'>'+x[1]+']'+x[3])
    return context_s.apply(lambda x: _check_to_flip(x,set(cosmic[sub_index])))

def postprocess_msigs(res: dict, cmap: pd.DataFrame, cosmic: pd.DataFrame, cosmic_index: str):
    """
    Post process ARD-NMF on mutational signatures.
    ------------------------
    Args:
        * res: results dictionary from ARD-NMF (see ardnmf function)
        * cmap: pandas dataframe mapping context96.num --> context96.word
        * cosmic: cosmic signatures dataframe
        * cosmic_index: feature index column in cosmic
            ** ex. in cosmic_v2, "Somatic Mutation Type" columns map to
                A[C>A]A, A[C>A]C, etc.

    Returns:
        * None, edits res dictionary directly
    """
    res["W"]["mut"] = _map_sigs(res["W"].join(cmap), cosmic)

    # Column names of NMF signatures & COSMIC References
    nmf_cols = ["S"+x for x in list(map(str, set(res["signatures"].max_id)))]
    ref_cols = list(cosmic.columns[cosmic.dtypes == 'float64'])

    # Create cosine similarity matrix
    X = res["W"].set_index("mut").join(cosmic.set_index(cosmic_index)).dropna(1).loc[:,nmf_cols+ref_cols]
    res["cosine"] = pd.DataFrame(cosine_similarity(X.T), index=X.columns, columns=X.columns).loc[ref_cols,nmf_cols]

    # Add assignments
    s_assign = dict(res["cosine"].idxmax())
    s_assign = {key:key+"-" + s_assign[key] for key in s_assign}

    # Update column names
    res["W"] = res["W"].rename(columns=s_assign)
    res["H"] = res["H"].rename(columns=s_assign)
    res["cosine"] = res["cosine"].rename(columns=s_assign)
