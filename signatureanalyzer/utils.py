import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
import h5py
from sklearn.metrics.pairwise import cosine_similarity
import pkg_resources
import re
import itertools
import os
import matplotlib.pyplot as plt

from sys import stdout ### GET rid of later
from .context import context_composite, context96, context1536, context78, context83, context_composite96

COMPL = {"A":"T","T":"A","G":"C","C":"G"}

# ---------------------------------
# IOUtils
# ---------------------------------
def file_loader(x):
    if x.endswith('.csv'):
        return pd.read_csv(x, index_col=0)
    elif x.endswith('.parquet'):
        return pd.read_parquet(x)
    else:
        return pd.read_csv(x, sep='\t', index_col=0)

# ---------------------------------
# NMF Utils
# ---------------------------------
def split_negatives(x: pd.DataFrame, tag: str = '_n', axis: int = 0):
    """
    Split dataframe into positive and negative components.
    --------------------
    Args:
        * x: pd.DataFrame input matrix
        * tag: string that will be added to the end of the negative variable names
        * axis: which axis to create a positive dimension for
            NOTE: this is for 2D numpy arrays

    Returns:
        * pd.DataFrame with new positive transformed matrix.
    """
    x_neg = -1 * x.copy()
    x_neg = x_neg.where(x_neg > 0, 0)

    if axis:
        x.columns = x.columns.astype(str)
        x_neg.columns = [x+'_n' for x in x.columns]
    else:
        x.index = x.index.astype(str)
        x_neg.index = [x+'_n' for x in x.index]

    return pd.concat([x.where(x > 0, 0), x_neg], axis=axis)

def l2fc(df: pd.DataFrame, center: str = 'median', axis: int = 1):
    """
    Log2 Fold-Change Input Dataframe
    -------------------------
    Args:
        * df: pd.DataFrame
        * center: center-metrix to compute log-fold change over
            ** 'median'
            ** 'mean'
        * axis: int axis to compute median and LFC across (i.e. samples)

    Returns:
        * pd.Dataframe: log2FC-transformed pandas dataframe
    """
    X = df.values

    if center == 'median':
        X_mid = np.median(X, axis)
    elif center == 'mean':
        X_mid = np.mean(X, axis)

    if axis==1:
        return pd.DataFrame(np.log2(X) - np.log2(X_mid)[:,np.newaxis], index=df.index, columns=df.columns)
    else:
        return pd.DataFrame(np.log2(X) - np.log2(X_mid)[np.newaxis], index=df.index, columns=df.columns)

def compute_phi(mu: float, var: float, beta: float):
    """
    Compute Phi
    ------------------------
    Compute the Dispersion parameter.
    """
    return var / (mu ** (2-beta))

def transfer_weights(W: pd.DataFrame, H: pd.DataFrame, active_thresh:float = 1e-2):
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
    W = W.copy()
    H = H.copy()

    # Active signatures
    nonzero_idx = (np.sum(H, axis=1) * np.sum(W, axis=0)) > active_thresh
    nsig = np.sum(nonzero_idx)

    # Raw matrices for active signatures
    W_active = W[:, nonzero_idx]
    H_active = H[nonzero_idx, :]
    
    # Normalize W and transfer weight to H matrix
    W_weight = np.sum(W_active, axis=0)
    W_final = W_active / W_weight
    H_final = W_weight[:, np.newaxis] * H_active

    return W_final, H_final, nsig, nonzero_idx

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
        Wnorm.iloc[:,j] *= H.sum(1).values[j]  # Multiple normalized signature contributions by their total mutation attribution to get total attribution per context
        Hnorm.iloc[j,:] *= W.sum(0).values[j]  # Multiply signature raw attributions by fraction of mutations per context
        
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
    Wnorm['max_norm'] = Wnorm.max(axis=1, skipna=True)
    
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
def load_reference_signatures(ref: str, verbose=True):
    """
    Load reference signatures.
    -------------------------
    Pre-processed Reference Mutational Signatures.
    """
    if ref == 'cosmic2':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/cosmic_v2/sa_cosmic2.tsv'), sep='\t').dropna(axis=1)
        reference_index = "Somatic Mutation Type"
    elif ref == 'cosmic3':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/cosmic_v3/sa_cosmic3_sbs.tsv'), sep='\t').dropna(axis=1)
        reference_index = "Somatic Mutation Type"
    elif ref == 'cosmic3_exome':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/cosmic_v3/sa_cosmic3_sbs_exome.tsv'), sep='\t').dropna(axis=1)
        reference_index = "Somatic Mutation Type"
    elif ref == 'cosmic3_DBS':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/cosmic_v3/sa_cosmic3_dbs.tsv'), sep='\t').dropna(axis=1)
        reference_index = "Somatic Mutation Type"
    elif ref == 'cosmic3_ID':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/cosmic_v3/sa_cosmic3_id.tsv'), sep='\t').dropna(axis=1)
        reference_index = "Mutation Type"
    elif ref == 'pcawg_SBS':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/PCAWG/sa_PCAWG_sbs.tsv'), sep='\t').dropna(axis=1)
        reference_index = 'Somatic Mutation Type'
    elif ref == 'pcawg_COMPOSITE':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/PCAWG/sa_PCAWG_composite.tsv'), sep='\t').dropna(axis=1)
        reference_index = 'Somatic Mutation Type'
    elif ref == 'pcawg_COMPOSITE96':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/PCAWG/sa_PCAWG_composite96.tsv'), sep='\t').dropna(axis=1)
        reference_index = 'Somatic Mutation Type'
    elif ref == 'pcawg_SBS_ID':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/PCAWG/sa_PCAWG_sbs_id.tsv'), sep='\t').dropna(axis=1)
        reference_index = 'Somatic Mutation Type'
    elif ref == 'pcawg_SBS96_ID':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/PCAWG/sa_PCAWG_sbs96_id.tsv'), sep='\t').dropna(axis=1)
        reference_index = 'Somatic Mutation Type'
    elif ref == 'polymerase_msi':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/POLE_MSI/POLE_MSI_1536SBS_ID.tsv'), sep='\t').dropna(axis=1)
        reference_index = 'Somatic Mutation Type'
    elif ref == 'polymerase_msi96':
        reference = pd.read_csv(pkg_resources.resource_filename('signatureanalyzer', 'ref/POLE_MSI/POLE_MSI_SBS96_ID.tsv'), sep='\t').dropna(axis=1)
        reference_index = 'Somatic Mutation Type'
    else:
        raise Exception("Not yet implemented for {}".format(ref))
    if verbose:
        print("   * Using {} signatures".format(ref))
    return reference, reference_index

def compl(seq: str, reverse: bool = False):
    """
    Gets the complement of a string
    Args:
        * seq: string (does not have to be base pair)
        * reverse: set to true to reverse seq

    Returns:
        * complement of seq
    """
    return ''.join([COMPL[x] if x in COMPL.keys() else x for x in (reversed(seq) if reverse else seq)])

def sbs_annotation_converter(x: str) -> str:
    """
    Eithers swaps from word -> arrow format for SBS or vice versa.
        word: (REF)(ALT)(LEFT)(RIGHT)
        arrow: (LEFT)[(REF)>(ALT)](RIGHT)
    """
    if '>' in x:
        return x[2]+x[4]+x[0]+x[6]
    else:
        return x[2]+'['+x[0]+'>'+x[1]+']'+x[3]

def sbs1536_annotation_converter(x: str) -> str:
    """
    Eithers swaps from word -> arrow format for 1536 SBS or vice versa.
        word: (REF)(ALT)(L-2)(L-1)(R+1)(R+2)
        arrow: (L-2)(L-1)[(REF)>(ALT)](R+1)(R+2)
    """
    if '>' in x:
        return x[3] + x[5] + x[:2] + x[7:9]
    else:
        return x[2:4] + '[' + x[0] + '>' + x[1] + ']' + x[4:6]
    
def _map_id_sigs(
    df: pd.DataFrame,
    ) -> pd.Series:
    """
    Map Insertion-Deletion Substitution Signatures.
    -----------------------
    Args:
        * df: pandas.core.frame.DataFrame with index to be mapped

    Returns:
        * pandas.core.series.Series with matching indices to input cosmic
    """
    def _convert_to_cosmic(x):
        i1 = 'DEL' if 'del' in x else 'INS'
        if x[0].isdigit():
            i2 = 'MH' if 'm' in x else 'repeats'
            i3 = re.search('[\d+]+', x).group()
        else:
            i2 = x[0]
            i3 = '1'
        i4 = re.search('[\d+]+$', x).group()
        if i1 == 'DEL' and i2 != 'MH':
            i4 = str(int(i4[0]) - 1) + i4[1:]
        return '_'.join([i1, i2, i3, i4])

    if df.index.name is None: df.index.name = 'index'
    df_idx = df.index.name

    context_s = df.reset_index()[df_idx]
    return context_s.apply(_convert_to_cosmic)

def _map_dbs_sigs(
    df: pd.DataFrame,
    cosmic_df: pd.DataFrame,
    sub_index: str = 'Substitution Type'
    ) -> pd.Series:
    """
    Map Doublet-Base Substitution Signatures.
    -----------------------
    Args:
        * df: pandas.core.frame.DataFrame with index to be mapped
        * cosmic_df: dataframe with Cosmic indices to map to
        * sub_index: substitution index - the column to map to in the cosmic dataframe

    Returns:
        * pandas.core.series.Series with matching indices to input cosmic
    """
    def _check_to_flip(x, ref):
        if x in ref:
            return x
        else:
            return compl(x[:2], reverse=True) + '>' + compl(x[3:], reverse=True)

    if df.index.name is None: df.index.name = 'index'
    df_idx = df.index.name

    context_s = df.reset_index()[df_idx]
    return context_s.apply(lambda x: _check_to_flip(x, set(cosmic_df[sub_index])))

def _map_sbs_sigs(
    df: pd.DataFrame,
    ref_df: pd.DataFrame,
    ref_type: str,
    sub_index: str = 'Substitution Type',
    ) -> pd.Series:
    """
    Map Single-Base Substitution Signatures.
    -----------------------
    Args:
        * df: pandas.core.frame.DataFrame with index to be mapped
        * ref_df: dataframe with reference indices to map to
        * sub_index: substitution index - the column to map to in the reference dataframe

    Returns:
        * pandas.core.series.Series with matching indices to input reference
    """
    if ref_type in ["pcawg_SBS", "pcawg_COMPOSITE", "pcawg_SBS_ID"]:
        def _check_to_flip(x, ref):
            if x[3:-3] in ref:
                return x
            else:
                return compl(x[7:], reverse=True) + '[' + compl(x[3]) + '>' + compl(x[5]) + ']' + compl(x[:2], reverse=True)
    else:
        def _check_to_flip(x, ref):
            if x[2:-2] in ref:
                return x
            else:
                return compl(x[6]) + '[' + compl(x[2]) + '>' + compl(x[4]) + ']' + compl(x[0])
            
    if df.index.name is None: df.index.name = 'index'
    df_idx = df.index.name

    if ">" not in df.index[0]:
        # Convert word format to arrow format
        if ref_type in ["pcawg_SBS","pcawg_COMPOSITE","pcawg_SBS_ID"]:
            context_s = df.reset_index()[df_idx].apply(sbs1536_annotation_converter)
        else:
            context_s = df.reset_index()[df_idx].apply(sbs_annotation_converter)
    else:
        # Already in arrow format
        context_s = df.reset_index()[df_idx]

    return context_s.apply(lambda x: _check_to_flip(x, set(ref_df[sub_index])))


def _map_composite_sigs(
    df: pd.DataFrame,
    ref_df: pd.DataFrame,
    ref_type: str,
    sub_index: str = 'Somatic Mutation Type'
    ) -> pd.Series:
    """
    Map composite signatures
    -----------------------
    Args:
        * df: pandas.core.frame.DataFrame with index to be mapped
    Returns:
        * pandas.core.series.Series with matching indices to input reference
    """
    
    if ref_type == 'pcawg_COMPOSITE':
        context_sbs_s = _map_sbs_sigs(df[df.index.isin(context1536)], ref_df.iloc[:1536], ref_type)
        context_dbs_s = _map_dbs_sigs(df[df.index.isin(context78)], ref_df.iloc[1536:1614])
    else:
        context_sbs_s = _map_sbs_sigs(df[df.index.isin(context96)], ref_df.iloc[:96], ref_type)
        context_dbs_s = _map_dbs_sigs(df[df.index.isin(context78)], ref_df.iloc[96:174])    
    
    context_id_s = df[df.index.isin(context83)].index.to_series()
    return pd.concat([context_sbs_s, context_dbs_s, context_id_s])
    

def _map_sbs_id_sigs(
    df: pd.DataFrame,
    ref_df: pd.DataFrame,
    ref_type: str,
    sub_index: str = 'Somatic Mutation Type',
    ) -> pd.Series:
    """
    Map signatures for SBS + ID spectra.
    -----------------------
    Args:
        * df: pandas.core.frame.DataFrame with index to be mapped
    Returns:
        * pandas.core.series.Series with matching indices to input reference
    """
    if ref_type == 'pcawg_SBS_ID':
        context_sbs_s = _map_sbs_sigs(df[df.index.isin(context1536)], ref_df.iloc[:1536], ref_type)
    else:
        context_sbs_s = _map_sbs_sigs(df[df.index.isin(context96)], ref_df.iloc[:96], ref_type)
    context_id_s = df[df.index.isin(context83)].index.to_series()
    return context_sbs_s.append(context_id_s)

def _map_polymerase96_id(
        df: pd.DataFrame,
        ref_df: pd.DataFrame,
        sub_index: str = 'Substitution Type',
        ) -> pd.Series:
    def _check_to_flip(x, ref):
        if x[2:-2] in ref:
            return x
        else:
            return compl(x[6]) + '[' + compl(x[2]) + '>' + compl(x[4]) + ']' + compl(x[0])
    if df.index.name is None: df.index.name = 'index'
    df_idx = df.index.name
    # Convert to arrow format. pole_msi96 starts with word format in spectra
    context_s = df.reset_index()[df_idx].map(lambda x: x if ('INS' in x or 'DEL' in x) else sbs_annotation_converter(x))
    # Reverse complement wherever necessary
    return context_s.apply(lambda x: _check_to_flip(x, set(ref_df[sub_index])) if '>' in x else x)
    
def postprocess_msigs(res: dict, ref: pd.DataFrame, ref_index: str, ref_type: str, minimum_similarity: float = 0.85):
    """
    Post process ARD-NMF on mutational signatures.
    ------------------------
    Args:
        * res: results dictionary from ARD-NMF (see ardnmf function)
        * ref: reference pd.DataFrmae
        * ref_index: feature index column in reference
            ** ex. in cosmic_v2, "Somatic Mutation Type" columns map to
                A[C>A]A, A[C>A]C, etc.
        * minimum_similarity: the minimum cosine similarity for mapping signature to reference name

    Returns:
        * None, edits res dictionary directly
    """
    # Annotate raw W matrix with mutation indices
    if ref_type in ['cosmic2','cosmic3','cosmic3_exome']:
        res["Wraw"]["mut"] = _map_sbs_sigs(res["Wraw"], ref, ref_type).values
    elif ref_type == 'cosmic3_DBS':
        res["Wraw"]["mut"] = _map_dbs_sigs(res["Wraw"], ref).values
    elif ref_type == 'cosmic3_ID':
        res["Wraw"]["mut"] = _map_id_sigs(res["Wraw"]).values  
    elif 'pcawg' in ref_type in ['pcawg_COMPOSITE', 'pcawg_COMPOSITE96', 'pcawg_SBS_ID', 'pcawg_SBS96_ID', 'pcawg_SBS']:
        # Map to PCAWG
        if ref_type in ['pcawg_COMPOSITE', 'pcawg_COMPOSITE96']:
            if ref_type == 'pcawg_COMPOSITE':
                res["Wraw"] = res["Wraw"].sort_index(key=lambda x: x.map(context_composite))
            else:
                res["Wraw"] = res["Wraw"].sort_index(key=lambda x: x.map(context_composite96))
            res["Wraw"]["mut"] = _map_composite_sigs(res["Wraw"], ref, ref_type).values
        elif ref_type == 'pcawg_SBS':
            res["Wraw"]["mut"] = _map_sbs_sigs(res["Wraw"], ref, ref_type).values
        elif ref_type in ['pcawg_SBS_ID','pcawg_SBS96_ID']:
            res["Wraw"]["mut"] = _map_sbs_id_sigs(res["Wraw"], ref, ref_type).values
            
        # load COSMIC 96 SBS and map
        cosmic_df_96, cosmic_idx_96 = load_reference_signatures("cosmic3", verbose=False)
        # Collapse 1536 to 96 if pentanucleotide context SBS
        if ref_type in ['pcawg_SBS','pcawg_COMPOSITE','pcawg_SBS_ID']: #['pcawg_SBS96_ID','pcawg_COMPOSITE96']:
            res["Wraw96"] = get96_from_1536(res["Wraw"].drop(columns=['mut'])[res["Wraw"].index.isin(context1536)])
            res["Wraw96"]["mut"] = _map_sbs_sigs(res["Wraw96"], cosmic_df_96, 'cosmic3_exome').values
        else:
            res["Wraw96"] = res["Wraw"][res["Wraw"].index.isin(context96)].copy()
            res["Wraw96"]["mut"] = _map_sbs_sigs(res["Wraw96"], cosmic_df_96, 'cosmic3_exome').values   
    elif ref_type in ['polymerase_msi', 'polymerase_msi96']:
        if ref_type == 'polymerase_msi96':
            res["Wraw"]["mut"] = _map_polymerase96_id(res["Wraw"], ref).values
        else:
            res["Wraw"]["mut"] = res["Wraw"].index
    else:
        raise Exception("Error: Invalid Reference Type (Not yet Implemented for {}".format(ref_type))
        
    # Column names of NMF signatures & References
    nmf_cols = list(res["signatures"].columns[res["signatures"].columns.str.match('S\d+')])
    ref_cols = list(ref.columns[ref.dtypes == 'float64'])
    if "pcawg" in ref_type:
        ref_cols_96 = list(cosmic_df_96.columns[cosmic_df_96.dtypes == 'float64'])
    
    # Create cosine similarity matrix
    X = res["Wraw"].set_index("mut").join(ref.set_index(ref_index)).dropna(axis=1).loc[:,nmf_cols+ref_cols]
    res["cosine"] = pd.DataFrame(cosine_similarity(X.T), index=X.columns, columns=X.columns).loc[ref_cols,nmf_cols]

    def map_sig_names(similarity_matrix,minimum_similarity):
        s_assign = dict(similarity_matrix.idxmax())
        s_max = dict(similarity_matrix.max())
        s_assign = {key:key+"-" + s_assign[key] if s_max[key] >= minimum_similarity else key+"-Unmatched" for key in s_assign}
        return(s_assign)

    # For PCAWG references, compute cosine similarity for COSMIC SBS as well
    if "pcawg" in ref_type:
        # Evaluate COSMIC cosine similarity
        X96 = res["Wraw96"].set_index("mut").join(cosmic_df_96.set_index(cosmic_idx_96)).dropna(axis=1).loc[:,nmf_cols+ref_cols_96]
        res["cosine_cosmic"] = pd.DataFrame(cosine_similarity(X96.T), index=X96.columns, columns=X96.columns).loc[ref_cols_96,nmf_cols]

        # Generate W96 matrix
        res["W96"] = res["Wraw96"].drop(columns='mut')
        res["W96"] = res["W96"].div(res["W96"].sum(0),1)
        res["W96"].columns = ['S' + str(x) for x in range(1,res["W96"].shape[1]+1)]
        
        # Get corresponding COSMIC signature name and rename
        s_assign96 = map_sig_names(res["cosine_cosmic"],minimum_similarity)
        res["cosine_cosmic"] = res["cosine_cosmic"].rename(columns=s_assign96)
        res["Wraw96"] = res["Wraw96"].rename(columns=s_assign96)
        res["W96"] = res["W96"].rename(columns=s_assign96)
        
            
    # Add assignments
    s_assign = map_sig_names(res["cosine"],minimum_similarity)

    # Update column names
    res["W"] = res["W"].rename(columns=s_assign)
    res["H"] = res["H"].rename(columns=s_assign)
    res["Wraw"] = res["Wraw"].rename(columns=s_assign)
    res["Hraw"] = res["Hraw"].T.rename(columns=s_assign)
    res["cosine"] = res["cosine"].rename(columns=s_assign)

def assign_signature_weights_to_maf(maf: pd.DataFrame, W: pd.DataFrame, H: pd.DataFrame):
    """
    Assign probabilities for each signature to each mutation in the maf
    ------------------------
    Args:
        * maf: input maf as a pd.DataFrame
        * W: W-matrix from NMF
        * H: H-matrix from NMF

    Returns:
        Maf with columns appended to the end for the weights of each signature
    """
    sig_columns = W.columns[W.columns.str.startswith('S')]
    W_prob = W.loc[maf[W.index.name], sig_columns].reset_index(drop=True)
    H_prob = H.loc[maf['sample'], sig_columns].reset_index(drop=True)
    signature_probabilities = W_prob * H_prob
    signature_probabilities /= np.sum(signature_probabilities.values, axis=1, keepdims=True)
    maf[sig_columns] = signature_probabilities
    return maf

# ---------------------------------
# Parsing Output H5 Utils
# ---------------------------------
def get_nruns_from_output(file: str) -> int:
    """
    Get number of NMF runs from an .h5 output file.
    ------------------------
    Args:
        * output file from NMF (.h5)

    Returns:
        * number of NMF runs (int)
    """
    with h5py.File(file) as f:
        return(max([int(x.split('run')[1]) for x in list(f.keys()) if 'run' in x])+1)

def get_nlogs_from_output(file: str) -> pd.DataFrame:
    """
    Returns a dataframe of the final output log from each run.
    ------------------------
    Args:
        * output file from NMF (.h5)

    Returns:
        * pd.DataFrame of reporting statistics for each run
    """
    n_runs = get_nruns_from_output(file)
    df = pd.concat([pd.read_hdf(file,"run{}/log".format(i)).reset_index().iloc[-1] for i in range(n_runs)],axis=1).T
    df.index = np.arange(n_runs)
    return df

# ---------------------------------
# Preprocessing input mafs
# ---------------------------------
def get_dnps_from_maf(maf: pd.DataFrame):
    """
    Get DNPs from a maf which has adjacent SNPs
    ________________________
    Args:
        * maf: maf DataFrame

    Returns:
        * pd.DataFrame of maf with only adjacent SNPs (annotated as DNPs)
    """
    sub_mafs = []
    for _, df in maf.loc[maf['Variant_Type'] == 'SNP'].groupby(['sample', 'Chromosome']):
        df = df.sort_values('Start_position')
        start_pos = np.array(df['Start_position'])
        pos_diff = np.diff(start_pos)
        idx = []
        if len(pos_diff) == 1 and pos_diff[0] == 1:
            idx.append(0)
        if len(pos_diff) >= 2 and pos_diff[0] == 1 and pos_diff[1] > 1:
            idx.append(0)
        idx.extend(np.flatnonzero((pos_diff[:-2] > 1) & (pos_diff[1:-1] == 1) & (pos_diff[2:] > 1)) + 1)
        if len(pos_diff) >= 2 and pos_diff[-1] == 1 and pos_diff[-2] > 1:
            idx.append(len(pos_diff) - 1)
        if idx:
            idx = np.array(idx)
            rows = df.iloc[idx][['Hugo_Symbol', 'Tumor_Sample_Barcode', 'sample',
                'Chromosome', 'Start_position', 'Reference_Allele', 'Tumor_Seq_Allele2']].reset_index(drop=True)
            rows_plus_one = df.iloc[idx + 1].reset_index()
            rows['Variant_Type'] = 'DNP'
            rows['End_position'] = rows['Start_position'] + 1
            rows['Reference_Allele'] = rows['Reference_Allele'] + rows_plus_one['Reference_Allele']
            rows['Tumor_Seq_Allele2'] = rows['Tumor_Seq_Allele2'] + rows_plus_one['Tumor_Seq_Allele2']
            sub_mafs.append(rows)
    return pd.concat(sub_mafs).reset_index(drop=True)

def get_true_snps_from_maf(maf: pd.DataFrame):
    """
    Get SNPs from a maf which has adjacent SNPs
    ________________________
    Args:
        * maf: maf DataFrame

    Returns:
        * pd.DataFrame of maf with adjacent SNPs filtered out
    """
    sub_mafs = []

    # Check if MAF includes End_position
    if not 'End_position' in list(maf):
        def get_endpos(row):
            if row['Variant_Type'] in ['SNP', 'DEL']:
                return int(row['Start_position']) + len(row['Reference_Allele']) - 1
            else:
                return int(row['Start_position']) + 1
        maf['End_position'] = maf.apply(get_endpos, axis=1)
    
    for _, df in maf.loc[maf['Variant_Type'].isin(['SNP','DEL',])].groupby(['sample', 'Chromosome']):
        df = df.sort_values('Start_position')
        start_pos = np.array(df['Start_position'])
        end_pos = np.array(df['End_position'])
        pos_diff = np.diff(start_pos)
        pos_diff_end = np.diff(end_pos)
        rem_idx = np.flatnonzero(pos_diff <= 1)
        rem_idx = np.concatenate([rem_idx, rem_idx + 1])
        rem_idx_end = np.flatnonzero(pos_diff_end <= 1)
        rem_idx_end = np.concatenate([rem_idx_end, rem_idx_end + 1])
        idx = np.delete(np.arange(len(start_pos)), np.append(rem_idx,rem_idx_end))
        if len(idx):
            rows = df.iloc[idx]
            sub_mafs.append(rows)
    temp = pd.concat(sub_mafs).reset_index(drop=True)
    return temp[temp['Variant_Type']=='SNP']


def get96_from_1536(W1536):
    """
    Convert 1536 W matrix to 96 context W matrix to extract COSMIC signatures
    ________________________
    Args:
        * W1536: 1536 context W matrix
    Returns:
        * pd.Dataframe of 96 context matrix
    """
    context96_df = pd.DataFrame(0, index=context96 , columns=W1536.columns)
    
    # Define conversion function based on format
    if ">" in W1536.index[0]:
        def convert(x):
            sbs = x[3] + x[5] + x[1] + x[7]
            if sbs[0] not in ['C','A']:
                sbs = compl(sbs[:2]) + compl(sbs[-1]) + compl(sbs[-2])
            return sbs
    else:
        def convert(x):
            sbs = x[:2] + x[3:5]
            if sbs[0] not in ['C','A']:
                sbs = compl(sbs[:2]) + compl(sbs[-1]) + compl(sbs[-2])
            return sbs
        
    # For each context in 96 SNV, sum all corresponding 1536 context rows
    for context in context96:
        context96_df.loc[context] = W1536[W1536.index.map(convert) == context].astype(np.float64).sum()
    return context96_df

def get_pole_pold_muts(maf: pd.DataFrame):
    """
    Prints sets of samples with POLE-exo mutation, POLD-exo mutation, and
    POLE-exo + POLD-exo mutation
    """
    pole_res = (268,471)
    pold_res = (304,517)
    pole = []
    pold = []
    if 'UniProt_AApos' in list(maf) or 'HGVSp_Short' in list(maf):
        if 'HGVSp_Short' in list(maf) and 'UniProt_AApos' not in list(maf):
            table = maf[maf['Hugo_Symbol'].isin(['POLE','POLD1'])][['Hugo_Symbol','Variant_Classification','Variant_Type','HGVSp_Short', 'Tumor_Sample_Barcode']].copy()
            table = table[table['Variant_Classification']=='Missense_Mutation']
            table.rename(columns={'HGVSp_Short':'UniProt_AApos'})
            get_pos = lambda x: int(''.join(c for c in x if c.isdigit()))
            table['UniProt_AApos'] = table['HGVSp_Short'].map(get_pos)
        else:
            table = maf[maf['Hugo_Symbol'].isin(['POLE','POLD1'])][['Hugo_Symbol','Variant_Classification','Variant_Type','UniProt_AApos', 'Tumor_Sample_Barcode']].copy()
            table = table[table['Variant_Classification'] == 'Missense_Mutation']
        for m in table.index:
            if table.loc[m,'Hugo_Symbol'] == 'POLE':
                if table.loc[m,'UniProt_AApos'] >= pole_res[0] and table.loc[m,'UniProt_AApos'] <= pole_res[1]:
                    pole.append(table.loc[m,'Tumor_Sample_Barcode'])
            else:
                if table.loc[m,'UniProt_AApos'] >= pole_res[0] and table.loc[m,'UniProt_AApos'] <= pole_res[1]:
                    pold.append(table.loc[m,'Tumor_Sample_Barcode'])
        stdout.write("POLE-exo* patients:\n{}\n".format(np.unique(pole)))
        stdout.write("POLD-exo* patients:\n{}\n".format(np.unique(pold)))
        stdout.write("POLE-exo* + POLD-exo* patients:\n{}\n".format(np.intersect1d(pole,pold)))
    else:
        stdout.write("Neither UniProt_AApos nor HGVSp_Short were found in the maf columns. Please try again with one of these columns")
    return np.unique(pole), np.unique(pold)

def plot_mutational_signatures(outdir, reference, k):
    from .plotting import k_dist
    from .plotting import signature_barplot, stacked_bar, signature_barplot_DBS, signature_barplot_ID, signature_barplot_composite, signature_barplot_sbs_id, signature_barplot_polymerase

    # Import plotting functions
    from .plotting import k_dist, signature_barplot, stacked_bar, signature_barplot_DBS, signature_barplot_ID, signature_barplot_composite, cosine_similarity_plot
    
    print("   * Saving report plots to {}".format(outdir))
    H = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "H")
    W = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "W")
    cosine = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "cosine")
    
    if reference == 'cosmic3_DBS':
        sys.stdout.write("Plotting Contributions Barplot:\n")
        _ = signature_barplot_DBS(W, contributions=np.sum(H))
    elif reference == 'cosmic3_ID':
        sys.stdout.write("Plotting Contributions Barplot:\n")
        _ = signature_barplot_ID(W, contributions=np.sum(H))
    elif reference == 'pcawg_SBS':
        #
        # COSMIC cosine similarity, COSMIC attribution stacked barplot
        # PCAWG contribution barplot
        #
        H96 = H.copy()
        W96 = pd.read_hdf(os.path.join(outdir, 'nmf_output.h5'), "W96")
        H96.columns = W96.columns.append(pd.Index(['max','max_id','max_norm']))
        # Plot 96 COSMIC cosine similarity
        sys.stdout.write("Plotting COSMIC Cosine Similarity:\n")
        _ = cosine_similarity_plot(pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "cosine_cosmic"))
        plt.savefig(os.path.join(outdir, "cosine_similarity_plot_96.pdf"), dpi=100, bbox_inches='tight')
        # Plot COSMIC Signature Attribution Stacked Barplot
        sys.stdout.write("Plotting COSMIC Attribution Barplot:\n")
        _ = stacked_bar(H96, 'cosmic3')
        plt.savefig(os.path.join(outdir, "signature_stacked_barplot_cosmic.pdf"), dpi=100, bbox_inches='tight')
        # Plot signature contribution barplot collapsed to 96 SBS
        sys.stdout.write("Plotting {} Contributions Barplot:\n".format(reference))
        _ = signature_barplot(W96, contributions=np.sum(H96))
    elif reference in ['pcawg_COMPOSITE','pcawg_COMPOSITE96']:
        #
        # COSMIC cosine similarity, COSMIC attribution stacked barplot, COSMIC contribution barplot
        # PCAWG contribution barplot
        #
        H96 = H.copy()
        W96 = pd.read_hdf(os.path.join(outdir, 'nmf_output.h5'), "W96")
        H96.columns = W96.columns.append(pd.Index(['max','max_id','max_norm']))        # COSMIC cosine similarity
        sys.stdout.write("Plotting COSMIC Cosine Similarity:\n")
        _ = cosine_similarity_plot(pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "cosine_cosmic"))
        plt.savefig(os.path.join(outdir, "cosine_similarity_plot_96.pdf"), dpi=100, bbox_inches='tight')
        # COSMIC attribution stacked barplot
        sys.stdout.write("Plotting COSMIC Attribution Barplot:\n")
        _ = stacked_bar(H96, 'cosmic3')
        plt.savefig(os.path.join(outdir,'signature_stacked_barplot_cosmic.pdf'), dpi=100, bbox_inches='tight')
        # COSMIC contribution barplot
        sys.stdout.write("Plotting COSMIC Contribution Barplot:\n")
        _ = signature_barplot(W96, contributions=np.sum(H96))
        plt.savefig(os.path.join(outdir, "signature_contributions_COSMIC.pdf"), dpi=100,bbox_inches='tight')
        # PCAWG contribution barplot
        if reference == 'pcawg_COMPOSITE':
            W_plot = pd.concat([get96_from_1536(W[W.index.isin(context1536)]),W[~W.index.isin(context1536)]])
        else:
            W_plot = W
        sys.stdout.write("Plotting {} Contributions Barplot:\n".format(reference))
        _ = signature_barplot_composite(W_plot, contributions=np.sum(H))
    elif reference in ['pcawg_SBS_ID', 'pcawg_SBS96_ID']:
        #
        # COSMIC cosine similarity, COSMIC attribution stacked barplot, COSMIC contribution barplot
        # PCAWG contribution barplot
        #
        H96 = H.copy()
        W96 = pd.read_hdf(os.path.join(outdir, 'nmf_output.h5'), "W96")
        H96.columns = W96.columns.append(pd.Index(['max','max_id','max_norm']))
        # COSMIC cosine similarity
        sys.stdout.write("Plotting COSMIC Cosine Similarity :\n")
        _ = cosine_similarity_plot(pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "cosine_cosmic"))
        plt.savefig(os.path.join(outdir, "cosine_similarity_plot_96.pdf"), dpi=100, bbox_inches='tight')
        # COSMIC attribution stacked barplot
        sys.stdout.write("Plotting COSMIC Attributions Barplot:\n")
        _ = stacked_bar(H96, 'cosmic3')
        plt.savefig(os.path.join(outdir,'signature_stacked_barplot_cosmic.pdf'), dpi=100, bbox_inches='tight')
        # COSMIC contribution barplot
        sys.stdout.write("Plotting COSMIC Contribution Barplot:\n")
        _ = signature_barplot(W96, contributions=np.sum(H96))
        plt.savefig(os.path.join(outdir, "signature_contributions_COSMIC.pdf"), dpi=100,bbox_inches='tight')
        # PCAWG contribution barplot
        if reference == 'pcawg_SBS_ID':
            W_plot = pd.concat([get96_from_1536(W[W.index.isin(context1536)]),W[~W.index.isin(context1536)]])
        else:
            W_plot = W
        sys.stdout.write("Plotting {} Contributions Barplot:\n".format(reference))
        _ = signature_barplot_sbs_id(W_plot, contributions=np.sum(H))
    elif reference in ['polymerase_msi', 'polymerase_msi96']:
        sys.stdout.write("Plotting Contributions Barplot:\n")
        W_plot = W
        if reference == 'polymerase_msi':
            W_plot = pd.concat([get96_from_1536(W[W.index.isin(context1536)]),W[~W.index.isin(context1536)]])
        _ = signature_barplot_polymerase(W_plot, contributions=np.sum(H))
    else:
        _ = signature_barplot(W, contributions=np.sum(H))
        
    # Plot signature contributions, attribution stacked barplot, K distribution, and cosine similarity
    plt.savefig(os.path.join(outdir, "signature_contributions.pdf"), dpi=100, bbox_inches='tight')
    sys.stdout.write("Plotting {} Attributions Barplot:\n".format(reference))
    _ = stacked_bar(H,reference)
    plt.savefig(os.path.join(outdir, "signature_stacked_barplot.pdf"), dpi=100, bbox_inches='tight')
    sys.stdout.write("Plotting K Histogram:\n")
    _ = k_dist(np.array(k, dtype=int))
    plt.savefig(os.path.join(outdir, "k_dist.pdf"), dpi=100, bbox_inches='tight')
    sys.stdout.write("Plotting {} Cosine Similarity:\n".format(reference))
    _ = cosine_similarity_plot(cosine)
    plt.savefig(os.path.join(outdir, "cosine_similarity_plot.pdf"), dpi=100, bbox_inches='tight')
    
