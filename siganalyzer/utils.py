import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
import h5py
from sklearn.metrics.pairwise import cosine_similarity
import pkg_resources

COMPL = {"A":"T","T":"A","G":"C","C":"G"}

# ---------------------------------
# IOUtils
# ---------------------------------
def file_loader(x):
    if x.endswith('.tsv') or x.endswith('.txt'):
        return pd.read_csv(x, sep='\t', index_col=0)
    elif x.endswith('.parquet'):
        return pd.read_parquet(x)
    else:
        return pd.read_csv(x, index_col=0)

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

    nonzero_idx = (np.sum(H, axis=1) * np.sum(W, axis=0)) > active_thresh

    W_active = W[:, nonzero_idx]
    H_active = H[nonzero_idx, :]
    nsig = np.sum(nonzero_idx)

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
def load_cosmic_signatures(cosmic: str):
    """
    Load cosmic signatures.
    -------------------------
    Pre-processed Cosmic Mutational Signatures.
    """
    if cosmic == 'cosmic2':
        print("   * Using {} signatures".format(cosmic))
        cosmic = pd.read_csv(pkg_resources.resource_filename('siganalyzer', 'ref/cosmic_v2/sa_cosmic2.tsv'), sep='\t').dropna(1)
        cosmic_index = "Somatic Mutation Type"
    elif cosmic == 'cosmic3':
        print("   * Using {} signatures".format(cosmic))
        cosmic = pd.read_csv(pkg_resources.resource_filename('siganalyzer', 'ref/cosmic_v3/sa_cosmic3_sbs.tsv'), sep='\t').dropna(1)
        cosmic_index = "Somatic Mutation Type"
    elif cosmic == 'cosmic3_exome':
        print("   * Using {} signatures".format(cosmic))
        cosmic = pd.read_csv(pkg_resources.resource_filename('siganalyzer', 'ref/cosmic_v3/sa_cosmic3_sbs_exome.tsv'), sep='\t').dropna(1)
        cosmic_index = "Somatic Mutation Type"
    else:
        raise Exception("Not yet implemented for {}".format(cosmic))

    return cosmic, cosmic_index

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

def postprocess_msigs(res: dict, cosmic: pd.DataFrame, cosmic_index: str):
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
    res["Wraw"]["mut"] = _map_sigs(res["Wraw"].reset_index(), cosmic).values

    # Column names of NMF signatures & COSMIC References
    nmf_cols = ["S"+x for x in list(map(str, set(res["signatures"].max_id)))]
    ref_cols = list(cosmic.columns[cosmic.dtypes == 'float64'])

    # Create cosine similarity matrix
    X = res["Wraw"].set_index("mut").join(cosmic.set_index(cosmic_index)).dropna(1).loc[:,nmf_cols+ref_cols]
    res["cosine"] = pd.DataFrame(cosine_similarity(X.T), index=X.columns, columns=X.columns).loc[ref_cols,nmf_cols]

    # Add assignments
    s_assign = dict(res["cosine"].idxmax())
    s_assign = {key:key+"-" + s_assign[key] for key in s_assign}

    # Update column names
    res["W"] = res["W"].rename(columns=s_assign)
    res["H"] = res["H"].rename(columns=s_assign)
    res["Wraw"] = res["Wraw"].rename(columns=s_assign)
    res["Hraw"] = res["Hraw"].T.rename(columns=s_assign)
    res["cosine"] = res["cosine"].rename(columns=s_assign)

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
    df = pd.concat([pd.read_hdf(file,"run{}/log".format(i)).reset_index().iloc[-1] for i in range(n_runs)],1).T
    df.index = np.arange(n_runs)
    return df

def get_dnps_from_maf(maf: pd.DataFrame):
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
    sub_mafs = []
    for _, df in maf.loc[maf['Variant_Type'] == 'SNP'].groupby(['sample', 'Chromosome']):
        df = df.sort_values('Start_position')
        start_pos = np.array(df['Start_position'])
        pos_diff = np.diff(start_pos)
        rem_idx = np.flatnonzero(pos_diff <= 1)
        rem_idx = np.concatenate([rem_idx, rem_idx + 1])
        idx = np.delete(np.arange(len(start_pos)), rem_idx)
        if len(idx):
            rows = df.iloc[idx]
            sub_mafs.append(rows)
    return pd.concat(sub_mafs).reset_index(drop=True)
