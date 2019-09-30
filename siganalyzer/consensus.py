import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
from .utils import get_nruns_from_output

def consensus_cluster(filepath: str):
    """
    Consensus clustering of ard-nmf results.
    -----------------------
    Args:
        * filepath: path to the output .5 file of ARD-NMF runs

    Returns:
        * pd.DataFrame: consensus matrix from results
        * pd.Series: assignment probability for selected cluster

    """
    niter = get_nruns_from_output(filepath)
    H_selected = pd.read_hdf(filepath, "H")

    x = np.vstack([pd.read_hdf(filepath, "run{}/H".format(i)).loc[:,'max_id'].values for i in range(niter)])
    consensus_matrix = np.vstack([(x[:,[y]] == x[:]).sum(0) for y in range(x.shape[1])])

    df = pd.DataFrame(consensus_matrix, index=H_selected.index, columns=H_selected.index)
    df = df.loc[H_selected.sort_values('max_id').index, H_selected.sort_values('max_id').index]

    assign_p = pd.concat([df.loc[
        H_selected[H_selected['max_id']==x].index,
        H_selected[H_selected['max_id']==x].index
    ].mean(1) for x in set(H_selected['max_id'])])

    assign_p.name = 'assignment'
    return df, assign_p
