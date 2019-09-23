import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import Union
import numpy as np
from .utils import compl

def plot_bar(H: pd.DataFrame, figsize: tuple = (8,8)):
    """
    Plot stacked barchart & normalized stacked barchart.
    --------------------------------------
    Args:
        * H: matrix output from NMF
        * figsize: size of figure (int,int)

    Returns:
        * figure

    Example usage:
        plot_bar(H)
    """
    H = H.iloc[:,:-3].copy()
    H['sum'] = H.sum(1)
    H = H.sort_values('sum', ascending=False)

    fig,axes = plt.subplots(2,1,figsize=figsize, sharex=True)

    H.iloc[:,:-1].plot(
        kind='bar',
        stacked=True,
        ax=axes[0],
        width=1.0,
        rasterized=True
    )

    axes[0].set_xticklabels([])
    axes[0].set_xticks([])
    axes[0].set_ylabel('Counts', fontsize=20)

    H_norm = H.iloc[:,:-1].div(H['sum'].values,axis=0)
    H_norm.plot(
        kind='bar',
        stacked=True,
        ax=axes[1],
        width=1.0,
        rasterized=True
    )

    axes[1].set_xticklabels([])
    axes[1].set_xticks([])
    axes[1].set_xlabel('Samples', fontsize=16)
    axes[1].set_ylabel('Fractions', fontsize=20)
    axes[1].get_legend().remove()

    return fig

def plot_k_dist(X: np.ndarray, figsize: tuple = (8,8)):
    """
    Selected signatures plot.
    --------------------------------------
    Args:
        * X: numbers to plot
        * figsize: size of figure (int,int)

    Returns:
        * fig

    Example usage:
        plot_k_dist(np.array(pd.read_hdf("output_nmf.h5","aggr").K))

    """
    fig,ax = plt.subplots(figsize=figsize)

    sns.countplot(X, ax=ax, linewidth=2, edgecolor='k', rasterized=True)
    ax.set_ylim(0,ax.get_ylim()[1]+int(ax.get_ylim()[1]*0.1))

    ax.set_title("Aggregate of ARD-NMF (n={})".format(len(X)), fontsize=20)
    ax.set_ylabel("Number of Iterations", fontsize=18)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=18)

    return fig

def plot_signatures(W: pd.DataFrame, contributions: Union[int, pd.Series] = 1):
    """
    Plots signatures from W-matrix
    --------------------------------------
    Args:
        * W: W-matrix
        * contributions: Series of total contributions, np.sum(H), from each
            signature if W is normalized; else, 1

    Returns:
        * fig

    Example usage:
        plot_signatures(W, np.sum(H))
    """
    W = W.copy()
    sig_columns = [c for c in W if c.startswith('S')]

    if isinstance(contributions, pd.Series):
        W = W[sig_columns] * contributions[sig_columns]
    else:
        W = W[sig_columns] * contributions

    n_sigs = len(sig_columns)

    change_map = {'CA': [], 'CG': [], 'CT': [], 'TA': [], 'TC': [], 'TG': []}
    for x in W.index:
        if x.startswith('A'):
            change_map[compl(x[:2])].insert(0, x)
        else:
            change_map[x[:2]].append(x)

    color_map = {'CA': 'cyan', 'CG': 'red', 'CT': 'yellow', 'TA': 'purple', 'TC': 'green', 'TG': 'blue'}
    context_label = ['-'.join(p) for p in itertools.product('ACGT', 'ACGT')]
    x_coords = range(16)
    fig, axes = plt.subplots(nrows=n_sigs, ncols=6, figsize=(20, 2 * n_sigs), sharex='col', sharey='row')
    for row, sig in enumerate(sig_columns):
        for col, chg in enumerate(['CA', 'CG', 'CT', 'TA', 'TC', 'TG']):
            ax = axes[row, col]
            bar_heights = W[sig].loc[change_map[chg]]
            ax.bar(x_coords, bar_heights, width=.95, linewidth=1.5, edgecolor='gray', color=color_map[chg], rasterized=True)
            ax.set_xlim(-.55, 15.55)
            if row == 0:
                ax.set_title('>'.join(chg), fontsize=18)
            if row < n_sigs - 1:
                ax.tick_params(axis='x', length=0)
            else:
                ax.set_xticks(x_coords)
                ax.set_xticklabels(context_label, fontfamily='monospace', rotation='vertical')
            if col > 0:
                ax.tick_params(axis='y', length=0)
            if col == 5:
                ax.text(1.05, .5, sig, fontsize=14, rotation=270, transform=ax.transAxes, verticalalignment='center')
    plt.subplots_adjust(wspace=.08, hspace=.15)
    plt.suptitle('Mutational Signatures', fontsize=24, horizontalalignment='right')
    fig.text(.08, .5, 'Contributions', rotation='vertical', verticalalignment='center', fontsize=20, fontweight='bold')
    fig.text(.51, .03, 'Motifs', horizontalalignment='center', fontsize=20, fontweight='bold')

    return fig

def plot_marker_heatmap(markers, H, figsize=(16,12)):
    """
    Plots signatures from W-matrix
    --------------------------------------
    """
    fig, ax = plt.subplots(figsize=figsize)
    cbar_ax = fig.add_axes([.91, 0.5, .025, .3])

    sns.heatmap(markers, ax=ax, cmap="YlGnBu", rasterized=True, cbar_ax=cbar_ax)
    v,c = np.unique(H['max_id'],return_counts=True)

    ax.vlines(np.cumsum(c), *ax.get_ylim())
    ax.set_xticks(np.cumsum(c)-c/2)
    ax.set_xticklabels(v, rotation=360,fontsize=14)

    ax.set_yticks(np.arange(markers.index.values.shape[0]))
    ax.set_yticklabels(markers.index.values, fontsize=5)

    ax.set_title('')
    ax.set_xlabel('NMF Signatures', fontsize=14)
    ax.set_ylabel('Genes', fontsize=14)
    cbar_ax.set_ylabel('Normalized Expression', fontsize=12)

    return fig
