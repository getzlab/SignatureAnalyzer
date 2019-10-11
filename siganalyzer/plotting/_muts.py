import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import Union
import numpy as np

from ..utils import compl
from ..spectra import context96, context78

def stacked_bar(H: pd.DataFrame, figsize: tuple = (8,8)):
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
    axes[1].set_ylim([0,1])

    return fig

def signature_barplot(W: pd.DataFrame, contributions: Union[int, pd.Series] = 1):
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
        signature_barplot(W, np.sum(H))
    """
    W = W.copy()
    for c in context96:
        if c not in W.index:
            W.loc[c] = 0
    W.sort_index(inplace=True)
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
    fig, axes = plt.subplots(nrows=n_sigs, ncols=6, figsize=(20, 2.5 * n_sigs), sharex='col', sharey='row')
    for row, sig in enumerate(sig_columns):
        for col, chg in enumerate(['CA', 'CG', 'CT', 'TA', 'TC', 'TG']):
            if n_sigs == 1:
                ax = axes[col]
            else:
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
    plt.suptitle('Mutational Signatures', y=1.18 - n_sigs * .06, fontsize=24, horizontalalignment='right')
    fig.text(.08, .5, 'Contributions', rotation='vertical', verticalalignment='center', fontsize=20, fontweight='bold')
    fig.text(.51, -.15 + n_sigs * .05, 'Motifs', horizontalalignment='center', fontsize=20, fontweight='bold')

    return fig

def signature_barplot_DBS(W, contributions):
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
        signature_barplot_DBS(W, np.sum(H))
    """
    W = W.copy()
    for c in context78:
        if c not in W.index:
            W.loc[c] = 0
    W.sort_index(inplace=True)
    sig_columns = [c for c in W if c.startswith('S')]
    if isinstance(contributions, pd.Series):
        W = W[sig_columns] * contributions[sig_columns]
    else:
        W = W[sig_columns] * contributions

    n_sigs = len(sig_columns)

    ref_map = {'AC': [], 'AT': [], 'CC': [], 'CG': [], 'CT': [], 'GC': [], 'TA': [], 'TC': [], 'TG': [], 'TT': []}
    for x in W.index:
        ref_map[x[:2]].append(x)
    x_coords = {ref: range(len(sigs)) for ref, sigs in ref_map.items()}

    color_map = {'AC': '#99CCFF', 'AT': '#0000FF', 'CC': '#CCFF99', 'CG': '#00FF00', 'CT': '#FF99CC',
                 'GC': '#FF0000', 'TA': '#FFCC99', 'TC': '#FF8000', 'TG': '#CC99FF', 'TT': '#8000FF'}
    fig, axes = plt.subplots(nrows=n_sigs, ncols=10, figsize=(20, 2.5 * n_sigs), sharex='col',
                             sharey='row', gridspec_kw={'width_ratios': (3, 2, 3, 2, 3, 2, 2, 3, 3, 3)})
    for row, sig in enumerate(sig_columns):
        for col, ref in enumerate(ref_map):
            if n_sigs == 1:
                ax = axes[col]
            else:
                ax = axes[row, col]
            bar_heights = W[sig].loc[ref_map[ref]]
            ax.bar(x_coords[ref], bar_heights, width=.95, linewidth=1.5, edgecolor='gray', color=color_map[ref],
                   rasterized=True)
            ax.set_xlim(-.55, x_coords[ref][-1] + .55)
            if row == 0:
                ax.set_title(ref)
            if row < n_sigs - 1:
                ax.tick_params(axis='x', length=0)
            else:
                xlabels = [x[3:] for x in ref_map[ref]]
                ax.set_xticks(x_coords[ref])
                ax.set_xticklabels(xlabels, fontfamily='monospace', rotation='vertical')
            if col > 0:
                ax.tick_params(axis='y', length=0)
            if col == 9:
                ax.text(1.05, .5, sig, fontsize=14, rotation=270, transform=ax.transAxes, verticalalignment='center')

    plt.subplots_adjust(wspace=.08, hspace=.15)
    plt.suptitle('Mutational Signatures', y=1.18 - n_sigs * .06, fontsize=24, horizontalalignment='right')
    fig.text(.08, .5, 'Contributions', rotation='vertical', verticalalignment='center', fontsize=20, fontweight='bold')
    fig.text(.51, -.15 + n_sigs * .05, 'Motifs', horizontalalignment='center', fontsize=20, fontweight='bold')

    return fig
