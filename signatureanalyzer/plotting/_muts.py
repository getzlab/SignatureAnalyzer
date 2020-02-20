import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import Union
import numpy as np
import re

from ..utils import compl, sbs_annotation_converter
from ..spectra import context96, context78, context83

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

def _map_sbs_sigs_back(df: pd.DataFrame) -> pd.Series:
    """
    Map Back Single-Base Substitution Signatures.
    -----------------------
    Args:
        * df: pandas.core.frame.DataFrame with index to be mapped

    Returns:
        * pandas.core.series.Series with matching indices to context96
    """
    def _check_to_flip(x, ref):
        if x in ref:
            return x
        else:
            return compl(x)

    if df.index.name is None: df.index.name = 'index'
    df_idx = df.index.name

    if ">" in df.index[0]:
        # Already in arrow format
        context_s = df.reset_index()[df_idx].apply(sbs_annotation_converter)
    else:
        # Already in word format
        context_s = df.reset_index()[df_idx]

    return context_s.apply(lambda x: _check_to_flip(x, context96.keys()))

def _map_id_sigs_back(df: pd.DataFrame) -> pd.Series:
    """
        Map Back Insertion-Deletion Signatures.
        -----------------------
        Args:
            * df: pandas.core.frame.DataFrame with index to be mapped

        Returns:
            * pandas.core.series.Series with matching indices to context83
        """
    if df.index.name is None: df.index.name = 'index'
    df_idx = df.index.name

    context_s = df.reset_index()[df_idx]

    def _convert_from_cosmic(x):
        if x in context83:
            return x
        i1, i2, i3, i4 = x.split('_')
        pre = i2 if i3 == '1' else i3
        main = i1.lower() + ('m' if i2 == 'MH' else '')
        if main == 'del':
            post = str(int(i4[0]) + 1) + i4[1:]
        else:
            post = i4
        return pre + main + post

    return context_s.apply(_convert_from_cosmic)

def signature_barplot(W: pd.DataFrame, contributions: Union[int, pd.Series] = 1):
    """
    Plots signatures from W-matrix for Single-Base Substitutions
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
    W.index = _map_sbs_sigs_back(W)

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

    context_label = []
    change_map = {'CA': [], 'CG': [], 'CT': [], 'TA': [], 'TC': [], 'TG': []}
    for p in itertools.product('ACGT', 'ACGT'):
        context = ''.join(p)
        compl_context = compl(context, reverse=True)
        context_label.append('-'.join(context))
        for key in change_map:
            if key.startswith('C'):
                change_map[key].append(key + context)
            else:
                change_map[key].append(compl(key) + compl_context)
    color_map = {'CA': 'cyan', 'CG': 'red', 'CT': 'yellow', 'TA': 'purple', 'TC': 'green', 'TG': 'blue'}

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
                if col == 0:
                    ax.text(51.2 / 16, 1.3, 'Mutational Signatures', transform=ax.transAxes,
                            horizontalalignment='center', fontsize=24)
            if row < n_sigs - 1:
                ax.tick_params(axis='x', length=0)
            else:
                ax.set_xticks(x_coords)
                ax.set_xticklabels(context_label, fontfamily='monospace', rotation='vertical')
                if col == 0:
                    ax.text(51.2 / 16, -.4, 'Motifs', transform=ax.transAxes, horizontalalignment='center', fontsize=20,
                            fontweight='bold')
            if col > 0:
                ax.tick_params(axis='y', length=0)
            if col == 5:
                ax.text(1.05, .5, sig, fontsize=14, rotation=270, transform=ax.transAxes, verticalalignment='center')

    plt.subplots_adjust(wspace=.08, hspace=.15)
    fig.text(.08, .5, 'Contributions', rotation='vertical', verticalalignment='center', fontsize=20, fontweight='bold')

    return fig

def signature_barplot_DBS(W, contributions):
    """
    Plots signatures from W-matrix for Doublet-Base Substitutions
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
                if col == 0:
                    ax.text(44.5 / 6, 1.2, 'Mutational Signatures', transform=ax.transAxes,
                            horizontalalignment='center', fontsize=24)
            if row < n_sigs - 1:
                ax.tick_params(axis='x', length=0)
            else:
                xlabels = [x[3:] for x in ref_map[ref]]
                ax.set_xticks(x_coords[ref])
                ax.set_xticklabels(xlabels, fontfamily='monospace', rotation='vertical')
                if col == 0:
                    ax.text(44.5 / 6, -.3, 'Motifs', transform=ax.transAxes, horizontalalignment='center', fontsize=20,
                            fontweight='bold')
            if col > 0:
                ax.tick_params(axis='y', length=0)
            if col == 9:
                ax.text(1.05, .5, sig, fontsize=14, rotation=270, transform=ax.transAxes, verticalalignment='center')

    plt.subplots_adjust(wspace=.08, hspace=.15)
    fig.text(.08, .5, 'Contributions', rotation='vertical', verticalalignment='center', fontsize=20, fontweight='bold')

    return fig

def signature_barplot_ID(W, contributions):
    """
    Plots signatures from W-matrix for Insertions-Deletions
    --------------------------------------
    Args:
        * W: W-matrix
        * contributions: Series of total contributions, np.sum(H), from each
            signature if W is normalized; else, 1

    Returns:
        * fig

    Example usage:
        signature_barplot_ID(W, np.sum(H))
    """
    W = W.copy()
    W.index = _map_id_sigs_back(W)
    for c in context83:
        if c not in W.index:
            W.loc[c] = 0
    W = W.loc[context83]
    sig_columns = [c for c in W if c.startswith('S')]
    if isinstance(contributions, pd.Series):
        W = W[sig_columns] * contributions[sig_columns]
    else:
        W = W[sig_columns] * contributions

    n_sigs = len(sig_columns)
    group_map = {'Cdel': [], 'Tdel': [], 'Cins': [], 'Tins': [],
                 '2del': [], '3del': [], '4del': [], '5+del': [],
                 '2ins': [], '3ins': [], '4ins': [], '5+ins': [],
                 '2delm': [], '3delm': [], '4delm': [], '5+delm': []}
    for x in W.index:
        group = re.search('.+?(?=[\d])', x).group(0)
        group_map[group].append(x)
    x_coords = {group: range(len(sigs)) for group, sigs in group_map.items()}

    color_map = {'Cdel': '#FFCC99', 'Tdel': '#FF8000', 'Cins': '#00FF00', 'Tins': '#00BB00',
                 '2del': '#FF99CC', '3del': '#FF3377', '4del': '#FF0000', '5+del': '#880000',
                 '2ins': '#99CCFF', '3ins': '#3377FF', '4ins': '#0000FF', '5+ins': '#000088',
                 '2delm': '#CC99FF', '3delm': '#9966FF', '4delm': '#8000FF', '5+delm': '#6000AA'}

    fig, axes = plt.subplots(nrows=n_sigs, ncols=16, figsize=(20, 2.5 * n_sigs), sharex='col',
                             sharey='row', gridspec_kw={'width_ratios': (6,) * 12 + (1, 2, 3, 5)})
    for row, sig in enumerate(sig_columns):
        for col, group in enumerate(group_map):
            if n_sigs == 1:
                ax = axes[col]
            else:
                ax = axes[row, col]
            bar_heights = W[sig].loc[group_map[group]]
            ax.bar(x_coords[group], bar_heights, width=.95, linewidth=1.5, edgecolor='gray', color=color_map[group],
                   rasterized=True)
            ax.set_xlim(-.55, x_coords[group][-1] + .55)
            if row == 0:
                ax.set_title(re.search('[\d+CT]+', group).group(0), color=color_map[group])
                if col == 0:
                    ax.text(44.5 / 6, 1.3, 'Mutational Signatures', transform=ax.transAxes,
                            horizontalalignment='center', fontsize=24)
                if group == 'Tdel':
                    ax.text(-.02, 1.16, '1bp deletions at repeats', fontsize=10, transform=ax.transAxes,
                            horizontalalignment='center', color=color_map[group])
                if group == 'Tins':
                    ax.text(-.02, 1.16, '1bp insertions at repeats', fontsize=10, transform=ax.transAxes,
                            horizontalalignment='center', color=color_map[group])
                if group == '4del':
                    ax.text(-.02, 1.16, '>1bp deletions at repeats', fontsize=10, transform=ax.transAxes,
                            horizontalalignment='center', color=color_map[group])
                if group == '4ins':
                    ax.text(-.02, 1.16, '>1bp insertions at repeats', fontsize=10, transform=ax.transAxes,
                            horizontalalignment='center', color=color_map[group])
                if group == '4delm':
                    ax.text(.8, 1.16, '>1bp deletions with microhomology', fontsize=10, transform=ax.transAxes,
                            horizontalalignment='center', color=color_map[group])
            if row < n_sigs - 1:
                ax.tick_params(axis='x', length=0)
            else:
                xlabels = [re.search('[\d+]+$', x).group(0) for x in group_map[group]]
                ax.set_xticks(x_coords[group])
                ax.set_xticklabels(xlabels, fontfamily='monospace')
                if col == 0:
                    ax.text(44.5 / 6, -.3, 'Motifs', transform=ax.transAxes, horizontalalignment='center', fontsize=20,
                            fontweight='bold')
            if col > 0:
                ax.tick_params(axis='y', length=0)
            if col == 15:
                ax.text(1.05, .5, sig, fontsize=14, rotation=270, transform=ax.transAxes, verticalalignment='center')

    plt.subplots_adjust(wspace=.08, hspace=.15)
    fig.text(.08, .5, 'Contributions', rotation='vertical', verticalalignment='center', fontsize=20, fontweight='bold')

    return fig
