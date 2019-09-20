import itertools
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from typing import Union

def plot_bar(H: pd.DataFrame, figsize: tuple(int,int) = (12,12)):
    """
    Plot stacked barchart & normalized version.

    Args:
        * H: matrix output from NMF
        * figsize: size of figure (int,int)

    Returns:
        * figure
    """
    H = H.iloc[:,:-3].copy()
    H['sum'] = H.sum(1)
    H = H.sort_values('sum', ascending=False)

    fig,axes = plt.subplots(2,1,figsize=figsize, sharex=True)

    H.iloc[:,:-1].plot(
        kind='bar',
        stacked=True,
        ax=axes[0],
        width=1.0
    )

    axes[0].set_xticklabels([])
    axes[0].set_xticks([])
    axes[0].set_ylabel('Counts', fontsize=20)

    H_norm = H.iloc[:,:-1].div(H['sum'].values,axis=0)
    H_norm.plot(
        kind='bar',
        stacked=True,
        ax=axes[1],
        width=1.0
    )

    axes[1].set_xticklabels([])
    axes[1].set_xticks([])
    axes[1].set_xlabel('Samples', fontsize=16)
    axes[1].set_ylabel('Fractions', fontsize=20)
    axes[1].get_legend().remove()

    return fig

def plot_k_dist(X: np.ndarray, figsize: tuple(int,int) = (8,8)):
    """
    Selected signatures plot.
    Args:
        X: numbers to plot
        figsize: size of figure (int,int)
    Returns:
        axis

    This may be called with:
        plot_k_dist(np.array([pd.read_hdf("output.h5","run{}/log".format(i)).K.iloc[-1] for i in range(250)]))

    """
    fig,ax = plt.subplots(figsize=figsize)

    sns.countplot(X, ax=ax, linewidth=2, edgecolor='k', rasterized=True)
    ax.set_ylim(0,ax.get_ylim()[1]+int(ax.get_ylim()[1]*0.1))

    ax.set_title("Aggregate of ARD-NMF (n={})".format(len(X)), fontsize=20)
    ax.set_ylabel("Number of Iterations", fontsize=18)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=18)

    return ax

def plot_signatures(W: pd.DataFrame, cohort: str, contributions: Union[int, pd.Series] = 1):
    """
    Plots signatures from W-matrix
    Args:
        W: W-matrix
        cohort: cohort name
        contributions: Series of total contributions from each signature if W is normalized, else 1
    """
    sig_columns = [c for c in W if c.startswith('S')]
    if isinstance(contributions, pd.Series):
        W = W[sig_columns] * contributions[sig_columns]
    else:
        W = W[sig_columns] * contributions
    n_sigs = len(sig_columns)
    change_map = {'CA': range(49, 65), 'CG': range(65, 81), 'CT': range(81, 97),
                  'TA': range(48, 32, -1), 'TC': range(32, 16, -1), 'TG': range(16, 0, -1)}
    color_map = {'CA': 'cyan', 'CG': 'red', 'CT': 'yellow', 'TA': 'purple', 'TC': 'green', 'TG': 'blue'}
    context_label = ['-'.join(p) for p in itertools.product('ACGT', 'ACGT')]
    x_coords = range(16)
    fig = plt.figure(figsize=(20, 2 * n_sigs))
    for row, sig in enumerate(sig_columns):
        ymax = 0
        ax_row = []
        for col, chg in enumerate(['CA', 'CG', 'CT', 'TA', 'TC', 'TG']):
            ax = fig.add_subplot(n_sigs, 6, row * 6 + col + 1)
            bar_heights = W[sig].loc[change_map[chg]]
            ax.bar(x_coords, bar_heights, width=.95, linewidth=1.5, edgecolor='gray', color=color_map[chg])
            ax.set_xlim(-.55, 15.55)
            ymax = max(ymax, ax.get_ylim()[1])
            if row == 0:
                ax.set_title('>'.join(chg), fontsize=18)
            if row < n_sigs - 1:
                ax.set_xticks([])
            else:
                ax.set_xticks(x_coords)
                ax.set_xticklabels(context_label, fontfamily='monospace')
                plt.xticks(rotation='vertical')
            if col > 0:
                ax.set_yticks([])
            if col == 5:
                ax.yaxis.set_label_position('right')
                ax.set_ylabel(sig, fontsize=18, rotation=270, labelpad=20)
            ax_row.append(ax)
        for a in ax_row:
            a.set_ylim(0, ymax)
    plt.subplots_adjust(wspace=.08, hspace=.15)
    plt.suptitle('Mutational Signatures in ' + cohort, fontsize=24, horizontalalignment='right')
    fig.text(.08, .5, 'Contributions', rotation='vertical', verticalalignment='center', fontsize=20, fontweight='bold')
    fig.text(.51, .03, 'Motifs', horizontalalignment='center', fontsize=20, fontweight='bold')
    with PdfPages('{}.signatures.pdf'.format(cohort)) as pdf:
        pdf.savefig(bbox_inches='tight')
    plt.savefig('{}.signatures.png'.format(cohort), bbox_inches='tight')
    plt.close()
