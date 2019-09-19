import itertools
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from typing import Union


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
