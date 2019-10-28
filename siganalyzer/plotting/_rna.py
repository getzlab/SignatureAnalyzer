import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import Union
import numpy as np
import pandas as pd

from ._utils import series_to_colors
from ._utils import color_list_to_matrix_and_cmap

def marker_heatmap(
    X: pd.DataFrame,
    signatures: pd.DataFrame,
    order_series: pd.Series,
    subset_genes: Union[pd.Series,None] = None,
    diff: float = 0.5,
    max_norm: float = 0.5,
    figsize: tuple = (16,12),
    cmap: str ="YlGnBu",
    display_y: bool = False
    ):
    """
    Plot marker map.
    -----------------------------
    Args:
        * X: pd.DataFrame of input sample x feature matrix
        * signatures: pd.DataFrame signatures output;
            this bundles information about the weightings of each feature (ex. gene) and
            what signature they map to
        * order_series: series of samples mapping to subgroups
            index: X.index
            values: subgrouping
        * subset_series: a pd.Series with the index as the gene name or ID that
            matches the marker matrix & has a "Subgroup" column for labeling
        * diff: difference of loading for called signature vs. rest
        * max_norm: strength of loading for called signature
        * figsize: size of figure
        * cmap: colormap for plot

    Returns:
        * plt.Figure
    """
    # Filter for marker genes
    signatures_filt = signatures[(signatures['diff'] > diff) & (signatures['max_norm'] > max_norm)]

    # Remove signatures with no marker genes associated
    order_series = order_series[order_series.isin(set(signatures_filt['max_id'].astype(int)))]

    # Filter X matrix
    sample_markers = X.loc[signatures_filt.index, order_series.sort_values().index]

    # Set horizontal lines
    hz_lines = np.unique(sample_markers.join(signatures_filt).loc[:,'max_id'].values, return_index=True)[1]

    fig, ax = plt.subplots(figsize=figsize)

    x0 = ax.get_position().x0
    x1 = ax.get_position().x1
    y0 = ax.get_position().y0
    y1 = ax.get_position().y1
    buf = y1*0.01

    cbar_ax = fig.add_axes([.91, 0.5, .025, .3])

    sns.heatmap(sample_markers, ax=ax, cmap=cmap, rasterized=True, cbar_ax=cbar_ax)
    v,c = np.unique(order_series, return_counts=True)

    # plot horizontal lines
    _c = np.cumsum(c)
    _ci = np.roll(_c,2)
    _ci[0] = 0
    _ci[1] = 0
    ax.hlines(hz_lines, _ci, _c, rasterized=True)

    # plot vertical lines
    _h = list(hz_lines)
    _h.append(ax.get_ylim()[0])
    ax.vlines(np.cumsum(c)[:-1], _h[:-2], _h[2:], rasterized=True)
    ax.vlines(np.cumsum(c)[:-1], *ax.get_ylim(), alpha=0.4, rasterized=True)

    # set ticks
    ax.set_xticks(np.cumsum(c)-c/2)
    ax.set_xticklabels(v, rotation=360,fontsize=14)
    ax.set_yticks(np.arange(sample_markers.index.values.shape[0]))


    # add gene markings
    if subset_genes is not None:
        ax.set_yticks([])
        ax.set_yticklabels([], rasterized=True)


        lax = fig.add_axes([x0-3*buf, y0, 2*buf, y1-y0])
        lax.set_xticks([])
        lax.set_yticks([])

        meta = sample_markers.drop(columns=list(sample_markers)).join(subset_genes).iloc[:,0]

        colors_conversion, meta_colormap = series_to_colors(meta)
        meta_colormap_inv = dict([[v,k] for k,v in meta_colormap.items()])
        meta_colormap_inv = {(k[0],k[1],k[2]):v for k,v in meta_colormap_inv.items()}
        cbar_lax = fig.add_axes([0.06, 0.68, 2*buf, .2])

        # Add heatmapping
        mat,cmap = color_list_to_matrix_and_cmap(colors_conversion)
        sns.heatmap(
            mat.T,
            cmap=cmap,
            ax=lax,
            yticklabels=False,
            xticklabels=False,
            cbar=True,
            cbar_ax=cbar_lax,
        )

        cb_ticks = [float(t.get_text().replace('−','-')) for t in cbar_lax.get_yticklabels()]

        color_value_mapping = dict()
        for v in np.unique(mat):
            color_code = list(cmap.__call__(v))
            color_code = tuple(color_code[:3])
            color_value_mapping[v] = meta_colormap_inv[color_code]

        cbar_lax.get_yaxis().set_ticks([])

        n_labels = len(list(color_value_mapping.keys()))
        vals = [x * ((n_labels)/(n_labels+1)) + 0.5 * ((n_labels)/(n_labels+1)) for x in list(color_value_mapping.keys())]

        cbar_lax.get_yaxis().set_ticks(vals)
        cbar_lax.get_yaxis().set_ticklabels(list(color_value_mapping.values()),)
        cbar_lax.yaxis.set_ticks_position('left')

        cbar_lax.set_frame_on(True)

        lax.set_ylabel('Marker Genes', fontsize=14)
        ax.set_ylabel("")

        for _, spine in lax.spines.items():
            spine.set_visible(True)

    else:
        if display_y:
            ax.set_yticklabels(sample_markers.index.values, fontsize=5, rasterized=True)
        else:
            ax.set_yticks([])
            ax.set_yticklabels([], rasterized=True)

        ax.set_ylabel('Genes', fontsize=14)

    ax.set_title('')
    ax.set_xlabel('NMF Signatures', fontsize=14)

    cbar_ax.set_ylabel('Normalized Expression', fontsize=12)

    [spine.set_visible(True) for _, spine in ax.spines.items()]

    return fig
