import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.cluster
import pandas as pd

def cosine_similarity_plot(cosine_similarity_matrix: pd.DataFrame):
    """
    Plots cosine similarity heatmap
    --------------------------------------
    Args:
        * cosine_similarity_matrix: DataFrame with cosine similarities of each component to each signature

    Returns:
        * fig

    Example usage:
        cosine_similarity_plot(pd.read_hdf('nmf_output.h5', 'cosine'))
    """

    linkage_mat = scipy.cluster.hierarchy.linkage(cosine_similarity_matrix, metric='euclidean')
    sig_order = scipy.cluster.hierarchy.dendrogram(linkage_mat, no_plot=True)['leaves']
    sorted_cosine_similarity_matrix = cosine_similarity_matrix.iloc[sig_order]

    fig, (ax, cbar_ax) = plt.subplots(figsize=(cosine_similarity_matrix.shape[1] + .5, cosine_similarity_matrix.shape[0] / 3),
                                      ncols=2, gridspec_kw={'width_ratios': [cosine_similarity_matrix.shape[1] * 2, 1]})
    sns.heatmap(sorted_cosine_similarity_matrix, vmin=.75, vmax=1, ax=ax, cbar_ax=cbar_ax, cmap='viridis', annot=True)
    ax.set_yticks(np.arange(cosine_similarity_matrix.shape[0]) + .5)
    ax.set_yticklabels(sorted_cosine_similarity_matrix.index)
    ymax, ymin = ax.get_ylim()
    ax.set_ylim(ymax + .5, ymin - .5)
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    cbar_ax.set_frame_on(True)

    return fig
