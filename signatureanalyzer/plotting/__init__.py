"""
Signature Analyzer plotting API.
"""
from ._utils import series_to_colors
from ._utils import color_list_to_matrix_and_cmap

from ._rna import marker_heatmap

from ._muts import signature_barplot
from ._muts import signature_barplot_DBS
from ._muts import signature_barplot_ID
from ._muts import stacked_bar

from ._nmf import k_dist
from ._nmf import consensus_matrix

from ._cosine import cosine_similarity_plot
