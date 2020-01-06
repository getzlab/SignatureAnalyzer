import pandas as pd
import os
from typing import Union
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

def series_to_colors(s: pd.Series, color_palette_map: str = 'husl', cdict: dict = None):
    """
    Convert a pandas series to a color_map.
    -----------------------
    Args:
        * s: pd.Series
        * color_palette_map: str for color map
        * cdict: dict of custom value - color mappings
            Ex. cdict = {'Low': 'Green', 'Intermediate':'Yellow', 'High': 'Red'}

    Returns:
        * pd.Series: series mapped to RGB values
        * color_map: dictionary mapping unique values in series to color
    """
    if cdict is None:
        color_labels = s.unique()
        rgb_values = sns.color_palette(color_palette_map, len(color_labels))
        color_map = dict(zip(color_labels, rgb_values))
    else:
        color_map = cdict

    for key in color_map:
        if isinstance(key, type(None)) or (isinstance(key, float) and np.isnan(key)):
            color_map[key] = (1.0,1.0,1.0,1.0)

    return s.map(color_map), color_map

def color_list_to_matrix_and_cmap(colors: list, order_dict: dict = None):
    """
    Stripped from Seaborn.
    -----------------------
    Turns a list of colors into a numpy matrix and matplotlib colormap
    These arguments can now be plotted using heatmap(matrix, cmap)
    and the provided colors will be plotted.

    Args:
        * colors : list of matplotlib colors
            Colors to label the rows or columns of a dataframe.
        * order_dict: explicit color ordering to provide
            Ex. order_dict = {'Green': 0, 'Yellow': 1, 'Red': 2}

    Returns:
        * matrix : numpy.array
            A numpy array of integer values, where each corresponds to a color
            from the originally provided list of colors
        * cmap : matplotlib.colors.ListedColormap
    """
    all_colors = set(colors)
    m = len(colors)
    colors = [colors]

    if order_dict is None:
        color_to_value = dict((col, i) for i, col in enumerate(all_colors))
    else:
        color_to_value = order_dict

    matrix = np.array([color_to_value[c] for color in colors for c in color])
    matrix = matrix.reshape((1,m))
    return matrix, mpl.colors.ListedColormap(color_to_value.keys())
