import pandas as pd
import os
from typing import Union
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

def series_to_colors(s: pd.Series, color_palette_map: str = 'husl'):
    """
    Convert a pandas series to a color_map.
    -----------------------
    Args:
        * s: pd.Series
        * color_palette_map: str for color map
    Returns:
        * pd.Series: series mapped to RGB values
    """
    color_labels = s.unique()
    rgb_values = sns.color_palette(color_palette_map, len(color_labels))
    color_map = dict(zip(color_labels, rgb_values))

    for key in color_map:
        if isinstance(key, type(None)) or (isinstance(key, float) and np.isnan(key)):
            color_map[key] = (1.0,1.0,1.0,1.0)

    return s.map(color_map), color_map

def color_list_to_matrix_and_cmap(colors: list, axis=0):
    """
    Stripped from Seaborn.
    -----------------------
    Turns a list of colors into a numpy matrix and matplotlib colormap
    These arguments can now be plotted using heatmap(matrix, cmap)
    and the provided colors will be plotted.

    Args:
        * colors : list of matplotlib colors
            Colors to label the rows or columns of a dataframe.
    Returns:
        * matrix : numpy.array
            A numpy array of integer values, where each corresponds to a color
            from the originally provided list of colors
        * cmap : matplotlib.colors.ListedColormap
    """
    all_colors = set(colors)
    m = len(colors)
    colors = [colors]

    color_to_value = dict((col, i) for i, col in enumerate(all_colors))
    matrix = np.array([color_to_value[c] for color in colors for c in color])

    matrix = matrix.reshape((1,m))
    return matrix, mpl.colors.ListedColormap(all_colors)
