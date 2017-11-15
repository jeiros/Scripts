import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob

# Plotting settings
# plt.style.use('notebook')


def figure_dims(width_pt, factor=0.45):
    """
    I copied this from here:
    https://www.archer.ac.uk/training/course-material/2014/07/SciPython_Cranfield/Slides/L04_matplotlib.pdf
    """
    WIDTH = width_pt  # Figure width in pt (usually from LaTeX)
    FACTOR = factor  # Fraction of the width you'd like the figure to occupy
    widthpt = WIDTH * FACTOR
    inperpt = 1.0 / 72.27
    golden_ratio = (np.sqrt(5) - 1.0) / 2.0  # because it looks good
    widthin = widthpt * inperpt
    heightin = widthin * golden_ratio
    figdims = [widthin, heightin]  # Dimensions as list
    return figdims
