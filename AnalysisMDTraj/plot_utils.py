import msmexplorer as msme
from msmexplorer.utils import msme_colors
from matplotlib.ticker import FuncFormatter
from matplotlib import pyplot as pp
from msmbuilder.io import load_meta
from traj_utils import split_trajs_by_type
import numpy as np


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


@msme_colors
def plot_tica_timescales(tica, meta, ax=None, color='cyan'):
    """
    Plot the timescales of a tica object

    :param tica: an msmbuilder tica object
    :param meta: an msmbuilder metadata object
    :param ax: a matplotlib axis
    :param color string: the color to plot

    :return ax: drawn matplotlib axis
    """
    if ax is None:
        ax = pp.gca()
    timestep = meta['step_ps'].unique()
    assert len(timestep) == 1, timestep
    timestep = float(timestep[0])  # ps
    to_us = (
        (1.0 / 1000)  # ps -> ns
        * (1.0 / 1000)  # ns -> us
        * (timestep / 1)  # steps -> ps
    )
    ax.hlines(
        tica.timescales_ * to_us,
        0, 1,
        color=color
    )
    ax.set_ylabel(r'Timescales / $\mathrm{\mu s}$', fontsize=18)
    ax.set_xticks([])
    ax.set_xlim((0, 1))
    return ax


def plot_ergodic_subspace(msm):
    pass

def plot_singletic_trajs(ttrajs, meta, system,
                         obs=(0, 1, 2), ylabels=None,
                         xlabel=None, title=None, figsize=None):
    """
    Plot each tIC vs. time on it's own axis and stack them vertically.
    By default, the first three tICs are plotted, but less (or more) can be
    chosen.

    :param ttrajs: dict, with keys as integers and values of np.arrays (tICA converted trajectories)
    :param meta: an msmbuilder metadata object
    :param system: str, the system inside 'meta' to plot
    :param obs: tuple (optional), the dimensions to plot
    :param ylabels: list or tuple (optional) of str, the labels of y axis
    :param xlabel: str (optional), the x label (to be shared)
    :param title: str (optional), the title of the plot
    :param figsize: tuple (optional), the figure dimensions
    :return axarr: an array of matplotlib axis
    """
    if ylabels is None:
        ylabels = ['tIC1', 'tIC2', 'tIC3']

    if len(obs) != len(ylabels):
        raise ValueError('Length of obs and ylabels is not equal.')

    def to_ns(x, pos):
        timestep = meta['step_ps'].unique()
        return "%d" % (x * timestep / 1000)


    # Get dictionary of specific sub system
    ttrajs_subtypes = split_trajs_by_type(ttrajs, meta)
    ttrajs_specific = ttrajs_subtypes[system]

    # Create the figure
    if figsize is None:
        figsize = figure_dims(1500)
    fig, axarr = pp.subplots(len(obs), 1, figsize=figsize, sharex=True, sharey=True)
    if title is not None:
        axarr[0].set(title=title)
    if xlabel is not None:
        axarr[-1].set(xlabel=xlabel)

    formatter = FuncFormatter(to_ns)
    indexes = meta[meta['type'] == system].index
    for j, ylabel in zip(range(len(obs)), ylabels):
        for index in indexes:
            ax = axarr[j]
            ax.plot(ttrajs_specific[index][:, obs[j]])
            ax.set_ylabel(ylabel)
            ax.xaxis.set_major_formatter(formatter)
    return axarr


def plot_overlayed_types(ttrajs, meta, obs=(0, 1), ax=None, stride=100,
                         xlabel=None, ylabel=None, plot_free_energy_kwargs=None, plot_kwargs=None):
    """
    Overlay each type of system inside the meta object onto the overall tICA
    free energy landscape.

    :param ttrajs: dict, with keys as integers and values of np.arrays (tICA converted trajectories)
    :param meta: an msmbuilder metadata object
    :param obs: tuple (optional), the dimensions to plot
    :param ax: matplotlib axis to plot in (optional)
    :param stride: int (optional, default=100), scatter every `stride` frames to avoid a cluttered plot
    :param xlabel: str (optional), the x label
    :param ylabel: str (optional), the y label
    :param plot_free_energy_kwargs: dict (optional), extra parameters to pass to msme.plot_free_energy
    :param plot_kwargs: dict (optional), extra parameters to pass to pyplot.plot
    :return ax: drawn matplotlib axis
    """
    if ax is None:
        ax = pp.gca()
    if plot_free_energy_kwargs is None:
        plot_free_energy_kwargs = {}
    if plot_kwargs is None:
        plot_kwargs = {}

    txx = np.concatenate(list(ttrajs.values()))
    ttrajs_subtypes = split_trajs_by_type(ttrajs, meta)
    msme.plot_free_energy(txx, obs=obs, ax=ax, **plot_free_energy_kwargs)

    for traj_id, traj_dict in ttrajs_subtypes.items():
        system_txx = np.concatenate(list(traj_dict.values()))
        ax.scatter(system_txx[::stride, obs[0]], system_txx[::stride, obs[1]], label=traj_id, **plot_kwargs)
    pp.legend()
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    return ax
