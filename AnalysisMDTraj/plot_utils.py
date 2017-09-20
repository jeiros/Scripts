import msmexplorer as msme
from msmexplorer.utils import msme_colors
from matplotlib.ticker import FuncFormatter
from matplotlib import pyplot as pp
from msmbuilder.io import load_meta
from traj_utils import split_trajs_by_type
import numpy as np
import pandas as pd
import seaborn as sns
from glob import glob
from natsort import natsorted, order_by_index, index_natsorted


def split(a, n):
    """
    Splits a list into approximately equally sized n chunks
     """
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def get_diff_pkls(glob_expression):
    """
    Load several pkl files which have the MMGBSA information for a run.
    :param glob_expression: str, A str representing a glob pattern for the pkl files to load
    :return concat_split_melted_dfs: list, A list with 4 sublists, each with a pd.DataFrame with MMGBSA data
    """
    pickle_list = natsorted(glob(glob_expression))
    # Load each pickle as a pd.DataFrame
    df_list = [pd.read_pickle(pkl) for pkl in pickle_list]
    # Add a column to the data frame that tells us which run this data belongs to
    # This will be used later on for plotting
    for df, dir_path in zip(df_list, pickle_list):
        run_name = dir_path.split('/')[-3].split('_')[0]
        df['run'] = run_name
    # Recover the 'TOTAL' value, which we are interested in, and identify by the 'run' column
    melted_dfs = [df.melt(id_vars='run', value_vars='TOTAL', value_name='MMGBSA (kcal/mol)')
                  for df in df_list]
    # Split the list into four almost-equally-sized chunks, for later on plotting in a 4x4 plot
    concat_split_melted_dfs = [pd.concat(x) for x in split(melted_dfs, 4)]
    return concat_split_melted_dfs


def violin_plot_mmgbsa_results(ligand_data, figsize=(12, 12), title=None):
    """
    Returns a 4x4 figure with the violin plot distributions of MMGBSA energy of every run
    :param ligand_data: list, A list with pd.DataFrames of the MMGBSA energies for each run
    :param figsize: tuple, The figure size (in inches)
    :param title: str, The title of the figure
    :return f: a matplotlib figure
    """
    if len(ligand_data) != 4:
        raise ValueError('ligand_data must be a list of length 4')
    f, ((ax0, ax1), (ax2, ax3)) = pp.subplots(2, 2, figsize=figsize, sharey=False)
    for ax, df in zip((ax0, ax1, ax2, ax3), ligand_data):
        sns.violinplot(x='run', y='MMGBSA (kcal/mol)', data=df, ax=ax)
        ax.set_ylabel(ylabel='MMGBSA (kcal/mol)', size=14)
        ax.set_xlabel(xlabel='Run', size=14)
    if title is not None:
        f.suptitle(title, size=22)
    return f


def bar_plot_mmgbsa_results(excel_file, sort=True, titles=None):
    """
    Load data from an Excel file with the summary of the MMGBSA results, in a sheet which has to be called "MMGBSA".
    Create a plot for each ligand that is found under the 'Ligand' column in the table.
    :param excel_file: str, Name of the Excel file with the data
    :param sort: bool, Whether to sort the plot by increasing MMGBSA values
    :param titles: list, Name for each of the plots (as many as there are ligands in the table)
    :return f_list: list, A list of matplotlib figures
    """
    sns.set_style('whitegrid')
    df = pd.read_excel(excel_file, sheetname="MMGBSA")
    df = df.reindex(index=order_by_index(df.index, index_natsorted(df.Run)))
    ligands = df.Ligand.unique()
    f_list = []
    if titles is None:
        titles = [None for _ in ligands]
    elif len(titles) != len(ligands):
        raise ValueError('len of ligands and titles is not equal.')
    for lig, title in zip(ligands, titles):
        lig_df = df[df.Ligand == lig]
        if sort:
            lig_df.sort_values(by='MMGBSA (mean)', inplace=True)
        ax = lig_df.plot(x="Run", y="MMGBSA (mean)", yerr='MMGBSA (Std)', kind='bar',
                         legend=False, figsize=figure_dims(1400), title=title)
        ax.set_ylabel(ylabel='MMGBSA (kcal/mol)', size=14)
        ax.set_xlabel(xlabel='Run', size=14)
        f = pp.gcf()
        f_list.append(f)
    return f_list


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
        figsize = figure_dims(1200)
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
