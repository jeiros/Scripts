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
import os


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
        ax.set_ylabel(ylabel=r'$\Delta$G binding (kcal/mol)', size=14)
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
    # sns.set_style('whitegrid')
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
        overall_mean = lig_df['MMGBSA (mean)'].mean()
        overall_std = lig_df['MMGBSA (mean)'].std()
        print("{} {:02f} {:02f}".format(lig, overall_mean, overall_std))
        xmin, xmax = ax.get_xlim()
        # Mean line
        ax.plot(
            [xmin, xmax], [overall_mean, overall_mean],
            linewidth=1.5,
            color='blue'
        )
        # Upper std bar
        ax.plot(
            [xmin, xmax],
            [overall_mean + overall_std, overall_mean + overall_std],
            linestyle='dashed',
            linewidth=1,
            color='blue'
        )
        # Lower std bar
        ax.plot(
            [xmin, xmax],
            [overall_mean - overall_std, overall_mean - overall_std],
            linestyle='dashed',
            linewidth=1,
            color='blue'
        )
        ax.set_ylim(top=0)
        ax.set_ylabel(ylabel=r'$\Delta$G binding (kcal/mol)', size=14)
        ax.set_xlabel(xlabel='Run', size=14)
        f = pp.gcf()
        f.tight_layout()
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


def plot_ergodic_subspace(msm, clusterer, obs=(0, 1), ax=None, alpha=1.,
                          color='blue', label=None, xlabel=None, ylabel=None,
                          scatter_kwargs=None):
    """
    Plot which cluster centers out of the clusterer object have been visited
    in the msm object.
    :param msm: A trained msmbuilder MSM object
    :param clusterer: A trained msmbuilder clusterer object
    :param obs: tuple, which dimensions to plot
    :param ax: a matplotlib.axes object
    :param alpha: float, transparency parameter for ax.scatter
    :param color: string, parameter for ax.scatter
    :param label: string, parameter for ax.scatter
    :param xlabel: string, parameter for ax.scatter
    :param ylabel: string, parameter for ax.scatter
    :param scatter_kwargs: dict, any other parameters for ax.scatter
    :return ax:  a matplotlib.axes object
    """
    if ax is None:
        ax = pp.gca()
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if scatter_kwargs is None:
        scatter_kwargs = {}

    prune = clusterer.cluster_centers_[:, obs]
    ax.scatter(
        prune[msm.state_labels_][:, 0],
        prune[msm.state_labels_][:, 1],
        color=color,
        label=label,
        alpha=alpha,
        **scatter_kwargs
    )

    return ax


def plot_singletic_trajs(ttrajs, meta, system, alpha=1,
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
            ax.plot(ttrajs_specific[index][:, obs[j]], alpha=alpha)
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
    pp.legend(loc='best')
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    return ax


def plot_microstates(msm, txx, clusterer, obs=(0, 1), eigenvector=1, ax=None,
                     clabel=None):
    """
    Taken from the msmbuilder template.
    Plot the microstate centers of an msm object on top of a grey hexbin
    tICA landscape. Color the microstates by the sign of the chosen eigenvector
    :param ax: matplotlib axis to plot in (optional)
    :param msm:
    :param txx:
    :param clusterer:
    :param obs:
    :param eigenvector:
    :param clabel:
    :return:
    """
    if ax is None:
        ax = pp.gca()

    txx = txx[:, obs]
    ax.hexbin(txx[:, 0], txx[:, 1],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    scale = 100 / np.max(msm.populations_)
    add_a_bit = 2
    prune = clusterer.cluster_centers_[:, obs]
    c = ax.scatter(prune[msm.state_labels_, 0],
                   prune[msm.state_labels_, 1],
                   s=scale * msm.populations_ + add_a_bit,
                   c=msm.left_eigenvectors_[:, eigenvector],
                   cmap='RdBu'
                   )
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 2", fontsize=16)
    pp.colorbar(c, label=clabel)
    return ax


def plot_src_sink(msm, clusterer, ev, txx, src, sink, clabel=None, title='', ax=None):
    """

    :param msm:
    :param clusterer:
    :param ev:
    :param txx:
    :param src:
    :param sink:
    :param clabel:
    :param title:
    :param ax:
    :return:
    """

    if ax is None:
        ax = pp.gca()
    ax.set_title(title)
    plot_microstates(msm=msm, eigenvector=ev, clabel=clabel, txx=txx, clusterer=clusterer, ax=ax)
    # source
    ax.scatter(
        clusterer.cluster_centers_[src][0],
        clusterer.cluster_centers_[src][1],
        marker='D',
        color='red',
        s=200,
    )
    # sink
    ax.scatter(
        clusterer.cluster_centers_[sink][0],
        clusterer.cluster_centers_[sink][1],
        marker='D',
        color='blue',
        s=200,
    )
    return ax


def plot_efhand_dists_src_sinks(src_glob, sink_glob, title=None, ax=None):
    if ax is None:
        ax = pp.gca()

    src = glob(src_glob)
    snk = glob(sink_glob)

    data = {}
    for f in src:
        atom = os.path.basename(f).split('-')[1].split('_')[0]
        data[atom] = [None, None]

    for f1, f2 in zip(src, snk):
        atom = os.path.basename(f1).split('-')[1].split('_')[0]
        src_mean = np.loadtxt(f1).mean(axis=0)[1]
        src_std = np.loadtxt(f1).std(axis=0)[1]
        snk_mean = np.loadtxt(f2).mean(axis=0)[1]
        snk_std = np.loadtxt(f2).std(axis=0)[1]

        data[atom][0] = (src_mean, src_std)
        data[atom][1] = (snk_mean, snk_std)

    # cpptraj names E32 of cTni as E280, rename that
    data['E32O'] = data.pop('E280O')
    data['E32OE1'] = data.pop('E280OE1')
    data['E32OE2'] = data.pop('E280OE2')

    # now do plot
    ax.set_title(title)
    n_feats_plot = len(data)
    xx = np.arange(1, n_feats_plot + 1)
    ax.errorbar(xx, [x[0][0] for x in data.values()],
                yerr=[x[0][1] for x in data.values()],
                label='source', linestyle='None', marker='s', color='red')
    ax.errorbar(xx, [x[1][0] for x in data.values()],
                yerr=[x[1][1] for x in data.values()],
                label='sink', linestyle='None', marker='*', color='blue')
    pp.legend()
    ax.set_xticks(xx)
    ax.set_xlim((0, n_feats_plot + 1))
    ax.set_xticklabels(
        list(data.keys())
    )
    for tick in ax.get_xticklabels():
        tick.set_rotation(60)
    ax.tick_params(labelsize=14)
    ax.set_ylabel("Distance (Å)", fontsize=16)
    return ax


def plot_cluster_centers(clusterer, centers, txx, ax=None, obs=(0, 1), from_clusterer=True, msm=None, add_bit=20):
    """
    Plots cluster centers as scatter points on top of a tICA landscape. Center IDs can be either in MSM labeling
    or directly from the clustering labelling.
    :param clusterer: a fit clusterer object, with .cluster_centers_ attribute
    :param centers: list of ints, the IDs of the cluster centers to plot
    :param txx: np.array of concatenated tIC trajs, shape = (n_frames, n_features)
    :param ax: a matplotlib axes object (optional)
    :param obs: tuple of ints, dimensions to plot (optional)
    :param from_clusterer: bool, are the centers id in clusterer indexing or msm?
        default True means IDs are from clusterer
    :param msm: a MarkovStateModel object which has been fit
    :param add_bit: control size of scatter plots
    :return ax: a matplotlib axes object
    :raises ValueError: if we select from_clusterer=False it means that the center IDs are in MSM-internal labeling,
        so we need an MSM object to map those back to the clusterer naming
    """
    if ax is None:
        pp.gca()
    ax = msme.plot_free_energy(txx, obs=obs, n_samples=5000,
                               gridsize=100, vmax=5.,
                               n_levels=8, cut=5, xlabel='tIC1',
                               ylabel='tIC2'
    )
    prune = clusterer.cluster_centers_[:, obs]
    if from_clusterer:
        chosen_centers = prune[centers]
        ax.scatter(chosen_centers[:, 0], chosen_centers[:, 1])
    else:
        if msm is None:
            raise ValueError('if from_clusterer is False please provide a fit MSM in the msm parameter')
        else:
            # Retrieve cluster centers from clusterer objects that have been used in this MSM
            centers_in_clusterer = []
            for k, v in msm.mapping_.items():
                if v in centers:
                    centers_in_clusterer.append(k)

            scale = 100 / np.max(msm.populations_)
            chosen_centers = prune[centers_in_clusterer]
            ax.scatter(chosen_centers[:, 0], chosen_centers[:, 1],
                       s=add_bit + (scale * msm.populations_),
            )

    return ax


def plot_tpt(msm, clusterer, txx, ev=1, ax=None, title=None, obs=(0, 1), num_paths=1):
    """
    Automatically plot a TPT plot of msmexplorer by selecting the microstates that have
    lowest (source) and highest (sink) value of the provided eigenvector
    :param msm: an MSM object
    :param clusterer: a clusterer object
    :param txx: np.array of concatenated tIC trajs, shape = (n_frames, n_features)
    :param ev: int, the eigenvector to plot the transition for
    :param ax: a matplotlib axes object (optional)
    :param title: str, the title
    :param obs: tuple of ints, dimensions to plot (optional)
    :param num_paths: int, the number of paths to plot
    :return ax:
    """
    if ax is None:
        pp.gca()
    prune = clusterer.cluster_centers_[:, obs]
    msm_states = prune[msm.state_labels_]
    pos = dict(zip(range(len(msm_states)), msm_states))
    w = (msm.left_eigenvectors_[:, ev] - msm.left_eigenvectors_[:, ev].min())
    w /= w.max()
    src, snk = get_source_sink(msm, clusterer=clusterer, eigenvector=ev)
    print('src', src, 'snk', snk)
    ax = msme.plot_free_energy(txx, obs=obs, n_samples=10000,
                               gridsize=100, vmax=5.,
                               n_levels=8, cut=5, xlabel='tIC1',
                               ylabel='tIC2', ax=ax)

    cmap = msme.utils.make_colormap(['rawdenim', 'lightgrey', 'pomegranate'])
    ax = msme.plot_tpaths(msm, [src], [snk], pos=pos, node_color=cmap(w),
                          alpha=.9, edge_color='black', ax=ax, num_paths=num_paths)
    if title is not None:
        ax.set_title(title)
    return ax