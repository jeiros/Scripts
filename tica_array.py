import msmexplorer as msme
from matplotlib.ticker import FuncFormatter
from matplotlib import pyplot as plt
import numpy as np


def plot_tic_array(tica_trajs, nrows, ncols, ts=0.2,
                   subplot_kwargs={}, trace_kwargs={}, free_energy_kwargs={}):
    '''
    Plot a nrows x cols array of the tIC projections, starting from the 1st one.

    Parameters
    ----------
    tica_trajs: list, or np.ndarray
        The tica transformed trajectory(ies)
    nrows: int
        Number of rows
    ncols: int
        Number of cols
    ts: float or int (default: 0.2)
        Timestep (in microseconds) between each frame in the trajectory
    subplot_kwargs: dict, optional
        Arguments to pass to plt.subplots
    trace_kwargs: dict, optional
        Arguments to pass to msme.plot_trace
    free_energy_kwargs: dict, optional
        Arguments to pass to msme.plot_free_energy

    Returns
    -------
    ax : matplotlib axis
        matplotlib figure axis

    '''

    def convert_to_mus(x, pos):
        'function for formatting the x axis of time trace plot'
        return x * ts

    if isinstance(tica_trajs, list):
        tica_trajs = np.concatenate(tica_trajs)
    elif not isinstance(tica_trajs, np.ndarray):
        raise ValuError('tica_trajs must be of type list or np.ndarray')

    fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, **subplot_kwargs)
    for i in range(nrows):
        for j in range(ncols):
            ax = ax_list[i][j]
            ax.grid(False)
            if i == j:
                msme.plot_trace(tica_trajs[:, i], ax=ax, **trace_kwargs)
                formatter = FuncFormatter(convert_to_mus)
                ax.xaxis.set_major_formatter(formatter)

                ax.grid(False, axis='y')
                ax.set_xlabel('Time ($\mu$s)')
                ax.sharey = False
            else:
                msme.plot_free_energy(tica_trajs, obs=(j, i), ax=ax,
                                      **free_energy_kwargs)
                # Bottom row
                if i == (nrows - 1):
                    ax.set_xlabel('tIC{}'.format(j + 1))
            # First column
            if j == 0:
                ax.set_ylabel('tIC{}'.format(i + 1))
    fig.tight_layout()
    return fig
