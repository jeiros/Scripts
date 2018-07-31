#!/usr/bin/env python

import mdtraj
import argparse
import msmexplorer as msme
from cycler import cycler
from itertools import cycle
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter
import seaborn as sns

parser = argparse.ArgumentParser(prog='rmsd.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''version1''')
parser.add_argument('prmtop', help="""Topology file matching Trajectories""")
parser.add_argument("Trajectories", help="""An indefinite amount of AMBER
                    trajectories""", nargs="+")
parser.add_argument('-st', '--stride', type=int, required=False, default=1)
parser.add_argument('-sl', '--select', type=str, required=False, default='all')
parser.add_argument('-o', '--out_file', type=str, required=False,
                    default='rmsd.pdf')


palette = cycler('color', sns.color_palette('colorblind', 10))
palette_cycled = cycle(palette)


def to_ns(x, pos):
    timestep = mdtraj.load_netcdf(args.Trajectories[0],
                                  args.prmtop, args.stride).timestep
    return '%d' % (int(x * timestep / 1000))


def plot_rsmd(traj_list, fout=None):
    rmsd_list = [mdtraj.rmsd(traj, traj, 0) * 10 for traj in traj_list]

    ax, side_ax = msme.plot_trace(rmsd_list[0], ylabel='RMSD (Å)', xlabel='Time (ns)',
                                  label=args.Trajectories[0][:-3],
                                  **next(palette_cycled))
    formatter = FuncFormatter(to_ns)
    ax.xaxis.set_major_formatter(formatter)
    if len(rmsd_list) > 1:
        for i, rmsd in enumerate(rmsd_list[1:]):
            msme.plot_trace(rmsd, ylabel='RMSD (Å)', xlabel='Time (ns)', ax=ax,
                            side_ax=side_ax,
                            label=args.Trajectories[i + 1][:-3],
                            **next(palette_cycled))

    if len(rmsd_list) > 5:
        ax.legend_.remove()
    sns.despine()
    f = plt.gcf()
    f.savefig(fout)


if __name__ == '__main__':
    args = parser.parse_args()
    print(args)
    top = mdtraj.load_prmtop(args.prmtop)
    atom_indices = top.select(args.select)
    trajs = [mdtraj.load(trj, top=args.prmtop, stride=args.stride,
                         atom_indices=atom_indices) for trj in args.Trajectories]
    print(trajs)
    plot_rsmd(trajs, args.out_file)
