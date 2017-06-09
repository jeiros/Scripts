#!/usr/bin/env python

import mdtraj
import argparse
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(prog='rmsd_map.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''version1''')

parser.add_argument('-t1', '--traj1', type=str, required=True)
parser.add_argument('-t2', '--traj2', type=str, required=True)
parser.add_argument('-p', '--prmtop', type=str, required=True)
parser.add_argument('-st', '--stride', type=int, required=False, default=1)
parser.add_argument('-o', '--out_file', type=str, required=False,
                    default='matrix.pdf')


def calculate_map(t1, t2):
    rmsd_matrix = np.empty((t1.n_frames, t2.n_frames))
    for i in range(t1.n_frames):
        rmsd_matrix[i, :] = mdtraj.rmsd(target=t2, reference=t1, frame=i)
    return rmsd_matrix * 10


def plot_matrix(mat, fout):
    f, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(data=mat, ax=ax, cmap='viridis')
    f.savefig(fout)

if __name__ == '__main__':
    args = parser.parse_args()
    print('Loading file 1: %s' % args.traj1)
    traj1 = mdtraj.load_netcdf(args.traj1, args.prmtop, args.stride)
    print('Traj1 has %d frames' % traj1.n_frames)
    print('Loading file 2: %s' % args.traj2)
    traj2 = mdtraj.load_netcdf(args.traj2, args.prmtop, args.stride)
    print('Traj2 has %d frames' % traj2.n_frames)
    print('Calculating and plotting RMSD map')
    plot_matrix(mat=calculate_map(traj1, traj2), fout=args.out_file)
    print('Done!\n')
