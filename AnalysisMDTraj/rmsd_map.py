#!/usr/bin/env python

import mdtraj
import argparse
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
import msmexplorer as msme

parser = argparse.ArgumentParser(prog='rmsd_map.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''version1''')

parser.add_argument('traj1', nargs='*')
parser.add_argument('-t2', '--traj2', required=False, nargs='*')
parser.add_argument('-p', '--prmtop', type=str, required=True)
parser.add_argument('-st', '--stride', type=int, required=False, default=1)
parser.add_argument('-sl', '--select', type=str, required=False, default='all')
parser.add_argument('-o', '--out_file', type=str, required=False,
                    default='matrix.pdf')
args = parser.parse_args()
# Timesteps which the trajes are saved at in disk
if args.traj2 is None:
    args.traj2 = args.traj1
tst1 = mdtraj.join(
    trajs=[
        mdtraj.load_frame(args.traj1[0], index=0, top=args.prmtop, atom_indices=[0]),
        mdtraj.load_frame(args.traj1[0], index=1, top=args.prmtop, atom_indices=[0]),
    ]
).timestep / 1000
tst2 = mdtraj.join(
    trajs=[
        mdtraj.load_frame(args.traj2[0], index=0, top=args.prmtop, atom_indices=[0]),
        mdtraj.load_frame(args.traj2[0], index=1, top=args.prmtop, atom_indices=[0]),
    ]
).timestep / 1000
print('tst1=%02fns\ttst2=%02fns' % (tst1, tst2))
# Timesteps which we have loaded with the trajs
TST1 = mdtraj.load(args.traj1[0], stride=args.stride,
                   top=args.prmtop, atom_indices=[0]).timestep / 1000
TST2 = mdtraj.load(args.traj2[0], stride=args.stride,
                   top=args.prmtop, atom_indices=[0]).timestep / 1000
print('TST1=%02fns\tTST2=%02fns' % (TST1, TST2))


def calculate_map(t1, t2):
    rmsd_matrix = np.empty((t1.n_frames, t2.n_frames))
    for i in range(t1.n_frames):
        rmsd_matrix[i, :] = mdtraj.rmsd(target=t2, reference=t1, frame=i)
    return rmsd_matrix * 10


def plot_matrix(mat, fout):
    f, ax = plt.subplots(figsize=(10, 10))
    if mat.shape[1] > 10:
        divisor_X = max(2, int(mat.shape[1] / 5))
    else:
        divisor_X = 1
    if mat.shape[0] > 10:
        divisor_Y = max(2, int(mat.shape[0] / 5))
    else:
        divisor_Y = 1
    # Only show multiples of divisor in the list
    x_ticks = list(map(lambda x: x if x % divisor_X == 0 else None, range(1, mat.shape[1] + 1)))
    y_ticks = list(map(lambda y: y if y % divisor_Y == 0 else None, range(1, mat.shape[0] + 1)))
    # Convert the elements of the list to strings
    x_ticks = list(map(lambda x: '{:0.0f}'.format(x * TST1) if type(x) == int else '', x_ticks))
    y_ticks = list(map(lambda y: '{:0.0f}'.format(y * TST2) if type(y) == int else '', y_ticks))

    sns.heatmap(data=mat, ax=ax, cmap='viridis',
                xticklabels=x_ticks, yticklabels=y_ticks,
                cbar_kws={'format': '%.1f',
                          'label': 'RMSD (Å)'})
    # Hide major ticks by setting their length to 0
    ax.tick_params(axis=u'both', which=u'both', length=0)
    plt.yticks(rotation='horizontal')
    ax.set(xlabel='Time (ns)', ylabel='Time (ns)')
    f.savefig(fout)

if __name__ == '__main__':
    args = parser.parse_args()
    print(args)
    print('Loading top')
    top = mdtraj.load_prmtop(args.prmtop)
    print('Selecting atoms')
    atom_indices = top.select(args.select)
    print('Loading file 1: %s' % args.traj1)
    traj1 = mdtraj.load(args.traj1, top=args.prmtop, stride=args.stride, atom_indices=atom_indices)
    print('Traj1 has %d frames' % traj1.n_frames)
    if args.traj2 is None:
        print('No secondary trajectory provided. Calculating self rmsd map.')
        traj2 = traj1
    else:
        print('Loading file 2: %s' % args.traj2)
        traj2 = mdtraj.load(args.traj2, top=args.prmtop, stride=args.stride, atom_indices=atom_indices)
        print('Traj2 has %d frames' % traj2.n_frames)
    print('Calculating and plotting RMSD map')
    m = calculate_map(traj1, traj2)
    plot_matrix(mat=m, fout=args.out_file)
    print('Done!\n')
