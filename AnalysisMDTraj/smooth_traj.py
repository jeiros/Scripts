#!/usr/bin/env python
import argparse
from msmbuilder.io import backup
import mdtraj
import os
parser = argparse.ArgumentParser(prog='smooth_traj.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''version1''')

parser.add_argument("top", help="A prmtop file", type=str)
parser.add_argument("traj", help="""A NetCDF MD trajectory""", nargs='+')

parser.add_argument("-w", "--width", help="The width of the filter", type=int,
                    default=10, required=False)
parser.add_argument("-n", "--name", help="The name of the output smoothed trajectory", type=str,
                    default="traj_smoothed.nc", required=False)
parser.add_argument("-r", "--ref", help="A reference structure to superpose to", type=str,
                    default=None, required=False)

if __name__ == '__main__':
    args = parser.parse_args()
    print(args)
    if len(args.traj) == 1:
        traj_name = os.path.basename(args.traj[0])[:-3]  # drop the .nc ending
        print('Loading traj...')
        traj = mdtraj.load(args.traj[0], top=args.top)
        print('Superposing...')
        atoms = traj.topology.select('name CA')
        if args.ref is None:
            traj.superpose(traj, 0, atom_indices=atoms)
        else:
            ref = mdtraj.load(args.ref)
            traj.superpose(ref, 0, atom_indices=atoms)
        if args.width > 1:
            print('Smoothing...')
            traj.smooth(args.width, inplace=True)
        backup(args.name)
        traj.save_netcdf(args.name)
    elif len(args.traj) > 1:
        print('Loading {} trajs as one...'.format(len(args.traj)))
        traj = mdtraj.load(args.traj, top=args.top)
        print('Superposing...')
        atoms = traj.topology.select('name CA')
        if args.ref is None:
            traj.superpose(traj, 0, atom_indices=atoms)
        else:
            ref = mdtraj.load(args.ref)
            traj.superpose(ref, 0, atom_indices=atoms)
        if args.width > 1:
            print('Smoothing...')
            traj.smooth(args.width, inplace=True)
        backup(args.name)
        traj.save_netcdf(args.name)
    print('Done!')
