#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(prog='smooth_traj.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''version1''')

parser.add_argument("top", help="A prmtop file", type=str)
parser.add_argument("traj", help="""A NetCDF MD trajectory""", nargs='+')
parser.add_argument("-n", "--name", help="The name of the output smoothed trajectoryaj", type=str,
                    default="traj_superposed.nc", required=False)
parser.add_argument("-r", "--ref", help="A reference structure to superpose to", type=str,
                    default=None, required=False)


if __name__ == '__main__':
    import mdtraj
    import os
    args = parser.parse_args()
    print(args)
    if len(args.traj) == 1:
        traj_name = os.path.basename(args.traj[0])[:-3]  # drop the .nc ending
        print('Loading traj...')
        traj = mdtraj.load(args.traj[0], top=args.top)
        print('Superposing...')
        # alpha carbons of the actin filament
        atoms = traj.topology.select('resid 1 to 6016 and name CA')
        if args.ref is None:
            traj.superpose(traj, 0, atom_indices=atoms)
        else:
            ref = mdtraj.load(args.ref)
            traj.superpose(ref, 0, atom_indices=atoms)
        traj.save_netcdf(''.join([traj_name, '_superposed.nc']))
    else:
        print('Loading {} trajs as one...'.format(len(args.traj)))
        traj = mdtraj.load(args.traj, top=args.top)
        print('Superposing...')
        # alpha carbons of the actin filament
        atoms = traj.topology.select('resid 1 to 6016 and name CA')
        if args.ref is None:
            traj.superpose(traj, 0, atom_indices=atoms)
        else:
            ref = mdtraj.load(args.ref)
            traj.superpose(ref, 0, atom_indices=atoms)
        traj.save_netcdf(args.name)
