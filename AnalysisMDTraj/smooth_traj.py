#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(prog='smooth_traj.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''version1''')

parser.add_argument("top", help="A prmtop file", type=str)
parser.add_argument("traj", help="""A NetCDF MD trajectory""", nargs='+')

parser.add_argument("-w", "--width", help="The width of the filter", type=int,
                    default=10, required=False)
parser.add_argument("-n", "--name", help="The name of the output smoothed trajectoryaj", type=str,
                    default="traj_smoothed.nc", required=False)

if __name__ == '__main__':
    import mdtraj
    from traj_utils import load_Trajs
    import os
    args = parser.parse_args()
    print(args)
    if len(args.traj) == 1:
        traj_name = os.path.basename(args.traj[0])[:-3]  # drop the .nc ending
        print('Loading traj...')
        traj = mdtraj.load(args.traj[0], top=args.top)
        print('Superposing...')
        traj.superpose(traj, 0)
        print('Smoothing...')
        smooth_traj = traj.smooth(args.width)
        smooth_traj.save_netcdf(''.join([traj_name, '_smoothed.nc']))
    elif len(args.traj) > 1:
        print('Loading {} trajs as one...'.format(len(args.traj)))
        traj = mdtraj.load(args.traj, top=args.top)
        print('Superposing...')
        traj.superpose(traj, 0)
        print('Smoothing...')
        smooth_traj = traj.smooth(args.width)
        smooth_traj.save_netcdf(args.name)
    print('Done!')
