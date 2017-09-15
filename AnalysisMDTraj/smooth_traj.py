#!/usr/bin/env python
import argparse
import mdtraj

parser = argparse.ArgumentParser(prog='smooth_traj.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''version1''')

parser.add_argument("traj", help="""A NetCDF MD trajectory""", type=str)
parser.add_argument("top", help="A prmtop file", type=str)
parser.add_argument("-w", "--width", help="The width of the filter", type=int,
                    default=10, required=False)

if __name__ == '__main__':
    args = parser.parse_args()
    traj_name = args.traj.split('.')[-2]
    traj = mdtraj.load(args.traj, top=args.top)
    smooth_traj = traj.smooth(args.width)
    smooth_traj.save_netcdf(''.join([traj_name, '_smoothed.nc']))
