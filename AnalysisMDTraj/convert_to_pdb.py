#!/usr/bin/env python
from glob import glob
from subprocess import call, PIPE
import os
import argparse
from traj_utils import write_cpptraj_script

parser = argparse.ArgumentParser(prog='convert_to_pdb.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''version0.1''')

parser.add_argument("Trajectories", help="""A glob expression of .nc file trajectories""", type=str)
parser.add_argument("Topology", help="""The topology .prmtop file that matches
                    the trajectory""", type=str)


if __name__ == '__main__':
    args = parser.parse_args()
    traj_fns = glob(args.Trajectories)
    top_fn = args.Topology
    if not os.path.exists('pdbs'):
        os.makedirs('pdbs')
    for traj in traj_fns:
        print('Converting {}'.format(traj))
        with open('cmds.cpptraj', 'w') as f:
            f.write(write_cpptraj_script(traj, top_fn))
        call(['cpptraj', '-i', 'cmds.cpptraj'], stdout=PIPE)
        if os.path.exists('cmds.cpptraj'):
            os.remove('cmds.cpptraj')
