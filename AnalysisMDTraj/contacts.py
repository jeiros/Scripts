#!/usr/bin/env python

import argparse
import sys
import time
from contact_utils import *
from traj_utils import load_Trajs_generator


parser = argparse.ArgumentParser(usage="""{} Trajs*.nc Topology.prmtop""".
                                 format(sys.argv[0]),
                                 epilog="""Load up a list of AMBER NetCDF
                                 trajectories and their corresponding topology
                                 with MDtraj. Calculate a heatmap between the
                                 two specified masks""")

parser.add_argument("Trajectories", help="""An indefinite amount of AMBER
                    trajectories""", nargs="+")

parser.add_argument("Topology", help="""The topology .prmtop file that matches
                    the trajectories""")

parser.add_argument("-s1", "--start1", help="""0-indexed value for the first
                    residue of Mask1""", type=int)

parser.add_argument("-e1", "--end1", help="""0-indexed value for the final
                    residue of Mask1""", type=int)

parser.add_argument("-s2", "--start2", help="""0-indexed value for the first
                    residue of Mask2""", type=int)

parser.add_argument("-e2", "--end2", help="""0-indexed value for the final
                    residue of Mask2""", type=int)
parser.add_argument("Map_type", help="""Type of calculation for the contact map.
                    Can be either mdtraj or cheng style.""", choices=['mdtraj',
                                                                      'cheng'])

parser.add_argument("-s", "--save", help="Save the plots as .png images",
                    action="store_true")

parser.add_argument("-st", "--stride", help="""Stride for the loading of the
                    trajectory. Must be a divisor of the chunk.
                    Default value is 1.""", default=1, type=int)

parser.add_argument("-ch", "--chunk", help="""Number of frames that will be
                    used by md.iterload to load up the trajectories. Must be
                    a multiplier of the stride.
                    Default is 100 frames.""", default=100, type=int)

parser.add_argument("-t", "--title", help="""Name of the image where
                    plot is stored. Default is PCA.""", default="PCA.png")

args = parser.parse_args()


def main():
    if args:
        print('\n', args, '\n')
        start = time.time()
        mask1, mask2, pairs = get_residuepairs(args.start1, args.end1,
                                               args.start2, args.end2)
        trjs = load_Trajs_generator(sorted(args.Trajectories),
                                    prmtop_file=args.Topology,
                                    stride=args.stride,
                                    chunk=args.chunk)
        print("Trajectories have been loaded after %.2f s.\n" %
              (time.time() - start))
        if args.Map_type == 'mdtraj':
            cmap = cmap_MDtraj(trjs, mask1, mask2, pairs)
            print("MDtraj style cmap has been calculated after %.2f s.\n" %
                  (time.time() - start))
        else:
            cmap = cmap_Cheng(trjs, mask1, mask2, pairs)
            print("Cheng style cmap has been calculated after %.2f s.\n" %
                  (time.time() - start))
        plot_heatmap(cmap, mask1, mask2)
        print("Total execution time: %.2f s." % (time.time() - start))

if __name__ == "__main__":
    main()
