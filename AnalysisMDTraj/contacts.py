#!/usr/bin/env python
import mdtraj as md
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser(usage="""{} Trajectories*.nc Topology.prmtop""".
                                 format(sys.argv[0]),
                                 epilog="""Load up a list of AMBER NetCDF 
                                 trajectories and their corresponding topology
                                 with MDtraj. Calculates contact maps""")

parser.add_argument("Trajectories", help="""An indefinite amount of AMBER
                    trajectories""", nargs="+")

parser.add_argument("Topology", help="""The topology .prmtop file that matches
                    the trajectories""")
parser.add_argument("-st", "--stride", help="""Stride for the loading of the
                    trajectory. Must be a divisor of the chunk.
                    Default value is 1.""", default=1, type=int)

parser.add_argument("-ch", "--chunk", help="""Number of frames that will be 
                    used by md.iterload to load up the trajectories. Must be
                    a multiplier of the stride.
                    Default is 100 frames.""", default=100, type=int)
args = parser.parse_args()


def load_Trajs_generator(traj_list, topology, stride, chunk):
    for traj in traj_list:
        for frag in md.iterload(traj, chunk=chunk, top=topology,
                                stride=stride):
            yield frag


def get_pairs(start1, end1, start2, end2, topology):
    top = md.load_prmtop(topology)
    mask1 = np.array([atom.index for atom in top.atoms
                     if atom.element.symbol is 'C' and
                     atom.residue.index in range(start1, end1)])
    mask2 = np.array([atom.index for atom in top.atoms
                     if atom.element.symbol is 'C' and
                     atom.residue.index in range(start2, end2)])
    pairs = top.select_pairs(mask1, mask2)
    return pairs, mask1, mask2


def distances(traj_generator, pairs):
    for chunk in traj_generator:
        distance_generator = md.compute_distances(chunk, pairs)
        yield distance_generator


def contact_counter(distance_generator, pairs):
    distance_list_counter = [0] * len(pairs)
    frame = 0
    for chunk in distance_generator:
        for i in range(0, len(chunk)):
            frame += 1
            for j in range(0, len(pairs)):
                print(frame, '\t', chunk[i][j])
                if chunk[i][j] < 5.4:
                    distance_list_counter[j] += 1
    newList = [x / frame for x in distance_list_counter]
    return(newList)





def main():
    print('\n', args, '\n')
    if args:
        traj_generator = load_Trajs_generator(traj_list=sorted(args.Trajectories),
                                              topology=args.Topology,
                                              stride=args.stride,
                                              chunk=args.chunk)
        pairs = get_pairs(start1=0, end1=89, start2=395, end2=412, topology=args.Topology)
        dist_generator = distances(traj_generator, pairs)
        contact_list = contact_counter(dist_generator, pairs)
        print(contact_list)



if __name__ == "__main__":
    main()


import numpy as np
import mdtraj as md
from glob import glob
filenames = sorted(glob("05*.nc"))
top_str = "repstr.c0_phosS1P_nowat.prmtop"
top_obj = md.load_prmtop(top_str)

trj_gen = load_Trajs_generator(filenames, top_str, stride=2, chunk=100)
pairs = get_pairs(0, 89, 395, 412, top_str)
dis_gen = distances(trj_gen, pairs)

