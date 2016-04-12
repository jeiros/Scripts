#!/usr/bin/env python

from traj_loading import load_Trajs_generator
import numpy as np
import mdtraj as md
import itertools
import seaborn as sns
from matplotlib import pyplot as plt
import argparse
import sys
import time

parser = argparse.ArgumentParser(usage="""{} Trajectories*.nc Topology.prmtop""".
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

parser.add_argument("-s", "--save", help="Save the plots as .eps images",
                    action="store_true")

parser.add_argument("-st", "--stride", help="""Stride for the loading of the
                    trajectory. Must be a divisor of the chunk.
                    Default value is 1.""", default=1, type=int)

parser.add_argument("-ch", "--chunk", help="""Number of frames that will be
                    used by md.iterload to load up the trajectories. Must be
                    a multiplier of the stride.
                    Default is 100 frames.""", default=100, type=int)

parser.add_argument("-t", "--title", help="""Name of the eps image where
                    plot is stored. Default is PCA.""", default="PCA.eps")

args = parser.parse_args()


def get_residuepairs(start1, end1, start2, end2):
    """
    Takes the beggining and end of two masks (0-indexed) and returns them as
    two lists, as well as the cartesian product between them.
    Parameters
    ----------
    start1: int
        Beggining of first mask
    end1: int
        End of first mask
    start2: int
        Beggining of second mask
    end2: int
        End of second mask
    Returns
    -------
    mask1: list of int
    mask2: list of int
    pairs: list of tuples
    """
    mask1 = list(range(start1, end1 + 1))
    mask2 = list(range(start2, end2 + 1))
    pairs = list(itertools.product(mask1, mask2))
    print('Mask1 has %d residues\n' % (end1 - start1 + 1))
    print('Mask2 has %d residues\n' % (end2 - start2 + 1))
    print('Number of residue-residue interactions that will be considered:%d\n'
          % len(pairs))
    return(mask1, mask2, pairs)


def cmap_MDtraj(traj_generator, mask1, mask2, pairs):
    print("Starting the cmap_MDtraj calculation...")
    frequency = np.zeros((len(mask1), len(mask2)))
    count = 0  # Keep the frame count
    for traj_chunk in traj_generator:
        count += traj_chunk.n_frames
        distances_inChunk = md.compute_contacts(traj_chunk, pairs, scheme='ca')
        column_sum = (distances_inChunk[0] <= 4).sum(0)  # Sum by column
        frequency += column_sum.reshape(len(mask1), len(mask2))  # Sum the partial result to frequency

    contact_frequency = frequency / count  # Total contact value for each residue-residue pair. From 0 to 1.
    print('Number of analyzed frames: %d\n' % count)
    print('Aggregate simulation time: %2.f ns\n' % (count * 0.02 * args.stride))
    return(contact_frequency)


def cmap_Cheng(traj_generator, mask1, mask2, pairs):
    print("Starting the cmap_Cheng calculation...")
    top = md.load_prmtop(args.Topology)
    frequency = np.zeros(len(pairs))
    count = 0  # Keep the frame count
    for traj_chunk in traj_generator:
        count += traj_chunk.n_frames
        index = 0  # To iterate through the residue-residue pair list
        for residue_pair in pairs:
            # Atom selection for each residue in the pair
            c_atoms_residue1 = top.select("resid %d and (type C)" %
                                          residue_pair[0])
            not_c_atoms_residue1 = top.select("resid %d and not type C" %
                                              residue_pair[0])
            c_atoms_residue2 = top.select("resid %d and type C" %
                                          residue_pair[1])
            not_c_atoms_residue2 = top.select("resid %d and not type C" %
                                              residue_pair[1])
            # Calculate all the possible distances between the C-C atoms and the
            # non C-C atoms. Results are stored in two np.arrays of shape:
            # (traj_chunk.n_frames, c_atoms_residue1*c_atoms_residue2)
            # (traj_chunk.n_frames, non_c_atoms_residue1*non_c_atoms_residue2)
            c_atoms_dist = md.compute_distances(traj_chunk,
                                                cartesianProduct(c_atoms_residue1,
                                                                 c_atoms_residue2))
            not_c_atoms_dist = md.compute_distances(traj_chunk,
                                                    cartesianProduct(not_c_atoms_residue1,
                                                                     not_c_atoms_residue2))
            # Implementation of the contact condition
            if (((c_atoms_dist <= 5.4).sum(1).any() > 0) and
               ((not_c_atoms_dist <= 4.6).sum(1).any() > 0)):
                frequency[index] += 1
            index += 1
    contact_frequency = (frequency/count).reshape(len(mask1), len(mask2))
    print('Number of analyzed frames: %d\n' % count)
    print('Aggregate simulation time: %2.f ns\n' % (count * 0.02 * args.stride))
    return(contact_frequency)


def plot_heatmap(contact_array, mask1, mask2):
    ax = sns.heatmap(contact_array,
                     vmin=0, vmax=1,
                     xticklabels=mask2,
                     yticklabels=mask1)
    ax.invert_yaxis()
    if args.save:
        plt.savefig(args.title, dpi=900, format='eps')
    else:
        sns.plt.show()


def cartesianProduct(x, y):
    """
    Implementation of cartesian product between two np.arrays
    Parameters
    ----------
    x, y: np.array
    Returns
    -------
    np.array with cartesian product

    """
    return(np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))]))


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
