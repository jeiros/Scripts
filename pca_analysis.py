#!/usr/bin/env python
"""
Trying to implement the following example
http://mdtraj.org/latest/examples/pca.html

Gives memory error for the moment
"""

import mdtraj as md
import sys
from glob import glob
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from itertools import combinations
import argparse

parser = argparse.ArgumentParser(usage="{} Trajectories*.nc Topology.prmtop".
                                 format(sys.argv[0]),
                                 epilog="Load up an AMBER trajectories and\
                                 their corresponding topology with MDtraj.\
                                 Plot the PCA done with sklearn's PCA class")

parser.add_argument("Trajectories", help="An indefinite amount of AMBER\
                    trajectories", nargs="+")
parser.add_argument("Topology", help="The topology .prmtop file that matches\
                    the trajectories")
parser.add_argument("-s", "--save", help="Save the plots as .png images",
                    action="store_true")
args = parser.parse_args()


def load_Trajs(names, topology):
    filenames = sorted(glob.glob(names))
    trajectories = md.load(filenames, top=topology)
    return trajectories


def pca_cartesian(filenames, topology):
    pca = PCA(n_components=2)
    first_frame = md.load_frame(filenames[0], 0, top=topology)
    pca_fits = []
    sim_frame = []
    for fragment in filenames:
        for chunk in md.iterload(fragment, chunk=100, top=topology):
            sim_frame.append(chunk.time)
            chunk_superposed = chunk.superpose(first_frame)
            pca_fits.append(pca.fit_transform(chunk.xyz.reshape(
                            chunk.n_frames, chunk.n_atoms * 3)))
    print(len(sim_frame), len(pca_fits))
    results = np.concatenate(pca_fits), np.concatenate(sim_frame)
    plt.scatter(results[0][:,0], results[0][:,1])

def pca_pwise_distance(trajectory):
    pca2 = PCA(n_components=2)
    atom_pairs = list(combinations(range(trajectory.n_atoms), 2))
    pairwise_distances = md.geometry.compute_distances(trajectory, atom_pairs)
    reduced_distances = pca2.fit_transform(pairwise_distances)
    plt.figure()
    plt.scatter(reduced_distances[:, 0], reduced_distances[:, 1],
                marker='x', c=trajectory.time)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('Pairwise distance PCA: alanine dipeptide')
    cbar = plt.colorbar()
    cbar.set_label('Time [ns]')
    if args.save:
        plt.savefig("test.png")
    else:
        plt.show()


def load_Trajectories(filenames, topology):
    pass


def PCA_hexbinplot(pca_data):
    pass


def rmsd(filenames, topology):
    rmsds = []
    first_frame = md.load_frame(filenames[0], 0, top=topology)
    for fragment in filenames:
        for chunk in md.iterload(fragment, chunk=100, top=topology):
            rmsds.append(md.rmsd(chunk, first_frame))
    return np.concatenate(rmsds)


def main():
    print('\n', args, '\n')
    if args:
        pca_cartesian(args.Trajectories, args.Topology)



if __name__ == "__main__":
    main()
