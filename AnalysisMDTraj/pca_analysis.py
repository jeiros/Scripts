#!/usr/bin/env python
"""
This scripts takes a list of NetCDF MD trajectories and their corresponding 
topology. Calculate a two-component PCA based on the pairwise distance between
alpha carbons of the protein bakcbone.
"""

import mdtraj as md
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import IncrementalPCA
import argparse

parser = argparse.ArgumentParser(usage="""{} Trajectories*.nc Topology.prmtop""".
                                 format(sys.argv[0]),
                                 epilog="""Load up a list of AMBER NetCDF 
                                 trajectories and their corresponding topology
                                 with MDtraj. Calculate a two-component PCA
                                 based on the pairwise distance between alpha
                                 carbons of the protein bakcbone. Plot a
                                 or scatter plot projection of the trajectories
                                 onto the PC space""")

parser.add_argument("Trajectories", help="""An indefinite amount of AMBER
                    trajectories""", nargs="+")

parser.add_argument("Topology", help="""The topology .prmtop file that matches
                    the trajectories""")

parser.add_argument("Plot_type", help="""The type of plot to be used. Can 
                    be a scatter plot or a hexbin plot.""",
                    choices=['scatter', 'hexbin'])

parser.add_argument("-s", "--save", help="Save the plots as .png images",
                    action="store_true")

parser.add_argument("-st", "--stride", help="""Stride for the loading of the
                    trajectory. Must be a divisor of the chunk.
                    Default value is 1.""", default=1, type=int)

parser.add_argument("-ch", "--chunk", help="""Number of frames that will be 
                    used by md.iterload to load up the trajectories. Must be
                    a multiplier of the stride.
                    Default is 100 frames.""", default=100, type=int)

parser.add_argument("-t", "--title", help="""Name of the png image where the PCA
                    plot is stored. Default is PCA.""", default="PCA")

args = parser.parse_args()


def load_Trajs(traj_list, topology, stride, chunk):
    """
    Iteratively loads a list of NetCDF files and returns them
    as a list of mdtraj.Trajectory objects

    Parameters
    ----------
    traj_list: list of str
            List with the names of trajectory files
    topology:  str
            Name of the topology file
    stride: int
            Frames to be used when loading the trajectories
    chunk:  int
            Number of frames to load at once from disk per iteration.
            If 0, load all.

    Returns
    -------
    list_chunks: list
            List of mdtraj.Trajectory objects, each of 'chunk' lenght
    """
    list_chunks = []
    for traj in traj_list:
        for frag in md.iterload(traj, chunk=chunk, top=topology,
                                stride=stride):
            list_chunks.append(frag)
    return(list_chunks)


def pca_pwise_distance(list_chunks, topology):
    """
    Takes a list of mdtraj.Trajectory objects and featurize them to backbone -
    Alpha Carbons pairwise distances. Perform 2 component Incremental 
    PCA on the featurized trajectory.

    Parameters
    ----------
    list_chunks: list of mdTraj.Trajectory objects
    topology: str
            Name of the Topology file

    Returns
    -------
    Y: np.array shape(frames, features)

    """
    pca = IncrementalPCA(n_components=2)
    top = md.load_prmtop(topology)
    ca_backbone = top.select("name CA")
    pairs = top.select_pairs(ca_backbone, ca_backbone)
    pair_distances = []
    for chunk in list_chunks:
        X = md.compute_distances(chunk, pairs)
        pair_distances.append(X)
    distance_array = np.concatenate(pair_distances)
    print("No. of data points: %d" % distance_array.shape[0])
    print("No. of features (pairwise distances): %d" % distance_array.shape[1])
    Y = pca.fit_transform(distance_array)
    return Y


# Plotting functions
def hex_plot(pca_array):
    PC1 = pca_array[:, 0]
    PC2 = pca_array[:, 1]
    plt.figure()
    plt.xlabel('PC1 (Å)')
    plt.ylabel('PC2 (Å)')
    plt.hexbin(x=PC1, y=PC2, bins='log', mincnt=1)
    cb = plt.colorbar()
    cb.set_label('log10(N)')
    if args.save:
        plt.savefig(args.title, dpi=600)
    else:
        plt.show()


def scatter_plot(pca_array):
    PC1 = pca_array[:, 0]
    PC2 = pca_array[:, 1]
    plt.figure()
    plt.xlabel('PC1 (Å)')
    plt.ylabel('PC2 (Å)')
    plt.scatter(x=PC1, y=PC2, marker='x')
    if args.save:
        plt.savefig(args.title, dpi=600)
    else:
        plt.show()


def main():
    print('\n', args, '\n')
    if args:
        list_chunks = load_Trajs(sorted(args.Trajectories), args.Topology,
                                 stride=args.stride, chunk=args.chunk)
        pca_array = pca_pwise_distance(list_chunks, args.Topology)
        if args.Plot_type == 'scatter':
            scatter_plot(pca_array)
        else:
            hex_plot(pca_array)


if __name__ == "__main__":
    main()
