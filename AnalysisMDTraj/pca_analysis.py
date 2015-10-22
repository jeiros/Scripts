#!/usr/bin/env python
"""
Trying to implement the following example
http://mdtraj.org/latest/examples/pca.html
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

parser = argparse.ArgumentParser(usage="""{} Trajectories*.nc Topology.prmtop""".
                                 format(sys.argv[0]),
                                 epilog="""Load up a list of AMBER NetCDF 
                                 trajectories and their corresponding topology
                                 with MDtraj.Calculate a two-component PCA
                                 based on the pairwise distance between alpha
                                 carbons of the protein bakcbone. Plot a hexbin
                                 projection of the trajectories onto the PC space""")

parser.add_argument("Trajectories", help="""An indefinite amount of AMBER
                    trajectories""", nargs="+")

parser.add_argument("Topology", help="""The topology .prmtop file that matches
                    the trajectories""")

parser.add_argument("-s", "--save", help="Save the plots as .png images",
                    action="store_true")

parser.add_argument("-st", "--stride", help="""Stride for the loading of the
                    trajectory. Must be a divisor of the chunk.
                    Default value is 1.""", default=1)

parser.add_argument("-ch", "--chunk", help="""Number of frames that will be 
                    used by md.iterload to load up the trajectories. Must be
                    a multiplier of the stride. 
                    Default is 50 frames.""", default=1)

args = parser.parse_args()


def load_Trajs(names, topology, stride=1, chunk = 50):
    list_chunks = []
    for file in names:
        for frag in md.iterload(file, chunk = chunk, top = topology,
                                stride=stride):
            list_chunks.append(frag)
    return(list_chunks)

def pca_pwise_distance(list_chunks):
    pca = PCA(n_components=2)
    topology = list_chunks[0].topology
    ca_backbone = topology.select("backbone and name CA")
    pair_distances = []

    no_of_pwise_distances = ca_backbone.shape[0]*(ca_backbone.shape[0] -1 )/ 2
    print("""The number of pairwise distances to be calculated is: %s"""
          % int(no_of_pwise_distances))

    for chunk in list_chunks:
        pairs = topology.select_pairs(ca_backbone, ca_backbone)
        X = md.compute_distances(chunk, pairs)
        pair_distances.append(X)
    distance_array = np.concatenate(pair_distances)
    print(distance_array.shape)
    Y = pca.fit_transform(distance_array)
    return Y

def hex_plot(pca_array):
    PC1 = pca_array[:,0]
    PC2 = pca_array[:,1]
    plt.figure()
    plt.xlabel('PC1 (nm)')
    plt.ylabel('PC2 (nm)')
    plt.hexbin(x=PC1, y=PC2, bins='log', mincnt=1)
    if args.save:
        plt.savefig("PCA.png", dpi=600)
    else:
        plt.show()


def main():
    print('\n', args, '\n')
    if args:
        list_chunks = load_Trajs(sorted(args.Trajectories), args.Topology,
                                 stride=int(args.stride), chunk=int(args.chunk))
        pca_array = pca_pwise_distance(list_chunks)
        hex_plot(pca_array)


if __name__ == "__main__":
    main()
