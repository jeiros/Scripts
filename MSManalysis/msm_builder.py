#!/usr/bin/env python

"""
Trying to apply MSM analysis on my trajectories
http://msmbuilder.org/latest/
This script takes trajectories and topology from AMBER and saves in a directory
different transicion matrices for a combination of lag times and cluster
centers
"""

import mdtraj as md
from msmbuilder.dataset import dataset
import msmbuilder.cluster as cluster
import msmbuilder.msm as msm
import glob
import numpy as np
import os
import argparse
import sys


parser = argparse.ArgumentParser(usage="\n{} Trajectories*.nc Topology.prmtop".
                                 format(sys.argv[0]),
                                 formatter_class=argparse.
                                 ArgumentDefaultsHelpFormatter,
                                 epilog="Load up AMBER trajectories and their\
                                 corresponding topology with MDtraj. Cluster \
                                 them using K means and save transition\
                                 matrices in a specified folder")
parser.add_argument("Trajectories", help="An indefinite amount of AMBER\
                    trajectories", nargs="+")
parser.add_argument("Topology", help="The topology .prmtop file that matches\
                    the trajectories")
parser.add_argument("-o", "--out_folder", default="msmResults", help="Name of\
                    the folder in which the results are stored", metavar="")

args = parser.parse_args()

lag_times = [1, 20, 50, 100, 200, 400]
cluster_list = [2**n for n in range(3,11)]  # 2^n with n(3..10)


def loadTrajs(names, top):
    """Lazy loading of all trajectories (they're only loaded once they are
     used). Helps memory usagefor large datasets."""
    ds = dataset(names, fmt="mdtraj", topology=top)
    return ds


def printTrajs(ds):
    for traj in ds:
        print(traj)

# loadTrajs2 is a slower version of loadTrajs


# def loadTrajs2(names, topology):
#     filenames = sorted(glob.glob(names))
#     trajectories = [md.load(filename, top=topology) for filename in filenames]
#     return trajectories


def LabelTrajectories(TrajectorySet, nClusters):
    """Instatiate an object of the KCenters class and call it's fit method on
    the TrajectorySet obtained through loadTrajs. Returns a labeled trajectory
    (A sequence of integers corresponding to the index of the conformational
    state occupied by the system at each time point on a trajectory. )"""
    kcenters = cluster.KCenters(nClusters, "rmsd")
    kcenters.fit(TrajectorySet)
    labels = kcenters.transform(TrajectorySet)
    return labels


def MakeSeveralClusters(TrajectorySet, cluster_list):
    """Returns a list with clustering method using different number of
    clusters 2^n with n ranging from 3 to 9 (cluster_list)"""
    final_list = []

    for clusters in cluster_list:
        final_list.append(LabelTrajectories(TrajectorySet, clusters))

    return final_list


def MakeMSM(LabeledTrajs):
    """Returns a matrix with i columns and j rows. The element i,j of the
    matrix contains the Transition Matrix of the MSM created with the
    amount of clustering centers corresponding to the i row and using
    the lag time corresponding to the j column."""

    transmat_matrix = [[0 for x in range(len(lag_times))]
                       for x in range(len(LabeledTrajs))]

    for i in range(len(LabeledTrajs)):
        for j in range(len(lag_times)):
            MSM = msm.MarkovStateModel(lag_times[j], verbose=False)
            parameters = MSM.fit(LabeledTrajs[i])
            transmat_matrix[i][j] = parameters.transmat_
    return transmat_matrix


def saveTransMatrices(transmat_matrix):
    """Retrieve each of the i x j Transition Matrices and store each one
    in different files in the results directory. The files are named
    Transmat_Nclusters_Nlagtime"""

    for i in range(len(cluster_list)):
        for j in range(len(lag_times)):
            item = transmat_matrix[i][j]
            if os.path.exists("./{}".format(args.out_folder)):
                np.savetxt("./{0}/Transmat_{1}clusters_{2}lagtime.dat".
                           format(args.out_folder, str(cluster_list[i]),
                                  str(lag_times[j])), item)
            else:
                os.mkdir("./{}".format(args.out_folder))
                np.savetxt("./{0}/Transmat_{1}clusters_{2}lagtime.dat".
                           format(args.out_folder, str(cluster_list[i]),
                                  str(lag_times[j])), item)


def main():
    if args:
        ds = loadTrajs(args.Trajectories, args.Topology)
        clusters = MakeSeveralClusters(ds, cluster_list)
        MarkovModel = MakeMSM(clusters)
        saveTransMatrices(MarkovModel)

if __name__ == "__main__":
    import sys
    main()
