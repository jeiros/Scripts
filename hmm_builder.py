#!/usr/bin/env python
import argparse
import sys

import mdtraj as md
from msm_builder import loadTrajs
from msmbuilder.featurizer import SuperposeFeaturizer
from msmbuilder.hmm import GaussianFusionHMM

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


def makeHMM(TrajectorySet, Topology):
    alpha_carbons = [a.index for a in Topology.atoms if a.name == 'CA']
    f = SuperposeFeaturizer(alpha_carbons, Topology)

    dataset = []
    for traj in TrajectorySet:
        dataset.append(f.partial_transform(traj))
    hmm = GaussianFusionHMM(n_states = 8, n_features = len(alpha_carbons))
    hmm.fit(dataset)
    print(hmm.timescales_)





def main():
    if args:
        ds = loadTrajs(args.Trajectories, args.Topology)
        top = md.load_prmtop(args.Topology)
        hmm = makeHMM(ds, top)

if __name__ == "__main__":
    import sys
    main()