#!/usr/bin/env python
import argparse
import sys
from glob import glob
import mdtraj as md
from msmbuilder.featurizer import SuperposeFeaturizer
from msmbuilder.hmm import GaussianHMM

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


def makeHMM(Trajectories, topology):
    top = md.load_prmtop(topology)
    alpha_carbons = [a.index for a in top.atoms if a.name == 'CA']
    filenames = sorted(glob(Trajectories))
    first_frame = md.load_frame(filenames[0], 0, top=top)

    f = SuperposeFeaturizer(alpha_carbons, first_frame)
    dataset = []
    for fragment in filenames:
            for chunk in md.iterload(fragment, chunk=100, top=top):
                dataset.append(f.partial_transform(chunk))
    hmm = GaussianHMM(n_states=8)
    hmm.fit(dataset)
    print(hmm.timescales_)
    return hmm


def main():
    if args:
        hmm = makeHMM(args.Trajectories, args.Topology)

if __name__ == "__main__":
    main()