#!/usr/bin/env python

"""
http://msmbuilder.org/latest/examples/hmm-and-msm.html
"""

import pdb
import os
import matplotlib
from matplotlib.pyplot import *
from msmbuilder.featurizer import SuperposeFeaturizer
from msmbuilder.example_datasets import AlanineDipeptide
from msmbuilder.hmm import GaussianHMM
from msmbuilder.cluster import KCenters
from msmbuilder.msm import MarkovStateModel



dataset = AlanineDipeptide().get()
trajectories = dataset.trajectories
topology = trajectories[0].topology

indices = [atom.index for atom in topology.atoms if atom.element.symbol in ['C', 'O', 'N']]
featurizer = SuperposeFeaturizer(indices, trajectories[0][0])
sequences = featurizer.transform(trajectories)

lag_times = [1, 10, 20, 30, 40]
hmm_ts0 = {}
hmm_ts1 = {}
n_states = [3, 5]


for n in n_states:
    pdb.set_trace()
    hmm_ts0[n] = []
    hmm_ts1[n] = []
    for lag_time in lag_times:
        strided_data = [s[i::lag_time] for s in sequences for i in range(lag_time)]
        hmm = GaussianHMM(n_states=n, n_init=1).fit(strided_data)
        timescales = hmm.timescales_ * lag_time
        hmm_ts0[n].append(timescales[0])
        hmm_ts1[n].append(timescales[1])
        print('n_states=%d\tlag_time=%d\ttimescales=%s' % (n, lag_time, timescales))
    print()

