#!/usr/bin/env python

"""
http://msmbuilder.org/latest/examples/hmm-and-msm.html
"""
import pdb
import os
import matplotlib
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
plt.style.use("ggplot")
from msmbuilder.featurizer import SuperposeFeaturizer
from msmbuilder.example_datasets import AlanineDipeptide
from msmbuilder.hmm import GaussianHMM
from msmbuilder.cluster import KCenters
from msmbuilder.msm import MarkovStateModel
from msmbuilder.dataset import dataset
import mdtraj as md
from glob import glob




filenames = sorted(glob("05_Prod_*.nc"))
topology = md.load_prmtop(glob("*nowat.prmtop")[0])

first_frame = md.load_frame(filenames[0], 0, top=topology)
indices = [atom.index for atom in topology.atoms if atom.element.symbol in ['C', 'O', 'N']]
featurizer = SuperposeFeaturizer(indices, first_frame)
sequences = []

for fragment in filenames:
    for chunk in md.iterload(fragment, chunk = 100, top = topology):
        sequences.append(featurizer.partial_transform(chunk))

lag_times = [1, 20, 50, 100, 200, 400]
hmm_ts0 = {}
hmm_ts1 = {}
n_states = [n for n in range(2,11,2)]


for n in n_states:
    hmm_ts0[n] = []
    hmm_ts1[n] = []
    for lag_time in lag_times:
        #pdb.set_trace()
        strided_data = [s[i::lag_time] for s in sequences for i in range(lag_time)]
        hmm = GaussianHMM(n_states=n, n_init=1).fit(strided_data)
        timescales = hmm.timescales_ * lag_time
        if len(timescales) > 1:
            hmm_ts0[n].append(timescales[0])
            hmm_ts1[n].append(timescales[1])
        else:
            hmm_ts0[n].append(timescales[0])
        print('n_states=%d\tlag_time=%d\ttimescales=%s' % (n, lag_time, timescales))
    print()

figure(figsize=(14,3))

for i, n in enumerate(n_states):
    subplot(1,len(n_states),1+i)
    plot(lag_times, hmm_ts0[n])
    plot(lag_times, hmm_ts1[n])
    if i == 0:
        ylabel('Relaxation Timescale')
    xlabel('Lag Time')
    title('%d states' % n)

show()

msmts0, msmts1 = {}, {}

n_states = [2**n for n in range(3,11)]

for n in n_states:
    msmts0[n] = []
    msmts1[n] = []
    for lag_time in lag_times:
        assignments = KCenters(n_clusters=n).fit_predict(sequences)
        msm = MarkovStateModel(lag_time=lag_time, verbose=False).fit(assignments)
        timescales = msm.timescales_
        msmts0[n].append(timescales[0])
        msmts1[n].append(timescales[1])
        print('n_states=%d\tlag_time=%d\ttimescales=%s' % (n, lag_time, timescales[0:2]))
    print()

figure(figsize=(14,3))

for i, n in enumerate(n_states):
    subplot(1,len(n_states),1+i)
    plot(lag_times, msmts0[n])
    plot(lag_times, msmts1[n])
    if i == 0:
        ylabel('Relaxation Timescale')
    xlabel('Lag Time')
    title('%d states' % n)

show()

