#!/usr/bin/env python
"""
HMM and MSM Timescales for Ala2

This example builds HMM and MSMs on the alanine_dipeptide dataset using varing lag times and numbers of states, and compares the relaxation timescales

Copied from <http://msmbuilder.org/latest/examples/hmm-and-msm.html>

"""

import os
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
plt.style.use("ggplot")
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~        HIDDEN MARKOV MODEL     ~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lag_times = [1, 10, 20, 30, 40]
hmm_ts0 = {}
hmm_ts1 = {}
n_states = [3, 5]

for n in n_states:
    hmm_ts0[n] = []
    hmm_ts1[n] = []
    for lag_time in lag_times:
        strided_data = [s[i::lag_time] for s in sequences for i in range(lag_time)]
        hmm = GaussianHMM(n_states=n, n_init=1).fit(strided_data)
        timescales = hmm.timescales_ * lag_time
        hmm_ts0[n].append(timescales[0])
        hmm_ts1[n].append(timescales[1])
        print('n_states=%d\tlag_time=%d\ttimescales=%s' % (n, lag_time, timescales))

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~         MARKOV STATE MODEL     ~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

msmts0, msmts1 = {}, {}
lag_times = [1, 10, 20, 30, 40]
n_states = [4, 8, 16, 32, 64]

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