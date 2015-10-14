#!/usr/bin/env python
"""
Bayesian Estimation of MSMs

< http://msmbuilder.org/latest/examples/bayesian-msm.html>
"""

import numpy as np
from matplotlib import pyplot as plt
plt.style.use("ggplot")
from mdtraj.utils import timing
from msmbuilder.example_datasets import load_doublewell
from msmbuilder.cluster import NDGrid
from msmbuilder.msm import BayesianMarkovStateModel, MarkovStateModel


trjs = load_doublewell(random_state=0)['trajectories']
plt.hist(np.concatenate(trjs), bins=50, log=True)
plt.ylabel('Frequency')
plt.show()

