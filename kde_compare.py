#!/usr/bin/env python

import numpy as np
from scipy import stats
import sys

p_dist = np.loadtxt(sys.argv[1], usecols=(1,))
wt_dist = np.loadtxt(sys.argv[2], usecols=(1,))


def integral(dist):
    if np.isnan(dist).any():
        dist = np.nan_to_num(dist)
    kde = stats.gaussian_kde(dist)
    integral = kde.integrate_box(4, 5.5)
    return(integral)

print(integral(p_dist)/integral(wt_dist))


import numpy as np
from scipy import stats
import sys
import matplotlib.pylab as plt


def filter_lines(f, stride):
    for i, line in enumerate(f):
        if i % stride and (i-1) % stride:
            yield line

with open("S69_distances.dat") as f:
    dist = np.loadtxt(filter_lines(f, 3),
                      usecols=(1,))

kde = stats.gaussian_kde(dist)
X = np.linspace(1, dist.shape[0], dist.shape[0])
plt.figure()
kdepdf = kde.evaluate(X)
plt.hist(dist, bins=200, normed=1)
plt.plot(X, kdepdf, label='kde', color='g')
plt.show()
