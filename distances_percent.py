#!/usr/bin/env python

# Run the script as ./script.py distances.dat 
# Produces a file distances.png


import numpy as np
import sys
import matplotlib.pyplot as plt

loaded_data = np.loadtxt(str(sys.argv[1]))


def percent(data, threshold):
    count = 0
    for i in range(0, data.shape[0]):
        if data[i][1] < threshold:
            count += 1
    return(count*100/data.shape[0])


def np_builder(thresholds_array, data_array):
    for i in thresholds_array:
        yield percent(data_array, i)


def print_res(data):
    nums = np.linspace(0, 10, 100)
    percent_array = np.fromiter(np_builder(nums, data), dtype=np.float)
    dist_percent_array = np.vstack((nums, percent_array))
    return(dist_percent_array)


def scatter_plot(array):
    x = array[0]
    y = array[1]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Distance to Ca (Å)')
    ax.set_ylabel('Cumulative %')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 100)
    ax.scatter(x, y)
    fig.savefig(str(sys.argv[1]).replace(".dat", "") + ".png", dpi=900)

scatter_plot(print_res(loaded_data))
