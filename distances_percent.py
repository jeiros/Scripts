#!/usr/bin/env python

# Run the script as ./script.py distances.dat
# Produces a file distances.png


import numpy as np
import sys
import matplotlib.pyplot as plt
import os

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


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def scatter_plot(array):
    x = array[0]
    y = array[1]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Distance to Ca (Å)')
    ax.set_ylabel('Cumulative %')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 100)
    ax.grid(True)
    ticklines = ax.get_xticklines() + ax.get_yticklines()
    gridlines = ax.get_xgridlines() + ax.get_ygridlines()
    for line in ticklines:
        line.set_linewidth(1)
    for line in gridlines:
        line.set_linestyle('-.')
    ax.scatter(x, y)
    plot_filename = str(sys.argv[1]).replace(".dat", "") + ".eps"
    if os.path.isfile(plot_filename):
        var = query_yes_no("File %s exists. Do you want to overwrite?:" % plot_filename, default=None)
        if not var:
            sys.exit("Do not overwrite. Exiting script...\n")
    fig.savefig(plot_filename, dpi=900)

scatter_plot(print_res(loaded_data))
