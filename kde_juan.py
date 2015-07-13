#!/usr/bin/env python
"""
Plots the PC1 vs PC2 projection using a kernel density estimator to smooth the
scatter plot
"""
import numpy as np
import pylab as pl
from scipy import stats
from scipy import interpolate
import matplotlib.pyplot as plt
import sys
import argparse


parser = argparse.ArgumentParser(usage="{} file1 [file2 [file3 [... ] ] ]".
                                 format(sys.argv[0]), epilog="Plot a series of\
                                 KDEplots for every evecs file passed as \
                                 argument. Forces the axis to be the same for \
                                 all of them for easier comparison")

parser.add_argument("files", help="An indefinite amount of files to plot",
                    nargs="+")

parser.add_argument("-s", "--save", help="Save the plots as .png images",
                    action="store_true")
parser.add_argument("-d", "--dots", help="Add the plots of the distribution \
                    to the plot.", action="store_true")
args = parser.parse_args()


def parse_data(files):
    """Returns a list called data that stores the X and Y values for each file.
    Find the minimum and maximum values of X and Y between all the files in
    order to always have the same edges"""
    data = []
    minx = 0
    miny = 0
    maxx = 0
    maxy = 0
    for file in args.files:
        f_colx = np.loadtxt(file, usecols=(1,), skiprows=1)
        f_coly = np.loadtxt(file, usecols=(2,), skiprows=1)
        minx = min(f_colx.min(), minx)
        maxx = max(f_colx.max(), maxx)
        miny = min(f_coly.min(), miny)
        maxy = max(f_coly.max(), maxy)
        data.append([f_colx, f_coly])

    # Added +-10 to the edges so the KDE is better seen (if not sometimes its
    # cut for the data that is in the limits of the plot)
    minx = minx - 10
    maxx = maxx + 10
    miny = miny - 10
    maxy = maxy + 10
    return (data, minx, maxx, miny, maxy)


def generate_kde(data):
    """ Returns a list called kernels with the kernel density estimator for
    each value stored in the data list"""
    kernels = []
    for i in range(0, len(args.files)):
        values = np.vstack((data[i][0], data[i][1]))
        kernels.append(stats.gaussian_kde(values))
    return kernels


def generate_shape(kernels, minx, maxx, miny, maxy):
    """Returns a list called shapes that stores the needed values to plot the
    shape of the KDE
    This is actually copied from here
    <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/
    scipy.stats.gaussian_kde.html>"""
    shapes = []
    for i in range(0, len(kernels)):
        X, Y = np.mgrid[minx:maxx:100j, miny:maxy:100j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        shape = X.shape
        sh = np.reshape(kernels[i](positions).T, shape)
        shapes.append(sh)
    return shapes


def generate_title_from_file(title):
    """This function removes the "myevecs" string and ".dat" from the files to
    create the title. The input files should look like
    "myevcs000-050ns_run1.dat" and the title that is created is 000-050ns_run
    """
    return title.replace("myevecs", "").replace(".dat", "")


def plot(shapes, data, minx, maxx, miny, maxy):
    figs = []
    for i in range(0, len(shapes)):
        figs.append(plt.figure())
        ax = figs[i].add_subplot(111)
        figs[i].suptitle(generate_title_from_file(sys.argv[i + 1]),
                         fontsize=14)
        ax.imshow(np.rot90(shapes[i]), cmap=plt.cm.gist_earth_r,
                  extent=[minx, maxx, miny, maxy])
        ax.set_xlabel('PCA1')
        ax.set_ylabel('PCA2')
        ax.set_xlim([minx, maxx])
        ax.set_ylim([miny, maxy])
        if args.dots:
            ax.plot(data[i][0], data[i][1], 'k.', markersize=2)
        if args.save:
            figs[i].savefig(sys.argv[i + 1].replace(".dat", "") + ".png")
        else:
            plt.show()


def main():
    print(args)
    if args.files:
        data, minx, maxx, miny, maxy = parse_data(args.files)
        kernels = generate_kde(data)
        shapes = generate_shape(kernels, minx, maxx, miny, maxy)
        plot(shapes, data, minx, maxx, miny, maxy)

if __name__ == "__main__":
    main()
