#!/usr/bin/env python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.cluster import KMeans
import msmexplorer as msme
import mdtraj
import argparse
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns

parser = argparse.ArgumentParser(prog='cluster_ligands.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''version1''')

parser.add_argument("Trajectories", help="""An indefinite amount of AMBER
                    trajectories""", nargs="+")
parser.add_argument('-p', '--prmtop', type=str, required=True)
parser.add_argument('-st', '--stride', type=int, required=False, default=1)
parser.add_argument('-l', '--ligand_selection', type=str, required=True)
parser.add_argument('-o', '--out_file', type=str, required=False,
                    default='cluster_ligands')


def figure_dims(width_pt, factor=0.45):
    """
    I copied this from here:
    https://www.archer.ac.uk/training/course-material/2014/07/SciPython_Cranfield/Slides/L04_matplotlib.pdf
    """
    WIDTH = width_pt  # Figure width in pt (usually from LaTeX)
    FACTOR = factor  # Fraction of the width you'd like the figure to occupy
    widthpt = WIDTH * FACTOR
    inperpt = 1.0 / 72.27
    golden_ratio = (np.sqrt(5) - 1.0) / 2.0  # because it looks good
    widthin = widthpt * inperpt
    heightin = widthin * golden_ratio
    figdims = [widthin, heightin]  # Dimensions as list
    return figdims


def plot_com_matrix(com_matrix):
    """
    Plots the projections of X, Y and Z of a center of mass position matrix
    Parameters
    ----------
    com_matrix: np.array of shape (frames, 3)

    Returns
    -------
    fig: plt.figure
    """
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    correspondance = {
        0: 'x',
        1: 'y',
        2: 'z'
    }
    # axes is a np.array of shape (3, 3)
    for i in range(3):
        for j in range(3):
            ax = axes[i][j]
            if i >= j:
                # we'll plot only the lower half
                if i == j:
                    sns.kdeplot(com_matrix[:, i], ax=ax)
                    ax.set(xlabel='%s (nm)' % correspondance[i],
                           ylabel='Density')
                else:
                    ax.set(xlabel='%s (nm)' % correspondance[j],
                           ylabel='%s (nm)' % correspondance[i])
                    msme.plot_free_energy(
                        com_matrix,
                        obs=(j, i),
                        ax=ax,
                        shade=True,
                        n_levels=5,
                        vmin=-1e-12,
                        cmap='viridis'
                    )
            else:
                ax.set_axis_off()
    return fig


def report_clusters(com_matrix, n_cluster_list=[2, 3, 4, 5, 6]):
    """
    Performs K means clustering on a center of mass position matrix.
    Reports the results in two plots for each number of clusters:
        1st Plot: Silhouette score summary
        2nd Plot: 3D representation of the center of mass positions and the 
            cluster labels that have been assigned.
    Parameters
    ----------
    com_matrix: np.array of shape (frames, 3)

    Returns
    -------
    fig_list: list of plt.figures, of length n_cluster_list
    """
    fig_list = []
    min_silhouette_value = 0
    for n_clusters in n_cluster_list:
        # Create a figure which will have two axis (1 row, 2 columns)
        fig = plt.figure(figsize=figure_dims(2500))
        ax1 = plt.subplot(121)  # Plot on the left will be simple 2D silhouette
        ax2 = plt.subplot(122, projection='3d')  # 3D plot on the right

        clusterer = KMeans(n_clusters=n_clusters)
        cluster_labels = clusterer.fit_predict(com_matrix)
        silhouette_avg = silhouette_score(com_matrix, cluster_labels)
        sample_silhouette_values = silhouette_samples(com_matrix, cluster_labels)

        # Silhouette values can go between -1 and 1, but here let's set the minimum
        # to whatever minimal value we have
        min_silhouette_value = min(min_silhouette_value, sample_silhouette_values.min())
        # Set axis values of first plot
        ax1.set_ylim([0, len(com_matrix) + (n_clusters + 1) * 10])
        ax1.set_xlim([min_silhouette_value, 1])
        y_lower = 10
        for i in range(n_clusters):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            ith_cluster_silhouette_values = \
                sample_silhouette_values[cluster_labels == i]

            ith_cluster_silhouette_values.sort()

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.spectral(float(i) / n_clusters)
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                              0, ith_cluster_silhouette_values,
                              facecolor=color, edgecolor=color, alpha=0.7)

            # Label the silhouette plots with their cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples

        ax1.set_title("Silhouette plot")
        ax1.set_xlabel("Silhouette coefficient")
        ax1.set_ylabel("Cluster label")

        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
        ax1.annotate('Avg. silhouette: %.2f' % silhouette_avg, xy=(0.75, 0.95), xycoords='axes fraction')
        ax1.set_yticks([])  # Clear the yaxis labels / ticks

        # 2nd Plot showing the actual clusters formed
        colors = cm.spectral(cluster_labels.astype(float) / n_clusters)
        ax2.scatter(com_matrix[:, 0], com_matrix[:, 1], com_matrix[:, 2],
                    marker='.', s=10, lw=0, alpha=0.2, c=colors, zorder=-1)

        # Labeling the clusters
        centers = clusterer.cluster_centers_
        # Plot the label of the cluster at its center on top of a white round dot
        for i, c in enumerate(centers):
            color = cm.spectral(float(i) / n_clusters)
            ax2.scatter(c[0], c[1], c[2], marker='o', alpha=1, s=200, color='white', lw=1, zorder=20)
            ax2.scatter(c[0], c[1], c[2], marker='$%d$' % i, alpha=1, s=100, color=color, zorder=20)

        ax2.set_xlabel("x (nm)")
        ax2.set_ylabel("y (nm)")
        ax2.set_zlabel("z (nm)")

        fig_list.append((fig, n_clusters))
    return fig_list


def plot_3d_time(com_matrix, time):
    fig = plt.figure(figsize=figure_dims(2500))
    ax = plt.subplot(111, projection='3d')  # 3D plot on the right
    # Set background color of 3D axis to white
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_zlabel("z (nm)")
    p = ax.scatter(com_matrix[:, 0], com_matrix[:, 1], com_matrix[:, 2],
                   c=time, cmap='viridis')
    cbar = fig.colorbar(p)
    cbar.set_label('Time (ns)')
    fig.savefig('%s_3d.pdf' % args.out_file)

    return ax


if __name__ == '__main__':
    args = parser.parse_args()
    print(args)
    top = mdtraj.load_prmtop(args.prmtop)
    traj = mdtraj.load([t for t in args.Trajectories], top=args.prmtop,
                       stride=args.stride)
    # Center all coordinates, makes center of geometry of the system (0, 0, 0)
    traj.center_coordinates()
    # Superpose trajectory onto first frame
    # This is exactly like the 'hold selection steady' command in Chimera
    traj.superpose(traj, 0)
    # Create separate trajectory for the ligand
    ligand_indices = top.select(args.ligand_selection)
    lig_traj = traj.atom_slice(ligand_indices)
    center_mass_ligand = mdtraj.compute_center_of_mass(lig_traj)

    # Plot projections
    f = plot_com_matrix(center_mass_ligand)
    f.savefig(args.out_file + '.pdf')
    # Plot clustering with Kmeans
    fig_list = report_clusters(center_mass_ligand)
    for fig in fig_list:
        f, n_clusters = fig[0], fig[1]
        f.savefig(args.out_file + str(n_clusters) + 'clusters_kmeans.pdf')

    # plot 3d alone
    time = np.linspace(0, traj.timestep * traj.n_frames / 1000, num=traj.n_frames)
    plot_3d_time(center_mass_ligand, time)
