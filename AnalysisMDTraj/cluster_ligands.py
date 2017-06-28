#!/usr/bin/env python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from sklearn.metrics import silhouette_score, silhouette_samples
import sklearn.cluster
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
parser.add_argument('-l', '--ligand_selection', type=str, required=False, default='all')
parser.add_argument('-o', '--out_file', type=str, required=False,
                    default='cluster_ligands')


def plot_com_matrix(com_matrix):
    f, axes = plt.subplots(3, 3, figsize=(15, 15))
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
                    ax.set(xlabel='%s (Å)' % correspondance[i],
                           ylabel='Density')
                else:
                    ax.set(xlabel='%s (Å)' % correspondance[j],
                           ylabel='%s (Å)' % correspondance[i])
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
    return f


def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)] * a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))


def report_clusters(X):
    fig_list = []
    avg_sils = []
    for n_clusters in [2, 3, 4, 5, 6]:
        fig = plt.figure()
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122, projection='3d')
        ax2.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax2.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax2.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        fig.set_size_inches(18, 7)
        ax1.set_xlim([-0.2, 1])
        ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])
        clusterer = sklearn.cluster.KMeans(n_clusters=n_clusters)
        cluster_labels = clusterer.fit_predict(X)
        silhouette_avg = silhouette_score(X, cluster_labels)
        avg_sils.append(silhouette_avg)
        sample_silhouette_values = silhouette_samples(X, cluster_labels)
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
        ax1.annotate('Average silhouette: %.2f' % silhouette_avg, xy=(0.80, 0.95), xycoords='axes fraction')
        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-0.2, -0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

        # 2nd Plot showing the actual clusters formed
        colors = cm.spectral(cluster_labels.astype(float) / n_clusters)
        ax2.scatter(X[:, 0], X[:, 1], X[:, 2], marker='.', s=30, lw=0, alpha=0.1,
                    c=colors)

        # Labeling the clusters
        centers = clusterer.cluster_centers_

        unique_colors = unique_rows(colors)
        for i, c in enumerate(zip(centers, unique_colors)):
            c, col = c[0], c[1]
            # Draw white circles at cluster centers
            ax2.scatter(c[0], c[1], c[2], marker='o', alpha=1, s=150, color='white')
            ax2.scatter(c[0], c[1], c[2], marker='$%d$' % i, alpha=1, s=100, color=col)

        ax2.set_xlabel("x (nm)")
        ax2.set_ylabel("y (nm)")
        ax2.set_zlabel("y (nm)")

        fig_list.append((fig, n_clusters))
    return fig_list


if __name__ == '__main__':
    args = parser.parse_args()
    print(args)
    top = mdtraj.load_prmtop(args.prmtop)
    traj = mdtraj.load([t for t in args.Trajectories], top=args.prmtop,
                       stride=args.stride)
    # Center all coordinates, makes center of geometry of the system (0, 0, 0)
    # traj.center_coordinates()
    # Superpose trajectory onto first frame
    # This is exactly like the 'hold selection steady' command in Chimera
    traj.superpose(traj, 0)
    # Create separate trajectory for the ligand
    ligand_indices = top.select(args.ligand_selection)
    lig_traj = traj.atom_slice(ligand_indices)
    center_mass_ligand = mdtraj.compute_center_of_mass(lig_traj)
    time = np.linspace(0, traj.timestep * traj.n_frames, num=traj.n_frames)
    f = plot_com_matrix(center_mass_ligand)
    f.savefig(args.out_file + '.pdf')
    fig_list = report_clusters(center_mass_ligand)
    for fig in fig_list:
        f, n_clusters = fig[0], fig[1]
        f.savefig(args.out_file + str(n_clusters) + '.pdf')
