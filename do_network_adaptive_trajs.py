#!/usr/bin/env python


import argparse
from adaptive_sampling_utils import network_from_simulations
import networkx as nx
from matplotlib import pyplot as plt
import seaborn as sns
from plot_utils import figure_dims
from fa2 import ForceAtlas2
from itertools import count
sns.set_style('ticks')
parser = argparse.ArgumentParser(prog='do_network_adaptive_trajs.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''
____________________________________________________________________________
| A program that generates a network from a set of folders following an    |
| adaptive sampling scheme                                                 |
----------------------------------------------------------------------------
''')

parser.add_argument('folders', type=str, help='A glob expression')
parser.add_argument('-t', '--title', type=str, required=False, default=None)
parser.add_argument('-o', '--outputfile', type=str, required=False, default='network_layout')
parser.add_argument('-l', '--labels', type=bool, required=False, default=False)


def main(args):
    G = network_from_simulations(args.folders, save=True)
    frc = ForceAtlas2(
        outboundAttractionDistribution=True,
        gravity=1,
    )
    positions = frc.forceatlas2_networkx_layout(G)
    degree = nx.degree(G)
    nodes = G.nodes()
    base_size = 100

    sizes = [
        (degree[node] + 1) * base_size for node in nodes
    ]

    groups = set(nx.get_node_attributes(G, 'epoch').values())
    mapping = dict(zip(sorted(groups), count()))
    colors = [mapping[G.node[n]['epoch']] for n in nodes]
    cmap = plt.cm.get_cmap('tab20c', len(groups))

    f, ax = plt.subplots(figsize=figure_dims(600, 0.9))

    ec = nx.draw_networkx_edges(G, positions, alpha=0.2, ax=ax)
    nc = nx.draw_networkx_nodes(G, positions, nodelist=nodes,
                                node_color=colors,
                                cmap=cmap, ax=ax, node_size=sizes)
    if args.labels:
        nx.draw_networkx_labels(G, positions, font_size=10, ax=ax)
    cbar = plt.colorbar(nc)
    if args.title is not None:
        ax.set_title(args.title)
    cbar.ax.set_ylabel('Epoch')
    plt.axis('off')
    if not args.outputfile.endswith('.pdf'):
        args.outputfile += '.pdf'
    f.savefig(args.outputfile)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
