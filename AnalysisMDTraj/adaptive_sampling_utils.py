from glob import glob
import networkx as nx
from natsort import natsorted
import re


class Simulation:
    """
    Class to retrieve data from a simulation named according to htmd adaptive's sampling
    """

    def __init__(self, name):
        """
        Initializer
        :param name: str
            Examples:
                e1s1_gen1
                e3s3_e1s1p0f88 <-  This means sim 3 of epoch3 spawned from the part 0, frame 88 of e1s1
        """
        self.name = name

    def _match(self):
        """
        Matches any occurence of the string e*s* where * can be any number of digits
        :return: list
        """
        return re.findall(r'(e\d*s\d)', self.name)

    @property
    def node(self):
        """
        The simulation itself
        :return: str
        """
        return self._match()[0]

    @property
    def parent(self):
        """
        Where the sim comes from (whatever comes after the first _ char)
        :return: str
        """
        matches = self._match()
        if len(matches) == 1:
            # Sim had no parent so it's epoch 1
            return None
        else:
            return matches[1]

    @property
    def epoch(self):
        """
        Find the first occurence of a digit, which is the epoch
        :return: int, the epoch this sim belongs to
        """
        return int(re.findall(r'[0-9]*', self.node)[1])


def network_from_simulations(fnames, save=False, gexf_name='network.xml'):
    """
    Build a directed graph from the structure of any given number of simulation files named according to HTMD conventions
    :param fnames: str, a glob expression of trajectory files
    :param save: bool, opt (default: False) whether to save the graph as a file loadable with Cytoscape
    :param gexf_name: str, opt (default: 'network.xml') filename of the graph file
    :return G: networkx.classes.digraph.DiGraph
    """
    fnames = natsorted(glob(fnames))
    clean_fnames = ['_'.join(x.split('/')[-1].split('_')[0:2])
                    for x in fnames]
    G = nx.DiGraph()
    for sim in clean_fnames:
        sim = Simulation(sim)
        G.add_node(sim.node, epoch=sim.epoch)
        if sim.parent is not None:
            G.add_edge(sim.parent, sim.node)

    if save:
        if not gexf_name.endswith('.xml'):
            gexf_name += '.xml'
        nx.write_graphml(G, path=gexf_name)
    return G
