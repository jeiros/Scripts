import argparse
import sys


class TrajectoryParser(object):
    """   """
    def __init__(self, trajectory_list, topology_file):
        self.trajectory_list = trajectory_list
        self.topology_file = topology_file

    def return_values():
        parser = argparse.ArgumentParser(usage="""{} Trajectories*.nc Topology.prmtop""".
                                 format(sys.argv[0]),
                                 epilog="""Load up a list of AMBER NetCDF 
                                 trajectories and their corresponding topology
                                 with MDtraj.""")
        args = parser.parse_args()
        print(args)
        return(args)


# parser = argparse.ArgumentParser(usage="""{} Trajectories*.nc Topology.prmtop""".
#                                  format(sys.argv[0]),
#                                  epilog="""Load up a list of AMBER NetCDF 
#                                  trajectories and their corresponding topology
#                                  with MDtraj. Calculates contact maps""")

# parser.add_argument("Trajectories", help="""An indefinite amount of AMBER
#                     trajectories""", nargs="+")

# parser.add_argument("Topology", help="""The topology .prmtop file that matches
#                     the trajectories""")
# parser.add_argument("-st", "--stride", help="""Stride for the loading of the
#                     trajectory. Must be a divisor of the chunk.
#                     Default value is 1.""", default=1, type=int)

# parser.add_argument("-ch", "--chunk", help="""Number of frames that will be 
#                     used by md.iterload to load up the trajectories. Must be
#                     a multiplier of the stride.
#                     Default is 100 frames.""", default=100, type=int)
# args = parser.parse_args()
