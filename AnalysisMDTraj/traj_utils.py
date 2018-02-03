import mdtraj
import numpy as np
from subprocess import call, PIPE
import matplotlib.pyplot as plt


def write_cpptraj_script(traj, top, frame1=1, frame2=1, outfile=None, write=True, run=False):
    """
    Create a cpptraj script to load specific range of frames from a trajectory and write them out to a file

    :param traj: str, Location in disk of trajectories to load
    :param top: str, Location in disk of the topology file
    :param frame1: int, The first frame to load
    :param frame2: int, The last frame to load
    :param outfile: str, Name (with file format extension) of the output trajectory
    :param write: bool, Whether to write the script to a file in disk
    :param run: bool, Whether to run the script after writing it to disk
    :return cmds: str, the string representing the cpptraj script
    """
    if run and not write:
        raise ValueError('Cannot call the script without writing it to disk')
    if outfile is None:
        outfile = 'pdbs/' + traj.split('.')[0] + '.pdb'
    commands = [
        'parm {}'.format(top),
        'trajin {} {} {}'.format(traj, frame1, frame2),
        'trajout {}'.format(outfile),
        'run'
    ]
    cmds = '\n'.join(commands)
    if write:
        with open('script.cpptraj', 'w') as f:
            f.write(cmds)
        if run:
            call(['cpptraj', '-i', 'script.cpptraj'], stdout=PIPE)

    return cmds


def load_Trajs(trajfiles_list, prmtop_file, stride=1, chunk=1000):
    """
    Iteratively loads a list of NetCDF files and returns them
    as a list of mdtraj.Trajectory objects

    Parameters
    ----------
    trajfiles_list: list of str
            List with the names of trajectory files
    prmtop_file:  str
            Name of the prmtop file
    stride: int
            Frames to be used when loading the trajectories
    chunk:  int
            Number of frames to load at once from disk per iteration.
            If 0, load all.

    Returns
    -------
    list_chunks: list
            List of mdtraj.Trajectory objects, each of 'chunk' lenght
    """
    list_chunks = []
    for traj in trajfiles_list:
        for frag in mdtraj.iterload(traj, chunk=chunk, top=prmtop_file,
                                    stride=stride):
            list_chunks.append(frag)
    return(list_chunks)


def load_Trajs_generator(trajfiles_list, prmtop_file, stride, chunk):
    """
    Iteratively loads a list of NetCDF files and returns them
    as an iterable of mdtraj.Trajectory objects
    Parameters
    ----------
    trajfiles_list: list of str
            List with the names of trajectory files
    prmtop_file:  str
            Name of the prmtop file
    stride: int
            Frames to be used when loading the trajectories
    chunk:  int
            Number of frames to load at once from disk per iteration.
            If 0, load all.
    Yields
    ------
    frag: mdtraj.Trajectory
    """
    try:
        for traj in trajfiles_list:
            for frag in mdtraj.iterload(traj, chunk=chunk, top=prmtop_file,
                                        stride=stride):
                yield frag
    except OSError:
        # User passed a single long trajectory as a string
        # so there's no need to iterate through it.
        for frag in mdtraj.iterload(trajfiles_list, chunk=chunk, top=prmtop_file,
                                    stride=stride):
            yield frag


def traj_list_to_dict(trajfiles_list, prmtop_file, stride=1):
    """
    Loads a list of trajs passed as a list of strings into a
    dictionary with keys as integers from 0
    """
    trajs_dict = {}
    for i, traj in enumerate(trajfiles_list):
        trajs_dict[i] = mdtraj.load(traj, top=prmtop_file, stride=stride)
    return trajs_dict


def split_trajs_by_type(traj_dict, meta):
    """
    Find the kind of types of simulations inside the meta object
    and build a dictionary that has them as keys. Then, build a dictionary
    of the trajs inside traj_dict that belong to each type.
    """

    if len(traj_dict) != len(meta):
        raise ValueError('Lengths of traj_dict and meta do not match.')

    type_set = set(meta['type'])
    # dict which stores each subtype dict of trajs
    type_dict = dict.fromkeys(type_set)

    for t in type_set:
        new_dict = {}
        for i, row in meta.iterrows():
            if row['type'] == t:
                new_dict[i] = traj_dict[i]
        type_dict[t] = new_dict
    return type_dict


def trim_centers_by_region(clusterer, x1=None, x2=None, y1=None, y2=None, obs=(0, 1)):
    """
    Find the cluster centers that fall within a user-defined region.

    :param clusterer: an msmbuilder cluster object
    :param x1: float The low limit of the x axis
    :param x2: float The high limit of the x axis
    :param y1: float The low limit of the y axis
    :param y2: float The high limit of the y axis
    :param obs: tuple, the dimensions to sample
    :return trimmed: np.array, Cluster centers that are within the region
    """
    if not hasattr(clusterer, 'cluster_centers_'):
        raise AttributeError('The provided clusterer object has no cluster_centers_ property.')
    centers = clusterer.cluster_centers_
    pruned = centers[:, obs]
    if x1 is None:
        x1 = np.min(pruned[:, 0])
    if y1 is None:
        y1 = np.min(pruned[:, 1])
    if x2 is None:
        x2 = np.max(pruned[:, 0])
    if y2 is None:
        y2 = np.max(pruned[:, 1])

    trimmed = centers[
        ((pruned[:, 0] > x1) & (pruned[:, 0] < x2)) &
        ((pruned[:, 1] > y1) & (pruned[:, 1] < y2))
    ]
    return trimmed


def cartesian_product(x, y):
    return np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])


def generate_traj_from_stateinds(inds, meta):
    for state_i, state_inds in enumerate(inds):
        traj = mdtraj.join(
            mdtraj.load_frame(meta.loc[traj_i]['traj_fn'], index=frame_i, top=meta.loc[traj_i]['top_fn'])
            for traj_i, frame_i in state_inds
        )
    traj.center_coordinates()
    traj.superpose(traj, 0)
    return traj


def load_in_vmd(dirname, inds):
    k = len(inds[0])
    templ = [
        '# Defaults',
        'mol default material AOChalky',
        'mol default representation NewCartoon',
        'color Display {Background} white',
        'axes location off',
    ]
    for i in range(k):
        templ += [
            '# State {}'.format(i),
            'mol new {}/{:03d}.pdb'.format(dirname, i),
            'mol rename top State-{}'.format(i),
            'mol modcolor 0 top ColorID {}'.format(i),
            'mol drawframes top 0 0:{k}'.format(k=k),
            'mol modselect 0 top resid 1 to 161',
            'mol modcolor 0 top ColorID 0',
            'mol addrep top',
            'mol modselect 1 top resid 162 to 248',
            'mol modcolor 1 top ColorID 7',
            'mol addrep top',
            'mol modselect 2 top resid 249 to 419',
            'mol modcolor 2 top  ColorID 1',
            'mol addrep top',
            'mol modselect 3 top not protein and not resname CAL',
            'mol modstyle 3 top Licorice',
            'mol addrep top',
            'mol modselect 4 top resname CAL',
            'mol modstyle 4 top VDW',
            'mol modcolor 4 top ColorID 6'
            '',
        ]
    return '\n'.join(templ)


def get_source_sink(msm, clusterer, eigenvector):
    """
    Get the source and sink of a given eigenvector, in cluster naming of clusterer object
    :param msm:
    :param clusterer:
    :param eigenvector:
    :return:
    """
    source_msm_naming = np.argmin(msm.left_eigenvectors_[:, eigenvector])
    sink_msm_naming = np.argmax(msm.left_eigenvectors_[:, eigenvector])

    source_clusterer_naming = msm.state_labels_[source_msm_naming]
    sink_clusterer_naming = msm.state_labels_[sink_msm_naming]

    assert msm.mapping_[source_clusterer_naming] == source_msm_naming
    assert msm.mapping_[sink_clusterer_naming] == sink_msm_naming

    return source_clusterer_naming, sink_clusterer_naming
