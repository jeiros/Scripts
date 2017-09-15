import mdtraj


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


def load_Trajs_generator(trajfiles_list, prmtop_file, stride=1, chunk=1000, verbose=False):
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
    for traj in trajfiles_list:
        if verbose:
            print("Loading {}".format(traj))
        for frag in mdtraj.iterload(traj, chunk=chunk, top=prmtop_file,
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
