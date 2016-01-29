import mdtraj as md


def load_Trajs(traj_list, topology, stride, chunk):
    """
    Iteratively loads a list of NetCDF files and returns them
    as a list of mdtraj.Trajectory objects

    Parameters
    ----------
    traj_list: list of str
            List with the names of trajectory files
    topology:  str
            Name of the topology file
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
    for traj in traj_list:
        for frag in md.iterload(traj, chunk=chunk, top=topology,
                                stride=stride):
            list_chunks.append(frag)
    return(list_chunks)


def load_Trajs_generator(traj_list, topology, stride, chunk):
    """
    Iteratively loads a list of NetCDF files and returnns them
    as an iterable of mdtraj.Trajectory objects
    Parameters
    ----------
    traj_list: list of str
            List with the names of trajectory files
    topology:  str
            Name of the topology file
    stride: int
            Frames to be used when loading the trajectories
    chunk:  int
            Number of frames to load at once from disk per iteration.
            If 0, load all.
    Yields
    ------
    frag: mdtraj.Trajectory
    """
    for traj in traj_list:
        for frag in md.iterload(traj, chunk=chunk, top=topology,
                                stride=stride):
            yield frag
