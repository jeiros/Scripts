import mdtraj


def load_Trajs(trajfiles_list, prmtop_file, stride, chunk):
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
    for traj in trajfiles_list:
        for frag in mdtraj.iterload(traj, chunk=chunk, top=prmtop_file,
                                    stride=stride):
            yield frag
