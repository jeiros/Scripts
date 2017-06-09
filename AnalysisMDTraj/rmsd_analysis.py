import mdtraj


def rmsd(filenames, topology):
    rmsds = []
    first_frame = mdtraj.load_frame(filenames[0], 0, top=topology)
    for fragment in filenames:
        for chunk in mdtraj.iterload(fragment, chunk=100, top=topology):
            rmsds.append(mdtraj.rmsd(chunk, first_frame))
    return np.concatenate(rmsds)


def pca(list_trajs):
    topology = list_trajs[0].topology
    ca_backbone = topology.select("backbone and name CA")
    pair_distances = []
    for chunk in list_trajs:
        pairs = topology.select_pairs(ca_backbone, ca_backbone)
        X = mdtraj.compute_distances(chunk, pairs)
        pair_distances.append(X)
    return np.concatenate(pair_distances)
