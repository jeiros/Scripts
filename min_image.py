#!/usr/bin/env python
from __future__ import print_function
import mdtraj as md
import time
import numpy as np
from copy import deepcopy
import itertools
import sys


def main():
    TOP = sys.argv[1]
    STRIDE = 10
    TRAJS = sys.argv[2:]

    print('Topology is %s' % TOP)
    print('Trajs are {}'.format(list(TRAJS)))

    t0 = time.time()
    T = md.load(TRAJS, top=TOP, stride=STRIDE)
    # Select alpha carbons
    T_gp = T.atom_slice(T.top.select("name CA"))
    # For use at the very end.
    T_gp0 = deepcopy(T_gp)
    print("%d frames were loaded using a stride of %d (%.1f ns)" % (T.n_frames, STRIDE, ((T.n_frames * T.timestep) / 1000)))
    print("Finished in %.1f seconds" % (time.time() - t0))

    t0 = time.time()
    print("Making images...")
    T_images = []
    vecs = []
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                if [i, j, k] == [0, 0, 0]:
                    continue
                vecs.append([i, j, k])
                T_copy = deepcopy(T_gp)
                displace = T.unitcell_vectors[:, 0, :] * i + T.unitcell_vectors[:, 1, :] * j + T.unitcell_vectors[:, 2, :] * k
                T_copy.xyz += displace[:, np.newaxis, :]
                T_images.append(T_copy)
    print("Done in %.3f seconds" % (time.time() - t0))

    t0 = time.time()
    print("Stacking images...")
    n = T_gp.n_atoms
    for T_i in T_images:
        T_gp = T_gp.stack(T_i)
    T_gp.unitcell_vectors *= 3
    print("Done in %.3f seconds" % (time.time() - t0))

    t0 = time.time()
    print("Measuring distances... (time consuming)")
    min_frames = []
    min_pairs = []
    min_dists = []
    for i in range(26):
        print("Working on lattice vector", vecs[i])
        sel1 = range(n)
        sel2 = list(np.arange(n) + (i + 1) * n)
        image_pairs = list(itertools.product(sel1, sel2))
        d = md.compute_distances(T_gp, image_pairs, periodic=False, opt=True)
        frame, pair = np.unravel_index(np.argmin(d), d.shape)
        min_frames.append(frame)
        min_pairs.append(image_pairs[pair])
        min_dists.append(d[frame, pair])
    print("Done in %.3f seconds" % (time.time() - t0))

    min_image = np.argmin(min_dists)
    i, j, k = vecs[min_image]
    T_copy = deepcopy(T)
    displace = T.unitcell_vectors[:, 0, :] * i + T.unitcell_vectors[:, 1, :] * j + T.unitcell_vectors[:, 2, :] * k
    T_copy.xyz += displace[:, np.newaxis, :]
    T = T.stack(T_copy)
    original_frame = min_frames[min_image] * 10 + 1
    print("Closest contact found at frame %i, lattice vector %s, atom pair %i-%i, distance %.3f Angstrom" % (original_frame, str(vecs[min_image]), min_pairs[min_image][0] % n, min_pairs[min_image][1] % n, 10 * min_dists[min_image]))
    pdb_filename = "min_image_f%i_%.1fA_%s.pdb" % (original_frame, 10 * min_dists[min_image], TRAJS[0][:-3])
    print('Saving %s' % pdb_filename)
    T[min_frames[min_image]].save_pdb(pdb_filename)

if __name__ == '__main__':
    try:
        main()
    except:
        sys.exit(0)
