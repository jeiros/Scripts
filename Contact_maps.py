
# coding: utf-8

# # Contact maps
# In this notebook we calculate the contact maps between two masks.

# In[1]:

from traj_loading import load_Trajs_generator
import numpy as np
import mdtraj as md
from glob import glob
import itertools


# Provide the topologies with glob expression.
# 
# 
# Provide the matching topology as a string.

# In[2]:

nc_file_list = sorted(glob(
        "/Users/je714/Troponin/IAN_Troponin/completehowarthcut/salted/ff14SB/run1/05_Prod_WTff14SB_000-050ns_run1.nc"))
topology = "/Users/je714/Troponin/IAN_Troponin/completehowarthcut/salted/ff14SB/run1/WT-ff14SB_clean.prmtop"


# Provide two masks to check the contacts. 0 indexed residue ranges only.
# 
# We generate a residue-residue pair list with `itertools.product()`.

# In[3]:

start1=0
end1=2
start2=3
end2=5
list1 = list(range(start1, end1 + 1))
list2 = list(range(start2, end2 + 1))
pairs = list(itertools.product(list1, list2))
print('Mask1 has %d residues\n' % (end1 - start1 + 1))
print('Mask2 has %d residues\n' % (end2 - start2 + 1))
print('Number of residue-residue interactions that will be considered: %d\n' % len(pairs))
print('Mask1:%s\nMask2:%s\nResidue pair list:%s\n' % (list1, list2, pairs))


# ## Contact implementation using `md.compute_contacts`.
# 
# It only allows for three types of contact definition. It gives much higher values than cpptraj unless really low values of distance are used as a cutoff.

# In[4]:

stride = 5
chunk = 100
trajs_generator= load_Trajs_generator(trajfiles_list=nc_file_list,
                                      prmtop_file=topology, stride=stride, chunk=chunk)
frequency = np.zeros(len(pairs))
count = 0 # Keep the frame count
for traj_chunk in trajs_generator:
    count += traj_chunk.n_frames
    distances_inChunk = md.compute_contacts(traj_chunk, pairs, scheme = 'closest-heavy')
    column_sum = (distances_inChunk[0] <= 4.6).sum(0).T # Sum by colum and transpose
    frequency += column_sum # Sum the partial result to frequency
contact_frequency = frequency / count # Total contact value for each residue-residue pair. From 0 to 1.
print('Number of analyzed frames: %d\n' % count)
print('Aggregate simulation time: %2.f ns\n' % (count * 0.02 * stride))
print(contact_frequency)


# ## Manual implementation of contacts
# Let's try to define a contact as Cheng et al. paper.
# 
# 
# A carbon-carbon distance of <= 5.4 Å and a distance between any other noncarbon atoms of <= 4.6 Å is defined as a *contact*.
# 
# We begin by adding a Cartesian product implementation for two numpy arrays. Extracted from [StackOverflow](http://stackoverflow.com/questions/11144513/numpy-cartesian-product-of-x-and-y-array-points-into-single-array-of-2d-points). We'll use it when computing the pairwise atoms to check the distance for in each residue-residue pair.

# In[5]:

def cartesianProduct_npArrays(x,y):
    return(np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))]))


# In[34]:

top = md.load_prmtop(topology)
stride = 5
chunk = 100
trajs_generator= load_Trajs_generator(trajfiles_list=nc_file_list,
                                      prmtop_file=topology, stride=stride, chunk=chunk)
frequency = np.zeros(len(pairs))
print(frequency)
count = 0 # Keep the frame count
for traj_chunk in trajs_generator:
    count += traj_chunk.n_frames
    for residue_pair in pairs:
        print(residue_pair)
        c_atoms_residue1 = top.select("resid %d and (type C)" % residue_pair[0])
        not_c_atoms_residue1 = top.select("resid %d and not type C" % residue_pair[0])
        c_atoms_residue2 = top.select("resid %d and type C" % residue_pair[1])
        not_c_atoms_residue2 = top.select("resid %d and not type C" % residue_pair[1])
        
        c_atoms_dist = md.compute_distances(traj_chunk,
                                            cartesianProduct_npArrays(c_atoms_residue1,
                                                                      c_atoms_residue2))
        not_c_atoms_dist = md.compute_distances(traj_chunk,
                                                cartesianProduct_npArrays(not_c_atoms_residue1,
                                                                          not_c_atoms_residue2))
        


# In[8]:

traj1 = md.load_netcdf(
    "/Users/je714/Troponin/IAN_Troponin/completehowarthcut/salted/ff14SB/run1/05_Prod_WTff14SB_000-050ns_run1.nc",
    top = top)


# In[10]:

c_atom_pairs = cartesianProduct_npArrays(top.select("resid 0 and type C"),
                                         top.select("resid 3 and type C"))
not_c_atom_pairs = cartesianProduct_npArrays(top.select("resid 0 and not type C"),
                                             top.select("resid 3 and not type C"))
c_distances = md.compute_distances(traj1, c_atom_pairs)
not_c_distances = md.compute_distances(traj1, not_c_atom_pairs)

