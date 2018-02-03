"""
Functions for the contacts script
"""

import numpy as np
import mdtraj as md
import itertools
import seaborn as sns
from matplotlib import pyplot as plt
from tqdm import tqdm
import pandas as pd
import datetime


class Region:

    def __init__(self, mask1, mask2, pair_list, name):
        self.dict = {
            'name': name,
            'mask1': mask1,
            'mask2': mask2,
            'pairs': pair_list
        }

    def __repr__(self):
        return """{name}
------
mask1: {mask1}
mask2: {mask2}
n_pairs : {n_pairs}
        """.format(
            name=self.dict['name'],
            mask1=self.dict['mask1'],
            mask2=self.dict['mask2'],
            n_pairs=len(self.dict['pairs'])
        )


def get_residuepairs(start1, end1, start2, end2):
    """
    Takes the beggining and end of two masks (0-indexed) and returns them as
    two lists, as well as the cartesian product between them.
    Parameters
    ----------
    start1: int
        Beggining of first mask
    end1: int
        End of first mask
    start2: int
        Beggining of second mask
    end2: int
        End of second mask
    Returns
    -------
    mask1: list of int
    mask2: list of int
    pairs: list of tuples
    """
    mask1 = list(range(start1, end1 + 1))
    mask2 = list(range(start2, end2 + 1))
    pairs = list(itertools.product(mask1, mask2))
    print('Mask1 has %d residues\n' % len(mask1))
    print('Mask2 has %d residues\n' % len(mask2))
    print('%d residue-residue interactions  will be considered\n' % len(pairs))
    return(mask1, mask2, pairs)


def cmap_MDtraj(traj_generator, mask1, mask2, pairs, distance_cutoff=7,
                scheme='ca'):
    """
    Returns a contact map between two masks. Contacts between two residues
    range between 0 (if never present) and 1 (if present in all frames).
    Parameters
    ----------
    traj_generator: generator
        MD trajectory as obtained through the load_Trajs_generator() function
    mask1, mask2, pairs:
        Two masks and the corresponding residue-residue list of tuples
        as obtained from the get_residuepairs() function
    distance_cutoff: int
        Distance value (in Å) to meet the criteria of contact.
    scheme: str
        Any of the three schemes allowed by md.compute_contacts
        ['ca', 'closest', 'closest-heavy']
    Returns
    -------
    contact_frequency: np.array of shape (len(mask1), len(mask2))
        Contact map between mask1 and mask2. Can be used directly as
        input by the sns.heatmap() function.
    """

    frequency = np.zeros((len(mask1), len(mask2)))  # store the partial sum of contacts
    frame_count = 0
    for traj_chunk in traj_generator:
        frame_count += traj_chunk.n_frames
        distances_in_chunk = md.compute_contacts(traj_chunk, pairs,
                                                 scheme=scheme)
        # Sum by column
        # Divide the cutoff by 10 as mdtraj uses nm instead of Å
        column_sum = (distances_in_chunk[0] <= distance_cutoff / 10).sum(0)
        # Sum the partial result to frequency
        frequency += column_sum.reshape(len(mask1), len(mask2))
    contact_frequency = frequency / frame_count
    # Total contact value for each residue-residue pair. From 0 to 1.
    return contact_frequency


def cmap_Cheng(traj_generator, mask1, mask2, pairs, topology):
    """
    Calculate contact maps between mask1 and mask2 as described in [1].

    Parameters
    ----------
    traj_generator: generator
        MD trajectory as obtained through the load_Trajs_generator() function
    mask1, mask2, pairs:
        Two masks and the corresponding residue-residue list of tuples
        as obtained from the get_residuepairs() function

    Returns
    -------
    contact_frequency: np.array of shape (len(mask1), len(mask2))
        Contact map between mask1 and mask2. Can be used directly as
        input by the sns.heatmap() function.
    References
    ----------
    [1] Cheng, Y. et al., 2014. Biophysical Journal, 107(7), pp.1675–1685.
    """
    print("Starting the cmap_Cheng calculation...")
    top = md.load_prmtop(topology)
    frequency = np.zeros(len(pairs))
    frame_count = 0
    for traj_chunk in tqdm(traj_generator):
        frame_count += traj_chunk.n_frames
        index = 0  # To iterate through the residue-residue pair list
        for residue_pair in pairs:
            # Atom selection for each residue in the pair
            c_atoms_residue1 = top.select("resid %d and (type C)" %
                                          residue_pair[0])
            c_atoms_residue2 = top.select("resid %d and type C" %
                                          residue_pair[1])

            not_c_atoms_residue1 = top.select("resid %d and not type C" %
                                              residue_pair[0])
            not_c_atoms_residue2 = top.select("resid %d and not type C" %
                                              residue_pair[1])
            # Calculate all the possible distances between the C-C atoms and
            # the non C-C atoms. Results are stored in two np.arrays of shape:
            # (traj_chunk.n_frames, c_atoms_residue1*c_atoms_residue2)
            # (traj_chunk.n_frames, non_c_atoms_residue1*non_c_atoms_residue2)
            c_atoms_dist = md.compute_distances(traj_chunk,
                                                cartesianProduct(c_atoms_residue1,
                                                                 c_atoms_residue2))
            not_c_atoms_dist = md.compute_distances(traj_chunk,
                                                    cartesianProduct(not_c_atoms_residue1,
                                                                     not_c_atoms_residue2))
            # Implementation of the contact condition
            if (((c_atoms_dist <= 5.4).sum(1).any() > 0) and
                    ((not_c_atoms_dist <= 4.6).sum(1).any() > 0)):
                frequency[index] += 1
            index += 1
    contact_frequency = (frequency / frame_count).reshape(len(mask1), len(mask2))
    print('Number of analyzed frames: %d\n' % frame_count)
    print('Aggregate simulation time: %2.f ns\n' % (frame_count * 0.02 * args.stride))
    return(contact_frequency)


def renumber_mask(mask):
    """
    Renumbers the mask to match the correct sequence of each
    subunit. Also takes into account the 0-indexed lists of
    Python.
    """
    if max(mask) > 160 and max(mask) < 249:
            # It's in cTnT
        new_mask = [x + 51 for x in mask]
    elif max(mask) >= 249:
        # It's in cTnI
        new_mask = [x - 247 for x in mask]
    else:
        # It's in cTnC
        new_mask = [x + 1 for x in mask]
    return(new_mask)


def plot_heatmap(contact_array, mask1, mask2, title=None, save=True,
                 x_label=None, y_label=None,
                 x_steps=True, y_steps=True,
                 min_value=0, max_value=1,
                 std_array=None):
    """
    Plot a single heatmap on a red color scale.
    """

    if std_array is not None:
        # Make figure bigger as the annotation takes a lot of space
        fig, ax = plt.subplots(figsize=(20, 20))
    else:
        fig, ax = plt.subplots(figsize=(10, 10))
    # Check what part of cTn mask1 is in and renumber accordingly
    new_mask1 = renumber_mask(mask1)
    new_mask2 = renumber_mask(mask2)

    # Convert to a pd.DataFrame so sns.heatmap() can read it's
    # row/column labels properly
    contact_df = pd.DataFrame(contact_array.T, index=new_mask2,
                              columns=new_mask1)

    if std_array is not None:
        ax = sns.heatmap(contact_df,
                         vmin=min_value,
                         vmax=max_value,
                         xticklabels=x_steps,
                         yticklabels=y_steps,
                         cmap='Reds',
                         linewidths=.5,
                         annot=std_array,
                         fmt='.2f',
                         annot_kws={'rotation': 'vertical'})
    else:
        ax = sns.heatmap(contact_df,
                         vmin=min_value,
                         vmax=max_value,
                         xticklabels=x_steps,
                         yticklabels=y_steps,
                         cmap='Reds',
                         linewidths=.5)

    plt.gca().invert_yaxis()

    plt.xticks(rotation='vertical')
    plt.yticks(rotation='horizontal')

    if title is not None:
        plt.title(title)
    if x_label is not None:
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)
    if save:
        if title is None:
            filename = datetime.date.today.isoformat()
        else:
            filename = title.replace(" ", "")

        filename += ".pdf"
        plt.savefig(filename, format='pdf')
    else:
        return ax


def plot_diffmap(contact_array, mask1, mask2, title=None, save=True,
                 x_label=None, y_label=None,
                 x_steps=True, y_steps=True,
                 std_array=None):
    """
    Plots a difference heatmap with a blue-grey-red diverging
    color scale.
    """
    if std_array is not None:
        # Make figure bigger as the annotation takes a lot of space
        fig, ax = plt.subplots(figsize=(20, 20))
    else:
        fig, ax = plt.subplots(figsize=(10, 10))

    # Check what part of cTn mask1 is in and renumber accordingly
    new_mask1 = renumber_mask(mask1)
    new_mask2 = renumber_mask(mask2)

    contact_df = pd.DataFrame(contact_array.T, index=new_mask2,
                              columns=new_mask1)

    # Diverging palette with light colour on the middle
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    if std_array is not None:
        ax = sns.heatmap(contact_df,
                         vmin=contact_array.min(),
                         vmax=contact_array.max(),
                         xticklabels=x_steps,
                         yticklabels=y_steps,
                         annot=std_array,
                         mask=std_array <= 0.01,  # Only values GE 0.01 will be displayed
                         fmt='.2f',
                         annot_kws={'rotation': 'vertical'},
                         cmap=cmap,
                         linewidths=.5)
    else:
        ax = sns.heatmap(contact_df,
                         vmin=contact_array.min(),
                         vmax=contact_array.max(),
                         xticklabels=x_steps,
                         yticklabels=y_steps,
                         cmap=cmap,
                         linewidths=.5)

    plt.gca().invert_yaxis()
    plt.xticks(rotation='vertical')
    plt.yticks(rotation='horizontal')

    if title is not None:
        plt.title(title)
    if x_label is not None:
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)
    if save:
        if title is None:
            filename = datetime.date.today.isoformat()
        else:
            filename = title.replace(" ", "")

        filename += ".pdf"
        plt.savefig(filename, format='pdf')
    else:
        return ax


def cartesianProduct(x, y):
    """
    Implementation of cartesian product between two np.arrays
    Parameters
    ----------
    x, y: np.array

    Returns
    -------
    np.array with cartesian product
    """
    return(np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))]))
