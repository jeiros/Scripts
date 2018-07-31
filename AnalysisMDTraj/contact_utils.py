"""
Functions for the contacts script
"""

import numpy as np
import mdtraj as md
import itertools
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
import datetime


def figure_dims(width_pt, factor=0.45):
    """
    I copied this from here:
    https://www.archer.ac.uk/training/course-material/2014/07/SciPython_Cranfield/Slides/L04_matplotlib.pdf
    """
    WIDTH = width_pt  # Figure width in pt (usually from LaTeX)
    FACTOR = factor  # Fraction of the width you'd like the figure to occupy
    widthpt = WIDTH * FACTOR
    inperpt = 1.0 / 72.27
    golden_ratio = (np.sqrt(5) - 1.0) / 2.0  # because it looks good
    widthin = widthpt * inperpt
    heightin = widthin * golden_ratio
    figdims = [widthin, heightin]  # Dimensions as list
    return figdims


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
{underline}
mask1: {mask1}
mask2: {mask2}
n_pairs : {n_pairs}\n""".format(
            name=self.dict['name'],
            underline='-' * len(self.dict['name']),
            mask1=self.dict['mask1'],
            mask2=self.dict['mask2'],
            n_pairs=len(self.dict['pairs'])
        )


def select_plotParams_fromTitle(title):
    if title == 'CcTnT - inhibitory peptide':
        y_step = True
        x_step = True
        y_lab = 'cTnI residue'
        x_lab = 'cTnT residue'
    if title == 'CcTnT - NcTnI':
        y_step = 2
        x_step = True
        y_lab = 'cTnI residue'
        x_lab = 'cTnT residue'
    if title == 'NcTnC - NcTnI':
        y_step = 2
        x_step = 4
        y_lab = 'cTnI residue'
        x_lab = 'cTnC residue'
    if title == 'NcTnC - switch peptide':
        y_step = True
        x_step = 4
        y_lab = 'cTnI residue'
        x_lab = 'cTnC residue'
    if title == 'NcTnI - inhibitory peptide':
        y_step = True
        x_step = 2
        y_lab = 'cTnI residue'
        x_lab = 'cTnI residue'
    if title == 'cTnC - inhibitory peptide':
        y_step = True
        x_step = 10
        y_lab = 'cTnI residue'
        x_lab = 'cTnC residue'
    if title == 'cTnC A-B - switch peptide':
        y_step = True
        x_step = 1
        y_lab = 'cTnI residue'
        x_lab = 'cTnC residue'

    return(y_step, x_step, y_lab, x_lab)


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


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
    for traj_chunk in traj_generator:
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


def renumber_mask(mask, top_fn=None):
    """
    Renumbers the mask to match the correct sequence of each
    subunit. Also takes into account the 0-indexed lists of
    Python.
    """
    if top_fn is not None:
        top = md.load_prmtop(top_fn)
    else:
        top = None
    if max(mask) > 160 and max(mask) < 249:
        # It's in cTnT
        if top is None:
            new_mask = [x + 51 for x in mask]
        else:
            new_mask = ['{} {}'.format(top.residue(x).name, x + 51) for x in mask]
    elif max(mask) >= 249:
        # It's in cTnI
        if top is None:
            new_mask = [x - 247 for x in mask]
        else:
            new_mask = ['{} {}'.format(top.residue(x).name, x - 247) for x in mask]
    else:
        # It's in cTnC
        if top is None:
            new_mask = [x + 1 for x in mask]
        else:
            new_mask = ['{} {}'.format(top.residue(x).name, x + 1) for x in mask]
    return(new_mask)


def plot_heatmap(contact_array, mask1, mask2, title=None, save=False,
                 x_label=None, y_label=None,
                 x_steps=True, y_steps=True,
                 min_value=0, max_value=1,
                 std_array=None, top_fn=None,
                 cbar_label='Contact frequency',
                 cmap='Reds'):
    """
    Plot a single heatmap on a red color scale.
    """
    fig, ax = plt.subplots(figsize=figure_dims(2000))
    # Check what part of cTn the masks are in and renumber accordingly
    new_mask1 = renumber_mask(mask1, top_fn=top_fn)
    new_mask2 = renumber_mask(mask2, top_fn=top_fn)

    # Convert to a pd.DataFrame so sns.heatmap() can read it's
    # row/column labels properly
    contact_df = pd.DataFrame(contact_array.T, index=new_mask2,
                              columns=new_mask1)

    if std_array is not None:
        ax = sns.heatmap(
            contact_df,
            vmin=min_value,
            vmax=max_value,
            xticklabels=x_steps,
            yticklabels=y_steps,
            cmap=cmap,
            linewidths=.5,
            annot=std_array,
            fmt='.2f',
            annot_kws={'rotation': 'vertical'},
            cbar_kws={'label': cbar_label}
        )
    else:
        ax = sns.heatmap(
            contact_df,
            vmin=min_value,
            vmax=max_value,
            xticklabels=x_steps,
            yticklabels=y_steps,
            cmap=cmap,
            linewidths=.5,
            cbar_kws={'label': cbar_label}
        )

    plt.xticks(rotation=45)
    plt.yticks(rotation='horizontal')

    if title is not None:
        plt.title(title)
    if x_label is not None:
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)

    fig.tight_layout()
    if save:
        if title is None:
            filename = datetime.date.today().isoformat()
        else:
            filename = title.replace(" ", "")

        filename += ".pdf"
        plt.savefig(filename, format='pdf')
    else:
        return ax


def plot_diffmap(contact_array, mask1, mask2, title=None, save=False,
                 x_label=None, y_label=None,
                 x_steps=True, y_steps=True,
                 top_fn=None, pval_arr=None,
                 cbar_label='$\Delta$ Contact frequency'):
    """
    Plots a difference heatmap with a blue-grey-red diverging
    color scale.
    """
    fig, ax = plt.subplots(figsize=figure_dims(2000))

    # Check what part of cTn the masks are in and renumber accordingly
    new_mask1 = renumber_mask(mask1, top_fn=top_fn)
    new_mask2 = renumber_mask(mask2, top_fn=top_fn)

    contact_df = pd.DataFrame(contact_array.T, index=new_mask2,
                              columns=new_mask1)

    # Diverging palette with light colour on the middle
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    maxval = max(abs(contact_array.min()), contact_array.max())

    if pval_arr is not None:
        ax = sns.heatmap(
            contact_df,
            vmin=-maxval,
            vmax=maxval,
            xticklabels=x_steps,
            yticklabels=y_steps,
            annot=pval_arr.T,
            annot_kws={'rotation': 'horizontal'},
            fmt='.2f',
            mask=pval_arr.T >= 0.01,

            cmap=cmap,
            linewidths=.5,
            cbar_kws={'label': cbar_label}
        )
    else:
        ax = sns.heatmap(
            contact_df,
            vmin=-maxval,
            vmax=maxval,
            xticklabels=x_steps,
            yticklabels=y_steps,
            cmap=cmap,
            linewidths=.5,
            cbar_kws={'label': cbar_label}
        )

    plt.xticks(rotation=45)
    plt.yticks(rotation='horizontal')

    if title is not None:
        plt.title(title)
    if x_label is not None:
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)

    fig.tight_layout()
    if save:
        if title is None:
            filename = datetime.date.today().isoformat()
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
