from cclib.io import ccread
import numpy as np
import glob
from cclib.parser import utils
import mdtraj
from matplotlib import pyplot as plt
import pandas as pd

symbols = utils.PeriodicTable().element[1:]  # drop first item which is None


def load_log_files(fnames):
    fname_list = glob.glob(fnames)
    print('Filenames are : ')
    [print(fname) for fname in fname_list]
    return [ccread(log) for log in fname_list]


def get_scf_energies(data_list, unit1='eV', unit2='kcal'):
    return [utils.convertor(d.scfenergies[0], unit1, unit2) for d in data_list]


def boltzmann_distribution(energy_dict, R=1.9872036e-3, T=300):
    return [np.exp(-(x[1] / (R * T)))for x in energy_dict.values()]


def plot_conformer_population(energies_list, R=1.9872036e-3, T=300, title=None):
    '''
    Parameters
    ----------
    energies_list: list, List of energies (could be obtained with get_scf_energies)
    R: float, the gas or boltzmann constant
    T: float, the temperature
    title: str (default=None), the title of the plot

    Returns
    -------
    ax: matplotlib.axis object
    ene_dict: dictionary, keys are ints with the conformer ID and the values
    are bidimensional tuple of the energy of said conformer and it's energy
    difference w.r.t the lowest energy conformer in the list
    '''
    # For each conformer, build a dictionary with it's energy value
    # and it's difference with respect to the lowest energy conformer
    ene_dict = dict.fromkeys(range(1, len(energies_list) + 1))
    for i, energy in enumerate(energies_list):
        ene_dict[i + 1] = (energy, energy - min(energies_list))
    # Boltzmann distribution
    exp_terms = boltzmann_distribution(ene_dict)
    # Do the plot
    f, ax = plt.subplots()
    ax.bar(range(1, len(energies_list) + 1), exp_terms / sum(exp_terms) * 100)
    ax.set_ylabel('Population %')
    ax.set_xlabel('Conformer ID')
    if title is not None:
        ax.set_title(title)
    return ax, ene_dict


def generate_traj(xyz_array, pdb_file):
    top = mdtraj.load_pdb(pdb_file).topology
    return mdtraj.Trajectory(xyz=xyz_array / 10, topology=top)


def compare_spe_energies(spe_gas_dict, spe_pcm_dict):
    df_gas = pd.DataFrame({
        'Conformer': range(1, len(spe_gas_dict) + 1),
        'Population': np.array(boltzmann_distribution(spe_gas_dict)) * 100 / np.array(boltzmann_distribution(spe_gas_dict)).sum(),
        'Method': 'Gas Phase'
    })
    df_pcm = pd.DataFrame({
        'Conformer': range(1, len(spe_pcm_dict) + 1),
        'Population': np.array(boltzmann_distribution(spe_pcm_dict)) * 100 / np.array(boltzmann_distribution(spe_pcm_dict)).sum(),
        'Method': 'PCM'
    })
    return pd.concat([df_gas, df_pcm])
