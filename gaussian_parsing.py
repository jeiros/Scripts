from cclib.io import ccread
import numpy as np
import glob
from cclib.parser import utils

symbols = [
    'Ac',
    'Ag',
    'Al',
    'Am',
    'Ar',
    'As',
    'At',
    'Au',
    'B',
    'Ba',
    'Be',
    'Bh',
    'Bi',
    'Bk',
    'Br',
    'C',
    'Ca',
    'Cd',
    'Ce',
    'Cf',
    'Cl',
    'Cm',
    'Cn',
    'Co',
    'Cr',
    'Cs',
    'Cu',
    'Db',
    'Ds',
    'Dy',
    'Er',
    'Es',
    'Eu',
    'F',
    'Fe',
    'Fl',
    'Fm',
    'Fr',
    'Ga',
    'Gd',
    'Ge',
    'H',
    'He',
    'Hf',
    'Hg',
    'Ho',
    'Hs',
    'I',
    'In',
    'Ir',
    'K',
    'Kr',
    'La',
    'Li',
    'Lr',
    'Lu',
    'Lv',
    'Md',
    'Mg',
    'Mn',
    'Mo',
    'Mt',
    'N',
    'Na',
    'Nb',
    'Nd',
    'Ne',
    'Ni',
    'No',
    'Np',
    'O',
    'Os',
    'P',
    'Pa',
    'Pb',
    'Pd',
    'Pm',
    'Po',
    'Pr',
    'Pt',
    'Pu',
    'Ra',
    'Rb',
    'Re',
    'Rf',
    'Rg',
    'Rh',
    'Rn',
    'Ru',
    'S',
    'Sb',
    'Sc',
    'Se',
    'Sg',
    'Si',
    'Sm',
    'Sn',
    'Sr',
    'Ta',
    'Tb',
    'Tc',
    'Te',
    'Th',
    'Ti',
    'Tl',
    'Tm',
    'U',
    'Uuo',
    'Uup',
    'Uus',
    'Uut',
    'V',
    'W',
    'Xe',
    'Y',
    'Yb',
    'Zn',
    'Zr'
]


def load_log_files(fnames):
    fname_list = glob.glob(fnames)
    return [ccread(log) for log in fname_list]


def get_scf_energies(data_list, unit1='eV', unit2='kcal'):
    return [utils.convertor(d.scfenergies[0], unit1, unit2) for d in data_list]


def boltzmann_distribution(energy_dict, R=1.9872036e-3, T=300):
    return [np.exp(-(x[1] / (R * T)))for x in energy_dict.values()]


def plot_conformer_population(energies_list, R=1.9872036e-3, T=300, title=None):
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
