#!/usr/bin/env python

from atomic_symbols import symbols
import argparse
parser = argparse.ArgumentParser(prog='create_Gaussian_scripts.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''
____________________________________________________________________________
| A program that generates multiple Gaussian scripts based on a template   |
| and an xyz file with several molecules                                   |
----------------------------------------------------------------------------
''')

parser.add_argument('-xyz', '--xyzfile', type=str, required=True)
parser.add_argument('-t', '--template', type=str, required=True)
parser.add_argument('-o', '--outputfile', type=str, required=True)


def get_xyz_confs(fname):
    """
    Parse the conformers in an .xyz file with several entries
    Parameters
    ----------
    fname: str, name of the .xyz file

    Returns
    -------
    molecules: list
        A list of length equal to the number of conformers in the xyz file
        Each entry of the list is a list of lenght equal to the number of atoms
        (It is assumed all entries in the xyz file are of the same molecule)
        The atom list is stored as a list of two-dimensional tuples:
            0 -> str atom symbol
            1 -> list of floats [x, y, z] coordinates of the atom
    """
    molecules = []
    with open(str(fname)) as f:
        # First line tells us how many atoms the molecule has
        n_atoms = int(f.readline())
        conformer = []
        for line in f:
            # Ignore empty lines
            if len(line.split()) > 0:
                if line.split()[0] in symbols:
                    split_line = line.split()
                    atom = split_line[0]
                    coords = [float(x) for x in split_line[1:]]
                    conformer.append((atom, coords))
                if len(conformer) == n_atoms:
                    molecules.append(conformer)
                    conformer = []

    return molecules


def read_template_file(fname):
    with open(fname, 'r') as f:
        pre_cmds = f.read()
    return pre_cmds


def write_Gaussian_script(out_fname, molecule):
    n_conformers = len(molecule)
    pre_cmds = read_template_file(args.template)
    for i, conformer in enumerate(molecule):
        with open(out_fname + '%02d.com' % (i + 1), 'w') as f:
            f.write(pre_cmds)
            for atom in conformer:
                f.write(' %s\t%.6f\t%.6f\t%.6f\n' % (atom[0],
                                                     atom[1][0],
                                                     atom[1][1],
                                                     atom[1][2]))
    print('Succesfully wrote %d files' % n_conformers)

if __name__ == '__main__':
    args = parser.parse_args()
    write_Gaussian_script(out_fname=args.outputfile,
                          molecule=get_xyz_confs(args.xyzfile))
