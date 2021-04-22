"""
Tool for solvating atomic structures written in .xyz format
Default solvent is Simple Point Charge (SPC) water
GROMACS-5.0 and ASE-3.9 or later required
Works with both orthogonal and non-orthogonal cells

author: Rasmus Kronberg
email: rasmus.kronberg@aalto.fi
"""

import os
import numpy as np
from ase.io import read, write
from argparse import ArgumentParser

from utils import angle, pbc, cleanup


def main():

    xyz = read(inp)
    print('Successfully read %s' % inp)

    # Try setting user-defined cell vectors
    if not lattice:
        print('No user-defined cell vectors, importing vectors from input')
    else:
        print('Setting user-defined cell vectors.')
        xyz.set_cell(lattice)

    # Write temporary .pdb file that GROMACS understands
    write('tmp.pdb', xyz)

    # Calculate box vector lengths and angles
    cell = xyz.get_cell()
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    alpha = angle(cell[1], cell[2])
    beta = angle(cell[0], cell[2])
    gamma = angle(cell[0], cell[1])

    # Solvate
    os.system('2>/dev/null 1>&2 gmx editconf -f tmp.pdb -o tmp.gro \
              -box %s %s %s -angles %s %s %s'
              % (a/10., b/10., c/10., alpha, beta, gamma))
    os.system('2>/dev/null 1>&2 gmx solvate -cp tmp.gro -cs spc216.gro \
              -o tmp_solv.pdb')

    solv = read('tmp_solv.pdb')
    elems = set(solv.get_chemical_symbols())

    # Wrap atoms in fractional coordinate space without
    # breaking water molecules
    solv = pbc(xyz, solv, cell)

    # Convert back to .xyz
    write('tmp_solv.xyz', solv)

    # Cleanup
    outfile = cleanup(inp, 'tmp_solv.xyz', elems)
    print('Solvated structure written to %s' % outfile)


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Script for solvating .xyz atomic coordinate file')
    parser.add_argument('-i', '--input', required=True, help='Input file')
    parser.add_argument('--cell', default=[], nargs='+', help='Cell vectors')
    args = parser.parse_args()
    inp = args.input
    lattice = args.cell
    main()
