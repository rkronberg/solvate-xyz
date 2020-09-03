'''
Utility functions for solvate.py

author: Rasmus Kronberg
email: rasmus.kronberg@aalto.fi
'''

import sys
import os
import numpy as np
from ase.io import write


def unitv(vector):

    # Returns the unit vector
    return vector/np.linalg.norm(vector)


def angle(u, v):

    # Returns the angle in degrees between vectors u and v
    uhat = unitv(u)
    vhat = unitv(v)
    return 360/(2*np.pi)*np.arccos(np.clip(np.dot(uhat, vhat), -1.0, 1.0))


def pbc(xyz, solv, cell):

    for atom in solv:
        # Original atoms and water oxygens
        if atom.index < len(xyz) or atom.symbol == 'O':
            frac = np.matmul(np.linalg.inv(cell.T), atom.position)
            frac = (frac + 1.0) % 1.0
            cart = np.matmul(cell.T, frac)
            shift = cart - atom.position
            atom.position = cart
        # Shift water hydrogens accordingly
        if atom.index >= len(xyz) and atom.symbol == 'O':
            solv[atom.index + 1].position += shift
            solv[atom.index + 2].position += shift

    return solv


def cleanup(file, elems):

    # Cleanup .xyz so that programs such as ase-gui understands the content
    # Fix header
    os.system("sed -i 's/Properties.*//' %s" % file)

    # Strip everything but symbols and coordinates
    outfile = '%s-solvated.xyz' % sys.argv[1].strip('.xyz')
    with open(outfile, 'w') as out:
        with open(file, 'r') as f:
            line = f.readline()
            while line != '':
                tmp = line.strip().split()
                if any(sym == tmp[0] for sym in elems):
                    string = '%3s %10.5f %10.5f %10.5f \n'
                    % (tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]))
                    out.write(string)
                else:
                    out.write(line)
                line = f.readline()

    os.system('rm -f tmp.pdb *tmp.gro* *tmp_solv.pdb* tmp_solv.xyz')

    return outfile
