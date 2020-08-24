"""
Tool for solvating atomic structures written in .xyz format
Default solvent is Simple Point Charge (SPC) water
GROMACS-5.0 and ASE-3.9 or later required
Works with both orthogonal and non-orthogonal cells

author: Rasmus Kronberg
email: rasmus.kronberg@aalto.fi
"""

import sys, os
import numpy as np
from ase.io import read, write
from utils import angle, pbc, cleanup

def main():

	# Read structure
	try:
		xyz = read(sys.argv[1])
		print('Successfully read %s' % sys.argv[1])
	except:
		print('Usage: python3 solvate.py your_input.xyz a1 a2 a3 b1 b2 b3 c1 c2 c3')
		print('Lattice vectors a, b and c are optional.')
		exit()

	# Try setting user-defined cell vectors
	try:
		lattice = np.array([float(x) for x in sys.argv[2:]]).reshape(3,3)
		print('Setting user-defined cell vectors.')
		xyz.set_cell(lattice)
	except:
		print('Incomplete user-defined cell vectors, importing cell vectors from input file')

	# Write temporary .pdb file that GROMACS understands
	write('tmp.pdb', xyz)

	# Calculate box vector lengths and angles
	cell = xyz.get_cell()
	a, b, c = np.linalg.norm(cell[0]), np.linalg.norm(cell[1]), np.linalg.norm(cell[2])
	alpha, beta, gamma = angle(cell[1], cell[2]), angle(cell[0], cell[2]), angle(cell[0], cell[1])

	# Solvate
	os.system('2>/dev/null 1>&2 gmx editconf -f tmp.pdb -o tmp.gro -box %s %s %s \
		-angles %s %s %s' % (a/10., b/10., c/10., alpha, beta, gamma))
	os.system('2>/dev/null 1>&2 gmx solvate -cp tmp.gro -cs spc216.gro -o tmp_solv.pdb')

	solv = read('tmp_solv.pdb')
	elems = set(solv.get_chemical_symbols())

	# Wrap atoms in fractional coordinate space without breaking water molecules
	solv = pbc(xyz, solv, cell)

	# Convert back to .xyz
	write('tmp_solv.xyz', solv)

	# Cleanup
	outfile = cleanup('tmp_solv.xyz', elems)
	print('Solvated structure written in %s' % outfile)
	
if __name__ == '__main__':
	main()
