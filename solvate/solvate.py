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

def unitv(vector):

	# Returns the unit vector
	return vector/np.linalg.norm(vector)

def angle(u, v):
	
	# Returns the angle in degrees between vectors u and v
	uhat = unitv(u)
	vhat = unitv(v)
	return 360/(2*np.pi)*np.arccos(np.clip(np.dot(uhat, vhat), -1.0, 1.0))

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
					string = '%3s %10.5f %10.5f %10.5f \n' %(tmp[0], float(tmp[1]), \
					float(tmp[2]), float(tmp[3]))
					out.write(string)
				else:
					out.write(line)
				line = f.readline()

	os.system('rm -f tmp.pdb *tmp.gro* *tmp_solv.pdb* tmp_solv.xyz')
	print('Solvated structure written in %s' % outfile)

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

	# Convert back to .xyz
	write('tmp_solv.xyz', solv)

	# Cleanup
	cleanup('tmp_solv.xyz', elems)
	
if __name__ == '__main__':
	main()
