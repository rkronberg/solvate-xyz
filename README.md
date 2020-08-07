# solvate-xyz
Script for solvating atomic structures written in .xyz format. Works with both orthogonal and non-orthogonal cells. Solvent structure is generated with GROMACS (v. 5.0 or later) using the Simple Point Charge (SPC) water model. ASE version 3.9 or later is required.

## Usage

```bash
$ python3 solvate.py your_input.xyz a1 a2 a3 b1 b2 b3 c1 c2 c3
```
The Cartesian cell vectors <img src="https://latex.codecogs.com/svg.latex?\mathbf{a}=[a_1,a_2,a_3]"/>, <img src="https://latex.codecogs.com/svg.latex?\mathbf{b}=[b_1,b_2,b_3]"/> and <img src="https://latex.codecogs.com/svg.latex?\mathbf{c}=[c_1,c_2,c_3]"/> are optional. If not specified, the code tries to read the cell vectors from the input file.
