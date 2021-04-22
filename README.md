# solvate-xyz
Script for solvating atomic structures written in .xyz format. Works with both orthogonal and non-orthogonal cells. Solvent structure is generated with GROMACS (v. 5.0 or later) using the Simple Point Charge (SPC) water model. ASE v. 3.9 or later is required.

## Usage

```bash
$ python3 solvate.py your_input.xyz a1 a2 a3 b1 b2 b3 c1 c2 c3
```
The Cartesian cell vectors ```[a1, a2, a3]```, ```[b1, b2, b3]``` and ```[c1, c2, c3]``` are optional. If not specified, the code tries to read the cell vectors from the input file.
