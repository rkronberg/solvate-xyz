# solvate-xyz
Script for solvating atomic structures written in .xyz format. Works with both orthogonal and non-orthogonal cells. Solvent structure is generated with GROMACS (v. 5.0 or later) using the Simple Point Charge (SPC) water model. ASE v. 3.9 or later is required.

## Usage

```console
$ python solvate.py --help
usage: solvate.py [-h] -i INPUT
                  [--cell CELL [CELL ...]]

Script for solvating .xyz atomic coordinate file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file
  --cell CELL [CELL ...]
                        Cell vectors
```

The Cartesian cell vectors are optional. If not specified, the code tries to read the cell vectors from the input file.
