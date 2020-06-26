# compute sc command

LAMMPS compute to get the spatialc correlation for potential energy per atom

## Installation

Put *.cpp *.h into LAMMPS/src/, and compule LAMMPS as usual.

## Synatax

Syntax follow LAMMPS compute rdf command.
Difference is that for each pair at each bin, six components are output.  They are

- Number of pair
- $\sum_{ij} e_i\cdot e_j$
- $\sum_i e_i$
- $\sum_j e_j$
- $\sum_i e_i^2$
- $\sum_j e_j^2$

## Output 

One can use fix ave/time to get the average correlation over different timestep.
See example for usage
