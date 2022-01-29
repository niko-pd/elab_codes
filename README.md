# elab_codes
Examples of codes that I used to elaborate the raw MD data for the paper:

M. Liao, P. Nicolini, L. Du, J. Yuan, S. Wang, H. Yu, J. Tang, P. Cheng, K. Watanabe, T. Taniguchi,
L. Gu, V.E.P. Claerbout, A. Silva, D. Kramer, T. Polcar, R. Yang, D. Shi, G. Zhang,
Ultra-low friction and edge-pinning effect in large-lattice-mismatch van der Waals heterostructures,
Nature Materials, 21, 47 (2022), doi: 10.1038/s41563-021-01058-4

To run the scripts, raw MD data from LAMMPS (not uploaded on github because they are large files) are needed.

The codes perform the following tasks:

1-friction: it reads the external force applied on a group of atoms in order to enforce the sliding motion
	      and it calculates average and standard deviations, together with the accumulated work done on the system;

2-res: it trasforms the average values of force calculated by 1-friction and transform into shear and
	 edge-pinning strengths as reported in the paper;

3-all: it reads the LAMMPS trajectory file containing the atomic positions and for selected groups of atoms
	 (e.g. bulk or edge atoms) it calculates: average atomic displacements (w.r.t. a reference configuration),
	 velocities, forces, potential and kinetic energies, root mean square displacements (w.r.t. a reference configuration);

4-distribution: it reads the LAMMPS trajectory file containing the atomic positions and for selected groups of atoms
	 	  (e.g. bulk or edge atoms) it calculates the distributions of root mean square displacements (w.r.t.
		  a reference configuration) and potential energies;

5-maps: it reads the LAMMPS trajectory file containing the atomic positions and for selected groups of atoms
	 (e.g. bulk or edge atoms) it generates the datasets used to produce the maps reported in the paper.

