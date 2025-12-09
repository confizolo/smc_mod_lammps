# SMC Tracking Analysis Utilities

This folder contains lightweight Python helpers used to post-process LAMMPS
trajectories produced by the SMC/bridging simulations. Two modules are
provided:

1. ``read_files.py``
	Minimal readers for LAMMPS ``.lammpstrj`` trajectory dumps and bond dump
	files. Each frame is returned as a dictionary with unwrapped/wrapped
	coordinates, particle indices, types, molecule identifiers, and the box
	lengths. The bond reader returns ``(atom1, atom2, type)`` arrays sorted by
	atom index to simplify comparisons between frames.

2. ``contact_dens.py``
	Tool for measuring intra- and inter-polymer SMC contact counts. Users can provide the trajectory path, the output directory, bead
	types to include, and frame slicing options.

3. ``environment.yml``
   Contains all the necessary libraries for running the scripts.

