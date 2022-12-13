# SMC Fix for LAMMPS

# Introduction 
The following repository contains code related to the development of a Fix for LAMMPS for loop extruding proteins (see [here](10.7554/eLife.14864), [here](10.1016/j.bpj.2021.11.015) and [here](10.1016/j.celrep.2016.04.085) for reference). Development started in 11/22 and based on LAMMPS(23 Jun 2022 - Update 2)(source [here](https://github.com/lammps/lammps/commit/88c8b6ec6feac6740d140393a0d409437f637f8b)). 
# What can Fix SMC do?
- Deploy multiple loop extruders
- Loop extrude along multiple directions
- Take into accounts border effects due to extruders interaction or polymer length in monodisperse melts
  
# Repository structure 
- data: contains ready-to-go data file for LAMMPS 
- exec: contains multiple bash and LAMMPS scripts for tests and simulations
- src: contain source files for fix_smc code
- test-analysis: contains analysis tools for tests

# How to perform analyses
The code in test-analysis requires scripts coming from a personal analysis library downloadable from [here](https://git.ecdf.ed.ac.uk/s2469797/analysis_scripts).

