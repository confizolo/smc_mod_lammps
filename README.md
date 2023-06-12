# SMC Fix for LAMMPS

# Introduction 
The following repository contains code related to the development of a Fix for LAMMPS for loop extruding proteins (see [here](10.7554/eLife.14864), [here](10.1016/j.bpj.2021.11.015) and [here](10.1016/j.celrep.2016.04.085) for reference). Development started in 11/22 and based on LAMMPS(23 Jun 2022 - Update 2)(source [here](https:-github.com/lammps/lammps/commit/88c8b6ec6feac6740d140393a0d409437f637f8b)). 
# What can Fix SMC do?
- Deploy multiple LEFs
- Loop extrude along multiple directions
- Take into accounts border effects due to LEFs interaction or polymer length in monodisperse melts
- Account for binding and unbinding of LEFs
- Stochastically move each feet to a maximum distance
- Block movements if encounters roadblocks

# Repository structure 
- data: contains ready-to-go data file for LAMMPS 
- exec: contains multiple bash and LAMMPS scripts for tests and simulations
- src: contain source files for fix_smc code
- test-analysis: contains analysis tools for tests

# How to run Fix SMC

The fix can be run inside a lammps script, by adding the following line:

fix "fixname" "considered beads" smc "par1" "par2" "par3" "par4" "par5" "par6" "par7" "par8" "par9" "par10" "par11" "par12" "par13" "par14" "par15" "par16"

1. nevery: Attempt the jump every nevery iteration
2. seed: random seed
3. prob: probability to attempt the jump
4. lpol: length of polymer(s)
5. poltype: give the form of the polymer (either linear or ring)
6. maxadir: maximum attempted movement of anchor (next attempted atom id: current anchor + adir)
7. maxhdir: minimum attempted movement of hinge (next attempted atom id: current hinge + hdir)
8. smcnum: number of deployed smcs
9. smctype: atom type of anchoring beads
10. smcbtype: bond type of anchoring beads after the first deployment
11. smcbitype: bond type of anchoring beads at the first deployment
12. cutoff: distance cutoff for attempted movements (Jump is accepted only if distance between new anchor and hinge is below the cutoff)
13. initmode: random or distributed according to uswe
14. kon: loading probability 
15. koff: unloading probability
16. blockbeads: type of beads that the extruder cannot grab, can be listed as an arbitrary long list (e.g.: 2 3 4 ...)


# Things to add
- Substitute type check with default LAMMPS mask