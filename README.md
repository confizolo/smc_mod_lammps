# SMC Fix for LAMMPS

# Introduction 
The following repository contains code related to the development of a Fix for LAMMPS for loop extruding proteins (see [here](10.7554/eLife.14864), [here](10.1016/j.bpj.2021.11.015) and [here](10.1016/j.celrep.2016.04.085) for reference). Development started in 11/22 and based on LAMMPS(23 Jun 2022 - Update 3)(source [here](https:-github.com/lammps/lammps/commit/88c8b6ec6feac6740d140393a0d409437f637f8b)).

The goal of this code is to introduce an easy and parallelised way to deploy and run loop extrusion on solution of polymers. To do so we developed a custom fix for LAMMPS that allows loop extrusion through moving bonds (see details [here](paper)).

# Repository structure 
- initfiles: contains sample datafiles to run tests
- exec: contains sample LAMMPS and bash scripts
- src: contain source files for fix_smc code

# How to run Fix SMC
The fix can be added to a regular LAMMPS script, by adding the following line with adequate parameters.

fix "fixname" "considered beads" smc "par1" "par2" "par3" "par4" "par5" "par6" "par7" "par8" "par9" "par10" "par11" "par12" "par13" "par14" "par15" "par16"

1. nevery: Attempt the jump every nevery iteration
2. seed: seed for random number draw
3. prob: probability to accept proposed movement
4. lpol: length of polymer(s) in solution
5. poltype: select shape of polymer(s) (either "linear" or "ring")
6. maxadir: maximum attempted movement of anchor (next attempted atom id: current anchor + adir)
7. maxhdir: minimum attempted movement of hinge (next attempted atom id: current hinge + hdir)
8. smcnum: number of deployed SMCs
9. smctype: atom type of beads representing SMCs' ends
10. smcbtype: SMCs' bond type after the first deployment
11. smcbitype: SMCs' bond type at the first deployment
12. cutoff: distance cutoff for attempted movements (these are is accepted only if distance between new anchor and hinge is below the cutoff)
13. initmode: define initialisation mode of extruders:
    1.  "random": deploys randomly the SMCs
    2.  "distributed": assign at least one SMC per polymer and then distribute remaining randomly
    3.  "full-distributed": distributes evenly SMCs over the polymers
14. kon: probability to load a free extruder every nevery step
15. koff: probability to unload an extruder every nevery step

# How to run Fix SMC Dump
An additional fix was developed to allow dumping of LEF positions (two ends). It can be exploited using the following command:

fix "fixname" "considered beads" dumpsmc "par1" "par2" "par3" "par4" 
1. nevery: Attempt the jump every nevery iteration
2. nsmc: number of deployed SMCs
3. dump_filename: name of file to use for dumping
4. fixname: name of the fix smc deployed in simulation

# Warnings

- At the current status SMCs can only be deployed all together at timestep 1 of simulation
- Simulations containing beads not belonging to polymers to extrude CAN'T be used given the depolyment algorithm for SMCs