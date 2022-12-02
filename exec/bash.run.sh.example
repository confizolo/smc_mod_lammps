#!/bin/bash
#
#SBATCH --job-name=maxdens.topo2.test
#SBATCH --output=maxdens.topo2.test.o
#
#SBATCH --ntasks=4
#SBATCH --partition=long
#SBATCH --time=168:00:00

#SBATCH --mem-per-cpu=512
LAMMPS=""
LAMMPSparallel="/Users/s2469797/Documents/lammps_testing/build/lmp"

lscript="test.simulation.lammps"

ncpu=4

cwd=`pwd`
 
#   r1=$(od -An -N2 -i /dev/urandom)
#   r2=$(od -An -N2 -i /dev/urandom)
#   echo "variable seedthermo0 equal ${r1}" >  parameters.dat
#   echo "variable seedthermo  equal ${r2}" >> parameters.dat
  
${LAMMPSparallel} -in ${lscript}
