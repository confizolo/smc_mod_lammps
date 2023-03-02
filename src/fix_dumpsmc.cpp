// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Filippo Conforto
------------------------------------------------------------------------- */

#include "fix_dumpsmc.h"

// The atom class provides access to various atom properties, including the atom type, molecular ID, position, velocity, and force
#include "atom.h"

// Came by default in fix_swap_atom
#include <cmath>
#include <cctype>
#include <cfloat>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "compute.h"
#include "modify.h"
#include "group.h"
#include "domain.h"
#include "region.h"
#include "random_park.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include<array> 

#include <fix_bond_history.h>

// Davide include them (with more)
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "citeme.h"
#include <stdlib.h>
#include <time.h>  /* Important for random number generator */
#include <sstream> // std::stringstream
#include <fstream> //ofstream
#include <algorithm>
#include <iterator>

#include "input.h"
#include "variable.h"
#include <iostream>
#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;

static
const char cite_fix_dumpsmc[] =
  "fix smc command:\n\n";

/* ---------------------------------------------------------------------- */

FixDUMPSMC::FixDUMPSMC(LAMMPS * lmp, int narg, char ** arg):
  Fix(lmp, narg, arg), connFix(nullptr){
    if (lmp -> citeme) lmp -> citeme -> add(cite_fix_dumpsmc);
 
    if ((narg != 7)) error -> all(FLERR, "Illegal fix dumpsmc command");

    nevery = utils::inumeric(FLERR, arg[3], false, lmp);
    if (nevery <= 0) error -> all(FLERR, "Illegal fix dumpsmc command");

    nsmc = utils::inumeric(FLERR, arg[4], false, lmp);
    if (nsmc <= 0) error -> all(FLERR, "Illegal fix dumpsmc command");

    // Flag to activate dump in restart file of fix smc structure
    restart_global = 1;

    dumpFile = utils::get_potential_file_path(arg[5]);

    connFixName = new char[static_cast < int > (sizeof(arg[6]) / sizeof(char))];
    std::copy(arg[6], arg[6] + static_cast < int > (sizeof(arg[6]) / sizeof(char)), connFixName);
    

  }

/* ---------------------------------------------------------------------- */

FixDUMPSMC::~FixDUMPSMC() {
}

/* ---------------------------------------------------------------------- */

int FixDUMPSMC::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDUMPSMC::init() {

  // Store the connected Fix ID pointer
  connFix = modify -> get_fix_by_id(connFixName);
  if (!connFix) error -> all(FLERR, "Illegal ausiliary Fix");
}

/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

void FixDUMPSMC::post_integrate() {
  if (update -> ntimestep == 1){
    dfile.open(dumpFile + ".txt", std::fstream::trunc | std::fstream::out);
    dfile.close();
  }
  if (update -> ntimestep % nevery) return;
  
  else{

    if (comm->me==0){      
      dfile.open(dumpFile + ".txt", std::ios_base::app);
      for (int i = 0; i < nsmc; i++)
        {
          dfile << update -> ntimestep << " " << i+1 << " " <<connFix->compute_array(i,0) << " " << connFix->compute_array(i,1) << std::endl;
        }
      dfile.close();
    }
  } 
}

/* ----------------------------------------------------------------------
   memory usage 
------------------------------------------------------------------------- */

double FixDUMPSMC::memory_usage() {
  double bytes = 2 * nsmc * sizeof(long);
  return bytes;
}

/*---------------------------------------------------------------------*/
/* Needed to write a restart file that can continue with the simulation*/
/*---------------------------------------------------------------------*/
void FixDUMPSMC::write_restart(FILE * fp) {
  // if ((debug)) utils::logmesg(lmp, "Writing restart for fix_smc \n");

  int rn = 0;
  long rlist[1];

  rlist[rn++] = static_cast < long > (update -> ntimestep);

  if (comm -> me == 0) {
    int size = rn * sizeof(long);
    fwrite( & size, sizeof(int), 1, fp);
    fwrite(rlist, sizeof(long), rn, fp);
  }

}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */
void FixDUMPSMC::restart(char * buf) {
  int rn = 0;
  long * rlist = (long * ) buf;

  bigint ntimestep_restart = static_cast < bigint > (rlist[rn++]);
  if (ntimestep_restart != update -> ntimestep)
    error -> all(FLERR, "Must not reset timestep when restarting fix smc");

  if (comm->me==0){ 

      ifile.open(dumpFile + ".txt");
      dfilerst.open(dumpFile + "temp" + ".txt", std::fstream::trunc | std::fstream::out);
      for (int t = 0; t < (int) ntimestep_restart/nevery; t++){
        for (int i = 0; i < nsmc; i++)
          {
            ifile >> temptime >> tempnum >> tempside1 >> tempside2;
            dfilerst << temptime << " " << tempnum << " " << tempside1<< " " << tempside2 << std::endl;
          }
      }
      ifile.close();
      dfilerst.close();

      ifile.open(dumpFile + "temp" + ".txt");
      dfile.open(dumpFile + ".txt", std::fstream::trunc | std::fstream::out);
      for (int t = 0; t < (int) ntimestep_restart/nevery; t++){
        for (int i = 0; i < nsmc; i++)
          {
            ifile >> temptime >> tempnum >> tempside1 >> tempside2;
            dfile << temptime << " " << tempnum << " " << tempside1<< " " << tempside2 << std::endl;
          }
      }

      ifile.close();
      dfile.close();

    }


}