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
   Contributing authors: Filippo Conforto (s2469797@ed.ac.uk)
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
  "fix dump smc command:\n\n";

/* ---------------------------------------------------------------------- */

/**
 * @brief Initializes the FixDUMPSMC fix for periodic dumping of SMC protein positions.
 * 
 * This constructor sets up the SMC position dump utility that periodically writes the locations
 * of SMC anchor and hinge positions to a text file for analysis and visualization. It parses
 * command-line arguments to configure dumping frequency, number of SMCs, output file path, and
 * references the parent FixSMC object containing the SMC data.
 * 
 * **Constructor Parameters** (after fix name and group):
 * - arg[3] (nevery): Integer specifying dump frequency (every nevery timesteps)
 * - arg[4] (nsmc): Number of SMC proteins to track
 * - arg[5] (dumpfile): Output filename prefix (without extension)
 * - arg[6] (fixname): Name of the FixSMC fix providing SMC position data
 * 
 * **Parsing Steps**:
 * 1. Register cite information if citation tracking enabled
 * 2. Validate argument count (must be exactly 7 arguments)
 * 3. Parse nevery: frequency of dump operations (must be > 0)
 * 4. Parse nsmc: number of SMCs to dump (must be > 0)
 * 5. Set restart_global=1 to enable SMC structure in restart file
 * 6. Store dumpFilestr: output filename (will become "filename.txt")
 * 7. Allocate and copy connFixName: stores FixSMC name for later retrieval
 * 
 * @param lmp [LAMMPS*] Pointer to LAMMPS instance
 * @param narg [int] Number of command-line arguments (should be 7)
 * @param arg [char**] Array of command-line arguments
 * @throws "Illegal fix dumpsmc command" - If narg != 7, nevery <= 0, or nsmc <= 0
 * @see FixDUMPSMC::init(), FixDUMPSMC::post_integrate()
 */
FixDUMPSMC::FixDUMPSMC(LAMMPS * lmp, int narg, char ** arg):
  Fix(lmp, narg, arg), connFix(nullptr){
    if (lmp -> citeme) lmp -> citeme -> add(cite_fix_dumpsmc);
 
    if ((narg != 7)) error -> all(FLERR, "Illegal fix dumpsmc command");

    nevery = utils::inumeric(FLERR, arg[3], false, lmp);
    if (nevery <= 0) error -> all(FLERR, "Illegal fix dumpsmc command");

    nsmc = utils::inumeric(FLERR, arg[4], false, lmp);
    if (nsmc <= 0) error -> all(FLERR, "Illegal fix dumpsmc command");

    restart_global = 1;
    dumpFilestr = arg[5];
    connFixName = new char[static_cast < int > (sizeof(arg[6]) / sizeof(char))];
    std::copy(arg[6], arg[6] + static_cast < int > (sizeof(arg[6]) / sizeof(char)), connFixName);
}

/**
 * @brief Destructor for FixDUMPSMC that cleans up allocated resources.
 * @see FixDUMPSMC::FixDUMPSMC()
 */
FixDUMPSMC::~FixDUMPSMC() {
}

/* ---------------------------------------------------------------------- */

int FixDUMPSMC::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

/**
 * @brief Initializes the dump fix by locating and caching the parent FixSMC pointer.
 * 
 * **Initialization Steps**:
 * 1. Call modify->get_fix_by_id(connFixName) to locate FixSMC by name
 * 2. Cache returned pointer in connFix member variable
 * 3. Validate that fix exists: if (!connFix), error "Illegal ausiliary Fix"
 * 
 * @see FixDUMPSMC::post_integrate(), FixDUMPSMC::connFix
 * @throws "Illegal ausiliary Fix" - If FixSMC with name connFixName not found
 */
void FixDUMPSMC::init() {

  // Store the connected Fix ID pointer
  connFix = modify -> get_fix_by_id(connFixName);
  if (!connFix) error -> all(FLERR, "Illegal ausiliary Fix");
}

/* -------------------------------------------*/
/**
 * @brief Periodically dumps SMC protein position data to a text file after integration.
 * 
 * This function is called after each timestep's integration and handles file I/O to record
 * SMC positions at specified intervals. It queries the parent FixSMC for current SMC anchor
 * and hinge positions and writes them to disk in a columnar text format. On the first timestep,
 * it creates/truncates the dump file; subsequently it appends data at intervals of nevery timesteps.
 * 
 * **Algorithm Steps**:
 * 1. On first timestep (update->ntimestep == 1):
 *    - Open dump file with truncate flag (erases previous content)
 *    - Close file immediately (creates empty file for later appending)
 * 2. Check dump frequency: if (update->ntimestep % nevery) return
 *    - Skip dumping if current timestep not a multiple of nevery
 * 3. If this is a dump timestep and processor rank is 0 (comm->me==0):
 *    - Open dump file in append mode (ios_base::app)
 *    - For each SMC i (0 to nsmc-1):
 *      * Query connFix->compute_array(i, 0): Get left end position (anchor)
 *      * Query connFix->compute_array(i, 1): Get right end position (hinge)
 *      * Write to file: "timestep  SMC_index  left_position  right_position\n"
 *    - Close file after all SMCs written
 * 4. Non-rank-0 processors skip file I/O (prevents concurrent writes)
 * 
 * **Output File Format**:
 * - Filename: "<dumpFilestr>.txt"
 * - Columns: timestep  SMC_number  left_end_bead_ID  right_end_bead_ID
 * - One line per SMC per dump timestep
 * - Space-separated text format (human-readable)
 * 
 * **Special Considerations**:
 * - Only rank-0 processor writes (avoids duplicate data in parallel runs)
 * - File created on timestep 1 to clear any previous simulation artifacts
 * - Append mode ensures multiple simulations can append to same file
 * - SMC indices are 1-indexed in output (i+1 instead of i)
 * 
 * @see FixDUMPSMC::init(), FixSMC::compute_array()
 */
void FixDUMPSMC::post_integrate() {

  // Opening the dump file on first iteration
  if (update -> ntimestep == 1){
    dfile.open(dumpFilestr + ".txt", std::fstream::trunc | std::fstream::out);
    dfile.close();
  }
  if (update -> ntimestep % nevery) return;
  
  // Writing dump information on selected file
  else{
    
    if (comm->me==0){      
      dfile.open(dumpFilestr + ".txt", std::ios_base::app);
      for (int i = 0; i < nsmc; i++)
        {
          // Writing on file timestep, number of LEF, left end position and right end position
          dfile << update -> ntimestep << " " << i+1 << " " <<connFix->compute_array(i,0) << " " << connFix->compute_array(i,1) << std::endl;
        }
      dfile.close();
    }
  } 
}

double FixDUMPSMC::memory_usage() {
  double bytes = 2 * nsmc * sizeof(long);
  return bytes;
}

/*-------------------------------------------*/
/*Add restart information in the restart file*/
/*-------------------------------------------*/
void FixDUMPSMC::write_restart(FILE * fp) {
  int rn = 0;
  long rlist[1];

  rlist[rn++] = static_cast < long > (update -> ntimestep);

  if (comm -> me == 0) {
    int size = rn * sizeof(long);
    fwrite( & size, sizeof(int), 1, fp);
    fwrite(rlist, sizeof(long), rn, fp);
  }

}

/*---------------------------------------------------*/
/*Use state info from restart file to restart the Fix*/
/*---------------------------------------------------*/
void FixDUMPSMC::restart(char * buf) {
  int rn = 0;
  long * rlist = (long * ) buf;

  bigint ntimestep_restart = static_cast < bigint > (rlist[rn++]);
  if (ntimestep_restart != update -> ntimestep)
    error -> all(FLERR, "Must not reset timestep when restarting fix smc");

  if (comm->me==0){ 

      // Copying the old dump file to a temporary one
      ifile.open(dumpFilestr + ".txt");
      dfilerst.open(dumpFilestr + "temp" + ".txt", std::fstream::trunc | std::fstream::out);
      for (int t = 0; t < (int) ntimestep_restart/nevery; t++){
        for (int i = 0; i < nsmc; i++)
          {
            ifile >> temptime >> tempnum >> tempside1 >> tempside2;
            dfilerst << temptime << " " << tempnum << " " << tempside1<< " " << tempside2 << std::endl;
          }
      }
      ifile.close();
      dfilerst.close();

      // Refilling the dump file until the restart time
      ifile.open(dumpFilestr + "temp" + ".txt");
      dfile.open(dumpFilestr + ".txt", std::fstream::trunc | std::fstream::out);
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