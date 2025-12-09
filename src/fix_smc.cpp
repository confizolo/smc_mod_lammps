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
   Contributing authors: Filippo Conforto (s2469797@ed.ac.uk), Choong Zheng Yang (zchoong001@e.ntu.edu.sg)
------------------------------------------------------------------------- */

#include "fix_smc.h"

// The atom class provides access to various atom properties, including the atom type, molecular ID, position, velocity, and force
#include "atom.h"

// Came by default in fix_swap_atom
#include <cmath>
#include <cctype>
#include <cfloat>
#include <cstring>
#include <iostream>
#include <fstream>
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
const char cite_fix_smc[] =
  "fix smc command:\n\n";

/* ---------------------------------------------------------------------- */

/**
 * @brief Initializes FixSMC for simulating Structural Maintenance of Chromosomes (SMC) proteins.
 * 
 * This constructor sets up a LAMMPS fix for simulating loop extrusion dynamics of SMC proteins
 * on DNA polymer chains. SMCs bind DNA at two positions (anchor and hinge) and can extrude loops
 * by moving one or both ends along the polymer. The fix supports multiple initialization modes,
 * bidirectional/monodirectional movement, probabilistic loading/unloading, and advanced features
 * like angle constraints and blocked bead types.
 * 
 * **Command Arguments** (arg[3:22+]):
 * - arg[3] (nevery): Update frequency (timesteps between extrusion attempts)
 * - arg[4] (seed): Random number seed
 * - arg[5] (prob): Acceptance probability for movements [0-1]
 * - arg[6] (lpol): Polymer length (beads per polymer)
 * - arg[7] (solsize): Total number of polymer beads
 * - arg[8] (poltype): "linear" or "ring" polymer topology
 * - arg[9] (dirmode): "bo" (bidirectional), "mo1" (anchor fixed), "mo2" (hinge fixed), "mora" (alternating)
 * - arg[10] (maxdir): Maximum movement distance per step
 * - arg[11] (smcnum): Number of SMC proteins
 * - arg[12] (smctype): Atom type for SMC ends
 * - arg[13] (smcbtype): Bond type for SMCs (after initial placement)
 * - arg[14] (smcbitype): Bond type for initial SMC placement
 * - arg[15] (cutoff): Distance cutoff for movement acceptance
 * - arg[16] (tancoff): Tangent product cutoff for angle constraint
 * - arg[17] (initmode): "random", "distributed", "full-distributed", or "fixed"
 * - arg[18] (kon): Loading probability [0-1]
 * - arg[19] (koff): Unloading probability [0-1]
 * - arg[20] (npatches): Number of patch atoms per bead
 * - arg[21] (debug): Enable debug CSV logging (0=off, 1=on)
 * - arg[22] (fixFname): Filename for initial positions (if initmode=fixed)
 * - arg[23+] (blockt): List of bead types SMCs cannot bind
 * 
 * **Initialization Steps**:
 * 1. Validate argument count (minimum 19 required)
 * 2. Parse nevery frequency and integration settings
 * 3. Parse random seed and movement probability
 * 4. Parse polymer parameters (length, size, topology)
 * 5. Parse movement mode and maximum distance
 * 6. Parse SMC parameters (number, types, bond types)
 * 7. Parse acceptance criteria (distance and angle cutoffs)
 * 8. Parse initialization mode and loading/unloading rates
 * 9. Allocate anchor/hinge position arrays
 * 10. Initialize based on initmode (random, distributed, or from file)
 * 11. Allocate availability list and blocked bead type array
 * 12. Create random number generator with processor-unique seed
 * 
 * @param lmp [LAMMPS*] Pointer to LAMMPS instance
 * @param narg [int] Number of command arguments (minimum 19)
 * @param arg [char**] Command argument array
 * @throws "Illegal fix smc command" - Invalid arguments or parameter ranges
 * @throws "Illegal position of hinge or anchor read from file" - Invalid initial positions
 * @see FixSMC::post_integrate(), FixSMC::load_smc(), FixSMC::place_smc()
 */
FixSMC::FixSMC(LAMMPS * lmp, int narg, char ** arg):
  Fix(lmp, narg, arg),
  anch(nullptr), hing(nullptr), smctype(0), smcbtype(0), smcnum(0), debug(0) {
    if (lmp -> citeme) lmp -> citeme -> add(cite_fix_smc);
    
    // Check on the number of arguments given to the fix
    if (narg < 19) error -> all(FLERR, "Illegal fix smc command");

    // Define attempt rate
    nevery = utils::inumeric(FLERR, arg[3], false, lmp);
    if (nevery <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Flag to activate dump in restart file
    restart_global = 1;

    // Activate flag for array and scalar return
    extscalar=1;
    extarray=1;
    scalar_flag = 1;
    array_flag = 1;

    // Define the seed using during generation of random numbers
    seed = utils::inumeric(FLERR, arg[4], false, lmp);

    // Define acceptance probability for movements
    prob = utils::numeric(FLERR, arg[5], false, lmp);
    if (prob < 0.0 || prob > 1.0)
      error -> all(FLERR, "Illegal fix topo2 command");

    // Define length of polymer(s) present in simulation to check boundary conditions
    lpol = utils::inumeric(FLERR, arg[6], false, lmp);
    if (lpol <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Define size of the solution in number of beads
    solsize = utils::inumeric(FLERR, arg[7], false, lmp);
    if (solsize < 0)
      error -> all(FLERR, "Illegal fix smc command, number of beads smaller than 0");

    // Define type of extrusion according to the kind of polymer(s) present in simulation
    if (strcmp(arg[8], "linear") == 0) {
      ring = 0;
    } else if (strcmp(arg[8], "ring") == 0) {
      ring = 1;
    } else {
      error -> all(FLERR, "Illegal fix smc command, indefinite polymer type");
    }
    
    
    // Define type of extrusion direction
    if (strcmp(arg[9], "bo") == 0) {
      dirmode = 0;
    } else if (strcmp(arg[9], "mo1") == 0) {
      dirmode = 1;
    } else if (strcmp(arg[9], "mo2") == 0) {
      dirmode = 2;
    } else if (strcmp(arg[9], "mora") == 0) {
      dirmode = 3;
    } else {
      error -> all(FLERR, "Illegal fix smc command, no definite direction");
    }

    // Define maximum movement of SMCs' hinge
    maxdir = utils::inumeric(FLERR, arg[10], false, lmp);
    
    if (maxdir < 0) error -> all(FLERR, "Illegal fix smc command, hinge and must not have same direction");

    // Define number of SMCs to deploy
    smcnum = utils::inumeric(FLERR, arg[11], false, lmp);
    if (smcnum <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Define atom type of SMCs' ends
    smctype = utils::inumeric(FLERR, arg[12], false, lmp);
    if (smctype <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Define deployment bond type of SMCs
    smcbtype = utils::inumeric(FLERR, arg[13], false, lmp);
    if (smcbtype <= 0) error -> all(FLERR, "Illegal fix smc command");

    //Define angle type for patches
    atype=2;

    // Define deployment bond type of SMCs after the first movement
    smcbitype = utils::inumeric(FLERR, arg[14], false, lmp);
    if (smcbitype <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Define distance cutoff to be checked before SMCs' movement
    cutoff = utils::numeric(FLERR, arg[15], false, lmp);
    if (cutoff < 0)
      error -> all(FLERR, "Illegal cutoff value");

    // Define distance cutoff to be checked before SMCs' movement
    tancoff = utils::numeric(FLERR, arg[16], false, lmp);
    if (abs(tancoff) > 1)
      error -> all(FLERR, "Illegal tangent cutoff value");

    initmode = 0;

    // Define init mode for SMCs' deployment 
    if (strcmp(arg[17], "random") == 0) {
      initmode = 0;
    } else if (strcmp(arg[17], "distributed") == 0) {
      initmode = 1;
    } else if (strcmp(arg[17], "full-distributed") == 0) {
      initmode = 2;
    } else if (strcmp(arg[17], "fixed") == 0) {
      initmode = 3;
    } else {
      error -> all(FLERR, "Illegal fix smc command, initmode not present");
    }

    // Define loading probability of one SMC on simulation polymers
    kon = utils::numeric(FLERR, arg[18], false, lmp);
    if ((kon <= 0) || (kon>1))
      error -> all(FLERR, "Illegal fix smc command, kon is out of the interval (0,1]");

    // Define unloading probability of one SMC on simulation polymers
    koff = utils::numeric(FLERR, arg[19], false, lmp);
    if ((koff < 0) || (koff>1))
      error -> all(FLERR, "Illegal fix smc command, koff is out of the interval [0,1]");

    // Define number of patches per bead
    npatches = utils::inumeric(FLERR, arg[20], false, lmp);
    if (npatches < 0)
      error -> all(FLERR, "Illegal fix smc command, npatches is smaller than 0");

    debug = utils::numeric(FLERR, arg[21], false, lmp);

    int argnum = 22;

    if (initmode == 3){  
      if (narg<23) error -> all(FLERR, "Illegal fix smc command, smc startfile absent");
      fixFname = arg[22];
      argnum += 1;
    }

    // Define number of type of beads to avoid
    nblockt = narg - argnum;

    blockt = new int[nblockt];

    // Set the type of beads to avoid
    for (int i = argnum; i < narg; i++) {
      blockt[i-argnum] = utils::inumeric(FLERR, arg[i], false, lmp);;
    }

    anch = new long[smcnum];
    hing = new long[smcnum];

    // Initialise SMC positions
    if (initmode == 3) {
      // Read from the text file
      std::ifstream smcfile(fixFname);
      for (int i = 0; i < smcnum; i++) {
        smcfile >> anch[i] >> hing[i];

        if ((hing[i] <= 0) || (anch[i] <= 0) || (hing[i] > solsize) || (anch[i] > solsize))
        error -> all(FLERR, "Illegal position of hinge or anchor read from file");
      }
    }
    else {
      for (int i = 0; i < smcnum; i++) {
        hing[i] = -1;
        anch[i] = -1;
      }
    }

    // Initialise list of available positions
    av_list = new long[atom->natoms];
    num_avl = solsize / (4*(npatches+1));

    // Initialise seed generator
    random_equal = new RanPark(lmp, seed);

    // copy = special list for one atom
    // size = ms^2 + ms is sufficient
    // b/c in rebuild_special_one() neighs of all 1-2s are added,
    //   then a dedup(), then neighs of all 1-3s are added, then final dedup()
    // this means intermediate size cannot exceed ms^2 + ms

    int maxspecial = atom->maxspecial;
    copy = new tagint[maxspecial*maxspecial + maxspecial];

    // To get a different random number every time the program is executed
    srand(time(NULL) * seed);
  }

/* ---------------------------------------------------------------------- */

/**
 * @brief Destructor for FixSMC that deallocates all dynamically allocated memory.
 * 
 * **Cleanup Steps**:
 * 1. Delete random number generator (RanPark)
 * 2. Delete anchor position array (anch)
 * 3. Delete hinge position array (hing)
 * 4. Delete blocked bead type array (blockt)
 * 5. Delete special list working buffer (copy)
 * 6. Destroy availability list (av_list) using LAMMPS memory management
 * 
 * 
 * @see FixSMC::FixSMC()
 */
FixSMC::~FixSMC() { 

  // Deleting pointers
  delete random_equal;
  delete anch;
  delete hing;
  delete blockt;
  delete [] copy;
  
  memory -> destroy(av_list);

}

/* ---------------------------------------------------------------------- */

int FixSMC::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSMC::init() {

  // Check if pair and bond styles are initialised
  if (force -> pair == nullptr || force -> bond == nullptr)
    error -> all(FLERR, "Fix smc requires pair and bond styles");

  if (debug) {
    debugfile.open("log_fix_smc.txt", std::fstream::trunc | std::fstream::out);
    debugfile<<"#timestep,nproc,smcnum,prev_anch,prev_anch_type,prev_hing,prev_hing_type,anch,anch_type,hing,hing_type,rej_prob,rej_lim,rej_type,rej_dist,rej_tang,dist,tang"<<std::endl;
    debugfile.close();
  }

}

/* ---------------------------------------------------------------------- */

/**
 * @brief Main SMC loop extrusion engine that drives dynamics and movement every nevery timesteps.
 * 
 * This is the core function called after each integration timestep to perform SMC extrusion
 * dynamics. It handles two phases: initialization (timestep 1) and dynamics (every nevery steps).
 * During initialization, it loads SMC proteins based on kon probability. During dynamics, it
 * performs SMC loading/unloading, attempts extrusion movements, validates geometry constraints,
 * and updates topology.
 * 
 * **Phase 1: Initialization (timestep == 1)**:
 * 1. For each SMC i (0 to smcnum-1):
 *    - Draw random number (MPI_Bcast from rank 0)
 *    - If random < kon: Call load_smc(i) to place SMC
 * 2. Place loaded SMCs with place_smc() using newsmc=true
 * 3. Debug output for initial placement
 * 
 * **Phase 2: Dynamics (every nevery timesteps)**:
 * 1. Return if maxdir==0 (no movement allowed)
 * 2. For each SMC i (0 to smcnum-1):
 *    - **Load/Unload decision**:
 *      * If unloaded (anch[i]<0 or hing[i]<0):
 *        - If random < kon: load_smc(i) and place_smc() with newsmc=true
 *      * If loaded: If random < koff: remove_smc() and mark as unloaded
 *    - **Skip if unloaded**: Continue to next SMC if not loaded
 *    - **Movement proposal**:
 *      * Draw random number; if > prob: reject (debug_post with flag=1)
 *      * Determine movement directions (hdir, adir) based on dirmode:
 *        - dirmode=0 (bo): Bidirectional movement
 *        - dirmode=1 (mo1): Only hinge moves
 *        - dirmode=2 (mo2): Only anchor moves
 *        - dirmode=3 (mora): Alternating direction per SMC index
 *      * Apply random magnitude: hdir, adir ∈ [-(maxdir), +maxdir]
 *    - **Boundary checks**: Adjust directions if movement exceeds polymer bounds
 *      * Handle ring vs linear topology differently
 *      * Zero out blocked directions
 *    - **Overlap checks**: Verify no collision with other SMCs
 *      * Skip movement if full overlap detected
 *      * Reduce movement dimensions if partial overlap
 *    - **Blocked bead checks**: Verify atoms at new positions not blocked type
 *      * Use atom->map() to get local indices
 *      * MPI_Bcast results to all processors
 *    - **Geometry validation**:
 *      * compute_xyz(): Get unwrapped coordinates of anchor/hinge beads
 *      * Calculate distance between new anchor and hinge: dist = |r_anch - r_hing|
 *      * Calculate tangent product: tan = (old_vec · new_vec) / (|old| |new|)
 *      * Reject if dist > cutoff
 *      * Reject if tan < tancoff (angle constraint)
 *    - **Movement execution** (if all checks pass):
 *      * remove_smc(): Restore old anchor/hinge types and remove bond
 *      * Update anch[i] and hing[i] arrays
 *      * place_smc(): Set new atom types and create bond with newsmc=false
 * 3. MPI_Barrier after each SMC to synchronize processors
 * 
 * **Movement Direction Modes**:
 * - "bo" (bidirectional): Both ends can move independently
 * - "mo1" (anchor fixed): Only hinge moves
 * - "mo2" (hinge fixed): Only anchor moves
 * - "mora" (alternating): First half of smcs have anchor move, second half have hinge move
 * 
 * **Boundary Handling**:
 * - Linear polymers: Stop movement at boundaries (adir/hdir = 0)
 * - Ring polymers: Wrap movement 
 * 
 * **Error Flags in debug_post()**:
 * - err1: Movement rejected by probability (rand > prob)
 * - err2: Collision with other SMC
 * - err3: Blocked bead type encountered
 * - err4: Distance constraint violated (dist > cutoff)
 * - err5: Angle constraint violated (tan < tancoff)
 * 
 * **MPI Coordination**:
 * - MPI_Bcast: Random numbers, loaded SMC positions
 * - MPI_Allreduce: Global reduction for blocked bead checks
 * - MPI_Barrier: Synchronize processor states
 * 
 * @see FixSMC::load_smc(), FixSMC::place_smc(), FixSMC::remove_smc()
 * @see FixSMC::compute_xyz(), FixSMC::debug_pre(), FixSMC::debug_post()
 * @see FixSMC::check_avl(), FixSMC::compile_avl_list()
 */
void FixSMC::post_integrate() {

  if (update -> ntimestep == 1) {

    double lrand;

    // Load SMCs according to kon probability
    for (int i = 0; i < smcnum; i++)
    {
      if (comm -> me == 0) lrand = random_equal -> uniform();
      MPI_Bcast( & lrand, 1, MPI_DOUBLE, 0, world);
      if (lrand < kon) load_smc(i);
    }

    MPI_Barrier(world);
    
    for (int i = 0; i < smcnum; i++) {
      // Place the loaded SMCs
      if ((anch[i]<0) || (hing[i]<0)) continue;
      debug_pre(i);
      place_smc(anch[i],hing[i], true);
      MPI_Barrier(world);
      debug_post(i,0,0,0,0,0,0,0);
    }

    return;

  } 
  
  // Propose movements with a fixed frequency
  else if (update -> ntimestep % nevery == 0){

    // If movement is disabled along both directions stop execution
    if (maxdir == 0) return;

    int mnew;
    int mannew;

    double olddist = 0;
    double dist = 0;
    double tan = 0;

    std::array<double, 3> xyzanch, xyzhing, xyzoldanch, xyzoldhing; 

    auto histories = modify -> get_fix_by_style("BOND_HISTORY");
    int n_histories = histories.size();
    
    int idnewhi;
    int idhi;
    int idnewan;
    int idan;

    double jrand = 0;
    double krand = 0;

    bool flag = 0;

    int tempadir = 0;
    int temphdir = 0;

    double rand;
    
    for (int i = 0; i < smcnum; i++) {

      // Barrier to check that each processor has defined correctly each smc
      MPI_Barrier(world);
    
      // Draw two random numbers for the unloading/loading
      double lrand;

      if (comm -> me == 0) lrand = random_equal -> uniform();

      MPI_Bcast( & lrand, 1, MPI_DOUBLE, 0, world);

      // Load/Unload SMC according to the kon/koff probability
      if ((anch[i]<0) || (hing[i]<0)){
        if (lrand < kon) {
          load_smc(i);
          debug_pre(i);
          place_smc(anch[i],hing[i],true);
          MPI_Barrier(world);
          debug_post(i,0,0,0,0,0,0,0);
        }
      }
      else{
        if (lrand < koff) {
          debug_pre(i);
          remove_smc(anch[i],hing[i]);
          MPI_Barrier(world);
          debug_post(i,0,0,0,0,0,0,0);
          anch[i]=-1;
          hing[i]=-1;
        }
      }

      MPI_Barrier(world);

      // Check if SMC are loaded
      if ((anch[i]<0) || (hing[i]<0)) continue;

      // Draw a random number for the jump attempt
      if (comm -> me == 0) rand = random_equal -> uniform();
      MPI_Bcast( & rand, 1, MPI_DOUBLE, 0, world);

      // First accept the move with a certain probability prob
      if (rand > prob) {
        debug_pre(i);
        debug_post(i,1,0,0,0,0,0,0);
        continue;
      }

      if (maxdir>1){
        if (comm -> me == 0) jrand = random_equal -> uniform();
        MPI_Bcast( & jrand, 1, MPI_DOUBLE, 0, world);
        if (comm -> me == 0) krand = random_equal -> uniform();
        MPI_Bcast( & krand, 1, MPI_DOUBLE, 0, world);
      }

      if (maxdir == 0){
        adir = 0;
        hdir = 0;
      }
      else if (dirmode==0){
        hdir = (int)(maxdir * krand) + (maxdir) / (abs(maxdir));
        adir = -((int)(maxdir * jrand) + (maxdir) / (abs(maxdir)));
      }
      else if (dirmode==1){
        hdir = (int)(maxdir * krand) + (maxdir) / (abs(maxdir));
        adir = 0;
      }
      else if (dirmode==2){
        hdir = 0;
        adir = -((int)(maxdir * jrand) + (maxdir) / (abs(maxdir)));
      }
      else if ((dirmode==3) && (i<smcnum/2)){
        hdir = (int)(maxdir * krand) + (maxdir) / (abs(maxdir));
        adir = 0;
      }
      else if ((dirmode==3) && (i>=smcnum/2)){
        hdir = 0;
        adir = -((int)(maxdir * jrand) + (maxdir) / (abs(maxdir)));
      }

      // Temporary direction if the smc is going towards the polymer end or another smc bead
      tempadir = adir;
      temphdir = hdir;

      // Check if the SMCs' ends are moving over the polymer ends
      if (hdir < 0) {
        if (((hing[i] + hdir) % (lpol) <= 0) || (((hing[i] + hdir) / (lpol)) != ((hing[i]) / (lpol))) || ((hing[i] + hdir) > solsize)){
          if (!ring) temphdir = 0;
          else temphdir = (hing[i] + hdir) % (lpol) - hing[i];
        }
      } else {
        if (((hing[i] + hdir) / (lpol)) != ((hing[i]) / (lpol)) || ((hing[i] + hdir) > solsize)){
          if (!ring) temphdir = 0;
          else temphdir = (hing[i] + hdir) % (lpol) - hing[i];
        }
      }
      if (adir < 0) {
        if (((anch[i] + adir) % (lpol) <= 0) || (((anch[i] + adir) / (lpol)) != ((anch[i]) / (lpol))) || ((anch[i] + adir) > solsize)){
          if (!ring) tempadir = 0;
          else tempadir = (anch[i] + hdir) % (lpol) - anch[i];
        }
      } else {
        if (((anch[i] + adir) / lpol) != ((anch[i]) / lpol) || ((anch[i] + adir) > solsize)){
          if (!ring) tempadir = 0;
          else tempadir = (anch[i] + hdir) % (lpol) - anch[i];
        }
      }

      // Skip movement if both the ends are not moving
      if ((tempadir == 0) && (temphdir == 0)) {
        debug_pre(i);
        debug_post(i,1,0,0,0,0,0,0);
        continue;
      }

      // Check if the new movement is forbidden because of superposition of SMCs
      flag = 0;
      for (int j = 0; j < smcnum; j++) {
        if (j == i){
          if (((anch[i] + tempadir) == (hing[i] + temphdir)) || ((anch[i] + tempadir) == hing[i]) || (anch[i] == (hing[i] + temphdir))){
            flag=1;
            break;
          }
        } 
        if ((((anch[i] + tempadir) == anch[j]) || ((anch[i] + tempadir) == hing[j])) && (((hing[i] + temphdir) == hing[j]) || ((hing[i] + temphdir) == anch[j]))) {
          flag = 1;
          break;
        } else if (((anch[i] + tempadir) == anch[j]) || ((anch[i] + tempadir) == hing[j])) tempadir = 0;
        else if (((hing[i] + temphdir) == hing[j]) || ((hing[i] + temphdir) == anch[j])) temphdir = 0;
      }

      // Interrup internal l_sample loop
      if (flag) {
        debug_pre(i);
        debug_post(i,0,1,0,0,0,0,0);
        continue;
      }
      
      // Finding local identifiers of the current hinge and anchor
      idan = atom -> map(map_to_beads(anch[i]));
      idhi = atom -> map(map_to_beads(hing[i]));
    
      // Finding local identifiers of the proposed new hinge and anchor
      idnewhi = atom -> map(map_to_beads(hing[i] + temphdir));
      idnewan = atom -> map(map_to_beads(anch[i] + tempadir));

      for (int j = 0; j < nblockt ; j++)
      {
          if ((((idnewan) >= 0) && (idnewan < atom -> nlocal)) && atom -> type[idnewan] == blockt[j]) {tempadir = 0; MPI_Bcast(&tempadir,1,MPI_INT,comm->me,world);}
          if ((((idnewhi) >= 0) && (idnewhi < atom -> nlocal)) && atom -> type[idnewhi] == blockt[j]) {temphdir = 0; MPI_Bcast(&temphdir,1,MPI_INT,comm->me,world);}
      }

      MPI_Barrier(world);
      if ((tempadir == 0) && (temphdir == 0)) {
        debug_pre(i);
        debug_post(i,0,0,1,0,0,0,0);
        continue;
      }

      xyzoldanch = compute_xyz(map_to_beads(anch[i]));
      xyzoldhing = compute_xyz(map_to_beads(hing[i]));
      xyzanch = compute_xyz(map_to_beads(anch[i]+tempadir));
      xyzhing = compute_xyz(map_to_beads(hing[i]+temphdir));

      olddist = 0;
      dist = 0;
      tan = 0;

      //Compute distance and tangent on old and new set of beads
      for (int k = 0; k < 3; k++) {
        olddist += (xyzoldanch[k] - xyzoldhing[k]) * (xyzoldanch[k] - xyzoldhing[k]);
        dist += (xyzanch[k] - xyzhing[k]) * (xyzanch[k] - xyzhing[k]);
        tan += (xyzoldanch[k] - xyzoldhing[k]) * (xyzanch[k] - xyzhing[k]);
      }

      olddist = sqrt(olddist);
      dist = sqrt(dist);
      tan /= (olddist*dist);

      // Barrier to check that each processor has defined correctly each smc
      MPI_Barrier(world);

      // Check if the distance is small enough to accept the movement
      if (dist > cutoff) {
        // Debug print if got rejected because of distance between beads
        debug_pre(i);
        debug_post(i,0,0,0,1,0,dist,tan);
        continue;
      }

      // Check if the product between the two tangents is smaller than a cutoff
      if (tan < tancoff) {
        // Debug print if got rejected because of tangent product value
        debug_pre(i);
        debug_post(i,0,0,0,0,1,dist,tan);
        continue;
      }

    // If found a suitable new pair of beads and move the bond
    remove_smc(anch[i], hing[i]);
    debug_pre(i);

    place_smc((anch[i] + tempadir), (hing[i] + temphdir), false);
    anch[i] += tempadir;
    hing[i] += temphdir;

    debug_post(i,0,0,0,0,0,dist,tan);

    // Barrier to check that each processor has defined correctly each smc
    MPI_Barrier(world);

    } 

    return;

  }

}

/*----------------------------------------------------*/

/**
 * @brief Converts logical bead index to actual atom ID accounting for patch atoms.
 * 
 * This function maps a logical polymer bead position i (used in SMC code) to the
 * actual LAMMPS atom ID. Since each bead can have attached patch atoms for visualization,
 * the mapping accounts for the gap created by patches.
 * 
 * **Mapping Formula**:
 * - atom_id = (i - 1) × (npatches + 1) + 1
 * - This shifts bead indices by (npatches + 1) positions to account for patch atoms
 * 
 * **Example** (npatches=1):
 * - Logical bead 1 → Atom 1 (main bead)
 * - Logical bead 2 → Atom 3 (skipping atom 2, which is patch for bead 1)
 * - Logical bead 3 → Atom 5 (skipping atom 4, which is patch for bead 2)
 * 
 * **Special Considerations**:
 * - Essential for correct SMC positioning when using patch atoms
 * - npatches determines spacing between logical beads in atom numbering
 * - Inverse function: logical_i = (atom_id - 1) / (npatches + 1) + 1
 * 
 * @param i [long] Logical polymer bead index
 * @return [long] Actual LAMMPS atom ID
 * @see FixSMC::place_smc(), FixSMC::remove_smc(), FixSMC::check_avl()
 */
long FixSMC::map_to_beads(long i){
  return (i-1)*(npatches+1) + 1;
}

/*----------------------------------------------------*/

/**
 * @brief Verifies if an anchor can be placed at bead position i (geometric validation).
 * 
 * This function checks whether a new SMC anchor can be placed at position i by verifying:
 * 1. Border conditions are satisfied
 * 2. No overlap with existing SMC positions (anchor, hinge, or adjacent beads)
 * 3. No overlap with blocked bead types
 * 
 * The function validates based on the anchor-hinge configuration where the hinge would be
 * placed at i + 2 relative to the anchor.
 * 
 * **Validation Steps**:
 * 1. Calculate hinge position: tmphing = i + 2
 * 2. Get local atom indices via atom->map():
 *    - idnewan: New anchor atom
 *    - idnewhi: New hinge atom
 *    - idnewmid: Center bead atom
 * 3. Calculate midpoint: mdbead = (i + tmphing)/2
 * 4. **Border check**: Return false if:
 *    - mdbead % lpol == 1 (near start of polymer)
 *    - mdbead == 1 (very start)
 *    - mdbead % lpol == 0 (near end of polymer)
 * 5. **SMC overlap check**: Return false if any existing SMC j occupies:
 *    - Same anchor/hinge as proposed (i or tmphing)
 *    - Neighbors of proposed positions (i±1, tmphing±1)
 * 6. **Blocked bead check**: For each blocked bead type:
 *    - Check anchor atom type against all blocked types
 *    - Check hinge atom type against all blocked types
 *    - Check midpoint atom type against all blocked types
 *    - MPI_Bcast results to ensure global consistency
 * 7. Return true if all checks pass
 * 
 * **Special Considerations**:
 * - Used during load_smc() to find valid starting positions
 * - Enforces minimum spacing from polymer boundaries
 * - Prevents SMC overlap or stacking
 * - Avoids binding to forbidden atom types (e.g., crosslinks)
 * - MPI coordination ensures all processors agree on validity
 * 
 * @param i [long] Logical anchor position to test
 * @return [bool] True if anchor can be placed, false if blocked
 * @see FixSMC::load_smc(), FixSMC::compile_avl_list(), FixSMC::check_avl()
 */
bool FixSMC::check_avl(long i){
  long tmphing;
  long mdbead;
  bool flag = 0;

  int mnew;
  int mannew;
  int midnew;

  int idnewhi;
  int idnewan;
  int idnewmid;

  // Define hinge position according to movement direction
  tmphing = i + 2 * maxdir / (abs(maxdir));

  idnewhi = atom -> map(map_to_beads(tmphing));
  idnewan = atom -> map(map_to_beads(i));
  idnewmid = atom -> map(map_to_beads((i+tmphing)/2));

  // Define SMC's center
  mdbead = (i + tmphing)/2;

  // Check if the smc is wrongly positioned (border conditions)
  if ((mdbead%lpol == 1) || (mdbead == 1) || (mdbead%lpol == 0)) {
    return 0;
  }

  flag = 0;
  // Check if we are moving over SMCs' beads
  for (int j = 0; j < smcnum; j++) {
    if ((anch[j]<0) || (hing[j]<0)) continue;
    if (((i == anch[j]) || (tmphing == hing[j])) || ((i == hing[j]) || (tmphing == anch[j])) || ((i + 1) == anch[j]) || ((i - 1) == anch[j]) || ((tmphing + 1) == anch[j]) || ((tmphing - 1) == anch[j]) || ((i + 1) == hing[j]) || ((i - 1) == hing[j]) || ((tmphing + 1) == hing[j]) || ((tmphing - 1) == hing[j])) {
      flag = 1;
      break;
    }
  }

  if (flag) return 0;
  
  // Check if we are moving over not available beads
  for (int j = 0; j < nblockt ; j++)
    {
      if ((((mannew = idnewan) >= 0) && (idnewan < atom -> nlocal)) && atom -> type[mannew] == blockt[j]) {
        flag=1; 
        MPI_Bcast(&flag,1,MPI_INT,comm->me,world); 
        break;
      }
      if ((((midnew = idnewmid) >= 0) && (idnewmid < atom -> nlocal)) && atom -> type[midnew] == blockt[j]) {
        flag=1; 
        MPI_Bcast(&flag,1,MPI_INT,comm->me,world); 
        break;
      }
      if ((((mnew = idnewhi) >= 0) && (idnewhi < atom -> nlocal)) && atom -> type[mnew] == blockt[j]) {
        flag=1; 
        MPI_Bcast(&flag,1,MPI_INT,comm->me,world); 
        break;
      }
    }

  if (flag) return 0;

  return 1;
}

/*----------------------------------------------------*/

/**
 * @brief Compiles a globally synchronized list of available bead positions for SMC loading.
 * 
 * This function builds a comprehensive list of all unoccupied bead positions that can
 * accommodate new SMC proteins. It combines local availability checks with MPI synchronization
 * to ensure all processors have the same availability information, which is critical for
 * reproducible SMC loading across parallel simulations.
 * 
 * **Algorithm Steps**:
 * 1. Initialize available positions count: num_avl = 0
 * 2. Initialize temporary availability array: temp_avl_list[atom->nlocal]
 * 3. **Compute MPI displacement arrays**:
 *    - MPI_Allgather to get rcounts: nlocal atoms per processor
 *    - Calculate displs: cumulative byte offsets for Allgatherv
 * 4. **Local availability check**:
 *    - For each bead i from 1 to solsize/(npatches+1):
 *      * If bead is local (atom->map returns valid index):
 *        - Call check_avl(i) for geometric validation
 *        - If valid: add to temp_avl_list, increment temp_num_avl
 * 5. **Global reduction**: MPI_Allreduce(temp_num_avl, num_avl, SUM)
 *    - Combines counts from all processors
 * 6. **Gather availability lists**:
 *    - MPI_Allgatherv to collect temp_avl_list from all processors
 *    - Results stored in av_list at offsets specified by displs
 * 7. **Sort results**: std::sort(av_list, av_list + atom->natoms)
 *    - Ensures consistent ordering for reproducible random selection
 * 
 * **Data Structures**:
 * - temp_avl_list: Local list of available positions (size: atom->nlocal)
 * - av_list: Global list of available positions (all processors, size: atom->natoms)
 * - rcounts: Number of atoms on each processor (size: comm->nprocs)
 * - displs: Byte offsets for Allgatherv communication (size: comm->nprocs)
 * 
 * **MPI Barriers**:
 * - Called before Allgatherv to ensure all processors reach same synchronization point
 * 
 * **Special Considerations**:
 * - Uses std::fill and std::sort from C++ standard library
 * - Allocates dynamic array: delete[] temp_avl_list at end
 * - Called at beginning of load_smc() and every nevery timesteps in post_integrate()
 * - Result av_list sorted to ensure deterministic random selection
 * - Performance: O(lpol) for local checks + O(lpol log lpol) for sorting + MPI communication
 * 
 * @see FixSMC::check_avl(), FixSMC::load_smc(), FixSMC::post_integrate()
 * @throws error->all if num_avl == 0 (insufficient space for SMCs)
 */
void FixSMC::compile_avl_list(){
  num_avl = 0;
  
  // Getting size of array to compute
  long temp_num_avl = 0;
  long * temp_avl_list = new long[atom->nlocal];

  // Filling the list with bead number solsize + 1
  std::fill(temp_avl_list, temp_avl_list + atom->nlocal, solsize/(npatches+1) + 1);

  // Define displacements for MPI parallelisation
  int displs[comm->nprocs];
  int rcounts[comm->nprocs];

  // Filling the list with bead number solsize + 1
  std::fill(displs, displs + comm->nprocs, 0);

  MPI_Allgather(&atom->nlocal, 1, MPI_INT, rcounts, 1, MPI_INT, world);

  for (int i = 1; i < comm->nprocs; i++)
  {
    displs[i] += rcounts[i-1] + displs[i-1];
  }

  for (int i = 1; i <= solsize/(npatches+1); i++)
  {
    if (((atom->map(map_to_beads(i))) >= 0) && (atom->map(map_to_beads(i)) < atom -> nlocal)){
      if (check_avl(i)) {
      temp_avl_list[temp_num_avl] = i;
      temp_num_avl++;
      }
    }
  }

  // Gather back available positions list
  MPI_Barrier(world);

  MPI_Allreduce(&temp_num_avl, &num_avl, 1, MPI_LONG, MPI_SUM, world);

  MPI_Allgatherv(temp_avl_list, atom->nlocal, MPI_LONG, av_list, rcounts, displs, MPI_LONG, world);

  // Sort the list to have a working availability list
  std::sort(av_list, av_list + atom->natoms);
  
  delete[] temp_avl_list;
}

/*----------------------------------------------------*/

/**
 * @brief Loads (initializes) an SMC protein at a valid position with specified initialization mode.
 * 
 * This function creates a new SMC protein by selecting an anchor
 * position and setting the corresponding hinge position. The selection strategy depends on
 * the initialization mode (initmode):
 * - **Mode 1 (distributed)**: One SMC per polymer, placed at random position within polymer
 * - **Mode 2 (full-distributed)**: SMCs distributed among polymers sequentially
 * - **Mode 3 (fixed)**: Load from input file (this function returns early)
 * - **Mode 4 (random)**: Random placement from available positions list
 * 
 * **Algorithm Steps**:
 * 1. **Compile availability list**: Call compile_avl_list() to get all valid positions
 * 2. **Check polymer count**: npol = solsize / (lpol × (npatches+1))
 * 3. **Select anchor position** based on initmode:
 *    - **Mode 1 + i < npol + timestep==1**:
 *      * Generate random position: anch[i] = random() × lpol + i × lpol
 *      * Verify validity: loop until check_avl(anch[i]) returns true
 *      * Ensures one SMC per polymer (first nevery timesteps)
 *    - **Mode 2 + timestep==1**:
 *      * Generate position across multiple polymers: anch[i] = random() × lpol + (i%npol) × lpol
 *      * Verify validity with loop
 *    - **Mode 3 + timestep==1**: Return early (load from restart file)
 *    - **Mode 4** (all other cases):
 *      * Select from available list: anch[i] = av_list[random() × num_avl]
 *      * Uses pre-compiled list for faster repeated calls
 * 4. **Validate availability**: Check num_avl > 0, error if no valid positions
 * 5. **Set hinge position**: hing[i] = anch[i] + 2 (always 2 beads from anchor)
 * 6. **Broadcast positions**:
 *    - MPI_Bcast(anch, smcnum, ..., 0, world) to distribute from rank-0
 *    - MPI_Bcast(hing, smcnum, ..., 0, world) to distribute from rank-0
 * 7. **Synchronize**: MPI_Barrier ensures all processors have updated SMC data
 * 
 * **MPI Coordination**:
 * - Random number generation happens only on rank-0 (comm->me == 0)
 * - All processors broadcast values to ensure consistency
 * - Barriers ensure synchronized state after loading
 * 
 * **Hinge Position Rules**:
 * - Always placed 2 beads away from anchor
 * - For bidirectional: hing = anch + 2 or anch - 2 (determined by movement)
 * - Ensures minimum SMC size and prevents overlap
 * 
 * **Special Considerations**:
 * - Called from post_integrate() during loading phase (when random() < kon)
 * - Only called on rank-0 for actual random selection (other processors receive via Bcast)
 * - Mode 3 early return: assumes SMCs already loaded from restart file
 * - Fails if num_avl == 0: all positions occupied or invalid
 * - Performance: O(check_avl) per attempt in modes 1-2, O(1) in mode 4
 * 
 * @param i [long] SMC index to load (0 to nsmc-1)
 * @see FixSMC::compile_avl_list(), FixSMC::check_avl(), FixSMC::post_integrate()
 * @throws error->all(FLERR, ...) if num_avl == 0 (insufficient space)
 */
void FixSMC::load_smc(long i) {
  compile_avl_list();

  int npol = solsize / (lpol*(npatches+1)); 

  // Distribute the SMCs in the system, placing at least one for each polymer
  if ((initmode == 1) && (i < npol) && (update -> ntimestep  == 1)) {
    do {
      if (comm->me == 0) anch[i] = static_cast < int > (random_equal -> uniform() * lpol + i * lpol);
      MPI_Bcast(anch, smcnum, MPI_LONG, 0, world);
      MPI_Barrier(world);
    } while (! check_avl(anch[i]));
  }
  // Distribute uniformly SMCs among the polymers
  else if ((initmode == 2) && (update -> ntimestep  == 1)) {
    do {
      if (comm->me == 0) anch[i] = static_cast < int > (random_equal -> uniform() * lpol + (i%npol) * lpol);
      MPI_Bcast(anch, smcnum, MPI_LONG, 0, world);
      MPI_Barrier(world);
    } while (! check_avl(anch[i]));
  }
  // Skip since loaded from file
  else if ((initmode == 3) && (update -> ntimestep  == 1)) {
    return;
  }
  // Distribute randomly extruders among polymers
  else{
    if (comm->me == 0) anch[i] = av_list[static_cast < int > (random_equal -> uniform() * num_avl)];
  }  

  if (num_avl == 0) error -> all(FLERR, "Not enough space for the smcs");

  hing[i] = anch[i] + 2;

  // Cast the chosen position to each processor
  MPI_Bcast(anch, smcnum, MPI_LONG, 0, world);
  MPI_Bcast(hing, smcnum, MPI_LONG, 0, world);

  MPI_Barrier(world);
}

/**
 * @brief Creates a bond between two atoms and updates the special neighbor list.
 * 
 * This function establishes a bond between two atoms
 * in the simulation. It handles both the bond list storage and the special neighbor list
 * (1-2, 1-3, 1-4 nearest neighbor relations) which are critical for excluded volume
 * interactions in LAMMPS.
 * 
 * **Algorithm Steps**:
 * 1. **Map atom tags to local indices**:
 *    - atom_map1 = atom->map(atom1)
 *    - atom_map2 = atom->map(atom2)
 *    - Local indices only valid if in range [0, nlocal)
 * 2. **Determine bond storage** (Newton's 3rd law symmetry):
 *    - If newton_bond = false OR atom1 < atom2:
 *      * Store bond in atom1's bond_atom and bond_type arrays
 *      * Increment atom1's num_bond counter
 *      * Increment global atom->nbonds counter
 *      * Set reneighboring flag for next timestep
 *    - Otherwise: Bond stored only by atom2 (not both)
 * 3. **Update special neighbor list** for atom1:
 *    - Get special list pointers: slist, n1, n2, n3
 *    - Check if atom2 already in special list (m from n1 to n3)
 *    - If found: remove duplicate (shift down and decrement counts)
 *    - If n3 = maxspecial: error (list full)
 *    - Insert atom2 at position n1 (make it a 1-2 neighbor)
 *    - Update counts: n1++, n2++, n3++
 * 4. **Error checking**:
 *    - Verify num_bond[atom_map1] < bond_per_atom
 *    - Verify n3 < maxspecial after insertion
 * 
 * **Special Considerations**:
 * - Only operates on local atoms (atom_map1 must be in [0, nlocal))
 * - Newton's 3rd law pairs handled by force->newton_bond setting
 * - Triggers reneighboring if new bond created (next_reneighbor = ntimestep)
 * - Special list must not exceed maxspecial entries
 * - Duplicate removal prevents re-adding existing bonds
 * - Used during SMC movement to connect anchor and hinge
 * - Called twice per movement: once for each direction (1→2 and 2→1)
 * 
 * @param atom1 [long] Tag of first atom in bond
 * @param atom2 [long] Tag of second atom in bond
 * @param btype [int] Bond type (smcbitype for new SMC, smcbtype for moving SMC)
 * @see FixSMC::place_smc(), FixSMC::post_integrate(), FixSMC::remove_bond()
 * @throws error->one if bond exceeds bonds_per_atom limit
 * @throws error->one if special list exceeds maxspecial limit
 */
void FixSMC::create_bond(long atom1, long atom2, int btype){

  int atom_map1 = atom->map(atom1);
  int atom_map2 = atom->map(atom2);

  tagint *slist;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  int i,j,k,m,n,ii,jj,inum,jnum,itype,jtype,n1,n2,n3,possible;
  int *type = atom->type;

  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  // if newton_bond is set, only store with atom2
  // if not newton_bond, store bond with both atom1 and atom2

  if (!force->newton_bond || atom1 < atom2) {
    if ((atom_map1 >= 0) && (atom_map1 < atom->nlocal)) {
      if (num_bond[atom_map1] == atom->bond_per_atom)
        error->one(FLERR,"New bond exceeded bonds per atom in fix bond/create");
      bond_type[atom_map1][num_bond[atom_map1]] = btype;
      bond_atom[atom_map1][num_bond[atom_map1]] = atom2;
      num_bond[atom_map1]++;
    }

    atom->nbonds++;

    next_reneighbor = update->ntimestep;
  }

  if ((atom_map1 >= 0) && (atom_map1 < atom->nlocal)){
    // add a 1-2 neighbor to special bond list
   
    slist = special[atom_map1];
    n1 = nspecial[atom_map1][0];
    n2 = nspecial[atom_map1][1];
    n3 = nspecial[atom_map1][2];
    for (m = n1; m < n3; m++)
      if (slist[m] == atom2) break;
    if (m < n3) {
      for (n = m; n < n3-1; n++) slist[n] = slist[n+1];
      n3--;
      if (m < n2) n2--;
    }
    if (n3 == atom->maxspecial)
      error->one(FLERR,
                "New bond exceeded special list size in fix smc");
    for (m = n3; m > n1; m--) slist[m] = slist[m-1];
    slist[n1] = atom2;
    nspecial[atom_map1][0] = n1+1;
    nspecial[atom_map1][1] = n2+1;
    nspecial[atom_map1][2] = n3+1;

  }

}

/**
 * @brief Removes a bond between two atoms and updates the special neighbor list.
 * 
 * This function breaks a chemical bond (or SMC cohesin connection) between two atoms.
 * It handles both removal from the bond list and cleanup from the special neighbor list.
 * Additionally, it manages any bond history information (for bond_history fix) that needs
 * to be cleaned up.
 * 
 * **Algorithm Steps**:
 * 1. **Identify bond history fixes**:
 *    - Query modify->get_fix_by_style("BOND_HISTORY") for all history managers
 *    - n_histories = count of such fixes
 * 2. **Determine bond storage location** (Newton's 3rd law):
 *    - If newton_bond = false OR atom1 < atom2:
 *      * Bond stored in atom1's lists; proceed to deletion
 *    - Otherwise: Bond stored elsewhere, skip this phase
 * 3. **Delete bond from atom1** (if applicable):
 *    - Map atom1 tag to local index
 *    - Search atom1's bond_atom array for atom2
 *    - When found at index m:
 *      * Shift all bonds after m down by one position
 *      * If bond_history present: shift history data accordingly
 *      * Decrement num_bond[atom1]
 *      * Decrement global atom->nbonds
 *      * Set reneighboring flag: next_reneighbor = ntimestep
 * 4. **Clean up bond history** (if applicable):
 *    - After shifting, delete history at position (num_bond-1)
 *    - Calls FixBondHistory->delete_history() for all history fixes
 * 5. **Remove atom2 from special neighbor list** of atom1:
 *    - Get special list pointers: slist, n1, n3
 *    - Search 1-2 neighbors (m from 0 to n1)
 *    - When atom2 found at position m:
 *      * Shift all atoms after m down by one position
 *      * Decrement all three counts: n1--, n2--, n3--
 * 6. **Global synchronization**:
 *    - Triggers reneighboring across all processors
 *  
 * @param atom1 [long] Tag of first atom in bond to remove
 * @param atom2 [long] Tag of second atom in bond to remove
 * @see FixSMC::remove_smc(), FixSMC::post_integrate(), FixSMC::create_bond()
*/
void FixSMC::remove_bond(long atom1, long atom2){
  
  int i,j,k,m,n,i1,i2,n1,n3,type;

  // find instances of bond history to delete data
  auto histories = modify->get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  // break bonds

  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  tagint *slist;

  int atom_map1 = atom->map(atom1);
  int atom_map2 = atom->map(atom2);

  // delete bond from atom1 if atom1 stores it

  if (!force->newton_bond || atom1 < atom2) {
    if ((atom_map1 >= 0) && (atom_map1 < atom->nlocal)) {
      for (m = 0; m < num_bond[atom_map1]; m++) {
        if (bond_atom[atom_map1][m] == atom2) {
          for (k = m; k < num_bond[atom_map1]-1; k++) {
            bond_atom[atom_map1][k] = bond_atom[atom_map1][k+1];
            bond_type[atom_map1][k] = bond_type[atom_map1][k+1];
            if (n_histories > 0)
              for (auto &ihistory: histories)
                dynamic_cast<FixBondHistory *>(ihistory)->shift_history(atom_map1,k,k+1);
          }
          if (n_histories > 0)
            for (auto &ihistory: histories)
              dynamic_cast<FixBondHistory *>(ihistory)->delete_history(atom_map1,num_bond[atom_map1]-1);
          num_bond[atom_map1]--;
          break;
        }
      }
    }

    atom->nbonds--;

    next_reneighbor = update->ntimestep;
  }

  // remove atom2 from special bond list for atom atom1

  if ((atom_map1 >= 0) && (atom_map1 < atom->nlocal)){

    slist = special[atom_map1];
    n1 = nspecial[atom_map1][0];
    for (m = 0; m < n1; m++)
      if (slist[m] == atom2) break;
    n3 = nspecial[atom_map1][2];
    for (; m < n3-1; m++) slist[m] = slist[m+1];
    nspecial[atom_map1][0]--;
    nspecial[atom_map1][1]--;
    nspecial[atom_map1][2]--;

  }
}


/**
 * @brief Visually and topologically marks an SMC by changing atom types and creating bonds.
 * 
 * This function transitions an SMC from an inactive state to an active state by:
 * 1. Changing the atom types of anchor and hinge beads (visual coloring for VMD)
 * 2. Creating bonds between anchor and hinge
 * 3. Updating the special neighbor list (1-2, 1-3, 1-4 relations)
 * 4. Recoloring patch atoms with distinct types for visualization
 * 
 * The function uses different bond types depending on whether this is a newly loaded SMC
 * (newsmc=true, uses smcbitype) or a moving SMC (newsmc=false, uses smcbtype).
 * 
 * **Algorithm Steps**:
 * 1. **Map logical positions to local atom indices**:
 *    - idhi = atom->map(map_to_beads(h)) - Hinge position
 *    - idan = atom->map(map_to_beads(a)) - Anchor position
 * 2. **Change anchor type**:
 *    - If anchor is local (idan in [0, nlocal)):
 *      * atom->type[idan] = smctype
 * 3. **Change hinge type**:
 *    - If hinge is local (idhi in [0, nlocal)):
 *      * atom->type[idhi] = smctype + 2
 * 4. **Create anchor-hinge bond**:
 *    - If newsmc (newly loaded):
 *      * create_bond(h, a, smcbitype) - Initial SMC bond
 *      * create_bond(a, h, smcbitype) - Reverse direction
 *    - Else (moving SMC):
 *      * create_bond(h, a, smcbtype) - Moving SMC bond
 *      * create_bond(a, h, smcbtype) - Reverse direction
 * 5. **Update special neighbor lists**:
 *    - Call update_topology(h, a) to rebuild 1-2, 1-3, 1-4 lists
 * 6. **Recolor patch atoms**:
 *    - For each patch c from 1 to npatches:
 *      * Anchor patch: type = smctype + 1
 *      * Hinge patch: type = smctype + 3
 * 
 * **Atom Type Encoding**:
 * - Inactive polymer bead: type = 1
 * - Active anchor: type = smctype
 * - Active hinge: type = smctype + 2
 * - Anchor patch: type = smctype + 1
 * - Hinge patch: type = smctype + 3
 * 
 * **Bond Types**:
 * - smcbitype: Initial SMC bond (newly loaded SMC)
 * - smcbtype: Movement bond (moving SMC updating its anchor/hinge)
 * - Allows tracking when bonds are created vs updated in post-processing
 * 
 * @param a [long] Logical anchor bead position
 * @param h [long] Logical hinge bead position
 * @param newsmc [bool] True if newly loaded SMC (uses smcbitype), false if moving (uses smcbtype)
 * @see FixSMC::remove_smc(), FixSMC::create_bond(), FixSMC::update_topology(), FixSMC::post_integrate()
 */
void FixSMC::place_smc(long a, long h, bool newsmc) {
  long mhi;
  long man;

  long idhi = atom -> map(map_to_beads(h));
  long idan = atom -> map(map_to_beads(a));

  // Changing type of the new anchor
  if (((man = idan) >= 0) && (idan < atom -> nlocal)) {
    atom -> type[man] = smctype;
  }

  // Create new bond between new hinge and anchor if not already present
  if (((mhi = idhi) >= 0) && (idhi < atom -> nlocal)) {
    // Changing type of new hing
    atom -> type[mhi] = smctype + 2;
  }

  // Creating new SMC bond
  if (newsmc){
    create_bond(map_to_beads(h), map_to_beads(a), smcbitype);
    create_bond(map_to_beads(a), map_to_beads(h), smcbitype);

  }
  else {
    create_bond(map_to_beads(h), map_to_beads(a), smcbtype);
    create_bond(map_to_beads(a), map_to_beads(h), smcbtype);
  }

  // Update special neighbors topology
  update_topology(map_to_beads(h), map_to_beads(a));

  // Recoloring patches
  for (int c = 1; c <= npatches; c++)
  {
    long idhipc = atom -> map(map_to_beads(h) + c);
    long idanpc = atom -> map(map_to_beads(a) + c);

    if (((man = idanpc) >= 0) && (idanpc < atom -> nlocal)) {
      atom -> type[man] = smctype + 1;
    }

    if (((mhi = idhipc) >= 0) && (idhipc < atom -> nlocal)) {
      atom -> type[mhi] = smctype + 3;
    }
  }

}

/**
 * @brief Restores SMC atoms to polymer state by removing bonds and resetting types.
 * 
 * This function transitions an SMC from an active state to an inactive state by:
 * 1. Resetting anchor and hinge atom types back to inactive polymer (type 1)
 * 2. Removing bonds between anchor and hinge
 * 3. Updating the special neighbor list (removing 1-2 relations)
 * 4. Restoring patch atom types to inactive polymer
 * 
 * This is the inverse operation of place_smc() and is called when an SMC unloads
 * from the polymer during the unloading phase of the SMC dynamics.
 * 
 * **Algorithm Steps**:
 * 1. **Map logical positions to local atom indices**:
 *    - idhi = atom->map(map_to_beads(h)) - Hinge position
 *    - idan = atom->map(map_to_beads(a)) - Anchor position
 * 2. **Reset anchor type**:
 *    - If anchor is local (idan in [0, nlocal)):
 *      * atom->type[idan] = 1 (back to polymer)
 * 3. **Reset hinge type**:
 *    - If hinge is local (idhi in [0, nlocal)):
 *      * atom->type[idhi] = 1 (back to polymer)
 * 4. **Remove anchor-hinge bond**:
 *    - remove_bond(h, a) - Hinge to anchor
 *    - remove_bond(a, h) - Anchor to hinge
 *    - Triggers compaction of bond arrays and special list
 * 5. **Update special neighbor lists**:
 *    - Call update_topology(h, a) to rebuild 1-2, 1-3, 1-4 lists
 * 6. **Reset patch atom types**:
 *    - For each patch c from 1 to npatches:
 *      * Anchor patch: type = 1 (back to polymer)
 *      * Hinge patch: type = 1 (back to polymer)
 * 
 * **Type Restoration**:
 * All atom types return to 1 (inactive polymer bead)
 * 
 * **Special Considerations**:
 * - Inverse of place_smc() function
 * - Only updates local atoms
 * - Bond removal triggers special list cleanup via remove_bond()
 * - Used during SMC unloading (when random() < koff)
 * - Called from post_integrate() in unloading phase
 * - Patch restoration helps reduce visualization clutter
 * - No MPI communication required (type changes are local)
 * - Must properly clean up bonds before type change to avoid inconsistencies
 * 
 * @param a [long] Logical anchor bead position
 * @param h [long] Logical hinge bead position
 * @see FixSMC::place_smc(), FixSMC::remove_bond(), FixSMC::update_topology(), FixSMC::post_integrate()
 */
void FixSMC::remove_smc(long a, long h) {
  long mhi;
  long man;

  long idhi = atom -> map(map_to_beads(h));
  long idan = atom -> map(map_to_beads(a));

  // Changing type of the old anchor
  if (((man = idan) >= 0) && (idan < atom -> nlocal)) {
    atom -> type[man] = 1;
  }

  // Changing type of the old hinge
  if (((mhi = idhi) >= 0) && (idhi < atom -> nlocal)) {
    atom -> type[mhi] = 1;
  }

  remove_bond(map_to_beads(h),map_to_beads(a));
  remove_bond(map_to_beads(a),map_to_beads(h));

  // Update special neighbors topology
  update_topology(map_to_beads(h), map_to_beads(a));

  // Recoloring patches
  for (int c = 1; c <= npatches; c++)
  {
    long idhipc = atom -> map(map_to_beads(h) + c);
    long idanpc = atom -> map(map_to_beads(a) + c);
    
    if (((man = idanpc) >= 0) && (idanpc < atom -> nlocal)) {
      atom -> type[man] = 2;
    }

    if (((mhi = idhipc) >= 0) && (idhipc < atom -> nlocal)) {
      atom -> type[mhi] = 2;
    }
  }

}

/*-----------------------------*/
/*Compute position of i-th bead*/
/*-----------------------------*/
std::array<double, 3> FixSMC::compute_xyz(long b){

    std::array<double,3> rtxyz, xyztemp;

    int count;
    int counts;

    double unwrap[3];

    count = 0;

    xyztemp[0] = 0;
    xyztemp[1] = 0;
    xyztemp[2] = 0;

    // Computing the distance between the proposed beads' ends using xyz values from different processors
    if ((atom->map(b) >= 0) && (atom->map(b) < (atom -> nlocal))) {
      domain -> unmap(atom -> x[atom->map(b)], atom -> image[atom->map(b)], unwrap);
      xyztemp[0] += unwrap[0];
      xyztemp[1] += unwrap[1];
      xyztemp[2] += unwrap[2];
      count += 1;
    }
    
    //MPI Barrier to avoid computational errors due to value collection
    MPI_Barrier(world);

    MPI_Allreduce(&xyztemp, &rtxyz, 3, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&count, &counts, 1, MPI_INT, MPI_SUM, world);

    for (int k = 0; k < 3; k++)
    {
      if (counts>0) rtxyz[k] /= counts;
      else rtxyz[k] = NAN;
    }

    return rtxyz;
}

/**
 * @brief Retrieves the MPI-synchronized, periodic-boundary-unwrapped position of a bead.
 * 
 * This function computes the Cartesian coordinates of a polymer bead (identified by its tag b)
 * across all processors in parallel. It handles the cases where:
 * 1. The bead is local to the current processor
 * 2. The bead is on a remote processor (ghost atom)
 * 3. Periodic boundary conditions have wrapped the bead's coordinates
 * 
 * The function uses MPI_Allreduce to gather position information from all processors and
 * computes an average position (which is correct since only one processor will have the bead
 * as a local atom).
 * 
 * **Algorithm Steps**:
 * 1. **Initialize arrays**:
 *    - xyztemp[3] = {0, 0, 0} (local contribution)
 *    - rtxyz[3] (return value, MPI result)
 *    - count = 0 (flag if bead is local)
 *    - counts (global count, should be 1)
 * 2. **Local check**: If bead tag b maps to local atom:
 *    - domain->unmap(): Convert wrapped coordinates to unwrapped using image flags
 *      * Accounts for periodic boundary crossings
 *      * Computes absolute position: x_abs = x_wrapped + image × box_length
 *    - Store unwrapped coordinates in xyztemp[3]
 *    - Set count = 1
 * 3. **MPI synchronization**: MPI_Barrier to ensure all processors ready
 * 4. **Global reduction**: 
 *    - MPI_Allreduce(xyztemp, rtxyz, 3, MPI_DOUBLE, MPI_SUM, world)
 *    - All processors send their contributions, one processor has valid data
 *    - MPI_Allreduce(count, counts, 1, MPI_INT, MPI_SUM, world)
 *    - Sums should equal 1 (bead found on exactly one processor)
 * 5. **Normalize**: Divide by counts to average (for robustness)
 *    - If counts > 0: rtxyz[k] /= counts
 *    - If counts == 0: rtxyz[k] = NAN (bead not found, error condition)
 * 6. **Return**: Array of three doubles (x, y, z coordinates)
 * 
 * **Periodic Boundary Handling**:
 * - domain->unmap() reconstructs absolute position accounting for box wrapping
 * - Essential for computing distances across periodic boundaries
 * - Allows correct distance calculations even when SMC spans box edge
 * 
 * **MPI Coordination**:
 * - Only one processor has bead as local atom (count=1)
 * - All other processors contribute zeros
 * - MPI_SUM reduces all zeros + one real position = the position
 * - Division by counts averages (only matters if counts=1, which it should)
 * 
 * **Special Considerations**:
 * - Expensive operation: barrier + 2 × Allreduce per bead
 * - Should cache results if same bead queried multiple times
 * - Used in post_integrate() to compute distances between SMC anchor and hinge
 * - Used in movement validation to check if proposed movement is within range
 * - Returns NAN if bead not found (should not happen in correct simulation)
 * - std::array<double,3> is C++11 feature for fixed-size array
 * 
 * @param b [long] Atom tag of the bead to locate
 * @return [std::array<double,3>] Unwrapped {x, y, z} coordinates, or {NAN, NAN, NAN} if not found
 * @see FixSMC::post_integrate(), domain->unmap()
 */
std::array<double, 3> FixSMC::compute_xyz(long b){

    std::array<double,3> rtxyz, xyztemp;

    int count;
    int counts;

    double unwrap[3];

    count = 0;

    xyztemp[0] = 0;
    xyztemp[1] = 0;
    xyztemp[2] = 0;

    // Computing the distance between the proposed beads' ends using xyz values from different processors
    if ((atom->map(b) >= 0) && (atom->map(b) < (atom -> nlocal))) {
      domain -> unmap(atom -> x[atom->map(b)], atom -> image[atom->map(b)], unwrap);
      xyztemp[0] += unwrap[0];
      xyztemp[1] += unwrap[1];
      xyztemp[2] += unwrap[2];
      count += 1;
    }
    
    //MPI Barrier to avoid computational errors due to value collection
    MPI_Barrier(world);

    MPI_Allreduce(&xyztemp, &rtxyz, 3, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&count, &counts, 1, MPI_INT, MPI_SUM, world);

    for (int k = 0; k < 3; k++)
    {
      if (counts>0) rtxyz[k] /= counts;
      else rtxyz[k] = NAN;
    }

    return rtxyz;
}

/**
 * @brief Logs SMC anchor and hinge atom types before movement/change (debugging utility).
 * 
 * This function writes a debug log entry recording the state of SMC i's anchor and hinge
 * atoms **before** any modification. It outputs:
 * - Timestep and processor rank
 * - SMC index i
 * - Anchor and hinge bead tags and atom types
 * 
 * The function is called before post_integrate() makes changes to SMC structure, allowing
 * post-processing analysis to track state transitions.
 * 
 * **Algorithm Steps**:
 * 1. **Check debug flag**: Only execute if debug == true
 * 2. **Map beads to local atoms**:
 *    - idhi = atom->map(map_to_beads(hing[i]))
 *    - idan = atom->map(map_to_beads(anch[i]))
 * 3. **Synchronize output across processors**:
 *    - For each processor nproc from 0 to nprocs-1:
 *      * Processor nproc writes if anchor is local (flag1 control)
 *      * MPI_Bcast(flag1, ..., nproc, world) - synchronize flags
 *      * All processors call MPI_Barrier
 *      * Processor nproc writes if hinge is local (flag2 control)
 *      * MPI_Bcast(flag2, ..., nproc, world) - synchronize flags
 *      * All processors call MPI_Barrier
 * 4. **Log format** (comma-separated CSV):
 *    - Anchor line: timestep,rank,smc_index,anchor_tag,anchor_type,
 *    - Hinge line (start): hinge_tag,hinge_type,
 *    - Post-movement continuation: ...see debug_post
 * 5. **File output**:
 *    - Opens "log_fix_smc.txt" in append mode
 *    - Each processor writes only its local data
 *    - Synchronization ensures no concurrent writes
 * 
 * **MPI Synchronization**:
 * - Strict processor-by-processor logging (nproc=0 then 1, then 2, etc.)
 * - Barriers prevent concurrent file I/O
 * - Broadcast flags ensure consistent state checking across MPI processes
 * 
 * **CSV Format**:
 * - Pre-change header (one line split across processes):
 *   * Timestep, processor rank, SMC index, anchor_tag, anchor_type, hinge_tag, hinge_type, [errors], [distance], [angle]
 * - Allows parsing with standard CSV tools
 * - One line per SMC movement per timestep (if debug==true)
 * 
 * **Overhead**:
 * - Multiple MPI barriers (2 × nprocs) per call
 * - File I/O per processor per call
 * - Significant slowdown: only use for debugging small systems
 * 
 * **Special Considerations**:
 * - Only called when debug == true (controlled by command-line argument)
 * - Called at beginning of SMC modification (before place/remove changes)
 * - Paired with debug_post() to track before/after states
 * - File grows unbounded (user must manage log file rotation)
 * - Synchronous MPI barriers required for consistent output order
 * 
 * @param i [int] SMC index to log
 * @see FixSMC::debug_post(), FixSMC::post_integrate()
 */
void FixSMC::debug_pre(int i){

  if (debug){
    int flag1 = 1;
    int flag2 = 1;
    
    long idhi = atom -> map(map_to_beads(hing[i]));
    long idan = atom -> map(map_to_beads(anch[i]));
    
    for (int nproc = 0; nproc < comm->nprocs; nproc++){
      if ((idan >= 0) && (idan < atom->nlocal) && (comm->me==nproc) && flag1) {
        debugfile.open("log_fix_smc.txt", std::ios_base::app);
        debugfile<<update -> ntimestep<<","<<comm->me<<","<<i<<","<<map_to_beads(anch[i])<<","<<atom->type[idan]<<",";
        debugfile.close();
        flag1 = 0;
      }
      MPI_Bcast(&flag1,1,MPI_INT,nproc,world);
      MPI_Barrier(world);

      if ((idhi >= 0) && (idhi < atom->nlocal) && (comm->me==nproc) && flag2) {
        debugfile.open("log_fix_smc.txt", std::ios_base::app);
        debugfile<<map_to_beads(hing[i])<<","<<atom->type[idhi]<<",";
        debugfile.close();
        flag2 = 0;
      }  

      MPI_Bcast(&flag2,1,MPI_INT,nproc,world);
      MPI_Barrier(world);
    }
  }

}

/**
 * @brief Logs SMC anchor and hinge atom types after movement/change plus error flags and metrics.
 * 
 * This function writes a debug log entry recording the **post-modification** state of SMC i's
 * anchor and hinge atoms. It is called **after** post_integrate() has modified the SMC structure,
 * allowing tracking of state transitions and error conditions.
 * 
 * **Algorithm Steps**:
 * 1. **Check debug flag**: Only execute if debug == true
 * 2. **Map beads to local atoms**:
 *    - idhi = atom->map(map_to_beads(hing[i]))
 *    - idan = atom->map(map_to_beads(anch[i]))
 * 3. **Synchronize output per processor** (same pattern as debug_pre):
 *    - For each processor nproc from 0 to nprocs-1:
 *      * Processor nproc logs anchor if local (flag1 control)
 *      * MPI_Bcast/Barrier synchronization
 *      * Processor nproc logs hinge + metrics if local (flag2 control)
 *      * MPI_Bcast/Barrier synchronization
 * 4. **Log format continuation** from debug_pre:
 *    - Anchor line: anchor_tag,anchor_type,
 *    - Hinge + metrics line: hinge_tag,hinge_type,err1,err2,err3,err4,err5,dist,tan\\n
 * 5. **Error flags** (boolean, written as 0 or 1):
 *    - err1: Probability rejection (random > kon/koff)
 *    - err2: Collision with another SMC
 *    - err3: Blocked bead type encountered
 *    - err4: Distance too large (SMC too stretched)
 *    - err5: Angle/tangent constraint violated
 * 6. **Metrics**:
 *    - dist: Distance between anchor and hinge (Euclidean)
 *    - tan: Tangent value (angle dot product or similar metric)
 *    - NAN values indicate distance could not be computed
 * 7. **File output**:
 *    - Continues append to "log_fix_smc.txt"
 *    - Ends line with std::endl (newline)
 *    - CSV format complete on newline
 * 
 * **CSV Line Format** (complete):
 * ```
 * timestep,rank,smc_index,anch_tag,anch_type_pre,hinge_tag,hinge_type_pre,
 * anch_tag,anch_type_post,hinge_tag,hinge_type_post,err1,err2,err3,err4,err5,dist,tan
 * ```
 * 
 * **Error Analysis in Post-Processing**:
 * - If all err flags = 0: Movement succeeded (types should differ pre/post)
 * - If err1 = 1: Movement rejected probabilistically (types unchanged)
 * - If err2-5 = 1: Movement rejected for physical reasons (types unchanged)
 * - If dist = NAN: Distance calculation failed (bead not found)
 * 
 * **MPI Synchronization**:
 * - Same strict processor-by-processor pattern as debug_pre
 * - 2 × nprocs barriers per call
 * - Ensures file I/O is serialized and consistent
 * 
 * **Overhead**:
 * - Multiple MPI barriers (2 × nprocs) per call
 * - File I/O and std::endl (flush after every entry)
 * - Significant slowdown: only use for small systems or occasional timesteps
 * 
 * **Special Considerations**:
 * - Only called when debug == true
 * - Paired with debug_pre() at start/end of SMC modifications
 * - Called even if movement fails (to log error reasons)
 * - std::endl flushes buffer (slow) but ensures data safety
 * - Log grows unbounded; user must manage rotation
 * - Error flags allow post-processing to identify rejection reasons
 * 
 * @param i [int] SMC index to log
 * @param err1 [bool] Probability rejection flag
 * @param err2 [bool] Collision with another SMC flag
 * @param err3 [bool] Blocked bead type flag
 * @param err4 [bool] Distance too large flag
 * @param err5 [bool] Angle/tangent constraint flag
 * @param dist [double] Distance between anchor and hinge (Euclidean)
 * @param tan [double] Tangent/angle metric value
 * @see FixSMC::debug_pre(), FixSMC::post_integrate()
 */
void FixSMC::debug_post(int i, bool err1, bool err2, bool err3, bool err4, bool err5, double dist, double tan){
  
  if (debug){
    int flag1 = 1;
    int flag2 = 1;

    long idhi = atom -> map(map_to_beads(hing[i]));
    long idan = atom -> map(map_to_beads(anch[i]));

    for (int nproc = 0; nproc < comm->nprocs; nproc++){
      if ((idan >= 0) && (idan < atom->nlocal) && (comm->me==nproc) && flag1) {
        debugfile.open("log_fix_smc.txt", std::ios_base::app);
        debugfile<<map_to_beads(anch[i])<<","<<atom->type[idan]<<",";
        debugfile.close();
        flag1 = 0;
      }
      MPI_Bcast(&flag1,1,MPI_INT,nproc,world);
      MPI_Barrier(world);

      if ((idhi >= 0) && (idhi < atom->nlocal) && (comm->me==nproc)  && flag2) {
        debugfile.open("log_fix_smc.txt", std::ios_base::app);
        debugfile<<map_to_beads(hing[i])<<","<<atom->type[idhi]<<","<<err1<<","<<err2<<","<<err3<<","<<err4<<","<<err5<<","<<dist<<","<<tan<<std::endl;
        debugfile.close();
        flag2 = 0;
      }

      MPI_Bcast(&flag2,1,MPI_INT,nproc,world);
      MPI_Barrier(world);
    }
  }
}

double FixSMC::memory_usage() {
  double bytes = 2 * smcnum * sizeof(long) + atom->natoms * sizeof(long);
  return bytes;
}

/*----------------------------------------------------*/

/**
 * @brief Serializes SMC state to restart file for checkpoint/restart capability.
 * 
 * **Restart File Format**:
 * 1. Size field: int (size in bytes of following data)
 *    - size = (3 + 2×smcnum) × sizeof(long)
 * 2. next_reneighbor: long (timestep of next neighbor list rebuild)
 *    - Critical for proper integration continuation
 * 3. ntimestep: long (simulation timestep when restart written)
 *    - Used for validation on restart
 * 4. smcnum: long (number of SMCs in simulation)
 *    - Checked against current smcnum on restart
 * 5. SMC data pairs (one per SMC):
 *    - anch[i]: long (anchor bead position)
 *    - hing[i]: long (hinge bead position)
 *    - Total: 2 × smcnum longs
 * 
 * @param fp [FILE*] Open file pointer for restart file (write mode)
 * @see FixSMC::restart(), LAMMPS restart mechanism
 */
void FixSMC::write_restart(FILE * fp) {

  int rn = 0;
  long rlist[3 + 2 * smcnum];

  rlist[rn++] = static_cast < long > (next_reneighbor);
  rlist[rn++] = static_cast < long > (update -> ntimestep);

  rlist[rn++] = static_cast < long > (smcnum);

  // Saving SMCs' ends positions
  for (int i = 0; i < smcnum; i++) {
    rlist[rn++] = anch[i];
    rlist[rn++] = hing[i];
  }

  if (comm -> me == 0) {
    int size = rn * sizeof(long);
    fwrite( & size, sizeof(int), 1, fp);
    fwrite(rlist, sizeof(long), rn, fp);
  }

}

/**
 * @brief Restores SMC state from restart file after restart/continuation.
 * 
 * **Algorithm Steps**:
 * 1. **Initialize buffer pointer**:
 *    - rn = 0 (index into buffer)
 *    - rlist = (long*)buf (cast buffer pointer)
 * 2. **Read next_reneighbor** (for integration continuity):
 *    - next_reneighbor = rlist[rn++]
 *    - Determines when neighbor list will be rebuilt
 *    - Critical for force/neighbor list consistency
 * 3. **Validate timestep consistency**:
 *    - ntimestep_restart = rlist[rn++]
 *    - Check: ntimestep_restart == update->ntimestep
 *    - Error if mismatch: "Must not reset timestep when restarting fix smc"
 *    - Prevents SMC data from being applied to wrong timestep
 * 4. **Validate SMC count**:
 *    - smcnum_rest = rlist[rn++]
 *    - Check: smcnum_rest == smcnum (from input command)
 *    - Error if mismatch: "Invalid restart, number of SMCs has changed!"
 *    - Ensures array sizes match (anch[] and hing[] allocation)
 * 5. **Restore SMC positions**:
 *    - For i=0 to smcnum-1:
 *      * anch[i] = rlist[rn++]
 *      * hing[i] = rlist[rn++]
 *    - Direct array assignment (no additional initialization)
 *    - Assumes anch[] and hing[] already allocated
 * 6. **Activate SMC atoms** (place_smc):
 *    - For i=0 to smcnum-1:
 *      * place_smc(anch[i], hing[i], false)
 *      * false: indicates moving SMC (not newly loaded)
 *      * Updates atom types, creates bonds, rebuilds topology
 * 7. **Synchronization**:
 *    - After loop: MPI_Barrier(world)
 *    - Ensures all processors have restored state
 * 
 * @param buf [char*] Buffer pointer to restart file data (from LAMMPS)
 * @see FixSMC::write_restart(), FixSMC::place_smc(), LAMMPS restart mechanism
 * @throws error->all if timestep mismatch with restart file
 * @throws error->all if SMC count mismatch with restart file
 */
void FixSMC::restart(char * buf) {

  int rn = 0;
  long * rlist = (long * ) buf;

  next_reneighbor = static_cast < bigint > (rlist[rn++]);

  bigint ntimestep_restart = static_cast < bigint > (rlist[rn++]);
  if (ntimestep_restart != update -> ntimestep)
    error -> all(FLERR, "Must not reset timestep when restarting fix smc");

  int smcnum_rest = rlist[rn++];
  if (smcnum_rest != smcnum)
    error -> all(FLERR, "Invalid restart, number of SMCs has changed!");

  // Loading SMCs' ends positions
  for (int j = 0; j < smcnum; j++) {
    anch[j] = static_cast < long > (rlist[rn++]);
    hing[j] = static_cast < long > (rlist[rn++]);
  }

}

/*---------------------*/
/*Return number of SMCs*/
/*---------------------*/
double FixSMC::compute_scalar() {
  return smcnum;
}

/*--------------------------------------------------------------------------------*/
/*Return left end position of i-th SMCs' left end if flag=0 or right end if flag=1*/
/*--------------------------------------------------------------------------------*/
double FixSMC::compute_array(int i, int flag) {
  int rflag;

  // Check which is the left end
  if (hdir != 0) {
    rflag = hdir / abs(hdir);
  } else {
    rflag = -adir / abs(adir);
  }

  rflag = (1 + rflag / (abs(rflag))) / 2;

  if (flag) {
    if (rflag) return hing[i];
    else return anch[i];
  } else {
    if (rflag) return anch[i];
    else return hing[i];
  }
}

/* ----------------------------------------------------------------------
   update special neighbors topology for atoms with id1 or id2
   rebuild special list for each atom that is influenced by id1 or id2
    influenced atoms are those that have id1 or id2 as a neighbor
    or have id1 or id2 as a neighbor of a neighbor
    rebuild_special_one() is called for each influenced atom
---------------------------------------------------------------------- */

void FixSMC::update_topology(int id1, int id2)
{
  int i,j,k,n,influence,influenced,found;
  tagint *slist;

  tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    influenced = 0;
    slist = special[i];

    influence = 0;
    if (tag[i] == id1 || tag[i] == id2) influence = 1;
    else {
      n = nspecial[i][2];
      found = 0;
      for (k = 0; k < n; k++)
        if (slist[k] == id1 || slist[k] == id2) found++;
      if (found == 2) influence = 1;
    }
    if (!influence) continue;
    influenced = 1;
  
    if (influenced) rebuild_special_one(i);
  
  }

}

void FixSMC::rebuild_special_one(int m)
{
  int i,j,n,n1,cn1,cn2,cn3;
  tagint *slist;

  tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  // existing 1-2 neighs of atom M

  slist = special[m];
  n1 = nspecial[m][0];
  cn1 = 0;
  for (i = 0; i < n1; i++)
    copy[cn1++] = slist[i];

  // new 1-3 neighs of atom M, based on 1-2 neighs of 1-2 neighs
  // exclude self
  // remove duplicates after adding all possible 1-3 neighs

  cn2 = cn1;
  for (i = 0; i < cn1; i++) {
    n = atom->map(copy[i]);
    if (n < 0)
      error->one(FLERR,"Fix smc needs ghost atoms from further away");
    slist = special[n];
    n1 = nspecial[n][0];
    for (j = 0; j < n1; j++)
      if (slist[j] != tag[m]) copy[cn2++] = slist[j];
  }

  cn2 = dedup(cn1,cn2,copy);
  if (cn2 > atom->maxspecial)
    error->one(FLERR,"Special list size exceeded in fix bond/create");

  // new 1-4 neighs of atom M, based on 1-2 neighs of 1-3 neighs
  // exclude self
  // remove duplicates after adding all possible 1-4 neighs

  cn3 = cn2;
  for (i = cn1; i < cn2; i++) {
    n = atom->map(copy[i]);
    if (n < 0)
      error->one(FLERR,"Fix smc needs ghost atoms from further away");
    slist = special[n];
    n1 = nspecial[n][0];
    for (j = 0; j < n1; j++)
      if (slist[j] != tag[m]) copy[cn3++] = slist[j];
  }

  cn3 = dedup(cn2,cn3,copy);
  if (cn3 > atom->maxspecial)
    error->one(FLERR,"Special list size exceeded in fix smc");

  // store new special list with atom M

  nspecial[m][0] = cn1;
  nspecial[m][1] = cn2;
  nspecial[m][2] = cn3;
  memcpy(special[m],copy,cn3*sizeof(int));
}

/* ----------------------------------------------------------------------
   remove all ID duplicates in copy from Nstart:Nstop-1
   compare to all previous values in copy
   return N decremented by any discarded duplicates
------------------------------------------------------------------------- */

int FixSMC::dedup(int nstart, int nstop, tagint *copy)
{
  int i;

  int m = nstart;
  while (m < nstop) {
    for (i = 0; i < m; i++)
      if (copy[i] == copy[m]) {
        copy[m] = copy[nstop-1];
        nstop--;
        break;
      }
    if (i == m) m++;
  }

  return nstop;
}