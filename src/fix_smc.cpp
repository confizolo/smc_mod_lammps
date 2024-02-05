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

FixSMC::FixSMC(LAMMPS * lmp, int narg, char ** arg):
  Fix(lmp, narg, arg),
  anch(nullptr), hing(nullptr), smctype(0), smcbtype(0), smcnum(0), debug(0) {
    if (lmp -> citeme) lmp -> citeme -> add(cite_fix_smc);
    // Number of arguments for the fix. The first three arguments are parsed by Fix base class constructor.
    // The rest are specific to this fix. 15 are mandatory
    // 4. nevery: Attempt the jump every nevery iteration
    // 5. seed: seed for random number draw
    // 6. prob: probability to accept proposed movement
    // 7. lpol: length of polymer(s) in solution
    // 8. solsize: number of atoms belonging to the solution of polymers
    // 9. poltype: select shape of polymer(s) (either "linear" or "ring")
    // 10. adir: attempted movement of anchor (next attempted atom id: current anchor + adir)
    // 11. hdir: attempted movement of hinge (next attempted atom id: current hinge + hdir)
    // 12. smcnum: number of deployed SMCs
    // 13. smctype: atom type of beads representing SMCs' ends
    // 14. smcbtype: SMCs' bond type after the first deployment
    // 15. smcbitype: SMCs' bond type at the first deployment
    // 16. cutoff: distance cutoff for attempted movements (these are accepted only if distance between new anchor and hinge is below the cutoff)
    // 17. tancoff: tangent cutoff for attempted movements (these are accepted only if scalar product of the tangents is smaller than this value)
    // 18. initmode: define initialisation mode of extruders:
    //    1.  "random": deploys randomly the SMCs
    //    2.  "distributed": assign at least one SMC per polymer and then distribute remaining randomly
    //    3.  "full-distributed": distributes evenly SMCs over the polymers
    //    4.  "fixed": assign SMCs' initial positions according to the filename defined at 23.
    // 19. kon: probability to load a free extruder every nevery step
    // 20. koff: probability to unload an extruder every nevery step
    // 21. patches: number of patches per bead to recolor (type associated is smctype + 1)
    // 22. debug: activate debug mode with detailed report on log_fix_smc.txt
    // 23. fixFname: file containing hinge and anchor position list separated by a space (optional)
    // 24. blockbeads: type of beads that the extruder cannot grab, can be listed as an arbitrary long list (e.g.: 2 3 4 ...) (optional)

    // Check on the number of arguments given to the fix
    if (narg < 19) error -> all(FLERR, "Illegal fix smc command");

    // Define attempt rate
    nevery = utils::inumeric(FLERR, arg[3], false, lmp);
    if (nevery <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Flag to activate dump in restart file
    restart_global = 1;

    // Activate flag for array and scalar return
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

    // Define movement of SMCs' anchor
    maxadir = utils::inumeric(FLERR, arg[9], false, lmp);

    // Define maximum movement of SMCs' hinge
    maxhdir = utils::inumeric(FLERR, arg[10], false, lmp);
    
    if (maxhdir * maxadir > 0) error -> all(FLERR, "Illegal fix smc command, hinge and must not have same direction");

    // Define number of SMCs to deploy
    smcnum = utils::inumeric(FLERR, arg[11], false, lmp);
    if (smcnum <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Define atom type of SMCs' ends
    smctype = utils::inumeric(FLERR, arg[12], false, lmp);
    if (smctype <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Define deployment bond type of SMCs
    smcbtype = utils::inumeric(FLERR, arg[13], false, lmp);
    if (smcbtype <= 0) error -> all(FLERR, "Illegal fix smc command");

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

    // To get a different random number every time the program is executed
    srand(time(NULL) * seed);
  }

/* ---------------------------------------------------------------------- */

FixSMC::~FixSMC() {

  // // Loop over the instantiated SMCs, remove bonds and change types
  // for (int i = 0; i < smcnum; i++) {
  //   remove_smc(anch[i], hing[i]);
  // }

  // Deleting pointers
  delete random_equal;
  delete anch;
  delete hing;
  delete blockt;
  
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

/* -------------------------------------------*/
/*Main extrusion code to run after integration*/
/*--------------------------------------------*/
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
    if ((maxhdir == 0) && (maxadir == 0)) return;

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
        MPI_Barrier(world);
        debug_post(i,1,0,0,0,0,0,0);
        continue;
      }

      if (maxadir>1){
        if (comm -> me == 0) jrand = random_equal -> uniform();
        MPI_Bcast( & jrand, 1, MPI_DOUBLE, 0, world);
      }

      if (maxhdir>1){
        if (comm -> me == 0) krand = random_equal -> uniform();
        MPI_Bcast( & krand, 1, MPI_DOUBLE, 0, world);
      }

      if (maxadir == 0){
        adir = 0;
      }
      else {
        adir = (int)(maxadir * jrand) + (maxadir) / (abs(maxadir));
      }

      if (maxhdir == 0){
        hdir = 0;
      }
      else if (maxhdir > 0) {
        hdir = (int)(maxhdir * krand) + (maxhdir) / (abs(maxhdir));
      }

      // Temporary direction if the smc is going towards the polymer end or another smc bead
      tempadir = adir;
      temphdir = hdir;

      // Check if the SMCs' ends are moving over the polymer ends
      if (hdir / (abs(hdir)) < 0) {
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
      if (adir / (abs(adir)) < 0) {
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
        MPI_Barrier(world);
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
        MPI_Barrier(world);
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
        MPI_Barrier(world);
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
        MPI_Barrier(world);
        debug_post(i,0,0,0,1,0,dist,tan);
        continue;
      }

      // Check if the product between the two tangents is smaller than a cutoff
      if (tan < tancoff) {
        // Debug print if got rejected because of tangent product value
        debug_pre(i);
        MPI_Barrier(world);
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
/*Map polymer bead index to real simulation bead index*/
/*----------------------------------------------------*/
long FixSMC::map_to_beads(long i){
  return (i-1)*(npatches+1) + 1;
}

/*----------------------------------------------------*/
/*Check if it is possible to place an anchor on bead i*/
/*----------------------------------------------------*/
bool FixSMC::check_avl(long i){
  long tmphing;
  long mdbead;
  bool flag = 0;

  int mnew;
  int mannew;

  int idnewhi;
  int idnewan;

  // Define hinge position according to movement direction
  if (maxhdir != 0) tmphing = i + 2 * maxhdir / (abs(maxhdir));
  else tmphing = i - 2 * (maxadir) / (abs(maxadir));

  idnewhi = atom -> map(map_to_beads(tmphing));
  idnewan = atom -> map(map_to_beads(i));
  
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
      if ((((mnew = idnewhi) >= 0) && (idnewhi < atom -> nlocal)) && atom -> type[mnew] == blockt[j]) {
        flag=1; 
        MPI_Bcast(&flag,1,MPI_INT,comm->me,world); 
        break;
      }
    }

  if (flag) return 0;

  return 1;
}

/*------------------------------------------------------*/
/*Compile a list of available position for SMCs' anchors*/
/*------------------------------------------------------*/
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

/*-------------*/
/*Load i-th SMC*/
/*-------------*/
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

  if (maxhdir != 0) hing[i] = anch[i] + 2;
  else hing[i] = anch[i] - 2 * maxadir / abs(maxadir);

  MPI_Barrier(world);

  // Cast the chosen position to each processor
  MPI_Bcast(anch, smcnum, MPI_LONG, 0, world);
  MPI_Bcast(hing, smcnum, MPI_LONG, 0, world);

  MPI_Barrier(world);
}

/*--------------*/
/*Place SMC anchor and hinge*/
/*--------------*/
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
    atom -> type[mhi] = smctype;

    bool create = 1;
    for (int ibond = 0; ibond < atom -> num_bond[mhi]; ibond++) {
      if ((atom -> bond_type[mhi][ibond] == smcbtype) || (atom -> bond_type[mhi][ibond] == smcbitype)) {
        create = 0;
      }
    }

    if (create) {
      // Creating new SMC bond
      if (atom -> num_bond[mhi] == atom -> bond_per_atom) error -> one(FLERR, "New bond exceeded bonds per atom limit of {} in create_bonds", atom -> bond_per_atom);
      if (newsmc){
        atom -> bond_type[mhi][atom -> num_bond[mhi]] = smcbitype;
      }
      else {
        atom -> bond_type[mhi][atom -> num_bond[mhi]] = smcbtype;
      }
      atom -> bond_atom[mhi][atom -> num_bond[mhi]] = map_to_beads(a);
      atom -> num_bond[mhi]++;
    }
  }

  // Recoloring patches
  for (int c = 1; c <= npatches; c++)
  {
    long idhipc = atom -> map(map_to_beads(h) + c);
    long idanpc = atom -> map(map_to_beads(a) + c);

    if (((man = idanpc) >= 0) && (idanpc < atom -> nlocal)) {
      atom -> type[man] = smctype + 1;
    }

    if (((mhi = idhipc) >= 0) && (idhipc < atom -> nlocal)) {
      atom -> type[mhi] = smctype + 1;
    }
  }

}

/*--------------*/
/*Remove SMC anchor and hinge*/
/*--------------*/
void FixSMC::remove_smc(long a, long h) {
  long mhi;
  long man;

  long idhi = atom -> map(map_to_beads(h));
  long idan = atom -> map(map_to_beads(a));

  // Get bond histories to apply bond changes
  auto histories = modify -> get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  // Changing type of the old anchor
  if (((man = idan) >= 0) && (idan < atom -> nlocal)) {
    atom -> type[man] = 1;
  }

  if (((mhi = idhi) >= 0) && (idhi < atom -> nlocal)) {
    // Resetting type
    atom -> type[mhi] = 1;

    // Deleting old SMC bond
    for (int ibond = 0; ibond < atom -> num_bond[mhi]; ibond++) {
      if ((atom -> bond_type[mhi][ibond] == smcbtype) || (atom -> bond_type[mhi][ibond] == smcbitype)) {
        atom -> bond_type[mhi][ibond] = atom -> bond_type[mhi][atom -> num_bond[mhi] - 1];
        atom -> bond_atom[mhi][ibond] = atom -> bond_atom[mhi][atom -> num_bond[mhi] - 1];

        if (n_histories > 0)
          for (auto & ihistory: histories) {
            dynamic_cast < FixBondHistory * > (ihistory) -> shift_history(mhi, ibond, atom -> num_bond[mhi] - 1);
            dynamic_cast < FixBondHistory * > (ihistory) -> delete_history(mhi, atom -> num_bond[mhi] - 1);
          }

        atom -> num_bond[mhi]--;
        break;
      }
    }
  }

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

/*-----------------------*/
/*Pre change debug print*/
/*-----------------------*/
void FixSMC::debug_pre(int i){
  int flag1 = 1;
  int flag2 = 1;
  
  long idhi = atom -> map(map_to_beads(hing[i]));
  long idan = atom -> map(map_to_beads(anch[i]));
  
  for (int nproc = 0; nproc < comm->nprocs; nproc++){
    if ((debug) && (idan >= 0) && (idan < atom->nlocal) && (comm->me==nproc) && flag1) {
      debugfile.open("log_fix_smc.txt", std::ios_base::app);
      debugfile<<update -> ntimestep<<","<<comm->me<<","<<i<<","<<map_to_beads(anch[i])<<","<<atom->type[idan]<<",";
      debugfile.close();
      flag1 = 0;
    }
    MPI_Bcast(&flag1,1,MPI_INT,nproc,world);
    MPI_Barrier(world);

    if ((debug) && (idhi >= 0) && (idhi < atom->nlocal) && (comm->me==nproc) && flag2) {
      debugfile.open("log_fix_smc.txt", std::ios_base::app);
      debugfile<<map_to_beads(hing[i])<<","<<atom->type[idhi]<<",";
      debugfile.close();
      flag2 = 0;
    }  

    MPI_Bcast(&flag2,1,MPI_INT,nproc,world);
    MPI_Barrier(world);
  }
}

/*-----------------------*/
/*Post change debug print*/
/*-----------------------*/
void FixSMC::debug_post(int i, bool err1, bool err2, bool err3, bool err4, bool err5, double dist, double tan){
  int flag1 = 1;
  int flag2 = 1;

  long idhi = atom -> map(map_to_beads(hing[i]));
  long idan = atom -> map(map_to_beads(anch[i]));

  for (int nproc = 0; nproc < comm->nprocs; nproc++){
    if ((debug) && (idan >= 0) && (idan < atom->nlocal) && (comm->me==nproc) && flag1) {
      debugfile.open("log_fix_smc.txt", std::ios_base::app);
      debugfile<<map_to_beads(anch[i])<<","<<atom->type[idan]<<",";
      debugfile.close();
      flag1 = 0;
    }
    MPI_Bcast(&flag1,1,MPI_INT,nproc,world);
    MPI_Barrier(world);

    if ((debug) && (idhi >= 0) && (idhi < atom->nlocal) && (comm->me==nproc)  && flag2) {
      debugfile.open("log_fix_smc.txt", std::ios_base::app);
      debugfile<<map_to_beads(hing[i])<<","<<atom->type[idhi]<<","<<err1<<","<<err2<<","<<err3<<","<<err4<<","<<err5<<","<<dist<<","<<tan<<std::endl;
      debugfile.close();
      flag2 = 0;
    }

    MPI_Bcast(&flag2,1,MPI_INT,nproc,world);
    MPI_Barrier(world);
  }
}
/*------------------------------------*/
/*Memory usage of hing and anch arrays*/
/*------------------------------------*/
double FixSMC::memory_usage() {
  double bytes = 2 * smcnum * sizeof(long) + atom->natoms * sizeof(long);
  return bytes;
}

/*-------------------------------------------*/
/*Add restart information in the restart file*/
/*-------------------------------------------*/
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

/*---------------------------------------------------*/
/*Use state info from restart file to restart the Fix*/
/*---------------------------------------------------*/
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
