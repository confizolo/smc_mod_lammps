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

#include "fix_smc.h"

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
    // 8. poltype: select shape of polymer(s) (either "linear" or "ring")
    // 9. adir: attempted movement of anchor (next attempted atom id: current anchor + adir)
    // 10. hdir: attempted movement of hinge (next attempted atom id: current hinge + hdir)
    // 11. smcnum: number of deployed SMCs
    // 12. smctype: atom type of beads representing SMCs' ends
    // 13. smcbtype: SMCs' bond type after the first deployment
    // 14. smcbitype: SMCs' bond type at the first deployment
    // 15. cutoff: distance cutoff for attempted movements (these are is accepted only if distance between new anchor and hinge is below the cutoff)
    // 16. initmode: define initialisation mode of extruders:
    //    1.  "random": deploys randomly the SMCs
    //    2.  "distributed": assign at least one SMC per polymer and then distribute remaining randomly
    //    3.  "full-distributed": distributes evenly SMCs over the polymers
    // 17. kon: probability to load a free extruder every nevery step
    // 18. koff: probability to unload an extruder every nevery step
    // 18.5 anchor position
    // 19. blockbeads: type of beads that the extruder cannot grab, can be listed as an arbitrary long list (e.g.: 2 3 4 ...)

    // Check on the number of arguments given to the fix
    if (narg < 18) error -> all(FLERR, "Illegal fix smc command");

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

    // Define lenght of polymer(s) present in simulation to check boundary conditions
    lpol = utils::inumeric(FLERR, arg[6], false, lmp);
    if (lpol <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Define type of extrusion according to the kind of polymer(s) present in simulation
    if (strcmp(arg[7], "linear") == 0) {
      ring = 0;
    } else if (strcmp(arg[7], "ring") == 0) {
      ring = 1;
    } else {
      error -> all(FLERR, "Illegal fix smc command, indefinite polymer type");
    }

    // Define movement of SMCs' anchor
    maxadir = utils::inumeric(FLERR, arg[8], false, lmp);

    // Define movement of SMCs' hinge
    maxhdir = utils::inumeric(FLERR, arg[9], false, lmp);
    //if (hdir * adir >= 0) error -> all(FLERR, "Illegal fix smc command, hinge and must not have same direction");

    // Define number of SMCs to deploy
    smcnum = utils::inumeric(FLERR, arg[10], false, lmp);
    if (smcnum <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Define atom type of SMCs' ends
    smctype = utils::inumeric(FLERR, arg[11], false, lmp);
    if (smctype <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Define deployment bond type of SMCs
    smcbtype = utils::inumeric(FLERR, arg[12], false, lmp);
    if (smcbtype <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Define deployment bond type of SMCs after the first movement
    smcbitype = utils::inumeric(FLERR, arg[13], false, lmp);
    if (smcbitype <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Define distance cutoff to be checked before SMCs' movement
    cutoff = utils::numeric(FLERR, arg[14], false, lmp);
    if (cutoff < 0)
      error -> all(FLERR, "Illegal fix smc command");

    initmode = 0;

    // Define init mode for SMCs' deployment 
    if (strcmp(arg[15], "random") == 0) {
      initmode = 0;
    } else if (strcmp(arg[15], "distributed") == 0) {
      initmode = 1;
    } else if (strcmp(arg[15], "full-distributed") == 0) {
      initmode = 2;
    } else if (strcmp(arg[15], "fixed") == 0) {
      initmode = 3;
    } else {
      error -> all(FLERR, "Illegal fix smc command, initmode not present");
    }

    // Define loading probability of one SMC on simulation polymers
    kon = utils::numeric(FLERR, arg[16], false, lmp);
    if ((kon <= 0) || (kon>1))
      error -> all(FLERR, "Illegal fix smc command, kon is out of the interval (0,1]");

    // Define unloading probability of one SMC on simulation polymers
    koff = utils::numeric(FLERR, arg[17], false, lmp);
    if ((koff < 0) || (koff>1))
      error -> all(FLERR, "Illegal fix smc command, koff is out of the interval [0,1]");

    // fixed_hinge_position = 200;
    fixed_hinge_position = utils::inumeric(FLERR, arg[18], false, lmp);

    nblockt = narg - 19;

    blockt = new int[nblockt];

    for (int i = 19; i < narg; i++) {
      blockt[i-19] = utils::inumeric(FLERR, arg[i], false, lmp);;
    }
    if ((debug)) utils::logmesg(lmp, "atom type --> {}, nblockt = {}", blockt[0], nblockt);


    xyzanch = nullptr;
    xyzhing = nullptr;

    anch = new long[smcnum];
    hing = new long[smcnum];

    // Initialise SMC positions
    for (int i = 0; i < smcnum; i++) {
      hing[i] = -1;
      anch[i] = -1;
    }

    // Initialise list of available positions
    av_list = new long[atom->natoms];
    num_avl = atom->natoms / 4;

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

  memory -> destroy(xyzanch);
  memory -> destroy(xyzhing);
  
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

    for (int i = 0; i < smcnum; i++)
    { 
      if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Anchors are" + std::to_string(anch[i]) + " " + std::to_string(hing[i]) + "\n");
    }

    for (int i = 0; i < smcnum; i++) {
      // Place the loaded SMCs
      if ((anch[i]<0) || (hing[i]<0)) continue;
      place_smc(anch[i],hing[i], true);
    }

    return;

  } 

  // Propose movements with a fixed frequency
  else if (update -> ntimestep % nevery == 0){

    // If movement is disabled along both directions stop execution
    if ((hdir == 0) && (adir == 0)) return;

    double * xyzanchtemp = nullptr;
    double * xyzhingtemp = nullptr;

    int * anchcount = nullptr;
    int * hingcount = nullptr;

    int * hingcounts = nullptr;
    int * anchcounts = nullptr;

    double unwrap[3];

    int mnew;
    int mannew;

    double dist = 0;

    auto histories = modify -> get_fix_by_style("BOND_HISTORY");
    int n_histories = histories.size();

    int idnewhi;
    int idhi;
    int idnewan;
    int idan;

    double jrand;
    double krand;

    bool flag = 0;
    bool l_flag = 0;

    int tempadir = 0;
    int temphdir = 0;

    double rand;

    for (int i = 0; i < smcnum; i++) {

      // if ((debug)) utils::logmesg(lmp, "hi\n");
      for (int l_sample = 0; l_sample < std::max(abs(maxadir), abs(maxhdir)); l_sample ++){


      if (comm -> me == 0) jrand = random_equal -> uniform();
      MPI_Bcast( & jrand, 1, MPI_DOUBLE, 0, world);

      if (comm -> me == 0) krand = random_equal -> uniform();
      MPI_Bcast( & krand, 1, MPI_DOUBLE, 0, world);

      // fudge by + 1 is required to not have a divide by zero later on...

      if (maxadir == 0){
        adir = 0;
      }
      else if(maxadir > 0){
        adir = round(maxadir * jrand) + 1;
      }
      else if (maxadir < 0){
        adir = round(maxadir * jrand) - 1;
      }

      if (maxhdir == 0){
        hdir = 0;
      }
      else if (maxhdir > 0) {
        hdir = round(maxhdir * krand) + 1;
      }
      else if (maxhdir < 0) {
        hdir = round(maxhdir * krand) - 1;
      }

      tempadir = adir;
      temphdir = hdir;

      if ((debug)) utils::logmesg(lmp, "lsample {} | adir: {}, hdir: {}\n",l_sample, adir, hdir);

      // Draw two random numbers for the unloading/loading
      double lrand;

      if (comm -> me == 0) lrand = random_equal -> uniform();

      MPI_Bcast( & lrand, 1, MPI_DOUBLE, 0, world);

      // Load/Unload SMC according to the kon/koff probability
      if ((anch[i]<0) || (hing[i]<0)){
        if (lrand < kon) {
          load_smc(i);
          place_smc(anch[i],hing[i],true);
        }
      }
      else{
        if (lrand < koff) {
          remove_smc(anch[i],hing[i]);
          anch[i]=-1;
          hing[i]=-1;
        }
      }

      MPI_Barrier(world);

      // Check if SMC are loaded
      // change to `continue` if removing the l_sample loop
      if ((anch[i]<0) || (hing[i]<0)) break;

      // Draw a random number for the jump attempt
      if (comm -> me == 0) rand = random_equal -> uniform();
      MPI_Bcast( & rand, 1, MPI_DOUBLE, 0, world);

      // First accept the move with a certain probability prob
      // change to `continue` if removing the l_sample loop
      if (rand > prob) {
        break;
      }

      // Temporary direction if the smc is going towards the polymer end or another smc bead

      // Check if the SMCs' ends are moving over the polymer ends
      if (hdir / abs(hdir) < 0) {
        if ((hing[i] + hdir) % lpol == 0){
          if (!ring) temphdir = 0;
          else temphdir = (lpol-1);
        }
      } else {
        if ((hing[i] + hdir) % lpol == 1){
          if (!ring) temphdir = 0;
          else temphdir = 1-lpol;
        }
      }
      if (adir / abs(adir) < 0) {
        if ((anch[i] + adir) % lpol == 0){
          if (!ring) tempadir = 0;
          else tempadir = (lpol-1);
        }
      } else {
        if ((anch[i] + adir) % lpol == 1){
          if (!ring) tempadir = 0;
          else tempadir = 1-lpol;
        }
      }

      // Skip movement if both the ends are not moving
      if ((tempadir == 0) && (temphdir == 0)) continue;

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

      // change to `continue` if removing the l_sample loop
      if (flag) {
        break;
      }
      
      if ((debug)) utils::logmesg(lmp, "SMC {} Current anchor {}, current hinge {}; next anchor {}, next hinge {} \n", i,anch[i],hing[i],anch[i]+tempadir,hing[i]+temphdir);

      // Finding local identifier of the atom
      idan = atom -> map(anch[i]);
      idhi = atom -> map(hing[i]);

      for (int j = 0; j < nblockt ; j++)
      {
          if ((((idnewan) >= 0) && (idnewan < atom -> nlocal)) && atom -> type[idnewan] == blockt[j]) {tempadir = 0; MPI_Bcast(&tempadir,1,MPI_INT,comm->me,world);}
          if ((((idnewhi) >= 0) && (idnewhi < atom -> nlocal)) && atom -> type[idnewhi] == blockt[j]) {temphdir = 0; MPI_Bcast(&temphdir,1,MPI_INT,comm->me,world);}
      }

      MPI_Barrier(world);
      if ((tempadir == 0) && (temphdir == 0)) continue;

      memory -> destroy(xyzanchtemp);
      memory -> create(xyzanchtemp, 3, "FixSMC::post_integrate()");
      memory -> destroy(xyzhingtemp);
      memory -> create(xyzhingtemp, 3, "FixSMC::post_integrate()");
      memory -> destroy(anchcount);
      memory -> create(anchcount, 1, "FixSMC::post_integrate()");
      memory -> destroy(hingcount);
      memory -> create(hingcount, 1, "FixSMC::post_integrate()");

      anchcount[0] = 0;
      hingcount[0] = 0;

      xyzanchtemp[0] = 0;
      xyzanchtemp[1] = 0;
      xyzanchtemp[2] = 0;

      xyzhingtemp[0] = 0;
      xyzhingtemp[1] = 0;
      xyzhingtemp[2] = 0;

      memory -> destroy(xyzanch);
      memory -> create(xyzanch, 3, "FixSMC::post_integrate()");
      memory -> destroy(xyzhing);
      memory -> create(xyzhing, 3, "FixSMC::post_integrate()");

      // Computing the distance between the proposed beads' ends using xyz values from different processors
      if (((mannew = idnewan) >= 0) && (idnewan < (atom -> nlocal))) {
        domain -> unmap(atom -> x[mannew], atom -> image[mannew], unwrap);
        xyzanchtemp[0] += unwrap[0];
        xyzanchtemp[1] += unwrap[1];
        xyzanchtemp[2] += unwrap[2];
        anchcount[0] += 1;
      }

      if (((mnew = idnewhi) >= 0) && (idnewhi < (atom -> nlocal))) {
        domain -> unmap(atom -> x[mnew], atom -> image[mnew], unwrap);

        xyzhingtemp[0] += unwrap[0];
        xyzhingtemp[1] += unwrap[1];
        xyzhingtemp[2] += unwrap[2];
        hingcount[0] += 1;
      }

      memory -> destroy(anchcounts);
      memory -> create(anchcounts, 1, "FixSMC::post_integrate()");
      memory -> destroy(hingcounts);
      memory -> create(hingcounts, 1, "FixSMC::post_integrate()");

      //MPI Barrier to avoid computational errors due to value collection
      MPI_Barrier(world);

      MPI_Allreduce(xyzanchtemp, xyzanch, 3, MPI_DOUBLE, MPI_SUM, world);
      MPI_Allreduce(xyzhingtemp, xyzhing, 3, MPI_DOUBLE, MPI_SUM, world);
      MPI_Allreduce(anchcount, anchcounts, 1, MPI_INT, MPI_SUM, world);
      MPI_Allreduce(hingcount, hingcounts, 1, MPI_INT, MPI_SUM, world);

      dist = 0;

      for (int k = 0; k < 3; k++) {
        xyzanch[k] = xyzanch[k] / anchcounts[0];
        xyzhing[k] = xyzhing[k] / hingcounts[0];
        dist += (xyzanch[k] - xyzhing[k]) * (xyzanch[k] - xyzhing[k]);
      }

      if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Number of counts is " + std::to_string(anchcounts[0]) + " Anchor " + std::to_string(hingcounts[0]) + " Hinge " + "\n");
      if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Proposed distance is " + std::to_string(sqrt(dist)) + "\n");

      // Check if the distance is small enough to accept the movement
      if (!(dist > cutoff * cutoff || (anchcounts[0] == 0) || (hingcounts[0] == 0))) {
        
        if ((temphdir!=0) || (tempadir!=0)){
          l_flag = 1;
          break;
        }
      }
      // Barrier to check that each processor has defined correctly each smc
      MPI_Barrier(world);

    } 
    // loop over samples
    if (l_flag){
      remove_smc(anch[i], hing[i]);
      place_smc(anch[i] + tempadir, hing[i] + temphdir, false);
      anch[i] += tempadir;
      hing[i] += temphdir;

    } 

    // Barrier to check that each processor has defined correctly each smc
    MPI_Barrier(world);

  } // loop over smc


    memory -> destroy(xyzanch);
    memory -> destroy(xyzhing);
    memory -> destroy(xyzanchtemp);
    memory -> destroy(xyzhingtemp);
    memory -> destroy(anchcount);
    memory -> destroy(hingcount);
    memory -> destroy(hingcounts);
    memory -> destroy(anchcounts);

    return;

  }

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
  if (hdir != 0) tmphing = i + 2 * hdir / abs(hdir);
  else tmphing = i - 2 * adir / abs(adir);

  idnewhi = atom -> map(tmphing);
  idnewan = atom -> map(tmphing);
  
  // Define SMC's center
  mdbead = (i + tmphing)/2;

  // Check if the smc is wrongly positioned (border conditions)
  if ((mdbead%lpol == 1) || (mdbead == 1) || (mdbead%lpol == 0)) {
    return 0;
  }

  flag = 0;
  // Check if we are superimposing other SMCs' beads
  for (int j = 0; j < smcnum; j++) {
    if ((anch[j]<0) || (hing[j]<0)) continue;
    if (((i == anch[j]) || (tmphing == hing[j])) || ((i == hing[j]) || (tmphing == anch[j])) || ((i + 1) == anch[j]) || ((i - 1) == anch[j]) || ((tmphing + 1) == anch[j]) || ((tmphing - 1) == anch[j]) || ((i + 1) == hing[j]) || ((i - 1) == hing[j]) || ((tmphing + 1) == hing[j]) || ((tmphing - 1) == hing[j])) {
      flag = 1;
      break;
    }
  }

  if (flag) return 0;

  for (int j = 0; j < nblockt ; j++)
    {
        if ((((mannew = idnewan) >= 0) && (idnewan < atom -> nlocal)) && atom -> type[mannew] == blockt[j]) {
          flag=1; 
          MPI_Bcast(&flag,1,MPI_INT,comm->me,world); 
          break;
        }
        if ((((mnew = idnewhi) >= 0) && (idnewhi < atom -> nlocal)) && atom -> type[mnew] == blockt[j]) {flag=1; MPI_Bcast(&flag,1,MPI_INT,comm->me,world); break;}
    }

  MPI_Barrier(world);

  if (flag) return 0;

  return 1;
}

/*------------------------------------------------------*/
/*Compile a list of available position for SMCs' anchors*/
/*------------------------------------------------------*/
void FixSMC::compile_avl_list(){
  num_avl = 0;

  // Scatter availability checking among the processors
  int atoms_per_rank = atom->natoms / comm->nprocs;

  // Getting the range of atoms to consider
  long min = comm->me * atoms_per_rank + 1;
  long max = comm->me * atoms_per_rank + atoms_per_rank;

  if (comm->me == comm->nprocs-1) max = atom->natoms;
  
  // Getting size of array to compute
  int size = max - min + 1;
  long temp_num_avl = 0;
  long * temp_avl_list = new long[size];

  // Filling the list with bead number atom->natoms + 1
  std::fill(temp_avl_list, temp_avl_list + size, atom->natoms + 1);

  // Define displacements for MPI parallelisation
  int displs[comm->nprocs];
  int rcounts[comm->nprocs];

  for (int i = 0; i < comm->nprocs; i++)
  {
    displs[i] = i * atoms_per_rank;
    rcounts[i] = atoms_per_rank;
  }

  rcounts[comm->nprocs - 1] = atom->natoms - atoms_per_rank * (comm->nprocs - 1);

  for (int i = min; i <= max  ; i++)
  {
    if (check_avl(i)) {
    temp_avl_list[temp_num_avl] = i;
    temp_num_avl++;
    }
  }
  
  // Gather back available positions list
  MPI_Barrier(world);

  MPI_Allreduce(&temp_num_avl, &num_avl, 1, MPI_LONG, MPI_SUM, world);

  MPI_Allgatherv(temp_avl_list, size, MPI_LONG, av_list, rcounts, displs, MPI_LONG, world);

  // Sort the list to have a working availability list
  std::sort(av_list, av_list + atom->natoms);

  delete temp_avl_list;
}

/*-------------*/
/*Load i-th SMC*/
/*-------------*/
void FixSMC::load_smc(long i) {
  compile_avl_list();
  int npol = atom->natoms / lpol;

  // Distribute the SMCs in the system, placing at least one for each polymer
  if ((initmode == 1) && (i * lpol < atom -> natoms) && (update -> ntimestep  == 1)) {
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
  // Distribute randomly extruders among polymers
  else if ((initmode == 3) && (update -> ntimestep  == 1)) {
    anch[i] = fixed_hinge_position; // this is so manual omg
  }

  else{
    anch[i] = av_list[static_cast < int > (random_equal -> uniform() * num_avl)];
  }  

  if (num_avl == 0) error -> all(FLERR, "Not enough space for the smcs");

  if (hdir != 0) hing[i] = anch[i] + 2 * hdir / abs(hdir);
  else hing[i] = anch[i] - 2 * adir / abs(adir);
  

  MPI_Barrier(world);
  // Cast the chosen position to each processor
  MPI_Bcast(anch, smcnum, MPI_LONG, 0, world);
  MPI_Bcast(hing, smcnum, MPI_LONG, 0, world);

}

void FixSMC::place_smc(long a, long h, bool newsmc) {
  long mhi;
  long man;

  long idhi = atom -> map(h);
  long idan = atom -> map(a);

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
      atom -> bond_atom[mhi][atom -> num_bond[mhi]] = a;
      atom -> num_bond[mhi]++;
    }
  }

}

void FixSMC::remove_smc(long a, long h) {
  long mhi;
  long man;

  long idhi = atom -> map(h);
  long idan = atom -> map(a);

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
  if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Writing restart for fix_smc \n");

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

  if ((debug)) utils::logmesg(lmp, "End of writing restart for fix_smc \n");

}

/*---------------------------------------------------*/
/*Use state info from restart file to restart the Fix*/
/*---------------------------------------------------*/
void FixSMC::restart(char * buf) {
  if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Reading restart for fix_smc \n");

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

  if ((comm->me==0) && (debug)) utils::logmesg(lmp, "End of reading restart for fix_smc \n");

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
