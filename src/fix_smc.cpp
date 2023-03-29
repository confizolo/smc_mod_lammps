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
  anch(nullptr), hing(nullptr), smctype(0), smcbtype(0), smcnum(0), debug(0), connFix(nullptr) {
    if (lmp -> citeme) lmp -> citeme -> add(cite_fix_smc);
    // Number of arguments for the fix. The first three arguments are parsed by Fix base class constructor.
    // The rest are specific to this fix. 11 are mandatory
    // 4. nevery: Attempt the jump every nevery iteration
    // 5. seed: random seed
    // 6. prob: probability to attempt the jump
    // 7. lpol: length of polymer(s)
    // 8. poltype: give the form of the polymer (either linear or ring)
    // 9. adir: attempted movement of anchor (next attempted atom id: current anchor + adir)
    // 10. hdir: attempted movement of hinge (next attempted atom id: current hinge + hdir)
    // 11. smcnum: number of deployed smcs
    // 12. smctype: atom type of anchoring beads
    // 13. smcbtype: bond type of anchoring beads after the first deployment
    // 14. smcbitype: bond type of anchoring beads at the first deployment
    // 15. cutoff: distance cutoff for attempted movements (Jump is accepted only if distance between new anchor and hinge is below the cutoff)
    // 16. initmode: random or distributed according to uswe
    // 17. kon: loading probability 
    // 18. koff: unloading probability
    // 19. FixID: Name of ID to get informations about 

    if ((narg != 18) && (narg != 19)) error -> all(FLERR, "Illegal fix smc command");

    nevery = utils::inumeric(FLERR, arg[3], false, lmp);
    if (nevery <= 0) error -> all(FLERR, "Illegal fix smc command");

    // Flag to activate dump in restart file of fix smc structure
    restart_global = 1;

    // Activate flag for array and scalar returning
    scalar_flag = 1;
    array_flag = 1;

    seed = utils::inumeric(FLERR, arg[4], false, lmp);

    prob = utils::numeric(FLERR, arg[5], false, lmp);
    if (prob < 0.0 || prob > 1.0)
      error -> all(FLERR, "Illegal fix topo2 command");

    lpol = utils::inumeric(FLERR, arg[6], false, lmp);
    if (lpol <= 0) error -> all(FLERR, "Illegal fix smc command");

    if (strcmp(arg[7], "linear") == 0) {
      ring = 0;
    } else if (strcmp(arg[7], "ring") == 0) {
      ring = 1;
    } else {
      error -> all(FLERR, "Illegal fix smc command, indefinite polymer type");
    }

    adir = utils::inumeric(FLERR, arg[8], false, lmp);

    hdir = utils::inumeric(FLERR, arg[9], false, lmp);
    if (hdir * adir >= 0) error -> all(FLERR, "Illegal fix smc command, hinge and must not have same direction");

    smcnum = utils::inumeric(FLERR, arg[10], false, lmp);
    if (smcnum <= 0) error -> all(FLERR, "Illegal fix smc command");

    smctype = utils::inumeric(FLERR, arg[11], false, lmp);
    if (smctype <= 0) error -> all(FLERR, "Illegal fix smc command");

    smcbtype = utils::inumeric(FLERR, arg[12], false, lmp);
    if (smcbtype <= 0) error -> all(FLERR, "Illegal fix smc command");

    smcbitype = utils::inumeric(FLERR, arg[13], false, lmp);
    if (smcbitype <= 0) error -> all(FLERR, "Illegal fix smc command");

    cutoff = utils::numeric(FLERR, arg[14], false, lmp);
    if (cutoff < 0)
      error -> all(FLERR, "Illegal fix smc command");

    initmode = 0;

    if (strcmp(arg[15], "random") == 0) {
      initmode = 0;
    } else if (strcmp(arg[15], "distributed") == 0) {
      initmode = 1;
    } else if (strcmp(arg[15], "full-distributed") == 0) {
      initmode = 2;
    } else {
      error -> all(FLERR, "Illegal fix smc command, initmode not present");
    }

    kon = utils::numeric(FLERR, arg[16], false, lmp);
    if ((kon <= 0) || (kon>1))
      error -> all(FLERR, "Illegal fix smc command, kon is out of the interval (0,1]");

    koff = utils::numeric(FLERR, arg[17], false, lmp);
    if ((koff < 0) || (koff>1))
      error -> all(FLERR, "Illegal fix smc command, koff is out of the interval [0,1]");

    if (narg == 19) {
      connFixName = new char[static_cast < int > (sizeof(arg[18]) / sizeof(char))];
      std::copy(arg[18], arg[18] + static_cast < int > (sizeof(arg[18]) / sizeof(char)), connFixName);
    } else {
      connFixName = new char[5];
      connFixName = "nofix";
    }

    xyzanch = nullptr;
    xyzhing = nullptr;

    anch = new long[smcnum];
    hing = new long[smcnum];

    for (int i = 0; i < smcnum; i++) {
      hing[i] = -1;
      anch[i] = -1;
    }

    av_list = new long[atom->natoms];
    num_avl = atom->natoms;

    random_equal = new RanPark(lmp, seed);

    // To get a different random number every time the program is executed
    srand(time(NULL) * seed);
  }

/* ---------------------------------------------------------------------- */

FixSMC::~FixSMC() {
  int m;
  int man;

  // Get bond histories to apply bond changes
  auto histories = modify -> get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  int idhi;
  int idan;

  // Loop over the SMCs instantiated and remove
  for (int i = 0; i < smcnum; i++) {
    remove_smc(anch[i], hing[i]);
  }

  delete random_equal;
  delete anch;
  delete hing;

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

  if (force -> pair == nullptr || force -> bond == nullptr)
    error -> all(FLERR, "Fix smc requires pair and bond styles");

  // Store the connected Fix ID pointer
  if (utils::strmatch(connFixName, "nofix") == 1) {
    connFix = nullptr;
  } else {
    connFix = modify -> get_fix_by_id(connFixName);
    if (!connFix) error -> all(FLERR, "Illegal ausiliary Fix");
  }

}

/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

void FixSMC::post_integrate() {

  if (update -> ntimestep == 1) {

    // Draw two random numbers for the unloading/loading
    double lrand;

    // Define a random anchor position inside a monodisperse system with L=lpol
    // Define randomly the position of the smc hinge and anchors along one of the polymers
    for (int i = 0; i < smcnum; i++)
    {
    if (comm -> me == 0) lrand = random_equal -> uniform();
    MPI_Bcast( & lrand, 1, MPI_DOUBLE, 0, world);

    if (lrand < kon) load_smc(i);
    }
    
    MPI_Barrier(world);

    // for (int i = 0; i < smcnum; i++)
    // { 
    //   if ((debug)) utils::logmesg(lmp, "Anchors are" + std::to_string(anch[i]) + " " + std::to_string(hing[i]) + "\n");
    // }

    // Create bonds according to chosen beads
    for (int i = 0; i < smcnum; i++) {
      // Check if smc are loaded
      if ((anch[i]<0) || (hing[i]<0)) continue;
      place_smc(anch[i],hing[i], true);
    }

    return;

  } 

  else {
    for (int i = 0; i < smcnum; i++) {

      // Draw two random numbers for the unloading/loading
      double lrand;

      if (comm -> me == 0) lrand = random_equal -> uniform();

      MPI_Bcast( & lrand, 1, MPI_DOUBLE, 0, world);

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
    }
  }

  if (update -> ntimestep % nevery == 0){

    // Return if the smcs are still
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

    int tempadir = adir;
    int temphdir = hdir;
    bool flag = 0;

    // for (int i = 0; i < smcnum; i++)
    // { 
    //   if ((debug)) utils::logmesg(lmp, "Anchors are" + std::to_string(anch[i]) + " " + std::to_string(hing[i]) + "\n");
    // }

    double rand;

    for (int i = 0; i < smcnum; i++) {
      
      // Check if smc are loaded
      if ((anch[i]<0) || (hing[i]<0)) continue;

      // Draw a random number for the jump attempt
      if (comm -> me == 0) rand = random_equal -> uniform();
      MPI_Bcast( & rand, 1, MPI_DOUBLE, 0, world);

      if (rand > prob) {
        continue;
      }

      // Temporary direction if the smc is going towards the polymer end or another smc bead
      tempadir = adir;
      temphdir = hdir;

      // Check if we are going to the polymer border on one side or on the other

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
      
      if (connFix) {
        // Check for conflicts with connected fix
        for (int k = 0; k < connFix -> compute_scalar(); k++) {
          if (((anch[i] + adir) >= connFix -> compute_array(k, 0)) && ((anch[i] + adir) < (connFix -> compute_array(k, 1)))) {
            tempadir = 0;
          }
          if (((hing[i] + hdir) >= connFix -> compute_array(k, 0)) && ((hing[i] + hdir) < (connFix -> compute_array(k, 1)))) {
            temphdir = 0;
          }
        }
      }

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

      if (flag) {
        continue;
      }

      // if ((debug)) utils::logmesg(lmp, "SMC {} Current anchor {}, current hinge {}; next anchor {}, next hinge {} \n", i,anch[i],hing[i],anch[i]+tempadir,hing[i]+temphdir);

      idnewhi = atom -> map(hing[i] + temphdir);
      idhi = atom -> map(hing[i]);
      idnewan = atom -> map(anch[i] + tempadir);
      idan = atom -> map(anch[i]);

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

      // Computing the distance between the new beads in a parallel way
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

      // if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Number of counts is " + std::to_string(anchcounts[0]) + " Anchor " + std::to_string(hingcounts[0]) + " Hinge " + "\n");
      // if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Proposed distance is " + std::to_string(sqrt(dist)) + "\n");

      // Check if the distance is small enough to run the jump
      if (!(dist > cutoff * cutoff || (anchcounts[0] == 0) || (hingcounts[0] == 0))) {
        
        if ((temphdir!=0) || (tempadir!=0)){
          remove_smc(anch[i], hing[i]);
          place_smc(anch[i] + tempadir, hing[i] + temphdir, false);
        }

        anch[i] += tempadir;
        hing[i] += temphdir;

      }
      // Barrier to check that each processor has defined correctly each smc
      MPI_Barrier(world);
    }

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

/* ----------------------------------------------------------------------
   memory usage of hing and anch arrays
------------------------------------------------------------------------- */

double FixSMC::memory_usage() {
  double bytes = 2 * smcnum * sizeof(long) + atom->natoms * sizeof(long);
  return bytes;
}

/*---------------------------------------------------------------------*/
/* Needed to write a restart file that can continue with the simulation*/
/*---------------------------------------------------------------------*/
void FixSMC::write_restart(FILE * fp) {
  // if ((debug)) utils::logmesg(lmp, "Writing restart for fix_smc \n");

  int rn = 0;
  long rlist[3 + 2 * smcnum];

  rlist[rn++] = static_cast < long > (next_reneighbor);
  rlist[rn++] = static_cast < long > (update -> ntimestep);

  rlist[rn++] = static_cast < long > (smcnum);

  // Saving SMCs positions
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

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */
void FixSMC::restart(char * buf) {
  if ((debug)) utils::logmesg(lmp, "Reading restart for fix_smc \n");

  int rn = 0;
  long * rlist = (long * ) buf;

  next_reneighbor = static_cast < bigint > (rlist[rn++]);

  bigint ntimestep_restart = static_cast < bigint > (rlist[rn++]);
  if (ntimestep_restart != update -> ntimestep)
    error -> all(FLERR, "Must not reset timestep when restarting fix smc");

  int smcnum_rest = rlist[rn++];
  if (smcnum_rest != smcnum)
    error -> all(FLERR, "Invalid restart, number of smcs has changed!");

  // Loading SMCs positions
  for (int j = 0; j < smcnum; j++) {
    anch[j] = static_cast < long > (rlist[rn++]);
    hing[j] = static_cast < long > (rlist[rn++]);
  }

  if ((debug)) utils::logmesg(lmp, "End of reading restart for fix_smc \n");

}

/*
Returns number of smcs
*/
double FixSMC::compute_scalar() {
  return smcnum;
}

/*
Returns position of i smc hinge or anchor depending on the flag 
*/
double FixSMC::compute_array(int i, int flag) {
  int rflag;
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
bool FixSMC::check_avl(long i){
  long tmphing;
  long mdbead;
  bool flag = 0;

  // Instantiate the bead according to the direction
  if (hdir != 0) tmphing = i + 2 * hdir / abs(hdir);
  else tmphing = i - 2 * adir / abs(adir);

  mdbead = (i + tmphing)/2;

  // Check if the smc is wrongly positioned (border conditions)
  if ((mdbead%lpol == 0) || (mdbead == 1)) {
    return 0;
  }

  flag = 0;
  // Check if we are superimposing other beads
  for (int j = 0; j < smcnum; j++) {
    if ((anch[j]<0) || (hing[j]<0)) continue;
    if (((i == anch[j]) || (tmphing == hing[j])) || ((i == hing[j]) || (tmphing == anch[j])) || ((i + 1) == anch[j]) || ((i - 1) == anch[j]) || ((tmphing + 1) == anch[j]) || ((tmphing - 1) == anch[j]) || ((i + 1) == hing[j]) || ((i - 1) == hing[j]) || ((tmphing + 1) == hing[j]) || ((tmphing - 1) == hing[j])) {
      flag = 1;
      break;
    }
  }

  if (flag) return 0;

  if (connFix) {
    // Check for conflicts with connected fix
    for (int k = 0; k < connFix -> compute_scalar(); k++) {
      if (connFix -> compute_array(k, 1) > atom -> natoms) continue;
      if ((i >= connFix -> compute_array(k, 0)) && (i < (connFix -> compute_array(k, 1)))) {
        flag = 1;
        break;
      }
      if ((tmphing >= connFix -> compute_array(k, 0)) && (tmphing < (connFix -> compute_array(k, 1)))) {
        flag = 1;
        break;
      }
    }
  }

  if (flag) return 0;

  return 1;
}

void FixSMC::compile_avl_list(){
  num_avl = 0;

  if (comm->me==0) {
    for (int i = 1; i <= atom->natoms; i++)
    {
      if (check_avl(i)) {
      av_list[num_avl] = i;
      num_avl++;
      }
    }
  }

  MPI_Bcast(av_list, atom->natoms, MPI_LONG, 0, world);
  MPI_Bcast( &num_avl, 1, MPI_DOUBLE, 0, world);
}

void FixSMC::load_smc(long i) {
  int npol = atom->natoms % lpol;
  if (num_avl == 0) error -> all(FLERR, "Not enough space for the smcs");
  if (comm->me==0) {
    if ((initmode == 1) && (i * lpol < atom -> natoms) && (update -> ntimestep  == 1)) {
      do {
      anch[i] = static_cast < int > (random_equal -> uniform() * lpol + i * lpol);
      } while (not check_avl(anch[i]));
    }
    else if ((initmode == 2) && (update -> ntimestep  == 1)) {
      do {
        anch[i] = static_cast < int > (random_equal -> uniform() * lpol + (i%npol) * lpol);
      } while (not check_avl(anch[i]));
    }
    else{
      compile_avl_list();
      anch[i] = av_list[static_cast < int > (random_equal -> uniform() * num_avl)];
    }  
  }

  if (hdir != 0) hing[i] = anch[i] + 2 * hdir / abs(hdir);
  else hing[i] = anch[i] - 2 * adir / abs(adir);

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
    // Changing type
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