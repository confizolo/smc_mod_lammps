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

#include "fix_smc.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "citeme.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "random_mars.h"
#include "update.h"
#include "create_bonds.h"
#include "delete_bonds.h"
#include <iostream>
#include <cmath>
#include <cstring>
#include <utils.h>

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_smc[] =
  "fix smc command:\n\n";

/* ---------------------------------------------------------------------- */

FixSMC::FixSMC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  anch(0),hing(0), smctype(0), smcbtype(0), debug(0), list(nullptr), random(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_smc);

  if (narg != 13) error->all(FLERR,"Illegal fix smc command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix smc command");

  next_reneighbor = -1;


  int seed = utils::inumeric(FLERR,arg[4],false,lmp);
  random = new RanMars(lmp,seed);
  
  prob = utils::numeric(FLERR, arg[5], false, lmp);
  if (prob < 0.0 || prob > 1.0)
            error->all(FLERR, "Illegal fix topo2 command");

  lpol = utils::inumeric(FLERR,arg[6],false,lmp);
  if (lpol <= 0) error->all(FLERR,"Illegal fix smc command");

  adir = utils::inumeric(FLERR,arg[7],false,lmp);

  hdir = utils::inumeric(FLERR,arg[8],false,lmp);
  if (hdir*adir >= 0) error->all(FLERR,"Illegal fix smc command, hinge and must not have same direction");

  smctype = utils::inumeric(FLERR,arg[9],false,lmp);
  if (smctype <= 0) error->all(FLERR,"Illegal fix smc command");

  smcbtype = utils::inumeric(FLERR,arg[10],false,lmp);
  if (smcbtype <= 0) error->all(FLERR,"Illegal fix smc command");

  smcbitype = utils::inumeric(FLERR,arg[11],false,lmp);
  if (smcbitype <= 0) error->all(FLERR,"Illegal fix smc command");

  cutoff = utils::numeric(FLERR, arg[12], false, lmp);
  if (cutoff <0)
            error->all(FLERR, "Illegal fix topo2 command");

  // To get a different random number every time the program is executed
  srand(time(NULL) * seed);
  
  xyzanch = nullptr;
  xyzhing = nullptr;
}

/* ---------------------------------------------------------------------- */

FixSMC::~FixSMC()
{
  int m;
  int man;

  auto histories = modify->get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  const int idhi = atom->map(hing);
  const int idan = atom->map(anch);

  if (((m = idhi) >= 0) && (hdir!=0)){

    atom->type[m]=1;

    // Deleting old SMC bond
    for (int ibond = 0; ibond < atom->num_bond[m]; ibond++) {
      if ((atom->bond_type[m][ibond] == smcbtype) || (atom->bond_type[m][ibond] == smcbitype)){
        atom->bond_type[m][ibond] = atom->bond_type[m][atom->num_bond[m]-1];
        atom->bond_atom[m][ibond] = atom->bond_atom[m][atom->num_bond[m]-1];

        if (n_histories > 0)
          for (auto &ihistory: histories) {
            dynamic_cast<FixBondHistory *>(ihistory)->shift_history(m,ibond,atom->num_bond[m]-1);
            dynamic_cast<FixBondHistory *>(ihistory)->delete_history(m,atom->num_bond[m]-1);
            }

        atom->num_bond[m]--;
        break;
      }
    }
  }

  if ((man = idan) >= 0) {
      atom->type[man] = 1;
  }

  delete random;
  memory->destroy(xyzanch);
  memory->destroy(xyzhing);
  memory->destroy(list);

}

/* ---------------------------------------------------------------------- */

int FixSMC::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSMC::init()
{


  if (force->pair == nullptr || force->bond == nullptr)
    error->all(FLERR,"Fix smc requires pair and bond styles");

}


void FixSMC::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

void FixSMC::post_integrate()
{
  int m;
  int man;

  if (update->ntimestep % nevery) return;

  else if (update->ntimestep == nevery){
  //Initialize a random SMC within 5 beads of distance

  // Define a random anchor position inside a monodisperse system with L=1000
  anch = static_cast<int> (random->uniform() * lpol);
    
  if (hdir!=0) hing = anch + 2*hdir/abs(hdir);
  else hing = anch - 2*adir/abs(adir);

  const int idhi = atom->map(hing);
  const int idan = atom->map(anch);

  // Change type to defined hinge and anchor beads if in processor
  if ((m = idhi) >= 0){

    atom->type[m] = smctype;
    // Creating new SMC bond
    if (atom->num_bond[m] == atom->bond_per_atom) error->one(FLERR, "New bond exceeded bonds per atom limit of {} in create_bonds", atom->bond_per_atom);
    atom->bond_type[m][atom->num_bond[m]] = smcbitype;
    atom->bond_atom[m][atom->num_bond[m]] = anch;
    atom->num_bond[m]++;
  }
  if ((man = idan) >= 0) {
      atom->type[man] = smctype;
  }
  return;
  }
  
  else{

  if (random->uniform() > prob) return;
  
  // Check if we are going to the polymer border
  if ((hing + hdir)%lpol == 0) hdir = 0;  
  if ((anch + adir)%lpol == 0) adir = 0;  

  if ((hdir==0) && (adir==0)) return;
  
  int mnew;
  int mannew;

  const int idnewhi = atom->map(hing+hdir);
  const int idhi = atom->map(hing);
  const int idnewan = atom->map(anch+adir);
  const int idan = atom->map(anch);

  double *xyzanchtemp = nullptr;
  double *xyzhingtemp = nullptr;

  int *anchcount= nullptr;
  int *hingcount= nullptr;

  memory->destroy(xyzanchtemp);
  memory->create(xyzanchtemp,3,"FixSMC::post_integrate()");
  memory->destroy(xyzhingtemp);
  memory->create(xyzhingtemp,3,"FixSMC::post_integrate()");
  memory->destroy(anchcount);
  memory->create(anchcount,1,"FixSMC::post_integrate()");
  memory->destroy(hingcount);
  memory->create(hingcount,1,"FixSMC::post_integrate()");

  anchcount[0] = 0;
  hingcount[0] = 0;

  xyzanchtemp[0] = 0;
  xyzanchtemp[1] = 0;
  xyzanchtemp[2] = 0;

  xyzhingtemp[0] = 0;
  xyzhingtemp[1] = 0;
  xyzhingtemp[2] = 0;

  memory->destroy(xyzanch);
  memory->create(xyzanch,3,"FixSMC::post_integrate()");
  memory->destroy(xyzhing);
  memory->create(xyzhing,3,"FixSMC::post_integrate()");
  
  double unwrap[3];


  if (((mannew = idnewan) >= 0) && (idnewan<(atom->nlocal))){
    domain->unmap(atom->x[mannew], atom->image[mannew], unwrap);
    xyzanchtemp[0] += unwrap[0];
    xyzanchtemp[1] += unwrap[1];
    xyzanchtemp[2] += unwrap[2];
    anchcount[0]++;
  }

  if (((mnew = idnewhi) >= 0) && (idnewhi<(atom->nlocal))){
    domain->unmap(atom->x[mnew], atom->image[mnew], unwrap);

    xyzhingtemp[0] += unwrap[0];
    xyzhingtemp[1] += unwrap[1];
    xyzhingtemp[2] += unwrap[2];
    hingcount[0]++;
  } 

  int* hingcounts = nullptr;
  int* anchcounts = nullptr;

  memory->destroy(anchcounts);
  memory->create(anchcounts,1,"FixSMC::post_integrate()");
  memory->destroy(hingcounts);
  memory->create(hingcounts,1,"FixSMC::post_integrate()");

  MPI_Allreduce(xyzanchtemp, xyzanch, 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(xyzhingtemp, xyzhing, 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(anchcount, anchcounts, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(hingcount, hingcounts, 1, MPI_DOUBLE, MPI_SUM, world);

  double dist = 0;

  for (int k=0; k<3 ; k++){
    xyzanch[k] = xyzanch[k]/anchcounts[0];
    xyzhing[k] = xyzhing[k]/hingcounts[0];
    dist += (xyzanch[k]-xyzhing[k])*(xyzanch[k]-xyzhing[k]);
  }
  if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Number of counts is " + std::to_string(anchcounts[0]) + " Anchor " + std::to_string(hingcounts[0]) + " Hinge " + "\n");
  if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Proposed distance is " + std::to_string(sqrt(dist)) + "\n");
  if (dist > cutoff*cutoff || (anchcounts[0]==0) || (hingcounts[0]==0)) return;

  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;

  auto histories = modify->get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  if (((m = idhi) >= 0) && (hdir!=0)){

    atom->type[m]=1;

    // Deleting old SMC bond
    for (int ibond = 0; ibond < atom->num_bond[m]; ibond++) {
      if ((bond_type[m][ibond] == smcbtype) || (bond_type[m][ibond] == smcbitype)){
        atom->bond_type[m][ibond] = atom->bond_type[m][num_bond[m]-1];
        atom->bond_atom[m][ibond] = atom->bond_atom[m][num_bond[m]-1];

        if (n_histories > 0)
          for (auto &ihistory: histories) {
            dynamic_cast<FixBondHistory *>(ihistory)->shift_history(m,ibond,atom->num_bond[m]-1);
            dynamic_cast<FixBondHistory *>(ihistory)->delete_history(m,atom->num_bond[m]-1);
            }

        atom->num_bond[m]--;
        break;
      }
    }
  }

  if (((man = idan) >= 0) && (adir!=0)){
    atom->type[man] = 1;
  }

  if (((mannew = idnewan) >= 0) && (adir!=0)){
    atom->type[mannew] = smctype;
  }

  if (((mnew = idnewhi) >= 0) && (hdir!=0)){
    atom->type[mnew]=smctype;

    // Creating new SMC bond
    if (num_bond[mnew] == atom->bond_per_atom) error->one(FLERR, "New bond exceeded bonds per atom limit of {} in create_bonds", atom->bond_per_atom);
    bond_type[mnew][num_bond[mnew]] = smcbtype;
    bond_atom[mnew][num_bond[mnew]] = (anch+adir);
    num_bond[mnew]++;

    }

    anch+=adir;
    hing+=hdir;

    memory->destroy(xyzanch);
    memory->destroy(xyzhing);
    memory->destroy(xyzanchtemp);
    memory->destroy(xyzhingtemp);
    memory->destroy(anchcount);
    memory->destroy(hingcount);
    memory->destroy(hingcounts);
    memory->destroy(anchcounts);

    return;
  } 
  
  }


/* ----------------------------------------------------------------------
   memory usage of alist
------------------------------------------------------------------------- */

double FixSMC::memory_usage()
{
  double bytes = 0;
  return bytes;
}


/***********************************************************************/
/* Needed to write a restart file that can continue with the simulation*/
/***********************************************************************/
void FixSMC::write_restart(FILE *fp)
{
    int n = 0;
    double list[4];
    list[n++] = ubuf(next_reneighbor).d;
    list[n++] = ubuf(update->ntimestep).d;
    list[n++] = ubuf(anch).d;
    list[n++] = ubuf(hing).d;

    if (comm->me == 0)
    {
        int size = n * sizeof(double);
        fwrite(&size, sizeof(int), 1, fp);
        fwrite(list, sizeof(double), n, fp);
    }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */
void FixSMC::restart(char *buf)
{
    int n = 0;
    double *list = (double *)buf;

    next_reneighbor = (bigint)ubuf(list[n++]).i;

    bigint ntimestep_restart = (bigint)ubuf(list[n++]).i;
    if (ntimestep_restart != update->ntimestep)
        error->all(FLERR, "Must not reset timestep when restarting fix smc");

    anch = (bigint)ubuf(list[n++]).i;

    hing = (bigint)ubuf(list[n++]).i;

}

