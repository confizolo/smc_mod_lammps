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

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_smc[] =
  "fix smc command:\n\n";

/* ---------------------------------------------------------------------- */

FixSMC::FixSMC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  anch(0),hing(0), smctype(0), smcbtype(0), active(0), debug(1), dir(0), type(nullptr), x(nullptr), list(nullptr), random(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_smc);

  if (narg != 7) error->all(FLERR,"Illegal fix smc command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix smc command");

  next_reneighbor = -1;

  // initialize Marsaglia RNG with processor-unique seed

  int seed = utils::inumeric(FLERR,arg[4],false,lmp);
  random = new RanMars(lmp,seed + comm->me);
  
  smctype = utils::inumeric(FLERR,arg[5],false,lmp);
  if (smctype <= 0) error->all(FLERR,"Illegal fix smc command");

  smcbtype = utils::inumeric(FLERR,arg[6],false,lmp);
  if (smcbtype <= 0) error->all(FLERR,"Illegal fix smc command");
}

/* ---------------------------------------------------------------------- */

FixSMC::~FixSMC()
{
  delete random;

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

  // pair and bonds must be defined
  // no dihedral or improper potentials allowed
  // special bonds must be 0 1 1

  if (force->pair == nullptr || force->bond == nullptr)
    error->all(FLERR,"Fix smc requires pair and bond styles");

  if (force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support fix smc");

  if (force->angle == nullptr && atom->nangles > 0 && comm->me == 0)
    error->warning(FLERR,"Fix smc will not preserve correct angle "
                   "topology because no angle_style is defined");

  // need a half neighbor list, built every Nevery steps

  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);

  //Initialize a random SMC within 5 beads of distance

  // Define a random anchor position inside a monodisperse system with L=1000

  anch = static_cast<int> (random->uniform() * atom->natoms);
  
  if (anch%1000 >= 995){hing = anch-5;}
  else { hing = anch+5;}

  // Assign positive direction to SMC (hinge will move to higher tag beads)
  dir = 1;
}

/* ---------------------------------------------------------------------- */

void FixSMC::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   look for and perform swaps
   NOTE: used to do this every pre_neighbor(), but think that is a bug
         b/c was doing it after exchange() and before neighbor->build()
         which is when neigh lists are actually out-of-date or even bogus,
         now do it based on user-specified Nevery, and trigger reneigh
         if any swaps performed, like fix bond/create
------------------------------------------------------------------------- */

void FixSMC::post_integrate()
{
  int i,j,inum,jnum;
  int inext,ibond, ibondtype;
  int *ilist,*jlist,*numneigh,**firstneigh;

  if (update->ntimestep % nevery) return;

  else if (update->ntimestep == nevery){

  // Change type to defined hinge and anchor beads if in processor
  for (int l = 0; l < atom->nlocal; l++)
  {

    if (atom->tag[l]==hing){
      atom->type[l] = smctype;
      // Creating new SMC bond
      if (atom->num_bond[l] == atom->bond_per_atom) error->one(FLERR, "New bond exceeded bonds per atom limit of {} in create_bonds", atom->bond_per_atom);
      atom->bond_type[l][atom->num_bond[l]] = smcbtype;
      atom->bond_atom[l][atom->num_bond[l]] = anch;
      atom->num_bond[l]++;
    }
    else if (atom->tag[l]==anch) {
      atom->type[l] = smctype;
    }

  }
    return;
  }
  
  else{
  // local ptrs to atom arrays

  tagint *tag = atom->tag;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int nlocal = atom->nlocal;
  int n = -1;

  auto histories = modify->get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  type = atom->type;
  x = atom->x;

  neighbor->build_one(list,1);
  inum = list->inum;
  ilist = list->ilist;
  
  // The hinge must be inside the current processor to have a successful bond/type change

  int active=0;
  for (int ii = 0; ii < inum; ii++)
  {
    i = ilist[ii];
    if (tag[i]==hing) {
      n = num_bond[i];
      for (ibond = 0; ibond < n; ibond++) {
        inext = atom->map(bond_atom[i][ibond]);
        if ((bond_type[i][ibond] == 1) || (tag[inext]>hing)){
        // The operation requires the atom to be on the same processor
        if (inext >= nlocal || inext < 0) return;

        active=1;

        atom->type[i]=1;
        atom->type[inext]=smctype;

        // Deleting old SMC bond
        for (ibond = 0; ibond < n; ibond++) {
          if (bond_type[i][ibond] == smcbtype){
            atom->bond_type[i][ibond] = atom->bond_type[i][n-1];
            atom->bond_atom[i][ibond] = atom->bond_atom[i][n-1];

            if (n_histories > 0)
              for (auto &ihistory: histories) {
                dynamic_cast<FixBondHistory *>(ihistory)->shift_history(i,ibond,n-1);
                dynamic_cast<FixBondHistory *>(ihistory)->delete_history(i,n-1);
                }

            atom->num_bond[i]--;
            break;
          }
        }

        // Creating new SMC bond
        if (num_bond[inext] == atom->bond_per_atom) error->one(FLERR, "New bond exceeded bonds per atom limit of {} in create_bonds", atom->bond_per_atom);
        bond_type[inext][num_bond[inext]] = smcbtype;
        bond_atom[inext][num_bond[inext]] = anch;
        num_bond[inext]++;

        hing = tag[inext];

        break;
        }
      }
  }
  }

  if (!active) return;


  // trigger immediate reneighboring if a movement occurred on one or more procs

  int active_any;
  MPI_Allreduce(&active,&active_any,1,MPI_INT,MPI_SUM,world);
  if (active_any) next_reneighbor = update->ntimestep;

  if (comm->me == 0 && screen) {
    if (debug) {
      fmt::print(screen,"  Current anchor    : {}"
                          "  Current hinge  : {}\n",
                  anch,hing);
    }
  }

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


