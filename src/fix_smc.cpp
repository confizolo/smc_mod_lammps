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
#include <cmath>
#include <cstring>
#include <utils.h>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_smc[] =
  "fix smc command:\n\n";

/* ---------------------------------------------------------------------- */

FixSMC::FixSMC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  anch(nullptr),hing(nullptr), smctype(0), smcbtype(0), smcnum(0), debug(0), random(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_smc);
  // Number of arguments for the fix. The first three arguments are parsed by Fix base class constructor.
  // The rest are specific to this fix. 11 are mandatory
  // 4. nevery: Attempt the jump every nevery iteration
  // 5. seed: random seed
  // 6. prob: probability to attempt the jump
  // 7. lpol: length of polymer(s)
  // 8. adir: attempted movement of anchor (next attempted atom id: current anchor + adir)
  // 9. hdir: attempted movement of hinge (next attempted atom id: current hinge + hdir)
  // 10. smcnum: number of deployed smcs
  // 11. smctype: atom type of anchoring beads
  // 12. smcbtype: bond type of anchoring beads after the first deployment
  // 13. smcbitype: bond type of anchoring beads at the first deployment
  // 14. cutoff: distance cutoff for attempted movements (Jump is accepted only if distance between new anchor and hinge is below the cutoff)
  // 15. FixID: Name of ID to get informations about 

  if ((narg != 14) && (narg!=15)) error->all(FLERR,"Illegal fix smc command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix smc command");

  // Flag to activate dump in restart file of fix smc structure
  restart_global = 1;

  // Activate flag for array and scalar returning
  scalar_flag = 1;
  array_flag = 1;

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

  smcnum = utils::inumeric(FLERR,arg[9],false,lmp);
  if (smcnum <= 0) error->all(FLERR,"Illegal fix smc command");

  smctype = utils::inumeric(FLERR,arg[10],false,lmp);
  if (smctype <= 0) error->all(FLERR,"Illegal fix smc command");

  smcbtype = utils::inumeric(FLERR,arg[11],false,lmp);
  if (smcbtype <= 0) error->all(FLERR,"Illegal fix smc command");

  smcbitype = utils::inumeric(FLERR,arg[12],false,lmp);
  if (smcbitype <= 0) error->all(FLERR,"Illegal fix smc command");

  cutoff = utils::numeric(FLERR, arg[13], false, lmp);
  if (cutoff <0)
            error->all(FLERR, "Illegal fix topo2 command");

  if (narg==15){
    connFixName = new char[static_cast<int>(sizeof(arg[14]) / sizeof(char))];
    std::copy(arg[14],arg[14]+static_cast<int>(sizeof(arg[14])/sizeof(char)),connFixName);
  }
  else
  {
        connFixName = new char[5];
        connFixName = "nofix";
  }

  // To get a different random number every time the program is executed
  srand(time(NULL) * seed);
  
  xyzanch = nullptr;
  xyzhing = nullptr;

  anch = new long[smcnum];
  hing = new long[smcnum];

}

/* ---------------------------------------------------------------------- */

FixSMC::~FixSMC()
{
  int m;
  int man;

  // Get bond histories to apply bond changes
  auto histories = modify->get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  int idhi;
  int idan;

  // Loop over the SMCs instantiated
  for (int i = 0; i < smcnum; i++)
  {
    idhi = atom->map(hing[i]);
    idan = atom->map(anch[i]);

    // Check if hinge is in the current processor and delete bonds/reset to original type
    if (((m = idhi) >= 0) && (idhi<atom->nlocal)){

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
    // Check if anchor is in the current processor and reset the type
    if (((man = idan) >= 0) && (idan<atom->nlocal)) {
        atom->type[man] = 1;
    }

  }
  
  delete random;
  delete anch;
  delete hing;
  delete connFixName;
  memory->destroy(xyzanch);
  memory->destroy(xyzhing);  

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

  // Store the connected Fix ID pointer
  if (strcmp(connFixName, "nofix") == 0){
    connFix == nullptr;
  }
  else{
    connFix = modify->get_fix_by_id(connFixName);
    if (!connFix) error->all(FLERR, "Illegal ausiliary Fix");
  }

  for (int i = 0; i < smcnum; i++)
  {
    hing[i] = atom->natoms+1;
    anch[i] = atom->natoms+1;
  }
  
}


/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

void FixSMC::post_integrate()
{
  int m;
  int man;

  if (update->ntimestep % nevery) return;

  else if (update->ntimestep == nevery){

  //Initialize a random SMC within 2 beads of distance
  int idhi;
  int idan;
  
  // Define a random anchor position inside a monodisperse system with L=lpol
  if (comm->me==0){

  // Define randomly the position of the smc hinge and anchors along one of the polymers
  int i=0;
  bool flag=0;

  while (i<smcnum)
  { 
    anch[i] = static_cast<int> (random->uniform() * atom->natoms);
      
    // Check if the smc is wrongly positioned (border conditions)
    if ((anch[i]+1)%lpol==0){anch[i]-=1;}
    if ((anch[i]-1)%lpol==1){anch[i]+=1;}
    if ((anch[i]+1)%lpol==1){anch[i]-=2;}
    if ((anch[i]-1)%lpol==0){anch[i]+=2;}

    // Instantiate the bead according to the direction
    if (hdir!=0) hing[i] = anch[i] + 2*hdir/abs(hdir);
    else hing[i] = anch[i] - 2*adir/abs(adir);

    flag =0;
    // Check if we are superimposing other beads
    for (int j = 0; j < i; j++)
    {
      if(((anch[i]==anch[j])||(hing[i]==hing[j])) || ((anch[i]==hing[j])||(hing[i]==anch[j]))){flag =1; break;}
    }

    if (connFix)
    {
        // Check for conflicts with connected fix
        for (int k = 0; k < connFix->compute_scalar(); k++)
        {
            if (connFix->compute_array(k, 1) > atom->natoms) continue;
            if (((anch[i]) >= connFix->compute_array(k, 0)) && ((anch[i]) < (connFix->compute_array(k, 1))))
            {
              flag =1;
              break;
            }
            if (((hing[i]) >= connFix->compute_array(k, 0)) && ((hing[i]) < (connFix->compute_array(k, 1))))
            {
              flag=1;
              break;
            }
        }
    }

    if(flag){continue;}

    i++;
  }


  }

  // Cast the chosen position to each processor
  MPI_Bcast(anch,smcnum,MPI_LONG,0,world);
  MPI_Bcast(hing,smcnum,MPI_LONG,0,world);

  MPI_Barrier(world);

  // for (int i = 0; i < smcnum; i++)
  // { 
  //   if ((debug)) utils::logmesg(lmp, "Anchors are" + std::to_string(anch[i]) + " " + std::to_string(hing[i]) + "\n");
  // }
    
  // Create bonds according to chosen beads
  for (int i = 0; i < smcnum; i++)
  {
    
    idhi = atom->map(hing[i]);
    idan = atom->map(anch[i]);

    // Change type to defined hinge and anchor beads if in processor
    if (((m = idhi) >= 0) && (idhi<atom->nlocal)){

      atom->type[m] = smctype;
      // Creating new SMC bond
      if (atom->num_bond[m] == atom->bond_per_atom) error->one(FLERR, "New bond exceeded bonds per atom limit of {} in create_bonds", atom->bond_per_atom);
      atom->bond_type[m][atom->num_bond[m]] = smcbitype;
      atom->bond_atom[m][atom->num_bond[m]] = anch[i];
      atom->num_bond[m]++;
    }
    // Change bead anchor type to smctype
    if (((man = idan) >= 0) && (idan<atom->nlocal)) {
        atom->type[man] = smctype;
    }
  }
  
  return;

  }
  
  else{

  // Draw a random number for the jump attempt
  if (random->uniform() > prob) return;
  
  // Return if the smcs are still
  if ((hdir==0) && (adir==0)) return;

  double *xyzanchtemp = nullptr;
  double *xyzhingtemp = nullptr;

  int *anchcount= nullptr;
  int *hingcount= nullptr;

  int* hingcounts = nullptr;
  int* anchcounts = nullptr;

  double unwrap[3];

  int mnew;
  int mannew;

  double dist = 0;

  auto histories = modify->get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  int idnewhi;
  int idhi;
  int idnewan;
  int idan; 

  int tempadir = adir;
  int temphdir = hdir;
  bool flag=0;
  
  // for (int i = 0; i < smcnum; i++)
  // { 
  //   if ((debug)) utils::logmesg(lmp, "Anchors are" + std::to_string(anch[i]) + " " + std::to_string(hing[i]) + "\n");
  // }
    
  for (int i = 0; i < smcnum; i++)
  {
    // Temporary direction if the smc is going towards the polymer end or another smc bead
    tempadir = adir;
    temphdir = hdir;

    // Check if we are going to the polymer border on one side or on the other
    if (hdir/abs(hdir) < 0){
      if ((hing[i] + hdir)%lpol == 0) temphdir = 0;
    }
    else{
      if ((hing[i] + hdir)%lpol == 1) temphdir = 0;
    }
    if (adir/abs(adir) < 0){
      if ((anch[i] + adir)%lpol == 0) tempadir = 0;
    }
    else{
      if ((anch[i] + adir)%lpol == 1) tempadir = 0;
    }

    if (connFix)
    {
        // Check for conflicts with connected fix
        for (int k = 0; k < connFix->compute_scalar(); k++)
        {
            if (((anch[i] + adir)>= connFix->compute_array(k, 0)) && ( (anch[i] + adir)< (connFix->compute_array(k, 1))))
            {
              tempadir = 0;
            }
            if (((hing[i] + hdir)>= connFix->compute_array(k, 0)) && ( (hing[i] + hdir)< (connFix->compute_array(k, 1))))
            {
              temphdir = 0;
            }
        }
    }

    if ((tempadir==0) && (temphdir==0)) continue;

    // Check if the new movement is forbidden because of superposition of SMCs
    flag=0;
    for (int j = 0; j < smcnum; j++)
    {
      if (j==i) continue;
      if((((anch[i]+tempadir)==anch[j])||((anch[i]+tempadir)==hing[j])) && (((hing[i]+temphdir)==hing[j])||((hing[i]+temphdir)==anch[j]))){flag =1; break;}
      else if (((anch[i]+tempadir)==anch[j])||((anch[i]+tempadir)==hing[j])) tempadir=0;
      else if (((hing[i]+temphdir)==hing[j])||((hing[i]+temphdir)==anch[j])) temphdir=0;
    }

    if(flag){continue;}

    // if ((debug)) utils::logmesg(lmp, "SMC {} Current anchor {}, current hinge {}; next anchor {}, next hinge {} \n", i,anch[i],hing[i],anch[i]+tempadir,hing[i]+temphdir);
    
    idnewhi = atom->map(hing[i]+temphdir);
    idhi = atom->map(hing[i]);
    idnewan = atom->map(anch[i]+tempadir);
    idan = atom->map(anch[i]);

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
    
    // Computing the distance between the new beads in a parallel way
    if (((mannew = idnewan) >= 0) && (idnewan<(atom->nlocal))){
      domain->unmap(atom->x[mannew], atom->image[mannew], unwrap);
      xyzanchtemp[0] += unwrap[0];
      xyzanchtemp[1] += unwrap[1];
      xyzanchtemp[2] += unwrap[2];
      anchcount[0]+=1;
    }

    if (((mnew = idnewhi) >= 0) && (idnewhi<(atom->nlocal))){
      domain->unmap(atom->x[mnew], atom->image[mnew], unwrap);

      xyzhingtemp[0] += unwrap[0];
      xyzhingtemp[1] += unwrap[1];
      xyzhingtemp[2] += unwrap[2];
      hingcount[0]+=1;
    } 
    
    memory->destroy(anchcounts);
    memory->create(anchcounts,1,"FixSMC::post_integrate()");
    memory->destroy(hingcounts);
    memory->create(hingcounts,1,"FixSMC::post_integrate()");

    //MPI Barrier to avoid computational errors due to value collection
    MPI_Barrier(world);

    MPI_Allreduce(xyzanchtemp, xyzanch, 3, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(xyzhingtemp, xyzhing, 3, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(anchcount, anchcounts, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(hingcount, hingcounts, 1, MPI_INT, MPI_SUM, world);

    dist = 0;

    for (int k=0; k<3 ; k++){
      xyzanch[k] = xyzanch[k]/anchcounts[0];
      xyzhing[k] = xyzhing[k]/hingcounts[0];
      dist += (xyzanch[k]-xyzhing[k])*(xyzanch[k]-xyzhing[k]);
    }

    // if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Number of counts is " + std::to_string(anchcounts[0]) + " Anchor " + std::to_string(hingcounts[0]) + " Hinge " + "\n");
    // if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Proposed distance is " + std::to_string(sqrt(dist)) + "\n");
    
    // Check if the distance is small enough to run the jump
    if (!(dist > cutoff*cutoff || (anchcounts[0]==0) || (hingcounts[0]==0))) {


    if (((m = idhi) >= 0) && (idhi < atom->nlocal)){
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

      // Move the bond to the new anchor even if the hinge is still if not already present
      if (temphdir==0){
      
      // if ((comm->me==0) && (debug)) utils::logmesg(lmp, "Changing bond keeping hinge fixed \n");

      bool create = 1;
      for (int ibond = 0; ibond < atom->num_bond[m]; ibond++) {
        if ((atom->bond_type[m][ibond] == smcbtype) || (atom->bond_type[m][ibond] == smcbitype)){
          create = 0;
        }
      }

      if (create) {
      // Creating new SMC bond
      if (atom->num_bond[m] == atom->bond_per_atom) error->one(FLERR, "New bond exceeded bonds per atom limit of {} in create_bonds", atom->bond_per_atom);
      atom->bond_type[m][atom->num_bond[m]] = smcbtype;
      atom->bond_atom[m][atom->num_bond[m]] = (anch[i]+tempadir);
      atom->num_bond[m]++;
      }

      }

      else atom->type[m]=1;

    }

    // Changing type of the old anchor
    if (((man = idan) >= 0) && (tempadir!=0) && (idan < atom->nlocal)){
      atom->type[man] = 1;
    }

    // Changing type of the new anchor
    if (((mannew = idnewan) >= 0) && (tempadir!=0) && (idnewan < atom->nlocal)){
      atom->type[mannew] = smctype;
    }

    // Create new bond between new hinge and anchor if not already present
    if (((mnew = idnewhi) >= 0) && (temphdir!=0) && (idnewhi < atom->nlocal)){
      atom->type[mnew]=smctype;
      
      bool create = 1;
      for (int ibond = 0; ibond < atom->num_bond[mnew]; ibond++) {
        if ((atom->bond_type[mnew][ibond] == smcbtype) || (atom->bond_type[mnew][ibond] == smcbitype)){
          create = 0;
        }
      }

      if (create) {
      // Creating new SMC bond
      if (atom->num_bond[mnew] == atom->bond_per_atom) error->one(FLERR, "New bond exceeded bonds per atom limit of {} in create_bonds", atom->bond_per_atom);
      atom->bond_type[mnew][atom->num_bond[mnew]] = smcbtype;
      atom->bond_atom[mnew][atom->num_bond[mnew]] = (anch[i]+tempadir);
      atom->num_bond[mnew]++;
      }
      }

      anch[i]+=tempadir;
      hing[i]+=temphdir;

      }
    // Barrier to check that each processor has defined correctly each smc
    MPI_Barrier(world);
  } 

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
   memory usage of hing and anch arrays
------------------------------------------------------------------------- */

double FixSMC::memory_usage()
{
  double bytes = 2 * smcnum *  sizeof(long);
  return bytes;
}


/*---------------------------------------------------------------------*/
/* Needed to write a restart file that can continue with the simulation*/
/*---------------------------------------------------------------------*/
void FixSMC::write_restart(FILE *fp)
{
    // if ((debug)) utils::logmesg(lmp, "Writing restart for fix_smc \n");

    int restart_n = 0;
    long restart_list[2+2*smcnum];

    restart_list[restart_n++] = static_cast<long>(next_reneighbor);
    restart_list[restart_n++] = static_cast<long>(update->ntimestep);
    
    // Saving SMCs positions
    for (int i = 0; i < smcnum; i++)
    {
      restart_list[restart_n++] = anch[i];
      restart_list[restart_n++] = hing[i];
    }

    if (comm->me == 0) {
    int size = restart_n * sizeof(long);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(restart_list, sizeof(long), restart_n, fp);
    }

    if ((debug)) utils::logmesg(lmp, "End of writing restart for fix_smc \n");

}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */
void FixSMC::restart(char *buf)
{
    if ((debug)) utils::logmesg(lmp, "Reading restart for fix_smc \n");

    int restart_n = 0;
    long *restart_list = (long *)buf;

    next_reneighbor = static_cast<bigint>(restart_list[restart_n++]);

    bigint ntimestep_restart = static_cast<bigint>(restart_list[restart_n++]);

    if (ntimestep_restart != update->ntimestep)
        error->all(FLERR, "Must not reset timestep when restarting fix smc");

    // Loading SMCs positions
    for (int j = 0; j < smcnum; j++)
    {
      anch[j] = static_cast<long>(restart_list[restart_n++]);     
      hing[j] = static_cast<long>(restart_list[restart_n++]);     
    }

  // Store the connected Fix ID pointer
  if (strcmp(connFixName, "nofix") == 0){
    connFix == nullptr;
  }
  else{
    connFix = modify->get_fix_by_id(connFixName);
    if (!connFix) error->all(FLERR, "Illegal ausiliary Fix");
  }

  if ((debug)) utils::logmesg(lmp, "End of reading restart for fix_smc \n");

}

/*
Returns number of smcs
*/
double FixSMC::compute_scalar(){
  return smcnum;
}
 
/*
Returns position of i smc hinge or anchor depending on the flag 
*/
double FixSMC::compute_array(int i, int flag){
  int rflag;   
  if (hdir!=0) {
    rflag = hdir/abs(hdir);
    }
    else {
    rflag = -adir/abs(adir);
    }

  rflag = (1+rflag/(abs(rflag)))/2;

  if (flag){
    if (rflag) return hing[i];
    else return anch[i];
  }
  else{
    if (rflag) return anch[i];
    else return hing[i];
  }
}