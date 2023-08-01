/* -*- c++ -*- ----------------------------------------------------------
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
#ifdef FIX_CLASS
// clang-format off
FixStyle(smc,FixSMC);
// clang-format on
#else

#ifndef LMP_FIX_SMC_H
#define LMP_FIX_SMC_H

#include "fix.h"
#include <string>
#include <fstream>

namespace LAMMPS_NS
{

   class FixSMC : public Fix
   {
   public:
      FixSMC(class LAMMPS *, int, char **);
      ~FixSMC() override;
      int setmask() override;
      void init() override;
      void post_integrate() override;
      double memory_usage() override;
      void write_restart(FILE *fp) override;
      void restart(char *) override;
      double compute_array(int, int) override;
      double compute_scalar() override;

      bool check_avl(long);
      void compile_avl_list();
      void load_smc(long);
      void place_smc(long, long, bool);
      void remove_smc(long, long);
      std::array<double, 3> compute_xyz(long);
      void debug_pre(int);
      void debug_post(int, bool, bool, bool, bool, bool, double, double);

   private:
      long *anch, *hing, *av_list;
      long num_avl;
      int seed, smctype, smcbtype, smcbitype, lpol, maxadir, maxhdir, adir, hdir, smcnum, initmode, ring, nblockt;
      int *blockt;
      double prob, cutoff, tancoff, kon, koff;
      bool debug;
      std::ofstream debugfile;
      class RanPark *random_equal;
      std::string fixFname;
   };

} // namespace LAMMPS_NS

#endif
#endif
