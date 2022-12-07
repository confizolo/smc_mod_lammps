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

#ifdef FIX_CLASS
// clang-format off
FixStyle(smc,FixSMC);
// clang-format on
#else

#ifndef LMP_FIX_SMC_H
#define LMP_FIX_SMC_H

#include "fix.h"

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

   private:
      long *anch, *hing;
      int smctype, smcbtype, smcbitype, lpol, adir, hdir, smcnum;
      double prob, cutoff;
      double *xyzanch, *xyzhing;
      bool debug;
      class RanMars *random;
   };

} // namespace LAMMPS_NS

#endif
#endif
