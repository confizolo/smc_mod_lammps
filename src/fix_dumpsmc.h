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
   Contributing authors: Filippo Conforto
------------------------------------------------------------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(dumpsmc,FixDUMPSMC);
// clang-format on
#else

#ifndef LMP_FIX_DUMPSMC_H
#define LMP_FIX_DUMPSMC_H

#include "fix.h"
#include <string>
#include <fstream>

namespace LAMMPS_NS
{

   class FixDUMPSMC : public Fix
   {
   public:
      FixDUMPSMC(class LAMMPS *, int, char **);
      ~FixDUMPSMC() override;
      int setmask() override;
      void init() override;
      void post_integrate() override;
      double memory_usage() override;
      void write_restart(FILE *fp) override;
      void restart(char *) override;

   private:
      long nevery;
      int nsmc;
      class Fix *connFix;
      char *connFixName;
      std::string dumpFile;
      std::ofstream dfile;
   };

} // namespace LAMMPS_NS

#endif
#endif
