/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(bfield,FixBfield)
// clang-format on

#else

#ifndef LMP_FIX_BFIELD_H
#define LMP_FIX_BFIELD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBfield : public Fix {

 public:
  FixBfield(class LAMMPS *, int, char **);
  ~FixBfield() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void initial_integrate(int) override;
  void post_integrate() override;
  void post_force(int) override;
  double memory_usage() override;
  double compute_scalar() override;
  double compute_vector(int) override;

 private:
  class Region *region;
  char *idregion;
  int varflag, iregion;
  char *xstr, *ystr, *zstr, *estr;
  int xvar, yvar, zvar, evar, xstyle, ystyle, zstyle, estyle;
  int maxatom;
  double dtf;
  double **v0, **fb;
  double B[3];
  double omega[3];
  double qBm2f;
  int bmuflag;
  int qflag;
  int force_flag;
  double fsum[4], fsum_all[4];
};
}    // namespace LAMMPS_NS

#endif
#endif
