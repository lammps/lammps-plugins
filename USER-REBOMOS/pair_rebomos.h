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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(rebomos,PairREBOMoS);
// clang-format on
#else

#ifndef LMP_PAIR_REBOMOS_H
#define LMP_PAIR_REBOMOS_H

#include "math_const.h"
#include "pair.h"
#include <cmath>

namespace LAMMPS_NS {

class PairREBOMoS : public Pair {
 public:
  PairREBOMoS(class LAMMPS *);
  ~PairREBOMoS() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

 protected:
  double cutljrebosq;    // cut for when to compute
                         // REBO neighs of ghost atoms

  double **lj1, **lj2, **lj3, **lj4;    // pre-computed LJ coeffs for M,S types
  double cut3rebo;                      // maximum distance for 3rd REBO neigh

  int maxlocal;             // size of numneigh, firstneigh arrays
  int pgsize;               // size of neighbor page
  int oneatom;              // max # of neighbors for one atom
  MyPage<int> *ipage;       // neighbor list pages
  int *REBO_numneigh;       // # of pair neighbors for each atom
  int **REBO_firstneigh;    // ptr to 1st neighbor of each atom

  double *closestdistsq;    // closest owned atom dist to each ghost
  double *nM, *nS;          // sum of weighting fns with REBO neighs

  double rcmin[2][2], rcmax[2][2], rcmaxsq[2][2], rcmaxp[2][2];
  double Q[2][2], alpha[2][2], A[2][2], BIJc[2][2], Beta[2][2];
  double b0[2], b1[2], b2[2], b3[2], b4[2], b5[2], b6[2];
  double bg0[2], bg1[2], bg2[2], bg3[2], bg4[2], bg5[2], bg6[2];
  double a0[2], a1[2], a2[2], a3[2];
  double rcLJmin[2][2], rcLJmax[2][2], rcLJmaxsq[2][2];
  double epsilon[2][2], sigma[2][2];

  void REBO_neigh();
  void FREBO(int, int);
  void FLJ(int, int);

  double bondorder(int, int, double *, double, double, double **, int);

  double gSpline(double, int, double *);
  double PijSpline(double, double, int, double *);

  void read_file(char *);

  void allocate();

  // ----------------------------------------------------------------------
  // S'(t) and S(t) cutoff functions
  // added to header for inlining
  // ----------------------------------------------------------------------

  /* ----------------------------------------------------------------------
     cutoff function Sprime
     return cutoff and dX = derivative
     no side effects
  ------------------------------------------------------------------------- */

  inline double Sp(double Xij, double Xmin, double Xmax, double &dX) const
  {
    double cutoff;

    const double t = (Xij - Xmin) / (Xmax - Xmin);
    if (t <= 0.0) {
      cutoff = 1.0;
      dX = 0.0;
    } else if (t >= 1.0) {
      cutoff = 0.0;
      dX = 0.0;
    } else {
      cutoff = 0.5 * (1.0 + cos(t * MathConst::MY_PI));
      dX = (-0.5 * MathConst::MY_PI * sin(t * MathConst::MY_PI)) / (Xmax - Xmin);
    }
    return cutoff;
  };
};
}    // namespace LAMMPS_NS

#endif
#endif
