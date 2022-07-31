/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(rebomos,PairREBOMoS)

#else

#ifndef LMP_PAIR_REBOMOS_H
#define LMP_PAIR_REBOMOS_H
#define PI 3.141592653

#include "pair.h"
#include "math.h"
#include "math_const.h"

namespace LAMMPS_NS {

class PairREBOMoS : public Pair {
 public:
  PairREBOMoS(class LAMMPS *);
  virtual ~PairREBOMoS();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();

 protected:
  int **pages;                     // neighbor list pages
  int *map;                        // 0 (M), 1 (S), or -1 (NULL) for each type

  int me;
  int ljflag;
  int maxlocal;                    // size of numneigh, firstneigh arrays
  int maxpage;                     // # of pages currently allocated
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom

  double cutlj;                    // user-specified LJ cutoff
  double cutljrebosq;              // cut for when to compute
                                   // REBO neighs of ghost atoms

  double **cutljsq;                // LJ cutoffs for M,S types
  double **lj1,**lj2,**lj3,**lj4;  // pre-computed LJ coeffs for M,S types
  double cut3rebo;                 // maximum distance for 3rd REBO neigh

  int *REBO_numneigh;              // # of pair neighbors for each atom
  int **REBO_firstneigh;           // ptr to 1st neighbor of each atom
  double *closestdistsq;           // closest owned atom dist to each ghost
  double *nM,*nS;                  // sum of weighting fns with REBO neighs

  double rcmin[2][2],rcmax[2][2],rcmaxsq[2][2],rcmaxp[2][2];
  double Q[2][2],alpha[2][2],A[2][2],BIJc[2][2],Beta[2][2];
  double b0[2],b1[2],b2[2],b3[2],b4[2],b5[2],b6[2];
  double bg0[2],bg1[2],bg2[2],bg3[2],bg4[2],bg5[2],bg6[2];
  double a0[2],a1[2],a2[2],a3[2];
  double rcLJmin[2][2],rcLJmax[2][2],rcLJmaxsq[2][2];
  double epsilon[2][2],sigma[2][2];

  void REBO_neigh();
  void FREBO(int, int);
  void FLJ(int, int);

  double bondorder(int, int, double *, double, double, double **, int);

  double gSpline(double, int, double *);
  double PijSpline(double, double, int, double *);

  void add_pages(int howmany = 1);
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

  inline double Sp(double Xij, double Xmin, double Xmax, double &dX) const {
    double cutoff;

    if (Xij <= Xmin){
        cutoff = 1.0;                         
        dX = 0.0;
    } else if (Xij >=Xmax){
      cutoff = 0.0;
      dX = 0.0;
    } else{
      cutoff = 0.5*(1.0 + cos(PI*(Xij - Xmin)/(Xmax - Xmin)));       
      dX = -0.5*PI*sin(PI*(Xij - Xmin)/(Xmax - Xmin))/(Xmax - Xmin);        
    }
   return cutoff;
  }; 
 };
}
  /* ----------------------------------------------------------------------
     LJ cutoff function Sp2
     return cutoff and dX = derivative
     no side effects
  ------------------------------------------------------------------------- */

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style REBOMoS requires atom IDs

This is a requirement to use the REBOMoS potential.

E: Pair style REBOMoS requires newton pair on

See the newton command.  This is a restriction to use the REBOMoS
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Neighbor list overflow, boost neigh_modify one or page

There are too many neighbors of a single atom.  Use the neigh_modify
command to increase the neighbor page size and the max number of
neighbors allowed for one atom.

E: Cannot open REBOMoS potential file %s

The specified REBOMoS potential file cannot be opened.  Check that the
path and name are correct.

*/
