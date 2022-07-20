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

PairStyle(aeam,PairAEAM)

#else

#ifndef LMP_PAIR_AEAM_H
#define LMP_PAIR_AEAM_H
#include "stdio.h"
#include "pair.h"

namespace LAMMPS_NS {

class PairAEAM : public Pair {
 public:
  double cutforcesq,cutmax;

  // per-atom arrays

  double *rho,*fp;

  // potentials as array data
	
	
  int nrhomax,nrmax;
  int nfrho,nrhor,nz2r;
  double **frho,**rhor,**z2r;
  int *nrz2r, *nrrho;
  double *drz2r, *drrho;
  int *type2frho,**type2rhor,**type2z2r;
  
  // potentials in spline form used for force computation

  double rdr1,rdrho,rdr2;
  double ***rhor_spline,***frho_spline,***z2r_spline;

  PairAEAM(class LAMMPS *);
  virtual ~PairAEAM();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
	
	


 protected:
  int nmax;                   // allocated size of per-atom arrays
  int *map;                   // which element each atom type maps to


  // potentials as file data

  struct Setfl {
    int nelements,nnonangular,nangular,nrhomax, nrmax;
	char **elements;
	char **nonangular;
	char **angular;
    double **cut, **dr;
	int *nrho;
    double *mass, *drho;
	int **nr;
    double **frho,***rhor,***z2r;
  };
  Setfl *setfl;
	


  void allocate();
  void array2spline();
  void interpolate(int, double, double *, double **);
  void grab(FILE *, int, double *);

  virtual void read_file(char *);
  virtual void file2array();
};

}

#endif
#endif






