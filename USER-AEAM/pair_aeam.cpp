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

#include "pair_aeam.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairAEAM::PairAEAM(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  manybody_flag = 1;
  unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);

  nmax = 0;
  rho = nullptr;
  fp = nullptr;
  setfl = nullptr;
  frho = nullptr;
  rhor = nullptr;
  z2r = nullptr;

  frho_spline = nullptr;
  rhor_spline = nullptr;
  z2r_spline = nullptr;

  // set comm size needed by this Pair

  comm_forward = 1;
  comm_reverse = 1;
  one_coeff = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairAEAM::~PairAEAM()
{
  memory->destroy(rho);
  memory->destroy(fp);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete[] type2frho;
    memory->destroy(type2rhor);
    memory->destroy(type2z2r);
    delete[] nrrho;
    delete[] drrho;
    delete[] nrz2r;
    delete[] drz2r;
  }

  if (setfl) {
    for (int i = 0; i < setfl->nelements; i++) delete[] setfl->elements[i];
    delete[] setfl->elements;
    delete[] setfl->nonangular;
    delete[] setfl->angular;
    delete[] setfl->mass;
    delete[] setfl->nrho;
    delete[] setfl->drho;
    memory->destroy(setfl->nr);
    memory->destroy(setfl->dr);
    memory->destroy(setfl->cut);
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    delete setfl;
  }

  memory->destroy(frho);
  memory->destroy(rhor);
  memory->destroy(z2r);

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);
}

/* ---------------------------------------------------------------------- */

void PairAEAM::compute(int eflag, int vflag)
{
  int i, j, k, ii, jj, kk, m, m1, m2;
  int itype, jtype, ktype;
  double xtmp, ytmp, ztmp, evdwl, fpair;
  double r1, r2, r3, p, p1, p2, recip, phip, phi, rsq1, rsq2, rsq3;
  double *coeff, delr1[3], delr2[3], delr3[3], fj[3], fk[3];
  int *ilist, *numneigh, **firstneigh, *jlist, inum, jnum;
  double fij, fik, cs, delcs, ftet, delcs2;
  double dcosij, dcosik, dcosjk;
  double dfij, dfik, dfijr;
  double poteam, Feam, F2b;
  double DFij, DFik, DFjk, FFij, FFik, FFjk;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  constexpr double THIRD = 1.0 / 3.0;
  constexpr double minrho = 0.0000000000001;

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(rho, nmax, "pair:rho");
    memory->create(fp, nmax, "pair:fp");
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  int nnonangular = setfl->nnonangular;
  double cuttempij, cuttempik, CutDec;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int itypemap, jtypemap, ktypemap;

  // zero out density
  if (newton_pair) {
    for (i = 0; i < nall; i++) rho[i] = 0.0;
  } else
    for (i = 0; i < nlocal; i++) rho[i] = 0.0;

  // rho = density at each atom
  // loop over neighbors of my nonangular dependent atoms
  // Calculation of density conteribution for all atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    itypemap = itype - 1;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      jtypemap = jtype - 1;
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;

      rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];
      r1 = sqrt(rsq1);

      if (itype > setfl->nnonangular && jtype > setfl->nnonangular)
        CutDec = 1.5;
      else
        CutDec = 0;

      cuttempij = setfl->cut[itypemap][jtypemap] - CutDec;

      if (r1 > cuttempij) continue;
      rdr1 = 1 / setfl->dr[itypemap][jtypemap];
      p1 = r1 * rdr1 + 1.0;
      m1 = static_cast<int>(p1);
      m1 = MIN(m1, setfl->nr[itypemap][jtypemap] - 1);

      p1 -= m1;
      p1 = MIN(p1, 1.0);
      coeff = rhor_spline[type2rhor[itype][jtype]][m1];
      fij = ((coeff[3] * p1 + coeff[4]) * p1 + coeff[5]) * p1 + coeff[6];
      if (itype <= setfl->nnonangular) {    //atom I is noangular
        rho[i] += fij;

      } else {    //atom I is angular
        for (kk = jj + 1; kk < jnum; kk++) {

          k = jlist[kk];
          k &= NEIGHMASK;
          ktype = type[k];
          ktypemap = ktype - 1;
          delr2[0] = x[k][0] - xtmp;
          delr2[1] = x[k][1] - ytmp;
          delr2[2] = x[k][2] - ztmp;

          if (ktype > setfl->nnonangular)
            CutDec = 1.5;    //atom I and K are angular
          else
            CutDec = 0;

          cuttempik = setfl->cut[itypemap][ktypemap] - CutDec;

          rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
          r2 = sqrt(rsq2);

          if (r2 > cuttempik) continue;

          delr3[0] = x[k][0] - x[j][0];
          delr3[1] = x[k][1] - x[j][1];
          delr3[2] = x[k][2] - x[j][2];

          rsq3 = delr3[0] * delr3[0] + delr3[1] * delr3[1] + delr3[2] * delr3[2];
          r3 = sqrt(rsq3);

          rdr2 = 1 / setfl->dr[itypemap][ktypemap];

          p2 = r2 * rdr2 + 1.0;
          m2 = static_cast<int>(p2);
          m2 = MIN(m2, setfl->nr[itypemap][ktypemap] - 1);
          p2 -= m2;
          p2 = MIN(p2, 1.0);
          coeff = rhor_spline[type2rhor[itype][ktype]][m2];
          fik = ((coeff[3] * p2 + coeff[4]) * p2 + coeff[5]) * p2 + coeff[6];
          cs = (rsq1 + rsq2 - rsq3) / (2 * r1 * r2);
          delcs = cs + THIRD;
          ftet = delcs * delcs;
          rho[i] += 2 * fij * fik * ftet;
        }
      }
    }
  }

  // communicate and sum densities

  if (newton_pair) comm->reverse_comm(this);

  double ni, deli, ci, Fptmp;

  // fp = derivative of embedding energy at each atom
  // poteam = embedding energy at each atom

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    itypemap = itype - 1;
    rdrho = 1 / setfl->drho[itypemap];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    if (itype <= nnonangular) {    //atom i is Me
      ni = 1;
      deli = 0;
      ci = 0;
    } else {    //atom i is Si
      ni = 0.5;
      deli = 1;
      ci = 2;
    }

    p = pow(rho[i], ni) * rdrho + 1.0;
    m = static_cast<int>(p);
    m = MAX(1, MIN(m, setfl->nrho[itypemap] - 1));
    p -= m;
    p = MIN(p, 1.0);
    coeff = frho_spline[type2frho[itype]][m];
    fp[i] = (coeff[0] * p + coeff[1]) * p + coeff[2];    //F'(sum(fij))
    if (eflag) {
      //F(rhoi) for all atoms F(sum(fij)) for Me F(2sum(fij.fik.ftet)^0.5 for Si
      poteam = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
      if (eflag_global) eng_vdwl += poteam;
      if (eflag_atom) {
        if (itype <= nnonangular) {
          eatom[i] += poteam;
        } else {
          eatom[i] += THIRD * poteam;
        }
      }
    }
  }

  // communicate derivative of embedding function

  comm->forward_comm(this);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    itypemap = itype - 1;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    if (itype <= nnonangular) {
      ni = 1;
      deli = 0;
      ci = 0;
    } else {
      ni = 0.5;
      deli = 1;
      ci = 2;
    }

    if (rho[i] > minrho)
      Fptmp = ni *
          pow(rho[i],
              (ni -
               1));    //for Me Fptmp=1  for Si Fptmp=0.5(2sum(fij.fik.ftet)))^(-0.5)=0.5/rho(i)
    else
      Fptmp = 0;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      jtypemap = jtype - 1;

      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;

      rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];
      r1 = sqrt(rsq1);

      if (r1 > setfl->cut[itypemap][jtypemap]) continue;

      rdr1 = 1 / setfl->dr[itypemap][jtypemap];
      p1 = r1 * rdr1 + 1.0;
      m1 = static_cast<int>(p1);
      m1 = MIN(m1, setfl->nr[itypemap][jtypemap] - 1);

      p1 -= m1;
      p1 = MIN(p1, 1.0);

      //fij density contribution of atom J on atom I
      //dfij derivative of density contribution of atom J on atom I
      //phi pair potential
      //phip derivative of pair potential
      coeff = rhor_spline[type2rhor[itype][jtype]][m1];
      fij = ((coeff[3] * p1 + coeff[4]) * p1 + coeff[5]) * p1 + coeff[6];
      dfij = (coeff[0] * p1 + coeff[1]) * p1 + coeff[2];
      coeff = z2r_spline[type2z2r[itype][jtype]][m1];
      phip = (coeff[0] * p1 + coeff[1]) * p1 + coeff[2];
      phi = ((coeff[3] * p1 + coeff[4]) * p1 + coeff[5]) * p1 + coeff[6];

      recip = 1 / r1;
      dfijr = dfij * recip;
      Feam = -(1 - deli) * Fptmp * fp[i] * dfijr;    // Me: -F'(sum(fij))*1*dfij/r1  Si: 0

      F2b = -phip * recip;
      fpair = Feam + 0.5 * F2b;

      f[i][0] -= delr1[0] * fpair;
      f[i][1] -= delr1[1] * fpair;
      f[i][2] -= delr1[2] * fpair;

      f[j][0] += delr1[0] * fpair;
      f[j][1] += delr1[1] * fpair;
      f[j][2] += delr1[2] * fpair;

      if (eflag) {
        evdwl = 0.5 * phi;
        if (eflag_global) eng_vdwl += evdwl;
        if (eflag_atom) eatom[i] += evdwl;
      }

      if (evflag)
        ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delr1[0], delr1[1], delr1[2]);

      if (itype <= setfl->nnonangular)
        continue;    //atom I is Au

      else {    //atom I is Si
        for (kk = jj + 1; kk < jnum; kk++) {

          k = jlist[kk];
          k &= NEIGHMASK;
          ktype = type[k];
          ktypemap = ktype - 1;
          delr2[0] = x[k][0] - xtmp;
          delr2[1] = x[k][1] - ytmp;
          delr2[2] = x[k][2] - ztmp;

          if (ktype > setfl->nnonangular)
            CutDec = 1.5;
          else
            CutDec = 0;

          cuttempik = setfl->cut[itypemap][ktypemap] - CutDec;

          rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
          r2 = sqrt(rsq2);

          if (r2 > cuttempik) continue;

          delr3[0] = x[k][0] - x[j][0];
          delr3[1] = x[k][1] - x[j][1];
          delr3[2] = x[k][2] - x[j][2];

          rsq3 = delr3[0] * delr3[0] + delr3[1] * delr3[1] + delr3[2] * delr3[2];
          r3 = sqrt(rsq3);

          rdr2 = 1 / setfl->dr[itypemap][ktypemap];

          p2 = r2 * rdr2 + 1.0;
          m2 = static_cast<int>(p2);
          m2 = MIN(m2, setfl->nr[itypemap][ktypemap] - 1);
          p2 -= m2;
          p2 = MIN(p2, 1.0);
          coeff = rhor_spline[type2rhor[itype][ktype]][m2];
          fik = ((coeff[3] * p2 + coeff[4]) * p2 + coeff[5]) * p2 + coeff[6];    //fik for all atoms
          dfik = (coeff[0] * p2 + coeff[1]) * p2 + coeff[2];    //dfik for all atoms

          cs = (rsq1 + rsq2 - rsq3) / (2 * r1 * r2);
          dcosij = 1 / r2 - cs / r1;
          dcosik = 1 / r1 - cs / r2;
          dcosjk = -r3 / (r1 * r2);
          delcs = (cs + THIRD);
          ftet = delcs * delcs;
          delcs2 = 2 * delcs;

          DFij = ci * (fik * dfij * ftet + fij * fik * delcs2 * dcosij);
          DFik = ci * (fij * dfik * ftet + fij * fik * delcs2 * dcosik);
          DFjk = ci * fij * fik * delcs2 * dcosjk;

          FFij = -Fptmp * fp[i] * DFij / r1;
          FFik = -Fptmp * fp[i] * DFik / r2;
          FFjk = -Fptmp * fp[i] * DFjk / r3;

          fj[0] = delr1[0] * FFij - delr3[0] * FFjk;
          fj[1] = delr1[1] * FFij - delr3[1] * FFjk;
          fj[2] = delr1[2] * FFij - delr3[2] * FFjk;

          fk[0] = delr2[0] * FFik + delr3[0] * FFjk;
          fk[1] = delr2[1] * FFik + delr3[1] * FFjk;
          fk[2] = delr2[2] * FFik + delr3[2] * FFjk;

          f[i][0] -= fj[0] + fk[0];
          f[i][1] -= fj[1] + fk[1];
          f[i][2] -= fj[2] + fk[2];
          f[j][0] += fj[0];
          f[j][1] += fj[1];
          f[j][2] += fj[2];
          f[k][0] += fk[0];
          f[k][1] += fk[1];
          f[k][2] += fk[2];

          if (evflag) ev_tally3(i, j, k, 0.0, 0.0, fj, fk, delr1, delr2);
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairAEAM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  delete[] map;
  map = new int[n + 1];
  for (int i = 1; i <= n; i++) map[i] = -1;

  type2frho = new int[n + 1];
  nrrho = new int[n * n];
  drrho = new double[n * n];
  nrz2r = new int[n * n];
  drz2r = new double[n * n];
  memory->create(type2rhor, n + 1, n + 1, "pair:type2rhor");
  memory->create(type2z2r, n + 1, n + 1, "pair:type2z2r");
}

/* ----------------------------------------------------------------------
 global settings 
 ------------------------------------------------------------------------- */

void PairAEAM::settings(int narg, char ** /*arg*/)
{

  if (narg != 0) error->all(FLERR, "Illegal pair_style command");
}
/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO funcfl file
------------------------------------------------------------------------- */

void PairAEAM::coeff(int narg, char **arg)
{
  int i, j;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes) error->all(FLERR, "Incorrect args for pair coefficients");

  if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
    error->all(FLERR, "Incorrect args for pair coefficients");

  // read AEAM setfl file

  if (setfl) {
    for (i = 0; i < setfl->nelements; i++) delete[] setfl->elements[i];
    delete[] setfl->elements;
    for (i = 0; i < setfl->nnonangular; i++) delete[] setfl->nonangular[i];
    delete[] setfl->nonangular;
    for (i = 0; i < setfl->nangular; i++) delete[] setfl->angular[i];
    delete[] setfl->angular;
    delete[] setfl->mass;
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    delete setfl;
  }

  setfl = new Setfl();
  read_file(arg[2]);

  // read args that map atom types to elements in potential file

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i], "NULL") == 0) {
      map[i - 2] = -1;
      continue;
    }
    for (j = 0; j < setfl->nelements; j++)
      if (strcmp(arg[i], setfl->elements[j]) == 0) break;
    if (j < setfl->nelements)
      map[i - 2] = j;
    else
      error->all(FLERR, "No matching element in AEAM potential file");
  }

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i], setfl->elements[i - 3]) != 0) {
      error->all(FLERR, "no matching atom order of input file and potential file");
    }
  }

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++) setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j

  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR, i, setfl->mass[map[i]]);
        count++;
      }
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAEAM::init_style()
{
  file2array();
  array2spline();

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAEAM::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
  cutmax = setfl->cut[i - 1][j - 1];
  cutforcesq = cutmax * cutmax;
  return cutmax;
}

/* ----------------------------------------------------------------------
 read a multi-element DYNAMO setfl file
 ------------------------------------------------------------------------- */

void PairAEAM::read_file(char *filename)
{
  Setfl *file = setfl;

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];

  if (me == 0) {
    fptr = fopen(filename, "r");
    if (fptr == nullptr) {
      char str[128];
      sprintf(str, "Cannot open AEAM potential file %s", filename);
      error->one(FLERR, str);
    }
  }

  /*-----------------------------------------------------------------------------
	 omit the first ten+1 lines(header of the aeam), read 
	 -----------------------------------------------------------------------------*/
  int n, i, j;
  int nheader1 = 12;

  if (me == 0) {
    for (i = 0; i < nheader1; i++) fgets(line, MAXLINE, fptr);
    n = strlen(line) + 1;
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, world);
  MPI_Bcast(line, n, MPI_CHAR, 0, world);

  sscanf(line, "%d %d %d", &file->nelements, &file->nnonangular, &file->nangular);
  int nwords = utils::count_words(line);
  if (nwords != file->nelements + 3 || file->nelements != file->nnonangular + file->nangular)
    error->all(FLERR, "Incorrect element names in AEAM potential file");

  char **words = new char *[file->nelements + 3];
  nwords = 0;
  char *first = strtok(line, " \t\n\r\f");
  while (words[nwords++] = strtok(NULL, " \t\n\r\f")) continue;

  file->elements = new char *[file->nelements];
  file->nonangular = new char *[file->nnonangular];
  file->angular = new char *[file->nangular];
  for (i = 0; i < file->nelements; i++) file->elements[i] = utils::strdup(words[i + 2]);

  delete[] words;

  file->mass = new double[file->nelements];
  file->nrho = new int[file->nelements];
  file->drho = new double[file->nelements];

  file->nrmax = 0;
  file->nrhomax = 0;

  for (i = 0; i < file->nelements; i++) {
    if (me == 0) {
      fgets(line, MAXLINE, fptr);
      sscanf(line, "%d %lg %lg", &file->nrho[i], &file->drho[i], &file->mass[i]);
    }

    MPI_Bcast(&file->nrho[i], 1, MPI_INT, 0, world);
    MPI_Bcast(&file->drho[i], 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&file->mass[i], 1, MPI_DOUBLE, 0, world);

    if (file->nrhomax < file->nrho[i]) file->nrhomax = file->nrho[i];
  }

  memory->create(file->nr, file->nelements, file->nelements, "pair:nr");
  memory->create(file->dr, file->nelements, file->nelements, "pair:dr");
  memory->create(file->cut, file->nelements, file->nelements, "pair:cut");

  for (i = 0; i < file->nelements; i++) {
    for (j = 0; j < file->nelements; j++) {
      if (me == 0) {
        fgets(line, MAXLINE, fptr);
        sscanf(line, "%d %lg %lg", &file->nr[i][j], &file->dr[i][j], &file->cut[i][j]);
      }

      MPI_Bcast(&file->nr[i][j], 1, MPI_INT, 0, world);
      MPI_Bcast(&file->dr[i][j], 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&file->cut[i][j], 1, MPI_DOUBLE, 0, world);
      if (file->nrmax < file->nr[i][j]) file->nrmax = file->nr[i][j];
    }
  }

  memory->create(file->frho, file->nelements, file->nrhomax + 1, "pair:frho");
  memory->create(file->rhor, file->nelements, file->nelements, file->nrmax + 1, "pair:rhor");
  memory->create(file->z2r, file->nelements, file->nelements, file->nrmax + 1, "pair:z2r");

  for (i = 0; i < file->nelements; i++) {
    if (me == 0) grab(fptr, file->nrho[i], &file->frho[i][1]);
    MPI_Bcast(&file->frho[i][1], file->nrho[i], MPI_DOUBLE, 0, world);
  }

  for (i = 0; i < file->nelements; i++)
    for (j = 0; j < file->nelements; j++) {
      if (me == 0) grab(fptr, file->nr[i][j], &file->rhor[i][j][1]);
      MPI_Bcast(&file->rhor[i][j][1], file->nr[i][j], MPI_DOUBLE, 0, world);
    }

  for (i = 0; i < file->nelements; i++)
    for (j = 0; j <= i; j++) {
      if (me == 0) grab(fptr, file->nr[i][j], &file->z2r[i][j][1]);
      MPI_Bcast(&file->z2r[i][j][1], file->nr[i][j], MPI_DOUBLE, 0, world);
    }

  if (me == 0) fclose(fptr);
}

/* ----------------------------------------------------------------------
 copy read-in setfl potential to standard array format
 ------------------------------------------------------------------------- */

void PairAEAM::file2array()
{
  int i, j, m, n;
  int ntypes = atom->ntypes;

  nrhomax = setfl->nrhomax;
  nrmax = setfl->nrmax;

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of setfl elements + 1 for zero array

  nfrho = setfl->nelements + 1;
  memory->destroy(frho);
  memory->create(frho, nfrho, nrhomax + 1, "pair:frho");

  // copy each element's frho to global frho

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= setfl->nrho[i]; m++) frho[i][m] = setfl->frho[i][m];

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrhomax; m++) frho[nfrho - 1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to element (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++) {
    if (map[i] >= 0)
      type2frho[i] = map[i];
    else
      type2frho[i] = nfrho - 1;
  }

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of setfl elements^2

  nrhor = setfl->nelements * setfl->nelements;
  memory->destroy(rhor);
  memory->create(rhor, nrhor, nrmax + 1, "pair:rhor");

  // copy each element's rhor to global rhor
  n = 0;

  for (i = 0; i < setfl->nelements; i++)
    for (j = 0; j < setfl->nelements; j++) {
      for (m = 1; m <= setfl->nr[i][j]; m++) rhor[n][m] = setfl->rhor[i][j][m];
      nrrho[n] = setfl->nr[i][j];
      drrho[n] = setfl->dr[i][j];
      n++;
    }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used
  int mapp = 0;
  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++) {
      type2rhor[i][j] = mapp;
      mapp += 1;
    }

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of setfl elements

  nz2r = setfl->nelements * (setfl->nelements + 1) / 2;
  memory->destroy(z2r);
  memory->create(z2r, nz2r, nrmax + 1, "pair:z2r");

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < setfl->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= setfl->nr[i][j]; m++) z2r[n][m] = setfl->z2r[i][j][m];
      nrz2r[n] = setfl->nr[i][j];
      drz2r[n] = setfl->dr[i][j];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow, icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairAEAM::array2spline()
{
  int nrhoint, nrint;
  double drhoint, drint;

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);

  memory->create(frho_spline, nfrho, nrhomax + 1, 7, "pair:frho");
  memory->create(rhor_spline, nrhor, nrmax + 1, 7, "pair:rhor");
  memory->create(z2r_spline, nz2r, nrmax + 1, 7, "pair:z2r");

  for (int i = 0; i < nfrho; i++) {
    if (i < (nfrho - 1)) {
      nrhoint = setfl->nrho[i];
      drhoint = setfl->drho[i];
    } else {    //for the last part of Frho which is all zero use for hybrid
      nrhoint = setfl->nrho[0];
      drhoint = setfl->drho[0];
    }
    interpolate(nrhoint, drhoint, frho[i], frho_spline[i]);
  }

  for (int i = 0; i < nrhor; i++) {
    nrint = nrrho[i];
    drint = drrho[i];
    interpolate(nrint, drint, rhor[i], rhor_spline[i]);
  }

  for (int i = 0; i < nz2r; i++) {
    nrint = nrz2r[i];
    drint = drz2r[i];
    interpolate(nrint, drint, z2r[i], z2r_spline[i]);
  }
}

/* ---------------------------------------------------------------------- */

void PairAEAM::interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6] - spline[1][6]);
  spline[n - 1][5] = 0.5 * (spline[n][6] - spline[n - 2][6]);
  spline[n][5] = spline[n][6] - spline[n - 1][6];

  for (int m = 3; m <= n - 2; m++)
    spline[m][5] =
        ((spline[m - 2][6] - spline[m + 2][6]) + 8.0 * (spline[m + 1][6] - spline[m - 1][6])) /
        12.0;

  for (int m = 1; m <= n - 1; m++) {
    spline[m][4] = 3.0 * (spline[m + 1][6] - spline[m][6]) - 2.0 * spline[m][5] - spline[m + 1][5];
    spline[m][3] = spline[m][5] + spline[m + 1][5] - 2.0 * (spline[m + 1][6] - spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5] / delta;
    spline[m][1] = 2.0 * spline[m][4] / delta;
    spline[m][0] = 3.0 * spline[m][3] / delta;
  }
}

/* ----------------------------------------------------------------------
   grab n values from file fp and put them in list
   values can be several to a line
   only called by proc 0
------------------------------------------------------------------------- */

void PairAEAM::grab(FILE *fptr, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line, MAXLINE, fptr);
    ptr = strtok(line, " \t\n\r\f");
    list[i++] = atof(ptr);
    while (ptr = strtok(NULL, " \t\n\r\f")) list[i++] = atof(ptr);
  }
}

/* ---------------------------------------------------------------------- */

int PairAEAM::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    buf[m++] = fp[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairAEAM::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) fp[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int PairAEAM::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = rho[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void PairAEAM::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double PairAEAM::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom * 6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}
