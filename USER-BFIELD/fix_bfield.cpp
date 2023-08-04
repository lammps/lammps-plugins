/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_bfield.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NONE, CONSTANT, EQUAL, ATOM };

using MathConst::MY_PI;
static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

FixBfield::FixBfield(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), region(nullptr), idregion(nullptr), xstr(nullptr), ystr(nullptr),
    zstr(nullptr), estr(nullptr), v0(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix bfield", error);

  dynamic_group_allow = 0;
  vector_flag = 1;
  scalar_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  extscalar = 1;

  energy_global_flag = 1;
  virial_global_flag = virial_peratom_flag = 0;

  xstyle = ystyle = zstyle = estyle = NONE;

  if (utils::strmatch(arg[3], "^v_")) {
    xstr = utils::strdup(arg[3] + 2);
  } else {
    B[0] = utils::numeric(FLERR, arg[3], false, lmp);
    xstyle = CONSTANT;
  }

  if (utils::strmatch(arg[4], "^v_")) {
    ystr = utils::strdup(arg[4] + 2);
  } else {
    B[1] = utils::numeric(FLERR, arg[4], false, lmp);
    ystyle = CONSTANT;
  }

  if (utils::strmatch(arg[5], "^v_")) {
    zstr = utils::strdup(arg[5] + 2);
  } else {
    B[2] = utils::numeric(FLERR, arg[5], false, lmp);
    zstyle = CONSTANT;
  }

  // optional args

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + " region", error);
      region = domain->get_region_by_id(arg[iarg + 1]);
      if (!region) error->all(FLERR, "Region {} for fix {} does not exist", arg[iarg + 1], style);
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown keyword for fix {} command: {}", style, arg[iarg]);
    }
  }

  force_flag = 0;
  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;

  maxatom = atom->nmax;
  memory->create(v0, maxatom, 3, "v0:bfield");
}

/* ---------------------------------------------------------------------- */

FixBfield::~FixBfield()
{
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] estr;
  delete[] idregion;
  memory->destroy(v0);
}

/* ---------------------------------------------------------------------- */

int FixBfield::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBfield::init()
{
  qflag = bmuflag = 0;
  // charges needed for lorentz force
  if (atom->q_flag) qflag = 1;
  //  if (atom->bmu_flag && atom->torque_flag) bmuflag = 1;
  //  if (!qflag && !bmuflag)
  //    error->all(FLERR,"Fix bfield requires atom attribute q or bmu");

  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR, "Variable {} for x-field in fix {} does not exist", xstr, style);
    if (input->variable->equalstyle(xvar))
      xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar))
      xstyle = ATOM;
    else
      error->all(FLERR, "Variable {} for x-field in fix {} is invalid style", xstr, style);
  }

  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR, "Variable {} for y-field in fix {} does not exist", ystr, style);
    if (input->variable->equalstyle(yvar))
      ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar))
      ystyle = ATOM;
    else
      error->all(FLERR, "Variable {} for y-field in fix {} is invalid style", ystr, style);
  }

  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR, "Variable {} for z-field in fix {} does not exist", zstr, style);
    if (input->variable->equalstyle(zvar))
      zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar))
      zstyle = ATOM;
    else
      error->all(FLERR, "Variable {} for z-field in fix {} is invalid style", zstr, style);
  }

  // set units
  // B-field input in Tesla (T) for all unit sets except LJ
  // LJ: B in tau*q/m
  if (strcmp(update->unit_style, "lj") == 0)
    qBm2f = 1;
  else if (strcmp(update->unit_style, "real") == 0)
    qBm2f = 1.60217646e-19 / 1.66054e-27 /
        1e15;    // coulomb per electron charge / kg per amu / fs per s
  else if (strcmp(update->unit_style, "metal") == 0)
    qBm2f = 1.60217646e-19 / 1.66054e-27 /
        1e12;    // coulomb per electron charge / kg per amu / ps per s
  else if (strcmp(update->unit_style, "si") == 0)
    qBm2f = 1 / 1 / 1;    // coulomb per coulomb / kg per kg / s per s
  else if (strcmp(update->unit_style, "cgs") == 0)
    qBm2f = 3.356e-10 / 1.66054e-24 / 1;    // coulomb per statcoulomb / g per amu / s per s
  else if (strcmp(update->unit_style, "electron") == 0)
    qBm2f = 1.60217646e-19 / 1.66054e-27 /
        1e15;    // coulomb per electron charge / kg per amu / fs per s
  else if (strcmp(update->unit_style, "micro") == 0)
    qBm2f = 1e-12 / 1.66054e-12 / 1e6;    // coulomb per picocoulomb / kg per picogram / ms per s
  else if (strcmp(update->unit_style, "nano") == 0)
    qBm2f = 1.60217646e-19 / 1.66054e-6 /
        1e9;    // coulomb per electron charge / kg per attogram / ns per s
  else
    error->all(FLERR, "Illegal units {} in fix {}", update->unit_style, style);

  dtf = 0.5 * update->dt * force->ftm2v;

  // order of pre_integrate fixes important
  // fix bfield needs current v(t)
  // so fix bfield should occur before other pre_integrate methods
  int bfield = 0;

  for (const auto &ifix : modify->get_fix_list()) {
    if (utils::strmatch(ifix->style, "^bfield") == 0) bfield = 1;
    if (ifix->time_integrate) {
      if (!bfield)
        error->all(FLERR, "fix bfield must be defined before NVE style time integrator");
      if (!utils::strmatch(ifix->style, "^nve"))
        error->all(FLERR, "fix bfield requires an NVE style integrator");
    }
  }

  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix {} does not exist", idregion, style);
  }

  // variable flags
  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    error->all(FLERR, "Fix bfield cannot use atom-style variables");
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else
    varflag = CONSTANT;

  if (qflag) {

    // x,v update is for limit of weak B-field
    // will not work well for strong fields
    double max_omega = 2 * MY_PI * SMALL / update->dt;
    double *q = atom->q;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    double dtfm, c1;

    if (rmass) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {

          dtfm = dtf / rmass[i];
          c1 = qBm2f * q[i] / rmass[i];
          omega[0] = c1 * B[0];
          omega[1] = c1 * B[1];
          omega[2] = c1 * B[2];

          if (fabs(omega[0]) > max_omega || omega[1] > max_omega ||
              omega[2] > max_omega)    // Spreiter Eq. 1
            error->warning(FLERR, "fix bfield does not support strong magnetic fields");
        }
    } else {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {

          dtfm = dtf / mass[type[i]];
          c1 = qBm2f * q[i] / mass[type[i]];
          omega[0] = c1 * B[0];
          omega[1] = c1 * B[1];
          omega[2] = c1 * B[2];

          if (fabs(omega[0]) > max_omega || omega[1] > max_omega ||
              omega[2] > max_omega)    // Spreiter Eq. 1
            error->warning(FLERR, "fix bfield does not support strong magnetic fields");
        }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBfield::setup(int vflag)
{
  if (strstr(update->integrate_style, "verlet")) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBfield::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
  save current velocites
  runs before update to 1/2 step velocity
------------------------------------------------------------------------- */

void FixBfield::initial_integrate(int vflag)
{
  if (!qflag) return;

  int nlocal = atom->nlocal;
  double **v = atom->v;
  double **f = atom->f;

  // reallocate v0 array if necessary
  if (nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(v0);
    memory->create(v0, maxatom, 3, "bfield:v0");
  }

  for (int i = 0; i < nlocal; i++) {
    v0[i][0] = v[i][0];
    v0[i][1] = v[i][1];
    v0[i][2] = v[i][2];
  }
}

/* ----------------------------------------------------------------------
   update v,x for b-field
   v,x are analytically integrated
   Spreiter and Walter, J. Comp. Phys. 102-119 (1999)
------------------------------------------------------------------------- */

void FixBfield::post_integrate()
{
  if (!qflag) return;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *q = atom->q;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  imageint *image = atom->image;
  double vx, vy, vz, c1;
  double fx, fy, fz;
  double dtfm;
  double dtv = update->dt;
  double dtv_omega0, dtv_omega1, dtv_omega2;
  double half_dtfm, half_dtv_omega0, half_dtv_omega1, half_dtv_omega2;
  double unwrap[3];

  // update region if necessary

  if (region) {
    region = domain->get_region_by_id(idregion);
    region->prematch();
  }

  // fsum[0] = "potential energy" for added force
  // fsum[123] = extra force added to atoms

  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
  force_flag = 0;

  // constant bfield

  if (rmass) {

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;

        dtfm = dtf / rmass[i];
        c1 = qBm2f * q[i] / rmass[i];

        omega[0] = c1 * B[0];
        omega[1] = c1 * B[1];
        omega[2] = c1 * B[2];
        vx = v0[i][0];
        vy = v0[i][1];
        vz = v0[i][2];

        //velocity update
        dtv_omega0 = dtv * omega[0];
        dtv_omega1 = dtv * omega[1];
        dtv_omega2 = dtv * omega[2];
        half_dtfm = 0.5 * dtfm;
        half_dtv_omega0 = 0.5 * dtv_omega0;
        half_dtv_omega1 = 0.5 * dtv_omega1;
        half_dtv_omega2 = 0.5 * dtv_omega2;

        // B0
        v[i][1] += dtv_omega0 * (vz + half_dtfm * f[i][1] - half_dtv_omega0 * vy);
        v[i][2] += -dtv_omega0 * (vy + half_dtfm * f[i][0] + half_dtv_omega0 * vz);
        // B1
        v[i][0] += -dtv_omega1 * (vz + half_dtfm * f[i][0] + half_dtv_omega1 * vx);
        v[i][2] += dtv_omega1 * (vx + half_dtfm * f[i][1] - half_dtv_omega1 * vz);
        // B2
        v[i][0] += dtv_omega2 * (vy + half_dtfm * f[i][1] - half_dtv_omega2 * vx);
        v[i][1] += -dtv_omega2 * (vx + half_dtfm * f[i][0] + half_dtv_omega2 * vy);

        // position update
        // B0
        x[i][1] += (dtv * half_dtv_omega0 * vz);
        x[i][2] += (-dtv * half_dtv_omega0 * vy);
        // B1
        x[i][0] += (-dtv * half_dtv_omega1 * vz);
        x[i][2] += (dtv * half_dtv_omega1 * vx);
        // B2
        x[i][0] += (dtv * half_dtv_omega2 * vy);
        x[i][1] += (-dtv * half_dtv_omega2 * vx);

        // estimate total magnetic force for analysis
        // not used for dynamics
        fx = q[i] * (vy * B[2] - vz * B[1]);
        fy = q[i] * (vz * B[0] - vx * B[2]);
        fz = q[i] * (vx * B[1] - vy * B[0]);
        domain->unmap(x[i], image[i], unwrap);
        fsum[0] -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
        fsum[1] += fx;
        fsum[2] += fy;
        fsum[3] += fz;
      }

  } else {

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
        dtfm = dtf / mass[type[i]];
        c1 = qBm2f * q[i] / mass[type[i]];
        //update for B-field
        omega[0] = c1 * B[0];
        omega[1] = c1 * B[1];
        omega[2] = c1 * B[2];

        vx = v0[i][0];
        vy = v0[i][1];
        vz = v0[i][2];

        //velocity update
        dtv_omega0 = dtv * omega[0];
        dtv_omega1 = dtv * omega[1];
        dtv_omega2 = dtv * omega[2];
        half_dtfm = 0.5 * dtfm;
        half_dtv_omega0 = 0.5 * dtv_omega0;
        half_dtv_omega1 = 0.5 * dtv_omega1;
        half_dtv_omega2 = 0.5 * dtv_omega2;

        // B0
        v[i][1] += dtv_omega0 * (vz + half_dtfm * f[i][1] - half_dtv_omega0 * vy);
        v[i][2] += -dtv_omega0 * (vy + half_dtfm * f[i][0] + half_dtv_omega0 * vz);
        // B1
        v[i][0] += -dtv_omega1 * (vz + half_dtfm * f[i][0] + half_dtv_omega1 * vx);
        v[i][2] += dtv_omega1 * (vx + half_dtfm * f[i][1] - half_dtv_omega1 * vz);
        // B2
        v[i][0] += dtv_omega2 * (vy + half_dtfm * f[i][1] - half_dtv_omega2 * vx);
        v[i][1] += -dtv_omega2 * (vx + half_dtfm * f[i][0] + half_dtv_omega2 * vy);

        // position update
        // B0
        x[i][1] += (dtv * half_dtv_omega0 * vz);
        x[i][2] += (-dtv * half_dtv_omega0 * vy);
        // B1
        x[i][0] += (-dtv * half_dtv_omega1 * vz);
        x[i][2] += (dtv * half_dtv_omega1 * vx);
        // B2
        x[i][0] += (dtv * half_dtv_omega2 * vy);
        x[i][1] += (-dtv * half_dtv_omega2 * vx);

        // compute total magnetic force for analysis
        // not used for dynamics
        fx = q[i] * (vy * B[2] - vz * B[1]);
        fy = q[i] * (vz * B[0] - vx * B[2]);
        fz = q[i] * (vx * B[1] - vy * B[0]);
        domain->unmap(x[i], image[i], unwrap);
        fsum[0] -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
        fsum[1] += fx;
        fsum[2] += fy;
        fsum[3] += fz;
      }
    /*
      C A/ps T = C Ang/ps N /(Amp m) = 100 N wychodzi fx, czyli podzielic na 100 zeby dostac w N
      eV/Ang jednostka sily w metal, = 1.6e-19 J/Ang = 1.6e-19 N m/1e-10 m= 1.6e-9 N, czyli 1N = (1eV/Ang)/1.6e-9
      czyli 0.01 u = (1ev/Ang)/1.6e-9, 1eV/ang = 1.6e-11 u        
    */
  }
}
/* ----------------------------------------------------------------------
    no force, torque = bmu cross B
------------------------------------------------------------------------- */

void FixBfield::post_force(int vflag)
{
  if (!bmuflag) return;

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  //  double **bmu = atom->bmu;
  double **t = atom->torque;
  double tx, ty, tz;

  //update region if necessary

  if (region) {
    region = domain->get_region_by_id(idregion);
    region->prematch();
  }

  // fsum[0] = "potential energy" for added force
  force_flag = 0;

  if (varflag == EQUAL) {
    modify->clearstep_compute();
    if (xstyle == EQUAL) B[0] = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) B[1] = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) B[2] = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);
  }

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixBfield::memory_usage()
{
  double bytes;
  bytes = atom->nmax * 3 * sizeof(double);    //v0
  return bytes;
}

/* ----------------------------------------------------------------------
   return energy added by fix
------------------------------------------------------------------------- */

double FixBfield::compute_scalar(void)
{
  if (force_flag == 0) {
    MPI_Allreduce(fsum, fsum_all, 4, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return fsum_all[0];
}

/* ----------------------------------------------------------------------
   return total extra force due to fix
------------------------------------------------------------------------- */

double FixBfield::compute_vector(int n)
{
  if (force_flag == 0) {
    MPI_Allreduce(fsum, fsum_all, 4, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return fsum_all[n + 1];
}
