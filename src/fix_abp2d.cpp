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

/* ----------------------------------------------------------------------
   Changed from fix_nve.cpp May 14th 2016 version. Guang Shi, July 2016
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Brownian Dynamics integrator. Euler Algorithm.
------------------------------------------------------------------------- */

#include <math.h>
#include "math_extra.h"
#include <stdio.h>
#include <string.h>
#include "fix_abp2d.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "utils.h"
#include "update.h"
#include "error.h"
#include "random_mars.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixABP2D::FixABP2D(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"abp2d") != 0 && narg <= 6)
    error->all(FLERR,"Illegal fix abp2d command");

  /* t_target = force->numeric(FLERR,arg[3]); // set temperature */
  /* t_period = force->numeric(FLERR,arg[4]); // same as t_period in fix_langevin_overdamp.cpp */
  /* seed = force->inumeric(FLERR,arg[5]); //seed for random number generator. integer */
  /* spf = force->numeric(FLERR,arg[6]); //seed for random number generator. integer */

  if (t_target <= 0.0) error->all(FLERR,"Fix abp2d temperature must be > 0.0");
  if (t_period <= 0.0) error->all(FLERR,"Fix abp2d period must be > 0.0");
  if (spf < 0) error->all(FLERR,"Fix abp2d self-propulsion must be > 0.0");
  if (seed <= 0) error->all(FLERR,"Illegal fix abp2d command");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixABP2D::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixABP2D::init()
{
  dtv = update->dt;  // timestep
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(2*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixABP2D::initial_integrate(int vflag)
{
  double dtfm;
  double randf;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        randf = sqrt(rmass[i]) * gfactor;
        x[i][0] += dtv * dtfm * (v[i][0]+randf*random->gaussian());
        x[i][1] += dtv * dtfm * (v[i][1]+randf*random->gaussian());
        x[i][2] = 0.0;
        v[i][2] += randf*random->gaussian();
        const double ang = v[i][2];
        v[i][0] = spf * cos(ang);
        v[i][1] = spf * sin(ang);
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        randf = sqrt(mass[type[i]]) * gfactor;
        x[i][0] += dtv * dtfm * (v[i][0]+randf*random->gaussian());
        x[i][1] += dtv * dtfm * (v[i][1]+randf*random->gaussian());
        x[i][2] = 0.0;
        v[i][2] += randf*random->gaussian();
        const double ang = v[i][2];
        v[i][0] = spf * cos(ang);
        v[i][1] = spf * sin(ang);
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixABP2D::reset_dt()
{
  dtv = update->dt;
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(2*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
}
