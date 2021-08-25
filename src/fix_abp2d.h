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

#ifdef FIX_CLASS

FixStyle(abp2d,FixABP2D)

#else

#ifndef LMP_FIX_BD_H
#define LMP_FIX_BD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixABP2D : public Fix {
 public:
  FixABP2D(class LAMMPS *, int, char **);
  virtual ~FixABP2D() {}
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void reset_dt();

 protected:
  double t_target,t_period;
  double dtv,dtf;
  double gfactor;
  double spf;
  int mass_require;

  class RanMars *random;
  int seed;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
