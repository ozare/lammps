/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(bop/atom,ComputeBopAtom)

#else

#ifndef LMP_COMPUTE_BOP_ATOM_H
#define LMP_COMPUTE_BOP_ATOM_H

#include "compute.h"

class SphericalHarmonics;

namespace LAMMPS_NS {

  class ComputeBopAtom : public Compute {
  public:
    ComputeBopAtom(class LAMMPS *, int, char **);
    ~ComputeBopAtom();
    void init();
    void init_list(int, class NeighList *);
    void compute_peratom();
    void set_arrays(int);
    double memory_usage();

  private:
    int nmax;

    double **invariants;
    char *id_fix;
    class FixStore *fix;
    class NeighList *list;

    bool has_cutoff;
    double cutsq;
    SphericalHarmonics *spherical_harmonics;
    int grid_quality;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute displace/atom fix ID

Self-explanatory.

*/
