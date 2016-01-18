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

ComputeStyle(uncertainty,ComputeUncertainty)

#else

#ifndef COMPUTE_UNCERTAINTY_H
#define COMPUTE_UNCERTAINTY_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeUncertainty : public Compute {
 public:
  ComputeUncertainty(class LAMMPS *, int, char **);
  ~ComputeUncertainty();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

  //user methods
  static double compute_uncertainty(double, double, double, double);

 private:
  int nmax,maxneigh,nnn;
  double *distsq;
  int *nearest;
  class NeighList *list;
  double *centro;

  void select(int, int, double *);
  void select2(int, int, double *, int *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute centro/atom requires a pair style be defined

This is because the computation of the centro-symmetry values
uses a pairwise neighbor list.

W: More than one compute centro/atom

It is not efficient to use compute centro/atom more than once.

*/
