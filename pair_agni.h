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

#ifdef PAIR_CLASS

PairStyle(agni,PairAgni)

#else

#ifndef LMP_PAIR_AGNI_H
#define LMP_PAIR_AGNI_H

#include "pair.h"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

using namespace std;

namespace LAMMPS_NS {

class PairAgni : public Pair {
 public:
  PairAgni(class LAMMPS *);
  virtual ~PairAgni();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void readUserFile();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);


 protected:
  double cut_global;
  double **cut; //DON'T TOUCH THIS, WILL MESS UP ATOMS!!!!!!!!
  //double **epsilon,**sigma;//NO CHANGE
  //double **lj1,**lj2,**lj3,**lj4,**offset;//NO CHANGE
  double *cut_respa;

  //user variables
  string inputFile;
  int nTrain;
  vector<double> Rc,sigma1,lambda,b;
  vector<string> elementList;
  vector< vector<double> > eta,yU,alpha;
  vector< vector < vector<double> > > xU;

  bool start;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

*/