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
   Contributing author: Michel Perez (U Lyon) for non-fcc lattices
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "compute_uncertainty.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "pair_agni.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

ComputeUncertainty::ComputeUncertainty(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute uncertainty command"); // arg == 3, ID, group-ID, style

  //if (strcmp(arg[3],"fcc") == 0) nnn = 12;
  ///else if (strcmp(arg[3],"bcc") == 0) nnn = 8;
  //else nnn = force->inumeric(FLERR,arg[3]);

  //if (nnn <= 0 || nnn % 2)
    //error->all(FLERR,"Illegal neighbor value for compute uncertainty command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;

  //USer defined
  epsilon = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeUncertainty::~ComputeUncertainty()
{
  //user defined
  memory->destroy(epsilon);
}

/* ---------------------------------------------------------------------- */

//end user methods
void ComputeUncertainty::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute uncertainty requires a pair style be defined");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"uncertainty") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute uncertainty");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeUncertainty::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeUncertainty::compute_uncertainty()
{
  memory->create(epsilon, list->inum, "uncertainty:epsilon");

  PairAgni *agni = (PairAgni *) force->pair_match("agni",1);

  for(int i = 0; i < list->inum; i++)
    epsilon[i] = agni->a[0] + agni->dMin[i]*(agni->a[1]) + agni->a[2]*pow(agni->dMin[i],2.0);

  vector_atom = epsilon; // vector_atom is parent array that send information to dump styles
}

void ComputeUncertainty::compute_peratom()
{
  if (atom->nlocal > nmax)
    nmax = atom->nmax;

  neighbor->build_one(list);

  compute_uncertainty();
}

/* ----------------------------------------------------------------------
   2 select routines from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   2nd routine sorts auxiliary array at same time
------------------------------------------------------------------------- */

#define SWAP(a,b)   tmp = a; a = b; b = tmp;
#define ISWAP(a,b) itmp = a; a = b; b = itmp;


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeUncertainty::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
