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

  nnn = 0;

  nmax = 0;
  centro = NULL;
  maxneigh = 0;
  distsq = NULL;
  nearest = NULL;

  //USer defined
  epsilon = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeUncertainty::~ComputeUncertainty()
{
  memory->destroy(centro);
  memory->destroy(distsq);
  memory->destroy(nearest);

  //user defined
  memory->destroy(epsilon);
}

/* ---------------------------------------------------------------------- */



//user methods
void ComputeUncertainty::compute_uncertainty()
{
  memory->create(epsilon, list->inum, "uncertainty:epsilon");

  PairAgni *agni = (PairAgni *) force->pair_match("agni",1);

  for(int i = 0; i < list->inum; i++)
    epsilon[i] = agni->a[0] + agni->dMin[i]*(agni->a[1]) + agni->a[2]*pow(agni->dMin[i],2.0);

  vector_atom = epsilon; // vector_atom is parent array that send information to dump styles
}

/* ------------------------------------------------------------------------*/



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

void ComputeUncertainty::compute_peratom()
{
  int i,j,k,ii,jj,kk,n,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,value;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow centro array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(centro);
    nmax = atom->nmax;
    memory->create(centro,nmax,"uncertainty:centro");
    vector_atom = centro;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // npairs = number of unique pairs

  int nhalf = nnn/2;
  int npairs = nnn * (nnn-1) / 2;
  double *pairs = new double[npairs];

  // compute centro-symmetry parameter for each atom in group
  // use full neighbor list

  double **x = atom->x;
  int *mask = atom->mask;
  double cutsq = force->pair->cutforce * force->pair->cutforce;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // insure distsq and nearest arrays are long enough

      if (jnum > maxneigh) {
        memory->destroy(distsq);
        memory->destroy(nearest);
        maxneigh = jnum;
        memory->create(distsq,maxneigh,"uncertainty:distsq");
        memory->create(nearest,maxneigh,"uncertainty:nearest");
      }

      // loop over list of all neighbors within force cutoff
      // distsq[] = distance sq to each
      // nearest[] = atom indices of neighbors

      n = 0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          distsq[n] = rsq;
          nearest[n++] = j;
        }
      }

      // if not nnn neighbors, centro = 0.0

      if (n < nnn) {
        centro[i] = 0.0;
        continue;
      }


      // R = Ri + Rj for each of npairs i,j pairs among nnn neighbors
      // pairs = squared length of each R

      n = 0;
      for (j = 0; j < nnn; j++) {
        jj = nearest[j];
        for (k = j+1; k < nnn; k++) {
          kk = nearest[k];
          delx = x[jj][0] + x[kk][0] - 2.0*xtmp;
          dely = x[jj][1] + x[kk][1] - 2.0*ytmp;
          delz = x[jj][2] + x[kk][2] - 2.0*ztmp;
          pairs[n++] = delx*delx + dely*dely + delz*delz;
        }
      }


      // centrosymmetry = sum of nhalf smallest squared values

      value = 0.0;
      for (j = 0; j < nhalf; j++) value += pairs[j];
      centro[i] = value;
    } else centro[i] = 0.0;
  }

  delete [] pairs;

  //user stuff
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
