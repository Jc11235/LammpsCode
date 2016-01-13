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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_agni.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;

/* ---------------------------------------------------------------------- */

PairAgni::PairAgni(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  writedata = 1;
  start = true;
}

/* ---------------------------------------------------------------------- */

PairAgni::~PairAgni()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
   // memory->destroy(epsilon);
    //memory->destroy(sigma);

    //memory->destroy(offset);

    //memory->destroy(nTrain);

  }
}

/* ---------------------------------------------------------------------- */

void PairAgni::compute(int eflag, int vflag)
{
  
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  double kx, ky, kz;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  if(start == true)
  {
    readUserFile(); //user input files
    start = false;
  }

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  //our stuff 
 

  for (ii = 0; ii < inum; ii++) {
    vector<double> Vx(eta.size()),Vy(eta.size()),Vz(eta.size()); 
    for(int n = 0; n < eta.size(); n++)
    {
      Vx.at(n) = 0;
      Vy.at(n) = 0;
      Vz.at(n) = 0;
    }

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) 
    {
      j = jlist[jj];

      factor_lj = 0;
      //j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

     // cout<<"dx: "<<delx<<" dy: "<<dely<<" dz: "<<delz<<endl;
      //our stuff
      double cF,wX,wY,wZ;

      //cout<<"Rc: "<<cutsq[itype][jtype]<<endl;
      cF = .5*(cos((M_PI*sqrt(rsq))/sqrt(cutsq[itype][jtype])) + 1.0);
      wX = delx/sqrt(rsq);
      wY = dely/sqrt(rsq);
      wZ = delz/sqrt(rsq);

      for(int n = 0; n < eta.size(); n++)
      {
        Vx.at(n) += wX*cF*exp(-(eta.at(n)*rsq)); 
        Vy.at(n) += wY*cF*exp(-(eta.at(n)*rsq));
        Vz.at(n) += wZ*cF*exp(-(eta.at(n)*rsq)); 
      }
    }

    //our stuff - Force prediction
    //cout<<"atom#: "<<i<<" x: "<<xtmp<<" y: "<<ytmp<<" z: "<<ztmp<<endl;
    for(int n = 0; n < nTrain; n++)
    {
      double kx = 0.0;
      double ky = 0.0;
      double kz = 0.0;

      for(int m = 0; m < eta.size(); m++)
      {        
        kx += pow((Vx.at(m) - xU[m][n]), 2.0);  
        ky += pow((Vy.at(m) - xU[m][n]), 2.0);
        kz += pow((Vz.at(m) - xU[m][n]), 2.0);
        //cout<<" xU: "<<xU[m][n]<<endl;
      }
      
      //cout<<" alpha: "<<alpha.at(n)<<"  y:"<<yU.at(n)<<endl;      
      
      //cout<<"alpha: "<<alpha.at(n)<<endl;
      f[i][0] += alpha.at(n)*exp(-kx/(2.0*pow(sigma1,2.0)));
      f[i][1] += alpha.at(n)*exp(-ky/(2.0*pow(sigma1,2.0)));
      f[i][2] += alpha.at(n)*exp(-kz/(2.0*pow(sigma1,2.0)));
      //cout<<setprecision(15)<<"fi0: "<<alpha.at(n)*exp(-kx/(2.0*pow(sigma1,2.0)))<<endl;
     // cout<<"kx = "<<kx<<" ky = "<<ky<<"kz = "<<kz<<endl;
    }
    //cout<<"fi0: "<<f[i][0]<<endl;
    f[i][0] += b;
    f[i][1] += b;
    f[i][2] += b;

    /*for(int m = 0; m < eta.size(); m++)
      {
        cout<<"Vx: " <<Vx[m]<<" Vy: "<<Vy[m]<<" Vz: "<<Vz[m]<<endl;
      } */  
    //cout<<"fx: " <<f[i][0]<<" fy: "<<f[i][1]<<" fz: "<<f[i][2]<<endl;

    //cout<<"\n\n";
  }
  //end of our stuff

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairAgni::compute_inner()
{
  
}

/* ---------------------------------------------------------------------- */

void PairAgni::compute_middle()
{
  
}

/* ---------------------------------------------------------------------- */

void PairAgni::compute_outer(int eflag, int vflag)
{
 
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairAgni::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  //memory->create(epsilon,n+1,n+1,"pair:epsilon");
 // memory->create(sigma,n+1,n+1,"pair:sigma");

  //memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairAgni::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairAgni::coeff(int narg, char **arg)
{
  inputFile = arg[2]; //sets the user input filename from in.eam


  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[3]); //changed from arg[2] to arg[3] since filename is now arg[2]
  double sigma_one = force->numeric(FLERR,arg[4]);//changed from arg[3] to arg[4] 

  double cut_one = cut_global;
  if (narg == 6) cut_one = force->numeric(FLERR,arg[5]); //changed from narg = 5 to narg = 6, arg[4] to arg[5]

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //epsilon[i][j] = epsilon_one;
      //sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAgni::init_style()
{
  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;

    if (respa == 0) 
      {
        irequest = neighbor->request(this,instance_me);
        neighbor->requests[irequest]->half = 0;
        neighbor->requests[irequest]->full = 1;
      }
    else if (respa == 1) {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    } 
    else {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 2;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respamiddle = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    }

  } else {
  //  cout<<"Reading full list"<<endl;
    irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
   // cout<<"Finished reading full list"<<endl;
  }

  // set rRESPA cutoffs

  if (strstr(update->integrate_style,"respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairAgni::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAgni::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    //epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               //sigma[i][i],sigma[j][j]);
    //sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  if (offset_flag) {
    //double ratio = sigma[i][j] / cut[i][j];
    //offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } //else offset[i][j] = 0.0;

  //offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && cut[i][j] < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    //double sig2 = sigma[i][j]*sigma[i][j];
    //double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    //etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
     // sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
   // ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
     // sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }

  return cut[i][j];
}
//reads in the user input file, reads values in as references to objects
void PairAgni::readUserFile()
{
  ifstream infile(inputFile.c_str());//input file

  string line = "";  //initializes line
  string varSet = ""; //determines where the xU values will start

  bool changeJ; //changes entered 0's to avoid erasing
  bool tempBreak; //boolean for deleting the first element of the item vector
  bool resizeXU = true; //resizes xU based on nTrain
  bool varStart = true; //temp boolean to start 
  
  double j; //temp variables for string to double values
  
  int xUStart = 1000000; //arbitrarily large value to ensure this only happens when it needs to
  int xUcount = 0;
  int count = 0; //line countre

  while(getline(infile,line))
  {    
    changeJ = false;
    tempBreak = false;
    varStart = true;
    varSet = "";

    if(line != "clear" && line.size() > 0)//allows the user to clear and ignores empty lines
    {
      //reads in the line and then splits it at the spaces
      vector<string> elements; //elements will become the strings after the split, but in vector form for easy accessing
      char delim = ' '; //deliminator tells getline where to split the string
      string item; //item is a null string variable that gets initialized to strings after being split
      stringstream ss(line); //ss is the sstream that stores line

      while(getline(ss,item,delim))
      {
        elements.push_back(item); //pushes the line to the elements vector
      }
      //erases verctor elements that correspond to spaces and initialize variables
      for(int i = 0; i < elements.size(); i++)
      {        
        j = atof(elements.at(i).c_str()); //converts a const *char to a float
        
        if(i == 0 && varStart == true)
        {
          varSet = elements.at(0); //stores var name 
          varStart = false; //ensures thi is only done once per line
        }

        //determines all different elements that are in the system
        if(varSet == "Dataset" && i == 0)
        {
          	for(int e = 0; e < elements.size(); e++)
          	{
          		string buf; // Have a buffer string
    			stringstream ss(elements.at(e)); // Insert the string into a stream

    			if(ss >> buf)
        			elementList.push_back(buf); //pushes the elements name to the element list
          	} 

          	elementList.erase(elementList.begin()); //gets rid of the word Dataset               	
        }
        if(elements.at(i) == "0.0") //changes the value of added 0's to avoid be erased
        {
          j = 1.0;
          changeJ = true;
        } 
        //if(i == 0 && count != 4 && tempBreak == false && count <= xUStart) //erases the first element
        if(i == 0 && count != 4 && tempBreak == false)
        {
            elements.erase(elements.begin());
            i--;//modifies the index to avoid segmentation faults
            tempBreak = true;
        }       
        else if(j == 0 && count != 4)//erases strings
        {         
            elements.erase(elements.begin());
            i--;//modifies the index to avoid segmentation faults              
        }
        else if(j != 0 && tempBreak == true) //initializes values
        {
          if(changeJ == true)//restores the original value of the changed 0's
            j = atof(elements.at(i).c_str());
          if(varSet == "Rc") //Rc
            Rc =  j; 
          else if(varSet == "eta") //eta
            eta.push_back(j);
          else if(varSet == "sigma")//sigma
            sigma1= j;
          else if(varSet == "lambda")//lambda
            lambda = j;
          else if(varSet == "b") //b
            b = j;
          else if(varSet == "n_train")//n-train
            nTrain = (int) j;
          else if(varSet == "endVar")
            xUStart = count + 1;  //sets when xU will start          
          else if(count > xUStart && count <= xUStart + nTrain)//xU,yU,alpha
          {             
            if(i < elements.size()-2)//xU
            {
              //resizes xU based on nTrain and eta.size()
              if(resizeXU == true)
              {
                xU.resize(eta.size());//resize number of columes in xU
                for(int l = 0; l < xU.size(); l++)
                {
                  xU[l].resize(nTrain);//resizes number of rows in xU
                }
                resizeXU = false;
              }              
              for(int k = 0; k < eta.size(); k++)
              {
                if(i == k)
                  xU[k].at(xUcount) = j; //adds j to the correct xU index
              }
            }
            else if(i == elements.size()-2)//yU
              yU.push_back(j);
            else if(i == elements.size()-1)//alpha
              alpha.push_back(j);
          }
        }
      }
      if(count > xUStart) //counts until nTrain
        xUcount++;
    }
    count++;
  }

  //send variables to nodes through MPI
  MPI_Bcast(&Rc,1,MPI_DOUBLE,0,world);//Rc
  MPI_Bcast(&sigma1,1,MPI_DOUBLE,0,world);//sigma
  MPI_Bcast(&lambda,1,MPI_DOUBLE,0,world);//lambda
  MPI_Bcast(&b,1,MPI_DOUBLE,0,world);//b
  MPI_Bcast(&nTrain,1,MPI_INT,0,world);//ntrain
  for(int k = 0; k < xU.size(); k++)
  {
  	for(int p = 0; p < xU.at(0).size(); p++)
  	{
		MPI_Bcast(&xU[k][p],1,MPI_DOUBLE,0,world);//xU
  	}
  }  
  for(int p = 0; p < yU.size(); p++)
  {
	MPI_Bcast(&yU[p],1,MPI_DOUBLE,0,world);//yU
  }  
  for(int k = 0; k < alpha.size(); k++)
  {
	MPI_Bcast(&alpha[k],1,MPI_DOUBLE,0,world);//alpha
  }
  for(int k = 0; k < eta.size(); k++)
  {
  	MPI_Bcast(&eta[k],1,MPI_DOUBLE,0,world);//eta
  } 

 	
}
/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairAgni::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
     
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairAgni::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {         
        }       
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairAgni::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairAgni::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairAgni::write_data(FILE *fp)
{
  //for (int i = 1; i <= atom->ntypes; i++)
    //fprintf(fp,"%d %g %g\n",i,epsilon[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairAgni::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      //fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut[i][j]);
    	int t = 0;
}

/* ---------------------------------------------------------------------- */

double PairAgni::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  /*double r2inv,r6inv,forcelj,philj;

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  fforce = factor_lj*forcelj*r2inv;

  philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
    offset[itype][jtype];
    */
  return 0;
}

/* ---------------------------------------------------------------------- */

void *PairAgni::extract(const char *str, int &dim)
{
  dim = 2;
  //if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  //if (strcmp(str,"sigma") == 0) return (void *) sigma;
  return NULL;
}
