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
   Contributing authors: Bin Xu
------------------------------------------------------------------------- */

#include "compute_sc.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeSC::ComputeSC(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  rdfpair(NULL), nrdfpair(NULL), ilo(NULL), ihi(NULL), jlo(NULL), jhi(NULL),
  hist(NULL), histall(NULL), typecount(NULL), icount(NULL), jcount(NULL),
  duplicates(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute rdf command");

  array_flag = 1;
  extarray = 0;

  nbin = force->inumeric(FLERR,arg[3]);
  if (nbin < 1) error->all(FLERR,"Illegal compute rdf command");

  // optional args
  // nargpair = # of pairwise args, starting at iarg = 4

  cutflag = 0;

  // Need per atom energy
  peatomflag = 1;

  // need to store timestep
  timeflag = 1;

  int iarg;
  for (iarg = 4; iarg < narg; iarg++)
    if (strcmp(arg[iarg],"cutoff") == 0) break;

  int nargpair = iarg - 4;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute sc command");
      cutoff_user = force->numeric(FLERR,arg[iarg+1]);
      if (cutoff_user <= 0.0) cutflag = 0;
      else cutflag = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal compute sc command");
  }

  // pairwise args
  //
  // nargpair is the number of arg about itype and jtype

  if (nargpair == 0) npairs = 1;  // if nargpair == 0, gr for all-all
  else {
    if (nargpair % 2) error->all(FLERR,"Illegal compute sc command");
    npairs = nargpair/2;
  }

  size_array_rows = nbin;
  //size_array_cols = 1 + 2*npairs;
  size_array_cols = 1 + 6*npairs;

  int ntypes = atom->ntypes;
  memory->create(rdfpair,npairs,ntypes+1,ntypes+1,"sc:rdfpair");
  memory->create(nrdfpair,ntypes+1,ntypes+1,"sc:nrdfpair");
  ilo = new int[npairs];
  ihi = new int[npairs];
  jlo = new int[npairs];
  jhi = new int[npairs];

  if (nargpair == 0) {
    ilo[0] = 1; ihi[0] = ntypes;
    jlo[0] = 1; jhi[0] = ntypes;
  } else {
    iarg = 4;
    for (int ipair = 0; ipair < npairs; ipair++) {
      force->bounds(FLERR,arg[iarg],atom->ntypes,ilo[ipair],ihi[ipair]);
      force->bounds(FLERR,arg[iarg+1],atom->ntypes,jlo[ipair],jhi[ipair]);
      if (ilo[ipair] > ihi[ipair] || jlo[ipair] > jhi[ipair])
        error->all(FLERR,"Illegal compute sc command");
      iarg += 2;
    }
  }

  int i,j;
  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      nrdfpair[i][j] = 0; // number of rdfpair for itype and jtype will be used

  int ihisto;
  for (int m = 0; m < npairs; m++)
    for (i = ilo[m]; i <= ihi[m]; i++)
      for (j = jlo[m]; j <= jhi[m]; j++) {
        ihisto = nrdfpair[i][j]++; 
        rdfpair[ihisto][i][j] = m; // know itype and jtype, iterate over range(nrdfpair[i][j]) to get index of rdfpair
      }

  memory->create(hist,npairs,nbin,"sc:hist");
  memory->create(histall,npairs,nbin,"sc:histall");
  memory->create(dot_sum,npairs,nbin,"sc:dot_sum");
  memory->create(dot_sum_all,npairs,nbin,"sc:dot_sum_all");
  memory->create(sum_i,npairs,nbin,"sc:sum_i");
  memory->create(sum_i_all,npairs,nbin,"sc:sum_i_all");
  memory->create(sum_j,npairs,nbin,"sc:sum_j");
  memory->create(sum_j_all,npairs,nbin,"sc:sum_j_all");
  memory->create(sum_sq_i,npairs,nbin,"sc:sum_sq_i");
  memory->create(sum_sq_i_all,npairs,nbin,"sc:sum_sq_i_all");
  memory->create(sum_sq_j,npairs,nbin,"sc:sum_sq_j");
  memory->create(sum_sq_j_all,npairs,nbin,"sc:sum_sq_j_all");
  memory->create(array,nbin,1+6*npairs,"sc:array");
  typecount = new int[ntypes+1];
  icount = new int[npairs];
  jcount = new int[npairs];
  duplicates = new int[npairs];

  dynamic = 0;
  natoms_old = 0;
}

/* ---------------------------------------------------------------------- */

ComputeSC::~ComputeSC()
{
  memory->destroy(rdfpair);
  memory->destroy(nrdfpair);
  delete [] ilo;
  delete [] ihi;
  delete [] jlo;
  delete [] jhi;
  memory->destroy(hist);
  memory->destroy(histall);
  memory->destroy(dot_sum);
  memory->destroy(dot_sum_all);
  memory->destroy(sum_i);
  memory->destroy(sum_i_all);
  memory->destroy(sum_j);
  memory->destroy(sum_j_all);
  memory->destroy(sum_sq_i);
  memory->destroy(sum_sq_i_all);
  memory->destroy(sum_sq_j);
  memory->destroy(sum_sq_j_all);
  memory->destroy(array);
  delete [] typecount;
  delete [] icount;
  delete [] jcount;
  delete [] duplicates;
}

/* ---------------------------------------------------------------------- */

void ComputeSC::init()
{

  if (!force->pair)
    error->all(FLERR,"Compute sc requires a pair style be defined ");

  if (cutflag) {
    double skin = neighbor->skin;
    mycutneigh = cutoff_user + skin;

    double cutghost;            // as computed by Neighbor and Comm
    if (force->pair)
      cutghost = MAX(force->pair->cutforce+skin,comm->cutghostuser);
    else
      cutghost = comm->cutghostuser;

    if (mycutneigh > cutghost)
      error->all(FLERR,"Compute sc cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");
    if (force->pair && mycutneigh < force->pair->cutforce + skin)
      if (comm->me == 0)
        error->warning(FLERR,"Compute sc cutoff less than neighbor cutoff - "
                       "forcing a needless neighbor list build");

    delr = cutoff_user / nbin;
  } else delr = force->pair->cutforce / nbin;

  delrinv = 1.0/delr;

  // set 1st column of output array to bin coords

  for (int i = 0; i < nbin; i++)
    array[i][0] = (i+0.5) * delr;

  // initialize normalization, finite size correction, and changing atom counts

  natoms_old = atom->natoms;
  dynamic = group->dynamic[igroup];
  if (dynamic_user) dynamic = 1;
  init_norm();

  // need an occasional half neighbor list
  // if user specified, request a cutoff = cutoff_user + skin
  // skin is included b/c Neighbor uses this value similar
  //   to its cutneighmax = force cutoff + skin
  // also, this NeighList may be used by this compute for multiple steps
  //   (until next reneighbor), so it needs to contain atoms further
  //   than cutoff_user apart, just like a normal neighbor list does

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
  if (cutflag) {
    neighbor->requests[irequest]->cut = 1;
    neighbor->requests[irequest]->cutoff = mycutneigh;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSC::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSC::init_norm()
{
  int i,j,m;

  // count atoms of each type that are also in group

  const int nlocal = atom->nlocal;
  const int ntypes = atom->ntypes;
  const int * const mask = atom->mask;
  const int * const type = atom->type;

  for (i = 1; i <= ntypes; i++) typecount[i] = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) typecount[type[i]]++;

  // icount = # of I atoms participating in I,J pairs for each histogram
  // jcount = # of J atoms participating in I,J pairs for each histogram
  // duplicates = # of atoms in both groups I and J for each histogram

  for (m = 0; m < npairs; m++) {
    icount[m] = 0;
    for (i = ilo[m]; i <= ihi[m]; i++) icount[m] += typecount[i];
    jcount[m] = 0;
    for (i = jlo[m]; i <= jhi[m]; i++) jcount[m] += typecount[i];
    duplicates[m] = 0;
    for (i = ilo[m]; i <= ihi[m]; i++)
      for (j = jlo[m]; j <= jhi[m]; j++)
        if (i == j) duplicates[m] += typecount[i];
  }

  int *scratch = new int[npairs];
  MPI_Allreduce(icount,scratch,npairs,MPI_INT,MPI_SUM,world);
  for (i = 0; i < npairs; i++) icount[i] = scratch[i];
  MPI_Allreduce(jcount,scratch,npairs,MPI_INT,MPI_SUM,world);
  for (i = 0; i < npairs; i++) jcount[i] = scratch[i];
  MPI_Allreduce(duplicates,scratch,npairs,MPI_INT,MPI_SUM,world);
  for (i = 0; i < npairs; i++) duplicates[i] = scratch[i];
  delete [] scratch;
}

/* ---------------------------------------------------------------------- */

void ComputeSC::compute_array()
{
  invoked_peratom = update->ntimestep;
  if (update->eflag_atom != invoked_peratom)
    error->all(FLERR,"Per-atom energy was not tallied on needed timestep");

  double * eatom = force->pair->eatom; // get peratom energy for pair potential

  int i,j,m,ii,jj,inum,jnum,itype,jtype,ipair,jpair,ibin,ihisto;
  double xtmp,ytmp,ztmp,delx,dely,delz,r;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double factor_lj,factor_coul;

  if (natoms_old != atom->natoms) {
    dynamic = 1;
    natoms_old = atom->natoms;
  }

  // if the number of atoms has changed or we have a dynamic group
  // or dynamic updates are requested (e.g. when changing atom types)
  // we need to recompute some normalization parameters

  if (dynamic) init_norm();

  invoked_array = update->ntimestep;

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero the histogram counts

  for (i = 0; i < npairs; i++)
    for (j = 0; j < nbin; j++){
      hist[i][j] = 0;
      dot_sum[i][j] = 0.;
      sum_i[i][j] = 0.;
      sum_j[i][j] = 0.;
      sum_sq_i[i][j] = 0.;
      sum_sq_j[i][j] = 0.;
    }

  for (i = 0; i < npairs; i++)
    for (j = 0; j < nbin; j++)
      dot_sum[i][j] = 0;

  // tally the RDF
  // both atom i and j must be in fix group
  // itype,jtype must have been specified by user
  // consider I,J as one interaction even if neighbor pair is stored on 2 procs
  // tally I,J pair each time I is central atom, and each time J is central

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double e_i,e_j;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    e_i = eatom[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      // if both weighting factors are 0, skip this pair
      // could be 0 and still be in neigh list for long-range Coulombics
      // want consistency with non-charged pairs which wouldn't be in list

      if (factor_lj == 0.0 && factor_coul == 0.0) continue;

      if (!(mask[j] & groupbit)) continue;
      jtype = type[j];
      ipair = nrdfpair[itype][jtype];
      jpair = nrdfpair[jtype][itype];
      if (!ipair && !jpair) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      r = sqrt(delx*delx + dely*dely + delz*delz);
      ibin = static_cast<int> (r*delrinv);
      if (ibin >= nbin) continue;
      e_j = eatom[j];

      for (ihisto = 0; ihisto < ipair; ihisto++) {
        m = rdfpair[ihisto][itype][jtype];
        hist[m][ibin] += 1.0;
	dot_sum[m][ibin] += e_i*e_j;
	sum_i[m][ibin] += e_i;
	sum_j[m][ibin] += e_j;
	sum_sq_i[m][ibin] += e_i*e_i;
	sum_sq_j[m][ibin] += e_j*e_j;
      }
      if (newton_pair || j < nlocal) {
        for (ihisto = 0; ihisto < jpair; ihisto++) {
          m = rdfpair[ihisto][jtype][itype];
          hist[m][ibin] += 1.0;
	  dot_sum[m][ibin] += e_i*e_j;
	  sum_i[m][ibin] += e_i;
	  sum_j[m][ibin] += e_j;
	  sum_sq_i[m][ibin] += e_i*e_i;
	  sum_sq_j[m][ibin] += e_j*e_j;
        }
      }
    }
  }

  // sum histograms across procs

  MPI_Allreduce(hist[0],histall[0],npairs*nbin,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(dot_sum[0],dot_sum_all[0],npairs*nbin,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(sum_i[0],sum_i_all[0],npairs*nbin,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(sum_j[0],sum_j_all[0],npairs*nbin,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(sum_sq_i[0],sum_sq_i_all[0],npairs*nbin,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(sum_sq_j[0],sum_sq_j_all[0],npairs*nbin,MPI_DOUBLE,MPI_SUM,world);

  // convert counts to g(r) and coord(r) and copy into output array
  // vfrac = fraction of volume in shell m
  // npairs = number of pairs, corrected for duplicates
  // duplicates = pairs in which both atoms are the same

  //double constant,vfrac,gr,ncoord,rlower,rupper,normfac;

  if (domain->dimension == 3) {
  //  constant = 4.0*MY_PI / (3.0*domain->xprd*domain->yprd*domain->zprd);

    for (m = 0; m < npairs; m++) {
      //normfac = (icount[m] > 0) ? static_cast<double>(jcount[m])
      //          - static_cast<double>(duplicates[m])/icount[m] : 0.0;
      //ncoord = 0.0;
      for (ibin = 0; ibin < nbin; ibin++) {
      //  rlower = ibin*delr;
      //  rupper = (ibin+1)*delr;
      //  vfrac = constant * (rupper*rupper*rupper - rlower*rlower*rlower);
      //  if (vfrac * normfac != 0.0)
      //    gr = histall[m][ibin] / (vfrac * normfac * icount[m]);
      //  else gr = 0.0;
      //  if (icount[m] != 0)
      //    ncoord += gr * vfrac * normfac;
        array[ibin][1+6*m] = histall[m][ibin];
        array[ibin][2+6*m] = dot_sum_all[m][ibin];
        array[ibin][3+6*m] = sum_i_all[m][ibin];
        array[ibin][4+6*m] = sum_j_all[m][ibin];
        array[ibin][5+6*m] = sum_sq_i_all[m][ibin];
        array[ibin][6+6*m] = sum_sq_j_all[m][ibin];
        //array[ibin][1+2*m] = gr;
        //array[ibin][2+2*m] = ncoord;
      }
    }

  } else {
    //constant = MY_PI / (domain->xprd*domain->yprd);

    for (m = 0; m < npairs; m++) {
      //ncoord = 0.0;
      //normfac = (icount[m] > 0) ? static_cast<double>(jcount[m])
      //          - static_cast<double>(duplicates[m])/icount[m] : 0.0;
      for (ibin = 0; ibin < nbin; ibin++) {
        //rlower = ibin*delr;
        //rupper = (ibin+1)*delr;
        //vfrac = constant * (rupper*rupper - rlower*rlower);
        //if (vfrac * normfac != 0.0)
        //  gr = histall[m][ibin] / (vfrac * normfac * icount[m]);
        //else gr = 0.0;
        //if (icount[m] != 0)
        //  ncoord += gr * vfrac * normfac;
        //array[ibin][1+2*m] = gr;
        //array[ibin][2+2*m] = ncoord;
        array[ibin][1+6*m] = histall[m][ibin];
        array[ibin][2+6*m] = dot_sum_all[m][ibin];
        array[ibin][3+6*m] = sum_i_all[m][ibin];
        array[ibin][4+6*m] = sum_j_all[m][ibin];
        array[ibin][5+6*m] = sum_sq_i_all[m][ibin];
        array[ibin][6+6*m] = sum_sq_j_all[m][ibin];
      }
    }
  }
}

