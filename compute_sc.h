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

ComputeStyle(sc,ComputeSC)

#else

#ifndef LMP_COMPUTE_SC_H
#define LMP_COMPUTE_SC_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSC : public Compute {
 public:
  ComputeSC(class LAMMPS *, int, char **);
  ~ComputeSC();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();

 private:
  int nbin;              // # of rdf bins
  int cutflag;           // user cutoff flag
  int npairs;            // # of rdf pairs
  double delr,delrinv;   // bin width and its inverse
  double cutoff_user;    // user-specified cutoff
  double mycutneigh;     // user-specified cutoff + neighbor skin
  int ***rdfpair;        // map 2 type pair to rdf pair for each histo
  int **nrdfpair;        // # of histograms for each type pair
  int *ilo,*ihi,*jlo,*jhi;
  double **hist;         // histogram bins
  double **histall;      // summed histogram bins across all procs
  double **dot_sum;         // dot sum bins
  double **dot_sum_all;      // summed dot sum bins across all procs
  double **sum_i;         // sum i bins
  double **sum_i_all;      // summed  sum i bins across all procs
  double **sum_j;         // sum j bins
  double **sum_j_all;      // summed  sum j bins across all procs
  double **sum_sq_i;         // sum sq i bins
  double **sum_sq_i_all;      // summed  sum sq i bins across all procs
  double **sum_sq_j;         // sum sq j bins
  double **sum_sq_j_all;      // summed  sum sq j bins across all procs

  int *typecount;
  int *icount,*jcount;
  int *duplicates;

  class NeighList *list; // half neighbor list
  void init_norm();
  bigint natoms_old;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute rdf requires a pair style be defined or cutoff specified

UNDOCUMENTED

E: Compure rdf cutoff exceeds ghost atom range - use comm_modify cutoff command

UNDOCUMENTED

W: Compute rdf cutoff less than neighbor cutoff - forcing a needless neighbor list build

UNDOCUMENTED

U: Compute rdf requires a pair style be defined

Self-explanatory.

*/
