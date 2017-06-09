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

/* FUDO| this should not be done, because we only want to use the derived classes
 */

#ifdef FIX_CLASS

FixStyle(iccs,FixICCS)    // This registers this fix class with LAMMPS.

#else

#ifndef LMP_FIX_ICCS_H
#define LMP_FIX_ICCS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixICCS : public Fix {
 public:
  FixICCS(class LAMMPS *, int, char **);
  ~FixICCS();
  void init();
  void reset_vectors();
  void pre_force(int);
  virtual void setup_pre_force(int);
  int modify_param(int narg, char **arg);

  int setmask();

  void add_vector(int);
  double *request_vector(int);

  int pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc);
  void unpack_forward_comm(int n, int first, double *buf);

 protected:
  int dim;                     // dimensionality of the thing to be aspc'ed
  bigint ndoftotal;           // total dof for entire problem

  int nvec;                   // local atomic dof = length of xvec
  double *qty;                // variables for atomic dof, as 1d vector

  int torqueflag,extraflag;

  int pair_compute_flag;            // 0 if pair->compute is skipped
  int kspace_compute_flag;          // 0 if kspace->compute is skipped

 private:
  int nvector;
  int *peratom;
  double **vectors;
  double bzr;
};

}

#endif
#endif
