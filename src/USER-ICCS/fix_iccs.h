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

  int pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc);
  void unpack_forward_comm(int n, int first, double *buf);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  void add_vector(int);
  double *request_vector(int);

 protected:
  int dim;                     // dimensionality of the thing to be aspc'ed
  bigint ndoftotal;           // total dof for entire problem

  int nvec;                   // local atomic dof = length of xvec
  double *qty;                // variables for atomic dof, as 1d vector

  int torqueflag,extraflag;

  int pair_compute_flag;            // 0 if pair->compute is skipped
  int kspace_compute_flag;          // 0 if kspace->compute is skipped

 private:
  char *id_ef, *id_diel, *id_area, *id_srfx, *id_srfy, *id_srfz;
  class Compute *c_ef;

  int nvector;
  int *peratom;
  double **vectors;
  double bzr;
  double damp, conv;
  int niter;
  int qinit;

  double bulk_perm;
  double *p_diel, *p_area, *p_srfx, *p_srfy, *p_srfz;
  double *contrast, *qprv, *qnxt;

  void calculate_contrast();
  void run();
  void iterate();
  int check_convergence();
  void calculate_charges_iccs();
  void update_charges();
  void initialize_charges();
};

}

#endif
#endif
