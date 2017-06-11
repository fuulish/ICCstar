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
   Contributing author: Frank Uhlig (ICP Stuttgart)
   Sources: Stefan Kesselheim, Marcello Sega, Christian Holm
            The ICC ⋆ Algorithm: A fast way to include dielectric
            boundary effects into molecular dynamics simulations
            https://arxiv.org/pdf/1003.1271.pdf
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_iccs.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "compute.h"
#include "neighbor.h"
#include "domain.h"
#include "modify.h"
#include "pair.h"
#include "output.h"
#include "thermo.h"
#include "timer.h"
#include "error.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define TWOPI 6.283185307179586
#define EPS 1.E-04

/* ---------------------------------------------------------------------- */

FixICCS::FixICCS(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg)
{

  if (narg < 10) error->all(FLERR,"Illegal fix ICCS command");

  //FUX| need input:
  //      - compute/efield
  //      - property names: 1 dielectric 3 components surface normal

  bulk_perm = force->numeric(FLERR,arg[3]);

  int n = strlen(arg[4]) + 1;
  id_ef = new char[n];
  strcpy(id_ef,arg[4]);

  n = strlen(arg[5]) + 1;
  id_diel = new char[n];
  strcpy(id_diel, arg[5]);

  n = strlen(arg[6]) + 1;
  id_area = new char[n];
  strcpy(id_area, arg[6]);

  n = strlen(arg[7]) + 1;
  id_srfx = new char[n];
  strcpy(id_srfx, arg[7]);

  n = strlen(arg[8]) + 1;
  id_srfy = new char[n];
  strcpy(id_srfy, arg[8]);

  n = strlen(arg[9]) + 1;
  id_srfz = new char[n];
  strcpy(id_srfz, arg[9]);

  //FU| memory-related
  // nvector = 0;
  // peratom = NULL;
  // vectors = NULL;
  // coeffs = NULL;

  //FUX| we need to communicate charges, make sure that happens --> see reax/c
  atom->add_callback(0);

  comm_forward = 1;

  damp = 0.9;
  conv = EPS;
  niter = 20;
}

/* ---------------------------------------------------------------------- */

FixICCS::~FixICCS()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  delete [] id_ef;
  delete [] id_diel;
  delete [] id_area;
  delete [] id_srfx;
  delete [] id_srfy;
  delete [] id_srfz;

  memory->destroy(contrast);

  memory->destroy(qprv);
  memory->destroy(qnxt);

  // delete locally stored data

  // memory->destroy(peratom);
  // for (int m = 0; m < nvector; m++) memory->destroy(vectors[m]);
  // memory->sfree(vectors);
  // memory->sfree(coeffs);
}

/* ---------------------------------------------------------------------- */

void FixICCS::init()
{
  //FUX| find the fixes and corresponding data that we care about
  //FUX| surface normals and dielectric contrast
  //
  //FUX| useful lines
  //     index[nvalue] = atom->find_custom(&arg[iarg][2],tmp);
  //     atom->dvector[index[j]][m] = atof(values[j+1]);

  //FUX| handle if compute id not found
  int ief = modify->find_compute(id_ef);
  if (ief < 0)
    error->all(FLERR,"Compute ID for fix efield/atom does not exist");

  c_ef = modify->compute[ief];

  //FUX| loopify; check if flag is 0/1 for integer/float
  int index, flag;
  
  index = atom->find_custom(id_diel, flag);
  p_diel = atom->dvector[index];
  
  index = atom->find_custom(id_area, flag);
  p_area = atom->dvector[index];
  
  index = atom->find_custom(id_srfx, flag);
  p_srfx = atom->dvector[index];
  
  index = atom->find_custom(id_srfy, flag);
  p_srfy = atom->dvector[index];

  index = atom->find_custom(id_srfz, flag);
  p_srfz = atom->dvector[index];

  //FUX| create memory for contrast
  int natoms = atom->natoms;
  memory->create(contrast,natoms,"iccs:contrast");

  memory->create(qprv,natoms,"iccs:qprv");
  memory->create(qnxt,natoms,"iccs:qnxt");

  calculate_contrast();

}

void FixICCS::reset_vectors()
{
}

int FixICCS::modify_param(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal fix_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"damp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      damp = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"conv") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      conv = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix_modify command");
  }

  return 2;
}

void FixICCS::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

void FixICCS::pre_force(int vflag)
{

  reset_vectors();
  run();
  //FUX| do the actual iterations and then forward communicate charges
  comm->forward_comm_fix(this);
}

int FixICCS::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

void FixICCS::run()
{
  //FUX|  keep local array with initial charges qprv
  //      --> necessary for convergence check
  //      
  //      keep local array with nxt charges qnxt

  int i;

  for( i=0; i<niter; i++ ) {
    iterate();
    check_convergence();
  }
}

void FixICCS::check_convergence()
{
}

void FixICCS::iterate()
{
}

void FixICCS::calculate_contrast()
{
  int i;
  int n = atom->natoms;

  for ( i=0; i<n; i++ ) {

    contrast[i] = bulk_perm / TWOPI * p_area[i];

    if ( p_diel[i] < 1 )
      contrast[i] = -1.;
    else
      contrast[i] = (bulk_perm - p_diel[i]) / (bulk_perm + p_diel[i]);

    // printf("CONTRAST %i %f\n", i, contrast[i]);
  }

}

// /* ---------------------------------------------------------------------- */
// 
// int FixICCS::pack_reverse_comm(int n, int first, double *buf)
// {
//   int i, m;
//   for(m = 0, i = first; m < n; m++, i++) buf[m] = q[i];
//   return n;
// }
// 
// /* ---------------------------------------------------------------------- */
// 
// void FixICCS::unpack_reverse_comm(int n, int *list, double *buf)
// {
//   for(int m = 0; m < n; m++) q[list[m]] += buf[m];
// }

/* ---------------------------------------------------------------------- */

int FixICCS::pack_forward_comm(int n, int *list, double *buf,
                                  int pbc_flag, int *pbc)
{
  int m;

  for(m = 0; m < n; m++) buf[m] = atom->q[list[m]];

  return n;
}

/* ---------------------------------------------------------------------- */

void FixICCS::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  for(m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
}

