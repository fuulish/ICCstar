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

/* ---------------------------------------------------------------------- */

FixICCS::FixICCS(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg)
{

  if (narg < 5) error->all(FLERR,"Illegal fix ICCS command");

  //FU| memory-related
  // nvector = 0;
  // peratom = NULL;
  // vectors = NULL;
  // coeffs = NULL;

  //FUX| we need to communicate charges, make sure that happens --> see reax/c
  atom->add_callback(0);

  comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixICCS::~FixICCS()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored data

  // memory->destroy(peratom);
  // for (int m = 0; m < nvector; m++) memory->destroy(vectors[m]);
  // memory->sfree(vectors);
  // memory->sfree(coeffs);
}

/* ---------------------------------------------------------------------- */

void FixICCS::init()
{
}

void FixICCS::reset_vectors()
{
}

int FixICCS::modify_param(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal fix_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"testing") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      // do something here
      // ...
      iarg += 2;
    } else error->all(FLERR,"Illegal fix_modify command");
  }

  return 2;
}

void FixICCS::setup_pre_force(int vflag)
{
    reset_vectors();
}

void FixICCS::pre_force(int vflag)
{
  //FUX| do the actual iterations and then forward communicate charges
  comm->forward_comm_fix(this);
}

int FixICCS::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
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

