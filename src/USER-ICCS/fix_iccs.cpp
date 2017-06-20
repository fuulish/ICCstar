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
#define EPS 1.E-02

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
  nvector = 0;
  peratom = NULL;
  vectors = NULL;
  // coeffs = NULL;

  //FUX| we need to communicate charges, make sure that happens --> see reax/c
  atom->add_callback(0);

  comm_forward = 1;

  damp = 0.9;
  conv = EPS;
  niter = 20;

  qinit = 0;
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

  // memory->destroy(contrast);

  memory->destroy(peratom);
  for (int m = 0; m < nvector; m++) memory->destroy(vectors[m]);
  memory->sfree(vectors);

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

  //FUX| create memory for contrast
  int natoms = atom->natoms;

  // memory->create(contrast,natoms+1,"iccs:contrast");
  // memory->create(qprv,natoms+1,"iccs:qprv");
  // memory->create(qnxt,natoms+1,"iccs:qnxt");

  // contrast = memory->create(contrast,natoms+1,"iccs:contrast");
  qprv = memory->create(qprv,natoms+1,"iccs:qprv");
  qnxt = memory->create(qnxt,natoms+1,"iccs:qnxt");

  add_vector(1);
  contrast = request_vector(0);

  //FUX| loopify; check if flag is 0/1 for integer/float
  reset_vectors();
  calculate_contrast();
}

void FixICCS::reset_vectors()
{
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

  contrast = request_vector(0);
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
    } else if (strcmp(arg[iarg],"niter") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      niter = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix_modify command");
  }

  return 2;
}

void FixICCS::setup_pre_force(int vflag)
{
  reset_vectors();
  initialize_charges();
  pre_force(vflag);
}

void FixICCS::pre_force(int vflag)
{

  reset_vectors();
  run();
  //FUX| do the actual iterations and then forward communicate charges
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
  //
  //FUX| what happens with non-zero charges???

  int i;
  int converged = 0;

  // calculate_contrast();

  for( i=0; i<niter; i++ ) {

    // post_scf_checks();

    iterate();
    converged = check_convergence();

    if( converged )
      break;
  }

  if( !(converged) )
    error->all(FLERR,"Convergence could not be achieved in maximum number of iterations");

  // if( converged)
  //   if( comm->me == 0 )
  //     printf("ICC* converged in %i iterations.\n", i);

  post_scf_checks();

}

void FixICCS::iterate()
{
  backup_charges();
  c_ef->compute_peratom();

  //FUX| this could become more complicated if fixes like ASPC/DRUDE SCF calculations get involved
  calculate_charges_iccs();
  update_charges();
}

void FixICCS::post_scf_checks()
{
  double qtot = 0.0;
  double qtotall = 0.0;

  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  double *q = atom->q;

  for( i=0; i<nlocal; ++i )
    if( mask[i] & groupbit ) {
      qtot += q[i];
      // printf("Q: %g\n", q[i]);
      // printf("C: %g\n", contrast[i]);
    }

  MPI_Allreduce(&qtot, &qtotall, 1, MPI_DOUBLE, MPI_SUM, world);

  // printf( "Total ICC* charge is %g\n", qtotall);

}

void FixICCS::backup_charges()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  double *q = atom->q;

  for( i=0; i<nlocal; ++i )
    if( mask[i] & groupbit )
      qprv[i] = q[i];
}

//FUX| calculate_charges()
//FUX| update_charges()

void FixICCS::calculate_charges_iccs()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  double **ef = c_ef->array_atom;

  // printf("NEW ITERATION\n");
  for( i=0; i<nlocal; ++i )
    if( mask[i] & groupbit )
      qnxt[i] = contrast[i] * ( ef[i][0]*p_srfx[i] + ef[i][1]*p_srfy[i] + ef[i][2]*p_srfz[i] );

}

void FixICCS::initialize_charges()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *q = atom->q;

  for( i=0; i<nlocal; ++i )
    if( mask[i] & groupbit )
      q[i] = 0.01 * ( (float) rand() / RAND_MAX - 0.5);

  comm->forward_comm_fix(this);
  force->kspace->qsum_qsq();
      
  qinit = 1;
}

void FixICCS::update_charges()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *q = atom->q;
  double onemdamp = 1. - damp;

  for( i=0; i<nlocal; ++i )
    if( mask[i] & groupbit )
      q[i] = onemdamp * qprv[i] + damp * qnxt[i];

  comm->forward_comm_fix(this);
  force->kspace->qsum_qsq();
}

int FixICCS::check_convergence()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *q = atom->q;

  int isnotconv = 0;
  int allisnotconv = 0;

  for( i=0; i<nlocal; ++i )
    if( mask[i] & groupbit )
      if( fabs( ( q[i] - qprv[i] ) / qprv[i] ) > conv )
        isnotconv += 1;

  //FUX| substitute MPI_SUM by MPI_MAX? and check for largest charge
  //FUX| avoids havhing to calculate the sum...
  MPI_Allreduce(&isnotconv,&allisnotconv,1,MPI_INT,MPI_SUM,world);

  if ( allisnotconv )
    return 0;
  else
    return 1;

}

void FixICCS::calculate_contrast()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  // double fpieps = 0.0030119505336064496;
  double fpieps = 0.0030114702717118813;

  for ( i=0; i<nlocal; i++ )
    if( mask[i] & groupbit ) {

      contrast[i] = bulk_perm / TWOPI * p_area[i] * fpieps;

      if ( p_diel[i] < 1 )
        contrast[i] *= -1.;
      else
        contrast[i] *= (bulk_perm - p_diel[i]) / (bulk_perm + p_diel[i]);

    }

}

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

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixICCS::memory_usage()
{
  double bytes = 0.0;
  for (int m = 0; m < nvector; m++)
    bytes += atom->nmax*peratom[m]*sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixICCS::grow_arrays(int nmax)
{
  for (int m = 0; m < nvector; m++)
    memory->grow(vectors[m],peratom[m]*nmax,"minimize:vector");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixICCS::copy_arrays(int i, int j, int delflag)
{
  int m,iper,nper,ni,nj;

  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper*i;
    nj = nper*j;
    for (iper = 0; iper < nper; iper++) vectors[m][nj++] = vectors[m][ni++];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixICCS::pack_exchange(int i, double *buf)
{
  int m,iper,nper,ni;

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper*i;
    for (iper = 0; iper < nper; iper++) buf[n++] = vectors[m][ni++];
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixICCS::unpack_exchange(int nlocal, double *buf)
{
  int m,iper,nper,ni;

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper*nlocal;
    for (iper = 0; iper < nper; iper++) vectors[m][ni++] = buf[n++];
  }
  return n;
}

/* ----------------------------------------------------------------------
   allocate/initialize memory for a new vector with N elements per atom
------------------------------------------------------------------------- */

void FixICCS::add_vector(int n)
{
  memory->grow(peratom,nvector+1,"minimize:peratom");
  peratom[nvector] = n;

  vectors = (double **)
    memory->srealloc(vectors,(nvector+1)*sizeof(double *),"minimize:vectors");
  memory->create(vectors[nvector],atom->nmax*n,"minimize:vector");

  int ntotal = n*atom->nlocal;
  for (int i = 0; i < ntotal; i++) vectors[nvector][i] = 0.0;
  nvector++;
}

/* ----------------------------------------------------------------------
   return a pointer to the Mth vector
------------------------------------------------------------------------- */

double *FixICCS::request_vector(int m)
{
  return vectors[m];
}
