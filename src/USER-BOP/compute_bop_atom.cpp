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

#include <math.h>
#include <string.h>

#include "compute_bop_atom.h"
#include "spherical_harmonics.h"

#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "fix_store.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "force.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"


//---------------------------------------- bop_calculator_constants.h ----------------------------------

// siatka interpolacyjna metody FSI
// liczba wezlow dla wielomianow Legendre'a
// liczba wezlow dla funkcji trygonometrycznych

const int BOP_CALCULATOR_N_KNOTS_LEGENDRE_0      = 36;
const int BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_0 = 36;

const int BOP_CALCULATOR_N_KNOTS_LEGENDRE_1      = 600;
const int BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_1 = 600;

const int BOP_CALCULATOR_N_KNOTS_LEGENDRE_2      = 1200;
const int BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_2 = 1200;

const int BOP_CALCULATOR_N_KNOTS_LEGENDRE_3      = 2400;
const int BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_3 = 2400;

const int BOP_CALCULATOR_N_KNOTS_LEGENDRE_4      = 4800;
const int BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_4 = 4800;

const int BOP_CALCULATOR_N_KNOTS_LEGENDRE_5      = 9600;
const int BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_5 = 9600;


// symbole 3j Wignera ( l l l ) ( m1 m2 m3 )
// zakodowane "na twardo"

//***************************
// dla l = 4:

// wchodzace z czynnikiem 1x
// 1) ( 0 0 0 )
const double wigner_3j_4_0_0_0 = 3.0 * sqrt( 2.0 / 1001.0 );

// wchodzace z czynnikiem 6x lub -6x
// 2) ( -1 0 1 )
const double wigner_3j_4_m1_0_1 = - 3.0 * sqrt( 1.0 / 2002.0 );
// 3) ( -2 0 2 )
const double wigner_3j_4_m2_0_2 = - ( 1.0 / 3.0 ) * sqrt( 11.0 / 182.0 );
// 4) ( -3 0 3 )
const double wigner_3j_4_m3_0_3 = sqrt( 7.0 / 286.0 );
// 5) ( -4 0 4 )
const double wigner_3j_4_m4_0_4 = ( 1.0 / 3.0 ) * sqrt( 14.0 / 143.0 );

// wchodzace z czynnikiem 12x
// 6) ( -4 1 3 )
const double wigner_3j_4_m4_1_3 = - ( 1.0 / 3.0 ) * sqrt( 35.0 / 143.0 );

// wchodzace z czynnikiem -12x
// 7) ( -3 1 2 )
const double wigner_3j_4_m3_1_2 = - ( 1.0 / 3.0 ) * sqrt( 5.0 / 143.0 );

// wchodzace z czynnikiem 6x
// 8) ( -4 2 2 )
const double wigner_3j_4_m4_2_2 = sqrt( 5.0 / 143.0 );
// 9) ( -2 1 1 )
const double wigner_3j_4_m2_1_1 = 2.0 * sqrt( 5.0 / 1001.0 );

//***************************
// dla l = 6:

// wchodzace z czynnikiem 1x
// 1) ( 0 0 0 )
const double wigner_3j_6_0_0_0 = - 20.0 * sqrt( 1.0 / 46189.0 );

// wchodzace z czynnikiem 6x lub -6x
// 2) ( -1 0 1 )
const double wigner_3j_6_m1_0_1 = 10.0 * sqrt( 1.0 / 46189.0 );
// 3) ( -2 0 2 )
const double wigner_3j_6_m2_0_2 = sqrt( 11.0 / 4199.0 );
// 4) ( -3 0 3 )
const double wigner_3j_6_m3_0_3 = - ( 43.0 / 2.0 ) * sqrt( 1.0 / 46189.0 );
// 5) ( -4 0 4 )
const double wigner_3j_6_m4_0_4 = 4.0 * sqrt( 1.0 / 46189.0 );
// 6) ( -5 0 5 )
const double wigner_3j_6_m5_0_5 = ( 5.0 / 2.0 ) * sqrt( 11.0 / 4199.0 );
// 7) ( -6 0 6 )
const double wigner_3j_6_m6_0_6 = sqrt( 11.0 / 4199.0 );

// wchodzace z czynnikiem 12x
// 8) ( -6 1 5 )
const double wigner_3j_6_m6_1_5 = - sqrt( 77.0 / 8398.0 );
// 9) ( -6 2 4 )
const double wigner_3j_6_m6_2_4 = sqrt( 70.0 / 4199.0 );
// 10) ( -4 1 3 )
const double wigner_3j_6_m4_1_3 = ( 5.0 / 2.0 ) * sqrt( 35.0 / 46189.0 );

// wchodzace z czynnikiem -12x
// 11) ( -5 1 4 )
const double wigner_3j_6_m5_1_4 = - ( 3.0 / 2.0 ) * sqrt( 21.0 / 4199.0 );
// 12) ( -5 2 3 )
const double wigner_3j_6_m5_2_3 = sqrt( 7.0 / 4199.0 );
// 13) ( -3 1 2 )
const double wigner_3j_6_m3_1_2 = 3.0 * sqrt( 21.0 / 92378.0 );

// wchodzace z czynnikiem 6x
// 14) ( -6 3 3 )
const double wigner_3j_6_m6_3_3 = - 2.0 * sqrt( 21.0 / 4199.0 );
// 15) ( -4 2 2 )
const double wigner_3j_6_m4_2_2 = - 6.0 * sqrt( 14.0 / 46189.0 );
// 16) ( -2 1 1 )
const double wigner_3j_6_m2_1_1 = - 2.0 * sqrt( 105.0 / 46189.0 );


// czynniki multiplikacji
// dla l = 4
const double mult_4_0  =   1.0 * wigner_3j_4_0_0_0;
const double mult_4_1  =  -6.0 * wigner_3j_4_m1_0_1;
const double mult_4_2  =   6.0 * wigner_3j_4_m2_0_2;
const double mult_4_3  =  -6.0 * wigner_3j_4_m3_0_3;
const double mult_4_4  =   6.0 * wigner_3j_4_m4_0_4;
const double mult_4_5  =  12.0 * wigner_3j_4_m4_1_3;
const double mult_4_6  = -12.0 * wigner_3j_4_m3_1_2;
const double mult_4_7  =   6.0 * wigner_3j_4_m4_2_2;
const double mult_4_8  =   6.0 * wigner_3j_4_m2_1_1;
// dla l = 6
const double mult_6_0  =   1.0 * wigner_3j_6_0_0_0;
const double mult_6_1  =  -6.0 * wigner_3j_6_m1_0_1;
const double mult_6_2  =   6.0 * wigner_3j_6_m2_0_2;
const double mult_6_3  =  -6.0 * wigner_3j_6_m3_0_3;
const double mult_6_4  =   6.0 * wigner_3j_6_m4_0_4;
const double mult_6_5  =  -6.0 * wigner_3j_6_m5_0_5;
const double mult_6_6  =   6.0 * wigner_3j_6_m6_0_6;
const double mult_6_7  =  12.0 * wigner_3j_6_m6_1_5;
const double mult_6_8  =  12.0 * wigner_3j_6_m6_2_4;
const double mult_6_9  =  12.0 * wigner_3j_6_m4_1_3;
const double mult_6_10 = -12.0 * wigner_3j_6_m5_1_4;
const double mult_6_11 = -12.0 * wigner_3j_6_m5_2_3;
const double mult_6_12 = -12.0 * wigner_3j_6_m3_1_2;
const double mult_6_13 =   6.0 * wigner_3j_6_m6_3_3;
const double mult_6_14 =   6.0 * wigner_3j_6_m4_2_2;
const double mult_6_15 =   6.0 * wigner_3j_6_m2_1_1;

const double norm_factor_l4_ = sqrt(4.0 * M_PI / 9.0);
const double norm_factor_l6_ = sqrt(4.0 * M_PI / 13.0);

//----------------------------------------END bop_calculator_constants.h ----------------------------------

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeBopAtom::ComputeBopAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
//  if (narg != 3) error->all(FLERR,"Illegal compute displace/atom command");

  peratom_flag = 1;
  size_peratom_cols = 4;
  create_attribute = 1;

  has_cutoff = false;
  cutsq = -1.0;
  grid_quality = 5;

  for (int a = 3; a < narg; ++a)
  {
    if (strcmp(arg[a], "cutoff") == 0)
    {
      if (a >= narg - 1)
      {
        error->all(FLERR, "Missing argument for 'cutoff' in compute bop/atom command");
      }
      a++;
      double cutoff = atof(arg[a]);
      has_cutoff = true;
      cutsq = cutoff * cutoff;
    }
    else if (strcmp(arg[a], "quality") == 0)
    {
      if (a >= narg - 1)
      {
        error->all(FLERR, "Missing argument for 'quality' in compute bop/atom command");
      }
      a++;
      grid_quality = atoi(arg[a]);
    }
  }
  if (!has_cutoff)
  {
    error->all(FLERR, "Missing required option 'cutoff' in compute bop/atom command");
  }


  if ( grid_quality == 0 )
    spherical_harmonics = new SphericalHarmonics(BOP_CALCULATOR_N_KNOTS_LEGENDRE_0, BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_0, memory);
  else if ( grid_quality == 1 )
    spherical_harmonics = new SphericalHarmonics(BOP_CALCULATOR_N_KNOTS_LEGENDRE_1, BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_1, memory);
  else if ( grid_quality == 2 )
    spherical_harmonics = new SphericalHarmonics(BOP_CALCULATOR_N_KNOTS_LEGENDRE_2, BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_2, memory);
  else if ( grid_quality == 3 )
    spherical_harmonics = new SphericalHarmonics(BOP_CALCULATOR_N_KNOTS_LEGENDRE_3, BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_3, memory);
  else if ( grid_quality == 4 )
    spherical_harmonics = new SphericalHarmonics(BOP_CALCULATOR_N_KNOTS_LEGENDRE_4, BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_4, memory);
  else if ( grid_quality == 5 )
    spherical_harmonics = new SphericalHarmonics(BOP_CALCULATOR_N_KNOTS_LEGENDRE_5, BOP_CALCULATOR_N_KNOTS_TRIGONOMETRIC_5, memory);
  else
  {
    error->all(FLERR, "Invalid 'quality' value in compute bop/atom command, supported values are 0-5");
  }

  nmax = 0;
  invariants = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBopAtom::~ComputeBopAtom()
{
  memory->destroy(invariants);
  delete spherical_harmonics;
}

/* ---------------------------------------------------------------------- */

void ComputeBopAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute bop/atom requires a pair style be defined");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"bop/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute bop/atom");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeBopAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeBopAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow local invariants array if necessary
  if (atom->nlocal > nmax) {
    memory->destroy(invariants);
    nmax = atom->nmax;
    memory->create(invariants,nmax,4,"bop/atom:bop");
    array_atom = invariants;
  }

  neighbor->build_one(list);

  double **x = atom->x;
  int *mask = atom->mask;

  for (int ii = 0; ii < list->inum; ii++) {
    int i = list->ilist[ii];

    for (int inv = 0; inv < 4; inv++)
    {
      invariants[i][inv] = 0.0;
    }

    if (mask[i] & groupbit) {
      double *xi = x[i];
      int* jlist = list->firstneigh[i];
      int jnum = list->numneigh[i];

      cmpx q_lm_vector[12];
      for (int c = 0; c < 12; ++c)
      {
        q_lm_vector[c].re = 0;
        q_lm_vector[c].im = 0;
      }

      // loop over list of all neighbors within force cutoff
      // distsq[] = distance sq to each
      // nearest[] = atom indices of neighbors

      int number_of_bonds = 0;
      double total_area = 0;
      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj] & NEIGHMASK;

        double* xj = x[j];
        double delx = x[i][0] - x[j][0];
        double dely = x[i][1] - x[j][1];
        double delz = x[i][2] - x[j][2];
        double rsq = delx*delx + dely*dely + delz*delz;
        if (has_cutoff && rsq < cutsq) {
          double r = sqrt(rsq);

          // The calculation of cos_theta and phi is copied from function bond_cart2sph in bond.h
          double cos_theta = delz / r;
          double phi;

          if ( delx == 0.0 )
          {
            if ( dely == 0.0 )
              phi = 0.0;
            else if ( dely > 0.0 )
              phi = HALF_OF_PI;
            else
              phi = THREE_HALFS_OF_PI;
          }
          else
          {
            phi = atan(dely / delx);
            if ( delx < 0.0 )
              phi += PI;
            else if ( dely < 0.0 )
              phi += TWO_PI;
          }

          double area = 1.0; // TODO: Calculate the area
          cmpx vec_Ylm_ij[12];
          if (grid_quality == 0)
          {
            spherical_harmonics->getSphHarmVecExact(cos_theta, phi, vec_Ylm_ij);
          }
          else
          {
            spherical_harmonics->getSphHarmVecFSI(cos_theta, phi, vec_Ylm_ij);
          }
          for (int lm = 0; lm < 12; lm++)
          {
            q_lm_vector[lm].re += area * vec_Ylm_ij[lm].re;
            q_lm_vector[lm].im += area * vec_Ylm_ij[lm].im;
          }
          total_area += area;
          number_of_bonds += 1;
        }
      }

      if (number_of_bonds == 0)
      {
        continue;
      }

      double modulus_squared[12];
      double q4_squared, q6_squared;
      double q4, q6;
      double w4, w6;
      double bracket_w4, bracket_w6;
      double inv_of_total_area = 1.0 / total_area;


// l = 4
// m = 0
      modulus_squared[0] = q_lm_vector[0].re * q_lm_vector[0].re;
      q4_squared  = modulus_squared[0];
// m = 1
      modulus_squared[1] = q_lm_vector[1].re * q_lm_vector[1].re + q_lm_vector[1].im * q_lm_vector[1].im;
      q4_squared += 2 * modulus_squared[1];
// m = 2
      modulus_squared[2] = q_lm_vector[2].re * q_lm_vector[2].re + q_lm_vector[2].im * q_lm_vector[2].im;
      q4_squared += 2 * modulus_squared[2];
// m = 3
      modulus_squared[3] = q_lm_vector[3].re * q_lm_vector[3].re + q_lm_vector[3].im * q_lm_vector[3].im;
      q4_squared += 2 * modulus_squared[3];
// m = 4
      modulus_squared[4] = q_lm_vector[4].re * q_lm_vector[4].re + q_lm_vector[4].im * q_lm_vector[4].im;
      q4_squared += 2 * modulus_squared[4];

// l = 6
// m = 0
      modulus_squared[5] = q_lm_vector[5].re * q_lm_vector[5].re;
      q6_squared  = modulus_squared[5];
// m = 1
      modulus_squared[6] = q_lm_vector[6].re * q_lm_vector[6].re + q_lm_vector[6].im * q_lm_vector[6].im;
      q6_squared += 2 * modulus_squared[6];
// m = 2
      modulus_squared[7] = q_lm_vector[7].re * q_lm_vector[7].re + q_lm_vector[7].im * q_lm_vector[7].im;
      q6_squared += 2 * modulus_squared[7];
// m = 3
      modulus_squared[8] = q_lm_vector[8].re * q_lm_vector[8].re + q_lm_vector[8].im * q_lm_vector[8].im;
      q6_squared += 2 * modulus_squared[8];
// m = 4
      modulus_squared[9] = q_lm_vector[9].re * q_lm_vector[9].re + q_lm_vector[9].im * q_lm_vector[9].im;
      q6_squared += 2 * modulus_squared[9];
// m = 5
      modulus_squared[10] = q_lm_vector[10].re * q_lm_vector[10].re + q_lm_vector[10].im * q_lm_vector[10].im;
      q6_squared += 2 * modulus_squared[10];
// m = 6
      modulus_squared[11] = q_lm_vector[11].re * q_lm_vector[11].re + q_lm_vector[11].im * q_lm_vector[11].im;
      q6_squared += 2 * modulus_squared[11];

// obliczamy q4 i q6
      q4 = sqrt(q4_squared);
      q6 = sqrt(q6_squared);

      invariants[i][0] = q4 * norm_factor_l4_ * inv_of_total_area;
      invariants[i][1] = q6 * norm_factor_l6_ * inv_of_total_area;

// obliczamy w4 i w6

// przyczynek do w4:
// Re(q_{4,0}) * [ + 1 x W3j(4 4 4 |  0 0 0) | q_{4, 0} |^2 +
//                 - 6 x W3j(4 4 4 | -1 0 1) | q_{4, 1} |^2 +
//                 + 6 x W3j(4 4 4 | -2 0 2) | q_{4, 2} |^2 +
//                 - 6 x W3j(4 4 4 | -3 0 3) | q_{4, 3} |^2 +
//                 + 6 x W3j(4 4 4 | -4 0 4) | q_{4, 4} |^2 ]

      bracket_w4 = mult_4_0 * modulus_squared[0]
                   + mult_4_1 * modulus_squared[1]
                   + mult_4_2 * modulus_squared[2]
                   + mult_4_3 * modulus_squared[3]
                   + mult_4_4 * modulus_squared[4];

// przyczynek do w6:
// Re(q_{6,0}) * [ + 1 x W3j(6 6 6 |  0 0 0) | q_{6, 0} |^2 +
//                 - 6 x W3j(6 6 6 | -1 0 1) | q_{6, 1} |^2 +
//                 + 6 x W3j(6 6 6 | -2 0 2) | q_{6, 2} |^2 +
//                 - 6 x W3j(6 6 6 | -3 0 3) | q_{6, 3} |^2 +
//                 + 6 x W3j(6 6 6 | -4 0 4) | q_{6, 4} |^2 +
//                 - 6 x W3j(6 6 6 | -5 0 5) | q_{6, 5} |^2 +
//                 + 6 x W3j(6 6 6 | -6 0 6) | q_{6, 6} |^2 ]

      bracket_w6 = mult_6_0 * modulus_squared[5]
                   + mult_6_1 * modulus_squared[6]
                   + mult_6_2 * modulus_squared[7]
                   + mult_6_3 * modulus_squared[8]
                   + mult_6_4 * modulus_squared[9]
                   + mult_6_5 * modulus_squared[10]
                   + mult_6_6 * modulus_squared[11];

      w4 = q_lm_vector[0].re * bracket_w4;
      w6 = q_lm_vector[5].re * bracket_w6;

// harmoniki sferyczne (indeksy)
// 0 - 4, 0   1 - 4, 1   2 - 4, 2   3 - 4, 3   4 - 4, 4
// 5 - 6, 0   6 - 6, 1   7 - 6, 2   8 - 6, 3   9 - 6, 4   10 - 6, 5   11 - 6, 6

// pozostale przyczynki do w4:
//  12 x W3j(4 4 4 | -4 1 3) Re( q_{4, 1} q_{4, 3} q_{4, 4}* )
      w4 += mult_4_5 * Cmpx_Re_Z1_Z2_Z3Star(q_lm_vector[1], q_lm_vector[3], q_lm_vector[4]);
// -12 x W3j(4 4 4 | -3 1 2) Re( q_{4, 1} q_{4, 2} q_{4, 3}* )
      w4 += mult_4_6 * Cmpx_Re_Z1_Z2_Z3Star(q_lm_vector[1], q_lm_vector[2], q_lm_vector[3]);
//   6 x W3j(4 4 4 | -4 2 2) Re( q_{4, 2} q_{4, 2} q_{4, 4}* )
      w4 += mult_4_7 * Cmpx_Re_Z1_Z1_Z2Star(q_lm_vector[2], q_lm_vector[4]);
//   6 x W3j(4 4 4 | -2 1 1) Re( q_{4, 1} q_{4, 1} q_{4, 2}* )
      w4 += mult_4_8 * Cmpx_Re_Z1_Z1_Z2Star(q_lm_vector[1], q_lm_vector[2]);

// pozostale przyczynki do w6:
//  12 x W3j(6 6 6 | -6 1 5) Re( q_{6, 1} q_{6, 5} q_{6, 6}* )
      w6 += mult_6_7  * Cmpx_Re_Z1_Z2_Z3Star(q_lm_vector[6], q_lm_vector[10], q_lm_vector[11]);
//  12 x W3j(6 6 6 | -6 2 4) Re( q_{6, 2} q_{6, 4} q_{6, 6}* )
      w6 += mult_6_8  * Cmpx_Re_Z1_Z2_Z3Star(q_lm_vector[7], q_lm_vector[9], q_lm_vector[11]);
//  12 x W3j(6 6 6 | -4 1 3) Re( q_{6, 1} q_{6, 3} q_{6, 4}* )
      w6 += mult_6_9  * Cmpx_Re_Z1_Z2_Z3Star(q_lm_vector[6], q_lm_vector[8], q_lm_vector[9]);
// -12 x W3j(6 6 6 | -5 1 4) Re( q_{6, 1} q_{6, 4} q_{6, 5}* )
      w6 += mult_6_10 * Cmpx_Re_Z1_Z2_Z3Star(q_lm_vector[6], q_lm_vector[9], q_lm_vector[10]);
// -12 x W3j(6 6 6 | -5 2 3) Re( q_{6, 2} q_{6, 3} q_{6, 5}* )
      w6 += mult_6_11 * Cmpx_Re_Z1_Z2_Z3Star(q_lm_vector[7], q_lm_vector[8], q_lm_vector[10]);
// -12 x W3j(6 6 6 | -3 1 2) Re( q_{6, 1} q_{6, 2} q_{6, 3}* )
      w6 += mult_6_12 * Cmpx_Re_Z1_Z2_Z3Star(q_lm_vector[6], q_lm_vector[7], q_lm_vector[8]);
//   6 x W3j(6 6 6 | -6 3 3) Re( q_{6, 3} q_{6, 3} q_{6, 6}* )
      w6 += mult_6_13 * Cmpx_Re_Z1_Z1_Z2Star(q_lm_vector[8], q_lm_vector[11]);
//   6 x W3j(6 6 6 | -4 2 2) Re( q_{6, 2} q_{6, 2} q_{6, 4}* )
      w6 += mult_6_14 * Cmpx_Re_Z1_Z1_Z2Star(q_lm_vector[7], q_lm_vector[9]);
//   6 x W3j(6 6 6 | -2 1 1) Re( q_{6, 1} q_{6, 1} q_{6, 2}* )
      w6 += mult_6_15 * Cmpx_Re_Z1_Z1_Z2Star(q_lm_vector[6], q_lm_vector[7]);

// obliczamy parametry W4 i W6
      invariants[i][2] = w4 / ( q4 * q4_squared );
      invariants[i][3] = w6 / ( q6 * q6_squared );
    }
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputeBopAtom::set_arrays(int i)
{
  double **xoriginal = fix->astore;
  double **x = atom->x;
  xoriginal[i][0] = x[i][0];
  xoriginal[i][1] = x[i][1];
  xoriginal[i][2] = x[i][2];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeBopAtom::memory_usage()
{
  double bytes = nmax*4 * sizeof(double);
  return bytes;
}
