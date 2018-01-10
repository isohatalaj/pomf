
#ifndef POMF_SMSIMU_H
#define POMF_SMSIMU_H

#include <cs.h>

#include "util.h"
#include "rgrid.h"
#include "smaux.h"

/**
 * Construct KFE left-hand side operator matrix, allocated into
 * `KFE_out`, the probability current operators, allocated into
 * `Js_out`, and the right-hand side vector, allocated into
 * `b_out`. Imposing reflecting boundary conditions on the boundaries
 * and normalisation to unity.
 */
int
pomf_smkfe(const rgrid_t *grid,
	   int dim, int wdim,
	   int (*sdefunc)(const double *x,
			  double *mu, double *sigma,
			  void *params),
	   void *params,
	   cs **KFE_out,
	   cs ***Js_out,
	   double **b_out);

/**
 * Change of variables for density and probability currents. Original
 * variables are x, new variables y, function psi maps x's to y's.
 * Supplied psi function does the mapping, plus evaluates the
 * jacobian, its inverse and determinant, and the hessians for each
 * component of the transform.
 */
int
pomf_covfj(rgrid_t *grid,
	   const double *fx, 
	   double **y_out, double *fy_out, double **Jy_out,
	   int (*psi)(const double *x, double *y,
		      double *jacobian,
		      double *inverse_jacobian,
		      double *det_jacobian,
		      double *hessians,
		      void *psi_params),
	   void *psi_params,
	   int wdim,
	   int (*musigma)(const double *x,
			  double *mu, double *sigma,
			  void *musigma_params),
	   void *musigma_params);

#endif
