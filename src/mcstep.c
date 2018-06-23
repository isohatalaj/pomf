
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "spgrid.h"

#include "pomf_model.h"
#include "pomf_mcsimu.h"

#define FDOUBLE " %25.15lf"
#define FPHOLDR " %25s"


int
main(void)
{
  int status;
  gsl_rng *rng = NULL;
  double *states = NULL;

  const int dim = POMF_DIM;
  /* const int wdim = POMF_WDIM; */
  const int ndata = POMF_MCSIMU_NDATA; 

  gsl_rng_env_setup();

  rng = gsl_rng_alloc(gsl_rng_default);
  if (rng == NULL) { status = GSL_EFAILED; goto exit; }

  /* Setup model parameters */
  pomf_params_t p;

  /* Parametrizing TFP process via its equilibrium distribution mean,
     stddev, and rate of mean reversal. */
  /* const double a_mean = 0.10; */
  /* const double a_stddev = 0.05; */
  /* const double a_revrate = 1.0; */
  /* p.a_bar       = a_mean; */
  /* p.gamma       = a_revrate; */
  /* p.sigma_a     = sqrt(2*a_revrate)*a_stddev; */
  /* p.sigma_c     = 0.15; */
  /* p.sigma_c_bar = 0.15; */
  /* p.kappa       = 0.000; */
  /* p.r_bar       = 0.015; */
  /* p.rho         = 0.070;   /\* Note: Model depends only on difference rho - rho_bar *\/ */
  /* p.rho_bar     = 0.050; */

  /* const double a_mean = 0.070; */
  /* const double a_stddev = 0.035; */
  /* const double a_revrate = 0.05; */
  /* p.a_bar       = a_mean; */
  /* p.gamma       = a_revrate; */
  /* p.sigma_a     = sqrt(2*a_revrate)*a_stddev; */
  /* p.sigma_c     = 0.100; */
  /* p.sigma_c_bar = 0.100; */
  /* p.kappa       = -0.707; */
  /* p.r_bar       = 0.010; */
  /* p.rho         = 0.070;   /\* Note: Model depends only on difference rho - rho_bar *\/ */
  /* p.rho_bar     = 0.050; */

  const double a_mean = 0.070;
  const double a_stddev = 0.035;
  const double a_revrate = 0.5;
  p.a_bar       = a_mean;
  p.gamma       = a_revrate;
  p.sigma_a     = sqrt(2*a_revrate)*a_stddev;
  p.sigma_c     = 0.025;
  p.sigma_c_bar = 0.025;
  p.kappa       = -0.707;
  p.r_bar       = 0.010;
  p.rho         = 0.070;   /* Note: Model depends only on difference rho - rho_bar */
  p.rho_bar     = 0.050;
  
  double range_est[2*POMF_DIM];
  const double sigma_a_tilde = pomf_sigma_a_tilde_star(&p);
  pomf_limitfunc(range_est, &p);

  fprintf(stderr, "## sigma_a_tilde = %lf\n", sigma_a_tilde);

  /* Sparse grid data structure */
  spgrid_t grid;
  double h[POMF_DIM];
  h[0] = (range_est[2*0 + 1] - range_est[2*0 + 0]) / 50;
  h[1] = (range_est[2*1 + 1] - range_est[2*1 + 0]) / 50;
  h[2] = (range_est[2*2 + 1] - range_est[2*2 + 0]) / 50;
  
  spgrid_init(&grid, dim, ndata, h);
  spgrid_set_bound(&grid, POMF_I_ETA, 0.0, 1.0);

  /* Simulation parameters */
  const double t_end = 1.0;
  const size_t n_copies = 100000;
  const size_t n_steps = 10;
  const size_t n_skip = 9;
  const double dt = t_end / n_steps;

  /* Simulation data */
  states = malloc(POMF_DIM*n_copies*sizeof(*states));
  if (states == NULL) { status = GSL_ENOMEM; goto exit; }
  
  /* Run simulation */
  pomf_mcsimu_datainit(states, n_copies, &p);

  status = pomf_mcsimu_step(rng, states, n_copies,
			    n_steps, n_skip, dt, &p, &grid);
  if (status) { goto exit; }

  /* Take marginals */
  int i;
  for (i = 0; i < dim; ++i)
    {
      spgrid_t mgrid;
      rgrid_t rgrid;
      double **data;
      
      spgrid_marginal(&grid, &mgrid, i);
      spgrid_to_rgrid(&mgrid, &rgrid, &data);

      rgrid_print_data(&rgrid, NULL, ndata + 1, data, NULL);
      printf("\n\n");

      spgrid_destroy(&mgrid);
      rgrid_destroy(&rgrid);
      Xfree_many(data, ndata + 1, free);
    }


 exit:
  /* Done, exit. */

  spgrid_destroy(&grid);

  if (states != NULL) free(states);
  if (rng != NULL) gsl_rng_free(rng);

  if (status) GSL_ERROR("Error encountered", status);
  return GSL_SUCCESS;
}


