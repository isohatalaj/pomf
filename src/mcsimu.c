
#include <stdio.h>
#include <math.h>

#include "pomf_mcsimu.h"
#include "pomf_model.h"

int
pomf_mcsimu_add_data(spgrid_t *grid,
		     const double *x,
		     const double *dx,
		     double r, double dt)
{
  int status = OK;
  double y[POMF_DIM + 2];
  double hx[POMF_DIM];

  spgrid_volelem(grid, x, hx);
  
  y[0] = dx[0] / (hx[1] * hx[2] * dt);
  y[1] = dx[1] / (hx[0] * hx[2] * dt);
  y[2] = dx[2] / (hx[0] * hx[1] * dt);
  y[3] = r;
  y[4] = 1.0;
  
  CHECK( spgrid_add_data(grid, x, y), FAILED );

 exit:

  return status;
}

void
pomf_mcsimu_datainit(double *states,
		     int n_copies,
		     pomf_params_t *p)
{
  size_t i_copy;

  for (i_copy = 0; i_copy < n_copies; ++i_copy)
    {
      states[POMF_DIM*i_copy + POMF_I_A] = p->a_bar;
      states[POMF_DIM*i_copy + POMF_I_DELTA] = 0.0;
      states[POMF_DIM*i_copy + POMF_I_ETA] = 0.25;
    }
}

int
pomf_mcsimu_step(gsl_rng *rng,
		 double *states,
		 const size_t n_copies,
		 const size_t n_steps,
		 const size_t n_skip,
		 const double dt,
		 pomf_params_t *p,
		 spgrid_t *grid)
{
  int status = OK;
  
  const double sqrt_dt = sqrt(dt);

  size_t i_step;
  size_t i_copy;

  for (i_step = 0; i_step < n_steps; ++i_step)
    {
      for (i_copy = 0; i_copy < n_copies; ++i_copy)
	{
	  double *x = &states[i_copy*POMF_DIM];
	  double dw[POMF_WDIM];
	  int i, j;

	  for (i = 0; i < POMF_WDIM; ++i)
	    dw[i] = gsl_ran_gaussian_ziggurat(rng, sqrt_dt);

	  double mu[POMF_DIM];
	  double sigma[POMF_DIM*POMF_WDIM];
	  pomf_sdefunc(x, mu, sigma, p);

	  double dx[POMF_DIM];
	  
	  for (i = 0; i < POMF_DIM; ++i)
	    {
	      dx[i] = mu[i]*dt;
	      for (j = 0; j < POMF_WDIM; ++j) dx[i] += sigma[i*POMF_WDIM + j]*dw[j];
	    }

	  /* Enforce [0,1] domain for eta by reflection, if finite time step
	   * attempts to jump out of the allowed range. */
	  if (x[POMF_I_ETA] + dx[POMF_I_ETA] < 0.0)
	    dx[POMF_I_ETA] = -2*x[POMF_I_ETA] - dx[POMF_I_ETA];

	  if (x[POMF_I_ETA] + dx[POMF_I_ETA] > 1.0)
	    dx[POMF_I_ETA] = 2*(1 - x[POMF_I_ETA]) - dx[POMF_I_ETA];

	  double r;
	  
	  if (grid && i_step >= n_skip)
	    {
	      r = pomf_r_star(x[POMF_I_A] - x[POMF_I_DELTA], x[POMF_I_ETA], p);
	      CHECK( pomf_mcsimu_add_data(grid, x, dx, r, dt), FAILED );
	    }

	  for (i = 0; i < POMF_DIM; ++i) x[i] += dx[i];

	  if (grid && i_step >= n_skip)
	    {
	      r = pomf_r_star(x[POMF_I_A] - x[POMF_I_DELTA], x[POMF_I_ETA], p);
	      CHECK( pomf_mcsimu_add_data(grid, x, dx, r, dt), FAILED );
	    }
	}
    }

 exit:

  return status;
}
