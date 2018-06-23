
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
  double y[POMF_MCSIMU_NDATA];
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

static int
basicstep(gsl_rng *rng,
	  const double *x,
	  double *dx,
	  double dt,
	  double sqrt_dt,
	  pomf_params_t *p)
{
  int status = OK;
  
  double dw[POMF_WDIM];
  int i, j;

  if (x[POMF_I_ETA] < 0.0 || x[POMF_I_ETA] > 1.0)
    {
      ERRMSG("eta variable out of bounds at start of euler step\n");
      ERRMSG("eta = %lf, a = %lf, delta = %lf\n",
	     x[POMF_I_ETA], x[POMF_I_A], x[POMF_I_DELTA]);
	     
      FAILWITH(FAILED);
    }
  
  for (i = 0; i < POMF_WDIM; ++i)
    dw[i] = gsl_ran_gaussian_ziggurat(rng, sqrt_dt);
  
  double mu[POMF_DIM];
  double sigma[POMF_DIM*POMF_WDIM];

  CHECK_( pomf_sdefunc(x, mu, sigma, p) );

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

  if (x[POMF_I_ETA] + dx[POMF_I_ETA] < 0.0 || x[POMF_I_ETA] + dx[POMF_I_ETA] > 1.0)
    {
      ERRMSG("Euler step out of bounds even after reflection, too large dt?\n");
      FAILWITH(FAILED);
    }

 exit:

  return status;
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
	  double dx[POMF_DIM];
	  int i;
	  
	  CHECK_( basicstep(rng, x, dx, dt, sqrt_dt, p) );

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

pomf_mcsimu_t *
pomf_mcsimu_make(size_t n_samples)
{
  int status = OK;
  pomf_mcsimu_t *self = NULL;

  CHECK_NULL( self = malloc(sizeof(*self)), OUT_OF_MEM );

  self->rng = NULL;
  self->n_samples = n_samples;
  self->samples = NULL;

  CHECK_NULL( self->rng = gsl_rng_alloc(gsl_rng_default), OUT_OF_MEM );
  CHECK_NULL( self->samples = malloc(self->n_samples*POMF_NOBS*sizeof(*self->samples)), OUT_OF_MEM );
  
 exit:
  if (status)
    {
      Xfree_(self, pomf_mcsimu_free);
      self = NULL;
    }

  return self;
}

void
pomf_mcsimu_free(pomf_mcsimu_t *self)
{
  if (self)
    {
      Xfree_(self->rng, gsl_rng_free);
      Xfree(self->samples);
      Xfree(self);
    }
}

int
pomf_mcsimu_init(pomf_mcsimu_t *self,
		 double a_tilde_init,
		 double r_init,
		 double delta_init_mean,
		 double delta_init_stdev,
		 pomf_params_t *p)
{
  int status = OK;
  int i_sample;

  for (i_sample = 0; i_sample < self->n_samples; ++i_sample)
    {
      if (delta_init_stdev >= 0.0)
	{
	  double delta_init = delta_init_mean
	    + gsl_ran_gaussian_ziggurat(self->rng, delta_init_stdev);

	  self->samples[POMF_NOBS*i_sample + POMF_IOBS_ATILDE] = a_tilde_init;
	  self->samples[POMF_NOBS*i_sample + POMF_IOBS_R] = r_init;
	  self->samples[POMF_NOBS*i_sample + POMF_IOBS_DELTA] = delta_init;
	}
      else
	{
	  self->samples[POMF_NOBS*i_sample + POMF_IOBS_ATILDE] = a_tilde_init;
	  self->samples[POMF_NOBS*i_sample + POMF_IOBS_R] = r_init;
	}
    }

  return status;
}

int
pomf_mcsimu_run(pomf_mcsimu_t *self,
		double t_end,
		double dt,
		pomf_params_t *p)
{
  int status = OK;
  int i_sample;

  for (i_sample = 0; i_sample < self->n_samples; ++i_sample)
    {
      const double sqrt_dt = sqrt(dt);
      
      double x[POMF_DIM], dx[POMF_DIM];
      double delta_init = self->samples[POMF_NOBS*i_sample + POMF_IOBS_DELTA];
      double a_tilde_init = self->samples[POMF_NOBS*i_sample + POMF_IOBS_ATILDE];
      double r_init = self->samples[POMF_NOBS*i_sample + POMF_IOBS_R];
      double a_init = a_tilde_init + delta_init;
      double eta_init = pomf_eta_star(a_tilde_init, r_init, p);

      x[POMF_I_ETA] = eta_init;
      x[POMF_I_A] = a_init;
      x[POMF_I_DELTA] = delta_init;

      double t = 0.0;

      while (t < t_end)
	{
	  if (t + dt <= t_end)
	    {
	      CHECK_( basicstep(self->rng, x, dx, dt, sqrt_dt, p) );
	      t += dt;
	    }
	  else
	    {
	      CHECK_( basicstep(self->rng, x, dx, t_end - t, sqrt(t_end - t), p) );
	      t = t_end;
	    }

	  int i;
	  for (i = 0; i < POMF_DIM; ++i) x[i] += dx[i];

	}

      double a_tilde = x[POMF_I_A] - x[POMF_I_DELTA];
      double r = pomf_r_star(a_tilde, x[POMF_I_ETA], p);
      double delta = x[POMF_I_DELTA];

      self->samples[POMF_NOBS*i_sample + POMF_IOBS_ATILDE] = a_tilde;
      self->samples[POMF_NOBS*i_sample + POMF_IOBS_R] = r;
      self->samples[POMF_NOBS*i_sample + POMF_IOBS_DELTA] = delta;
    }

 exit:

  return status;
}

double *
pomf_mcsimu_get_samples(pomf_mcsimu_t *self)
{
  return self->samples;
}
