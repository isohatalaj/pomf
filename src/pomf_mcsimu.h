
#ifndef POMF_MCSIMU_H
#define POMF_MCSIMU_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "util.h"
#include "spgrid.h"

#include "pomf_model.h"

#define POMF_MCSIMU_NDATA (POMF_DIM + 2)

int
pomf_mcsimu_add_data(spgrid_t *grid,
		     const double *x, const double *dx,
		     double r, double dt);

void
pomf_mcsimu_datainit(double *states,
		     int n_copies,
		     pomf_params_t *p);

int
pomf_mcsimu_step(gsl_rng *rng,
		 double *states,
		 const size_t n_copies,
		 const size_t n_steps,
		 const size_t n_skip,
		 const double dt,
		 pomf_params_t *p,
		 spgrid_t *grid);


typedef struct {
  gsl_rng *rng;
  size_t n_samples;
  double *samples;
} pomf_mcsimu_t;

pomf_mcsimu_t *
pomf_mcsimu_make(size_t n_samples);

void
pomf_mcsimu_free(pomf_mcsimu_t *self);

double *
pomf_mcsimu_get_samples(pomf_mcsimu_t *self);

/**
 * Initialize `pomf_mcsimu_t` object for a simulation run with given
 * initial values of `a_tilde`, `r`, and an initial, Gaussian
 * distribution for the `delta` variable. If `delta_init_stdev` is
 * negative, then `delta` component of the state is not reset. Number
 * of samples to generate is stored in the `pomf_mcsimu_t` structure,
 * and is fixed at the time the object is created. Returns zero for
 * success.
 */
int
pomf_mcsimu_init(pomf_mcsimu_t *self,
		 double a_tilde_init,
		 double r_init,
		 double delta_init_mean,
		 double delta_init_stdev,
		 pomf_params_t *p);

/**
 * Run a simulation of the model SDEs with initial state set by
 * `pomf_mcsimu_init` or left by a previous call to `pomf_mcsimu_run`.
 * Simulation end time and time step length are `t_end` and `dt`
 * respectively. Values of the process at time `t_end` are then stored
 * in the array `samples` in `pomf_mcsimu_t` structure. Returns zero
 * for success.
 */
int
pomf_mcsimu_run(pomf_mcsimu_t *self,
		double t_end,
		double dt,
		pomf_params_t *p);


#endif
