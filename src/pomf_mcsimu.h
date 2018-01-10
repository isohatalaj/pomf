
#ifndef POMF_MCSIMU_H
#define POMF_MCSIMU_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "util.h"
#include "spgrid.h"

#include "pomf_model.h"

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

#endif
