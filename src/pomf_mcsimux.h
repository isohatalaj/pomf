/**
 * Extended MC simulations of the model equations. Instead of just
 * simulating the main macro variables of the model, these routines
 * also compute e.g. the output growth predicted by the model.
 *
 */

#ifndef POMF_MCSIMUEXT_H
#define POMF_MCSIMUEXT_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "util.h"
#include "spgrid.h"

/* #include "pomf_model.h" */


/**
 * Number of elements in the state array used in the simulations.
 */
#define POMF_XDIM 6

/**
 * Number of independent noise sources
 */
#define POMF_XWDIM 2

/**
 * Indices to the state array used in the simulations.
 */
#define POMF_XSIM_IA 0
#define POMF_XSIM_IDELTA 1
#define POMF_XSIM_IETA 2
#define POMF_XSIM_INTOT 3
#define POMF_XSIM_IY 4
#define POMF_XSIM_IYTOT 5

/**
 * Number of variables to record from a simulation per time point.
 * Currently: a_tilde, r, delta, Ntot, y, ytot, where Ntot, y, ytot
 * are scaled so that their initial values are all ones.
 */
#define POMF_XELEMS 6

/**
 * Indices to elements of the sample data arrays.
 */
#define POMF_XIATILDE 0
#define POMF_XIR 1
#define POMF_XIDELTA 2
#define POMF_XINTOT 3
#define POMF_XIY 4
#define POMF_XIYTOT 5

/* These definitions export the above constants into the simulation
 * library, so that they are accessible e.g. in Python code. */
extern const int pomf_xdim;
extern const int pomf_xwdim;
extern const int pomf_xsim_ia;
extern const int pomf_xsim_idelta;
extern const int pomf_xsim_ieta;
extern const int pomf_xsim_intot;
extern const int pomf_xsim_iy;
extern const int pomf_xsim_iytot;
extern const int pomf_xelems;
extern const int pomf_xiatilde;
extern const int pomf_xir;
extern const int pomf_xidelta;
extern const int pomf_xintot;
extern const int pomf_xiy;
extern const int pomf_xiytot;

/* REPRODUCED IN PYTHON CODE -- ANY CHANGES HERE MUST BE ALSO
 * IMMEDIATELY MIRRORED IN THE PYTHON SOURCES AS WELL. 
 *
 * Also, the first elements should match those of the structure
 * pomf_params_t, since at the moment, for a few function calls,
 * pomf_params_t and pomf_paramx_t are being used interchangably for
 * function calls that only depend on the main model parameters.
 * 
 **/
typedef struct {
  /* Reduced set of parameters. */
  double a_bar;
  double gamma;
  double sigma_c;
  double sigma_c_bar;
  double sigma_a;
  double kappa;
  double r_bar;
  double rho;
  double rho_bar;

  double sin_theta;
  double sin_theta_bar;

  double a_tilde_init;
  double r_init; 

  /* Parameters for simulating the extra variables. */
  double N0_over_y0;
  double Ntot0_over_ytot0;
  
} pomf_paramx_t;

/**
 * Caches frequently used constants that are derived from the static
 * parameters.
 */
typedef struct {
  double delta_equil_stdev;
  double gamma_delta;
  double s;
  double s_bar;
  double sigma_a_tilde;
  double kappa_prime;
  double sigma_c_inv;
  double sigma_c_bar_inv;
} pomf_pcachex_t;

void
pomf_mcsimux_cache(pomf_pcachex_t *pc,
		   const pomf_paramx_t *px);

typedef struct {
  gsl_rng *rng;
  size_t n_samples;
  double *samples;
  pomf_paramx_t px;
  pomf_pcachex_t pc;
} pomf_mcsimux_t;

double
pomf_r_starx_(const double at, const double eta,
	      const pomf_paramx_t *px,
	      const pomf_pcachex_t *pc);

double
pomf_r_starx(const double at, const double eta,
	     const pomf_paramx_t *px);

double
pomf_eta_starx_(const double at, const double r,
		const pomf_paramx_t *px,
		const pomf_pcachex_t *pc);

double
pomf_eta_starx(const double at, const double r,
	       const pomf_paramx_t *px);

double
pomf_delta_equil_stdevx(const pomf_paramx_t *px);

double
pomf_sigma_a_tilde_starx(const pomf_paramx_t *px);

pomf_mcsimux_t *
pomf_mcsimux_make(size_t n_samples);

void
pomf_mcsimux_free(pomf_mcsimux_t *self);
		  
int
pomf_mcsimux_init(pomf_mcsimux_t *self,
		  pomf_paramx_t *px);

int
pomf_mcsimux_init_wobs(pomf_mcsimux_t *self,
		       double r,
		       double y,
		       double ytot);

int
pomf_mcsimux_run(pomf_mcsimux_t *self,
		 double t_end,
		 double dt);

double *
pomf_mcsimux_get_samples(pomf_mcsimux_t *self);

#endif

