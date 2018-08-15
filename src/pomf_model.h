
#ifndef POMF_MODEL_H
#define POMF_MODEL_H

#define POMF_DIM 3
#define POMF_DIM2 (POMF_DIM*POMF_DIM)
#define POMF_WDIM 3
#define POMF_NOBS 3

#define POMF_I_A 0
#define POMF_I_DELTA 1
#define POMF_I_ETA 2
#define POMF_I_XI 2

/* Indices of observable quantities in (observable) state vectors. */
#define POMF_IOBS_ATILDE 0
#define POMF_IOBS_R 1
#define POMF_IOBS_DELTA 2

/* Retrieve constants POMF_DIM, POMF_WDIM, etc. Used by the Python
   interface. */
extern const int pomf_dim;
extern const int pomf_wdim;
extern const int pomf_nobs;
extern const int pomf_iobs_atilde;
extern const int pomf_iobs_r;
extern const int pomf_iobs_delta;

/**
 * Model parameters. This struct is mirrored in the Python code.
 */
typedef struct {
  double a_bar;
  double gamma;
  double sigma_c;
  double sigma_c_bar;
  double sigma_a;
  double kappa;
  double r_bar;
  double rho;
  double rho_bar;

  /* For diagonalizing transform, used only in the third alt. change
   * of variables version and should be set using the appropriate
   * function.
   */
  double x_star[POMF_DIM];
  double U[POMF_DIM2];
} pomf_params_t;

double
pomf_r_star(const double at, const double eta,
	    const pomf_params_t *p);

double
pomf_eta_star(const double at, const double r,
	      const pomf_params_t *p);

double
pomf_sigma_a_tilde_star(const pomf_params_t *p);


/**
 * Compute the standard deviation of the equilibrium, unconditional
 * distribution of `delta := a - a_tilde`. The mean of the
 * distribution is zero.
 */
double
pomf_delta_equil_stdev(const pomf_params_t *p);

void
pomf_limitfunc(double bs[2*POMF_DIM], void *params);

int
pomf_limitfunc2(double bs[2*POMF_DIM], void *params);

void
pomf_limitfunc_alt(double bs[2*POMF_DIM], void *params);

int
pomf_limitfunc_alt2(double bs[2*POMF_DIM], void *params);

int
pomf_limitfunc_alt3(double bs[2*POMF_DIM], void *params);

int
pomf_Eta(const double *x, double *y,
	 double *jacobian,
	 double *inv_jacobian,
	 double *det_jacobian,
	 double *hessians,
	 void *params);

/**
 * Compute the drift and diffusion parameters for total net worth
 * Ntot, firm output y, and total output ytot.
 */
int
pomf_ext_musigma(const double a, const double eta, const double delta,
		 const double Ntot,
		 const double N0_over_y0,
		 const double Ntot0_over_ytot0,
		 double *mu_Ntot, double *sigma_Ntot,
		 double *mu_y, double *sigma_y,
		 double *mu_ytot, double *sigma_ytot,
		 const pomf_params_t *p);

int
pomf_sdefunc(const double x[POMF_DIM], 
	     double mu_out[POMF_DIM],
	     double sigma_out[POMF_DIM*POMF_WDIM],
	     void *params);

int
pomf_sdefunc2(const double x[POMF_DIM], 
	      double mu_out[POMF_DIM],
	      double sigma_out[POMF_DIM*POMF_WDIM],
	      double sigma_diff[POMF_DIM*POMF_DIM*POMF_WDIM],
	      void *params);

int
pomf_sdefunc_alt(const double x[POMF_DIM], 
		 double mu_out[POMF_DIM],
		 double sigma_out[POMF_DIM*POMF_WDIM],
		 void *params);

int
pomf_sdefunc_alt2(const double x[POMF_DIM], 
		  double mu_out[POMF_DIM],
		  double sigma_out[POMF_DIM*POMF_WDIM],
		  double sigma_diff[POMF_DIM*POMF_DIM*POMF_WDIM],
		  void *params);

/**
 * alt2, but with a linear change of variables x = Uy + x*
 */
int
pomf_sdefunc_alt3(const double x[POMF_DIM], 
		  double mu_out[POMF_DIM],
		  double sigma_out[POMF_DIM*POMF_WDIM],
		  double sigma_diff[POMF_DIM*POMF_DIM*POMF_WDIM],
		  void *params);

int
pomf_setup_alt3(const double x_star[POMF_DIM],
		pomf_params_t *p);


#endif
