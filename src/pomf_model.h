
#ifndef POMF_MODEL_H
#define POMF_MODEL_H

#define POMF_DIM 3
#define POMF_DIM2 (POMF_DIM*POMF_DIM)
#define POMF_WDIM 3

#define POMF_I_A 0
#define POMF_I_DELTA 1
#define POMF_I_ETA 2
#define POMF_I_XI 2

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
} pomf_params_t;

double
pomf_r_star(const double at, const double eta,
	    const pomf_params_t *p);

double
pomf_sigma_a_tilde_star(const pomf_params_t *p);

void
pomf_limitfunc(double bs[2*POMF_DIM], void *params);

void
pomf_limitfunc_alt(double bs[2*POMF_DIM], void *params);

int
pomf_Eta(const double *x, double *y,
	 double *jacobian,
	 double *inv_jacobian,
	 double *det_jacobian,
	 double *hessians,
	 void *params);

int
pomf_sdefunc(const double x[POMF_DIM], 
	     double mu_out[POMF_DIM],
	     double sigma_out[POMF_DIM*POMF_WDIM],
	     void *params);

int
pomf_sdefunc_alt(const double x[POMF_DIM], 
		 double mu_out[POMF_DIM],
		 double sigma_out[POMF_DIM*POMF_WDIM],
		 void *params);

#endif
