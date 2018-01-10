
#include <math.h>

#include "pomf_model.h"


double
pomf_r_star(const double at, const double eta,
	    const pomf_params_t *p)
{
  const double eta_prime = 1.0 - eta;
  const double s = 1.0/(p->sigma_c*p->sigma_c);
  const double s_bar = 1.0/(p->sigma_c_bar*p->sigma_c_bar);

  return
    ((at*s - 1)*eta + (p->r_bar*s_bar)*eta_prime) /
    (s*eta + s_bar*eta_prime);
}

double
pomf_sigma_a_tilde_star(const pomf_params_t *p)
{
  const double kappa = p->kappa;
  const double kappa_prime = sqrt(1 - kappa*kappa);
  
  const double t0 = p->gamma*p->sigma_c;
  const double t1 = kappa*p->sigma_a + p->gamma*p->sigma_c;
  const double t2 = p->sigma_a*kappa_prime;
  
  return -t0 + sqrt(t1*t1 + t2*t2); 
}

void
pomf_limitfunc(double bs[2*POMF_DIM], void *params)
{
  pomf_params_t *p = params;
  
  const double kappa = p->kappa;
  const double sigma_a_tilde = pomf_sigma_a_tilde_star(p);
  const double sat_over_sc   = sigma_a_tilde / p->sigma_c;

  const double a_mean = p->a_bar;
  const double a_revrate = p->gamma;
  const double a_stddev = p->sigma_a/sqrt(2*a_revrate);

  const double delta_mean = 0.0;  
  const double delta_revrate = (p->gamma + sat_over_sc);
  const double sigma_delta
    = sqrt(-2*p->sigma_a*kappa*sigma_a_tilde + sigma_a_tilde*sigma_a_tilde + p->sigma_a*p->sigma_a);
  const double delta_stddev = sigma_delta/sqrt(2*delta_revrate);

  bs[2*POMF_I_A + 0] = a_mean - 3*a_stddev;
  bs[2*POMF_I_A + 1] = a_mean + 3*a_stddev;
  bs[2*POMF_I_DELTA + 0] = delta_mean - 3*delta_stddev;
  bs[2*POMF_I_DELTA + 1] = delta_mean + 3*delta_stddev;
  bs[2*POMF_I_ETA + 0] = 0.001;
  bs[2*POMF_I_ETA + 1] = 0.999;
}

double
eta2xi(double eta)
{
  return log(eta/(1 - eta));
}

double
xi2eta(double xi)
{
  return 1.0 / (1 + exp(-xi));
}

/* xi -> eta with jacobian and friends */
int
pomf_Eta(const double *x, double *y,
	 double *jacobian,
	 double *inv_jacobian,
	 double *det_jacobian,
	 double *hessians,
	 void *psi_params)
{
  const double xi = x[POMF_I_XI];
  const double a = x[POMF_I_A];
  const double delta = x[POMF_I_DELTA];

  const double zeta = exp(-xi);
  const double eta = 1.0 / (1 + zeta);

  y[POMF_I_ETA] = eta;
  y[POMF_I_A] = a;
  y[POMF_I_DELTA] = delta;

  *det_jacobian = zeta*eta*eta;
    
  jacobian[POMF_I_ETA*POMF_DIM + POMF_I_ETA  ] = *det_jacobian;
  jacobian[POMF_I_ETA*POMF_DIM + POMF_I_A    ] = 0.0;
  jacobian[POMF_I_ETA*POMF_DIM + POMF_I_DELTA] = 0.0;

  jacobian[POMF_I_A*POMF_DIM + POMF_I_ETA  ] = 0.0;
  jacobian[POMF_I_A*POMF_DIM + POMF_I_A    ] = 1.0;
  jacobian[POMF_I_A*POMF_DIM + POMF_I_DELTA] = 0.0;

  jacobian[POMF_I_DELTA*POMF_DIM + POMF_I_ETA  ] = 0.0;
  jacobian[POMF_I_DELTA*POMF_DIM + POMF_I_A    ] = 0.0;
  jacobian[POMF_I_DELTA*POMF_DIM + POMF_I_DELTA] = 1.0;

  inv_jacobian[POMF_I_ETA*POMF_DIM + POMF_I_ETA  ] = 1.0 / *det_jacobian;
  inv_jacobian[POMF_I_ETA*POMF_DIM + POMF_I_A    ] = 0.0;
  inv_jacobian[POMF_I_ETA*POMF_DIM + POMF_I_DELTA] = 0.0;

  inv_jacobian[POMF_I_A*POMF_DIM + POMF_I_ETA  ] = 0.0;
  inv_jacobian[POMF_I_A*POMF_DIM + POMF_I_A    ] = 1.0;
  inv_jacobian[POMF_I_A*POMF_DIM + POMF_I_DELTA] = 0.0;

  inv_jacobian[POMF_I_DELTA*POMF_DIM + POMF_I_ETA  ] = 0.0;
  inv_jacobian[POMF_I_DELTA*POMF_DIM + POMF_I_A    ] = 0.0;
  inv_jacobian[POMF_I_DELTA*POMF_DIM + POMF_I_DELTA] = 1.0;
  
  int i;
  for (i = 0; i < POMF_DIM*POMF_DIM2; ++i) hessians[i] = 0.0;

  double *hess_eta = hessians + POMF_I_ETA*POMF_DIM2;

  hess_eta[POMF_I_ETA*POMF_DIM + POMF_I_ETA] = *det_jacobian * (zeta - 1) * eta;

  return 0;
}

void
pomf_limitfunc_alt(double bs[POMF_DIM], void *params)
{
  pomf_limitfunc(bs, params);
  bs[2*POMF_I_XI + 0] = eta2xi(bs[2*POMF_I_ETA + 0]);
  bs[2*POMF_I_XI + 1] = eta2xi(bs[2*POMF_I_ETA + 1]);
}
  
int
pomf_sdefunc(const double x[POMF_DIM], 
	     double mu_out[POMF_DIM],
	     double sigma_out[POMF_DIM*POMF_WDIM],
	     void *params)
{
  pomf_params_t *p = params;

  /* Some stuff here could be precomputed, but let's just keep it
     simple for now. */
  
  const double kappa = p->kappa;
  const double kappa_prime = sqrt(1 - kappa*kappa);

  const double a = x[POMF_I_A];
  const double delta = x[POMF_I_DELTA]; 
  const double a_tilde = a - delta;
  const double eta = x[POMF_I_ETA];
  const double eta_prime = 1 - eta;
  
  const double r = pomf_r_star(a_tilde, eta, p);
  const double h_h_prime = eta*(1.0 - eta);

  const double sh_a_c        =        (a - r) / p->sigma_c;
  const double sh_a_tilde_c  =  (a_tilde - r) / p->sigma_c;
  const double sh_r_c_bar    = (r - p->r_bar) / p->sigma_c_bar;
  const double r_diff        =  r - p->r_bar;
  const double rho_diff      = p->rho - p->rho_bar;
  const double sigma_a_tilde = pomf_sigma_a_tilde_star(p);
  const double sat_over_sc   = sigma_a_tilde / p->sigma_c;
  
  double mu_h = r_diff - rho_diff + sh_a_c*sh_a_tilde_c - sh_r_c_bar*sh_r_c_bar
    - (sh_a_tilde_c - sh_r_c_bar)*(sh_a_tilde_c*eta + sh_r_c_bar*eta_prime);
  double sigma_h = sh_a_tilde_c - sh_r_c_bar;

  mu_out[POMF_I_A] = p->gamma*(p->a_bar - a);
  mu_out[POMF_I_DELTA] = -(p->gamma + sat_over_sc)*delta;
  mu_out[POMF_I_ETA] = h_h_prime*mu_h;
  
  sigma_out[POMF_I_A*POMF_WDIM + 0] = p->sigma_a*kappa;
  sigma_out[POMF_I_A*POMF_WDIM + 1] = p->sigma_a*kappa_prime;
  sigma_out[POMF_I_A*POMF_WDIM + 2] = 0.0;
  
  sigma_out[POMF_I_DELTA*POMF_WDIM + 0] = p->sigma_a*kappa - sigma_a_tilde;
  sigma_out[POMF_I_DELTA*POMF_WDIM + 1] = p->sigma_a*kappa_prime;
  sigma_out[POMF_I_DELTA*POMF_WDIM + 2] = 0.0;  
  
  sigma_out[POMF_I_ETA*POMF_WDIM + 0] = h_h_prime*sigma_h;
  sigma_out[POMF_I_ETA*POMF_WDIM + 1] = 0.0;
  sigma_out[POMF_I_ETA*POMF_WDIM + 2] = 0.0;
  
  return 0;
}

int
pomf_sdefunc_alt(const double x[POMF_DIM], 
		 double mu_out[POMF_DIM],
		 double sigma_out[POMF_DIM*POMF_WDIM],
		 void *params)
{
  pomf_params_t *p = params;

  /* Some stuff here could be precomputed, but let's just keep it
     simple for now. */
  
  const double kappa = p->kappa;
  const double kappa_prime = sqrt(1 - kappa*kappa);

  const double a = x[POMF_I_A];
  const double delta = x[POMF_I_DELTA]; 
  const double a_tilde = a - delta;
  const double eta = xi2eta(x[POMF_I_XI]);
  const double eta_prime = 1 - eta;
  
  const double r = pomf_r_star(a_tilde, eta, p);

  const double sh_a_c        =        (a - r) / p->sigma_c;
  const double sh_a_tilde_c  =  (a_tilde - r) / p->sigma_c;
  const double sh_r_c_bar    = (r - p->r_bar) / p->sigma_c_bar;
  const double r_diff        =  r - p->r_bar;
  const double rho_diff      = p->rho - p->rho_bar;
  const double sigma_a_tilde = pomf_sigma_a_tilde_star(p);
  const double sat_over_sc   = sigma_a_tilde / p->sigma_c;
  
  double mu_h = r_diff - rho_diff + sh_a_c*sh_a_tilde_c - sh_r_c_bar*sh_r_c_bar
    - (sh_a_tilde_c - sh_r_c_bar)*(sh_a_tilde_c*eta + sh_r_c_bar*eta_prime);
  double sigma_h = sh_a_tilde_c - sh_r_c_bar;

  mu_out[POMF_I_A] = p->gamma*(p->a_bar - a);
  mu_out[POMF_I_DELTA] = -(p->gamma + sat_over_sc)*delta;
  mu_out[POMF_I_XI] = mu_h - 0.5*(1 - 2*eta)*sigma_h*sigma_h;
  
  sigma_out[POMF_I_A*POMF_WDIM + 0] = p->sigma_a*kappa;
  sigma_out[POMF_I_A*POMF_WDIM + 1] = p->sigma_a*kappa_prime;
  sigma_out[POMF_I_A*POMF_WDIM + 2] = 0.0;
  
  sigma_out[POMF_I_DELTA*POMF_WDIM + 0] = p->sigma_a*kappa - sigma_a_tilde;
  sigma_out[POMF_I_DELTA*POMF_WDIM + 1] = p->sigma_a*kappa_prime;
  sigma_out[POMF_I_DELTA*POMF_WDIM + 2] = 0.0;  
  
  sigma_out[POMF_I_ETA*POMF_WDIM + 0] = sigma_h;
  sigma_out[POMF_I_ETA*POMF_WDIM + 1] = 0.0;
  sigma_out[POMF_I_ETA*POMF_WDIM + 2] = 0.0;
  
  return 0;
}

