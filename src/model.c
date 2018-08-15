
#include <math.h>
#include <string.h>

#include "util.h"
#include "pomf_model.h"

const int pomf_dim = POMF_DIM;
const int pomf_wdim = POMF_WDIM;
const int pomf_nobs = POMF_NOBS;
const int pomf_iobs_atilde = POMF_IOBS_ATILDE;
const int pomf_iobs_r = POMF_IOBS_R;
const int pomf_iobs_delta = POMF_IOBS_DELTA;

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
pomf_eta_star(const double at, const double r,
	      const pomf_params_t *p)
{
  const double s = 1.0/(p->sigma_c*p->sigma_c);
  const double s_bar = 1.0/(p->sigma_c_bar*p->sigma_c_bar);

  return 
    (p->r_bar - r)*s_bar /
    (r*(s - s_bar) - (at*s - 1 - p->r_bar*s_bar));
}

double
pomf_delta_equil_stdev(const pomf_params_t *p)
{
  const double sigma_a = p->sigma_a;
  const double sigma_a_tilde = pomf_sigma_a_tilde_star(p);
  const double kappa = p->kappa;
  const double sigma_c = p->sigma_c;
  const double sat_over_sc = sigma_a_tilde / sigma_c;
  const double gamma = p->gamma;

  const double gamma_delta = gamma + sat_over_sc;
  const double sigma_delta_sqr = sigma_a*sigma_a
    - 2*kappa*sigma_a*sigma_a_tilde + sigma_a_tilde*sigma_a_tilde;

  return sqrt(0.5*sigma_delta_sqr/gamma_delta);
}


double
pomf_r_star_diff(const double at, const double eta,
		 double *dr_dat,
		 double *dr_deta,
		 const pomf_params_t *p)
{
  const double eta_prime = 1.0 - eta;
  const double s = 1.0/(p->sigma_c*p->sigma_c);
  const double s_bar = 1.0/(p->sigma_c_bar*p->sigma_c_bar);

  const double P = (at*s - 1)*eta + (p->r_bar*s_bar)*eta_prime;
  const double Q = s*eta + s_bar*eta_prime;

  *dr_dat = s*eta / Q;
  *dr_deta = ((at*s - 1 - p->r_bar*s_bar)*Q - P*(s - s_bar)) / (Q*Q);

  return P / Q;
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
		 const pomf_params_t *p)
{
  int status = OK;

  const double N = eta*Ntot;
  const double Nb = (1 - eta)*Ntot;
  const double a_tilde = a - delta;
  
  const double rho = p->rho;
  const double r = pomf_r_star(a_tilde, eta, p);
  const double r_bar = p->r_bar;
  const double eta_prime = 1 - eta;
  const double sh_a_c        =       (a - r) / p->sigma_c;
  const double sh_a_tilde_c  = (a_tilde - r) / p->sigma_c;
  const double sh_r_c_bar    =   (r - r_bar) / p->sigma_c_bar;

  const double mu_N = r - rho + sh_a_c*sh_a_tilde_c;
  const double sigma_N = sh_a_tilde_c;
  const double mu_Nb = r_bar - r + sh_r_c_bar*sh_r_c_bar;
  const double sigma_Nb = sh_r_c_bar;

  *mu_Ntot = mu_N*eta + mu_Nb*eta_prime;
  *sigma_Ntot = sigma_N*eta + sigma_Nb*eta_prime;

  const double ma = a * sh_a_tilde_c / p->sigma_c;
  const double mb = r * sh_r_c_bar / p->sigma_c_bar;
  const double sa = sh_a_tilde_c;
  const double sb = sh_r_c_bar;

  const double Ns = N * N0_over_y0;
  const double Ntots = Ntot * Ntot0_over_ytot0;
  
  *mu_y = ma * Ns;
  *sigma_y = sa * Ns;

  *mu_ytot = (ma*eta + mb*eta_prime) * Ntots;
  *sigma_ytot = (sa*eta + sb*eta_prime) * Ntots;

  return status;  
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

int
pomf_limitfunc2(double bs[2*POMF_DIM], void *params)
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
  bs[2*POMF_I_ETA + 0] = 0.0;
  bs[2*POMF_I_ETA + 1] = 1.0;

  return OK;
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
pomf_limitfunc_alt2(double bs[POMF_DIM], void *params)
{
  pomf_limitfunc2(bs, params);
  bs[2*POMF_I_XI + 0] = eta2xi(bs[2*POMF_I_ETA + 0]);
  bs[2*POMF_I_XI + 1] = eta2xi(bs[2*POMF_I_ETA + 1]);

  return OK;
}

int
pomf_sdefunc(const double x[POMF_DIM], 
	     double mu_out[POMF_DIM],
	     double sigma_out[POMF_DIM*POMF_WDIM],
	     void *params)
{
  return pomf_sdefunc2(x, mu_out, sigma_out, NULL, params);
}

int
pomf_sdefunc2(const double x[POMF_DIM], 
	      double mu_out[POMF_DIM],
	      double sigma_out[POMF_DIM*POMF_WDIM],
	      double sigma_diff[POMF_DIM*POMF_DIM*POMF_WDIM],
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

  double dr_dat;
  double dr_deta;
  
  const double r = pomf_r_star_diff(a_tilde, eta, &dr_dat, &dr_deta, p);
  const double h_h_prime = eta*(1.0 - eta);

  const double dr_da = dr_dat;
  const double dr_ddelta = -dr_dat;

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

  const double dsh_a_tilde_c_da = (1.0 - dr_da) / p->sigma_c;
  const double dsh_a_tilde_c_ddelta = (-1.0 - dr_ddelta) / p->sigma_c;
  const double dsh_a_tilde_c_deta = -dr_deta / p->sigma_c;
  
  const double dsh_r_c_bar_da = dr_da / p->sigma_c_bar;
  const double dsh_r_c_bar_ddelta = dr_ddelta / p->sigma_c_bar;
  const double dsh_r_c_bar_deta = dr_deta / p->sigma_c_bar;

  const double dsigma_h_da = dsh_a_tilde_c_da - dsh_r_c_bar_da;
  const double dsigma_h_ddelta = dsh_a_tilde_c_ddelta - dsh_r_c_bar_ddelta;
  const double dsigma_h_deta = dsh_a_tilde_c_deta - dsh_r_c_bar_deta;

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

  if (sigma_diff != NULL)
    {
      /* a derivative */
      double *dsigma_da = &sigma_diff[POMF_I_A*POMF_WDIM*POMF_DIM];
      
      dsigma_da[POMF_I_A*POMF_WDIM + 0] = 0.0;
      dsigma_da[POMF_I_A*POMF_WDIM + 1] = 0.0;
      dsigma_da[POMF_I_A*POMF_WDIM + 2] = 0.0;
      
      dsigma_da[POMF_I_DELTA*POMF_WDIM + 0] = 0.0;
      dsigma_da[POMF_I_DELTA*POMF_WDIM + 1] = 0.0;
      dsigma_da[POMF_I_DELTA*POMF_WDIM + 2] = 0.0;  
      
      dsigma_da[POMF_I_ETA*POMF_WDIM + 0] = h_h_prime*dsigma_h_da;
      dsigma_da[POMF_I_ETA*POMF_WDIM + 1] = 0.0;
      dsigma_da[POMF_I_ETA*POMF_WDIM + 2] = 0.0;
  
      /* delta derivative */
      double *dsigma_ddelta = &sigma_diff[POMF_I_DELTA*POMF_WDIM*POMF_DIM];
      
      dsigma_ddelta[POMF_I_A*POMF_WDIM + 0] = 0.0;
      dsigma_ddelta[POMF_I_A*POMF_WDIM + 1] = 0.0;
      dsigma_ddelta[POMF_I_A*POMF_WDIM + 2] = 0.0;
      
      dsigma_ddelta[POMF_I_DELTA*POMF_WDIM + 0] = 0.0;
      dsigma_ddelta[POMF_I_DELTA*POMF_WDIM + 1] = 0.0;
      dsigma_ddelta[POMF_I_DELTA*POMF_WDIM + 2] = 0.0;  
      
      dsigma_ddelta[POMF_I_ETA*POMF_WDIM + 0] = h_h_prime*dsigma_h_ddelta;
      dsigma_ddelta[POMF_I_ETA*POMF_WDIM + 1] = 0.0;
      dsigma_ddelta[POMF_I_ETA*POMF_WDIM + 2] = 0.0;
      
      /* eta derivative */
      double *dsigma_deta = &sigma_diff[POMF_I_ETA*POMF_WDIM*POMF_DIM];
      
      dsigma_deta[POMF_I_A*POMF_WDIM + 0] = 0.0;
      dsigma_deta[POMF_I_A*POMF_WDIM + 1] = 0.0;
      dsigma_deta[POMF_I_A*POMF_WDIM + 2] = 0.0;
      
      dsigma_deta[POMF_I_DELTA*POMF_WDIM + 0] = 0.0;
      dsigma_deta[POMF_I_DELTA*POMF_WDIM + 1] = 0.0;
      dsigma_deta[POMF_I_DELTA*POMF_WDIM + 2] = 0.0;  
      
      dsigma_deta[POMF_I_ETA*POMF_WDIM + 0] = h_h_prime*dsigma_h_deta + (1 - 2*eta)*sigma_h;
      dsigma_deta[POMF_I_ETA*POMF_WDIM + 1] = 0.0;
      dsigma_deta[POMF_I_ETA*POMF_WDIM + 2] = 0.0;
    }
  
  return 0;
}

int
pomf_sdefunc_alt(const double x[POMF_DIM], 
		 double mu_out[POMF_DIM],
		 double sigma_out[POMF_DIM*POMF_WDIM],
		 void *params)
{
  return pomf_sdefunc_alt2(x, mu_out, sigma_out, NULL, params);
}

int
pomf_sdefunc_alt2(const double x[POMF_DIM], 
		  double mu_out[POMF_DIM],
		  double sigma_out[POMF_DIM*POMF_WDIM],
		  double sigma_diff[POMF_DIM*POMF_DIM*POMF_WDIM],
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

  double dr_dat;
  double dr_deta;
  const double r = pomf_r_star_diff(a_tilde, eta, &dr_dat, &dr_deta, p);
  const double dr_da = dr_dat;
  const double dr_ddelta = -dr_dat;

  const double sh_a_c        =        (a - r) / p->sigma_c;
  const double sh_a_tilde_c  =  (a_tilde - r) / p->sigma_c;
  const double sh_r_c_bar    = (r - p->r_bar) / p->sigma_c_bar;
  const double r_diff        =  r - p->r_bar;
  const double rho_diff      = p->rho - p->rho_bar;
  const double sigma_a_tilde = pomf_sigma_a_tilde_star(p);
  const double sat_over_sc   = sigma_a_tilde / p->sigma_c;

  const double dsh_a_tilde_c_da = (1.0 - dr_da) / p->sigma_c;
  const double dsh_a_tilde_c_ddelta = (-1.0 - dr_ddelta) / p->sigma_c;
  const double dsh_a_tilde_c_deta = -dr_deta / p->sigma_c;
  
  const double dsh_r_c_bar_da = dr_da / p->sigma_c_bar;
  const double dsh_r_c_bar_ddelta = dr_ddelta / p->sigma_c_bar;
  const double dsh_r_c_bar_deta = dr_deta / p->sigma_c_bar;

  const double dsigma_h_da = dsh_a_tilde_c_da - dsh_r_c_bar_da;
  const double dsigma_h_ddelta = dsh_a_tilde_c_ddelta - dsh_r_c_bar_ddelta;
  const double dsigma_h_deta = dsh_a_tilde_c_deta - dsh_r_c_bar_deta;

  const double deta_dxi = eta*eta_prime;  

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
  
  sigma_out[POMF_I_XI*POMF_WDIM + 0] = sigma_h;
  sigma_out[POMF_I_XI*POMF_WDIM + 1] = 0.0;
  sigma_out[POMF_I_XI*POMF_WDIM + 2] = 0.0;

  if (sigma_diff != NULL)
    {
      /* a derivative */
      double *dsigma_da = &sigma_diff[POMF_I_A*POMF_WDIM*POMF_DIM];
      
      dsigma_da[POMF_I_A*POMF_WDIM + 0] = 0.0;
      dsigma_da[POMF_I_A*POMF_WDIM + 1] = 0.0;
      dsigma_da[POMF_I_A*POMF_WDIM + 2] = 0.0;
      
      dsigma_da[POMF_I_DELTA*POMF_WDIM + 0] = 0.0;
      dsigma_da[POMF_I_DELTA*POMF_WDIM + 1] = 0.0;
      dsigma_da[POMF_I_DELTA*POMF_WDIM + 2] = 0.0;  
      
      dsigma_da[POMF_I_XI*POMF_WDIM + 0] = dsigma_h_da;
      dsigma_da[POMF_I_XI*POMF_WDIM + 1] = 0.0;
      dsigma_da[POMF_I_XI*POMF_WDIM + 2] = 0.0;
      
      /* delta derivative */
      double *dsigma_ddelta = &sigma_diff[POMF_I_DELTA*POMF_WDIM*POMF_DIM];
      
      dsigma_ddelta[POMF_I_A*POMF_WDIM + 0] = 0.0;
      dsigma_ddelta[POMF_I_A*POMF_WDIM + 1] = 0.0;
      dsigma_ddelta[POMF_I_A*POMF_WDIM + 2] = 0.0;
      
      dsigma_ddelta[POMF_I_DELTA*POMF_WDIM + 0] = 0.0;
      dsigma_ddelta[POMF_I_DELTA*POMF_WDIM + 1] = 0.0;
      dsigma_ddelta[POMF_I_DELTA*POMF_WDIM + 2] = 0.0;  
      
      dsigma_ddelta[POMF_I_XI*POMF_WDIM + 0] = dsigma_h_ddelta;
      dsigma_ddelta[POMF_I_XI*POMF_WDIM + 1] = 0.0;
      dsigma_ddelta[POMF_I_XI*POMF_WDIM + 2] = 0.0;
      
      /* eta derivative */
      double *dsigma_deta = &sigma_diff[POMF_I_XI*POMF_WDIM*POMF_DIM];
      
      dsigma_deta[POMF_I_A*POMF_WDIM + 0] = 0.0;
      dsigma_deta[POMF_I_A*POMF_WDIM + 1] = 0.0;
      dsigma_deta[POMF_I_A*POMF_WDIM + 2] = 0.0;
      
      dsigma_deta[POMF_I_DELTA*POMF_WDIM + 0] = 0.0;
      dsigma_deta[POMF_I_DELTA*POMF_WDIM + 1] = 0.0;
      dsigma_deta[POMF_I_DELTA*POMF_WDIM + 2] = 0.0;  
      
      dsigma_deta[POMF_I_XI*POMF_WDIM + 0] = dsigma_h_deta * deta_dxi;
      dsigma_deta[POMF_I_XI*POMF_WDIM + 1] = 0.0;
      dsigma_deta[POMF_I_XI*POMF_WDIM + 2] = 0.0;
    }
  
  return 0;
}

void
mmadd(const double *m, int n, int transp, int howmany,
      const double *x, int xstride, int xnext,
      const double *a, int astride, int anext,
      double *y, int ystride, int ynext)
{
  int i, j, k;

  for (k = 0; k < howmany; ++k)
    {
      const double *xthis = x + k*xnext;
      const double *athis = a + k*anext;
      double *ythis = y + k*ynext;

      for (i = 0; i < n; ++i)
	{
	  double s = 0.0;
	  for (j = 0; j < n; ++j)
	    {
	      if (transp)
		s += m[j*n + i]*xthis[j*xstride];
	      else
		s += m[i*n + j]*xthis[j*xstride];
	    }

	  if (a)
	    ythis[i*ystride] = s + athis[i*astride];
	  else
	    ythis[i*ystride] = s;
	}
    }
}

void
dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);

int
pomf_setup_alt3(const double x_star[POMF_DIM],
		pomf_params_t *p)
{
  int status = OK;

  double mu_x[POMF_DIM];
  double sigma_x[POMF_DIM*POMF_WDIM];

  CHECK_( pomf_sdefunc_alt2(x_star, mu_x, sigma_x, NULL, p) );

  memcpy(p->x_star, x_star, POMF_DIM*sizeof(*x_star));

  /* Compute the diffusion matrix. */
  double D[POMF_DIM2];
  int i, j, k;
  for (i = 0; i < POMF_DIM; ++i)
    {
      for (j = i; j < POMF_DIM; ++j)
	{
	  double s = 0.0;
	  for (k = 0; k < POMF_WDIM; ++k)
	    {
	      s += sigma_x[i*POMF_WDIM + k]*sigma_x[j*POMF_WDIM + k];
	    }

	  D[i*POMF_DIM + j] = D[j*POMF_DIM + i] = 0.5 * s;
	}
    }

  memcpy(p->U, D, POMF_DIM2*sizeof(double));

  /* Diagonalize the diffusion matrix */
  char jobz = 'V';
  char uplo = 'U';
  int n = POMF_DIM;
  int lda = POMF_DIM;
  double w[POMF_DIM];
  double work[3*POMF_DIM-1];
  int lwork = 3*POMF_DIM-1;
  int info;
  
  dsyev_(&jobz, &uplo, &n, p->U, &lda, w, work, &lwork, &info);
  CHECK( info, FAILED );

  fprintf(stderr, "U = \n");
  for (i = 0; i < POMF_DIM; ++i)
    {
      for (j = 0; j < POMF_DIM; ++j)
  	{
	  fprintf(stderr, "%25.15lf ", p->U[i*POMF_DIM + j]);
  	}
      fprintf(stderr, "\n");
    }

  double DUT[POMF_DIM2];
  double UDUT[POMF_DIM2];

  mmadd(D, POMF_DIM, 0, POMF_DIM,
	p->U, 1, POMF_DIM, 
	NULL, 0, 0,
	DUT, POMF_DIM, 1);

  mmadd(p->U, POMF_DIM, 0, POMF_DIM,
	DUT, POMF_DIM, 1,
	NULL, 0, 0,
	UDUT, POMF_DIM, 1);

  fprintf(stderr, "U D UT = \n");
  for (i = 0; i < POMF_DIM; ++i)
    {
      for (j = 0; j < POMF_DIM; ++j)
  	{
	  fprintf(stderr, "%25.15lf ", UDUT[i*POMF_DIM + j]);
  	}
      fprintf(stderr, "\n");
    }

  fprintf(stderr, "\nW = ");
  for (j = 0; j < POMF_DIM; ++j)
    {
      fprintf(stderr, "%25.15lf ", w[j]);
    }
  fprintf(stderr, "\n");

  
  for (i = 0; i < POMF_DIM; ++i)
    {
      for (j = i + 1; j < POMF_DIM; ++j)
  	{
  	  double t = p->U[i*POMF_DIM + j];
  	  p->U[i*POMF_DIM + j] = p->U[j*POMF_DIM + i];
  	  p->U[j*POMF_DIM + i] = t;
  	}
    }

  double y[POMF_DIM] = {0.1, 0.0, 0.0};
  double mu_y[POMF_DIM];
  double sigma_y[POMF_DIM*POMF_WDIM];
  double sigma_y_diff[POMF_DIM*POMF_DIM*POMF_WDIM];
  pomf_sdefunc_alt3(y, mu_y, sigma_y, sigma_y_diff, p);

  fprintf(stderr, "\nD(y) = \n");
  for (i = 0; i < POMF_DIM; ++i)
    {
      for (j = 0; j < POMF_DIM; ++j)
	{
	  double s = 0.0;
	  for (k = 0; k < POMF_WDIM; ++k)
	    {
	      s += sigma_y[i*POMF_WDIM + k]*sigma_y[j*POMF_WDIM + k];
	    }
	  fprintf(stderr, "%25.15lf ", 0.5*s);
	}
      fprintf(stderr, "\n");
    }
  
   status = 1;
    
 exit:

  return status;  
}

int
pomf_sdefunc_alt3(const double y[POMF_DIM], 
		  double mu_out[POMF_DIM],
		  double sigma_out[POMF_DIM*POMF_WDIM],
		  double sigma_diff[POMF_DIM*POMF_DIM*POMF_WDIM],
		  void *params)
{
  int status = OK;
  pomf_params_t *p = params;

  double x[POMF_DIM];
  double mu_x[POMF_DIM];
  double sigma_x[POMF_DIM*POMF_WDIM];
  double sigma_x_diff[POMF_DIM*POMF_DIM*POMF_WDIM];

  mmadd(p->U, POMF_DIM, 0, 1,
	y, 1, 0,
	p->x_star, 1, 0,
	x, 1, 0);

  CHECK_( pomf_sdefunc_alt2(x, mu_x, sigma_x, sigma_x_diff, params) );

  mmadd(p->U, POMF_DIM, 1, 1,
	mu_x, 1, 0,
	NULL, 0, 0,
	mu_out, 1, 0);

  mmadd(p->U, POMF_DIM, 1, POMF_WDIM,
	sigma_x, POMF_WDIM, 1,
	NULL, 0, 0,
	sigma_out, POMF_WDIM, 1);

  int i, j, k, l, m;
  for (i = 0; i < POMF_DIM2*POMF_WDIM; ++i) sigma_diff[i] = 0.0;
  for (l = 0; l < POMF_DIM; ++l) /* derivative along y-space axis l */
    {
      for (i = 0; i < POMF_DIM; ++i) /* row of derivative mtx */
	{
	  for (j = 0; j < POMF_WDIM; ++j) /* col of derivative mtx */
	    {
	      for (k = 0; k < POMF_DIM; ++k) /* UT left-multipl. */
		{
		  for (m = 0; m < POMF_DIM; ++m) /* x-space derivative m */
		    {
		      sigma_diff[l*POMF_DIM*POMF_WDIM + i*POMF_WDIM + j]
			+=
			p->U[k*POMF_DIM + i]
			* sigma_x_diff[m*POMF_DIM*POMF_WDIM + k*POMF_WDIM + j]
			* p->U[m*POMF_DIM + l];
		    }
		}
	    }
	}
    }
  
 exit:

  return status;
}

int
pomf_limitfunc_alt3(double bs[2*POMF_DIM], void *params)
{
  pomf_params_t *p = params;
  
  pomf_limitfunc2(bs, params);
  bs[2*POMF_I_XI + 0] = eta2xi(bs[2*POMF_I_ETA + 0]);
  bs[2*POMF_I_XI + 1] = eta2xi(bs[2*POMF_I_ETA + 1]);

  int i;
  for (i = 0; i < POMF_DIM; ++i)
    {
      bs[2*i + 0] = -3.0;
      bs[2*i + 1] = +3.0;

      /* bs[2*i + 0] -= p->x_star[i]; */
      /* bs[2*i + 1] -= p->x_star[i]; */
    }
  /* x = Uy + x* */

  bs[2*0 + 0] = -0.5;
  bs[2*0 + 1] = +0.5;

  bs[2*1 + 0] = -0.5;
  bs[2*1 + 1] = +0.5;

  bs[2*2 + 0] = -2.0;
  bs[2*2 + 1] = +2.0;

  return OK;  
}

