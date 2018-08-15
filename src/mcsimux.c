
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pomf_mcsimux.h"
#include "pomf_model.h"

const int pomf_xdim = POMF_XDIM;
const int pomf_xwdim = POMF_XWDIM;
const int pomf_xsim_ia = POMF_XSIM_IA;
const int pomf_xsim_idelta = POMF_XSIM_IDELTA;
const int pomf_xsim_ieta = POMF_XSIM_IETA;
const int pomf_xsim_intot = POMF_XSIM_INTOT;
const int pomf_xsim_iy = POMF_XSIM_IY;
const int pomf_xsim_iytot = POMF_XSIM_IYTOT;

const int pomf_xelems = POMF_XELEMS;
const int pomf_xiatilde = POMF_XIATILDE;
const int pomf_xir = POMF_XIR;
const int pomf_xidelta = POMF_XIDELTA;
const int pomf_xintot = POMF_XINTOT;
const int pomf_xiy = POMF_XIY;
const int pomf_xiytot = POMF_XIYTOT;


void
pomf_mcsimux_cache(pomf_pcachex_t *pc,
		   const pomf_paramx_t *px)
{
  const double kappa = px->kappa;
  const double kappa_prime = sqrt(1 - kappa*kappa);

  /* sigma_a_tilde code */
  const double t0 = px->gamma*px->sin_theta*px->sigma_c;
  const double t1 = kappa*px->sigma_a + px->gamma*px->sin_theta*px->sigma_c;
  const double t2 = px->sigma_a*kappa_prime;
  
  const double sigma_a_tilde = -t0 + sqrt(t1*t1 + t2*t2);

  /* delta_equil_stdev code */
  const double sat_over_sc = sigma_a_tilde / (px->sin_theta*px->sigma_c);
  const double gamma_delta = px->gamma + sat_over_sc;
  const double sigma_delta_sqr = px->sigma_a*px->sigma_a
    - 2*kappa*px->sigma_a*sigma_a_tilde + sigma_a_tilde*sigma_a_tilde;
  const double delta_equil_stdev = sqrt(0.5*sigma_delta_sqr/gamma_delta);

  /* Fill the cache structure */
  pc->delta_equil_stdev = delta_equil_stdev;
  pc->gamma_delta = gamma_delta;
  pc->s = 1.0/(px->sigma_c*px->sigma_c);
  pc->s_bar = 1.0/(px->sigma_c_bar*px->sigma_c_bar);
  pc->sigma_a_tilde = sigma_a_tilde;
  pc->kappa_prime = sqrt(1.0 - px->kappa*px->kappa);
  pc->sigma_c_inv = 1.0/px->sigma_c;
  pc->sigma_c_bar_inv = 1.0/px->sigma_c_bar;
}

pomf_mcsimux_t *
pomf_mcsimux_make(size_t n_samples)
{
  int status = OK;
  pomf_mcsimux_t *self = NULL;

  CHECK_NULL( self = malloc(sizeof(*self)), OUT_OF_MEM );

  self->rng = NULL;
  self->n_samples = n_samples;
  self->samples = NULL;

  const size_t ssize = self->n_samples*POMF_XELEMS;

  CHECK_NULL( self->rng = gsl_rng_alloc(gsl_rng_default), OUT_OF_MEM );
  CHECK_NULL( self->samples = malloc(ssize*sizeof(*self->samples)), OUT_OF_MEM );
  
 exit:
  if (status)
    {
      Xfree_(self, pomf_mcsimux_free);
      self = NULL;
    }

  return self;
}

void
pomf_mcsimux_free(pomf_mcsimux_t *self)
{
  if (self)
    {
      Xfree_(self->rng, gsl_rng_free);
      Xfree(self->samples);
      Xfree(self);
    }
}

int
pomf_mcsimux_init(pomf_mcsimux_t *self,
		  pomf_paramx_t *px)
{
  int status = OK;

  memcpy(&self->px, px, sizeof(*px));
  pomf_mcsimux_cache(&self->pc, px);

  const double delta_init_mean = 0.0;
  const double delta_init_stdev = self->pc.delta_equil_stdev;

  int i_sample;

  for (i_sample = 0; i_sample < self->n_samples; ++i_sample)
    {
      double delta_init = delta_init_mean
	+ gsl_ran_gaussian_ziggurat(self->rng, delta_init_stdev);

      self->samples[POMF_XELEMS*i_sample + POMF_XIATILDE] = px->a_tilde_init;
      self->samples[POMF_XELEMS*i_sample + POMF_XIR] = px->r_init;
      self->samples[POMF_XELEMS*i_sample + POMF_XIDELTA] = delta_init;
      self->samples[POMF_XELEMS*i_sample + POMF_XINTOT] = 1.0;
      self->samples[POMF_XELEMS*i_sample + POMF_XIY] = 1.0;
      self->samples[POMF_XELEMS*i_sample + POMF_XIYTOT] = 1.0;
    }

  return status;
}

int
pomf_mcsimux_init_wobs(pomf_mcsimux_t *self, double r, double y, double ytot)
{
  int status = OK;
  int i_sample;

  for (i_sample = 0; i_sample < self->n_samples; ++i_sample)
    {
      self->samples[POMF_XELEMS*i_sample + POMF_XIR] = r;
      self->samples[POMF_XELEMS*i_sample + POMF_XIY] = y;
      self->samples[POMF_XELEMS*i_sample + POMF_XIYTOT] = ytot;
    }

  return status;
}

double
pomf_r_starx_(const double at, const double eta,
	      const pomf_paramx_t *px,
	      const pomf_pcachex_t *pc)
{
  const double eta_prime = 1.0 - eta;
  return
    ((at*pc->s - 1)*eta + px->r_bar*pc->s_bar*eta_prime) /
    (pc->s*eta + pc->s_bar*eta_prime);
}

double
pomf_r_starx(const double at, const double eta,
	     const pomf_paramx_t *px)
{
  pomf_pcachex_t pc;
  pomf_mcsimux_cache(&pc, px);
  return pomf_r_starx_(at, eta, px, &pc);
}

double
pomf_eta_starx_(const double at, const double r,
		const pomf_paramx_t *px,
		const pomf_pcachex_t *pc)
{
  return 
    (px->r_bar - r)*pc->s_bar /
    (r*(pc->s - pc->s_bar) - (at*pc->s - 1 - px->r_bar*pc->s_bar));
}

double
pomf_eta_starx(const double at, const double r,
	       const pomf_paramx_t *px)
{
  pomf_pcachex_t pc;
  pomf_mcsimux_cache(&pc, px);
  return pomf_eta_starx_(at, r, px, &pc);
}

double
pomf_delta_equil_stdevx(const pomf_paramx_t *px)
{
  pomf_pcachex_t pc;
  pomf_mcsimux_cache(&pc, px);
  return pc.delta_equil_stdev;
}

double
pomf_sigma_a_tilde_starx(const pomf_paramx_t *px)
{
  pomf_pcachex_t pc;
  pomf_mcsimux_cache(&pc, px);
  return pc.sigma_a_tilde;
}

/**
 * Compute the drift and diffusion for the extended simulation system.
 */
int
pomf_musigmax(const double x[POMF_XDIM],
	      double mu[POMF_XDIM],
	      double sigma[POMF_XDIM][POMF_XWDIM],
	      const pomf_paramx_t *px,
	      const pomf_pcachex_t *pc)
{
  int status = OK;

  const double a = x[POMF_XSIM_IA];
  const double delta = x[POMF_XSIM_IDELTA];
  const double eta = x[POMF_XSIM_IETA];
  const double Ntot = x[POMF_XSIM_INTOT];
  const double y = x[POMF_XSIM_IY];
  const double ytot = x[POMF_XSIM_IYTOT];

  const double at = a - delta;
  const double eta_prime = 1.0 - eta;
  const double N = eta * Ntot;
  const double N_bar = eta_prime * Ntot;
  const double eta_eta_prime = eta * eta_prime;

  const double r = pomf_r_starx_(at, eta, px, pc);
  const double r_diff = r - px->r_bar;
  const double rho_diff = px->rho - px->rho_bar;

  const double sh_a = (a - r) * pc->sigma_c_inv;
  const double sh_at = (at - r) * pc->sigma_c_inv;
  const double sh_r = r_diff * pc->sigma_c_bar_inv;

  const double sh_at_st = sh_at * px->sin_theta;
  const double sh_r_st = sh_r * px->sin_theta_bar;
  const double sigma_eta = sh_at_st - sh_r_st;

  const double mu_N = r - px->rho + sh_at*sh_a;
  const double sigma_N = sh_at_st;
  const double mu_N_bar = r - px->rho_bar + sh_r*sh_r;
  const double sigma_N_bar = sh_r_st;

  const double mu_Y = a * pc->sigma_c_inv * sh_at * N;
  const double mu_Y_bar = r * pc->sigma_c_bar_inv * sh_r * N_bar;
  const double sigma_Y = sh_at_st * N;
  const double sigma_Y_bar = sh_r_st * N_bar;

  mu[POMF_XSIM_IA] = px->gamma*(px->a_bar - a);
  mu[POMF_XSIM_IDELTA] = -pc->gamma_delta*delta;
  mu[POMF_XSIM_IETA] = eta_eta_prime *
    (r_diff - rho_diff + sh_a * sh_at - sh_r*sh_r
     - sigma_eta * (sh_at_st*eta + sh_r_st*eta_prime));

  mu[POMF_XSIM_INTOT] = (mu_N*eta + mu_N_bar*eta_prime) * Ntot;
  mu[POMF_XSIM_IY] = mu_Y * px->N0_over_y0;
  mu[POMF_XSIM_IYTOT] = (mu_Y + mu_Y_bar) * px->Ntot0_over_ytot0;

  sigma[POMF_XSIM_IA][0] = px->kappa*px->sigma_a;
  sigma[POMF_XSIM_IA][1] = pc->kappa_prime*px->sigma_a;
  sigma[POMF_XSIM_IDELTA][0] = sigma[POMF_XSIM_IA][0] - pc->sigma_a_tilde;
  sigma[POMF_XSIM_IDELTA][1] = sigma[POMF_XSIM_IA][1];
  sigma[POMF_XSIM_IETA][0] = eta_eta_prime * sigma_eta;
  sigma[POMF_XSIM_IETA][1] = 0.0;

  sigma[POMF_XSIM_INTOT][0] = (sigma_N*eta + sigma_N_bar*eta_prime) * Ntot;
  sigma[POMF_XSIM_INTOT][1] = 0.0;
  sigma[POMF_XSIM_IY][0] = sigma_Y * px->N0_over_y0;
  sigma[POMF_XSIM_IY][1] = 0.0;
  sigma[POMF_XSIM_IYTOT][0] = (sigma_Y + sigma_Y_bar) * px->Ntot0_over_ytot0;
  sigma[POMF_XSIM_IYTOT][1] = 0.0;

 exit:
  return status;
}

static int
extstep(pomf_mcsimux_t *self,
	const double *x,
	double *dx,
	double dt,
	double sqrt_dt)
{
  int status = OK;

  pomf_paramx_t *px = &self->px;
  pomf_pcachex_t *pc = &self->pc;
  gsl_rng *rng = self->rng;

  double dw[POMF_XWDIM];
  int i, j;

  if (x[POMF_XSIM_IETA] < 0.0 || x[POMF_XSIM_IETA] > 1.0)
    {
      ERRMSG("eta variable out of bounds at start of euler step\n");
      ERRMSG("eta = %lf, a = %lf, delta = %lf\n",
	     x[POMF_XSIM_IETA], x[POMF_XSIM_IA], x[POMF_XSIM_IDELTA]);
	     
      FAILWITH(FAILED);
    }
  
  for (i = 0; i < POMF_XWDIM; ++i)
    dw[i] = gsl_ran_gaussian_ziggurat(rng, sqrt_dt);

  double mu[POMF_XDIM];
  double sigma[POMF_XDIM][POMF_XWDIM];

  CHECK_( pomf_musigmax(x, mu, sigma, px, pc) );

  for (i = 0; i < POMF_XDIM; ++i)
    {
      dx[i] = mu[i]*dt;
      for (j = 0; j < POMF_XWDIM; ++j) dx[i] += sigma[i][j]*dw[j];
    }

  /* Enforce [0,1] domain for eta by reflection, if finite time step
   * attempts to jump out of the allowed range. */
  if (x[POMF_XSIM_IETA] + dx[POMF_XSIM_IETA] < 0.0)
    dx[POMF_XSIM_IETA] = -2*x[POMF_XSIM_IETA] - dx[POMF_XSIM_IETA];

  if (x[POMF_XSIM_IETA] + dx[POMF_XSIM_IETA] > 1.0)
    dx[POMF_XSIM_IETA] = 2*(1 - x[POMF_XSIM_IETA]) - dx[POMF_XSIM_IETA];

  if (x[POMF_XSIM_IETA] + dx[POMF_XSIM_IETA] < 0.0 ||
      x[POMF_XSIM_IETA] + dx[POMF_XSIM_IETA] > 1.0)
    {
      ERRMSG("Euler step out of bounds even after reflection, too large dt?\n");
      FAILWITH(FAILED);
    }

 exit:

  return status;
}

int
pomf_mcsimux_run(pomf_mcsimux_t *self,
		 double t_end,
		 double dt)
{
  int status = OK;
  int i_sample;
  
  const double sqrt_dt = sqrt(dt);

  for (i_sample = 0; i_sample < self->n_samples; ++i_sample)
    {
      double x[POMF_XDIM], dx[POMF_XDIM];

      double delta_init = self->samples[POMF_XELEMS*i_sample + POMF_XIDELTA];
      double a_tilde_init = self->samples[POMF_XELEMS*i_sample + POMF_XIATILDE];
      double r_init = self->samples[POMF_XELEMS*i_sample + POMF_XIR];
      double a_init = a_tilde_init + delta_init;
      double eta_init = pomf_eta_starx_(a_tilde_init, r_init, &self->px, &self->pc);

      double Ntot_init = self->samples[POMF_XELEMS*i_sample + POMF_XINTOT];
      double y_init = self->samples[POMF_XELEMS*i_sample + POMF_XIY];
      double ytot_init = self->samples[POMF_XELEMS*i_sample + POMF_XIYTOT];

      if (eta_init < 0.0) eta_init = 0.0;
      else if (eta_init > 1.0) eta_init = 1.0;

      x[POMF_XSIM_IETA] = eta_init;
      x[POMF_XSIM_IA] = a_init;
      x[POMF_XSIM_IDELTA] = delta_init;
      x[POMF_XSIM_INTOT] = Ntot_init;
      x[POMF_XSIM_IY] = y_init;
      x[POMF_XSIM_IYTOT] = ytot_init;

      double t = 0.0;

      while (t < t_end)
	{
	  if (t + dt <= t_end)
	    {
	      CHECK_( extstep(self, x, dx, dt, sqrt_dt) );
	      t += dt;
	    }
	  else
	    {
	      CHECK_( extstep(self, x, dx, t_end - t, sqrt(t_end - t)) );
	      t = t_end;
	    }

	  int i;
	  for (i = 0; i < POMF_XDIM; ++i) x[i] += dx[i];
	}
      
      /* ****** */

      double a_tilde = x[POMF_XSIM_IA] - x[POMF_XSIM_IDELTA];
      double r = pomf_r_starx_(a_tilde, x[POMF_XSIM_IETA], &self->px, &self->pc);
      double delta = x[POMF_XSIM_IDELTA];

      self->samples[POMF_XELEMS*i_sample + POMF_XIATILDE] = a_tilde;
      self->samples[POMF_XELEMS*i_sample + POMF_XIR] = r;
      self->samples[POMF_XELEMS*i_sample + POMF_XIDELTA] = delta;
      self->samples[POMF_XELEMS*i_sample + POMF_XINTOT] = x[POMF_XSIM_INTOT];
      self->samples[POMF_XELEMS*i_sample + POMF_XIY] = x[POMF_XSIM_IY];
      self->samples[POMF_XELEMS*i_sample + POMF_XIYTOT] = x[POMF_XSIM_IYTOT];
    }

 exit:

  return status;
}

double *
pomf_mcsimux_get_samples(pomf_mcsimux_t *self)
{
  return self->samples;
}
