
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pomf_mcsimux.h"

int
main()
{
  int status = OK;

  pomf_mcsimux_t *mcx = NULL;
  const int n_samples = 1;
  const double dt = 0.001;
  const double t_unit = 1.0;
  const int n_steps = 50;
  int i_step;

  pomf_paramx_t p;
    
  const double a_mean = 0.15;
  const double a_stddev = a_mean / 3.0;
  const double a_revrate = 0.1;
  p.a_bar         = a_mean;
  p.gamma         = a_revrate;
  p.sigma_a       = sqrt(2*a_revrate)*a_stddev;
  p.sigma_c       = a_mean * 2.0;
  p.sigma_c_bar   = a_mean * 2.0;
  p.kappa         = 0.0;
  p.r_bar         = 0.010;
  p.rho           = 0.040;
  p.rho_bar       = 0.030;
  p.a_tilde_init  = p.a_bar;
  p.sin_theta     = 0.1;
  p.sin_theta_bar = 0.1;

  const double r_init_0 = pomf_r_starx(p.a_tilde_init, 0.0, &p);
  const double r_init_1 = pomf_r_starx(p.a_tilde_init, 1.0, &p);

  p.r_init       = 0.5*(r_init_0 + r_init_1);
  p.N0_over_y0       = 1.0;
  p.Ntot0_over_ytot0 = 1.0;

  mcx = pomf_mcsimux_make(n_samples);

  CHECK_( pomf_mcsimux_init(mcx, &p) );

  for (i_step = 1; i_step <= n_steps; ++i_step)
    {
      CHECK_( pomf_mcsimux_run(mcx, t_unit, dt) );

      printf(" %25.15lf", i_step*t_unit);
      int i_sample;
      for (i_sample = 0; i_sample < n_samples; ++i_sample)
	{
	  printf(" %25.15lf %25.15lf %25.15lf %25.15le %25.15le %25.15le",
		 mcx->samples[i_sample*POMF_XELEMS + POMF_XIATILDE],   /* 2 */
		 mcx->samples[i_sample*POMF_XELEMS + POMF_XIR],        /* 3 */
		 mcx->samples[i_sample*POMF_XELEMS + POMF_XIDELTA],    /* 4 */
		 mcx->samples[i_sample*POMF_XELEMS + POMF_XINTOT],     /* 5 */
		 mcx->samples[i_sample*POMF_XELEMS + POMF_XIY],        /* 6 */
		 mcx->samples[i_sample*POMF_XELEMS + POMF_XIYTOT]);    /* 7 */
	}
      printf(" ?\n");
      /*
      for (i_sample = 0; i_sample < n_samples; ++i_sample)
	{
	  mcx->samples[i_sample*POMF_XELEMS + POMF_XIY] = 0.0;
	  mcx->samples[i_sample*POMF_XELEMS + POMF_XIYTOT] = 0.0;
	}
      */
    }

 exit:
  if (status) fprintf(stderr, "There were errors.\n");

  if (mcx) pomf_mcsimux_free(mcx);
}
