
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "pomf_model.h"

#include "fps.h"
#include "umfsolve.h"


double
rng_fun(double m, double s, void *state)
{
  gsl_rng *rng = state;
  return m + gsl_ran_gaussian_ziggurat(rng, s);
}

int
main()
{
  int status = OK;

  /* Model parameters --------------------------------------------- */
  pomf_params_t p;

  const double a_mean = 0.1;
  const double a_stddev = 0.05;
  const double a_revrate = 0.05;
  p.a_bar       = a_mean;
  p.gamma       = a_revrate;
  p.sigma_a     = sqrt(2*a_revrate)*a_stddev;
  p.sigma_c     = 0.15;
  p.sigma_c_bar = 0.15;
  p.kappa       = 0.00;
  p.r_bar       = 0.02;
  p.rho         = 0.07;
  p.rho_bar     = 0.06;

  /* const double a_mean = 0.070; */
  /* const double a_stddev = 0.035; */
  /* const double a_revrate = 0.5; */
  /* p.a_bar       = a_mean; */
  /* p.gamma       = a_revrate; */
  /* p.sigma_a     = sqrt(2*a_revrate)*a_stddev; */
  /* p.sigma_c     = 0.025; */
  /* p.sigma_c_bar = 0.025; */
  /* p.kappa       = 0.500; */
  /* p.r_bar       = 0.010; */
  /* p.rho         = 0.070;   /\* Note: Model depends only on difference rho - rho_bar *\/ */
  /* p.rho_bar     = 0.050; */
  
  int ns[POMF_DIM] = {41, 41, 41};
  // kfe_sde_t sde = { POMF_DIM, POMF_WDIM, pomf_sdefunc_alt3, pomf_limitfunc_alt3 };
  kfe_sde_t sde = { POMF_DIM, POMF_WDIM, pomf_sdefunc2, pomf_limitfunc2 };
  kfe_problem_t *prob = NULL;
  kfe_solver_t solver = kfe_solver_nil;
  cs *kfe = NULL;
  cs *ident = NULL;
  cs *bweuler = NULL;
  umfsolver_t umfsolver = UMFSOLVER_NIL;

  double x_init[POMF_DIM];
  /* x_init[0] = 0.0; */
  /* x_init[1] = 0.0; */
  /* x_init[2] = 0.0; // eta2xi(0.5); */
  x_init[POMF_I_A] = a_mean;
  x_init[POMF_I_DELTA] = 0.1;
  x_init[POMF_I_XI] = 0.5; // eta2xi(0.5);

  double *data_init = NULL, *data = NULL, *data_temp = NULL;
  const double dt = 0.1;
  const int niter = 100;
  const int n_samples = 100000;
  int iter;
  gsl_rng *rng = NULL;
  kfe_mcsolver_t mcsolver = kfe_mcsolver_nil;

  // CHECK_( pomf_setup_alt3(x_init, &p) );

  /* p.x_star[POMF_I_A] = 0.0; */
  /* p.x_star[POMF_I_DELTA] = 0.0; */
  /* p.x_star[POMF_I_XI] = 0.0; */

  /* p.U[0*POMF_DIM + 0] = 1.0; */
  /* p.U[0*POMF_DIM + 1] = 0.0; */
  /* p.U[0*POMF_DIM + 2] = 0.0; */
  /* p.U[1*POMF_DIM + 0] = 0.0; */
  /* p.U[1*POMF_DIM + 1] = 1.0; */
  /* p.U[1*POMF_DIM + 2] = 0.0; */
  /* p.U[2*POMF_DIM + 0] = 0.0; */
  /* p.U[2*POMF_DIM + 1] = 0.0; */
  /* p.U[2*POMF_DIM + 2] = 1.0; */

  CHECK_NULL( prob = kfe_problem_from_sde(&sde, &p), OUT_OF_MEM );

  gsl_rng_env_setup();

  CHECK_NULL( rng = gsl_rng_alloc(gsl_rng_default), OUT_OF_MEM );

  CHECK_( kfe_solver_init(&solver, prob, ns) );
  CHECK_( kfe_mcsolver_init(&mcsolver, &sde, n_samples, dt, rng_fun, rng) );

  CHECK_NULL( kfe = kfe_solver_make_KFE_matrix(&solver),     FAILED );
  CHECK_NULL( ident = eye(kfe->n, SPALLOC_CCS),              OUT_OF_MEM );
  CHECK_NULL( bweuler = cs_add(ident, kfe, 1.0, -dt),        OUT_OF_MEM );
  CHECK_NULL( data = malloc(kfe->n*sizeof(*data)),           OUT_OF_MEM );
  CHECK_NULL( data_init = malloc(kfe->n*sizeof(*data_init)), OUT_OF_MEM );
  CHECK_NULL( data_temp = malloc(kfe->n*sizeof(*data_temp)), OUT_OF_MEM );

  CHECK_( kfe_solver_set_init_state(&solver, x_init, data_init) );
  CHECK_( kfe_mcsolver_init_state_pars(&mcsolver, x_init, &p) );

  CHECK_( umfsolver_init(&umfsolver, bweuler) );

  for (iter = 1; iter <= niter; ++iter)
    {
      /* printf("# iter = %d, t = %lf\n", iter, iter*dt); */
      CHECK_( umfsolve_(&umfsolver, bweuler, data_init, data) );
      /* if ((iter - 1) % 10 == 0) */
      /* 	rgrid_print_data(&solver.rgrid, NULL, 1, &data, NULL); */
      /* printf("\n\n"); */

      CHECK_( kfe_mcsolver_step(&mcsolver) );

      int i;
      double s = 0.0;
      for (i = 0; i < kfe->n; ++i) s += data[i];
      fprintf(stderr, "# iter = %d, sum - 1 = %le\n", iter, s - 1.0);

      memcpy(data_init, data, kfe->n*sizeof(*data_init));
    }

  const int axis = 0;
  
  // rgrid_print_data(&solver.rgrid, NULL, 1, &data, NULL);
  {
    rgrid_t mgrid = rgrid_nil;
    double **mdata = NULL;
    rgrid_marginals(&solver.rgrid, NULL, 1, &data, axis, &mgrid, &mdata);
    rgrid_print_data(&mgrid, NULL, 1, mdata, NULL);
    rgrid_destroy(&mgrid);
    Xfree_many(mdata, 1, free);
  }

  printf("\n\n");

  int i;
  for (i = 0; i < solver.rgrid.m; ++i) data[i] = 0.0;
  for (i = 0; i < n_samples; ++i)
    {
      int ix[POMF_DIM];
      CHECK_( rgrid_nearest_i(&solver.rgrid, &mcsolver.data[POMF_DIM*i], ix) );
      int r = rgrid_rindex(&solver.rgrid, ix);
      data[r] += 1.0/n_samples;
    }

  {
    rgrid_t mgrid = rgrid_nil;
    double **mdata = NULL;
    rgrid_marginals(&solver.rgrid, NULL, 1, &data, axis, &mgrid, &mdata);
    rgrid_print_data(&mgrid, NULL, 1, mdata, NULL);
    rgrid_destroy(&mgrid);
    Xfree_many(mdata, 1, free);
  }
  
  // rgrid_print_data(&solver.rgrid, NULL, 1, &data, NULL);

 exit:
  Xfree(data);
  Xfree(data_init);
  Xfree(data_temp);
  Xfree_(kfe, cs_spfree);
  Xfree_(ident, cs_spfree);
  Xfree_(bweuler, cs_spfree);
  Xfree_(rng, gsl_rng_free);
  umfsolver_destroy(&umfsolver);
  kfe_solver_destroy(&solver);
  kfe_mcsolver_destroy(&mcsolver);
  kfe_problem_free(prob);

  return OK;
}
