
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pomf_smsimu.h"
#include "pomf_mcsimu.h"
#include "pomf_model.h"

int
main()
{
  int i, j;
  int status = OK;

  rgrid_t grid = rgrid_nil;

  cs *RHS_op = NULL, *LHS_op = NULL, **Jx_ops = NULL;
  double *rhs = NULL;
  double *fx = NULL;
  double **Jxs = NULL;

  double **y = NULL;
  double *fy = NULL;
  double **Jys = NULL;

  pomf_params_t p;

  /* Model parameters --------------------------------------------- */

  const double a_mean = 0.070;
  const double a_stddev = 0.035;
  const double a_revrate = 0.5;
  p.a_bar       = a_mean;
  p.gamma       = a_revrate;
  p.sigma_a     = sqrt(2*a_revrate)*a_stddev;
  p.sigma_c     = 0.025;
  p.sigma_c_bar = 0.025;
  p.kappa       = -0.707;
  //p.kappa       = 0.0;
  p.r_bar       = 0.010;
  p.rho         = 0.070;   /* Note: Model depends only on difference rho - rho_bar */
  p.rho_bar     = 0.050;
  
  /* Solution domain ---------------------------------------------- */
  const int ns[POMF_DIM] = {21, 21, 61};
  double bs[2*POMF_DIM];

  const double t_end = 0.1;
  const int n_step = 1;
  const double dt = t_end / n_step;
  
  /* Initialization ----------------------------------------------- */

  void (*lfunc)(double bs[POMF_DIM], void *params);
  int (*sfunc)(const double x[POMF_DIM], 
	       double mu_out[POMF_DIM],
	       double sigma_out[POMF_DIM*POMF_WDIM],
	       void *params);
  
  pomf_limitfunc(bs, &p);
  
  CHECK( rgrid_init(&grid, POMF_DIM, ns, bs), FAILED );

  fprintf(stderr, "# Init KFE... "); fflush(stderr);

  /* CHECK( pomf_smkfe(&grid, */
  /* 		    POMF_DIM, POMF_WDIM, &pomf_sdefunc, &p, */
  /* 		    &KFE_op, &Jx_ops, NULL), */
  /* 	 FAILED ); */

  CHECK( pomf_smkfe_timedep(&grid,
			    POMF_DIM, POMF_WDIM, &pomf_sdefunc, &p,
			    dt,
			    &LHS_op, &RHS_op, &Jx_ops),
	 FAILED );

  fprintf(stderr, "Done\n"); fflush(stderr);
  
  CHECK_NULL( fx = malloc(grid.m*sizeof(*fx)), OUT_OF_MEM );

  for (i = 0; i < grid.m; ++i) fx[i] = i == grid.m / 2 ? 1.0 : 0.0;

  /* Time-step the solution */

  for (i = 0; i < n_step; ++i)
    {
      fprintf(stderr, "# Time-step %d/%d... ", i + 1, n_step); fflush(stderr);

      CHECK( pomf_smsimu_step(&grid, dt, LHS_op, RHS_op, fx), FAILED );

      fprintf(stderr, "Done\n"); fflush(stderr);
    }
  
  /**/

  CHECK_NULL( Jxs = calloc(POMF_DIM, sizeof(*Jxs)), OUT_OF_MEM );

  for (i = 0; i < POMF_DIM; ++i)
    {
      CHECK_NULL( Jxs[i] = malloc(grid.m*sizeof(*Jxs[i])),
		  OUT_OF_MEM );
      for (j = 0; j < grid.m; ++j) (Jxs[i])[j] = 0.0;
      
      cs_gaxpy(Jx_ops[i], fx, Jxs[i]);
    }

  printf("# FULL DISTRIBUTION AND PROBABILITY CURRENTS\n");
  double *pJxs[POMF_DIM + 1]; // = {fy, Jys[0], Jys[1], Jys[2]};

  pJxs[0] = fx;
  pJxs[1] = Jxs[0];
  pJxs[2] = Jxs[1];
  pJxs[3] = Jxs[2];
  
  const char *headers[POMF_DIM + 1] = {"pdf", "J0", "J1", "J2"};
  rgrid_print_data(&grid, NULL, POMF_DIM + 1, pJxs, headers);

  /* Print rgrid_marginals */
  for (i = 0; i < POMF_DIM; ++i)
    {
      rgrid_t mgrid = rgrid_nil;
      double **mpJxs = NULL;
      rgrid_marginals(&grid, NULL, POMF_DIM + 1, pJxs, i, &mgrid, &mpJxs);

      double *mnodes[POMF_DIM - 1];

      printf("\n\n# MARGINALS ON AXIS %d\n", i);
      rgrid_print_data(&mgrid, NULL, POMF_DIM + 1, mpJxs, headers);

      rgrid_destroy(&mgrid);
      Xfree_many(mpJxs, POMF_DIM + 1, free);
    }
  
 exit:
  rgrid_destroy(&grid);

  Xfree_(LHS_op, cs_spfree);
  Xfree_(RHS_op, cs_spfree);
  Xfree(fx);
  Xfree_many(Jx_ops, POMF_DIM, cs_spfree);
  Xfree_many(Jxs, POMF_DIM, free);
  Xfree(fy);
  Xfree_many(Jys, POMF_DIM, free);
  Xfree(rhs);
  Xfree_many(y, POMF_DIM, free);

  return status;
}



