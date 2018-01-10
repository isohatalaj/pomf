
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pomf_smsimu.h"
#include "pomf_model.h"

int
main()
{
  int i, j;
  int status = OK;

  rgrid_t grid = rgrid_nil;

  cs *KFE_op = NULL, **Jx_ops = NULL;
  double *rhs = NULL;
  double *fx = NULL;
  double **Jxs = NULL;

  double **y = NULL;
  double *fy = NULL;
  double **Jys = NULL;

  pomf_params_t p;

  /* Model parameters --------------------------------------------- */
  const double a_mean = 0.1;
  const double a_stddev = 0.05;
  const double a_revrate = 0.5;
  p.a_bar       = a_mean;
  p.gamma       = a_revrate;
  p.sigma_a     = sqrt(2*a_revrate)*a_stddev;
  p.sigma_c     = 0.15;
  p.sigma_c_bar = 0.15;
  p.kappa       = 0.00;
  p.r_bar       = 0.02;
  p.rho         = 0.07;
  p.rho_bar     = 0.06;

  const int use_cov = 1;

  /* Solution domain ---------------------------------------------- */
  const int ns[POMF_DIM] = {31, 31, 41};
  double bs[2*POMF_DIM];
  
  /* Solver parameters -------------------------------------------- */
  const double solver_tol = 1e-9;
  const int solver_lu_order = 1;
  
  /* Initialization ----------------------------------------------- */

  void (*lfunc)(double bs[POMF_DIM], void *params);
  int (*sfunc)(const double x[POMF_DIM], 
	       double mu_out[POMF_DIM],
	       double sigma_out[POMF_DIM*POMF_WDIM],
	       void *params);
  
  if (use_cov)
    {
      lfunc = &pomf_limitfunc_alt;
      sfunc = &pomf_sdefunc_alt;
    }
  else
    {
      lfunc = &pomf_limitfunc;
      sfunc = &pomf_sdefunc;
    }
  
  lfunc(bs, &p);
  
  status = rgrid_init(&grid, POMF_DIM, ns, bs);
  if (status) goto exit;

  status = pomf_smkfe(&grid,
		      POMF_DIM, POMF_WDIM, sfunc, &p,
		      &KFE_op, &Jx_ops, &rhs);
  if (status) goto exit;

  CHECK_NULL( fx = malloc(grid.m*sizeof(*fx)), OUT_OF_MEM );
  memcpy(fx, rhs, grid.m*sizeof(*fx));

  CHECK( !cs_lusol(solver_lu_order, KFE_op, fx,
		   solver_tol),
	 FAILED );

  CHECK_NULL( Jxs = malloc(POMF_DIM*sizeof(*Jxs)), OUT_OF_MEM );

  for (i = 0; i < POMF_DIM; ++i)
    {
      CHECK_NULL( Jxs[i] = malloc(grid.m*sizeof(*Jxs[i])),
		  OUT_OF_MEM );
      for (j = 0; j < grid.m; ++j) (Jxs[i])[j] = 0.0;
      
      cs_gaxpy(Jx_ops[i], fx, Jxs[i]);
    }

  /* Change of variables */
  if (use_cov)
    {
      CHECK_NULL( y = calloc(POMF_DIM, sizeof(*y)), OUT_OF_MEM );
      CHECK_NULL( fy = malloc(grid.m*sizeof(*fy)), OUT_OF_MEM );
      CHECK_NULL( Jys = calloc(POMF_DIM, sizeof(*Jys)), OUT_OF_MEM );
      
      for (i = 0; i < POMF_DIM; ++i)
	{
	  CHECK_NULL( y[i] = malloc(grid.m*sizeof(*y[i])), OUT_OF_MEM );
	  CHECK_NULL( Jys[i] = malloc(grid.m*sizeof(*Jys[i])), OUT_OF_MEM );
	}
      
      CHECK( pomf_covfj(&grid, fx,
			y, fy, Jys,
			pomf_Eta, NULL,
			POMF_WDIM, pomf_sdefunc_alt, &p),
	     FAILED );
    }
  /* End change of variables */

  printf("# FULL DISTRIBUTION AND PROBABILITY CURRENTS\n");
  double *pJxs[POMF_DIM + 1]; // = {fy, Jys[0], Jys[1], Jys[2]};
  double **nodes;

  if (use_cov)
    {
      nodes = y;
      pJxs[0] = fy;
      pJxs[1] = Jys[0];
      pJxs[2] = Jys[1];
      pJxs[3] = Jys[2];
    }
  else
    {
      nodes = NULL;
      pJxs[0] = fx;
      pJxs[1] = Jxs[0];
      pJxs[2] = Jxs[1];
      pJxs[3] = Jxs[2];
    }
  
  const char *headers[POMF_DIM + 1] = {"pdf", "J0", "J1", "J2"};
  rgrid_print_data(&grid, nodes, POMF_DIM + 1, pJxs, headers);

  /* Print rgrid_marginals */
  for (i = 0; i < POMF_DIM; ++i)
    {
      rgrid_t mgrid = rgrid_nil;
      double **mpJxs = NULL;
      rgrid_marginals(&grid, nodes ? nodes[i] : NULL, POMF_DIM + 1, pJxs, i, &mgrid, &mpJxs);

      double *mnodes[POMF_DIM - 1];
      if (use_cov)
	{
	  if (i == 0) { mnodes[0] = nodes[1]; mnodes[1] = nodes[2]; }
	  if (i == 1) { mnodes[0] = nodes[0]; mnodes[1] = nodes[2]; }
	  if (i == 2) { mnodes[0] = nodes[0]; mnodes[1] = nodes[1]; }
	}

      printf("\n\n# MARGINALS ON AXIS %d\n", i);
      rgrid_print_data(&mgrid, nodes ? mnodes : NULL, POMF_DIM + 1, mpJxs, headers);

      rgrid_destroy(&mgrid);
      Xfree_many(mpJxs, POMF_DIM + 1, free);
    }
  
 exit:
  rgrid_destroy(&grid);

  Xfree_(KFE_op, cs_spfree);
  Xfree(fx);
  Xfree_many(Jx_ops, POMF_DIM, cs_spfree);
  Xfree_many(Jxs, POMF_DIM, free);
  Xfree(fy);
  Xfree_many(Jys, POMF_DIM, free);
  Xfree(rhs);
  Xfree_many(y, POMF_DIM, free);

  return status;
}


