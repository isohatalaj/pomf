
#include <string.h>

#include "pomf_smsimu.h"

/**
 * Helper function: Given a `dim` times `wdim` matrix `s` with
 * elements as the amplitudes of independent normal distributed random
 * variables, compute the correlation matrix `c`, `c = s * st / 2`
 * where `st` is the transpose of `s`.
 */
static void
ccovar(const int dim, const int wdim, const double *s, double *c)
{
  int i, j, k;

  for (i = 0; i < dim; ++i)
    {
      for (j = i; j < dim; ++j)
	{
	  double d = 0.0;

	  for (k = 0; k < wdim; ++k)
	    {
	      d += s[wdim*i + k] * s[wdim*j + k];
	    }

	  c[dim*j + i] = c[dim*i + j] = d;
	}
    }
}

int
pomf_smkfe(const rgrid_t *grid,
	   int dim, int wdim,
	   int (*sdefunc)(const double *x,
			  double *mu, double *sigma,
			  void *params),
	   void *params,
	   cs **KFE_out,
	   cs ***Js_out,
	   double **b_out)
{
  int i, j;
  int status = OK;
  double *work = NULL;
  int *iwork = NULL;
  const int m = grid->m;

  /* Sparse matrix allocated sizes could perhaps be tweaked, nzmax's
   * are pretty much just guessed here. */

  /* Drift and diffusion data, allocating `dim` elements for drift and
   * `dim*dim` for the diffusion matrices. Need one of such for each
   * datapoint in the grid. */
  const int ystride = dim*(1 + dim);
  double *y = NULL;  /* Vector of drift and diffusion data */
  double *b = NULL;  /* To hold the KF equation rhs vector, once we're
			done */

  cs *KFE = NULL;  /* To hold the KF equation lhs operator, once we're
		      done */
  cs **Js = NULL;  /* To hold the prob. current operators, once we're
		      done */
  cs **Ks = NULL;  /* Used in constructing boundary conditions
		      (selector matrices for rows of Js to be injected
		      into KFE) */
  cs **Ds = NULL;  /* Derivative operators, one for each dimension. */

  /* Temporary matrices */
  cs *T = NULL, *S = NULL, *U = NULL;

  /* Temporary storage for evaluating drift and diffusion data.  */
  const int work_size = dim + dim + dim*wdim + dim*dim;

  if (grid->dim != dim) { status = INVALID; goto exit; }
  
  CHECK_NULL( y = malloc(ystride*m*sizeof(*y)),       OUT_OF_MEM );
  CHECK_NULL( b = malloc(m*sizeof(*b)),               OUT_OF_MEM );
  CHECK_NULL( Js = calloc(dim, sizeof(*Js)),          OUT_OF_MEM );
  CHECK_NULL( Ks = calloc(dim, sizeof(*Ks)),          OUT_OF_MEM );
  CHECK_NULL( Ds = calloc(dim, sizeof(*Ds)),          OUT_OF_MEM );
  CHECK_NULL( iwork = malloc(dim*sizeof(*iwork)),     OUT_OF_MEM );
  CHECK_NULL( work = malloc(work_size*sizeof(*work)), OUT_OF_MEM );
  
  double *x = work;
  double *mu = x + dim;
  double *sigma = mu + dim;
  double *diffuse = sigma + dim*wdim;

  /* Allocate derivative operators */
  for (i = 0; i < dim; ++i)
    CHECK_NULL( Ds[i] = deriv(grid, i, SPALLOC_CCS), OUT_OF_MEM );

  /* Precompute the drift diffusion data. */
  i = 0;
  rgrid_start(grid, iwork, x);
  do
    {
      double *ddi = y + i*ystride;

      CHECK( sdefunc(x, mu, sigma, params), FAILED );

      ccovar(dim, wdim, sigma, diffuse);

      memcpy(ddi, mu, dim*sizeof(*ddi));
      memcpy(ddi + dim, diffuse, dim*dim*sizeof(*ddi));

      ++i;
    }
  while (rgrid_step(grid, iwork, x));
  
  /* First set the current operators to just the drift component,
   * adding the diffusion terms in subsequent loops. */
  for (i = 0; i < dim; ++i)
    CHECK_NULL( Js[i] = diag(m, y + i, ystride, SPALLOC_CCS), FAILED );


  /* Now add the diffusion terms. Here using sparse add/multiplies to
   * do some of the dirty work, but really, we could just splice in
   * the routines for constructing differentiation matrices, and
   * prolly get better performance. Now just trying to write something
   * that works and is easy to debug! */
  for (i = 0; i < dim; ++i)
    {
      for (j = 0; j < dim; ++j)
	{
	  CHECK_NULL( S = diag(m, y + dim + i*dim + j,
			       ystride, SPALLOC_CCS),   OUT_OF_MEM );
	  CHECK_NULL( T = cs_multiply(Ds[j], S),        FAILED );
	  CHECK_NULL( U = cs_add(Js[i], T, 1.0, -0.5),  FAILED );

	  cs_spfree(S); S = NULL;
	  cs_spfree(T); T = NULL;
	  cs_spfree(Js[i]); Js[i] = U; U = NULL;
	}
    }

  /* The current operators are ready, now take the divergence of the
   * vector field that is their range. */
  CHECK_NULL( KFE = cs_multiply(Ds[0], Js[0]), FAILED );
  
  for (i = 1; i < dim; ++i)
    {
      CHECK_NULL( S = cs_multiply(Ds[i], Js[i]),  FAILED );
      CHECK_NULL( T = cs_add(KFE, S, 1.0, 1.0),   FAILED );

      cs_spfree(S); S = NULL;
      cs_spfree(KFE); KFE = T; T = NULL;
    }
  
  /* OK, now do the boundary condition overrides. */
  
  /* Ks[i] selects the rows where Js[i] is included. T will be a
   * diagonal matrix with zero elements on rows that are excluded from
   * KFE in favour of a boundary condition, other diagonal elements
   * are unity => Left-multiply by T zeros rows corresponding to
   * boundary points; adding to T*KFE the Ks[i] weighted Js[i]
   * matrices fills these rows with the appropriate boundary
   * conditions. */
  CHECK_NULL( T = cs_spalloc(m, m, 4*m,
			     SPALLOC_VALUES, SPALLOC_TRIPLET),
	      OUT_OF_MEM );

  for (i = 0; i < dim; ++i)
    CHECK_NULL( Ks[i] = cs_spalloc(m, m, 4*m,
				   SPALLOC_VALUES, SPALLOC_TRIPLET),
		OUT_OF_MEM );


  /* Compute the index for the normalisation condition. */
  int m_norm = 0; // m / 2;
  
  j = 0;
  rgrid_start(grid, iwork, NULL);
  do
    {
      int n_boundary = 0;
      for (i = 0; i < dim; ++i)
	{
	  if (iwork[i] == 0)
	    CHECK( !cs_entry(Ks[i], j, j, +1.0), OUT_OF_MEM ); 
	  else if (iwork[i] == grid->ns[i] - 1)
	    CHECK( !cs_entry(Ks[i], j, j, -1.0), OUT_OF_MEM );
	  else
	    continue;

	  n_boundary++;
	}

      if (n_boundary)
	{
	  /* Don't override boundary condition row with the
	     normalisation condition, defer to next row. Does this
	     matter? Perhaps not, but let's just be safe. */

	  if (j == m_norm) m_norm++;
	}
      else if (j != m_norm)
	{
	  /* Rows that impose boundaries or normalisation are left to
	   * zero in matrix T, all others are unity. */

	  CHECK( !cs_entry(T, j, j, 1.0), OUT_OF_MEM );
	}
      
      ++j;
    }
  while (rgrid_step(grid, iwork, NULL));

  CHECK( compress_(&T), OUT_OF_MEM );

  for (i = 0; i < dim; ++i)
    CHECK( compress_(&Ks[i]), OUT_OF_MEM );

  CHECK_NULL( S = cs_multiply(T, KFE), OUT_OF_MEM );

  cs_spfree(KFE); KFE = S; S = NULL;
  cs_spfree(T); T = NULL;

  /* Now KFE has boundary rows zeroed. Now add in the boundary
   * conditions to fill these rows.  */

  /* Start with normalisation */
  CHECK_NULL( T = cs_spalloc(m, m, m,
			     SPALLOC_VALUES, SPALLOC_TRIPLET),
	      OUT_OF_MEM );

  double dV = 1.0;
  for (i = 0; i < dim; ++i) dV *= grid->hs[i];

  for (i = 0; i < m; ++i)
    {
      b[i] = i == m_norm ? 1.0 : 0.0;
      CHECK( !cs_entry(T, m_norm, i, dV), FAILED );
    }

  CHECK( compress_(&T), OUT_OF_MEM );
  CHECK_NULL( S = cs_add(KFE, T, 1.0, 1.0), FAILED );

  cs_spfree(KFE); KFE = S; S = NULL;
  cs_spfree(T); T = NULL;

  /* Now boundary currents */
  for (i = 0; i < dim; ++i)
    {
      CHECK_NULL( T = cs_multiply(Ks[i], Js[i]), FAILED );
      CHECK_NULL( S = cs_add(KFE, T, 1.0, 1.0),  FAILED );

      cs_spfree(KFE); KFE = S; S = NULL;
      cs_spfree(T); T = NULL;
    }

 exit:
  Xfree(work);
  Xfree(iwork);
  Xfree(y);
  Xfree_(S, cs_spfree);
  Xfree_(T, cs_spfree);
  Xfree_(U, cs_spfree);

  Xfree_many(Ds, dim, cs_spfree);
  Xfree_many(Ks, dim, cs_spfree);

  if (status == OK)
    {
      *KFE_out = KFE;
      *Js_out = Js;
      *b_out = b;
    }
  else
    {
      *KFE_out = NULL;
      *Js_out = NULL;
      *b_out = NULL;

      Xfree_(KFE, cs_spfree);
      Xfree_many(Js, dim, cs_spfree);
      Xfree(b);
    }

  return status;
}

int
pomf_covfj(rgrid_t *grid,
	   const double *fx, 
	   double **y_out, double *fy_out, double **Jy_out,
	   int (*psi)(const double *x, double *y,
		      double *jacobian,
		      double *inverse_jacobian,
		      double *det_jacobian,
		      double *hessians,
		      void *psi_params),
	   void *psi_params,
	   int wdim,
	   int (*musigma)(const double *x,
			  double *mu, double *sigma,
			  void *musigma_params),
	   void *musigma_params)
{
  int status = OK;
  int i, j, k, l;
  const int dim = grid->dim;
  const int dim2 = dim*dim;
  const int m = grid->m;

  const int wstride = dim2 /* inv_jac */ + dim /* mu_y */ + dim2 /* diffusion_y */;
  const int w_inv_jac = 0;
  const int w_mu_y = dim2;
  const int w_dfs_y = dim2 + dim;
  double *work = NULL;

  double *x = NULL, *y = NULL, *jac = NULL, *inv_jac = NULL, *hess = NULL;
  double *mu = NULL, *sigma = NULL, *diffusion = NULL;
  double *mu_y = NULL, *diffusion_y = NULL;
  double det_jac;
  int *mi = NULL;

  cs **Dxs = NULL, **Dys = NULL, **Jys = NULL;
  cs *S = NULL, *T = NULL, *U = NULL;

  CHECK_NULL( work = malloc(wstride*m*sizeof(*work)),  OUT_OF_MEM );
  CHECK_NULL( mi = malloc(dim*sizeof(*mi)),            OUT_OF_MEM );
  CHECK_NULL( x = malloc(dim*sizeof(*x)),              OUT_OF_MEM );
  CHECK_NULL( y = malloc(dim*sizeof(*y)),              OUT_OF_MEM );
  CHECK_NULL( jac = malloc(dim2*sizeof(*jac)),         OUT_OF_MEM );
  CHECK_NULL( inv_jac = malloc(dim2*sizeof(*inv_jac)), OUT_OF_MEM );
  CHECK_NULL( hess = malloc(dim2*dim*sizeof(*hess)),   OUT_OF_MEM );
  CHECK_NULL( mu = malloc(dim*sizeof(*mu)),                OUT_OF_MEM );
  CHECK_NULL( sigma = malloc(dim*wdim*sizeof(*sigma)),     OUT_OF_MEM );
  CHECK_NULL( diffusion = malloc(dim2*sizeof(*diffusion)), OUT_OF_MEM );
  CHECK_NULL( mu_y = malloc(dim*sizeof(*mu_y)),                OUT_OF_MEM );
  CHECK_NULL( diffusion_y = malloc(dim2*sizeof(*diffusion_y)), OUT_OF_MEM );

  CHECK_NULL( Dxs = calloc(dim, sizeof(*Dxs)), OUT_OF_MEM );
  CHECK_NULL( Dys = calloc(dim, sizeof(*Dys)), OUT_OF_MEM ); 
  CHECK_NULL( Jys = calloc(dim, sizeof(*Jys)), OUT_OF_MEM ); 
  
  l = 0;
  rgrid_start(grid, mi, x);
  do
    {
      CHECK( psi(x, y, jac, inv_jac, &det_jac, hess, psi_params), FAILED );
      CHECK( musigma(x, mu, sigma, musigma_params),               FAILED );
      ccovar(dim, wdim, sigma, diffusion);

      const double fy = fx[l] / det_jac;

      /* density */

      fy_out[l] = fy;

      /* y-coordinate */

      for (i = 0; i < dim; ++i) (y_out[i])[l] = y[i];

      /* y-drift */
      /* Ito: mu_y_i = mu_x_j do_x_j f_i + D_jk do_x_kj f_i */
      for (i = 0; i < dim; ++i) 
	{
	  double *hi = hess + i*dim2;
	  double sum = 0.0; 

	  for (j = 0; j < dim; ++j)
	    {
	      sum += mu[j]*jac[i*dim + j];
	      for (k = 0; k < dim; ++k)
		sum += 0.5*diffusion[k*dim + j]*hi[k*dim + j];
	    }

	  mu_y[i] = sum;
	}

      /* y-diffusion */
      /* Ito: sigma_y_ij = sigma_kj do_x_k f_i  
	      D_y_ij = 1/2 sigma_ik * sigma_jk
	             = 1/2 sigma_k1k do_x_k1 f_i sigma_k2k do_x_k2 f_j
	             = do_x_k1 f_i (1/2 sigma_k1k sigma_k2k) do_x_k2 f_j
	             = do_x_k1 f_i D_x_k1k2 do_x_k2 f_j
       */
      
      for (i = 0; i < dim; ++i)
	for (j = 0; j < dim; ++j)
	  {
	    double sum = 0.0;
	    int k1, k2;

	    for (k1 = 0; k1 < dim; ++k1)
	      for (k2 = 0; k2 < dim; ++k2)
		{
		  sum += jac[i*dim + k1] * diffusion[k1*dim + k2] * jac[j*dim + k2];
		}

	    diffusion_y[i*dim + j] = sum;
	  }
      
      memcpy(work + l*wstride + w_inv_jac, inv_jac, dim2*sizeof(*work));
      memcpy(work + l*wstride + w_mu_y, mu_y, dim*sizeof(*work));
      memcpy(work + l*wstride + w_dfs_y, diffusion_y, dim2*sizeof(*work));
      
      l++;
    }
  while (rgrid_step(grid, mi, x));

  /* x derivatives */
  for (i = 0; i < dim; ++i)
    CHECK_NULL( Dxs[i] = deriv(grid, i, SPALLOC_CCS), OUT_OF_MEM );

  /* y derivatives */
  for (i = 0; i < dim; ++i)
    {
      CHECK_NULL( S = diag(m, work + i*dim + 0, wstride, SPALLOC_CCS), FAILED);
      CHECK_NULL( Dys[i] = cs_multiply(S, Dxs[i]), FAILED );
      
      for (j = 1; j < dim; ++j)
	{
	  const double *inv_jac = work + w_inv_jac;

	  CHECK_NULL( S = diag(m, inv_jac + i*dim + j, wstride, SPALLOC_CCS), FAILED);
	  CHECK_NULL( T = cs_multiply(S, Dxs[i]), FAILED );
	  CHECK_NULL( U = cs_add(Dys[i], T, 1.0, 1.0), FAILED );

	  cs_spfree(S); S = NULL;
	  cs_spfree(T); T = NULL;
	  cs_spfree(Dys[i]); Dys[i] = U; U = NULL;
	}
    }

  /* y prob currents, drift part */
  for (i = 0; i < dim; ++i)
    CHECK_NULL( Jys[i] = diag(m, work + w_mu_y + i, wstride, SPALLOC_CCS), FAILED );

  /* y prob currents, diffusion part */
  for (i = 0; i < dim; ++i)
    {
      for (j = 0; j < dim; ++j)
	{
	  CHECK_NULL( S = diag(m, work + w_dfs_y + i*dim + j,
			       wstride, SPALLOC_CCS),   OUT_OF_MEM );
	  CHECK_NULL( T = cs_multiply(Dys[j], S),       FAILED );
	  CHECK_NULL( U = cs_add(Jys[i], T, 1.0, -0.5), FAILED );

	  cs_spfree(S); S = NULL;
	  cs_spfree(T); T = NULL;
	  cs_spfree(Jys[i]); Jys[i] = U; U = NULL;
	}
    }
  
  /* apply the operators */
  for (i = 0; i < dim; ++i)
    {
      double *Jyi = Jy_out[i];
      for (j = 0; j < m; ++j) Jyi[j] = 0.0;

      cs_gaxpy(Jys[i], fy_out, Jyi);
    }

 exit:
  Xfree_(S, cs_spfree);
  Xfree_(T, cs_spfree);
  Xfree_(U, cs_spfree);
  Xfree_many(Dxs, dim, cs_spfree);
  Xfree_many(Dys, dim, cs_spfree);
  Xfree_many(Jys, dim, cs_spfree);
  Xfree(diffusion_y);
  Xfree(mu_y);
  Xfree(diffusion);
  Xfree(sigma);
  Xfree(mu);
  Xfree(jac);
  Xfree(inv_jac);
  Xfree(hess);
  Xfree(x);
  Xfree(y);
  Xfree(mi);
  Xfree(work);

  return status;
}
