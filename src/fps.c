
#include "fps.h"
#include <string.h>
#include <math.h>
#include <float.h>

#define TOLFCT 1

const kfe_sde_t kfe_sde_nil = KFE_SDE_NIL;
const kfe_mcsolver_t kfe_mcsolver_nil = KFE_MCSOLVER_NIL;


int
kfe_mcsolver_init(kfe_mcsolver_t *self,
		  const kfe_sde_t *sde,
		  int n_samples, double dt,
		  kfe_nrng_t *rng,
		  void *rng_state)
{
  int status = OK;
  *self = kfe_mcsolver_nil;

  CHECK_NULL( self->data = malloc(sde->n*n_samples*sizeof(*self->data)),
	      OUT_OF_MEM );
  CHECK_NULL( self->bounds = malloc(2*sde->n*sizeof(*self->bounds)),
	      OUT_OF_MEM );
  memcpy(&self->sde, sde, sizeof(self->sde));
  self->n_samples = n_samples;
  self->dt = dt;
  self->t = 0;
  self->params = NULL;
  self->rng = rng;
  self->rng_state = rng_state;

 exit:
  if (status) kfe_mcsolver_destroy(self);

  return status;
}

void
kfe_mcsolver_destroy(kfe_mcsolver_t *self)
{
  if (self == NULL) return;
  Xfree(self->data);
  Xfree(self->bounds);
}


int
kfe_mcsolver_init_state_pars(kfe_mcsolver_t *self,
			     const double *init_state,
			     void *params)
{
  int status = OK;
  int i;

  for (i = 0; i < self->n_samples; ++i)
    {
      memcpy(&self->data[i*self->sde.n], init_state,
	     self->sde.n*sizeof(*self->data));
    }

  CHECK_( self->sde.eval_bounds(self->bounds, params) );

  self->params = params;
  self->t = 0.0;

 exit:

  return status;
}

int
kfe_mcsolver_step(kfe_mcsolver_t *self)
{
  int status = OK;
  const int n_samples = self->n_samples;
  const int n_x = self->sde.n;
  const int n_w = self->sde.k;

  const double dt = self->dt;
  const double sqrt_dt = sqrt(dt);

  double *dw = NULL, *b = NULL, *s = NULL, *dsdx = NULL;

  CHECK_NULL( dw = malloc(n_w*sizeof(*dw)), OUT_OF_MEM );
  CHECK_NULL( b = malloc(n_x*sizeof(*b)), OUT_OF_MEM );
  CHECK_NULL( s = malloc(n_w*n_x*sizeof(*s)), OUT_OF_MEM );
  CHECK_NULL( dsdx = malloc(n_w*n_x*n_x*sizeof(*dsdx)), OUT_OF_MEM );

  int i, j, k;
  for (i = 0; i < n_samples; ++i)
    {
      double *x = &self->data[i*n_x];

      for (j = 0; j < n_w; ++j)
	{
	  dw[j] = self->rng(0.0, sqrt_dt, self->rng_state);
	}

      CHECK_( self->sde.eval(x, b, s, dsdx, self->params) );

      for (j = 0; j < n_x; ++j)
      	{
      	  double dx = b[j]*dt;

      	  for (k = 0; k < n_w; ++k)
      	    {
      	      dx += s[j*n_w + k]*dw[k];
      	    }

      	  /* x' = x + dx
      	   * If x' > x1, then x'' = x1 - (x' - x1) = 2*x1 - x'
      	   * If x' < x0, then x'' = x0 + (x0 - x') = 2*x0 - x'
      	   */

      	  x[j] += dx;

      	  const double x_lo = self->bounds[2*j + 0];
      	  const double x_hi = self->bounds[2*j + 1];

      	  if (x[j] < x_lo)
      	    {
      	      x[j] = 2*x_lo - x[j];
      	    }
      	  else if (x[j] > x_hi)
      	    {
      	      x[j] = 2*x_hi - x[j];
      	    }
      	}

    }

  self->t += dt;

 exit:
  Xfree(dw);
  Xfree(b);
  Xfree(s);
  Xfree(dsdx);

  return status;
}


int kfe_gamma_size(int n) { return 2*n*(n-1); }
int kfe_gamma_stride(int n) { return n*(n-1)/2; }

int
kfe_gamma_index(int n, int l, int i, int j)
{
  const int stride = kfe_gamma_stride(n);
  if (i > j)
    {
      if (l == 1) { l = 2; }
      else if (l == 2) { l = 1; }
      int temp = i;
      i = j;
      j = temp;
    }

  return l*stride + i*(n-1) - (i-1)*i/2 + j;  
}

const kfe_solver_t kfe_solver_nil = KFE_SOLVER_NIL;

void
kfe_solver_destroy(kfe_solver_t *self)
{
  Xfree(self->gamma);
  Xfree(self->delta);
  Xfree(self->iwork);
  Xfree(self->work);
  rgrid_destroy(&self->rgrid);
}

int
kfe_solver_at_bound(kfe_solver_t *self,
		    const double *x,
		    int i)
{
  if (x[i] + 0.25*self->rgrid.hs[i] > self->rgrid.bs[2*i + 1])
    return 1;

  if (x[i] - 0.25*self->rgrid.hs[i] < self->rgrid.bs[2*i + 0])
    return -1;

  return 0;
}

int
kfe_solver_eval_D_(kfe_solver_t *self,
		   const double *x,
		   int i, double hi,
		   int j, double hj,
		   double *D,
		   double *x_temp)
{
  const int n_dim = self->prob->n;
  memcpy(x_temp, x, n_dim*sizeof(*x_temp));
  x_temp[i] += hi;
  x_temp[j] += hj;

  int k;
  for (k = 0; k < n_dim; ++k)
    {
      if (x_temp[k] < self->prob->bounds[2*k + 0] - self->rgrid.hs[k]/4.0 ||
      	  x_temp[k] > self->prob->bounds[2*k + 1] + self->rgrid.hs[k]/4.0)
  	{
  	  for (k = 0; k < n_dim*n_dim; ++k) D[k] = 0.0;
  	  return OK;
  	}
    }

  return self->prob->D(x_temp, D, self->prob->params);
}

int
kfe_solver_eval_K_(kfe_solver_t *self,
		   const double *x,
		   int i, double hi,
		   double *K,
		   double *x_temp)
{
  const int n_dim = self->prob->n;
  memcpy(x_temp, x, self->prob->n*sizeof(*x_temp));
  x_temp[i] += hi;

  int j;
  for (j = 0; j < n_dim; ++j)
    {
      if (x_temp[j] < self->prob->bounds[2*j + 0] - self->rgrid.hs[j]/4.0 ||
      	  x_temp[j] > self->prob->bounds[2*j + 1] + self->rgrid.hs[j]/4.0)
	{
	  for (j = 0; j < n_dim; ++j) K[j] = 0.0;
	  return OK;
	}
    }
  

  return self->prob->K(x_temp, K, self->prob->params);
}

int
kfe_solver_eval_delta_std(kfe_solver_t *self,
			  int i,
			  const double *x,
			  int h_axis, double h_len,
			  double *delta_i)
{
  if (FORCE_05_DELTA)
    {
      *delta_i = 0.5;
      return OK;
    }
  
  int status = OK;
  const int n_dim = self->prob->n;

  const int may_fail = MAY_FAIL;

  double *D = NULL, *K = NULL, *x_local = NULL, *x_temp = NULL;

  CHECK_NULL( D = malloc(n_dim*n_dim*sizeof(*D)), OUT_OF_MEM );
  CHECK_NULL( K = malloc(n_dim*sizeof(*K)), OUT_OF_MEM );
  CHECK_NULL( x_local = malloc(n_dim*sizeof(*x_local)), OUT_OF_MEM );
  CHECK_NULL( x_temp = malloc(n_dim*sizeof(*x_temp)), OUT_OF_MEM );
  
  memcpy(x_local, x, n_dim*sizeof(*x_local));
  x_local[h_axis] += h_len;

  /* if (kfe_solver_at_bound(self, x_local, i)) */
  /*   { */
  /*     *delta_i = 1.0; */
  /*     goto exit; */
  /*   } */
  
  CHECK_( kfe_solver_eval_K_(self, x_local, i, 0.5*self->rgrid.hs[i], K, x_temp) );
  const double Ki = K[i];
  
  CHECK_( kfe_solver_eval_D_(self, x_local, i, 0.5*self->rgrid.hs[i], i, 0.0, D, x_temp) );
  const double Dii = D[i*n_dim + i];

  double Ca, Cb;
  Ca = Cb = Dii/self->rgrid.hs[i]; 
  int j;
  for (j = 0; j < n_dim; ++j)
    {
      if (j == i) continue;
      int k = i*n_dim + j;
      
      double gamma[4], gamma_mj[4];
      
      CHECK_( kfe_solver_eval_gamma_std(self, i, j, x_local,
					0, 0.0, 0, 0.0,
					gamma) );
      CHECK_( kfe_solver_eval_gamma_std(self, i, j, x_local,
					j, -self->rgrid.hs[j], 0, 0.0,
					gamma_mj) );
      
      CHECK_( kfe_solver_eval_D_(self, x_local,
				 i, 0.5*self->rgrid.hs[i],
				 j, 0.0,
				 D, x_temp) );
      const double Dij_pi2 = D[k];
      
      CHECK_( kfe_solver_eval_D_(self, x_local,
				 i, 0.0,
				 j, 0.5*self->rgrid.hs[j],
				 D, x_temp) );
      const double Dij_pj2 = D[k];
      
      CHECK_( kfe_solver_eval_D_(self, x_local,
				 i, 0.0,
				 j, -0.5*self->rgrid.hs[j],
				 D, x_temp) );
      const double Dij_mj2 = D[k];
      
      CHECK_( kfe_solver_eval_D_(self, x_local,
				 i, self->rgrid.hs[i],
				 j, 0.5*self->rgrid.hs[j],
				 D, x_temp) );
      const double Dij_pi_pj2 = D[k];
      
      CHECK_( kfe_solver_eval_D_(self, x_local,
				 i, self->rgrid.hs[i],
				 j, -0.5*self->rgrid.hs[j],
				 D, x_temp) );
      const double Dij_pi_mj2 = D[k];

      
      Ca += +((Dij_pi2 + Dij_pj2)*gamma[1] - (Dij_pi2 + Dij_mj2)*gamma_mj[3])
	/ self->rgrid.hs[j];
      Cb += -((Dij_pi2 + Dij_pi_pj2)*gamma[0] - (Dij_pi2 + Dij_pi_mj2)*gamma_mj[2])
	/ self->rgrid.hs[j];
    }
  
  double Kabs = fabs(Ki);
  double delta;
  
  if (Ki < 0.0)
    {
      const double dhi = MIN(1.0 + Ca/Kabs, Cb/Kabs);
      if (dhi < 0.0)
	{
	  ERRMSG("FAILED EVALUATING delta_i, i = %d, reason:\n", i);
	  ERRMSG("Cannot satisfy DIAG. trans. rate +ivity: delta upper bnd. -ive\n");
	  ERRMSG("dhi = %lf\n", dhi);	  
	  if (may_fail) FAILWITH(FAILED);

	  if (ALLOW_NEGTRANSR)
	    delta = 0.0;
	  else
	    delta = dhi;
	}
      else
	{
	  /* mean */
	  // delta = MIN(dhi, 1.0) / 2.0; 

	  /* aggressive */
	  delta = MIN(0.50, dhi);
	}
    }
  else if (Ki > 0.0)
    {
      const double dlo = MAX(1.0 - Ca/Kabs, -Cb/Kabs);
      if (dlo > 1.0)
	{
	  ERRMSG("FAILED EVALUATING delta_i, i = %d, reason:\n", i);
	  ERRMSG("Cannot satisfy DIAG. trans. rate +ivity: delta lower bnd > 1.0\n");
	  ERRMSG("dlo = %lf\n", dlo);
	  ERRMSG("Ca = %lf, Cb = %lf\n", Ca, Cb);
	  int l;
	  for (l = 0; l < n_dim; ++l)
	    ERRMSG("x[%d] = %lf\n", l, x_local[l]);
	  if (may_fail) FAILWITH(FAILED);

	  if (ALLOW_NEGTRANSR)
	    delta = 1.0;
	  else
	    delta = dlo;
	}
      else
	{
	  /* mean */
	  // delta = (1.0 + MAX(dlo, 0.0)) / 2.0;
	  
	  /* aggressive */
	  delta = MAX(0.5, dlo);
	  // delta = MAX(0.5, delta);
	}

      // fprintf(stderr, "# delta[%d] = %lf (K > 0), dlo = %lf\n", i, delta, dlo);
      // delta = MAX(dlo, 0.0);
    }
  else
    {
      if (Ca < 0.0 || Cb < 0.0)
	{
	  ERRMSG("FAILED EVALUATING delta_i, i = %d, reason:\n", i);
	  ERRMSG("Cannot satisfy DIAG. trans. rate +ivity: K == 0 and Ca or Cb < 0.0\n");
	  ERRMSG("Ca = %lf, Cb = %lf\n", Ca, Cb);
	  int l;
	  for (l = 0; l < n_dim; ++l)
	    ERRMSG("x[%d] = %lf\n", l, x_local[l]);
	  
	  if (may_fail) FAILWITH(FAILED);
	}
      delta = 0.5;
    }

  // CHECK( delta < 0.0 || delta > 1.0, FAILED );
  *delta_i = delta;

  
 exit:
  Xfree(D);
  Xfree(K);
  Xfree(x_temp);
  Xfree(x_local);
  
  return status;
}

int
kfe_solver_eval_gamma_std(kfe_solver_t *self,
			  int i, int j,
			  const double *x,
			  int h_axis, double h_len,
			  int h_axis_2, double h_len_2,
			  double gamma_ij[4])
{
  int status = OK;
  double Dd0, Dd1, Dd2, Dd3;
  const int n_dim = self->prob->n;
  int k = i*n_dim + j;

  double *D = NULL, *x_local = NULL, *x_temp = NULL;

  CHECK_NULL( D = malloc(n_dim*n_dim*sizeof(*D)), OUT_OF_MEM );
  CHECK_NULL( x_local = malloc(n_dim*sizeof(*x_local)), OUT_OF_MEM );
  CHECK_NULL( x_temp = malloc(n_dim*sizeof(*x_temp)), OUT_OF_MEM );
  
  memcpy(x_local, x, n_dim*sizeof(*x_local));
  x_local[h_axis] += h_len;
  x_local[h_axis_2] += h_len_2;
  
  /* Gamma 3 */
  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, 0.5*self->rgrid.hs[i],
			     j, 0.0,
			     D, x_temp) );
  Dd3 = D[k];
  
  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, 0.0,
			     j, 0.5*self->rgrid.hs[j],
			     D, x_temp) );
  Dd3 += D[k];
  
  /* Gamma 0 */
  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, 0.5*self->rgrid.hs[i],
			     j, 1.0*self->rgrid.hs[j],
			     D, x_temp) );
  Dd0 = D[k];
  
  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, 1.0*self->rgrid.hs[i],
			     j, 0.5*self->rgrid.hs[j],
			     D, x_temp) );
  Dd0 += D[k];
  
  
  /* Gamma 1 */
  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, 0.5*self->rgrid.hs[i],
			     j, 1.0*self->rgrid.hs[j],
			     D, x_temp) );
  Dd1 = D[k];
  
  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, 0.0,
			     j, 0.5*self->rgrid.hs[j],
			     D, x_temp) );
  Dd1 += D[k];
  
  
  /* Gamma 2 */
  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, 1.0*self->rgrid.hs[i],
			     j, 0.5*self->rgrid.hs[j],
			     D, x_temp) );
  Dd2 = D[k];
  
  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, 0.5*self->rgrid.hs[i],
			     j, 0.0,
			     D, x_temp) );
  Dd2 += D[k];
  
  double gamma[4];
  gamma[0] = Dd0 >= 0.0 ? 1.0 : 0.0; gamma[1] = Dd1 <= 0.0 ? 1.0 : 0.0;
  gamma[3] = Dd3 >= 0.0 ? 1.0 : 0.0; gamma[2] = Dd2 <= 0.0 ? 1.0 : 0.0;
  /* gamma[0] = MAX(Dd0, 0.0); gamma[1] = -MIN(Dd1, 0.0); */
  /* gamma[3] = MAX(Dd3, 0.0); gamma[2] = -MIN(Dd2, 0.0); */

  int bi = kfe_solver_at_bound(self, x_local, i);
  int bj = kfe_solver_at_bound(self, x_local, j);
  
  if (bi == 1)
    {
      gamma[1] = 0.0; gamma[3] = 0.0;
    }
  if (bj == 1)
    {
      gamma[2] = 0.0; gamma[3] = 0.0;
    }
  
  double norm = gamma[0] + gamma[1] + gamma[2] + gamma[3];
  
  if (norm <= 0.0)
    {
      ERRMSG("Cannot satisfy OFF-DIAG. trans. rate +ivity\n");
      ERRMSG("Debug info: \n");
      FAILWITH(FAILED);
    }
  
  gamma_ij[0] = gamma[0] / norm;
  gamma_ij[1] = gamma[1] / norm;
  gamma_ij[2] = gamma[2] / norm;
  gamma_ij[3] = gamma[3] / norm;

 exit:
  Xfree(D);
  Xfree(x_temp);
  Xfree(x_local);
  
  return status;
}

int
kfe_solver_eval_beta_plus(kfe_solver_t *self,
			  const double *x,
			  int i,
			  int h_axis, double h_len,
			  double *beta_i)
{
  int status = OK;
  const int n_dim = self->prob->n;
  const double hinv = 1.0/self->rgrid.hs[i];

  double *D = NULL, *K = NULL, *x_temp = NULL, *x_local = NULL;

  CHECK_NULL( D = malloc(n_dim*n_dim*sizeof(*D)), OUT_OF_MEM );
  CHECK_NULL( K = malloc(n_dim*sizeof(*K)), OUT_OF_MEM );
  CHECK_NULL( x_temp = malloc(n_dim*sizeof(*x_temp)), OUT_OF_MEM );
  CHECK_NULL( x_local = malloc(n_dim*sizeof(*x_local)), OUT_OF_MEM );

  memcpy(x_local, x, n_dim*sizeof(*x_local));
  x_local[h_axis] += h_len;
  
  CHECK_( kfe_solver_eval_K_(self, x_local,
			     i, 0.5*self->rgrid.hs[i], K, x_temp) );
  const double Kia = K[i];

  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, 0.5*self->rgrid.hs[i],
			     0, 0.0, D, x_temp) );
  const double Dii = D[i*n_dim + i];

  double delta_i;
  CHECK_( kfe_solver_eval_delta_std(self, i, x_local, 0, 0.0, &delta_i) );

  // double S = -Kia*0.5*hinv + Dii*hinv*hinv;
  double S = -Kia*(1.0 - delta_i)*hinv + Dii*hinv*hinv;
  double scale = fabs(Kia*(1.0 - delta_i)*hinv) + fabs(Dii*hinv*hinv);
  double Od = 0.0;
  
  int j;
  for (j = 0; j < n_dim; ++j)
    {
      if (j == i) continue;
      
      const int k = i*n_dim + j;

      double Dij_pi2, Dij_pj2, Dij_mj2;

      CHECK_( kfe_solver_eval_D_(self, x_local,
				 i, +0.5*self->rgrid.hs[i],
				 0, 0.0, D, x_temp) );
      Dij_pi2 = D[k];

      CHECK_( kfe_solver_eval_D_(self, x_local,
				 j, +0.5*self->rgrid.hs[j],
				 0, 0.0, D, x_temp) );
      Dij_pj2 = D[k];

      CHECK_( kfe_solver_eval_D_(self, x_local,
				 j, -0.5*self->rgrid.hs[j],
				 0, 0.0, D, x_temp) );
      Dij_mj2 = D[k];

      double gamma[4];
      double gamma_mj[4];
      CHECK_( kfe_solver_eval_gamma_std(self, i, j, x_local,
					0, 0.0, 0, 0.0,
					gamma) );
      CHECK_( kfe_solver_eval_gamma_std(self, i, j, x_local,
					j, -self->rgrid.hs[j], 0, 0.0,
					gamma_mj) );

      const double hijinv = hinv / self->rgrid.hs[j];
      Od += hijinv*((Dij_pi2 + Dij_pj2)*gamma[1]
		    - (Dij_pi2 + Dij_mj2)*gamma_mj[3]);
    }

  S += Od;
  scale += fabs(Od);

  // if (S < 0.0)
  if (S < 0.0 && !(fabs(S) <= TOLFCT*scale*sqrt(DBL_EPSILON)))
    {
      ERRMSG("FAILED: ENCOUNTERED NEGATIVE TRANSITION RATE\n");
      ERRMSG("beta_i+ = %le < 0.0, i = %d\n", S, i);
      ERRMSG("scaled: %le\n", S/scale, i);
      ERRMSG("delta_i = %lf\n", delta_i);
      int l;
      for (l = 0; l < n_dim; ++l)
	ERRMSG("x[%d] = %lf\n", l, x_local[l]);
      
      if (!ALLOW_NEGTRANSR) FAILWITH(FAILED);
    }
  *beta_i = S;

 exit:
  Xfree(D);
  Xfree(K);
  Xfree(x_temp);
  Xfree(x_local);
  
  return status;
}

int
kfe_solver_eval_beta_minus(kfe_solver_t *self,
			   const double *x,
			   int i,
			   int h_axis, double h_len,
			   double *beta_i)
{
  int status = OK;
  const int n_dim = self->prob->n;
  const double hinv = 1.0/self->rgrid.hs[i];

  double *D = NULL, *K = NULL, *x_temp = NULL, *x_local = NULL;

  CHECK_NULL( D = malloc(n_dim*n_dim*sizeof(*D)), OUT_OF_MEM );
  CHECK_NULL( K = malloc(n_dim*sizeof(*K)), OUT_OF_MEM );
  CHECK_NULL( x_temp = malloc(n_dim*sizeof(*x_temp)), OUT_OF_MEM );
  CHECK_NULL( x_local = malloc(n_dim*sizeof(*x_local)), OUT_OF_MEM );

  memcpy(x_local, x, n_dim*sizeof(*x_local));
  x_local[h_axis] += h_len;
  
  CHECK_( kfe_solver_eval_K_(self, x_local, i, -0.5*self->rgrid.hs[i], K, x_temp) );
  const double Kia = K[i];

  CHECK_( kfe_solver_eval_D_(self, x_local, i, -0.5*self->rgrid.hs[i], 0, 0.0, D, x_temp) );
  const double Dii = D[i*n_dim + i];

  double delta_i;
  CHECK_( kfe_solver_eval_delta_std(self, i, x_local, i, -self->rgrid.hs[i], &delta_i) );
  CHECK( !isfinite(delta_i), FAILED );

  // double S = Kia*0.5*hinv + Dii*hinv*hinv;
  double S = Kia*delta_i*hinv + Dii*hinv*hinv;
  double scale = fabs(Kia*delta_i*hinv) + fabs(Dii*hinv*hinv);
  double Od = 0.0;

  int j;
  for (j = 0; j < n_dim; ++j)
    {
      const int k = i*n_dim + j;
      
      if (j == i) continue;

      double Dij_mi2, Dij_pj2, Dij_mj2;

      CHECK_( kfe_solver_eval_D_(self, x_local,
				 i, -0.5*self->rgrid.hs[i],
				 0, 0.0, D, x_temp) );
      CHECK( !isfinite(D[k]), FAILED );
      Dij_mi2 = D[k];

      CHECK_( kfe_solver_eval_D_(self, x_local,
				 j, +0.5*self->rgrid.hs[j],
				 0, 0.0, D, x_temp) );
      CHECK( !isfinite(D[k]), FAILED );
      Dij_pj2 = D[k];

      CHECK_( kfe_solver_eval_D_(self, x_local,
				 j, -0.5*self->rgrid.hs[j],
				 0, 0.0, D, x_temp) );
      CHECK( !isfinite(D[k]), FAILED );
      Dij_mj2 = D[k];

      double gamma_mi[4];
      double gamma_mi_mj[4]; 

      CHECK_( kfe_solver_eval_gamma_std(self, i, j, x_local,
					i, -self->rgrid.hs[i],
					0, 0.0,
					gamma_mi) );
      CHECK_( kfe_solver_eval_gamma_std(self, i, j, x_local,
					i, -self->rgrid.hs[i],
					j, -self->rgrid.hs[j],
					gamma_mi_mj) );

      CHECK( !isfinite(gamma_mi[0]), FAILED );
      CHECK( !isfinite(gamma_mi_mj[2]), FAILED );
      
      const double hijinv = hinv / self->rgrid.hs[j];
      Od -= hijinv*((Dij_mi2 + Dij_pj2)*gamma_mi[0]
		    - (Dij_mi2 + Dij_mj2)*gamma_mi_mj[2]);
    }

  S += Od;
  scale += fabs(Od);
  
  // if (S < 0.0)
  if (S < 0.0 && !(fabs(S) <= TOLFCT*scale*sqrt(DBL_EPSILON)))
    {
      ERRMSG("beta_i- = %lf < 0.0, i = %d\n", S, i);
      ERRMSG("delta_i = %lf\n", delta_i);
      int l;
      for (l = 0; l < n_dim; ++l)
	ERRMSG("x[%d] = %lf\n", l, x_local[l]);

      if (!ALLOW_NEGTRANSR) FAILWITH(FAILED);
    }
  *beta_i = S;

 exit:
  Xfree(D);
  Xfree(K);
  Xfree(x_temp);
  Xfree(x_local);
  
  return status;
}


int
kfe_solver_eval_theta_S_plus(kfe_solver_t *self,
			     const double *x,
			     int i, int j,
			     int h_axis, double h_len,
			     int h_axis_2, double h_len_2,
			     double *theta_S)
{
  int status = OK;
  const int n_dim = self->prob->n;
  const double hijinv = 1.0 / (self->rgrid.hs[i]*self->rgrid.hs[j]);

  double *D = NULL, *x_local = NULL, *x_temp = NULL;

  CHECK_NULL( D = malloc(n_dim*n_dim*sizeof(*D)), OUT_OF_MEM );
  CHECK_NULL( x_local = malloc(n_dim*sizeof(*x_local)), OUT_OF_MEM );
  CHECK_NULL( x_temp = malloc(n_dim*sizeof(*x_temp)), OUT_OF_MEM );

  memcpy(x_local, x, n_dim*sizeof(*x_local));
  x_local[h_axis] += h_len;
  x_local[h_axis_2] += h_len_2;

  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, 0.5*self->rgrid.hs[i],
			     0, 0.0, D, x_temp) );
  const double Da = D[i*n_dim + j];

  CHECK_( kfe_solver_eval_D_(self, x_local,
			     j, 0.5*self->rgrid.hs[j],
			     0, 0.0, D, x_temp) );
  const double Db = D[i*n_dim + j];

  double gamma[4];
  CHECK_( kfe_solver_eval_gamma_std(self, i, j, x_local,
				    0, 0.0, 0, 0.0,
				    gamma) );
  
  *theta_S = hijinv*(Da + Db)*gamma[3];
  CHECK( *theta_S < 0.0, FAILED );

 exit:
  Xfree(D);
  Xfree(x_local);
  Xfree(x_temp);

  return status;
}

int
kfe_solver_eval_theta_S_minus(kfe_solver_t *self,
			      const double *x,
			      int i, int j,
			      int h_axis, double h_len,
			      int h_axis_2, double h_len_2,
			      double *theta_S)
{
  int status = OK;
  const int n_dim = self->prob->n;
  const double hijinv = 1.0 / (self->rgrid.hs[i]*self->rgrid.hs[j]);

  double *D = NULL, *x_local = NULL, *x_temp = NULL;

  CHECK_NULL( D = malloc(n_dim*n_dim*sizeof(*D)), OUT_OF_MEM );
  CHECK_NULL( x_local = malloc(n_dim*sizeof(*x_local)), OUT_OF_MEM );
  CHECK_NULL( x_temp = malloc(n_dim*sizeof(*x_temp)), OUT_OF_MEM );

  memcpy(x_local, x, n_dim*sizeof(*x_local));
  x_local[h_axis] += h_len;
  x_local[h_axis_2] += h_len_2;

  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, -0.5*self->rgrid.hs[i],
			     0, 0.0, D, x_temp) );
  const double Da = D[i*n_dim + j];

  CHECK_( kfe_solver_eval_D_(self, x_local,
			     j, -0.5*self->rgrid.hs[j],
			     0, 0.0, D, x_temp) );
  const double Db = D[i*n_dim + j];

  double gamma[4];
  CHECK_( kfe_solver_eval_gamma_std(self, i, j, x_local,
				    i, -self->rgrid.hs[i],
				    j, -self->rgrid.hs[j],
				    gamma) );

  *theta_S = hijinv*(Da + Db)*gamma[0];
  CHECK( *theta_S < 0.0, FAILED );

 exit:
  Xfree(D);
  Xfree(x_local);
  Xfree(x_temp);

  return status;
}

int
kfe_solver_eval_theta_A(kfe_solver_t *self,
			const double *x,
			int i, int j,
			int h_axis, double h_len,
			int h_axis_2, double h_len_2,
			double *theta_A)
{
  int status = OK;
  const int n_dim = self->prob->n;
  const double hijinv = 1.0 / (self->rgrid.hs[i]*self->rgrid.hs[j]);

  double *D = NULL, *x_local = NULL, *x_temp = NULL;

  CHECK_NULL( D = malloc(n_dim*n_dim*sizeof(*D)), OUT_OF_MEM );
  CHECK_NULL( x_local = malloc(n_dim*sizeof(*x_local)), OUT_OF_MEM );
  CHECK_NULL( x_temp = malloc(n_dim*sizeof(*x_temp)), OUT_OF_MEM );

  memcpy(x_local, x, n_dim*sizeof(*x_local));
  x_local[h_axis] += h_len;
  x_local[h_axis_2] += h_len_2;

  CHECK_( kfe_solver_eval_D_(self, x_local,
			     i, 0.5*self->rgrid.hs[i],
			     0, 0.0, D, x_temp) );
  const double Da = D[i*n_dim + j];

  CHECK_( kfe_solver_eval_D_(self, x_local,
			     j, -0.5*self->rgrid.hs[j],
			     0, 0.0, D, x_temp) );
  const double Db = D[i*n_dim + j];

  double gamma[4];
  CHECK_( kfe_solver_eval_gamma_std(self, i, j, x_local,
				    0, 0.0,
				    j, -self->rgrid.hs[j],
				    gamma) );
  *theta_A = -hijinv*(Da + Db)*gamma[1];
  CHECK( *theta_A < 0.0, FAILED );

 exit:
  Xfree(D);
  Xfree(x_local);
  Xfree(x_temp);

  return status;
}

int kfe_solver_iwork_size(int n_dim) { return n_dim; }
int kfe_solver_work_size(int n_dim) { return 2*n_dim + n_dim*n_dim; }

int
kfe_solver_init(kfe_solver_t *self,
		kfe_problem_t *prob,
		const int *ns)
{
  int status = OK;
  *self = kfe_solver_nil;

  self->prob = prob;

  CHECK_( rgrid_init(&self->rgrid, prob->n, ns, prob->bounds) );

  const int n_grid = self->rgrid.m;
  const int n_gamma = kfe_gamma_size(prob->n);
  const int n_dim = prob->n;
  const int gamma_stride = kfe_gamma_stride(n_dim);
  
  CHECK_NULL( self->gamma = malloc(n_grid*n_gamma*sizeof(*self->gamma)),
	      OUT_OF_MEM );
  CHECK_NULL( self->delta = malloc(n_grid*sizeof(*self->delta)),
	      OUT_OF_MEM );
  CHECK_NULL( self->iwork = malloc(kfe_solver_iwork_size(n_dim)*sizeof(*self->iwork)),
	      OUT_OF_MEM );
  CHECK_NULL( self->work = malloc(kfe_solver_work_size(n_dim)*sizeof(*self->work)),
	      OUT_OF_MEM );

  int *grid_i = self->iwork;
  double *x = self->work;
  double *x_temp = self->work + n_dim;
  double *D = self->work + 2*n_dim;
  double *K = self->work + 2*n_dim + n_dim*n_dim;

  /* /\* GAMMA LOOP *\/ */
  /* rgrid_start(&self->rgrid, grid_i, x); */
  /* do */
  /*   { */
  /*     int i, j; */

  /*     /\* Second derivative axes *\/ */
  /*     for (i = 0; i < n_dim; ++i) */
  /* 	{ */
  /* 	  for (j = i + 1; j < n_dim; ++j) */
  /* 	    { */
  /* 	      double gamma[4]; */

  /* 	      CHECK_( kfe_solver_eval_gamma_std(self, i, j, x, 0, 0.0, gamma, */
  /* 						D, x_temp, K) ); */
	      
  /* 	      int offset = n_gamma*rgrid_rindex(&self->rgrid, grid_i); */
  /* 	      self->gamma[offset + kfe_gamma_index(n_dim, 0, i, j)] = gamma[0]; */
  /* 	      self->gamma[offset + kfe_gamma_index(n_dim, 1, i, j)] = gamma[1]; */
  /* 	      self->gamma[offset + kfe_gamma_index(n_dim, 2, i, j)] = gamma[2]; */
  /* 	      self->gamma[offset + kfe_gamma_index(n_dim, 3, i, j)] = gamma[3]; */
  /* 	    } */
  /* 	} */
      
  /*   } */
  /* while (rgrid_step(&self->rgrid, grid_i, x)); */

  /* /\* DELTA LOOP *\/ */
  /* rgrid_start(&self->rgrid, grid_i, x); */
  /* do */
  /*   { */
  /*     int i, j; */

  /*     for (i = 0; i < n_dim; ++i) */
  /* 	{ */
  /* 	  xxxxx; */
  /* 	} */
      
  /*   } */
  /* while (rgrid_step(&self->rgrid, grid_i, x)); */

  
 exit:
  if (status)
    {
      kfe_solver_destroy(self);
      *self = kfe_solver_nil;
    }

  return status;
}

int
kfe_solver_eval_Gamma(kfe_solver_t *self,
		      const double *x,
		      double *Gamma)
{
  int status = OK;
  const int n_dim = self->prob->n;
  int i, j;
  double S = 0;

  for (i = 0; i < n_dim; ++i)
    {
      const double hi = self->rgrid.hs[i];
      double beta_plus = 0.0, beta_minus = 0.0;

      if (kfe_solver_at_bound(self, x, i) != -1)
	CHECK_( kfe_solver_eval_beta_plus(self, x, i, i, -hi, &beta_plus) );

      if (kfe_solver_at_bound(self, x, i) != +1)
	CHECK_( kfe_solver_eval_beta_minus(self, x, i, i, +hi, &beta_minus) );

      CHECK( !isfinite(beta_plus), FAILED );
      CHECK( !isfinite(beta_minus), FAILED );

      S += beta_plus + beta_minus;
      
      for (j = i + 1; j < n_dim; ++j)
	{
	  const double hj = self->rgrid.hs[j];	  
	  double
	    theta_S_plus = 0.0, theta_S_minus = 0.0,
	    theta_A = 0.0, theta_A_prime = 0.0;

	  if (kfe_solver_at_bound(self, x, i) != -1 &&
	      kfe_solver_at_bound(self, x, j) != -1)
	    CHECK_( kfe_solver_eval_theta_S_plus(self, x, i, j,
						 i, -hi,
						 j, -hj,
						 &theta_S_plus) );
	  
	  if (kfe_solver_at_bound(self, x, i) != +1 &&
	      kfe_solver_at_bound(self, x, j) != +1)
	    CHECK_( kfe_solver_eval_theta_S_minus(self, x, i, j,
						  i, +hi,
						  j, +hj,
						  &theta_S_minus) );
	  
	  if (kfe_solver_at_bound(self, x, i) != -1 &&
	      kfe_solver_at_bound(self, x, j) != +1)
	    CHECK_( kfe_solver_eval_theta_A(self, x, i, j,
					    i, -hi,
					    j, +hj,
					    &theta_A) );

	  if (kfe_solver_at_bound(self, x, i) != +1 &&
	      kfe_solver_at_bound(self, x, j) != -1)
	    CHECK_( kfe_solver_eval_theta_A(self, x, j, i,
					    i, +hi,
					    j, -hj,
					    &theta_A_prime) );

	  CHECK( !isfinite(theta_S_plus), FAILED );
	  CHECK( !isfinite(theta_S_minus), FAILED );
	  CHECK( !isfinite(theta_A), FAILED );
	  CHECK( !isfinite(theta_A_prime), FAILED );

	  S += theta_S_plus + theta_S_minus + theta_A + theta_A_prime;
	}
    }

  *Gamma = S;

 exit:

  return status;
}

void
shift_(int n, const int *i, int *i_prime,
       int axis, int dir)
{
  memcpy(i_prime, i, n*sizeof(*i_prime));
  i_prime[axis] += dir;
}

void
shift2_(int n, const int *i, int *i_prime,
	int axis, int dir,
	int axis_2, int dir_2)
{
  memcpy(i_prime, i, n*sizeof(*i_prime));
  i_prime[axis] += dir;
  i_prime[axis_2] += dir_2;
}

cs *
kfe_solver_make_KFE_matrix(kfe_solver_t *self)
{
  int status = OK;
  const int n_dim = self->prob->n;
  const int n_tot = self->rgrid.m;

  double *x = NULL;
  int *i_grid = NULL, *j_grid = NULL;
  cs *A = NULL;

  CHECK_NULL( A = cs_spalloc(n_tot, n_tot, n_dim*n_tot,
			     SPALLOC_VALUES, SPALLOC_TRIPLET),
	      OUT_OF_MEM );
  CHECK_NULL( i_grid = malloc(n_dim*sizeof(*i_grid)),
	      OUT_OF_MEM );
  CHECK_NULL( j_grid = malloc(n_dim*sizeof(*j_grid)),
	      OUT_OF_MEM );
  CHECK_NULL( x = malloc(n_dim*sizeof(*x)),
	      OUT_OF_MEM );
  
  rgrid_start(&self->rgrid, i_grid, x);
  do
    {
      int row = rgrid_rindex(&self->rgrid, i_grid);
      int col;

      double Gamma;
      CHECK_( kfe_solver_eval_Gamma(self, x, &Gamma) ); 
      CHECK( !isfinite(Gamma), FAILED );
      CHECK( !cs_entry(A, row, row, -Gamma), OUT_OF_MEM );

      int i;
      for (i = 0; i < n_dim; ++i)
	{
	  int j;

	  if (i_grid[i] < self->rgrid.ns[i] - 1)
	    {
	      double beta_plus;
	      CHECK_( kfe_solver_eval_beta_plus(self, x, i, 0, 0.0, &beta_plus) );

	      shift_(n_dim, i_grid, j_grid, i, +1);
	      col = rgrid_rindex(&self->rgrid, j_grid);
	      CHECK( !(col < n_tot && col >= 0), FAILED );
	      CHECK( !cs_entry(A, row, col, beta_plus), OUT_OF_MEM );
	    }
	  else
	    {
	      // CHECK( beta_plus == 0.0, FAILED );
	    }

	  if (i_grid[i] > 0)
	    {
	      double beta_minus;
	      CHECK_( kfe_solver_eval_beta_minus(self, x, i, 0, 0.0, &beta_minus) );

	      shift_(n_dim, i_grid, j_grid, i, -1);
	      col = rgrid_rindex(&self->rgrid, j_grid);
	      CHECK( !(col < n_tot && col >= 0), FAILED );
	      CHECK( !cs_entry(A, row, col, beta_minus), OUT_OF_MEM );
	    }
	  else 
	    {
	      // CHECK( beta_minus == 0.0, FAILED );
	    }

	  
	  for (j = i + 1; j < n_dim; ++j)
	    {
	      /* Sanity check for gamma coefficients. */
	      double gamma_ij[4];
	      double gamma_ji[4];

	      CHECK_( kfe_solver_eval_gamma_std(self, i, j, x,
						0, 0.0, 0, 0.0,
						gamma_ij) );
	      CHECK_( kfe_solver_eval_gamma_std(self, j, i, x,
						0, 0.0, 0, 0.0,
						gamma_ji) );
	      CHECK( gamma_ij[0] != gamma_ji[0], FAILED );
	      CHECK( gamma_ij[3] != gamma_ji[3], FAILED );
	      CHECK( gamma_ij[1] != gamma_ji[2], FAILED );
	      CHECK( gamma_ij[2] != gamma_ji[1], FAILED );

	      double theta_S_plus, theta_S_minus, theta_A, theta_A_prime;

	      CHECK_( kfe_solver_eval_theta_S_plus(self, x, i, j, 0, 0.0, 0, 0.0, &theta_S_plus) );
	      CHECK_( kfe_solver_eval_theta_S_minus(self, x, i, j, 0, 0.0, 0, 0.0, &theta_S_minus) );
	      CHECK_( kfe_solver_eval_theta_A(self, x, i, j, 0, 0.0, 0, 0.0, &theta_A) );
	      CHECK_( kfe_solver_eval_theta_A(self, x, j, i, 0, 0.0, 0, 0.0, &theta_A_prime) );

	      if (i_grid[i] < self->rgrid.ns[i] - 1 && i_grid[j] < self->rgrid.ns[j] - 1)
		{
		  shift2_(n_dim, i_grid, j_grid, i, +1, j, +1);
		  col = rgrid_rindex(&self->rgrid, j_grid);
		  CHECK( !(col < n_tot && col >= 0), FAILED );		  
		  CHECK( !cs_entry(A, row, col, theta_S_plus), OUT_OF_MEM );
		}
	      else 
		{
		  // CHECK( theta_S_plus == 0.0, FAILED );
		}

	      if (i_grid[i] > 0 && i_grid[j] > 0)
		{
		  shift2_(n_dim, i_grid, j_grid, i, -1, j, -1);
		  col = rgrid_rindex(&self->rgrid, j_grid);
		  CHECK( !(col < n_tot && col >= 0), FAILED );		  
		  CHECK( !cs_entry(A, row, col, theta_S_minus), OUT_OF_MEM );
		}
	      else 
		{
		  // CHECK( theta_S_minus == 0.0, FAILED );
		}

	      if (i_grid[i] < self->rgrid.ns[i] - 1 && i_grid[j] > 0)
		{
		  shift2_(n_dim, i_grid, j_grid, i, +1, j, -1);
		  col = rgrid_rindex(&self->rgrid, j_grid);
		  CHECK( !(col < n_tot && col >= 0), FAILED );		  
		  CHECK( !cs_entry(A, row, col, theta_A), OUT_OF_MEM );
		}
	      else 
		{
		  // CHECK( theta_A == 0.0, FAILED );
		}

	      if (i_grid[i] > 0 && i_grid[j] < self->rgrid.ns[j] - 1)
		{
		  shift2_(n_dim, i_grid, j_grid, i, -1, j, +1);
		  col = rgrid_rindex(&self->rgrid, j_grid);
		  CHECK( !(col < n_tot && col >= 0), FAILED );		  
		  CHECK( !cs_entry(A, row, col, theta_A_prime), OUT_OF_MEM );
		}
	      else 
		{
		  // CHECK( theta_A_prime == 0.0, FAILED );
		}
	    }
	}
    }
  while (rgrid_step(&self->rgrid, i_grid, x));

  CHECK_( compress_(&A) );
  
 exit:
  if (status)
    {
      Xfree_(A, cs_spfree); A = NULL;
    }
  Xfree(x);
  Xfree(i_grid);
  Xfree(j_grid);

  return A;
}

int
kfe_solver_set_init_state(kfe_solver_t *self,
			  const double *x_init,
			  double *data)
{
  int status = OK;
  const int n_tot = self->rgrid.m;
  const int n_dim = self->prob->n;
  int *grid_i = NULL;

  CHECK_NULL( grid_i = malloc(n_dim*sizeof(*grid_i)), OUT_OF_MEM );

  int i;
  for (i = 0; i < n_tot; ++i)
    {
      data[i] = 0.0;
    }

  rgrid_nearest_i(&self->rgrid, x_init, grid_i);
  i = rgrid_rindex(&self->rgrid, grid_i);
  data[i] = 1.0;

 exit:
  Xfree(grid_i);

  return status;
}

typedef struct {
  kfe_sde_t sde;
  void *sde_params;
  double *s;
  double *dsdx;
  double *b;
} wrapper_t;

/*
 * Ki = bi - (dj Dij)
 */

int
eval_sde_D(const double *x, double *D, void *params)
{
  int status = OK;
  wrapper_t *p = params;
  const int n = p->sde.n;
  const int k = p->sde.k;

  CHECK_( p->sde.eval(x, p->b, p->s, p->dsdx, p->sde_params) );

  int i, j;
  for (i = 0; i < n; ++i)
    {
      for (j = i; j < n; ++j)
	{
	  int l;
	  double s = 0.0;

	  for (l = 0; l < k; ++l)
	    {
	      s += p->s[i*k + l] * p->s[j*k + l];
	    }

	  D[i*n + j] = D[j*n + i] = 0.5 * s;

	  CHECK( !isfinite(s), FAILED );
	}
    }

 exit:
  return status;
}

int
eval_sde_K(const double *x, double *K, void *params)
{
  int status = OK;
  wrapper_t *p = params;
  const int n = p->sde.n;
  const int k = p->sde.k;

  /* Ki = bi - dj Dij  
   * dj Dij = 1/2 dj (sil sjl)
   *        = 1/2 ((dj sil) sjl + sil (dj sjl))  
   */

  CHECK_( p->sde.eval(x, p->b, p->s, p->dsdx, p->sde_params) );

  int i, j, l;
  for (i = 0; i < n; ++i)
    {
      double s = 0.0;      

      for (j = 0; j < n; ++j)
	{
	  double *dj_s = p->dsdx + j*n*k;

	  for (l = 0; l < k; ++l)
	    {
	      s += dj_s[i*k + l] * p->s[j*k + l];
	      s += p->s[i*k + l] * dj_s[j*k + l];
	    }
	}

      K[i] = p->b[i] - 0.5*s;

      CHECK( !isfinite(K[i]), FAILED );
    }

 exit:
  return status;
}

kfe_problem_t *
kfe_problem_from_sde(kfe_sde_t *sde,
		     void *params)
{
  int status = OK;
  kfe_problem_t *self = NULL;
  wrapper_t *wrapper;
  const int n = sde->n;
  const int k = sde->k;
  void *ptr = NULL;
  char *block;
  const int block_size
    = sizeof(*self) + sizeof(*wrapper)
    + 2*n*sizeof(*self->bounds)
    + n*sizeof(*wrapper->b)
    + n*k*sizeof(*wrapper->s)
    + n*n*k*sizeof(*wrapper->dsdx);

  CHECK_NULL( ptr = malloc(block_size), OUT_OF_MEM );
  block = ptr;
  self          = (void*) block; block += sizeof(*self);
  wrapper       = (void*) block; block += sizeof(*wrapper);
  self->bounds  = (void*) block; block += 2*n*sizeof(*self->bounds);
  wrapper->b    = (void*) block; block += n*sizeof(*wrapper->b);
  wrapper->s    = (void*) block; block += n*k*sizeof(*wrapper->s);
  wrapper->dsdx = (void*) block; 

  memcpy(&wrapper->sde, sde, sizeof(wrapper->sde));
  wrapper->sde_params = params;

  CHECK_( sde->eval_bounds(self->bounds, params) );
  int i;
  for (i = 0; i < 2*n; ++i) CHECK( !isfinite(self->bounds[i]), FAILED );

  self->n = sde->n;
  self->D = eval_sde_D;
  self->K = eval_sde_K;
  self->params = wrapper;
  
 exit:
  if (status) { Xfree(ptr); self = NULL; }
  
  return self;
}

void
kfe_problem_free(kfe_problem_t *self)
{
  Xfree(self);
}

