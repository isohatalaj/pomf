
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rgrid.h"


const rgrid_t rgrid_nil = RGRID_NIL;

void
rgrid_destroy(rgrid_t *self)
{
  if (self == NULL) return;
  if (self->ns) free(self->ns);
  if (self->bs) free(self->bs);
  if (self->hs) free(self->hs);
  if (self->strides) free(self->strides);
}

int
rgrid_rindex(const rgrid_t *self, int *i)
{
  int j = self->strides[0]*i[0];
  int k;

  for (k = 1; k < self->dim; ++k)
    j += self->strides[k]*i[k];

  return j;
}

int
rgrid_init(rgrid_t *self, int dim, const int *ns, const double *bs)
{
  int i;  
  int status = OK;

  *self = rgrid_nil;

  CHECK_NULL( self->ns = malloc(dim*sizeof(*self->ns)),           OUT_OF_MEM );
  CHECK_NULL( self->bs = malloc(2*dim*sizeof(*self->bs)),         OUT_OF_MEM );
  CHECK_NULL( self->hs = malloc(dim*sizeof(*self->hs)),           OUT_OF_MEM );
  CHECK_NULL( self->strides = malloc(dim*sizeof(*self->strides)), OUT_OF_MEM );

  self->dim = dim;
  memcpy(self->ns, ns, dim*sizeof(*self->ns));
  memcpy(self->bs, bs, 2*dim*sizeof(*self->bs));
  self->m = iprod(dim, ns);

  for (i = 0; i < dim; ++i)
    {
      if (ns[i] <= 1) { status = INVALID; goto exit; }

      /* if this needs to be modified, be advised that spgrid code
	 assumes the following formula for step length */
      self->hs[i] = (bs[2*i + 1] - bs[2*i + 0]) / (ns[i] - 1);

      CHECK( !isfinite(self->hs[i]), FAILED );
      CHECK( self->hs[i] <= 0.0, FAILED );

      self->strides[dim - i - 1] = i == 0 ? 1 : ns[dim - i]*self->strides[dim - i];
    }

 exit:
  if (status) rgrid_destroy(self);
  return status;
}

void
rgrid_start(const rgrid_t *self, int *i, double *x)
{
  int k;
  for (k = 0; k < self->dim; ++k) i[k] = 0;

  if (x)
    {
      for (k = 0; k < self->dim; ++k) x[k] = self->bs[2*k + 0];
    }
}

int
rgrid_step(const rgrid_t *self, int *i, double *x)
{
  int done = 1;
  int k;
  for (k = self->dim - 1; k >= 0; --k)
    {
      if (i[k] == self->ns[k] - 1)
	{
	  i[k] = 0;
	}
      else
	{
	  i[k]++;
	  done = 0;
	  break;
	}
    }

  if (x)
    {
      for (k = 0; k < self->dim; ++k)
	{
	  x[k] = self->hs[k]*i[k] + self->bs[2*k + 0];
	}
    }
  
  return !done;
}

int
rgrid_print_data(const rgrid_t *grid, double **nodes, int n, double **y,
		 const char **headers)
{
  int i, j, k;

  int *mi = malloc(grid->dim*sizeof(*mi));
  double *x = malloc(grid->dim*sizeof(*x));

  if (mi == NULL || x == NULL) { return OUT_OF_MEM; }

  i = 0;
  rgrid_start(grid, mi, x);
  do
    {
      if (headers && i % 32 == 0)
	{
	  int col = 1;
	  
	  for (k = 0; k < grid->dim; ++k, ++col)
	    printf("#%5d%17s%-3d", col, "x", k);
      
	  for (j = 0; j < n; ++j, ++col)
	    printf("#%5d%20s", col, headers[j]);

	  printf("\n");
	}

      for (k = 0; k < grid->dim; ++k)
	printf(" %25.15le", nodes ? (nodes[k])[i] : x[k]);
      
      for (j = 0; j < n; ++j)
	printf(" %25.15le", (y[j])[i]);

      printf("\n");

      if (grid->dim == 2 && mi[1] == grid->ns[1] - 1)
	printf("\n");

      i++;
    }
  while (rgrid_step(grid, mi, x));

  return 0;
}

int
rgrid_marginals(const rgrid_t *grid, const double *x,
		int n, double **ys, int axis,
		rgrid_t *mgrid, double ***us_out)
{
  int status = OK;
  int i, j, k;
  const int dim = grid->dim;
  const int mdim = dim - 1;

  double *mbs = NULL;
  int *mns = NULL;

  double **us = NULL;

  *mgrid = rgrid_nil;

  CHECK_NULL( mbs = malloc(2*mdim*sizeof(*mbs)), OUT_OF_MEM );
  CHECK_NULL( mns = malloc(mdim*sizeof(*mns)),   OUT_OF_MEM );
  CHECK_NULL( us = calloc(n, sizeof(*us)),       OUT_OF_MEM );

  for (i = 0; i < mdim; ++i)
    {
      j = i < axis ? i : i + 1;
      mbs[2*i + 0] = grid->bs[2*j + 0];
      mbs[2*i + 1] = grid->bs[2*j + 1];
      mns[i] = grid->ns[j];
    }

  CHECK( rgrid_init(mgrid, mdim, mns, mbs), FAILED );

  for (i = 0; i < n; ++i)
    CHECK_NULL( us[i] = malloc(mgrid->m*sizeof(*us[i])),
		OUT_OF_MEM );
  
  int *mi = mns; /* Repurpose this array. */
  k = 0;

  rgrid_start(mgrid, mi, NULL);
  do
    {
      i = 0;
      for (j = 0; j < mdim; ++j)
	{
	  i += grid->strides[j >= axis ? j + 1 : j]*mi[j];
	}

      int l;
      for (l = 0; l < n; ++l)
	{
	  if (x == NULL)
	    {
	      double sum;

	      sum = 0.5 * (ys[l])[i + grid->strides[axis]*0];
	      for (j = 1; j < grid->ns[axis] - 1; ++j)
		{
		  sum += (ys[l])[i + grid->strides[axis]*j];
		}
	      sum += 0.5*(ys[l])[i + grid->strides[axis]*j];

	      (us[l])[k] = sum*grid->hs[axis];
	    }
	  else
	    {
	      double sum = 0;

	      for (j = 1; j < grid->ns[axis]; ++j)
		{
		  const double y0 = (ys[l])[i + grid->strides[axis]*(j - 1)];
		  const double y1 = (ys[l])[i + grid->strides[axis]*j];
		  const double x0 = x[i + grid->strides[axis]*(j - 1)];
		  const double x1 = x[i + grid->strides[axis]*j];
		  const double dx = x1 - x0;
		  
		  sum += 0.5*(y0 + y1)*dx;
		}

	      (us[l])[k] = sum;
	    }
	}
      k++;
    }
  while (rgrid_step(mgrid, mi, NULL));

 exit:
  Xfree(mbs);
  Xfree(mns);

  if (status == OK)
    {
      *us_out = us;
    }
  else
    {
      rgrid_destroy(mgrid);
      Xfree_many(us, n, free);
    }

  return status;
}

int
rgrid_nearest_i(const rgrid_t *self, const double *x, int *i)
{
  int k;
  for (k = 0; k < self->dim; ++k)
    {
      if (x[k] <= self->bs[2*k+0])
	i[k] = 0;
      else if (x[k] >= self->bs[2*k+1])
	i[k] = self->ns[k] - 1;

      const double h_ = (self->bs[2*k+1] - self->bs[2*k+0]) / self->ns[k];
      i[k] = floor((x[k] - self->bs[2*k+0]) / h_);
    }

  return OK;
}

double
rgrid_volel(const rgrid_t *self)
{
  return dprod(self->dim, self->hs);
}
