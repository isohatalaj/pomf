
#include "smaux.h"

int
compress_(cs **m)
{
  cs *m_ = cs_compress(*m);
  if (m_ == NULL) return OUT_OF_MEM;

  cs_spfree(*m);
  *m = m_;

  return OK;
}

cs *
diag(int n, const double *data, int stride, int triplet)
{
  int i, status = OK;
  cs *D = NULL;

  CHECK_NULL( D = cs_spalloc(n, n, n,
			     SPALLOC_VALUES, SPALLOC_TRIPLET),
	      OUT_OF_MEM );

  for (i = 0; i < n; ++i)
    {
      CHECK( !cs_entry(D, i, i, data[i*stride]), OUT_OF_MEM );
    }

  if (!triplet)
    {
      CHECK( compress_(&D), OUT_OF_MEM );
    }

 exit:
  if (status)
    {
      Xfree_(D, cs_spfree);
      D = NULL;
    }

  return D;
}

cs *
eye(int n, int triplet)
{
  int i, status = OK;
  cs *D = NULL;

  CHECK_NULL( D = cs_spalloc(n, n, n,
			 SPALLOC_VALUES, SPALLOC_TRIPLET),
	      OUT_OF_MEM );

  for (i = 0; i < n; ++i)
    {
      CHECK( !cs_entry(D, i, i, 1.0), OUT_OF_MEM );
    }

  if (!triplet)
    {
      CHECK( compress_(&D), OUT_OF_MEM );
    }

 exit:
  if (status)
    {
      Xfree_(D, cs_spfree);
      D = NULL;
    }

  return D;
}

cs *
deriv(const rgrid_t *grid, int axis, int triplet)
{
  int status = OK;
  cs *d = NULL;
  int i, j, k;
  const int m = grid->m;

  const int astride = grid->strides[axis];
  const int an = grid->ns[axis];
  const double hinv0 = 1.0 / grid->hs[axis];

  CHECK_NULL( d = cs_spalloc(m, m, 2*m,
			     SPALLOC_VALUES, SPALLOC_TRIPLET),
	      OUT_OF_MEM );

  /* j is the index on the diff'd axis, k counts up to each increment
   * of j, that is, k = 0...strides[axis]-1 */

  for (i = 0, j = 0, k = 0; i < m; ++i)
    {
      if (j == 0)
	{
	  CHECK( !cs_entry(d, i, i + 0*astride, -(3.0/2.0) * hinv0), FAILED );
	  CHECK( !cs_entry(d, i, i + 1*astride,     2.0    * hinv0), FAILED );
	  CHECK( !cs_entry(d, i, i + 2*astride, -(1.0/2.0) * hinv0), FAILED );
	}
      else if (j == an - 1)
	{
	  CHECK( !cs_entry(d, i, i - 0*astride, +(3.0/2.0) * hinv0), FAILED );
	  CHECK( !cs_entry(d, i, i - 1*astride,    -2.0    * hinv0), FAILED );
	  CHECK( !cs_entry(d, i, i - 2*astride, +(1.0/2.0) * hinv0), FAILED );
	}
      else
	{
	  CHECK( !cs_entry(d, i, i + 1*astride, +0.5 * hinv0), FAILED );
	  CHECK( !cs_entry(d, i, i - 1*astride, -0.5 * hinv0), FAILED );
	}

      if (k == astride - 1)
	{
	  if (j == an - 1) j = 0;
	  else ++j;
	  k = 0;
	}
      else
	{
	  ++k;
	}
    }

  if (!triplet)
    {
      CHECK( compress_(&d), OUT_OF_MEM );
    }
  
 exit:
  if (status)
    {
      Xfree_(d, cs_spfree);
      d = NULL;
    }
    
  return d;
}

static int
ixcmp(const void *p1, const void *p2)
{
  return *(int*)p1 - *(int*)p2;
}

int
rowsort(cs *A)
{
  int status = OK;
  int n = A->n, m = A->m;
  struct { int i; double x; } *ix = NULL;

  /* Inefficient, since we're allocating a temp vector which really
     shouldn't be necessary. I just want this working quick...
   */
  CHECK_NULL( ix = malloc(m*sizeof(*ix)), OUT_OF_MEM );

  int i, j, k;
  for (j = 0; j < n; ++j)
    {
      for (k = 0; k < A->p[j+1] - A->p[j]; ++k)
	{
	  ix[k].i = A->i[k + A->p[j]];
	  ix[k].x = A->x[k + A->p[j]];
	}

      qsort(ix, k, sizeof(*ix), ixcmp);

      for (k = 0; k < A->p[j+1] - A->p[j]; ++k)
	{
	  A->i[k + A->p[j]] = ix[k].i;
	  A->x[k + A->p[j]] = ix[k].x;
	}
    }

 exit:
  Xfree(ix);

  return status;
}
