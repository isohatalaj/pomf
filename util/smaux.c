
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
      double iprev, inext, hinv;

      if (j == 0)
	{
	  iprev = i;
	  inext = i + astride;
	  hinv = hinv0;
	}
      else if (j == an - 1)
	{
	  iprev = i - astride;
	  inext = i;
	  hinv = hinv0;
	}
      else
	{
	  iprev = i - astride;
	  inext = i + astride;
	  hinv = 0.5*hinv0;
	}

      CHECK( !cs_entry(d, i, inext, +hinv), FAILED );
      CHECK( !cs_entry(d, i, iprev, -hinv), FAILED );
      
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
