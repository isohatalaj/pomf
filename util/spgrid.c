
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "spgrid.h"

const spgrid_t spgrid_nil = {0, 0, NULL, NULL, NULL, NULL, 0};


static int
togrid(const double x, const double h)
{
  double ix;
  ix = ceil(x/h - 0.5);
  return ix;
}

static void
coord(const spgrid_t *self, const int *mi, double *x)
{
  int i;
  for (i = 0; i < self->dim; ++i)
    x[i] = mi[i]*self->h[i];
}

double
spgrid_volelem(const spgrid_t *grid, const double *x, double *hx)
{
  double dv = 1.0;
  int i;
  
  if (grid->bounds)
    {
      for (i = 0; i < grid->dim; ++i)
	{
	  double hi;
	  
	  if (grid->bounds[i] == NULL)
	    {
	      hi = grid->h[i];	  
	    }
	  else if (x[i] + 0.5*grid->h[i] > grid->bounds[i][1])
	    {
	      hi = grid->bounds[i][1] - (x[i] - 0.5*grid->h[i]);
	    }
	  else if (x[i] - 0.5*grid->h[i] < grid->bounds[i][0])
	    {
	      hi = x[i] + 0.5*grid->h[i] - grid->bounds[i][0];
	    }
	  else
	    {
	      hi = grid->h[i];
	    }

	  if (hx) hx[i] = hi;
	  
	  dv *= hi;		  
	}
    }
  else
    {
      for (i = 0; i < grid->dim; ++i)
	{
	  dv *= grid->h[i];
	  if (hx) hx[i] = grid->h[i];
	}
    }

  return dv;
}

void
spgrid_datum_free(spgrid_datum_t *self)
{
  if (self == NULL) return;
  Xfree(self->grid_index);
  Xfree(self->payload);
}

spgrid_datum_t *
spgrid_datum_new(int dim, int psize,
		 const int *grid_index, 
		 int dcount,
		 const double *payload)
{
  int status = OK;
  spgrid_datum_t *self = NULL;

  CHECK_NULL( self = malloc(sizeof(*self)), OUT_OF_MEM );
  CHECK_NULL( self->grid_index = malloc(dim*sizeof(*self->grid_index)), OUT_OF_MEM );
  CHECK_NULL( self->payload = malloc(psize*sizeof(*self->payload)), OUT_OF_MEM );

  self->dcount = dcount;
  memcpy(self->grid_index, grid_index, dim*sizeof(*self->grid_index));
  memcpy(self->payload, payload, psize*sizeof(*self->payload));
	      
 exit:
  if (status)
    {
      spgrid_datum_free(self);
      return NULL;
    }

  return self;
}

int
spgrid_init(spgrid_t *grid,
	    int dim, int psize,
	    const double *h)
{
  int status = OK;
  *grid = spgrid_nil;
  
  grid->dim = dim;
  grid->psize = psize;
  CHECK_NULL( grid->h = malloc(dim*sizeof(*grid->h)), OUT_OF_MEM );
  CHECK_NULL( grid->mi = malloc(dim*sizeof(*grid->mi)), OUT_OF_MEM );
  memcpy(grid->h, h, dim*sizeof(*h));
  grid->data = NULL;

 exit:
  if (status)
    {
      spgrid_destroy(grid);
      *grid = spgrid_nil;
    }

  return status;
}

void
spgrid_destroy(spgrid_t *grid)
{
  spgrid_datum_t *p, *tmp;

  HASH_ITER(hh, grid->data, p, tmp)
    {
      HASH_DEL(grid->data, p);
      spgrid_datum_free(p);
    }

  Xfree_many(grid->bounds, grid->dim, free);
  Xfree(grid->h);
  Xfree(grid->mi);

  *grid = spgrid_nil;
}

int
spgrid_set_bound(spgrid_t *grid,
		 int axis,
		 double min,
		 double max)
{
  int status = OK;

  if (grid->bounds == NULL)
    {
      CHECK_NULL( grid->bounds = calloc(grid->dim, sizeof(*grid->bounds)),
		  OUT_OF_MEM );
    }

  if (grid->bounds[axis] == NULL)
    {
      CHECK_NULL( grid->bounds[axis] = malloc(2*sizeof(*grid->bounds[axis])),
		  OUT_OF_MEM );
    }

  grid->bounds[axis][0] = min;
  grid->bounds[axis][1] = max;

 exit:

  return status;
}


int
spgrid_isize(const spgrid_t *grid)
{
  return grid->dim*sizeof(*grid->mi);
}

int
spgrid_add_data_(spgrid_t *grid,
		 const int *index,
		 int dcount,
		 const double *y)
{
  int status = OK;
  int i;
  spgrid_datum_t *new_datum = NULL;

  spgrid_datum_t *p;
  HASH_FIND(hh, grid->data, index, spgrid_isize(grid), p);

  if (p == NULL)
    {
      CHECK_NULL( new_datum = spgrid_datum_new(grid->dim,
					       grid->psize,
					       index, dcount, y),
		  OUT_OF_MEM );
      
      HASH_ADD_KEYPTR(hh, grid->data, index, spgrid_isize(grid), new_datum);

      new_datum = NULL;
    }
  else
    {
      p->dcount += dcount;
      for (i = 0; i < grid->psize; ++i)
	p->payload[i] += y[i];
    }
  
 exit:
  Xfree_(new_datum, spgrid_datum_free);
  
  return status;
}

int
spgrid_add_data(spgrid_t *grid,
		const double *x,
		const double *y)
{
  int i;

  if (grid->bounds)
    {
      for (i = 0; i < grid->dim; ++i)
	{
	  if (grid->bounds[i] == NULL) continue;
	  if (grid->bounds[i][0] > x[i] ||
	      grid->bounds[i][1] < x[i])
	    {
	      grid->dcount_out++;
	      return OK;
	    }
	}
    }

  for (i = 0; i < grid->dim; ++i)
    grid->mi[i] = togrid(x[i], grid->h[i]);

  return spgrid_add_data_(grid, grid->mi, 1, y);
}

/**
 * Given a sparse density grid `grid`, compute the marginal grid
 * `mgrid` by summing over the data along axis number
 * `axis`. Initializes the `mgrid` input to a new `spgrid` object.
 */
int
spgrid_marginal(spgrid_t *grid,
		spgrid_t *mgrid,
		int axis)
{
  int status = OK;
  int i;
  spgrid_datum_t *p, *tmp;
  double *mh = NULL;

  const int mdim = grid->dim - 1;
  const int psize = grid->psize; 

  CHECK_NULL( mh = malloc(mdim*sizeof(*mh)), OUT_OF_MEM );
  for (i = 0; i < mdim; ++i) mh[i] = grid->h[i < axis ? i : i + 1];

  CHECK( spgrid_init(mgrid, mdim, psize, mh), FAILED );

  if (grid->bounds)
    {
      for (i = 0; i < grid->dim; ++i)
	{
	  if (grid->bounds[i] == NULL) continue;
	  if (i == axis) continue;
	  CHECK( spgrid_set_bound(mgrid, i > axis ? i - 1 : i,
				  grid->bounds[i][0],
				  grid->bounds[i][1]),
		 FAILED );
	}
    }

  HASH_ITER(hh, grid->data, p, tmp)
    {
      for (i = 0; i < mdim; ++i)
	mgrid->mi[i] = p->grid_index[i < axis ? i : i + 1];

      CHECK( spgrid_add_data_(mgrid, mgrid->mi, p->dcount, p->payload),
	     FAILED );
    }  

 exit:
  Xfree(mh);

  return status;
}

int
spgrid_to_rgrid(const spgrid_t *grid,
		rgrid_t *rgrid,
		double ***data_out)
{
  int status = OK;
  int i, j;
  spgrid_datum_t *p, *tmp;
  int *minis = NULL, *maxis = NULL, *ns = NULL;
  double *bs = NULL, *x = NULL;
  double **data = NULL;
  int sum = 0;

  *rgrid = rgrid_nil;

  CHECK_NULL( minis = malloc(grid->dim*sizeof(*minis)), OUT_OF_MEM );
  CHECK_NULL( maxis = malloc(grid->dim*sizeof(*maxis)), OUT_OF_MEM );
  CHECK_NULL( ns = malloc(grid->dim*sizeof(*ns)),       OUT_OF_MEM );
  CHECK_NULL( bs = malloc(2*grid->dim*sizeof(*bs)),     OUT_OF_MEM );
  CHECK_NULL( x = malloc(grid->dim*sizeof(*x)),         OUT_OF_MEM );
    
  int first = 1;
  HASH_ITER(hh, grid->data, p, tmp)
    {
      if (first)
	{
	  for (i = 0; i < grid->dim; ++i)
	    maxis[i] = minis[i] = p->grid_index[i];
	  first = 0;
	}
      else
	{
	  for (i = 0; i < grid->dim; ++i)
	    {
	      if (minis[i] > p->grid_index[i]) minis[i] = p->grid_index[i];
	      if (maxis[i] < p->grid_index[i]) maxis[i] = p->grid_index[i];
	    }
	}

      sum += p->dcount;
    }

  for (i = 0; i < grid->dim; ++i)
    {
      bs[2*i + 0] = minis[i] * grid->h[i];
      bs[2*i + 1] = maxis[i] * grid->h[i];
      ns[i] = maxis[i] - minis[i] + 1;
    }

  CHECK( rgrid_init(rgrid, grid->dim, ns, bs), FAILED );

  CHECK_NULL( data = calloc(grid->psize + 1, sizeof(*data)), OUT_OF_MEM );
  for (i = 0; i < grid->psize + 1; ++i)
    CHECK_NULL( data[i] = malloc(rgrid->m*sizeof(*data[i])), OUT_OF_MEM );

  for (i = 0; i < rgrid->m; ++i)
    for (j = 0; j < grid->psize + 1; ++j)
      (data[j])[i] = 0.0;


  HASH_ITER(hh, grid->data, p, tmp)
    {
      coord(grid, p->grid_index, x);      
      double norm = sum * spgrid_volelem(grid, x, NULL);
      
      for (j = 0; j < grid->dim; ++j)
	grid->mi[j] = p->grid_index[j] - minis[j];

      i = rgrid_rindex(rgrid, grid->mi);
      
      data[0][i] = p->dcount / norm;
      
      for (j = 1; j < grid->psize + 1; ++j)
	data[j][i] = p->payload[j - 1] / sum;

      /* data[j][i] / data[0][i] = payload[j']*dv/p->dcount = */
    }  

 exit:
  
  if (status == OK)
    {
      *data_out = data;
    }
  else
    {
      rgrid_destroy(rgrid);
      Xfree_many(data, grid->psize + 1, free);
    }
  
  Xfree(x);
  Xfree(bs);
  Xfree(ns);
  Xfree(minis);
  Xfree(maxis);

  return status;
}

int
spgrid_dump(const spgrid_t *grid, FILE *file)
{
  int i;
  spgrid_datum_t *p, *tmp;
  
  HASH_ITER(hh, grid->data, p, tmp)
    {
      for (i = 0; i < grid->dim; ++i)
	fprintf(file, " %8d", p->grid_index[i]);

      fprintf(file, " %25.15le", (double) p->dcount);

      for (i = 0; i < grid->dim; ++i)
	fprintf(file, " %25.15le", p->grid_index[i]*grid->h[i]);

      for (i = 0; i < grid->psize; ++i)
	fprintf(file, " %25.15le", p->payload[i]);

      fprintf(file, "\n");
    }  

  return OK;
}
