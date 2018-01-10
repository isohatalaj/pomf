
#ifndef SPGRID_H
#define SPGRID_H

#include <stdio.h>

#include <uthash.h>

#include "util.h"
#include "rgrid.h"


typedef struct {
  int *grid_index;
  int dcount;
  double *payload;

  UT_hash_handle hh;
} spgrid_datum_t;

typedef struct {
  int dim;               /*< Number of grid dimensions */
  int psize;             /*< Length of payload vector */
  double *h;             /*< Grid step length, size `dim` */
  int *mi;               /*< Work array for use as hash key */
  spgrid_datum_t *data;  /*< Data hash array */
  double **bounds;
  int dcount_out;
} spgrid_t;

extern const spgrid_t spgrid_nil;

void
spgrid_datum_free(spgrid_datum_t *self);

spgrid_datum_t *
spgrid_datum_new(int dim, int psize,
		 const int *grid_index, 
		 int dcount,
		 const double *payload);
int
spgrid_init(spgrid_t *grid,
	    int dim, int psize,
	    const double *h);


/**
 * Set a boundary. Any datapoint added into the grid will be skipped
 * if outside given bounds, and further, density calculations will
 * crop the grid to the boundary (for proper evaluation of the density
 * near grid boundaries). Should be set before adding any data into
 * the grid, otherwise undefined behaviour.
 */
int
spgrid_set_bound(spgrid_t *grid,
		 int axis,
		 double min,
		 double max);

void
spgrid_destroy(spgrid_t *grid);

/**
 * Compute the size of a volume element containing the point `x`. If
 * no boundaries are set, this is just the product of grid spacings,
 * otherwise takes into account the possible clipping of grid boxes by
 * boundary planes. If `hx` is not null, the containing boxside
 * lengths are copied in here.
 */
double
spgrid_volelem(const spgrid_t *grid, const double *x, double *hx);

int
spgrid_isize(const spgrid_t *grid);

int
spgrid_add_data(spgrid_t *grid,
		const double *x,
		const double *y);

/**
 * Given a sparse density grid `grid`, compute the marginal grid
 * `mgrid` by summing over the data along axis number
 * `axis`. Initializes the `mgrid` input to a new `spgrid` object.
 */
int
spgrid_marginal(spgrid_t *grid,
		spgrid_t *mgrid,
		int axis);

int
spgrid_to_rgrid(const spgrid_t *grid, rgrid_t *rgrid,
		double ***data_out);

int
spgrid_dump(const spgrid_t *grid, FILE *file);

#endif
