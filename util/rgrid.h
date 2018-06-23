
#ifndef RGRID_H
#define RGRID_H

#include "util.h"

/**
 * Regular multidimensional grid object. Defines grid dimension,
 * numbers of samples along each dimension, and their start and end
 * points. Grid data is not stored in these structures.
 */
typedef struct {
  int dim;
  int *ns;
  double *bs;

  double *hs;
  int *strides;
  int m;		/*< Total number of points in the grid. */
} rgrid_t;

/**
 * Null grid for initializing stack allocated grid objects. Safe pass
 * into `rgrid_destroy`.
 */
#define RGRID_NIL {0, NULL, NULL, NULL, NULL, 0}
extern const rgrid_t rgrid_nil;

/**
 * Deallocate data associated with a grid object. Frees only objects
 * contained in the argument.
 */
void
rgrid_destroy(rgrid_t *self);

/**
 * Setup a multidimensional regular grid with constant step
 * size. Dimension is `dim`, numbers of sample points along each
 * dimension are in the `dim` sized array `ns`, and the bounds are in
 * the `2*dim` sized array `bs`, stored as a sequence of lower bound,
 * upper bound pairs.  The object should be released using
 * `rgrid_destroy` once done using it.
 */
int
rgrid_init(rgrid_t *self, int dim, const int *ns, const double *bs);

/**
 * Helper function: Initialize a vector of grid indices for looping
 * through grid points (sets all elements to zero). If pointer `x` is
 * not `NULL`, also the grid locations are computed (the lower bounds
 * along each dimension).
 */

void
rgrid_start(const rgrid_t *self, int *i, double *x);

/**
 * Helper function: Given a vector of grid indices, advance them one
 * step forward (last index fastest varying), and return 0 when rolled
 * back to beginning. If pointer `x` is not `NULL`, the corresponding
 * grid locations are computed into that array.
 */
int
rgrid_step(const rgrid_t *self, int *i, double *x);

/**
 * Compute index from a multi-index `i`.
 */
int
rgrid_rindex(const rgrid_t *self, int *i);

/**
 * Compute the size of a volume element.
 */
double
rgrid_volel(const rgrid_t *self);

/**
 * Find an approximating grid index to given point. Point selected so
 * that uniformly distributed x will give uniformly distributed i.
 */
int
rgrid_nearest_i(const rgrid_t *self, const double *x, int *i);

/**
 * Print data evaluated on a grid. Parameter `grid` holds the grid
 * description, `n` is the number of arrays to print, `y` is pointer
 * to the array of arrays holding the data to be printed, `headers` is
 * a list of strings to be used as headers in the output (can be
 * NULL). If `x` is not null, then this array is used as the grid
 * coordinates rather than the regular grid points implied by the grid
 * object.
 */
int
rgrid_print_data(const rgrid_t *grid, double **x, int n, double **y,
		 const char **headers);

/**
 * Construct marginal data grids. Given grid spec `grid`, and `n`
 * arrays of data evaluated on that grid, stored in `ys`, integrate
 * along axis `axis` and construct new grid `mgrid` with dimension
 * `axis` dropped out and results of integration allocated to arrays
 * pointed by `us_out`. If `xs` is non-null, this array is used for
 * the grid coordinates on the integrated axis instead of the regular
 * grid implied by the `grid` argument (array covers the whole grid,
 * i.e. it has `grid.m` elements). The output array and grid should be
 * completely uninitialized, this function will do all the allocation;
 * deallocation is caller's responsibility.
 */
int
rgrid_marginals(const rgrid_t *grid, const double *x, int n, double **ys, int axis,
		rgrid_t *mgrid, double ***us_out);


/* Example memory layout for 3x3x3 grid. The last coordinate varies
 * fastest:
 *
 *  ------
 *      z0
 *    y0z1
 *      z2
 *      --
 *      z0
 *  x0y1z1
 *      z2
 *      --
 *      z0
 *    y2z1
 *      z2
 *  ------
 *      z0
 *    y0z1
 *      z2
 *      --
 *      z0
 *  x1y1z1
 *      z2
 *      --
 *      z0
 *    y2z1
 *      z2
 *  ------
 *      z0
 *    y0z1
 *      z2
 *      --
 *      z0
 *  x2y1z1
 *      z2
 *      --
 *      z0
 *    y2z1
 *      z2
 *  ------
 *
 */

#endif
