
#ifndef SMAUX_H
#define SMAUX_H

/* Auxilliary functions for sparse matrices, such as constructing
 * diagonal and derivative matrices.*/

#include <cs.h>

#include "util.h"
#include "rgrid.h"

#define SPALLOC_VALUES 1
#define SPALLOC_TRIPLET 1
#define SPALLOC_CCS 0

/**
 * Helper function: Create reformatted copy of a triplet format sparse
 * matrix, and deallocate the old version. Essentially wraps 
 * `cs_compress` function to convert triple to compressed format.
 */
int
compress_(cs **m);

/**
 * Construct a diagonal sparse matrix from a given vector of
 * values. If `triplet == 1` then result is in triplet format,
 * otherwise it is converted to compressed column.
 */
cs *
diag(int n, const double *data, int stride, int triplet);

/**
 * Construct an identity sparse matrix. If `triplet == 1` then result
 * is in triplet format, otherwise it is converted to compressed
 * column.
 */
cs *
eye(int n, int triplet);

/**
 * Construct differentiation matrix along axis `axis` given grid
 * specification `grid`. If parameter `triplet == 1`, then output in
 * triplet format, otherwise converted to compressed column.  
 */
cs *
deriv(const rgrid_t *grid, int axis, int triplet);

/**
 * Rearrange `cs` column compressed matrix so that rows are stored in
 * ascending order.
 */
int
rowsort(cs *A);

#endif
