
#ifndef UTIL_H
#define UTIL_H

#define OK 0
#define OUT_OF_MEM 1
#define INVALID 2
#define FAILED 3
#define NO_CONVERGENCE 4

/**
 * Check `expr` (e.g. a `malloc`) for false result and set status
 * flag and goto exit in case of failure. Assumes `status` is defined
 * as `int` in the current context and `exit` label marks the cleanup
 * and return part of the function.
 */
#define CHECK(expr, errnum) do {		\
    if (expr) { status = errnum; goto exit; }	\
  } while (0)


/**
 * Check `expr` (e.g. a `malloc`) for `NULL` result and set status
 * flag and goto exit in case of failure. Assumes `status` is defined
 * as `int` in the current context and `exit` label marks the cleanup
 * and return part of the function.
 */
#define CHECK_NULL(expr, errnum) do {		\
    typeof(expr) _cn_expr = (expr);		\
    if (_cn_expr == NULL)			\
      { status = (errnum); goto exit; }		\
  } while (0)

/**
 * Free macro that checks for `NULL` pointers and uses custom free
 * function.
 */
#define Xfree(ptr) do {				\
    typeof(ptr) _xf_ptr = (ptr);		\
    if (_xf_ptr) free(_xf_ptr);			\
  } while (0)

/**
 * Free macro that checks for `NULL` pointers and uses custom free
 * function.
 */
#define Xfree_(ptr, freef) do {			\
    typeof(ptr) _xf_ptr = (ptr);		\
    typeof(freef) *_xf_freef = &(freef);	\
    if (_xf_ptr) _xf_freef(_xf_ptr);		\
  } while (0)
  
/**
 * Helper macro for deallocating arrays of pointers. Takes an array
 * `ptrs` with `n` elements and calls the function `freef` on each of
 * the elements, provided they are non-null. Finally `ptrs` is
 * deallocated using `free`.
 */
#define Xfree_many(ptrs, n, freef) do {				\
    typeof(ptrs) _xf_ptrs = (ptrs);				\
    size_t _xf_n = (n), _xf_i;					\
    typeof(freef) *_xf_freef = &(freef);			\
    if (_xf_ptrs)						\
      {								\
	for (_xf_i = 0; _xf_i < _xf_n; ++_xf_i)			\
	  if (_xf_ptrs[_xf_i]) _xf_freef(_xf_ptrs[_xf_i]);	\
	free(_xf_ptrs);						\
      }								\
  } while (0)							\

/**
 * Helper function, compute the product of `n` integers stored in
 * array `is`.
 */
int
iprod(int n, const int *is);


#endif
