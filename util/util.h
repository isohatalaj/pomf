
#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>

#define OK 0
#define OUT_OF_MEM 1
#define INVALID 2
#define FAILED 3
#define NO_CONVERGENCE 4

#define ERRMSG(format, ...) do {				\
    const char *fname = __FILE__;				\
    const char *fshort = fname;					\
    while (*fshort != '\0') fshort++;				\
    while ((fshort-1) != fname && *(fshort-1) != '/') fshort--;	\
    fprintf(stderr, "# ERROR @ %s(%d): " format,		\
	    fshort, __LINE__, ##__VA_ARGS__);			\
  } while (0)

#define FAILWITH(code) do {			\
    status = (code);				\
    goto exit;					\
  } while (0)

/**
 * Check `expr` (e.g. a `malloc`) for false result and set status
 * flag and goto exit in case of failure. Assumes `status` is defined
 * as `int` in the current context and `exit` label marks the cleanup
 * and return part of the function.
 */
#define CHECK(expr, errnum) do {				\
    if (expr)							\
      {								\
	ERRMSG("CHECK FAILED: %s\n", #expr);			\
	status = errnum;					\
	goto exit;						\
      }								\
  } while (0)

/**
 * Same as `CHECK` but use result of `expr` as `errnum`.
 */
#define CHECK_(expr) do {					\
    if ((status = (expr)))					\
      {								\
	ERRMSG("CHECK FAILED: %s\n", #expr);			\
	goto exit;						\
      }								\
  } while (0)

/**
 * Check `expr` (e.g. a `malloc`) for `NULL` result and set status
 * flag and goto exit in case of failure. Assumes `status` is defined
 * as `int` in the current context and `exit` label marks the cleanup
 * and return part of the function.
 */
#define CHECK_NULL(expr, errnum) \
  CHECK((expr) == NULL, errnum)



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
 * Safe maximum macro.
 */
#define MAX(a, b) ({				\
      typeof(a) _a = (a);			\
      typeof(b) _b = (b);			\
      _a > _b ? _a : _b; })

/**
 * Safe minimum macro.
 */
#define MIN(a, b) ({				\
      typeof(a) _a = (a);			\
      typeof(b) _b = (b);			\
      _a < _b ? _a : _b; })


/**
 * Helper function, compute the product of `n` integers stored in
 * array `is`.
 */
int
iprod(int n, const int *is);

/**
 * Helper function, compute the product of `n` doubles stored in
 * array `xs`.
 */
int
dprod(int n, const double *xs);


#endif

