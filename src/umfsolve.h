
#ifndef UMFSOLVE_H
#define UMFSOLVE_H

#include <cs.h>

#include "util.h"
#include "umfpack.h"
#include "smaux.h"

typedef struct {
  void *Numeric;
  void *Symbolic;
  double Info[UMFPACK_INFO];
  double Control[UMFPACK_CONTROL];
} umfsolver_t;

#define UMFSOLVER_NIL {NULL, NULL, {0}, {0}}

int
umfsolver_init(umfsolver_t *self, cs *A);

void
umfsolver_destroy(umfsolver_t *self);

int
umfsolve_(umfsolver_t *self, cs *A, double *b, double *x);

int
umfsolve(cs *A, double *b, double *x);

#endif
