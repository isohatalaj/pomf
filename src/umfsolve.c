
#include "umfsolve.h"

int
umfsolve(cs *A, double *b, double *x)
{
  int status = OK;

  int i;
  const int n = A->n;
  void *Numeric = NULL, *Symbolic = NULL;
  double Info[UMFPACK_INFO], Control[UMFPACK_CONTROL];

  CHECK( A->n != A->m, INVALID );
  CHECK( rowsort(A), FAILED );

  umfpack_di_defaults(Control);

  // fprintf(stderr, "# SYMBOLIC\n"); fflush(stderr);

  CHECK( umfpack_di_symbolic(n, n, A->p, A->i, A->x,
			     &Symbolic, Control, Info) < 0,
	 FAILED );

  // fprintf(stderr, "# NUMERIC\n"); fflush(stderr);

  CHECK( umfpack_di_numeric(A->p, A->i, A->x, Symbolic, &Numeric,
			    Control, Info) < 0,
	 FAILED );

  // fprintf(stderr, "# SOLVE\n"); fflush(stderr);

  CHECK( umfpack_di_solve(UMFPACK_A, A->p, A->i, A->x, x, b, Numeric,
			  Control, Info) < 0,
	 FAILED );

 exit:
  if (Numeric) umfpack_di_free_numeric(&Numeric);
  if (Symbolic) umfpack_di_free_symbolic(&Symbolic);

  return status;
}

int
umfsolver_init(umfsolver_t *self, cs *A)
{
  int status = OK;
  const int n = A->n;
  
  self->Numeric = NULL;
  self->Symbolic = NULL;

  CHECK( A->n != A->m, INVALID );
  CHECK( rowsort(A), FAILED );

  umfpack_di_defaults(self->Control);

  CHECK( umfpack_di_symbolic(n, n, A->p, A->i, A->x,
			     &self->Symbolic,
			     self->Control,
			     self->Info) < 0,
	 FAILED );

  CHECK( umfpack_di_numeric(A->p, A->i, A->x,
			    self->Symbolic, &self->Numeric,
			    self->Control, self->Info) < 0,
	 FAILED );

 exit:
  if (status) umfsolver_destroy(self);

  return status;  
}

void
umfsolver_destroy(umfsolver_t *self)
{
  if (self == NULL) return;
  if (self->Numeric) umfpack_di_free_numeric(&self->Numeric);
  if (self->Symbolic) umfpack_di_free_symbolic(&self->Symbolic);
}

int
umfsolve_(umfsolver_t *self, cs *A, double *b, double *x)
{
  int status = OK;
  
  CHECK( umfpack_di_solve(UMFPACK_A, A->p, A->i, A->x, x, b,
			  self->Numeric, 
			  self->Control, self->Info) < 0,
	 FAILED );

 exit:

  return status;
}
