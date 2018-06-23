
#include "util.h"

int
iprod(int n, const int *is)
{
  int k, p = 1;
  for (k = 0; k < n; ++k) p *= is[k];
  return p;
}

int
dprod(int n, const double *xs)
{
  int k;
  double p = 1.0;
  for (k = 0; k < n; ++k) p *= xs[k];
  return p;
}
