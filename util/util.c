
#include "util.h"

int
iprod(int n, const int *is)
{
  int k, p = 1;
  for (k = 0; k < n; ++k) p *= is[k];
  return p;
}
