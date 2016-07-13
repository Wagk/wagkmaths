#ifndef interpolating_poly_inl__
#define interpolating_poly_inl__

#include "interpolating_poly.h"

template<typename T>
static double NevilleRecursive(size_t i, size_t j, double x, const Plot2D<T>& range)
{
  if (i == j) return range.y[i];
  else
  {
    double xj = range.x[j];
    double xi = range.x[i];

    double pijm1 = NevilleRecursive(i, j - 1, x, range);
    double pip1j = NevilleRecursive(i + 1, j, x, range);

    return ((xj - x) * pijm1 + (x - xi) * pip1j) / (xj - xi);
  }
}

template<typename T>
double NevilleMethod(const Plot2D<T>& range, double x)
{
  return NevilleRecursive(0, range.Size() - 1, x, range);
}


#endif

