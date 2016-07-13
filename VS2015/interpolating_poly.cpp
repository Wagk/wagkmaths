#include "interpolating_poly.h"
#include <cassert>

static double NevilleRecursive(size_t i, size_t j, double x, const DataPlot& range)
{
  if (i == j) return range[i].second;
  else
  {
    double xj = range[j].first;
    double xi = range[i].first;

    double pijm1 = NevilleRecursive(i, j - 1, x, range);
    double pip1j = NevilleRecursive(i + 1, j, x, range);

    return ((xj - x) * pijm1 + (x - xi) * pip1j) / (xj - xi);
  }
}

double NevilleMethod(const DataPlot& range, double x)
{
  return NevilleRecursive(0, range.size() - 1, x, range);
}

/*!************************************************************

**************************************************************/
double CubicSpline(double x, double ym1pp, double ypp, double yim1, double yi, double xim1, double xi)
{
  double ximxim1 = xi - xim1;
  double xmxim1 = x - xim1;
  double ximx = xi - x;

  double v1 = (ypp   * (xmxim1)*((xmxim1 * xmxim1) - (ximxim1 * ximxim1))) / (6 * ximxim1);
  double v2 = (ym1pp * (ximx  )*((ximx   * ximx  ) - (ximxim1 * ximxim1))) / (6 * ximxim1);
  double v3 = (yi * xmxim1 + yim1 * ximx) / ximxim1;

  return v1 + v2 + v3;
}

/*!************************************************************

**************************************************************/
LagrangeBasisPolynomial::LagrangeBasisPolynomial(const ValueArray& mx /*= ValueArray()*/)
  : m_x(mx)
{

}

/*!************************************************************

**************************************************************/
double LagrangeBasisPolynomial::operator()(double x, size_t index)
{
  double val{1};
  double xj = m_x[index];

  for (size_t i = 0; i < m_x.size(); ++i)
  {
    if (i == index) continue;

    double xm = m_x[i];

    val *= (x - xm) / (xj - xm);
  }

  return val;
}
