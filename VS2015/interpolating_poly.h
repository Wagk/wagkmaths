#ifndef interpolating_poly_h__
#define interpolating_poly_h__

#include "math_defines.h"
#include "plot.h"

/*!************************************************************
	FullName	:NevilleMethod
	Returns		:double
	Parameter	:const DataPlot & range
	Parameter	:double x
	Brief:
				
	Assumes:
	Consider:
	Note:
**************************************************************/
double NevilleMethod(const DataPlot& range, double x);

template<typename T>
double NevilleMethod(const Plot2D<T>& range, double x);

/*!************************************************************

**************************************************************/
struct LagrangeBasisPolynomial
{
public:

  using ValueArray = std::vector<double>;

  LagrangeBasisPolynomial(const ValueArray& mx = ValueArray());

  double operator()(double x, size_t i);

private:

  const ValueArray& m_x;

};

/*!************************************************************
	FullName	:CubicSpline
	Returns		:double
	Parameter	:double x       the input to the function
	Parameter	:double ym1pp   the i - 1th second derivative
	Parameter	:double ypp     the ith second derivative
	Parameter	:double yim1    the i - 1th y-value
	Parameter	:double yi      the ith y-value
	Parameter	:double xim1    the i-1th x-value
	Parameter	:double xi      the ith x-value
	Brief:
				
	Assumes:
	Consider:
	Note:
**************************************************************/
double CubicSpline(double x, double ym1pp, double ypp, double yim1, 
  double yi, double xim1, double xi);

#include "interpolating_poly.inl"

#endif