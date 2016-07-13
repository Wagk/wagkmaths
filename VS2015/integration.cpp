#include "integration.h"

#include "interpolating_poly.h"

/*!************************************************************
  FullName	:RecursiveTrapezoidHelperBase
  Returns		:double
  Parameter	:const Plot2D<double> & data
  Parameter	:double start
  Parameter	:double end
  Brief:

  Assumes:
  Consider:
  Note:
**************************************************************/
static double RecursiveTrapezoidHelperBase(
  const Plot2D<double>& data,
  double start, double end)
{
  double coeff = (end - start) / 2.0;
  double fa = NevilleMethod(data, start);
  double fb = NevilleMethod(data, end);

  return coeff * (fa + fb);
}

/*!************************************************************
  FullName	:RecursiveTrapezoidHelper
  Returns		:double
  Parameter	:const Plot2D<double> & data
  Parameter	:double start
  Parameter	:double end
  Parameter	:double previter
  Parameter	:size_t iter
  Brief:

  Assumes:
  Consider:
  Note:
**************************************************************/
static double RecursiveTrapezoidHelper(
  const Plot2D<double>& data,
  double start, double end, double previter, size_t iter)
{

  double subdiv = std::pow(2.0, iter - 1);
  double coeff = (end - start) / subdiv;

  double sum = 0.0;

  size_t stopcount = std::pow(2, iter - 2);
  for (size_t i = 1; i <= stopcount; ++i)
  {

    double a = (subdiv - 2 * i + 1) * start;
    double b = (2 * i - 1) * end;

    double arg = (a + b) / subdiv;

    sum += NevilleMethod(data, arg);

  }

  double prev = previter * 0.5;
  return prev + coeff * sum;

}

/*!************************************************************

**************************************************************/
double RecursiveTrapezoidMethod(const Plot2D<double>& data, double start, double end, size_t iters)
{
  if (iters == 0) return 0.0;
  if (iters == 1) return RecursiveTrapezoidHelperBase(data, start, end);

  double previter = RecursiveTrapezoidHelperBase(data, start, end);
  for (size_t i = 2; i < iters; ++i)
  {
    previter = RecursiveTrapezoidHelper(data, start, end, previter, i);
  }

  return previter;
}

/*!************************************************************

**************************************************************/
double RombergIntegration(const Plot2D<double>& data, double start, double end, size_t degrees)
{
  return 0.0;
}

/*!************************************************************

**************************************************************/
double CompositeTrapezoidMethod(double start, double end, const std::vector<double>& data)
{
  double h = (end - start) / data.size();

  double startnend = (h / 2) * (data.front() + data.back());

  double middlin = 0.0;
  for (size_t i = 1; i < data.size() - 1; ++i)
  {
    middlin += data[i];
  }
  middlin *= h;

  return startnend + middlin;
}

