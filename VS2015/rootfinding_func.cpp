#include "rootfinding_func.h"

double Bisection(double begin, double end, double epsilon, std::function<Polynomial> func)
{
  double e = epsilon;
  double a = begin;
  double b = end;
  double c = 0;

  while ((b - a) > (2 * e))
  {
    c = (a + b) / 2;

    if (func(a) * func(c) > 0)
    {
      a = c;
    }
    else
    {
      b = c;
    }
  }

  return (a + b) / 2;
}

double NewtonRaphson(double start, double epsilon, std::function<Polynomial> func, std::function<Polynomial> derivative)
{
  double e = epsilon;
  double x1 = start;
  double x2 = x1;

  do
  {
    x2 = x1 - func(x1) / derivative(x1);

    if (std::abs(x2 - x1) <= epsilon)
    {
      return x1;
    }
    else
    {
      x1 = x2;
    }
  } while (true);
}
