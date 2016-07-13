#ifndef rootfinding_func_h__
#define rootfinding_func_h__

#include "math_defines.h"

double Bisection(double begin, double end, double epsilon, std::function<Polynomial> func);
double NewtonRaphson(double start, double epsilon, std::function<Polynomial> func, std::function<Polynomial> derivative);

#endif



