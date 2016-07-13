#ifndef integration_h__
#define integration_h__

#include "interpolating_poly.h"

//first find a way to generate a function out of the data
double RecursiveTrapezoidMethod(const Plot2D<double>& data, double start, double end, size_t iters);

double RombergIntegration(const Plot2D<double>& data, double start, double end, size_t degrees);

//data starts at `start`, and ends at `end`
double CompositeTrapezoidMethod(double start, double end, const std::vector<double>& data);



#endif

