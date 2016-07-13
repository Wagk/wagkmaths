#ifndef math_defines_h__
#define math_defines_h__

#include <vector>
#include <functional>

using Polynomial = double(double);
using Roots = std::vector<double>;
using Values = std::vector<double>;
using Range = std::pair<double, double>;
using Ranges = std::vector<Range>;
using DDCoefficients = std::vector<double>;
using DataPlot = std::vector<std::pair<double, double>>; //first is x, second is y

template<typename T>
using DataSet = std::vector<T>;

const double PI = 3.141592653589793238463;


#endif

