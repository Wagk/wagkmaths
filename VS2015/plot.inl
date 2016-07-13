#ifndef plot_inl__
#define plot_inl__

#include "plot.h"
#include "guarantees.h"

/*!************************************************************

**************************************************************/
template<typename T /*= double*/>
Plot2D<T>::Plot2D(const std::initializer_list<std::initializer_list<T>>& pairplot)
{
  T arr[2];
  for (const auto& subplot : pairplot)
  {
    size_t i = 0;
    for (const auto& elem : subplot)
    {
      arr[i] = elem;
      ++i;
    }
    x.push_back(arr[0]);
    y.push_back(arr[1]);
  }

  ENSURE(x.size() == y.size());
}

/*!************************************************************

**************************************************************/
template<typename T /*= double*/>
Plot2D<T>::Plot2D(const std::initializer_list<T>& xplot, const std::initializer_list<T>& yplot)
  : x(xplot)
  , y(yplot)
{
  ENSURE(x.size() == y.size());
}

/*!************************************************************

**************************************************************/
template<typename T /*= double*/>
Plot2D<T>::Plot2D(const std::vector<T>& xplot, const std::vector<T>& yplot)
  : x(xplot)
  , y(yplot)
{
  ENSURE(x.size() == y.size());
}

/*!************************************************************

**************************************************************/
template<typename T /*= double*/>
T Plot2D<T>::y_avg() const
{
  T avg{};
  for (const auto& elem : x)
  {
    avg += elem;
  }
  return avg / x.size();
}

/*!************************************************************

**************************************************************/
template<typename T /*= double*/>
T Plot2D<T>::x_avg() const
{
  T avg{};
  for (const auto& elem : y)
  {
    avg += elem;
  }
  return avg / x.size();
}

/*!************************************************************

**************************************************************/
template<typename T /*= double*/>
T Plot2D<T>::sum_each_y() const
{
  return sum_each_y([](double y) {return y; });
}

/*!************************************************************

**************************************************************/
template<typename T /*= double*/>
T Plot2D<T>::sum_each_x() const
{
  return sum_each_x([](double x) {return x; });
}

/*!************************************************************

**************************************************************/
template<typename T>
template<typename Functor>
T Plot2D<T>::sum_each_xy_pair(const Functor& func) const
{
  ENSURE(x.size() == y.size());

  T val{};

  for (size_t i = 0; i < x.size(); ++i)
  {
    val += func(x[i], y[i]);
  }

  return val;
}

/*!************************************************************

**************************************************************/
template<typename T>
template<typename Functor>
T Plot2D<T>::sum_each_x(const Functor& func) const
{
  T val{};

  for (size_t i = 0; i < x.size(); ++i) 
  {
    val += func(x[i]);
  }

  return val;
}

/*!************************************************************

**************************************************************/
template<typename T>
template<typename Functor>
T Plot2D<T>::sum_each_y(const Functor& func) const
{
  T val{};

  for (size_t i = 0; i < x.size(); ++i)
  {
    val += func(y[i]);
  }

  return val;
}

/*!************************************************************

**************************************************************/
template<typename T /*= double*/>
typename Plot2D<T>::SizeType Plot2D<T>::Size() const
{
  ENSURE(x.size() == y.size());

  return x.size();
}

#endif

