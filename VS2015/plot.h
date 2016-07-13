#ifndef plot_h__
#define plot_h__

#include <vector>
#include <array>


/*
 *  I dream that once day it will be generalised to N dimensional
 **/

template<unsigned N, typename T = double>
class Plot
{
public:
private:

  std::array<std::vector<T>, N> m_plotlist;

};

template<typename T = double>
class Plot2D
{
public:

  using SizeType = size_t;

  Plot2D() = default;
  Plot2D(const std::vector<T>& xplot, const std::vector<T>& yplot);
  Plot2D(const std::initializer_list<T>& xplot, const std::initializer_list<T>& yplot);
  Plot2D(const std::initializer_list<std::initializer_list<T>>& pairplot);

  std::vector<T> x;
  std::vector<T> y;

  T x_avg() const;
  T y_avg() const;

  SizeType Size() const;

  //TODO: implement this someday
  //std::pair<T, T> operator[](size_t index) const;


  //functors have to follow this format: T foo(T xi, T yi);
  template<typename Functor>
  T sum_each_xy_pair(const Functor& func) const;
  //functors have to follow this format: T foo(T xi);
  template < typename Functor>
  T sum_each_x(const Functor& func) const;
  //functors have to follow this format: T foo(T xi);
  template < typename Functor>
  T sum_each_y(const Functor& func) const;

  //these just sum the things
  T sum_each_x() const;
  T sum_each_y() const;
};






#include "plot.inl"

#endif

