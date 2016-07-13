#ifndef vector_inl__
#define vector_inl__

#include "vector.h"

/*!************************************************************

**************************************************************/
template<unsigned D, typename T /*= float*/>
Vector<D, T>::Vector(const std::initializer_list<T>& list)
{
  size_t i = 0;
  for (const auto& elem : list)
  {
    m[i] = elem;
    ++i;
  }
}

/*!************************************************************

**************************************************************/
template<unsigned D, typename T /*= float*/>
const T& Vector<D, T>::operator[](size_t index) const
{
  return m[index];
}

/*!************************************************************

**************************************************************/
template<unsigned D, typename T /*= float*/>
T& Vector<D, T>::operator[](size_t index)
{
  return m[index];
}

/*!************************************************************

**************************************************************/
template<unsigned D, typename T /*= float*/>
Vector<D, T>::Vector()
  : m{0}
{

}

/*!************************************************************

**************************************************************/
template<unsigned D, typename T, typename ReturnType /*= T*/>
ReturnType operator*(const Vector<D, T>& v1 , const Vector<D, T>& v2)
{
  ReturnType ret{};

  for (decltype(D) i = 0; i < D; ++i)
  {
    ret += (v1[i] * v2[i]);
  }

  return ret;
}

/*!************************************************************

**************************************************************/
template<unsigned D, typename T, typename InputT>
Vector<D, T>
operator*(const Vector<D, T>& vec, InputT scale)
{
  Vector<D, T> temp = vec;
  temp *= scale;
  return temp;
}

/*!************************************************************

**************************************************************/
template<unsigned D, typename T, typename InputT>
Vector<D, T> operator*(InputT scale, const Vector<D, T>& vec)
{
  return vec * scale;
}

/*!************************************************************

**************************************************************/
template<unsigned D, typename T, typename InputT>
Vector<D, T>& operator*=(Vector<D, T>& vec, InputT scale)
{
  for (size_t i = 0; i < D; ++i)
  {
    vec[i] *= scale;
  }

  return vec;
}

/*!************************************************************

**************************************************************/
template<unsigned D, typename T>
std::ostream& operator<<(std::ostream& os, const Vector<D, T>& vec)
{
  os << "[ ";
  for (size_t i = 0; i < D; ++i)
  {
    os << vec[i] << " ";
  }
  os << "]";

  return os;
}



#endif

