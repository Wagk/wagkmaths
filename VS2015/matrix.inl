#ifndef linear_inl__
#define linear_inl__

#include <array>
#include <type_traits>
#include "matrix.h"


/*!************************************************************

**************************************************************/
template<unsigned R, unsigned C /*= R*/, typename T /*= float*/>
Matrix<R, C, T>::Matrix(const std::initializer_list<std::initializer_list<T>>& list)
  : m()
{
  size_t i = 0; size_t j = 0;
  for (const auto& sublist : list)
  {
    for(const auto& elem : sublist)
    {
      m[i][j] = elem;
      ++j;
    }

    ++i;
    j = 0;
  }
}


/*!************************************************************

**************************************************************/
template<unsigned R, unsigned C /*= R*/, typename T /*= float*/>
const T& Matrix<R, C, T>::MatrixRow::operator[](unsigned i) const
{
  return arr[i];
}

/*!************************************************************

**************************************************************/
template<unsigned R, unsigned C /*= R*/, typename T /*= float*/>
T& Matrix<R, C, T>::MatrixRow::operator[](unsigned i)
{
  return arr[i];
}


/*!************************************************************

**************************************************************/
template<unsigned R, unsigned C /*= R*/, typename T /*= float*/>
Matrix<R, C, T>::MatrixRow::MatrixRow(Array& mat)
  : arr(mat)
{

}

/*!************************************************************

**************************************************************/
template<unsigned R, unsigned C /*= R*/, typename T /*= float*/>
typename const Matrix<R, C, T>::MatrixRow Matrix<R, C, T>::operator[](unsigned i) const
{
  return MatrixRow(const_cast<Matrix<R, C, T>* const>(this)->m[i]);
}

/*!************************************************************

**************************************************************/
template<unsigned R, unsigned C /*= R*/, typename T /*= float*/>
typename Matrix<R, C, T>::MatrixRow Matrix<R, C, T>::operator[](unsigned i)
{
  return MatrixRow(m[i]);
}

/*!************************************************************

**************************************************************/
template<unsigned R, unsigned C /*= R*/, typename T /*= float*/>
Matrix<R, C, T>::Matrix()
{
  SetIdentity(*this);
}

/*!************************************************************

**************************************************************/
template<unsigned R, unsigned C, typename T>
void SetIdentity(Matrix<R, C, T>& mat)
{
  for (unsigned i = 0; i < C; ++i)
  {
    for (unsigned j = 0; j < R; ++j)
    {
      if (i == j)
      {
        mat.m[i][j] = static_cast<T>(1);
      }
      else
      {
        mat.m[i][j] = static_cast<T>(0);
      }
    }
  }
}

/*!************************************************************

**************************************************************/
template<unsigned R, unsigned C, typename T>
Matrix<C, R, T> Transpose(const Matrix<R, C, T>& mat)
{
  return Matrix<C, R, T> mat;
}

/*!************************************************************

**************************************************************/
template<unsigned R, unsigned C /*= R*/, typename T /*= float*/>
Matrix<R, C, T>::Matrix(const Array2D& list)
  : m(list)
{

}

/*!************************************************************

**************************************************************/
template<unsigned R, typename T>
std::array<T, R> operator*(const std::array<T, R>& arr, T scalar)
{
  static_assert(std::is_floating_point<T>::value == true, "not real value!");

  std::array<T, R> ret = arr;

  for (auto& elem : ret)
  {
    elem *= scalar;
  }

  return ret;
}

/*!************************************************************

**************************************************************/
template<unsigned R, typename T>
std::array<T, R> operator+(const std::array<T, R>& arr, const std::array<T, R>& addby)
{
  std::array<T, R> ret = arr;

  for (unsigned i = 0; i < R; ++i)
  {
    ret[i] += addby[i];
  }

  return ret;
}

/*!************************************************************

**************************************************************/
template<unsigned R, typename T>
std::array<T, R> operator-(const std::array<T, R>& arr, const std::array<T, R>& subby)
{
  std::array<T, R> ret = arr;

  for (unsigned i = 0; i < R; ++i)
  {
    ret[i] -= subby[i];
  }

  return ret;
}


/*!************************************************************

**************************************************************/
template<unsigned R, typename T>
typename Vector<R, T> operator*(const Matrix<R, R, T>& mat, const Vector<R, T>& vec)
{
  Vector<R, T> ret;

  for (unsigned i = 0; i < R; ++i)
  {
    T value{};

    for (unsigned j = 0; j < R; ++j)
    {
      value += mat[i][j] * vec[j];
    }

    ret[i] = value;
  }

  return ret;
}

/*!************************************************************

**************************************************************/
template<unsigned R, unsigned C, typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<R, C, T>& mat)
{
  os << "[" << std::endl;
  for (size_t i = 0; i < R; ++i)
  {
    os << "\t";
    for (size_t j = 0; j < C; ++j)
    {
      os << mat[i][j] << " ";
    }
    os << std::endl;
  }
  os << "]";

  return os;
}


#endif

