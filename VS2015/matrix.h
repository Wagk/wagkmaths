#ifndef linear_h__
#define linear_h__

#include <array>
#include "vector.h"

#include <iostream>


/*
 *  The matrix is laid out:
 
    [0][0] [0][1] [0][2] [0][3]
    [1][0] [1][1] [1][2] [1][3]
    [2][0] [2][1] [2][2] [2][3]
    [3][0] [3][1] [3][2] [3][3]

 **/
template<unsigned R, unsigned C = R, typename T = double>
struct Matrix
{
public:

  using Array = std::array<T, C>;
  using Array2D = std::array<std::array<T, C>, R>;

  struct MatrixRow
  {
    MatrixRow(Array&);

    T& operator[](unsigned i);
    const T& operator[](unsigned i) const;

    Array& arr;
  };


  Matrix();
  Matrix(const Array2D& arr);
  Matrix(const std::initializer_list<std::initializer_list<T>>& list);
  Matrix& operator=(const Matrix&) = default;

  MatrixRow operator[](unsigned i);
  const MatrixRow operator[](unsigned i) const;

  std::array<std::array<T, C>, R> m; // m[R][C]

};

template<unsigned R, typename T>
Vector<R, T> operator*(const Matrix<R, R, T>& mat, const Vector<R, T>& vec);

using Mat2 = Matrix<2, 2, float>;
using Mat3 = Matrix<3, 3, float>;

template<unsigned R, unsigned C, typename T>
void SetIdentity(Matrix<R, C, T>& mat);

template<unsigned R, unsigned C, typename T>
Matrix<C, R, T> Transpose(const Matrix<R, C, T>& mat);

template<unsigned R, typename T>
std::array<T, R> operator*(const std::array<T, R>& arr, T scalar);

template<unsigned R, typename T>
std::array<T, R> operator+(const std::array<T, R>& arr, const std::array<T, R>& addby);

template<unsigned R, typename T>
std::array<T, R> operator-(const std::array<T, R>& arr, const std::array<T, R>& subby);

template<unsigned R, unsigned C, typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<R, C, T>& mat);

//row swapping should work with std::swap

#include "matrix.inl"

#endif
