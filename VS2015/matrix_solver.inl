#ifndef matrix_solver_inl__
#define matrix_solver_inl__

#include "matrix.h"
#include <cassert>
#include <type_traits>
#include <limits>

template<typename T>
static bool NonZero(T val)
{
  static_assert(std::is_floating_point<T>::value == true, "not real number!");

  return std::abs(val) < std::numeric_limits<T>::epsilon() == false;
}

/*!************************************************************

**************************************************************/
template<unsigned R, typename T /*= float*/>
Matrix<R, R, T> GaussJordan<R, T>::FindInverse(const Matrix<R, R, T>& in)
{
  AugmentedMat pair = std::make_pair(in, Matrix<R, R, T>());

  for (unsigned i = 0; i < R; ++i)
  {
    EliminateColumn(i, i, pair);
  }

  return pair.second;
}


/*!************************************************************

**************************************************************/
template<unsigned R, typename T /*= float*/>
void GaussJordan<R, T>::EliminateColumn(const unsigned col, const unsigned row, AugmentedMat& mat)
{
  //find a row in the column that is nonzero
  unsigned nonzero_row = 0;
  for (unsigned i = row; i < R; ++i)
  {
    if (NonZero(mat.first.m[i][col]) == true)
    {
      nonzero_row = i;
      break;
    }
  }
  EXPECT(nonzero_row < R);

  //once found, swap it to row we want
  SwapRow(mat, row, nonzero_row);

  //find the scalar that will normalize it to 1;
  T scalar = mat.first.m[nonzero_row][col];

  ScaleRow(mat, nonzero_row, 1/scalar);

  //for every row that isn't our row, scale and subtract a la Gauss Jordan
  for (unsigned i = 0; i < R; ++i)
  {
    if (i == nonzero_row || NonZero(mat.first.m[i][col]) == false) continue;

    T scale = mat.first.m[i][col];

    EXPECT(NonZero(scale));

    ScaleRow(mat, nonzero_row, scale);

    SubRow(mat, i, nonzero_row);

    ScaleRow(mat, nonzero_row, 1 / scale);
  }
  //column should be entirely zero except for ours, which should be 1
}

/*!************************************************************

**************************************************************/
template<unsigned R, typename T /*= float*/>
void GaussJordan<R, T>::AddRow(AugmentedMat& mat, unsigned target, unsigned i)
{
  mat.first.m[target] = mat.first.m[target] + mat.first.m[i];
  mat.second.m[target] = mat.second.m[target] + mat.second.m[i];
}

/*!************************************************************

**************************************************************/
template<unsigned R, typename T /*= float*/>
void GaussJordan<R, T>::SubRow(AugmentedMat& mat, unsigned target, unsigned i)
{
  mat.first.m[target] = mat.first.m[target] - mat.first.m[i];
  mat.second.m[target] = mat.second.m[target] - mat.second.m[i];
}

/*!************************************************************

**************************************************************/
template<unsigned R, typename T /*= float*/>
void GaussJordan<R, T>::SwapRow(AugmentedMat& mat, unsigned i, unsigned j)
{
  std::swap(mat.first.m[i], mat.first.m[j]);
  std::swap(mat.second.m[i], mat.second.m[j]);
}

/*!************************************************************

**************************************************************/
template<unsigned R, typename T /*= float*/>
void GaussJordan<R, T>::ScaleRow(AugmentedMat& mat, unsigned i, T scalar)
{
  static_assert(std::is_floating_point<T>::value == true, "not real number!");

  mat.first.m[i] = mat.first.m[i] * scalar;
  mat.second.m[i] = mat.second.m[i] * scalar;
}

#endif

