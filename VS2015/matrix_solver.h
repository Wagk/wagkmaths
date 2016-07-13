#ifndef matrix_solver_h__
#define matrix_solver_h__

#include "matrix.h"

template<unsigned R, typename T = double>
class GaussJordan
{
public:

  using AugmentedMat = std::pair<Matrix<R, R, T>, Matrix<R, R, T>>;

  static Matrix<R, R, T> FindInverse(const Matrix<R, R, T>& mat);

private:

  /*!************************************************************
  	FullName	:GaussJordan<T>::EliminateColumn
  	Returns		:void
  	Parameter	:unsigned col, the column to zero out
  	Parameter	:unsigned i, the row where the i is 1;
  	Brief:
  				zeroes out the column using gauss-jordan operations
  	Assumes:
  	Consider:
  	Note:
  **************************************************************/
  static void EliminateColumn(const unsigned col, const unsigned i, AugmentedMat& mat);

  static void ScaleRow(AugmentedMat& mat, unsigned i, T scalar);
  static void SwapRow(AugmentedMat& mat, unsigned i, unsigned j);
  static void AddRow(AugmentedMat& mat, unsigned target, unsigned i);
  static void SubRow(AugmentedMat& mat, unsigned target, unsigned i);

};



#include "matrix_solver.inl"

#endif

