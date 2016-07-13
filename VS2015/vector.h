#ifndef vector_h__
#define vector_h__

#include <array>

template<unsigned D, typename T = double>
struct Vector
{
public:

  using ArrayType = std::array<T, D>;
  static const unsigned Size = D;

  Vector();
  Vector(const std::initializer_list<T>& list);

  T& operator[](size_t index);
  const T& operator[](size_t index) const;

  ArrayType m;
};

using Vec4 = Vector<4, float>;
using Vec3 = Vector<3, float>;
using Vec2 = Vector<2, float>;

template<unsigned D, typename T, typename ReturnType = T>
ReturnType operator*(const Vector<D, T>& v1, const Vector<D, T>& v2);

template<unsigned D, typename T, typename InputT>
Vector<D, T> operator*(const Vector<D, T>& vec, InputT scale);

template<unsigned D, typename T, typename InputT>
Vector<D, T> operator*(InputT scale, const Vector<D, T>& vec);

template<unsigned D, typename T, typename InputT>
Vector<D, T>& operator*=(Vector<D, T>& vec, InputT scale);

template<unsigned D, typename T>
std::ostream& operator<<(std::ostream& os, const Vector<D, T>& vec);


#include "vector.inl"

#endif




