#ifndef fft_h__
#define fft_h__

#include <vector>

namespace Sig //signals processing
{

  template<typename T>
  struct Signal
  {

    Signal(const std::vector& rhs);

    std::vector<T> m_samples;

  };

  template<typename T>
  Signal<T> dft(const Signal<T>& rhs);
  template<typename T>
  Signal<T> inv_dft(const Signal<T>& rhs);
  template<typename T>
  Signal<T> fft(const Signal<T>& rhs);


} //


#endif