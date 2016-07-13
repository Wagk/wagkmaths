#ifndef differentiation_h__
#define differentiation_h__

template<typename Iterator, typename T = double>
T SecondOrderFirstForwardFF(Iterator iter, T h)
{
  return (-*(iter + 2) + 4 * *(iter + 1) - 3 * *iter) / (2 * h);
}

template<typename Iterator, typename T = double>
T ThirdOrderSecondForwardFF(Iterator iter, T h)
{
  return (-*(iter + 3) + 4 * *(iter + 2) - 5 * *(iter + 1) + 2 * *iter) / (h * h);
}

template<typename Iterator, typename T = double>
T SecondOrderFirstBackwardFF(Iterator iter, T h)
{
  return (3 * *iter - 4 * *(iter - 1) + *(iter - 2)) / (2 * h);
}

template<typename Iterator, typename T = double>
T ThirdOrderSecondBackwardFF(Iterator iter, T h)
{
  return (2 * *iter - 5 * *(iter - 1) + 4 * *(iter - 2) - *(iter - 3)) / (h * h);
}

template<typename Iterator, typename T = double>
T FirstOrderFirstCentralFF(Iterator iter, T h)
{
  return (*(iter + 1) - *(iter - 1)) / (2 * h);
}

template<typename Iterator, typename T = double>
T FirstOrderSecondCentralFF(Iterator iter, T h)
{
  return (*(iter + 1) - 2 * *iter + *(iter - 1)) / (h * h);
}



#endif

