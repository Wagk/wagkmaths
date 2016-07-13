#include "questions.h"
#include <cassert>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>

#include "interpolating_poly.h"
#include "matrix.h"
#include "matrix_solver.h"
#include "vector.h"
#include <algorithm>

#include "plot.h"

#include "differentiation.h"
#include "integration.h"

#include "math_defines.h"

/*!************************************************************
  FullName	:Chebyshev_first_kind
  Returns		:double
  Parameter	:double x
  Parameter	:size_t degree
  Brief:
        This is the chebyshev polynomial of the first kind,
        represented in trigonometric terms
  Assumes:
  Consider:
  Note:
**************************************************************/
double Chebyshev_first_kind(double x, size_t n)
{
  return std::cos(n * std::acos(x));
}

/*!************************************************************
  FullName	:Chebyshev_second_kind
  Returns		:double
  Parameter	:double x
  Parameter	:size_t degree
  Brief:
        This is the chebyshev polynomial of the second kind,
        represented in trigonometric terms
  Assumes:
  Consider:
  Note:
**************************************************************/
double Chebyshev_second_kind(double x, size_t n)
{
  double acos_x = std::acos(x);

  ENSURE(NonZero(acos_x) == true);

  return std::sin((n + 1) * acos_x) / acos_x;
}

/*!************************************************************
  FullName	:Chebyshev_first_derivative
  Returns		:double
  Parameter	:double x
  Parameter	:size_t degree
  Brief:
        n \cdot (second_kind)_{n-1}(x)
  Assumes:
  Consider:
  Note:
**************************************************************/
double Chebyshev_first_derivative(double x, size_t n)
{
  if (n == 0) return 0.0; //double check this
  else if (std::abs(x) == 1.0)
  {
    return std::pow(1, n);
  }
  else
  {
    return n * Chebyshev_second_kind(x, n - 1);
  }
}

/*!************************************************************
  FullName	:Chebyshev_second_derivative
  Returns		:double
  Parameter	:double x
  Parameter	:size_t degree
  Brief:
        \frac{n}{x^2 - 1}(n * T_n(x) - x * U_{n-1}(x))
  Assumes:
  Consider: We can't accept 1 or -1, since it will cause a
            division by 0
  Note:
**************************************************************/
double Chebyshev_second_derivative(double x, size_t n)
{
  if (n == 0) return 0.0;
  else if (x == 1.0)
  {
    return (n * n * n * n - n * n) / (3.0);
  }
  else if (x == -1.0)
  {
    return std::pow(-1.0, n) * ((n * n * n * n - n * n) / (3.0));
  }
  else
  {
    double denom = x * x - 1;
    double coeff = n / denom;
    double first_kind = n * Chebyshev_first_kind(x, n);
    double second_kind = x * Chebyshev_second_kind(x, n - 1);

    return coeff * (first_kind - second_kind);
  }
}

/*!************************************************************
  FullName	:Legendre
  Returns		:double
  Parameter	:double x
  Parameter	:size_t degree
  Brief:

  Assumes:
  Consider:
  Note:
**************************************************************/
double Legendre(double x, size_t degree)
{
  if (degree == 0) return 1.0;
  else if (degree == 1) return x;
  else
  {
    size_t n = 2;
    double pn;
    double pn_1 = x;
    double pn_2 = 1.0;

    while (n <= degree)
    {
      pn = ((2 * n - 1) * x * pn_1 - (n - 1) * pn_2) / n;
      pn_2 = pn_1;
      pn_1 = pn;
      ++n;
    }

    return pn;
  }
}

namespace P1
{
  /*!************************************************************
    FullName	:P1::Q1
    Returns		:Roots
    Parameter	:const QuestionInfo & qinfo
    Brief:

    Assumes:
    Consider:
    Note:
  **************************************************************/
  Roots Q1(const QuestionInfo& qinfo)
  {
    const auto& degree = qinfo.degree;
    std::function<Polynomial> chebyshev = [degree](double x)->double
    {
      return std::cos(degree * std::acos(x));
    };

    //we know the range is between [-1,1] inclusive
    Ranges ranges = FindRanges2(qinfo.begin, qinfo.end, qinfo.degree, chebyshev);

    Roots roots;
    for (const auto& elem : ranges)
    {
      double root = Bisection(elem.first, elem.second, qinfo.precision, chebyshev);
      roots.push_back(root);
    }

    return roots;
  }

  /*!************************************************************
    FullName	:P1::Q2
    Returns		:Roots
    Parameter	:const QuestionInfo & qinfo
    Brief:

    Assumes:
    Consider:
    Note:
  **************************************************************/
  Roots Q2(const QuestionInfo& qinfo)
  {
    // P_{n+1}(x) = ((2n+1) x P_n(x) - n P_{n-1}(x)) / (n + 1)qinfo
    const auto& degree = qinfo.degree;
    std::function<Polynomial> legendre = [degree](double x)->double
    {
      if (degree == 0) return 1.0;
      else if (degree == 1) return x;
      else
      {
        size_t n = 2;
        double pn;
        double pn_1 = x;
        double pn_2 = 1.0;

        while (n <= degree)
        {
          pn = ((2 * n - 1) * x * pn_1 - (n - 1) * pn_2) / n;
          pn_2 = pn_1;
          pn_1 = pn;
          ++n;
        }

        return pn;
      }
    };

    //we know the range is between [-1,1] inclusive
    Ranges ranges = FindRanges2(qinfo.begin, qinfo.end, qinfo.degree, legendre);

    Roots roots;
    for (const auto& elem : ranges)
    {
      double root = Bisection(elem.first, elem.second, qinfo.precision, legendre);
      roots.push_back(root);
    }

    return roots;
  }

  /*
    We need to:
      find a positive/negative pair of inputs such that we can perform bracketing


      This one is broken, use FindRanges2 instead
  */
  /*!************************************************************
    FullName	:P1::FindRanges
    Returns		:Ranges
    Parameter	:double begin
    Parameter	:double end
    Parameter	:size_t roots
    Parameter	:std::function<Polynomial> func
    Brief:

    Assumes:
    Consider:
    Note:
  **************************************************************/
  Ranges FindRanges(double begin, double end, size_t roots, std::function<Polynomial> func)
  {
    Ranges valid_ranges; //we return this

    Ranges range1; Ranges* from_range = &range1;
    Ranges range2; Ranges* to_range = &range2;

    double interval = (end - begin) / roots;
    double range_begin = begin;
    double range_end = begin + interval;

    //initialize the first range
    for (size_t i = 0; i < roots; ++i)
    {
      Range range;

      range.first = range_begin;
      range.second = range_end;

      range_begin = range_end;
      range_end += interval;

      from_range->push_back(range);
    }

    /*
      Sift through the ranges.
      if a range validates IVT (it intersects), we toss it into valid_ranges
      otherwise we break it into 2 and toss it into to_range.
    */
    while (valid_ranges.size() < roots)
    {
      for (const auto& elem : *from_range)
      {
        //if IVT succeeds
        if (func(elem.first) * func(elem.second) < 0.0)
        {
          valid_ranges.push_back(elem);
        }
        else
        {
          //otherwise cut range into half and hope for the best
          Range first_half;
          Range second_half;

          double midpoint = (elem.first + elem.second) / 2.0;

          first_half.first = elem.first;
          first_half.second = midpoint;

          second_half.first = midpoint;
          second_half.second = elem.second;

          to_range->push_back(first_half);
          to_range->push_back(second_half);
        }
      }

      //cleanup, swap pointers
      from_range->clear();
      std::swap(from_range, to_range);
    }

    return valid_ranges;
  }

  /*!************************************************************
    FullName	:P1::FindRanges2
    Returns		:Ranges
    Parameter	:double begin
    Parameter	:double end
    Parameter	:size_t roots
    Parameter	:std::function<Polynomial> func
    Brief:

    Assumes:
    Consider:
    Note:
  **************************************************************/
  Ranges FindRanges2(double begin, double end, size_t roots, std::function<Polynomial> func)
  {
    Ranges valid_ranges;

    double interval = (end - begin) / (roots * roots);
    double chaser = begin;
    double leader = chaser + interval;

    while (valid_ranges.size() < roots && leader <= end)
    {
      double chase_func = func(chaser);
      double lead_func = func(leader);

      if (chase_func * lead_func < 0.0)
      {
        Range range;
        range.first = chaser;
        range.second = leader;

        valid_ranges.push_back(range);

        chaser = leader;
      }

      leader += interval;
    }

    assert(valid_ranges.size() == roots);

    return valid_ranges;
  }

  /*!************************************************************
    FullName	:P1::Project1
    Returns		:void
    Brief:

    Assumes:
    Consider:
    Note:
  **************************************************************/
  void Project()
  {
    std::cout.precision(12);

    bool loop = true;
    while (loop == true)
    {
      size_t degree = 0;
      size_t question = 0;

      std::cout << "Assumptions: Range [-1,1], epsilon 0.000000000001" << std::endl;

      do {
        std::cout << "Question[1, 2]: " << std::endl;
        std::cin >> question;
      } while (question != 1 && question != 2);
      do {
        std::cout << "Degree[1, 16]: " << std::endl;
        std::cin >> degree;
      } while (degree <= 0 || degree > 16);

      QuestionInfo qinfo;
      qinfo.degree = degree;
      qinfo.begin = -1.0;
      qinfo.end = 1.0;
      qinfo.precision = 0.000000000001; //10^(-12)

      Roots roots;
      switch (question)
      {
      case 1:
        roots = Q1(qinfo);
        break;
      case 2:
        roots = Q2(qinfo);
        break;

      }

      std::cout << std::endl;
      for (double root : roots)
      {
        std::cout << std::scientific << std::setw(20) << root << std::endl;
      }
      std::cout << std::endl;

      std::string str;
      do
      {
        std::cout << "Again? [y/n]" << std::endl;
        std::cin >> str;
      } while (str != "y" && str != "n");

      if (str == "y")
      {
        loop = true;
      }
      else
      {
        loop = false;
      }

      std::cout << std::endl;
    }
  }
}

namespace P2
{
  /*!************************************************************
    FullName	:P2::Q1
    Returns		:Values
    Parameter	:const QuestionInfo & qinfo
    Brief:
          TODO: Verify that this works
    Assumes:
    Consider:
    Note:
  **************************************************************/
  Values Q1(const QuestionInfo& qinfo)
  {
    //populate the chebyshev polynomial
    std::vector<double> values; values.reserve(qinfo.n);

    //find inputs of chebyshev such that output evaluates to 1
    auto n = qinfo.n;
    auto a_cheb = [n](int k) {  return std::cos((k * PI) / n); };

    for (int i = 0; i <= qinfo.n; ++i)
    {
      values.push_back(a_cheb(i));
    }

    std::reverse(values.begin(), values.end());

    //then find out the lagrange basis
    LagrangeBasisPolynomial lagbasis(values);

    std::vector<double> lagrange;

    for (size_t i = 0; i < values.size(); ++i)
    {
      lagrange.push_back(lagbasis(qinfo.y, i));
    }

    return lagrange;
  }

  /*!************************************************************
    FullName	:P2::Q2
    Returns		:void
    Parameter	:const QuestionInfo & qinfo
    Brief:

    Assumes:
    Consider:
    Note:
  **************************************************************/
  double Q2(const QuestionInfo& qinfo)
  {
    DataPlot plot = {
      {-1.0, Legendre(-1.0, qinfo.n)},
      {-0.5, Legendre(-0.5, qinfo.n)},
      {+0.0, Legendre(+0.0, qinfo.n)},
      {+0.5, Legendre(+0.5, qinfo.n)},
      {+1.0, Legendre(+1.0, qinfo.n)} };

    //since the points are hardcoded and uniformly spaced, we
    //can make some pretty nice assumptions/optimizations with code

    Matrix<3, 3, double> mat3;
    mat3[0][0] = 2.0; mat3[0][1] = 0.5; mat3[0][2] = 0.0;
    mat3[1][0] = 0.5; mat3[1][1] = 2.0; mat3[1][2] = 0.5;
    mat3[2][0] = 0.0; mat3[2][1] = 0.5; mat3[2][2] = 2.0;

    GaussJordan<3, double> gjs; (gjs); //to remove warnings

    Matrix<3, 3, double> inv3 = gjs.FindInverse(mat3);

    Vector<3, double> tosolve;

    for (size_t i = 0, j = 1; i < Vec3::Size; ++i, ++j)
    {
      double yip1 = plot[j + 1].second;
      double yi = plot[j].second;
      double yim1 = plot[j - 1].second;

      double xip1 = plot[j + 1].first;
      double xi = plot[j].first;
      double xim1 = plot[j - 1].first;

      double eval = 6 * ((yip1 - yi) / (xip1 - xi) - (yi - yim1) / (xi - xim1));

      tosolve[i] = eval;
    }

    Vector<3, double> solved = inv3 * tosolve;

    //we (should) now have the second derivatives for points [1,5]
    //determine which interval does it belong to

    if (plot[0].first <= qinfo.y && qinfo.y <= plot[1].first)
    {
      return CubicSpline(qinfo.y,
        0.0, solved[0],
        plot[0].second, plot[1].second,
        plot[0].first, plot[1].first);
    }
    if (plot[1].first <= qinfo.y && qinfo.y <= plot[2].first)
    {
      return CubicSpline(qinfo.y,
        solved[0], solved[1],
        plot[1].second, plot[2].second,
        plot[1].first, plot[2].first);
    }
    if (plot[2].first <= qinfo.y && qinfo.y <= plot[3].first)
    {
      return CubicSpline(qinfo.y,
        solved[1], solved[2],
        plot[2].second, plot[3].second,
        plot[2].first, plot[3].first);
    }
    if (plot[3].first <= qinfo.y && qinfo.y <= plot[4].first)
    {
      return CubicSpline(qinfo.y,
        solved[2], 0.0,
        plot[3].second, plot[4].second,
        plot[3].first, plot[4].first);
    }

    ABORT();

    return double{};
  }

  /*!************************************************************
    FullName	:P2::Project
    Returns		:void
    Brief:

    Assumes:
    Consider:
    Note:
  **************************************************************/
  void Project()
  {
    std::cout.precision(12);

    bool loop = true;
    while (loop == true)
    {
      size_t degree = 0;
      double y = 0.f;
      size_t question = 0;

      std::cout << "Assumptions: 0 < n <= 16, -1 <= y <= 1" << std::endl;

      do {
        std::cout << "Question[1, 2]: " << std::endl;
        std::cin >> question;
      } while (question != 1 && question != 2);
      do {
        std::cout << "n(0, 16]: " << std::endl;
        std::cin >> degree;
      } while (degree <= 0 || degree > 16);
      do {
        std::cout << "y[-1, 1]: " << std::endl;
        std::cin >> y;
      } while (y < -1 || y > 1);

      QuestionInfo qinfo(static_cast<int>(degree), y);


      switch (question)
      {
      case 1:
      {
        std::vector<double> values = Q1(qinfo);

        std::cout << std::endl;
        for (double elem : values)
        {
          std::cout << std::scientific << std::setw(20) << elem << std::endl;
        }
        std::cout << std::endl;
      }
      break;
      case 2:
      {
        double eval = Q2(qinfo);

        std::cout << std::endl;
        std::cout << std::scientific << std::setw(20) << eval << std::endl;
        std::cout << std::endl;
      }
      break;

      }

      std::string str;
      do
      {
        std::cout << "Again? [y/n]" << std::endl;
        std::cin >> str;
      } while (str != "y" && str != "n");

      if (str == "y")
      {
        loop = true;
      }
      else
      {
        loop = false;
      }

      std::cout << std::endl;
    }
  }

  /*!************************************************************

  **************************************************************/
  QuestionInfo::QuestionInfo(int _i, double _y)
    : n(_i), y(_y)
  {
    ENSURE(0 < n && n <= 16);
    ENSURE(-1.0 <= y && y <= 1.0);
  }

}

namespace P3
{
  /*!************************************************************
    Yeah, we should puff this up a bit more
  **************************************************************/
  void Q1(const Q1QuestionInfo& qinfo)
  {
    std::cout << std::endl;
    std::cout << "T'n(y)  = " << Chebyshev_first_derivative(qinfo.y, qinfo.n) << std::endl;
    std::cout << "T''n(y) = " << Chebyshev_second_derivative(qinfo.y, qinfo.n) << std::endl;
    std::cout << std::endl;
  }

  /*!************************************************************

  **************************************************************/
  void Q2(const Q2QuestionInfo& qinfo)
  {
    //load the x
    std::vector<double> x = { -1.0, -0.5, 0.0, 0.5, 1.0 }; //hard coded yeaaaah!
    double h = 0.5;

    //load the y
    std::vector<double> y = x;
    std::transform(y.begin(), y.end(), y.begin(), [qinfo](double x) {return Legendre(x, qinfo.n); });

    //find out the P'n
    std::vector<double> yp = y;
    yp[0] = SecondOrderFirstForwardFF(&yp[0], h);
    for (size_t i = 1; i < yp.size() - 1; ++i)
    {
      yp[i] = FirstOrderFirstCentralFF(&yp[i], h);
    }
    yp[4] = SecondOrderFirstBackwardFF(&yp[4], h);

    //find out the P''n
    std::vector<double> ypp = y;
    ypp[0] = ThirdOrderSecondForwardFF(&ypp[0], h);
    for (size_t i = 1; i < ypp.size() - 1; ++i)
    {
      ypp[i] = FirstOrderSecondCentralFF(&ypp[i], h);
    }
    ypp[4] = ThirdOrderSecondBackwardFF(&ypp[4], h);

    //load the plots
    Plot2D<> first_deriv = { x, yp };
    Plot2D<> second_deriv = { x, ypp };

    double fp = NevilleMethod(first_deriv, qinfo.y);
    double fpp = NevilleMethod(second_deriv, qinfo.y);

    std::cout << std::endl;
    std::cout << "P'n(y)  = " << fp << std::endl;
    std::cout << "P''n(y) = " << fpp << std::endl;
    std::cout << std::endl;
  }

  /*!************************************************************

  **************************************************************/
  void Project()
  {
    std::cout.precision(12);

    bool loop = true;
    while (loop == true)
    {
      size_t degree = 0;
      double y = 0.f;
      size_t question = 0;

      do {
        std::cout << "Question[1, 2]: " << std::endl;
        std::cin >> question;
      } while (question != 1 && question != 2);

      switch (question)
      {
      case 1:
      {
        do {
          std::cout << "n(0, 16]: " << std::endl;
          std::cin >> degree;
        } while (degree <= 0 || degree > 16);
        do {
          std::cout << "y[-1, 1]: " << std::endl;
          std::cin >> y;
        } while (y < -1 || y > 1);

        Q1(Q1QuestionInfo(degree, y));
      }
      break;
      case 2:
      {
        do {
          std::cout << "n[6, 16]: " << std::endl;
          std::cin >> degree;
        } while (degree < 6 || degree > 16);
        do {
          std::cout << "y[-1, 1]: " << std::endl;
          std::cin >> y;
        } while (y < -1 || y > 1);

        Q2(Q2QuestionInfo(degree, y));
      }
      break;
      default:
        std::cout << "Unknown question!" << std::endl;
      }

      std::string str;
      do
      {
        std::cout << "Again? [y/n]" << std::endl;
        std::cin >> str;
      } while (str != "y" && str != "n");

      if (str == "y")
      {
        loop = true;
      }
      else
      {
        loop = false;
      }

      std::cout << std::endl;
    };
  }

  /*!************************************************************

  **************************************************************/
  Q1QuestionInfo::Q1QuestionInfo(size_t degree, double argument)
    : n(degree)
    , y(argument)
  {
    //as specified by the assignment
    EXPECT(0 < n && n <= 16);
    EXPECT(-1 <= y && y <= 1);
  }

  /*!************************************************************

  **************************************************************/
  Q2QuestionInfo::Q2QuestionInfo(size_t degree, double argument)
    : n(degree)
    , y(argument)
  {
    //as specified by the assignment
    EXPECT(6 <= n && n <= 16);
    EXPECT(-1 <= y && y <= 1);
  }

} //namespace P3 //for project 3

/*!************************************************************

**************************************************************/
P4::Q1QuestionInfo::Q1QuestionInfo(size_t degree, const std::vector<double>& y)
  : degree(degree)
  , data(y)
{
  EXPECT(1 <= degree && degree <= 16);
  EXPECT(y.size() == degree);
}

/*!************************************************************

**************************************************************/
void P4::Q1(const Q1QuestionInfo& qinfo)
{

  std::cout.precision(12);

  std::vector<double> chebyshev_x(qinfo.degree);
  for (size_t i = 1; i <= qinfo.degree; ++i)
  {
    double arg = ((2.0 * i - 1) / (2.0 * qinfo.degree));
    double node = std::cos(arg * PI);

    chebyshev_x[i - 1] = node;
  }

  //we got the roots, now assemble into the plot2d
  Plot2D<double> plot(chebyshev_x, qinfo.data);

  double quadrature = plot.sum_each_xy_pair([](double x, double y) {
    double denom = std::sqrt(1 - x * x);
    return y / denom;
  });

  quadrature *= PI / qinfo.degree;

  std::cout << "Evaluation: " << quadrature << std::endl;

}

/*!************************************************************

**************************************************************/
P4::Q2QuestionInfo::Q2QuestionInfo(size_t d, size_t m)
  : degree(d)
  , numnodes(m)
{
  EXPECT(1 <= degree && degree <= 16);
  EXPECT(degree <= numnodes && numnodes <= 2 * degree);
}

/*!************************************************************

**************************************************************/
void P4::Q2(const Q2QuestionInfo& qinfo)
{

  std::cout.precision(12);

  //populate the x values
  std::vector<double> legendre_x(qinfo.numnodes);

  const double delta = 2.0 / (qinfo.numnodes - 1);
  double sum = 0.0;
  for (size_t i = 0; i < qinfo.numnodes; ++i)
  {
    double value = -1.0 + sum;
    legendre_x[i] = value;
    sum += delta;
  }

  //populate the y values
  std::vector<double> legendre_y(qinfo.numnodes);
  for (size_t i = 0; i < qinfo.numnodes; ++i)
  {
    double result = Legendre(legendre_x[i], qinfo.degree);
    legendre_y[i] = result;
  }

  //evaluate
  Plot2D<double> data(legendre_x, legendre_y);
  //double finish = RecursiveTrapezoidMethod(data, -1.0, 1.0, /*arbitrary*/5);
  double finish = CompositeTrapezoidMethod(-1.0, 1.0, legendre_y);

  std::cout << finish << std::endl;

}

/*!************************************************************

**************************************************************/
void P4::Project()
{
  std::cout.precision(12);

  bool loop = true;
  while (loop == true)
  {
    size_t degree = 0;
    double y = 0.f;
    size_t question = 0;

    do {
      std::cout << "Question[1, 2]: " << std::endl;
      std::cin >> question;
    } while (question != 1 && question != 2);

    switch (question)
    {
    case 1:
    {
      do {
        std::cout << "Degree[1, 16]: " << std::endl;
        std::cin >> degree;
      } while (degree < 1 || 16 < degree);

      std::vector<double> y_values(degree);

      for (size_t i = 0; i < degree; ++i)
      {
        std::cout << "Node[" << i + 1 << "/" << degree << "]: ";
        std::cin >> y_values[i];
      }

      Q1(Q1QuestionInfo(degree, y_values));
    }
    break;
    case 2:
    {
      size_t m = 0;

      do {
        std::cout << "Degree[1, 16]: " << std::endl;
        std::cin >> degree;
      } while (degree < 1 || 16 < degree);;
      do {
        std::cout << "Nodes[" << degree << "," << 2 * degree << "]: " << std::endl;
        std::cin >> m;
      } while (m < degree || degree * 2 < m);

      Q2(Q2QuestionInfo(degree, m));
    }
    break;
    default:
      std::cout << "Unknown question!" << std::endl;
    }

    std::string str;
    do
    {
      std::cout << "Again? [y/n]" << std::endl;
      std::cin >> str;
    } while (str != "y" && str != "n");

    if (str == "y")
    {
      loop = true;
    }
    else
    {
      loop = false;
    }

    std::cout << std::endl;
  };
}


