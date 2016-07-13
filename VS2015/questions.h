#ifndef questions_h__
#define questions_h__

#include "rootfinding_func.h"
#include "guarantees.h"
#include "plot.h"

namespace P1
{
  struct QuestionInfo
  {
    double begin;
    double end;
    double precision;
    size_t degree;
  };

  Roots Q1(const QuestionInfo& qinfo);
  Roots Q2(const QuestionInfo& qinfo);
  void Project();

  Ranges FindRanges(double begin, double end, size_t roots, std::function<Polynomial> func);
  Ranges FindRanges2(double begin, double end, size_t roots, std::function<Polynomial> func);
}

namespace P2
{
  struct QuestionInfo
  {
    QuestionInfo(int i, double y);

    int n;
    double y;
  };

  Values Q1(const QuestionInfo& qinfo);

  void Project();
}

namespace P3
{

  struct Q1QuestionInfo
  {
    Q1QuestionInfo(size_t degree, double argument);
    size_t n;
    double y;
  };

  struct Q2QuestionInfo
  {
    Q2QuestionInfo(size_t degree, double argument);
    size_t n;
    double y;
  };

  void Q1(const Q1QuestionInfo& qinfo);
  void Q2(const Q2QuestionInfo& qinfo);

  void Project();

} //namespace P3

namespace P4
{
  struct Q1QuestionInfo
  {
    Q1QuestionInfo(size_t degree, const std::vector<double>& y);
    size_t degree;
    std::vector<double> data;
  };

  struct Q2QuestionInfo
  {
    Q2QuestionInfo(size_t degree, size_t numnodes);
    size_t degree;
    size_t numnodes;
  };

  void Q1(const Q1QuestionInfo& info);
  void Q2(const Q2QuestionInfo& info);

  void Project();
}

#endif