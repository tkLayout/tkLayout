/*
 * math_functions.h
 *
 *  Created on: 3. 11. 2016
 */
#ifndef INCLUDE_MATH_FUNCTIONS_H_
#define INCLUDE_MATH_FUNCTIONS_H_

#include <math.h>
#include <utility>

/**
 * Various global math helper functions, originally saved in global_funcs file
 */

//! Signum functions
template<typename ArgType> inline int signum(const ArgType& x) {
  return (x > ArgType(0)) - (x < ArgType(0));
}

//! Round function with given precision
template<int Precision, typename ArgRetType> inline ArgRetType roundprec(const ArgRetType& x) {
  static const float p = pow(10., Precision);
  return floor(x * p + 0.5) / p;
}

//! TODO: Document
template<int Magnification, typename ArgType> inline int mapint(const ArgType& x) { // magnification is the number of decimal digits to preserve
  static const float p = pow(10., Magnification);
  return floor(x * p + 0.5);
}

//! TODO: Document
template<int Miniaturization, typename RetType> inline RetType unmapint(int x) { // miniaturization is the number of decimal digits to restore (from those that were preserved)
  static const float p = pow(10., Miniaturization);
  return (RetType)x / p;
}

template<class I> struct RangePair : std::pair<I, I> {
  RangePair(const std::pair<I, I>& p) : std::pair<I, I>(p) {}
  I begin() const { return this->first; }
  I end()   const { return this->second; }
};

template<class I> inline RangePair<I> pair2range(const std::pair<I, I>& p) { return RangePair<I>(p); }

//! Get max of 2 numbers
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

//! Get min of 2 numbers
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

//! Get max value from a given vector of values
//! Example: double maxZ = maxget(vec.begin(), vec.end(), [](Module* m) { return m->maxZ(); }); // gets maxZ from a vector of modules
template<class I, class UnaryOperation> inline auto maxget(I begin, I end, UnaryOperation op) -> decltype(op(*begin)) {
  auto max = op(*begin++);
  for (auto it = begin; it != end; ++it) max = MAX(max, op(*it));
  return max;
}

//! Get min value from a given vector of values
template<class I, class UnaryOperation> inline auto minget(I begin, I end, UnaryOperation op) -> decltype(op(*begin)) {
  auto min = op(*begin++);
  for (auto it = begin; it != end; ++it) min = MIN(min, op(*it));
  return min;
}

// TODO: Document
template<class I, class MemFn> inline auto maxget2(I begin, I end, MemFn fn) -> typename std::remove_reference<decltype((*begin.*fn)())>::type {
  auto max = (*begin.*fn)();
  for (auto it = begin+1; it != end; ++it) max = MAX(max, (*it.*fn)());
  return max;
}

// TODO: Document
template<class I, class MemFn> inline auto minget2(I begin, I end, MemFn fn) -> typename std::remove_reference<decltype((*begin.*fn)())>::type {
  auto min = (*begin.*fn)();
  for (auto it = begin+1; it != end; ++it) min = MIN(min, (*it.*fn)());
  return min;
}

#endif /* INCLUDE_MATH_FUNCTIONS_H_ */
