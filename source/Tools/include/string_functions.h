/*
 * string_functions.h
 *
 *  Created on: 3. 11. 2016
 */

#ifndef INCLUDE_STRING_FUNCTIONS_H_
#define INCLUDE_STRING_FUNCTIONS_H_

//Include
#include <algorithm>
#include <sstream>
#include <string>

#include "StringConverter.h"

/**
 * Various global string helper functions, originally saved in global_funcs file
 */

//! Split string to individual sub-strings using a given separator
std::vector<std::string> split(const std::string& str, const std::string& seps = " \t\n", bool keepEmpty = false);

//! Split string to individual sub-strings using a given separator & then transform the sub-strings to required type
template<class T> std::vector<T> split(const std::string& str, const std::string& seps = " \t\n", bool keepEmpty = false) {

  auto vs = split(str, seps, keepEmpty);
  std::vector<T> vt(vs.size());
  std::transform(vs.begin(), vs.end(), vt.begin(), [](std::string s) { return str2any<T>(s); });
  return vt;
}

//! TODO: Document
template<class T, class I> std::string join(I begin, I end, const std::string& sep) {
  std::stringstream ss;
  std::copy(begin, end, std::ostream_iterator<T>(ss, sep.c_str()));
  return ss.str().substr(0, ss.str().length()-1);
}

//! TODO: Document
template<template<class> class T, class U> std::string join(const T<U>& vec, const std::string& sep) { return join<U>(vec.begin(), vec.end(), sep); }

std::string ltrim(std::string str);
std::string rtrim(std::string str);
std::string trim(std::string str);
std::string lctrim(std::string str, const std::string& chars);
std::string rctrim(std::string str, const std::string& chars);
std::string ctrim(std::string str , const std::string& chars);

#endif /* INCLUDE_STRING_FUNCTIONS_H_ */
