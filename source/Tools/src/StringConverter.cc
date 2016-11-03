/*
 * StringConverter.cc
 *
 *  Created on: 3. 11. 2016
 */
#include "StringConverter.h"

//
// Data container for enums
//
template<typename T> const std::vector<std::string> EnumTraits<T>::data = {};

//
// String to any, if any is string
//
template<> std::string StringConverter<_NOT_STRING_ENUM>::str2any<std::string>(const std::string& from) { return from; };

//
// String to any, if any is boolean
//
template<> bool StringConverter<_NOT_STRING_ENUM>::str2any<bool>(const std::string& from) {
  static std::map<std::string, bool> boolstr = { {"true", true}, {"TRUE", true}, {"True", true}, {"T", true}, {"1", true},
                                                 {"false", false}, {"FALSE", false}, {"False", false}, {"F", false}, {"0", false} };
  return boolstr.at(from);
}
