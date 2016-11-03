/*
 * StringConverter.h
 *
 *  Created on: 3. 11. 2016
 */
#ifndef INCLUDE_STRINGCONVERTER_H_
#define INCLUDE_STRINGCONVERTER_H_

// Include
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <iterator>
#include <type_traits>

// Define constructs for treatment of enum types
template<typename T> struct EnumTraits {
  static const std::vector<std::string> data;
};

#define define_enum_strings(E) template<> const std::vector<std::string> EnumTraits<E>::data
#define _STRING_ENUM     true
#define _NOT_STRING_ENUM false

//! String converter in 2 variants: for enum templates and non-enum templates (int, double, ...)
template<bool> class StringConverter {};

/**
 * @class StringConverter<false>
 * @brief String converter class for templates that are not of enum type. Class provides template2str or
 * str2template functionality. In the latter case, the template is output with required precision.
 */
template<> class StringConverter<_NOT_STRING_ENUM> {

 public:

  //! Any to string conversion with given precision at the output if template number
  template<typename ArgType> static std::string any2str(const ArgType& from, int precision = -1) {
    std::stringstream to("");
    if (precision > -1) {
      to.precision(precision);
      to.setf(std::ios::fixed, std::ios::floatfield);
    }
    to << from;
    return to.str();
  }

  //! String to string conversion (null conversion)
  static std::string any2str(const std::string& from) { return from; }

  //! Bool to string conversion
  static std::string any2str(bool from)               { return from == true ? "true" : "false"; }

  //! String to anything convesrion
  template<typename RetType> static RetType str2any(const std::string& from) {
    std::stringstream ssfrom(from);
    RetType to;
    ssfrom >> to;
    return to;
  }

 private:

  //! Private constructor, this class is meant to be called via a global function
  StringConverter();

}; // Class

/**
 * @class StringConverter<false>
 * @brief String converter class for templates that are of enum type. Class provides template2str or
 * str2template functionality. In the latter case, the template is output with required precision.
 */
template<> class StringConverter<_STRING_ENUM> {

public:

  //! Any to string conversion with given precision at the output if template number
  template<typename T> static std::string any2str(const T& from, int) {
    return EnumTraits<T>::data[static_cast<typename std::underlying_type<T>::type>(from)];
  }

  //! String to anything convesrion
  template<typename T> static T str2any(const std::string& from) {

    static auto begin = std::begin(EnumTraits<T>::data);
    static auto end   = std::end(EnumTraits<T>::data); // + EnumTraits<T>::data_size;
    return static_cast<T>(std::distance(begin, std::find(begin, end, from)));
    // TODO : throw exception if string not found
  }

private:

  //! Private constructor, this class is meant to be called via a global function
  StringConverter();
}; // Class

//! Global function converting any template to string with given precision (if template is a number)
template<typename T> std::string any2str(const T& from, int precision = -1) {
  return StringConverter<std::is_enum<T>::value>::any2str(from, precision);
};

//! Global function converting any string to template
template<typename T> T str2any(const std::string& from) {
  return StringConverter<std::is_enum<T>::value>::template str2any<T>(from);
};

#endif /* INCLUDE_STRINGCONVERTER_H_ */
