#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>

#include <global_funcs.hh>

template<typename T> const std::vector<std::string> EnumTraits<T>::data = {};

template<> std::string StringConverter<_NOT_STRING_ENUM>::str2any<std::string>(const std::string& from) { return from; };

/**
 * @brief Converts a string representation of a boolean to a bool value
 * 
 * This specialization provides a LUT in the form of a static std::map for boolean configurations.
 * 
 * Supported formats:
 * - Lowercase: "true", "false"
 * - Uppercase: "TRUE", "FALSE", "T", "F"
 * - Titlecase: "True", "False"
 * - Numeric:   "1", "0"
 * 
 * @param from The string to be converted.
 * @return The boolean representation of the string.
 * @throws std::out_of_range If the input string is not found in the lookup map.
 */
template<> bool StringConverter<_NOT_STRING_ENUM>::str2any<bool>(const std::string& from) {
  static std::map<std::string, bool> boolstr = { {"true", true}, {"TRUE", true}, {"True", true}, {"T", true}, {"1", true},
                                                 {"false", false}, {"FALSE", false}, {"False", false}, {"F", false}, {"0", false} };
  // Throws an exception for invalid inputs (e.g., "yes", "on")
  return boolstr.at(from); 
}

std::vector<std::string> split(const std::string& str, 
                               const std::string& seps /* = STD_WHITESPACE */, 
                               bool keepEmpty /*= false*/) {
    std::vector<std::string> tokens;
    // boost::is_any_of creates a predicate from the separator string
    // boost::token_compress_on handles the "keepEmpty" logic
    boost::split(tokens, str, boost::is_any_of(seps), keepEmpty ? boost::token_compress_off : boost::token_compress_on);
    return tokens;
}
