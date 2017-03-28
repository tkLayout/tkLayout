#include <global_funcs.hh>

template<typename T> const std::vector<std::string> EnumTraits<T>::data = {};

template<> std::string StringConverter<_NOT_STRING_ENUM>::str2any<std::string>(const std::string& from) { return from; };

template<> bool StringConverter<_NOT_STRING_ENUM>::str2any<bool>(const std::string& from) {
  static std::map<std::string, bool> boolstr = { {"true", true}, {"TRUE", true}, {"True", true}, {"T", true}, {"1", true},
                                                 {"false", false}, {"FALSE", false}, {"False", false}, {"F", false}, {"0", false} };
  return boolstr.at(from); 
}


std::vector<std::string> split(const std::string& str, const std::string& seps, bool keepEmpty) {
  std::vector<std::string> tokens;
  std::string token;
  size_t start = 0, end = 0;
  while ((end = str.find_first_of(seps, start)) != std::string::npos) {
    token = str.substr(start, end - start);
    if (!token.empty() || keepEmpty) tokens.push_back(token);
    start = end + 1;
  }
  token = str.substr(start);
  if (!token.empty() || keepEmpty) tokens.push_back(token);
  return tokens;
}


std::string ltrim(std::string str) {
  return str.erase(0, str.find_first_not_of(" \t\n"));
}

std::string rtrim(std::string str) {
  return str.erase(str.find_last_not_of(" \t\n")+1);
}

std::string trim(std::string str) {
  return ltrim(rtrim(str));
}

std::string lctrim(std::string str, const std::string& chars) {
  return str.erase(0, str.find_first_not_of(chars));
}

std::string rctrim(std::string str, const std::string& chars) {
  return str.erase(str.find_last_not_of(chars)+1);
}

std::string ctrim(std::string str, const std::string& chars) {
  return lctrim(rctrim(str, chars), chars);
}
