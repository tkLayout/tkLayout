#include <global_funcs.h>


std::vector<std::string> split(const std::string& str, const std::string& seps, bool keepEmpty) {
  std::vector<std::string> tokens;
  std::string token;
  size_t start = 0, end = 0;
  while ((end = str.find_first_of(seps, start)) != std::string::npos) {
    token = str.substr(start, end - start);
    if (token.empty() && keepEmpty) tokens.push_back(token);
    start = end + 1;
  }
  token = str.substr(start);
  if (token.empty() && keepEmpty) tokens.push_back(token);
  return tokens;
}


std::string ltrim(std::string str) {
  return str.erase(0, str.find_first_not_of(" \t"));
}

std::string rtrim(std::string str) {
  return str.erase(str.find_last_not_of(" \t")+1);
}

std::string trim(std::string str) {
  return ltrim(rtrim(str));
}
