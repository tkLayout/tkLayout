
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Property.hh"
#include "global_funcs.hh"

std::function<int()> noDefault() {
  return []() {
    throw std::logic_error("Tried to get value from an unset property");
    return 0;
  };
}
std::function<bool()> cacheIf(const bool &flag) {
  return [&flag]() { return flag; };
}

std::set<string> PropertyObject::globalMatchedProperties_;
std::set<string> PropertyObject::globalUnmatchedProperties_;
