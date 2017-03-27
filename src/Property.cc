
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <sstream>

#include "global_funcs.hh"
#include "Property.hh"



std::function<int()> noDefault() { return [](){ throw std::logic_error("Tried to get value from an unset property"); return 0; }; }
std::function<bool()> cacheIf(const bool& flag) { return [&flag]() { return flag; }; }

std::set<string> PropertyObject::globalMatchedProperties_;
std::set<string> PropertyObject::globalUnmatchedProperties_;
