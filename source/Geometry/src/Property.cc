
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Property.h"

std::function<int()> noDefault() { return [](){ throw std::logic_error("Tried to get value from an unset property"); return 0; }; }
std::function<bool()> cacheIf(const bool& flag) { return [&flag]() { return flag; }; }

std::set<std::string> PropertyObject::globalMatchedProperties_;
std::set<std::string> PropertyObject::globalUnmatchedProperties_;
