#include <boost/json.hpp>
namespace json = boost::json;

template<typename T>
json::value to_json_value(const T& v)            { return v; }              // numbers, bool
template<>
json::value to_json_value(const std::string& v)  { return json::value(v); }
