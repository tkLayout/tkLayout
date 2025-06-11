// PropertyJsonHelpers.hh
#pragma once

#include <boost/json.hpp>
#include <string>
#include <cstdint>
#include "Property.hh"

namespace json = boost::json;

//------------------------------------------------------------------------------
// 1. Small helpers: convert scalars to boost::json::value
//------------------------------------------------------------------------------

template<typename T>
inline json::value to_json_value(const T& v) { return v; }

template<>
inline json::value to_json_value<std::string>(const std::string& v) { return json::value(v); }

//------------------------------------------------------------------------------
// 2. Dump all run‑time properties from a Parsable‑derived object
//------------------------------------------------------------------------------

template<typename T>
struct has_propertyMap
{
private:
    template<typename U>
    static auto test(int) -> decltype(std::declval<U>().propertyMap(), std::true_type());
    template<typename>
    static std::false_type test(...);
public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

template<typename T>
inline std::enable_if_t<has_propertyMap<T>::value, json::object>
dump_properties(const T& objConst)
{
    json::object out;
    // Dump properties from objConst.propertyMap() if available
    const auto& propMap = objConst.propertyMap(); // string → PropertyBaseBase*
    for (const auto& pair : propMap) {
        const auto& name = pair.first;
        const auto* base = pair.second;
        if (auto p = dynamic_cast<const PropertyBase<int>*>(base))            out[name] = to_json_value((*p)());
        else if (auto p = dynamic_cast<const PropertyBase<double>*>(base))    out[name] = to_json_value((*p)());
        else if (auto p = dynamic_cast<const PropertyBase<std::string>*>(base)) out[name] = to_json_value((*p)());
        else if (auto p = dynamic_cast<const PropertyBase<bool>*>(base))      out[name] = to_json_value((*p)());
        else if (auto p = dynamic_cast<const PropertyBase<std::vector<int>>*>(base)) out[name] = to_json_value((*p)());
        else if (auto p = dynamic_cast<const PropertyBase<std::vector<double>>*>(base)) out[name] = to_json_value((*p)());
        else if (auto p = dynamic_cast<const PropertyBase<std::vector<std::string>>*>(base)) out[name] = to_json_value((*p)());
    }
    return out;
}

// Overload for types without propertyMap()
template<typename T>
inline std::enable_if_t<!has_propertyMap<T>::value, json::object>
dump_properties(const T&)
{
    return json::object{};
}
