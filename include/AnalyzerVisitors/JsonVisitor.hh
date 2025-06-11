// JsonVisitor.hh
#pragma once

#include <boost/json.hpp>
#include "Tracker.hh"
#include "Barrel.hh"
#include "Endcap.hh"
#include "Layer.hh"
#include "Ring.hh"
#include "DetectorModule.hh"
#include "GeometricModule.hh"
#include "Utilities/PropertyJsonHelpers.hh"

namespace json = boost::json;

class JsonVisitor
{
public:
    json::object build(const Tracker &tracker);
private:
    json::object visit_tracker(const Tracker &t);
    json::object visit_barrel(const Barrel &b);
    json::object visit_endcap(const Endcap &e);
    json::object visit_disk(const Disk &d);
    json::object visit_layer(const Layer &l);
    json::object visit_rodpairs(const RodPair& r);
    json::object visit_ring(const Ring &r);
    json::object visit_module(const DetectorModule &m);
    // json::object visit_geometric_module(const GeometricModule &m);
    // json::object visit_support(const SupportStructure &s);
    boost::json::object merge_json(const boost::json::object& a, const boost::json::object& b);
};
