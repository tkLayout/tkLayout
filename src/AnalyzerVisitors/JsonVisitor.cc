// AnalyzerVisitors/JsonVisitor.cc

#include <boost/json.hpp>
#include "Tracker.hh"
#include "Barrel.hh"
#include "Endcap.hh"
#include "Layer.hh"
#include "Ring.hh"
#include "DetectorModule.hh"
#include "GeometricModule.hh"
#include "Utilities/PropertyJsonHelpers.hh"
#include "AnalyzerVisitors/JsonVisitor.hh"

namespace json = boost::json;

boost::json::object JsonVisitor::merge_json(const boost::json::object& a, const boost::json::object& b) {
    boost::json::object result = a;
    for (auto& kv : b)
        result.insert_or_assign(kv.key(), kv.value());
    return result;
}

json::object JsonVisitor::build(const Tracker *tracker_ptr, const Tracker *pixel_ptr)
{
    json::object outer, inner;
    
    if (tracker_ptr != nullptr) outer = visit_tracker(*tracker_ptr);
    if (pixel_ptr != nullptr) inner = visit_tracker(*pixel_ptr);
    
    merge_json(outer, inner);
    return outer;

    // The following code is an example on how to use section names in the JSON output.
    // auto add_sections = [&sections](const Tracker& t) {
    //     for (const auto& b : t.barrels()) {
    //         json::object section;
    //         section["type"] = "barrel";
    //         json::array layer_names;
    //         for (const auto& l : b.layers())
    //             layer_names.emplace_back("L" + std::to_string(l.layerNumber()));
    //         section["layers"] = std::move(layer_names);
    //         sections[b.myid()] = std::move(section);
    //     }
    // 
    //     for (const auto& e : t.endcaps()) {
    //         json::object section;
    //         section["type"] = "endcap";
    //         json::array disk_names;
    //         for (const auto& d : e.disks())
    //             disk_names.emplace_back("D" + std::to_string(d.diskNumber()));
    //         section["disks"] = std::move(disk_names);
    //         sections[e.myid()] = std::move(section);
    //     }
    // };

}

json::object JsonVisitor::visit_tracker(const Tracker &t)
{
    json::object jo;
    json::array barrels;
    for (const auto &b : t.barrels())
        barrels.push_back(visit_barrel(b));
    jo["barrels"] = std::move(barrels);

    json::array endcaps;
    for (const auto &e : t.endcaps())
        endcaps.push_back(visit_endcap(e));
    jo["endcaps"] = std::move(endcaps);

    return jo;
}

json::object JsonVisitor::visit_barrel(const Barrel &b)
{
    json::object jo;
    // jo["properties"] = dump_properties(b);

    json::array layers;
    for (const auto &l : b.layers())
        layers.push_back(visit_layer(l));
    jo["layers"] = std::move(layers);
    return jo;
}

json::object JsonVisitor::visit_endcap(const Endcap &e)
{
    json::object jo;
    // jo["properties"] = dump_properties(e);

    json::array disks;
    for (const auto &d : e.disks())
        disks.push_back(visit_disk(d));
    jo["disks"] = std::move(disks);
    return jo;
}

json::object JsonVisitor::visit_disk(const Disk &d)
{
    json::object jo;
    json::array rings;

    for (const auto &ring : d.rings())
    {
        rings.push_back(visit_ring(ring));
    }
    jo["name"] = "ANAME"; // Placeholder for disk name; replace with actual name
    jo["rings"] = std::move(rings);
    return jo;
}

json::object JsonVisitor::visit_layer(const Layer &l)
{
    json::object jo;
    // jo["properties"] = dump_properties(l);

    json::array rods;
    for (const auto &r : l.rods())
        rods.push_back(visit_rodpairs(r));
    jo["name"] = "ANAME"; // Placeholder for disk name; replace with actual name
    jo["rods"] = std::move(rods);
    return jo;
}

json::object JsonVisitor::visit_rodpairs(const RodPair& r) {
    json::object jo;
    // jo["properties"] = dump_properties(r);

    json::array modules;
    for (const auto& mod : r.modules().first)
        modules.push_back(visit_module(mod));
    for (const auto& mod : r.modules().second)
        modules.push_back(visit_module(mod));

    jo["modules"] = std::move(modules);
    return jo;
}

json::object JsonVisitor::visit_ring(const Ring &r)
{
    json::object jo;
    // jo["properties"] = dump_properties(r);

    json::array mods;
    for (const auto &m : r.modules())
        mods.push_back(visit_module(m));
    jo["modules"] = std::move(mods);
    return jo;
}

json::object JsonVisitor::visit_module(const DetectorModule &m)
{
    json::object jo;
    // jo["properties"] = dump_properties(m);

    const auto &c = m.center();
    json::object pos;
    pos["x"] = c.X();
    pos["y"] = c.Y();
    pos["z"] = c.Z();
    jo["position"] = std::move(pos);

    return jo;
}
