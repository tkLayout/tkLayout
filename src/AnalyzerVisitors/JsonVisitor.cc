// AnalyzerVisitors/JsonVisitor.cc

#include <boost/json.hpp>
#include "Tracker.hh"
#include "Barrel.hh"
#include "Endcap.hh"
#include "Layer.hh"
#include "Ring.hh"
#include "DetectorModule.hh"
#include "GeometricModule.hh"
#include "AnalyzerVisitors/JsonVisitor.hh"

namespace json = boost::json;

json::object JsonVisitor::build(const Tracker *t, const Tracker *p)
{
    json::array trackers;
    if (t) trackers.push_back(visit_tracker(*t));
    if (p) trackers.push_back(visit_tracker(*p));

    return json::object{
        {"Trackers", std::move(trackers)}
    };
}

json::object JsonVisitor::visit_tracker(const Tracker &t)
{
    json::object jo;
    jo["type"] = "Tracker";
    jo["myid"] = t.myid();
    
    json::array barrels;
    for (const auto &b : t.barrels())
        barrels.push_back(visit_barrel(b));
    jo["Barrels"] = std::move(barrels);
    
    json::array endcaps;
    for (const auto &e : t.endcaps())
        endcaps.push_back(visit_endcap(e));
    jo["Endcaps"] = std::move(endcaps);

    return jo;
}

json::object JsonVisitor::visit_barrel(const Barrel &b)
{
    json::object jo;
    jo["type"] = "Barrel";
    jo["myid"] = b.myid();
    json::array layers;
    for (const auto &l : b.layers())
        layers.push_back(visit_layer(l));
    jo["Layers"] = std::move(layers);
    return jo;
}

json::object JsonVisitor::visit_endcap(const Endcap &e)
{
    json::object jo;
    jo["type"] = "Endcap";
    jo["myid"] = e.myid();
    // Add disks to the endcap
    json::array disks;
    for (const auto &d : e.disks())
        disks.push_back(visit_disk(d));
    jo["Disks"] = std::move(disks);
    return jo;
}

json::object JsonVisitor::visit_disk(const Disk &d)
{
    json::object jo;
    jo["type"] = "Disk";
    jo["myid"] = d.myid(); // myid() returns the disk number
    json::array rings;

    for (const auto &ring : d.rings())
        rings.push_back(visit_ring(ring));
    jo["Rings"] = std::move(rings);
    return jo;
}

json::object JsonVisitor::visit_layer(const Layer &l)
{
    json::object jo;
    jo["type"] = "Layer";
    jo["myid"] = l.myid();
    json::array ladderPairs;
    for (const auto &r : l.rods())
        ladderPairs.push_back(visit_rodpairs(r));
    jo["ladderPairs"] = std::move(ladderPairs);
    return jo;
}

json::object JsonVisitor::visit_rodpairs(const RodPair& r) {
    json::object jo;
    jo["type"] = "LadderPair";
    jo["myid"] = r.myid();
    json::array zPlusModules;
    json::array zMinusModules;
    for (const auto& mod : r.modules().first)
        zPlusModules.push_back(visit_module(mod));
    for (const auto& mod : r.modules().second)
        zMinusModules.push_back(visit_module(mod));

    jo["zPlusModules"] = std::move(zPlusModules);
    jo["zMinusModules"] = std::move(zMinusModules);
    return jo;
}

json::object JsonVisitor::visit_ring(const Ring &r)
{
    json::object jo;
    jo["type"] = "Ring";
    jo["myid"] = r.myid();
    json::array mods;
    for (const auto &m : r.modules())
        mods.push_back(visit_module(m));
    jo["Modules"] = std::move(mods);
    return jo;
}

json::object JsonVisitor::visit_module(const DetectorModule &m)
{
    json::object jo;
    jo["type"] = "Module";
    jo["myid"] = m.myid();
    jo["subType"] = m.moduleSubType();
    jo["detId"] = m.myDetId();
    
    jo["moduleType"] = m.moduleType();
    jo["dsDistance"] = m.dsDistance();

    const auto &c = m.center();
    json::object pos;
    pos["x"] = c.X();
    pos["y"] = c.Y();
    pos["z"] = c.Z();
    jo["position"] = std::move(pos);

    return jo;
}
