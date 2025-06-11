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
    json::value build(const Tracker &tk)
    {
        json::object jo;
        jo["properties"] = dump_properties(tk);

        json::array barrels;
        for (const auto &b : tk.barrels())
            barrels.push_back(visit_barrel(b));

        json::array endcaps;
        for (const auto &e : tk.endcaps())
            endcaps.push_back(visit_endcap(e));

        jo["barrels"] = std::move(barrels);
        jo["endcaps"] = std::move(endcaps);
        return jo;
    }

private:
    json::object visit_barrel(const Barrel &b)
    {
        json::object jo;
        jo["properties"] = dump_properties(b);

        json::array layers;
        for (const auto &l : b.layers())
            layers.push_back(visit_layer(l));
        jo["layers"] = std::move(layers);
        return jo;
    }

    json::object visit_endcap(const Endcap &e)
    {
        json::object jo;
        jo["properties"] = dump_properties(e);

        json::array disks;
        for (const auto &d : e.disks())
            disks.push_back(visit_disk(d));
        jo["disks"] = std::move(disks);
        return jo;
    }

    json::object visit_disk(const Disk &d)
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

    json::object visit_layer(const Layer &l)
    {
        json::object jo;
        jo["properties"] = dump_properties(l);

        json::array rods;
        for (const auto &r : l.rods())
            rods.push_back(visit_rodpairs(r));
        jo["rods"] = std::move(rods);
        return jo;
    }

    json::object visit_rodpairs(const RodPair& r) {
        json::object jo;
        jo["properties"] = dump_properties(r);

        json::array modules;
        for (const auto& mod : r.modules().first)
            modules.push_back(visit_module(mod));
        for (const auto& mod : r.modules().second)
            modules.push_back(visit_module(mod));

        jo["modules"] = std::move(modules);
        return jo;
    }


    json::object visit_ring(const Ring &r)
    {
        json::object jo;
        jo["properties"] = dump_properties(r);

        json::array mods;
        for (const auto &m : r.modules())
            mods.push_back(visit_module(m));
        jo["modules"] = std::move(mods);
        return jo;
    }

    json::object visit_module(const DetectorModule &m)
    {
        json::object jo;
        jo["properties"] = dump_properties(m);

        const auto &c = m.center();
        json::object pos;
        pos["x"] = c.X();
        pos["y"] = c.Y();
        pos["z"] = c.Z();
        jo["position"] = std::move(pos);

        return jo;
    }
};
