/**
 * @file Extractor.cc
 * @brief
 */

#include <Extractor.h>
namespace insur {
    //public
    void Extractor::analyse(MaterialTable& mt, MaterialBudget& mb, std::vector<Element>& e,
            std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
            std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::vector<SpecParInfo>& t) {
        Tracker& tr = mb.getTracker();
        InactiveSurfaces& is = mb.getInactiveSurfaces();
        std::vector<std::vector<ModuleCap> >& bc = mb.getBarrelModuleCaps();
        std::vector<std::vector<ModuleCap> >& ec = mb.getEndcapModuleCaps();
        //int layer;
        e.clear();
        c.clear();
        l.clear();
        s.clear();
        p.clear();
        a.clear();
        t.clear();
        // inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        AlgoInfo alg;
        SpecParInfo spec;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        pos.trans.dz = 0.0;
        pos.rot.name = "";
        pos.rot.phix = 0.0;
        pos.rot.phiy = 0.0;
        pos.rot.phiz = 0.0;
        pos.rot.thetax = 0.0;
        pos.rot.thetay = 0.0;
        pos.rot.thetaz = 0.0;
        // define top-level volume and hierarchy root
        shape.name_tag = xml_tracker;
        shape.type = tb;
        shape.dx = 0.0;
        shape.dxx = 0.0;
        shape.dy = 0.0;
        shape.dz = max_length;
        shape.rmin = inner_radius;
        shape.rmax = outer_radius + volume_width;
        s.push_back(shape);
        logic.name_tag = shape.name_tag;
        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
        logic.material_tag = xml_material_air;
        l.push_back(logic);
        spec.name = xml_full_tracker + xml_par_tail;
        spec.parameter.first = xml_tkddd_structure;
        spec.parameter.second = xml_full_tracker;
        spec.partselectors.push_back(xml_tracker);
        t.push_back(spec);
        spec.partselectors.clear();
        // define top-level barrel volume container (polycone)
        shape.type = pc;
        shape.name_tag = xml_tob;
        // create polycone info
        analyseBarrelContainer(tr, shape.rzup, shape.rzdown);
        // continue polycone bookkeeping
        s.push_back(shape);
        logic.name_tag = shape.name_tag;
        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
        l.push_back(logic);
        pos.parent_tag = xml_fileident + ":" + xml_tracker;
        pos.child_tag = logic.shape_tag;
        p.push_back(pos);
        spec.name = xml_tob_subdet + xml_par_tail;
        spec.parameter.first = xml_tkddd_structure;
        spec.parameter.second = xml_tob;
        spec.partselectors.push_back(xml_tob);
        t.push_back(spec);
        spec.partselectors.clear();
        // define top-level endcap volume containers (polycones)
        analyseBackwardEndcapContainer(tr, shape.rzup, shape.rzdown);
        if (!(shape.rzup.empty() || shape.rzdown.empty())) {
            shape.name_tag = xml_tidb;
            s.push_back(shape);
            logic.name_tag = shape.name_tag;
            logic.shape_tag = xml_fileident + ":" + logic.name_tag;
            l.push_back(logic);
            pos.parent_tag = xml_fileident + ":" + xml_tracker;
            pos.child_tag = logic.shape_tag;
            p.push_back(pos);
            spec.partselectors.push_back(xml_tidb);
        }
        analyseForwardEndcapContainer(tr, shape.rzup, shape.rzdown);
        if (!(shape.rzup.empty() || shape.rzdown.empty())) {
            shape.name_tag = xml_tidf;
            s.push_back(shape);
            logic.name_tag = shape.name_tag;
            logic.shape_tag = xml_fileident + ":" + logic.name_tag;
            l.push_back(logic);
            pos.parent_tag = xml_fileident + ":" + xml_tracker;
            pos.child_tag = logic.shape_tag;
            p.push_back(pos);
            spec.partselectors.push_back(xml_tidf);
        }
        if (!spec.partselectors.empty()) {
            spec.name = xml_tid_subdet + xml_par_tail;
            spec.parameter.first = xml_tkddd_structure;
            spec.parameter.second = xml_tid;
            t.push_back(spec);
        }
        // translate entries in mt to elementary materials
        analyseElements(mt, e);
        // analyse barrel
        analyseLayers(mt, bc, tr, c, l, s, p, a, t);
        // analyse endcaps
        analyseDiscs(mt, ec, tr, c, l, s, p, a, t);
        // barrel services
        analyseBarrelServices(is, c, l, s, p, t);
        // endcap services
        analyseEndcapServices(is, c, l, s, p, t);
        // supports
        analyseSupports(is, c, l, s, p, t);
    }
    //protected
    void Extractor::analyseElements(MaterialTable&mattab, std::vector<Element>& elems) {
        for (unsigned int i = 0; i < mattab.rowCount(); i++) {
            Element e;
            MaterialRow& r = mattab.getMaterial(i);
            e.tag = r.tag;
            e.density = r.density;
            e.atomic_weight = pow((r.ilength / 35.), 3); // magic!
            e.atomic_number = Z(r.rlength, e.atomic_weight);
            elems.push_back(e);
        }
    }
    
    void Extractor::analyseBarrelContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
            std::vector<std::pair<double, double> >& down) {
        bool is_short, previous_short = false;
        std::pair<double, double> rz;
        double rmax = 0.0, zmax = 0.0, zmin = 0.0;
        unsigned int layer, n_of_layers = t.getBarrelLayers()->size();
        std::vector<Layer*>* bl = t.getBarrelLayers();
        up.clear();
        down.clear();
        for (layer = 0; layer < n_of_layers; layer++) {
            is_short = (bl->at(layer)->getMinZ() > 0) || (bl->at(layer)->getMaxZ() < 0);
            if (is_short) {
                //short layer on z- side
                if (bl->at(layer)->getMaxZ() < 0) {
                    //indices 0, 1
                    if ((layer == 0) || ((layer == 1) && (previous_short))) {
                        rz.first = bl->at(layer)->getMinRho();
                        rz.second = bl->at(layer)->getMinZ();
                        up.push_back(rz);
                    }
                    else {
                        //new barrel reached
                        if (bl->at(layer)->getMinZ() != zmin) {
                            //new layer sticks out compared to old layer
                            if (bl->at(layer)->getMinZ() > zmin) rz.first = bl->at(layer)->getMinRho();
                            //old layer sticks out compared to new layer
                            else rz.first = rmax;
                            rz.second = zmin;
                            up.push_back(rz);
                            rz.second = bl->at(layer)->getMinZ();
                            up.push_back(rz);
                        }
                    }
                    //indices size - 2, size - 1
                    if ((layer == n_of_layers - 1) || ((layer == n_of_layers - 2) && (previous_short))) {
                        rz.first = bl->at(layer)->getMaxRho();
                        rz.second = bl->at(layer)->getMinZ();
                        up.push_back(rz);
                    }
                }
                //short layer on z+ side
                else {
                    //indices 0, 1
                    if ((layer == 0) || ((layer == 1) && (previous_short))) {
                        rz.first = bl->at(layer)->getMinRho();
                        rz.second = bl->at(layer)->getMaxZ();
                        down.push_back(rz);
                    }
                    else {
                        //new barrel reached
                        if (bl->at(layer)->getMaxZ() != zmax) {
                            //new layer sticks out compared to old layer
                            if (bl->at(layer)->getMaxZ() > zmax) rz.first = bl->at(layer)->getMinRho();
                            //old layer sticks out compared to new layer
                            else rz.first = rmax;
                            rz.second = zmax;
                            down.push_back(rz);
                            rz.second = bl->at(layer)->getMaxZ();
                            down.push_back(rz);
                        }
                    }
                    //indices size - 2, size - 1
                    if ((layer == n_of_layers - 1) || ((layer == n_of_layers - 2) && (previous_short))) {
                        rz.first = bl->at(layer)->getMaxRho();
                        rz.second = bl->at(layer)->getMinZ();
                        down.push_back(rz);
                    }
                }
            }
            //regular layer across z=0
            else {
                //index 0
                if (layer == 0) {
                    rz.first = bl->at(layer)->getMinRho();
                    rz.second = bl->at(layer)->getMinZ();
                    up.push_back(rz);
                    rz.second = bl->at(layer)->getMaxZ();
                    down.push_back(rz);
                }
                else {
                    //new barrel reached
                    if (bl->at(layer)->getMaxZ() != zmax) {
                        //new layer sticks out compared to old layer
                        if (bl->at(layer)->getMaxZ() > zmax) rz.first = bl->at(layer)->getMinRho();
                        //old layer sticks out compared to new layer
                        else rz.first = rmax;
                        rz.second = zmin;
                        up.push_back(rz);
                        rz.second = zmax;
                        down.push_back(rz);
                        rz.second = bl->at(layer)->getMinZ();
                        up.push_back(rz);
                        rz.second = bl->at(layer)->getMaxZ();
                        down.push_back(rz);
                    }
                }
                //index size - 1
                if (layer == n_of_layers - 1) {
                    rz.first = bl->at(layer)->getMaxRho();
                    rz.second = bl->at(layer)->getMinZ();
                    up.push_back(rz);
                    rz.second = bl->at(layer)->getMaxZ();
                    down.push_back(rz);
                }
            }
            rmax = bl->at(layer)->getMaxRho();
            if (bl->at(layer)->getMinZ() < 0) zmin = bl->at(layer)->getMinZ();
            if (bl->at(layer)->getMaxZ() > 0) zmax = bl->at(layer)->getMaxZ();
            previous_short = is_short;
        }
    }

    void Extractor::analyseBackwardEndcapContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
            std::vector<std::pair<double, double> >& down) {
        std::vector<Layer*>* e = t.getEndcapLayers();
        int first, last, guard;
        first = 0;
        guard = e->size();
        last = first;
        while (last < guard) {
            if (e->at(last)->getMaxZ() > 0) break;
            last++;
        }
        analyseEndcapContainer(*e, first, last, up, down);
    }
    
    void Extractor::analyseForwardEndcapContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
            std::vector<std::pair<double, double> >& down) {
        std::vector<Layer*>* e = t.getEndcapLayers();
        int first, last;
        first = 0;
        last = e->size();
        while (first < last) {
            if (e->at(first)->getMaxZ() > 0) break;
            first++;
        }
        analyseEndcapContainer(*e, first, last, up, down);
    }
    
    void Extractor::analyseEndcapContainer(std::vector<Layer*>& el, int start, int stop,
            std::vector<std::pair<double, double> >& up, std::vector<std::pair<double, double> >& down) {
        std::pair<double, double> rz;
        double rmin = 0.0, rmax = 0.0, zmax = 0.0;
        up.clear();
        down.clear();
        for (int i = start; i < stop; i++) {
            // special treatment for first disc
            if (i == start) {
                rmin = el.at(i)->getMinRho();
                rmax = el.at(i)->getMaxRho();
                rz.first = rmax;
                rz.second = el.at(i)->getMinZ();
                up.push_back(rz);
                rz.first = rmin;
                down.push_back(rz);
            }
            // disc beyond the first
            else {
                // endcap change larger->smaller
                if (rmax > el.at(i)->getMaxRho()) {
                    rz.second = zmax;
                    rz.first = rmax;
                    up.push_back(rz);
                    rz.first = rmin;
                    down.push_back(rz);
                    rmax = el.at(i)->getMaxRho();
                    rmin = el.at(i)->getMinRho();
                    rz.first = rmax;
                    up.push_back(rz);
                    rz.first = rmin;
                    down.push_back(rz);
                }
                // endcap change smaller->larger
                if (rmax < el.at(i)->getMaxRho()) {
                    rz.second = el.at(i)->getMinZ();
                    rz.first = rmax;
                    up.push_back(rz);
                    rz.first = rmin;
                    down.push_back(rz);
                    rmax = el.at(i)->getMaxRho();
                    rmin = el.at(i)->getMinRho();
                    rz.first = rmax;
                    up.push_back(rz);
                    rz.second = rmin;
                    down.push_back(rz);
                }
            }
            zmax = el.at(i)->getMaxZ();
            // special treatment for last disc
            if (i == stop - 1) {
                rz.first = rmax;
                rz.second = zmax;
                up.push_back(rz);
                rz.first = rmin;
                down.push_back(rz);
            }
        }
    }
    
    void Extractor::analyseLayers(MaterialTable& mt, std::vector<std::vector<ModuleCap> >& bc, Tracker& tr,
            std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
            std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::vector<SpecParInfo>& t) {
        int layer;
        std::vector<std::vector<ModuleCap> >::iterator oiter, oguard;
        std::vector<ModuleCap>::iterator iiter, iguard;
        // container inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        AlgoInfo alg;
        SpecParInfo apv, lspec, rspec, mspec;
        shape.dxx = 0.0;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        pos.trans.dz = 0.0;
        pos.rot.name = "";
        pos.rot.phix = 0.0;
        pos.rot.phiy = 0.0;
        pos.rot.phiz = 0.0;
        pos.rot.thetax = 0.0;
        pos.rot.thetay = 0.0;
        pos.rot.thetaz = 0.0;
        lspec.name = xml_subdet_layer + xml_par_tail;
        lspec.parameter.first = xml_tkddd_structure;
        lspec.parameter.second = xml_det_layer;
        rspec.name = xml_subdet_rod + xml_par_tail;
        rspec.parameter.first = xml_tkddd_structure;
        rspec.parameter.second = xml_det_rod;
        mspec.name = xml_subdet_tobdet + xml_par_tail;
        mspec.parameter.first = xml_tkddd_structure;
        mspec.parameter.second = xml_det_tobdet;
        // b_mod: one composite for every module position on rod
        // s and l: one entry for every module position on rod (box), one for every layer (tube), rods TBD
        // p: one entry for every layer (two for short layers), two modules, one wafer and active for each ring on rod
        // a: rods within layer (twice in case of a short layer)
        layer = 1;
        alg.name = xml_tobalgo;
        oguard = bc.end();
        // barrel caps layer loop
        for (oiter = bc.begin(); oiter != oguard; oiter++) {
            double rmin = tr.getBarrelLayers()->at(layer - 1)->getMinRho();
            double rmax = tr.getBarrelLayers()->at(layer - 1)->getMaxRho();
            double zmin = tr.getBarrelLayers()->at(layer - 1)->getMinZ();
            double zmax = tr.getBarrelLayers()->at(layer - 1)->getMaxZ();
            double deltar = findDeltaR(tr.getBarrelLayers()->at(layer - 1)->getModuleVector()->begin(),
                    tr.getBarrelLayers()->at(layer - 1)->getModuleVector()->end(), (rmin + rmax) / 2.0);
            double ds, dt;
            int segs;
            if (deltar == 0.0) continue;
            bool is_short = (zmax < 0.0) || (zmin > 0.0);
            bool is_relevant = !is_short || (zmin > 0.0);
            if (is_relevant) {
                shape.type = bx;
                shape.rmin = 0.0;
                shape.rmax = 0.0;
                std::set<int> rings;
                std::ostringstream lname, rname, pconverter;
                lname << xml_layer << layer;
                rname << xml_rod << layer;
                iguard = oiter->end();
                // module caps loop
                for (iiter = oiter->begin(); iiter != iguard; iiter++) {
                    if (rings.find(iiter->getModule().getRing()) == rings.end()) {
                        segs = iiter->getModule().getNSegments();
                        std::vector<ModuleCap>::iterator partner;
                        std::ostringstream matname, shapename, specname;
                        // module composite material
                        matname << xml_base_actcomp << "L" << layer << "P" << iiter->getModule().getRing();
                        c.push_back(createComposite(matname.str(), compositeDensity(*iiter, true), *iiter, true));
                        // module box
                        shapename << iiter->getModule().getRing() << lname.str();
                        shape.name_tag = xml_barrel_module + shapename.str();
                        shape.dx = iiter->getModule().getModuleThickness() / 2.0;
                        shape.dy = iiter->getModule().getArea() / iiter->getModule().getHeight() / 2.0;
                        shape.dz = iiter->getModule().getHeight() / 2.0;
                        s.push_back(shape);
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                        logic.material_tag = xml_fileident + ":" + matname.str();
                        l.push_back(logic);
                        pos.child_tag = logic.shape_tag;
                        if ((iiter->getModule().getMeanPoint().Rho() > (rmax - deltar / 2.0))
                                || ((iiter->getModule().getMeanPoint().Rho() < ((rmin + rmax) / 2.0))
                                && (iiter->getModule().getMeanPoint().Rho() > (rmin + deltar / 2.0)))) pos.trans.dx = deltar / 2.0 - shape.dx;
                        else pos.trans.dx = shape.dx - deltar / 2.0;
                        if (is_short) {
                            pos.parent_tag = xml_fileident + ":" + rname.str() + xml_plus;
                            pos.trans.dz = iiter->getModule().getMinZ() - ((zmax + zmin) / 2.0) + shape.dz;
                            p.push_back(pos);
                            pos.parent_tag = xml_fileident + ":" + rname.str() + xml_minus;
                            pos.trans.dz = -pos.trans.dz;
                            pos.copy = 2;
                            p.push_back(pos);
                            pos.copy = 1;
                        }
                        else {
                            pos.parent_tag = xml_fileident + ":" + rname.str();
                            partner = findPartnerModule(iiter, iguard, iiter->getModule().getRing());
                            if (iiter->getModule().getMeanPoint().Z() > 0) {
                                pos.trans.dz = iiter->getModule().getMaxZ() - shape.dz;
                                p.push_back(pos);
                                if (partner != iguard) {
                                    if ((partner->getModule().getMeanPoint().Rho() > (rmax - deltar / 2.0))
                                            || ((partner->getModule().getMeanPoint().Rho() < ((rmin + rmax) / 2.0))
                                            && (partner->getModule().getMeanPoint().Rho() > (rmin + deltar / 2.0)))) pos.trans.dx = deltar / 2.0 - shape.dx;
                                    else pos.trans.dx = shape.dx - deltar / 2.0;
                                    pos.trans.dz = partner->getModule().getMaxZ() - shape.dz;
                                    pos.copy = 2;
                                    p.push_back(pos);
                                    pos.copy = 1;
                                }
                            }
                            else {
                                pos.trans.dz = iiter->getModule().getMaxZ() - shape.dz;
                                pos.copy = 2;
                                p.push_back(pos);
                                pos.copy = 1;
                                if (partner != iguard) {
                                    if ((partner->getModule().getMeanPoint().Rho() > (rmax - deltar / 2.0))
                                            || ((partner->getModule().getMeanPoint().Rho() < ((rmin + rmax) / 2.0))
                                            && (partner->getModule().getMeanPoint().Rho() > (rmin + deltar / 2.0)))) pos.trans.dx = deltar / 2.0 - shape.dx;
                                    else pos.trans.dx = shape.dx - deltar / 2.0;
                                    pos.trans.dz = partner->getModule().getMaxZ() - shape.dz;
                                    p.push_back(pos);
                                }
                            }
                        }
                        // wafer
                        shape.name_tag = xml_barrel_module + shapename.str() + xml_base_waf;
                        shape.dx = calculateSensorThickness(*iiter, mt) / 2.0;
                        shape.dz = shape.dz / (double)(segs);
                        s.push_back(shape);
                        pos.parent_tag = logic.shape_tag;
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                        l.push_back(logic);
                        pos.child_tag = logic.shape_tag;
                        pos.trans.dx = shape.dx - iiter->getModule().getModuleThickness() / 2.0;
                        for (int i = 0; i < segs; i++) {
                            pos.copy = i + 1;
                            pos.trans.dz = (1 + 2 * i) * shape.dz - iiter->getModule().getHeight() / 2.0;
                            p.push_back(pos);
                        }
                        // active surface
                        shape.name_tag = xml_barrel_module + shapename.str() + xml_base_act;
                        s.push_back(shape);
                        pos.parent_tag = logic.shape_tag;
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                        logic.material_tag = xml_fileident + ":" + xml_sensor_silicon;
                        l.push_back(logic);
                        pos.child_tag = logic.shape_tag;
                        pos.copy = 1;
                        pos.trans.dz = 0.0;
                        p.push_back(pos);
                        // topology
                        mspec.partselectors.push_back(logic.name_tag);
                        specname << xml_apv_head << (iiter->getModule().getNStripsAcross() / 128) << xml_par_tail;
                        int id = findSpecParIndex(t, specname.str());
                        if (id >= 0) t.at(id).partselectors.push_back(logic.name_tag);
                        else {
                            apv.partselectors.clear();
                            apv.name = specname.str();
                            apv.parameter.first = xml_apv_number;
                            specname.str("");
                            specname << (iiter->getModule().getNStripsAcross() / 128);
                            apv.parameter.second = specname.str();
                            apv.partselectors.push_back(logic.name_tag);
                            t.push_back(apv);
                        }
                        rings.insert(iiter->getModule().getRing());
                    }
                }
                // rod(s)
                shape.name_tag = rname.str();
                if (is_short) shape.name_tag = shape.name_tag + xml_plus;
                dt = shape.dx;
                shape.dx = deltar / 2.0;
                if (is_short) shape.dz = (zmax - zmin) / 2.0;
                else shape.dz = zmax;
                s.push_back(shape);
                logic.name_tag = shape.name_tag;
                logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                logic.material_tag = xml_material_air;
                l.push_back(logic);
                rspec.partselectors.push_back(logic.name_tag);
                pconverter << logic.shape_tag;
                if (is_short) {
                    shape.name_tag = rname.str() + xml_minus;
                    s.push_back(shape);
                    logic.name_tag = shape.name_tag;
                    logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                    l.push_back(logic);
                    rspec.partselectors.push_back(logic.name_tag);
                }
                ds = fromRim(rmax, shape.dy);
                // layer
                shape.type = tb;
                shape.dx = 0.0;
                shape.dy = 0.0;
                pos.trans.dx = 0.0;
                pos.trans.dz = 0.0;
                shape.name_tag = lname.str();
                //if (is_short) shape.name_tag = shape.name_tag + xml_plus;
                shape.rmin = rmin;
                shape.rmax = rmax;
                //if (is_short) shape.dz = (zmax - zmin) / 2.0;
                //else shape.dz = zmax;
                shape.dz = zmax; //
                s.push_back(shape);
                logic.name_tag = shape.name_tag;
                logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                l.push_back(logic);
                pos.parent_tag = xml_fileident + ":" + xml_tob;
                pos.child_tag = logic.shape_tag;
                //if (is_short) pos.trans.dz = zmin + (zmax - zmin) / 2.0;
                p.push_back(pos);
                lspec.partselectors.push_back(logic.name_tag);
                // modules in rod algorithm(s)
                alg.parent = logic.shape_tag;
                alg.parameters.push_back(stringParam(xml_childparam, pconverter.str()));
                pconverter.str("");
                pconverter << (tr.getBarrelLayers()->at(layer - 1)->getTilt() + 90) << "*deg";
                alg.parameters.push_back(numericParam(xml_tilt, pconverter.str()));
                pconverter.str("");
                pconverter << tr.getBarrelLayers()->at(layer - 1)->getStartAngle();
                alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
                pconverter.str("");
                alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
                pconverter << (rmin + deltar / 2.0) << "*mm";
                alg.parameters.push_back(numericParam(xml_radiusin, pconverter.str()));
                pconverter.str("");
                pconverter << (rmax - ds - deltar / 2.0 - 2.0 * dt) << "*mm";
                alg.parameters.push_back(numericParam(xml_radiusout, pconverter.str()));
                pconverter.str("");
                if (is_short) {
                    pconverter << (zmin + (zmax - zmin) / 2.0) << "*mm";
                    alg.parameters.push_back(numericParam(xml_zposition, pconverter.str()));
                    pconverter.str("");
                }
                else alg.parameters.push_back(numericParam(xml_zposition, "0.0*mm"));
                pconverter << static_cast<BarrelLayer*>(tr.getBarrelLayers()->at(layer - 1))->getRods();
                alg.parameters.push_back(numericParam(xml_number, pconverter.str()));
                alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
                alg.parameters.push_back(numericParam(xml_incrcopyno, "1"));
                a.push_back(alg);
                // extras for short layers
                if (is_short) {
                    //shape.name_tag = lname.str() + xml_minus;
                    //s.push_back(shape);
                    //logic.name_tag = shape.name_tag;
                    //logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                    //l.push_back(logic);
                    //pos.child_tag = logic.shape_tag;
                    //pos.trans.dz = -(zmin + (zmax - zmin) / 2.0);
                    //p.push_back(pos);
                    //alg.parent = logic.shape_tag;
                    pconverter.str("");
                    pconverter << xml_fileident << ":" << rname.str() << xml_minus;
                    alg.parameters.front() = stringParam(xml_childparam, pconverter.str());
                    pconverter.str("");//
                    pconverter << -(zmin + (zmax - zmin) / 2.0) << "*mm";//
                    alg.parameters.at(6) = numericParam(xml_zposition, pconverter.str());//
                    a.push_back(alg);
                    //lspec.partselectors.push_back(logic.name_tag);
                }
                alg.parameters.clear();
            }
            layer++;
        }
        if (!lspec.partselectors.empty()) t.push_back(lspec);
        if (!rspec.partselectors.empty()) t.push_back(rspec);
        if (!mspec.partselectors.empty()) t.push_back(mspec);
    }
    
    void Extractor::analyseDiscs(MaterialTable& mt, std::vector<std::vector<ModuleCap> >& ec, Tracker& tr,
            std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
            std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::vector<SpecParInfo>& t) {
        int layer;
        std::vector<std::vector<ModuleCap> >::iterator oiter, oguard;
        std::vector<ModuleCap>::iterator iiter, iguard;
        // container inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        AlgoInfo alg;
        SpecParInfo apv, dspec, rspec, mspec;
        shape.dxx = 0.0;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        pos.trans.dz = 0.0;
        pos.rot.name = "";
        pos.rot.phix = 0.0;
        pos.rot.phiy = 0.0;
        pos.rot.phiz = 0.0;
        pos.rot.thetax = 0.0;
        pos.rot.thetay = 0.0;
        pos.rot.thetaz = 0.0;
        dspec.name = xml_subdet_wheel + xml_par_tail;
        dspec.parameter.first = xml_tkddd_structure;
        dspec.parameter.second = xml_det_wheel;
        rspec.name = xml_subdet_ring + xml_par_tail;
        rspec.parameter.first = xml_tkddd_structure;
        rspec.parameter.second = xml_det_ring;
        mspec.name = xml_subdet_tiddet + xml_par_tail;
        mspec.parameter.first = xml_tkddd_structure;
        mspec.parameter.second = xml_det_tiddet;
        // e_mod: one composite for every ring
        // s and l: one entry for every ring module, one for every ring, one for every disc
        // p: one entry for every disc, one for every ring, one module, wafer and active per ring
        // a: two per ring with modules inside ring
        layer = 1;
        alg.name = xml_ecalgo;
        oguard = ec.end();
        // endcap caps layer loop
        for (oiter = ec.begin(); oiter != oguard; oiter++) {
            std::set<int> ridx;
            std::map<int, RingInfo> rinfo;
            double rmin = tr.getEndcapLayers()->at(layer - 1)->getMinRho();
            double rmax = tr.getEndcapLayers()->at(layer - 1)->getMaxRho();
            double zmax = tr.getEndcapLayers()->at(layer - 1)->getMaxZ();
            double zmin = tr.getEndcapLayers()->at(layer - 1)->getMinZ();
            bool plus = tr.getEndcapLayers()->at(layer - 1)->getMaxZ() > 0;
            std::ostringstream dname, pconverter;
            dname << xml_disc << layer;
            shape.type = tp;
            shape.rmin = 0.0;
            shape.rmax = 0.0;
            pos.trans.dz = 0.0;
            iguard = oiter->end();
            for (iiter = oiter->begin(); iiter != iguard; iiter++) {
                if (ridx.find(iiter->getModule().getRing()) == ridx.end()) {
                    ridx.insert(iiter->getModule().getRing());
                    std::ostringstream matname, rname, mname, specname;
                    matname << xml_base_actcomp << "D" << layer << "R" << iiter->getModule().getRing();
                    c.push_back(createComposite(matname.str(), compositeDensity(*iiter, true), *iiter, true));
                    rname << xml_ring << iiter->getModule().getRing() << dname.str();
                    mname << xml_endcap_module << iiter->getModule().getRing() << dname.str();
                    // collect ring info
                    RingInfo ri;
                    ri.name = rname.str();
                    ri.childname = mname.str();
                    ri.fw = (iiter->getModule().getMeanPoint().Z() < (zmin + zmax) / 2.0);
                    ri.modules = static_cast<EndcapLayer*>(tr.getEndcapLayers()->at(layer - 1))->getModulesOnRing().at(iiter->getModule().getRing() - 1);
                    ri.rin = iiter->getModule().getMinRho();
                    ri.rout = iiter->getModule().getMaxRho();
                    ri.rmid = iiter->getModule().getMeanPoint().Rho();
                    ri.mthk = iiter->getModule().getModuleThickness();
                    ri.phi = iiter->getModule().getMeanPoint().Phi();
                    rinfo.insert(std::pair<int, RingInfo>(iiter->getModule().getRing(), ri));
                    // module trapezoid
                    shape.name_tag = mname.str();
                    shape.dx = static_cast<EndcapModule&>(iiter->getModule()).getWidthLo() / 2.0;
                    shape.dxx = static_cast<EndcapModule&>(iiter->getModule()).getWidthHi() / 2.0;
                    shape.dy = iiter->getModule().getHeight() / 2.0;
                    shape.dz = iiter->getModule().getModuleThickness() / 2.0;
                    s.push_back(shape);
                    logic.name_tag = shape.name_tag;
                    logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                    logic.material_tag = xml_fileident + ":" + matname.str();
                    l.push_back(logic);
                    // wafer
                    pos.parent_tag = logic.shape_tag;
                    shape.name_tag = mname.str() + xml_base_waf;
                    shape.dz = calculateSensorThickness(*iiter, mt) / 2.0;
                    s.push_back(shape);
                    logic.name_tag = shape.name_tag;
                    logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                    logic.material_tag = xml_material_air;
                    l.push_back(logic);
                    pos.child_tag = logic.shape_tag;
                    if (iiter->getModule().getMaxZ() > 0) pos.trans.dz = shape.dz - iiter->getModule().getModuleThickness() / 2.0;
                    else pos.trans.dz = iiter->getModule().getModuleThickness() / 2.0 - shape.dz;
                    p.push_back(pos);
                    pos.trans.dz = 0.0;
                    // active surface
                    pos.parent_tag = logic.shape_tag;
                    shape.name_tag = mname.str() + xml_base_act;
                    s.push_back(shape);
                    logic.name_tag = shape.name_tag;
                    logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                    logic.material_tag = xml_fileident + ":" + xml_sensor_silicon;
                    l.push_back(logic);
                    pos.child_tag = logic.shape_tag;
                    p.push_back(pos);
                    // topology
                    mspec.partselectors.push_back(logic.name_tag);
                    specname << xml_apv_head << (iiter->getModule().getNStripsAcross() / 128) << xml_par_tail;
                    int id = findSpecParIndex(t, specname.str());
                    if (id >= 0) t.at(id).partselectors.push_back(logic.name_tag);
                    else {
                        apv.partselectors.clear();
                        apv.name = specname.str();
                        apv.parameter.first = xml_apv_number;
                        apv.parameter.second = iiter->getModule().getNStripsAcross() / 128;
                        apv.partselectors.push_back(logic.name_tag);
                        t.push_back(apv);
                    }
                }
            }
            // rings
            shape.type = tb;
            shape.dx = 0.0;
            shape.dxx = 0.0;
            shape.dy = 0.0;
            shape.dz = findDeltaZ(tr.getEndcapLayers()->at(layer - 1)->getModuleVector()->begin(),
                    tr.getEndcapLayers()->at(layer - 1)->getModuleVector()->end(), (zmin + zmax) / 2.0) / 2.0;
            pos.parent_tag = xml_fileident + ":" + dname.str();
            std::set<int>::const_iterator siter, sguard = ridx.end();
            for (siter = ridx.begin(); siter != sguard; siter++) {
                if (rinfo[*siter].modules > 0) {
                    shape.name_tag = rinfo[*siter].name;
                    shape.rmin = rinfo[*siter].rin;
                    shape.rmax = rinfo[*siter].rout;
                    s.push_back(shape);
                    logic.name_tag = shape.name_tag;
                    logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                    logic.material_tag = xml_material_air;
                    l.push_back(logic);
                    pos.child_tag = logic.shape_tag;
                    if (rinfo[*siter].fw) pos.trans.dz = (zmin - zmax) / 2.0 + shape.dz;
                    else pos.trans.dz = (zmax - zmin) / 2.0 - shape.dz;
                    p.push_back(pos);
                    rspec.partselectors.push_back(logic.name_tag);
                    alg.parent = logic.shape_tag;
                    alg.parameters.push_back(stringParam(xml_childparam, xml_fileident + ":" + rinfo[*siter].childname));
                    pconverter << (rinfo[*siter].modules / 2);
                    alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
                    pconverter.str("");
                    alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
                    alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
                    alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
                    pconverter << rinfo[*siter].phi;
                    alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
                    pconverter.str("");
                    pconverter << rinfo[*siter].rmid;
                    alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
                    pconverter.str("");
                    alg.parameters.push_back(vectorParam(0, 0, shape.dz - rinfo[*siter].mthk / 2.0));
                    a.push_back(alg);
                    alg.parameters.clear();
                    alg.parameters.push_back(stringParam(xml_childparam, xml_fileident + ":" + rinfo[*siter].childname));
                    pconverter << (rinfo[*siter].modules / 2);
                    alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
                    pconverter.str("");
                    alg.parameters.push_back(numericParam(xml_startcopyno, "2"));
                    alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
                    alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
                    pconverter << (rinfo[*siter].phi + 2 * PI / (double)(rinfo[*siter].modules));
                    alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
                    pconverter.str("");
                    pconverter << rinfo[*siter].rmid;
                    alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
                    pconverter.str("");
                    alg.parameters.push_back(vectorParam(0, 0, rinfo[*siter].mthk / 2.0 - shape.dz));
                    a.push_back(alg);
                    alg.parameters.clear();
                }
            }
            pos.trans.dz = 0.0;
            //disc
            shape.name_tag = dname.str();
            shape.rmin = rmin;
            shape.rmax = rmax;
            shape.dz = (zmax - zmin) / 2.0;
            s.push_back(shape);
            logic.name_tag = shape.name_tag;
            logic.shape_tag = xml_fileident + ":" + logic.name_tag;
            logic.material_tag = xml_material_air;
            l.push_back(logic);
            pos.parent_tag = xml_fileident + ":";
            if (plus) pos.parent_tag = pos.parent_tag + xml_tidf;
            else pos.parent_tag = pos.parent_tag + xml_tidb;
            pos.child_tag = logic.shape_tag;
            pos.trans.dz = (zmax + zmin) / 2.0;
            p.push_back(pos);
            dspec.partselectors.push_back(logic.name_tag);
            layer++;
        }
        if (!dspec.partselectors.empty()) t.push_back(dspec);
        if (!rspec.partselectors.empty()) t.push_back(rspec);
        if (!mspec.partselectors.empty()) t.push_back(mspec);
    }
    
    void Extractor::analyseBarrelServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
            std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<SpecParInfo>& t) {
        // container inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        shape.type = tb;
        shape.dx = 0.0;
        shape.dxx = 0.0;
        shape.dy = 0.0;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        pos.rot.name = "";
        pos.rot.phix = 0.0;
        pos.rot.phiy = 0.0;
        pos.rot.phiz = 0.0;
        pos.rot.thetax = 0.0;
        pos.rot.thetay = 0.0;
        pos.rot.thetaz = 0.0;
        // b_ser: one composite for every service volume on the z+ side
        // s, l and p: one entry per service volume
        std::vector<InactiveElement>::iterator iter, guard;
        std::vector<InactiveElement>& bs = is.getBarrelServices();
        guard = bs.end();
        for (iter = bs.begin(); iter != guard; iter++) {
            std::ostringstream matname, shapename;
            matname << xml_base_serfcomp << iter->getCategory() << "R" << (int)(iter->getInnerRadius()) << "dZ" << (int)(iter->getZLength());
            shapename << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(iter->getZOffset());
            if ((iter->getZOffset() + iter->getZLength()) > 0) c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
            shape.name_tag = shapename.str();
            shape.dz = iter->getZLength() / 2.0;
            shape.rmin = iter->getInnerRadius();
            shape.rmax = shape.rmin + iter->getRWidth();
            s.push_back(shape);
            logic.name_tag = shapename.str();
            logic.shape_tag = xml_fileident + ":" + shapename.str();
            logic.material_tag = xml_fileident + ":" + matname.str();
            l.push_back(logic);
            pos.parent_tag = xml_fileident + ":" + xml_tracker;
            pos.child_tag = logic.shape_tag;
            pos.trans.dz = iter->getZOffset() + shape.dz;
            p.push_back(pos);
        }
    }
    
    void Extractor::analyseEndcapServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
            std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<SpecParInfo>& t) {
        // container inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        shape.type = tb;
        shape.dx = 0.0;
        shape.dxx = 0.0;
        shape.dy = 0.0;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        pos.rot.name = "";
        pos.rot.phix = 0.0;
        pos.rot.phiy = 0.0;
        pos.rot.phiz = 0.0;
        pos.rot.thetax = 0.0;
        pos.rot.thetay = 0.0;
        pos.rot.thetaz = 0.0;
        // e_ser: one composite for every service volume on the z+ side
        // s, l and p: one entry per service volume
        std::vector<InactiveElement>::iterator iter, guard;
        std::vector<InactiveElement>& es = is.getEndcapServices();
        guard = es.end();
        for (iter = es.begin(); iter != guard; iter++) {
            std::ostringstream matname, shapename;
            matname << xml_base_serfcomp << iter->getCategory() << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0));
            shapename << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0));
            if ((iter->getZOffset() + iter->getZLength()) > 0) c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
            shape.name_tag = shapename.str();
            shape.dz = iter->getZLength() / 2.0;
            shape.rmin = iter->getInnerRadius();
            shape.rmax = shape.rmin + iter->getRWidth();
            s.push_back(shape);
            logic.name_tag = shapename.str();
            logic.shape_tag = xml_fileident + ":" + shapename.str();
            logic.material_tag = xml_fileident + ":" + matname.str();
            l.push_back(logic);
            pos.parent_tag = xml_fileident + ":" + xml_tracker;
            pos.child_tag = logic.shape_tag;
            pos.trans.dz = iter->getZOffset() + shape.dz;
            p.push_back(pos);
        }
    }
    
    void Extractor::analyseSupports(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
            std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<SpecParInfo>& t) {
        // container inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        shape.type = tb;
        shape.dx = 0.0;
        shape.dxx = 0.0;
        shape.dy = 0.0;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        pos.rot.name = "";
        pos.rot.phix = 0.0;
        pos.rot.phiy = 0.0;
        pos.rot.phiz = 0.0;
        pos.rot.thetax = 0.0;
        pos.rot.thetay = 0.0;
        pos.rot.thetaz = 0.0;
        // b_sup, e_sup, o_sup, t_sup, u_sup: one composite per category
        // l, s and p: one entry per support part
        std::set<MaterialProperties::Category> found;
        std::set<MaterialProperties::Category>::iterator fres;
        std::vector<InactiveElement>::iterator iter, guard;
        std::vector<InactiveElement>& sp = is.getSupports();
        guard = sp.end();
        for (iter = sp.begin(); iter != guard; iter++) {
            std::ostringstream matname, shapename;
            matname << xml_base_lazycomp << iter->getCategory();
            shapename << xml_base_lazy << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset()));
            fres = found.find(iter->getCategory());
            if (fres == found.end()) {
                c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
                found.insert(iter->getCategory());
            }
            shape.name_tag = shapename.str();
            shape.dz = iter->getZLength() / 2.0;
            shape.rmin = iter->getInnerRadius();
            shape.rmax = shape.rmin + iter->getRWidth();
            s.push_back(shape);
            logic.name_tag = shapename.str();
            logic.shape_tag = xml_fileident + ":" + shapename.str();
            logic.material_tag = xml_fileident + ":" + matname.str();
            l.push_back(logic);
            pos.parent_tag = xml_fileident + ":" + xml_tracker;
            pos.child_tag = logic.shape_tag;
            if ((iter->getCategory() == MaterialProperties::o_sup) ||
                    (iter->getCategory() == MaterialProperties::t_sup)) pos.trans.dz = 0.0;
            else pos.trans.dz = iter->getZOffset() + shape.dz;
            p.push_back(pos);
        }
    }
    //private
    Composite Extractor::createComposite(std::string name, double density, MaterialProperties& mp, bool nosensors) {
        Composite comp;
        comp.name = name;
        comp.density = density;
        comp.method = wt;
        double m = 0.0;
        for (unsigned int i = 0; i < mp.localMassCount(); i++) {
            if (!nosensors || (mp.getLocalTag(i).compare(xml_sensor_silicon) != 0)) {
                std::pair<std::string, double> p;
                p.first = mp.getLocalTag(i);
                p.second = mp.getLocalMass(i);
                comp.elements.push_back(p);
                m = m + mp.getLocalMass(i);
            }
        }
        for (unsigned int i = 0; i < mp.exitingMassCount(); i++) {
            if (!nosensors || (mp.getExitingTag(i).compare(xml_sensor_silicon) != 0)) {
                std::pair<std::string, double> p;
                p.first = mp.getExitingTag(i);
                p.second = mp.getExitingMass(i);
                bool found = false;
                std::vector<std::pair<std::string, double> >::iterator iter, guard = comp.elements.end();
                for (iter = comp.elements.begin(); iter != guard; iter++) {
                    if (iter->first == p.first) {
                        found = true;
                        break;
                    }
                }
                if (found) iter->second = iter->second + p.second;
                else comp.elements.push_back(p);
                m = m + mp.getExitingMass(i);
            }
        }
        for (unsigned int i = 0; i < comp.elements.size(); i++)
            comp.elements.at(i).second = comp.elements.at(i).second / m;
        return comp;
    }
    
    std::vector<ModuleCap>::iterator Extractor::findPartnerModule(std::vector<ModuleCap>::iterator i,
            std::vector<ModuleCap>::iterator g, int ponrod, bool find_first) {
        std::vector<ModuleCap>::iterator res = i;
        if (i != g) {
            bool plus = false;
            if (!find_first) plus = i->getModule().getMeanPoint().Z() > 0;
            while (res != g) {
                if (res->getModule().getRing() == ponrod) {
                    if (find_first) break;
                    else {
                        if((plus && (res->getModule().getMeanPoint().Z() < 0))
                                || (!plus && (res->getModule().getMeanPoint().Z() > 0))) break;
                    }
                }
                res++;
            }
        }
        return res;
    }
    
    double Extractor::findDeltaR(std::vector<Module*>::iterator start,
            std::vector<Module*>::iterator stop, double middle) {
        std::vector<Module*>::iterator iter, mod1, mod2;
        double dr = 0.0;
        iter = start;
        mod1 = stop;
        mod2 = stop;
        for (iter = start; iter != stop; iter++) {
            if ((*iter)->getMeanPoint().Rho() > middle) {
                mod1 = iter;
                break;
            }
        }
        for (iter = mod1; iter != stop; iter++) {
            if ((*iter)->getMeanPoint().Rho() > middle) {
                if ((*iter)->getMeanPoint().Rho() < (*mod1)->getMeanPoint().Rho()) {
                    mod2 = iter;
                    break;
                }
                else if (!((*iter)->getMeanPoint().Rho() == (*mod1)->getMeanPoint().Rho())) {
                    mod2 = mod1;
                    mod1 = iter;
                    break;
                }
            }
        }
        dr = (*mod1)->getMinRho() - (*mod2)->getMinRho() + (*mod1)->getModuleThickness();
        return dr;
    }
    
    double Extractor::findDeltaZ(std::vector<Module*>::iterator start,
            std::vector<Module*>::iterator stop, double middle) {
        std::vector<Module*>::iterator iter, mod1, mod2;
        double dz = 0.0;
        iter = start;
        mod1 = stop;
        mod2 = stop;
        for (iter = start; iter != stop; iter++) {
            if ((*iter)->getMinZ() > middle) {
                mod1 = iter;
                break;
            }
        }
        for (iter = mod1; iter != stop; iter++) {
            if ((*iter)->getMinZ() > middle) {
                if ((*iter)->getMinZ() < (*mod1)->getMinZ()) {
                    mod2 = iter;
                    break;
                }
                else if (!((*iter)->getMinZ() == (*mod1)->getMinZ())) {
                    mod2 = mod1;
                    mod1 = iter;
                    break;
                }
            }
        }
        dz = (*mod1)->getMaxZ() - (*mod2)->getMinZ();
        return dz;
    }
    
    int Extractor::findSpecParIndex(std::vector<SpecParInfo>& specs, std::string label) {
        int idx = 0, size = (int)(specs.size());
        while (idx < size) {
            if (specs.at(idx).name.compare(label) == 0) return idx;
            idx++;
        }
        return -1;
    }
    
    double Extractor::calculateSensorThickness(ModuleCap& mc, MaterialTable& mt) {
        double t = 0.0;
        double m = 0.0, d = 0.0;
        for (unsigned int i = 0; i < mc.localMassCount(); i++) {
            if (mc.getLocalTag(i).compare(xml_sensor_silicon) == 0) m = m + mc.getLocalMass(i);
        }
        for (unsigned int i = 0; i < mc.exitingMassCount(); i++) {
            if (mc.getExitingTag(i).compare(xml_sensor_silicon) == 0) m = m + mc.getExitingMass(i);
        }
        try { d = mt.getMaterial(xml_sensor_silicon).density; }
        catch (std::exception& e) { return 0.0; }
        t = 1000 * m / (d * mc.getSurface());
        return t;
    }
    
    std::string Extractor::stringParam(std::string name, std::string value) {
        std::string res;
        res = xml_algorithm_string + name + xml_algorithm_value + value + xml_general_endline;
        return res;
    }
    
    std::string Extractor::numericParam(std::string name, std::string value) {
        std::string res;
        res = xml_algorithm_numeric + name + xml_algorithm_value + value + xml_general_endline;
        return res;
    }
    
    std::string Extractor::vectorParam(double x, double y, double z) {
        std::ostringstream res;
        res << xml_algorithm_vector_open << x << "," << y << "," << z << xml_algorithm_vector_close;
        return res.str();
    }
    
    double Extractor::compositeDensity(ModuleCap& mc, bool nosensors) {
        double d = mc.getSurface() * mc.getModule().getModuleThickness();
        if (nosensors) {
            double m = 0.0;
            for (unsigned int i = 0; i < mc.localMassCount(); i++) {
                if (mc.getLocalTag(i).compare(xml_sensor_silicon) != 0) m = m + mc.getLocalMass(i);
            }
            for (unsigned int i = 0; i < mc.exitingMassCount(); i++) {
                if (mc.getExitingTag(i).compare(xml_sensor_silicon) != 0) m = m + mc.getExitingMass(i);
            }
            d = 1000 * m / d;
        }
        else d = 1000 * mc.getTotalMass() / d;
        return d;
    }
    
    double Extractor::compositeDensity(InactiveElement& ie) {
        double d = ie.getRWidth() + ie.getInnerRadius();
        d = d * d - ie.getInnerRadius() * ie.getInnerRadius();
        d = 1000 * ie.getTotalMass() / (PI * ie.getZLength() * d);
        return d;
    }
    
    double Extractor::fromRim(double r, double w) {
        double s = asin(w / r);
        s = 1 - cos(s);
        s = s * r;
        return s;
    }
                
    int Extractor::Z(double x0, double A) {
        double d = 4 - 4 * (1.0 - 181.0 * A / x0);
        if (d > 0) return floor((sqrt(d) - 2.0) / 2.0 + 0.5);
        else return -1;
    }
}
