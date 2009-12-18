/**
 * @file tk2CMSSW.cc
 * @brief
 */

#include <tk2CMSSW.h>
namespace insur {
    // public
    void tk2CMSSW::translate(MaterialTable& mt, MaterialBudget& mb, std::string outsubdir) {
        std::string outpath = default_xmlpath;
        if (outsubdir.empty()) outpath = outpath + "/" + default_xml;
        else {
            if (outsubdir.at(0) == '/') outpath = outpath + outsubdir;
            else outpath = outpath + "/" + outsubdir;
        }
        if(outpath.at(outpath.size() - 1) != '/') outpath = outpath + "/";
        // analyse tracker system and build up collection of elements, composites, hierarchy, shapes, position and topology
        analyse(mt, mb, elements, composites, logic, shapes, positions, algos, specs);
        // translate collected information to XML and write to buffers
        std::ostringstream gbuffer, tbuffer;
        gbuffer << xml_preamble << xml_const_section;
        tbuffer << xml_preamble;
        materialSection(xml_trackerfile, elements, composites, gbuffer);
        logicalPartSection(logic, xml_trackerfile, gbuffer);
        solidSection(shapes, xml_trackerfile, gbuffer);
        posPartSection(positions, algos, xml_trackerfile, gbuffer);
        specParSection(specs, xml_specpars_label, tbuffer);
        gbuffer << xml_defclose;
        tbuffer << xml_defclose;
        // write contents of buffer to top-level file
        bfs::remove_all(outpath.c_str());
        bfs::create_directory(outpath);
        std::ofstream goutstream((outpath + xml_trackerfile).c_str());
        goutstream << gbuffer.str() << std::endl;
        goutstream.close();
        std::cout << "CMSSW XML output has been written to " << outpath << xml_trackerfile << std::endl;
        std::ofstream toutstream((outpath + xml_topologyfile).c_str());
        toutstream << tbuffer.str() << std::endl;
        toutstream.close();
        std::cout << "CMSSW topology output has been written to " << outpath << xml_topologyfile << std::endl;
    }
    
    // protected
    void tk2CMSSW::materialSection(std::string name , std::vector<Element>& e, std::vector<Composite>& c, std::ostringstream& stream) {
        stream << xml_material_section_open << name << xml_general_inter;
        for (unsigned int i = 0; i < e.size(); i++) elementaryMaterial(e.at(i).tag, e.at(i).density, e.at(i).atomic_number, e.at(i).atomic_weight, stream);
        for (unsigned int i = 0; i < c.size(); i++) compositeMaterial(c.at(i).name, c.at(i).density, c.at(i).method, c.at(i).elements, stream);
        stream << xml_material_section_close;
    }
    
    void tk2CMSSW::logicalPartSection(std::vector<LogicalInfo>& l, std::string label, std::ostringstream& stream) {
        std::vector<LogicalInfo>::const_iterator iter, guard = l.end();
        stream << xml_logical_part_section_open << label << xml_general_inter;
        for (iter = l.begin(); iter != guard; iter++) logicalPart(iter->name_tag, iter->shape_tag, iter->material_tag, stream);
        stream << xml_logical_part_section_close;
    }
    
    void tk2CMSSW::solidSection(std::vector<ShapeInfo>& s, std::string label, std::ostringstream& stream) {
        stream << xml_solid_section_open << label << xml_general_inter;
        for (unsigned int i = 0; i < s.size(); i++) {
            switch (s.at(i).type) {
                case bx : box(s.at(i).name_tag, s.at(i).dx, s.at(i).dy, s.at(i).dz, stream);
                break;
                case tp : trapezoid(s.at(i).name_tag, s.at(i).dx, s.at(i).dxx, s.at(i).dy, s.at(i).dz, stream);
                break;
                case tb : tubs(s.at(i).name_tag, s.at(i).rmin, s.at(i).rmax, s.at(i).dz, stream);
                break;
                case pc : polycone(s.at(i).name_tag, s.at(i).rzup, s.at(i).rzdown, stream);
                break;
                default: std::cerr << "solidSection(): unknown shape type found. Using box." << std::endl;
                box(s.at(i).name_tag, s.at(i).dx, s.at(i).dy, s.at(i).dz, stream);
            }
        }
        stream << xml_solid_section_close;
    }
    
    void tk2CMSSW::posPartSection(std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::string label, std::ostringstream& stream) {
        std::vector<PosInfo>::iterator piter, pguard = p.end();
        std::vector<AlgoInfo>::iterator aiter, aguard = a.end();
        stream << xml_pos_part_section_open << label << xml_general_inter;
        for (piter = p.begin(); piter != pguard; piter++) posPart(piter->parent_tag, piter->child_tag, piter->rot, piter->trans, piter->copy, stream);
        for (aiter = a.begin(); aiter != aguard; aiter++) algorithm(aiter->name, aiter->parent, aiter->parameters, stream);
        stream << xml_pos_part_section_close;
    }
    
    void tk2CMSSW::specParSection(std::vector<SpecPar>& t, std::string label, std::ostringstream& stream) {
        std::vector<SpecPar>::iterator titer, tguard = t.end();
        stream << xml_spec_part_section_open << label << xml_general_inter;
        for (titer = t.begin(); titer != tguard; titer++) specPar(titer->name, titer->parameter, titer->partselectors, stream);
        stream << xml_spec_part_section_close;
    }
    
    void tk2CMSSW::algorithm(std::string name, std::string parent,
            std::vector<std::string>& params, std::ostringstream& stream) {
        stream << xml_algorithm_open << name << xml_algorithm_parent << parent << xml_general_endline;
        for (unsigned int i = 0; i < params.size(); i++) stream << params.at(i);
        stream << xml_algorithm_close;
    }
    
    void tk2CMSSW::elementaryMaterial(std::string tag, double density, int a_number,
            double a_weight, std::ostringstream& stream) {
        stream << xml_elementary_material_open << tag << xml_elementary_material_first_inter << tag;
        stream << xml_elementary_material_second_inter << a_number << xml_elementary_material_third_inter;
        stream << a_weight << xml_elementary_material_fourth_inter << density;
        stream << xml_elementary_material_close;
    }
    
    void tk2CMSSW::compositeMaterial(std::string name, double density, CompType method,
            std::vector<std::pair<std::string, double> >& es, std::ostringstream& stream) {
        stream << xml_composite_material_open << name << xml_composite_material_first_inter;
        stream << density << xml_composite_material_second_inter ;
        switch (method) {
            case wt : stream << "mixture by weight";
            break;
            case vl : stream << "mixture by volume";
            break;
            case ap : stream << "mixture by atomic proportion";
            break;
            default: std::cerr << "tk2CMSSW::compositeMaterial(): unknown method identifier for composite material. Using mixture by weight." << std::endl;
            stream << "mixture by weight";
        }
        stream << xml_general_inter;
        for (unsigned int i = 0; i < es.size(); i++) {
            stream << xml_material_fraction_open << es.at(i).second << xml_material_fraction_inter;
            stream << xml_fileident << ":" << es.at(i).first << xml_material_fraction_close;
        }
        stream << xml_composite_material_close;
    }
    
    void tk2CMSSW::logicalPart(std::string name, std::string solid, std::string material, std::ostringstream& stream) {
        stream << xml_logical_part_open << name << xml_logical_part_first_inter << solid;
        stream << xml_logical_part_second_inter << material << xml_logical_part_close;
    }
    
    void tk2CMSSW::box(std::string name, double dx, double dy, double dz, std::ostringstream& stream) {
        stream << xml_box_open << name << xml_box_first_inter << dx << xml_box_second_inter << dy;
        stream << xml_box_third_inter << dz << xml_box_close;
    }
    
    void tk2CMSSW::trapezoid(std::string name, double dx, double dxx, double dy, double dz, std::ostringstream& stream) {
        stream << xml_trapezoid_open << name << xml_trapezoid_first_inter << dx;
        stream << xml_trapezoid_second_inter << dxx << xml_trapezoid_third_inter << dy;
        stream << xml_trapezoid_fourth_inter << dy << xml_trapezoid_fifth_inter << dz;
        stream << xml_trapezoid_close;
    }
    
    void tk2CMSSW::tubs(std::string name, double rmin, double rmax, double dz, std::ostringstream& stream) {
        stream << xml_tubs_open << name << xml_tubs_first_inter << rmin << xml_tubs_second_inter << rmax;
        stream << xml_tubs_third_inter << dz << xml_tubs_close;
    }
    
    void tk2CMSSW::polycone(std::string name, std::vector<std::pair<double, double> >& rzu,
            std::vector<std::pair<double, double> >& rzd, std::ostringstream& stream) {
        stream << xml_polycone_open << name << xml_polycone_inter;
        for (unsigned int i = 0; i < rzu.size(); i++) {
            stream << xml_rzpoint_open << rzu.at(i).first << xml_rzpoint_inter << rzu.at(i).second << xml_rzpoint_close;
        }
        for (unsigned int i = rzd.size(); i > 0; i--) {
            stream << xml_rzpoint_open << rzd.at(i - 1).first << xml_rzpoint_inter << rzd.at(i - 1).second << xml_rzpoint_close;
        }
        stream << xml_polycone_close;
    }
    
    void tk2CMSSW::posPart(std::string parent, std::string child, Rotation& rot, Translation& trans, int copy, std::ostringstream& stream) {
        stream << xml_pos_part_open << copy << xml_pos_part_first_inter << parent;
        stream << xml_pos_part_second_inter << child << xml_general_endline;
        if (!rot.name.empty()) rotation(rot.name, rot.phix, rot.phiy, rot.phiz, rot.thetax, rot.thetay, rot.thetaz, stream);
        if (!(trans.dx == 0.0 && trans.dy == 0.0 && trans.dz == 0.0)) translation(trans.dx, trans.dy, trans.dz, stream);
        stream << xml_pos_part_close;
    }
    
    void tk2CMSSW::rotation(std::string name, double phix, double phiy, double phiz,
            double thetax, double thetay, double thetaz, std::ostringstream& stream) {
        stream << xml_rotation_open << name << xml_rotation_first_inter << phix << xml_rotation_second_inter << phiy;
        stream << xml_rotation_third_inter << phiz << xml_rotation_fourth_inter << thetax << xml_rotation_fifth_inter;
        stream << thetay << xml_rotation_sixth_inter << thetaz << xml_rotation_close;
    }
    
    void tk2CMSSW::translation(double x, double y, double z, std::ostringstream& stream) {
        stream << xml_translation_open << x << xml_translation_first_inter << y << xml_translation_second_inter << z;
        stream << xml_translation_close;
    }
    
    std::string tk2CMSSW::stringParam(std::string name, std::string value) {
        std::string res;
        res = xml_algorithm_string + name + xml_algorithm_value + value + xml_general_endline;
        return res;
    }
    
    std::string tk2CMSSW::numericParam(std::string name, std::string value) {
        std::string res;
        res = xml_algorithm_numeric + name + xml_algorithm_value + value + xml_general_endline;
        return res;
    }
    
    std::string tk2CMSSW::vectorParam(double x, double y, double z) {
        std::ostringstream res;
        res << xml_algorithm_vector_open << x << "," << y << "," << z << "," << xml_algorithm_vector_close;
        return res.str();
    }
    
    void tk2CMSSW::specPar(std::string name, std::pair<std::string, std::string> param,
            std::vector<std::string>& partsel, std::ostringstream& stream) {
        stream << xml_spec_par_open << name << xml_general_inter;
        for (unsigned i = 0; i < partsel.size(); i++) {
            stream << xml_spec_par_selector << partsel.at(i) << xml_general_endline;
        }
        stream << xml_spec_par_parameter_first << param.first << xml_spec_par_parameter_second;
        stream << param.second << xml_spec_par_close;
    }
    
    // private
    void tk2CMSSW::analyse(MaterialTable& mt, MaterialBudget& mb, std::vector<Element>& e,
            std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
            std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::vector<SpecPar>& t) {
        //TODO: finish endcaps and transplant code after TIDF/TIDB volumes to analyseEndcapServices()
        Tracker& tr = mb.getTracker();
        InactiveSurfaces& is = mb.getInactiveSurfaces();
        std::vector<std::vector<ModuleCap> >& bc = mb.getBarrelModuleCaps();
        std::vector<std::vector<ModuleCap> >& ec = mb.getEndcapModuleCaps();
        std::vector<std::vector<ModuleCap> >::iterator oiter, oguard;
        std::vector<ModuleCap>::iterator iiter, iguard;
        //std::vector<InactiveElement>::iterator iter, guard;
        int layer;
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
        SpecPar spec;
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
        // define top-level barrel volume container (polycone)
        shape.type = pc;
        shape.name_tag = xml_tob;
        // create polycone info
        analyseBarrelContainer(tr, shape.rzup, shape.rzdown);
        s.push_back(shape);
        logic.name_tag = shape.name_tag;
        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
        l.push_back(logic);
        pos.parent_tag = xml_fileident + ":" + xml_tracker;
        pos.child_tag = logic.shape_tag;
        p.push_back(pos);
        // define top-level endcap volume containers (polycones)
        /**shape.name_tag = xml_tidf;
         * analyseForwardEndcapContainer(tr, shape.rzup, shape.rzdown);
         * s.push_back(shape);
         * logic.name_tag = shape.name_tag;
         * logic.shape_tag = xml_fileident + ":" + logic.name_tag;
         * l.push_back(logic);
         * pos.parent_tag = xml_fileident + ":" + xml_tracker;
         * pos.child_tag = logic.shape_tag;
         * p.push_back(pos);
         * shape.name_tag = xml_tidb;
         * analyseBackwardEndcapContainer(tr, shape.rzup, shape.rzdown);
         * s.push_back(shape);*/
        /* logic.name_tag = shape.name_tag;
         * logic.shape_tag = xml_fileident + ":" + logic.name_tag;
         * l.push_back(logic);
         * pos.parent_tag = xml_fileident + ":" + xml_tracker;
         * pos.child_tag = logic.shape_tag;
         * p.push_back(pos);*/
        shape.rzup.clear();
        shape.rzdown.clear();
        // translate entries in mt to elementary materials
        analyseElements(mt, e);
        // analyse barrel
        analyseLayers(bc, tr, c, l, s, p, a, t);
        // e_mod: one composite for every ring
        // s and l: one entry for every ring module, one for every ring, one for every disc
        // p: one entry for every disc, one for every ring, one module, wafer and active per ring
        // a: two per ring with modules inside ring
        layer = 1;
        alg.name = xml_ecalgo;
        oguard = ec.end();
        // endcap caps layer loop
        /* for (oiter = ec.begin(); oiter != oguard; oiter++) {
         * std::set<int> ridx;
         * std::map<int, RingInfo> rinfo;
         * double rmin = tr.getEndcapLayers()->at(layer - 1)->getMinRho();
         * double rmax = tr.getEndcapLayers()->at(layer - 1)->getMaxRho();
         * double zmax = tr.getEndcapLayers()->at(layer - 1)->getMaxZ();
         * double zmin = tr.getEndcapLayers()->at(layer - 1)->getMinZ();
         * std::ostringstream dname, pconverter;
         * dname << xml_disc << layer;
         * shape.type = tp;
         * shape.rmin = 0.0;
         * shape.rmax = 0.0;
         * pos.trans.dz = 0.0;
         * iguard = oiter->end();
         * for (iiter = oiter->begin(); iiter != iguard; iiter++) {
         * if (ridx.find(iiter->getModule().getRing()) == ridx.end()) {
         * ridx.insert(iiter->getModule().getRing());
         * std::ostringstream matname, rname, mname, specname;
         * matname << xml_base_actcomp << "D" << layer << "R" << iiter->getModule().getRing();
         * c.push_back(createComposite(matname.str(), compositeDensity(*iiter), *iiter));
         * rname << xml_ring << iiter->getModule().getRing() << dname.str();
         * mname << xml_endcap_module << iiter->getModule().getRing() << dname.str();
         * // collect ring info
         * RingInfo ri;
         * ri.name = rname.str();
         * ri.childname = mname.str();
         * ri.fw = (iiter->getModule().getMeanPoint().Z() < (zmin + zmax) / 2.0);
         * ri.modules = static_cast<EndcapLayer*>(tr.getEndcapLayers()->at(layer - 1))->getModulesOnRing().at(iiter->getModule().getRing() - 1);
         * ri.rin = iiter->getModule().getMinRho();
         * ri.rout = iiter->getModule().getMaxRho();
         * ri.rmid = iiter->getModule().getMeanPoint().Rho();
         * ri.mthk = iiter->getModule().getThickness();
         * ri.phi = iiter->getModule().getMeanPoint().Phi();
         * rinfo.insert(std::pair<int, RingInfo>(iiter->getModule().getRing(), ri));
         * // module trapezoid
         * shape.name_tag = mname.str();
         * shape.dx = static_cast<EndcapModule&>(iiter->getModule()).getWidthLo() / 2.0;
         * shape.dxx = static_cast<EndcapModule&>(iiter->getModule()).getWidthHi() / 2.0;
         * shape.dy = iiter->getModule().getHeight() / 2.0;
         * shape.dz = iiter->getModule().getThickness() / 2.0;
         * s.push_back(shape);
         * logic.name_tag = shape.name_tag;
         * logic.shape_tag = xml_fileident + ":" + logic.name_tag;
         * logic.material_tag = xml_material_air;
         * l.push_back(logic);
         * // wafer
         * pos.parent_tag = logic.shape_tag;
         * shape.name_tag = mname.str() + xml_base_waf;
         * s.push_back(shape);
         * logic.name_tag = shape.name_tag;
         * logic.shape_tag = xml_fileident + ":" + logic.name_tag;
         * logic.material_tag = xml_material_air;
         * l.push_back(logic);
         * pos.child_tag = logic.shape_tag;
         * p.push_back(pos);
         * // active surface
         * pos.parent_tag = logic.shape_tag;
         * shape.name_tag = mname.str() + xml_base_act;
         * s.push_back(shape);
         * logic.name_tag = shape.name_tag;
         * logic.shape_tag = xml_fileident + ":" + logic.name_tag;
         * logic.material_tag = xml_fileident + ":" + matname.str();
         * l.push_back(logic);
         * pos.child_tag = logic.shape_tag;
         * p.push_back(pos);
         * // topology
         * specname << xml_apv_head << (iiter->getModule().getNStripsAcross() / 128) << xml_apv_tail;
         * int id = findSpecParIndex(t, specname.str());
         * if (id >= 0) t.at(id).partselectors.push_back(logic.name_tag);
         * else {
         * spec.partselectors.clear();
         * spec.name = specname.str();
         * spec.parameter.first = xml_apv_number;
         * spec.parameter.second = iiter->getModule().getNStripsAcross() / 128;
         * spec.partselectors.push_back(logic.name_tag);
         * t.push_back(spec);
         * }
         * }
         * }
         * // rings
         * shape.type = tb;
         * shape.dx = 0.0;
         * shape.dxx = 0.0;
         * shape.dy = 0.0;
         * shape.dz = findDeltaZ(tr.getEndcapLayers()->at(layer - 1)->getModuleVector()->begin(),
         * tr.getEndcapLayers()->at(layer - 1)->getModuleVector()->end(), (zmin + zmax) / 2.0) / 2.0;
         * pos.parent_tag = xml_fileident + ":" + dname.str();
         * std::set<int>::const_iterator siter, sguard = ridx.end();
         * for (siter = ridx.begin(); siter != sguard; siter++) {
         * shape.name_tag = rinfo[*siter].name;
         * shape.rmin = rinfo[*siter].rin;
         * shape.rmax = rinfo[*siter].rout;
         * s.push_back(shape);
         * logic.name_tag = shape.name_tag;
         * logic.shape_tag = xml_fileident + ":" + logic.name_tag;
         * logic.material_tag = xml_material_air;
         * l.push_back(logic);
         * pos.child_tag = logic.shape_tag;
         * if (rinfo[*siter].fw) pos.trans.dz = (zmin - zmax) / 2.0 + shape.dz;
         * else pos.trans.dz = (zmax - zmin) / 2.0 - shape.dz;
         * p.push_back(pos);
         * alg.parent = logic.shape_tag;
         * alg.parameters.push_back(stringParam(xml_childparam, rinfo[*siter].childname));
         * pconverter << (rinfo[*siter].modules / 2);
         * alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
         * pconverter.str("");
         * alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
         * alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
         * alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
         * pconverter << rinfo[*siter].phi;
         * alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
         * pconverter.str("");
         * pconverter << rinfo[*siter].rmid;
         * alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
         * pconverter.str("");
         * alg.parameters.push_back(vectorParam(0, 0, (zmin - zmax) / 2.0 + rinfo[*siter].mthk / 2.0));
         * a.push_back(alg);
         * alg.parameters.clear();
         * alg.parameters.push_back(stringParam(xml_childparam, rinfo[*siter].childname));
         * pconverter << (rinfo[*siter].modules / 2);
         * alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
         * pconverter.str("");
         * alg.parameters.push_back(numericParam(xml_startcopyno, "2"));
         * alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
         * alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
         * pconverter << (rinfo[*siter].phi + 2 * PI / (double)(rinfo[*siter].modules));
         * alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
         * pconverter.str("");
         * pconverter << rinfo[*siter].rmid;
         * alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
         * pconverter.str("");
         * alg.parameters.push_back(vectorParam(0, 0, (zmax - zmin) / 2.0 - rinfo[*siter].mthk / 2.0));
         * a.push_back(alg);
         * alg.parameters.clear();
         * }
         * pos.trans.dz = 0.0;
         * //disc
         * shape.name_tag = dname.str();
         * shape.rmin = rmin;
         * shape.rmax = rmax;
         * shape.dz = (zmax - zmin) / 2.0;
         * s.push_back(shape);
         * logic.name_tag = shape.name_tag;
         * logic.shape_tag = xml_fileident + ":" + logic.name_tag;
         * logic.material_tag = xml_material_air;
         * l.push_back(logic);
         * pos.parent_tag = xml_fileident + ":" + xml_tracker;
         * pos.child_tag = logic.shape_tag;
         * pos.trans.dz = (zmax + zmin) / 2.0;
         * p.push_back(pos);
         * layer++;
         * }*/
        // barrel services
        analyseBarrelServices(is, c, l, s, p, t);
        // endcap services
        analyseEndcapServices(is, c, l, s, p, t);
        // supports
        analyseSupports(is, c, l, s, p, t);
    }
    
    void tk2CMSSW::analyseElements(MaterialTable&mattab, std::vector<Element>& elems) {
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
    
    void tk2CMSSW::analyseBarrelContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
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
    
    void tk2CMSSW::analyseForwardEndcapContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
            std::vector<std::pair<double, double> >& down) {
        //TODO
        //get endcap layers
        //find first and last forward disc (iterators start = discs->begin, stop = start => condition)
        // call analyseEndcapContainer() with those iterators and pts
    }
    
    void tk2CMSSW::analyseBackwardEndcapContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
            std::vector<std::pair<double, double> >& down) {
        //TODO
        //get endcap layers
        //find first and last backward disc (iterators start = discs->begin => condition, stop = discs->end)
        // call analyseEndcapContainer() with those iterators and pts
    }
    
    void tk2CMSSW::analyseEndcapContainer(std::vector<Module*>::iterator i, std::vector<Module*>::iterator g,
            std::vector<std::pair<double, double> >& up, std::vector<std::pair<double, double> >& down) {
        if (i != g) {
            up.clear();
            down.clear();
            //fill up pts according to notes using iterators as boundaries
        }
    }
    
    void tk2CMSSW::analyseLayers(std::vector<std::vector<ModuleCap> >& bc, Tracker& tr, std::vector<Composite>& c,
            std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<PosInfo>& p,
            std::vector<AlgoInfo>& a, std::vector<SpecPar>& t) {
        int layer;
        std::vector<std::vector<ModuleCap> >::iterator oiter, oguard;
        std::vector<ModuleCap>::iterator iiter, iguard;
        // container inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        AlgoInfo alg;
        SpecPar spec;
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
                        std::vector<ModuleCap>::iterator partner;
                        std::ostringstream matname, shapename, specname;
                        // module composite material
                        matname << xml_base_actcomp << "L" << layer << "P" << iiter->getModule().getRing();
                        c.push_back(createComposite(matname.str(), compositeDensity(*iiter), *iiter));
                        // module box
                        shapename << iiter->getModule().getRing() << lname.str();
                        shape.name_tag = xml_barrel_module + shapename.str();
                        shape.dx = iiter->getModule().getThickness() / 2.0;
                        shape.dy = iiter->getModule().getArea() / iiter->getModule().getHeight() / 2.0;
                        shape.dz = iiter->getModule().getHeight() / 2.0;
                        s.push_back(shape);
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                        logic.material_tag = xml_material_air;
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
                        s.push_back(shape);
                        pos.parent_tag = logic.shape_tag;
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                        l.push_back(logic);
                        pos.child_tag = logic.shape_tag;
                        pos.copy = 1;
                        pos.trans.dx = 0.0;
                        pos.trans.dz = 0.0;
                        p.push_back(pos);
                        // active surface
                        shape.name_tag = xml_barrel_module + shapename.str() + xml_base_act;
                        s.push_back(shape);
                        pos.parent_tag = logic.shape_tag;
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                        logic.material_tag = xml_fileident + ":" + matname.str();
                        l.push_back(logic);
                        pos.child_tag = logic.shape_tag;
                        p.push_back(pos);
                        // topology
                        specname << xml_apv_head << (iiter->getModule().getNStripsAcross() / 128) << xml_apv_tail;
                        int id = findSpecParIndex(t, specname.str());
                        if (id >= 0) t.at(id).partselectors.push_back(logic.name_tag);
                        else {
                            spec.partselectors.clear();
                            spec.name = specname.str();
                            spec.parameter.first = xml_apv_number;
                            specname.str("");
                            specname << (iiter->getModule().getNStripsAcross() / 128);
                            spec.parameter.second = specname.str();
                            spec.partselectors.push_back(logic.name_tag);
                            t.push_back(spec);
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
                pconverter << logic.shape_tag;
                if (is_short) {
                    shape.name_tag = rname.str() + xml_minus;
                    s.push_back(shape);
                    logic.name_tag = shape.name_tag;
                    logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                    l.push_back(logic);
                }
                ds = fromRim(rmax, shape.dy);
                // layer
                shape.type = tb;
                shape.dx = 0.0;
                shape.dy = 0.0;
                pos.trans.dx = 0.0;
                pos.trans.dz = 0.0;
                shape.name_tag = lname.str();
                if (is_short) shape.name_tag = shape.name_tag + xml_plus;
                shape.rmin = rmin;
                shape.rmax = rmax;
                if (is_short) shape.dz = (zmax - zmin) / 2.0;
                else shape.dz = zmax;
                s.push_back(shape);
                logic.name_tag = shape.name_tag;
                logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                l.push_back(logic);
                pos.parent_tag = xml_fileident + ":" + xml_tob;
                //pos.parent_tag = xml_fileident + ":" + xml_tracker;
                pos.child_tag = logic.shape_tag;
                if (is_short) pos.trans.dz = zmin + (zmax - zmin) / 2.0;
                p.push_back(pos);
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
                alg.parameters.push_back(numericParam(xml_zposition, "0.0*mm"));
                pconverter << static_cast<BarrelLayer*>(tr.getBarrelLayers()->at(layer - 1))->getRods();
                alg.parameters.push_back(numericParam(xml_number, pconverter.str()));
                alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
                alg.parameters.push_back(numericParam(xml_incrcopyno, "1"));
                a.push_back(alg);
                if (is_short) {
                    shape.name_tag = lname.str() + xml_minus;
                    s.push_back(shape);
                    logic.name_tag = shape.name_tag;
                    logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                    l.push_back(logic);
                    pos.child_tag = logic.shape_tag;
                    pos.trans.dz = -(zmin + (zmax - zmin) / 2.0);
                    p.push_back(pos);
                    alg.parent = logic.shape_tag;
                    pconverter.str("");
                    pconverter << xml_fileident << ":" << rname.str() << xml_minus;
                    alg.parameters.front() = stringParam(xml_childparam, pconverter.str());
                    a.push_back(alg);
                }
                alg.parameters.clear();
            }
            layer++;
        }
    }
    
    //analyseEndcaps()
    
    void tk2CMSSW::analyseBarrelServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
            std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<SpecPar>& t) {
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
    
    void tk2CMSSW::analyseEndcapServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
            std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<SpecPar>& t) {
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
    
    void tk2CMSSW::analyseSupports(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
            std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<SpecPar>& t) {
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
            if (iter->getCategory() != MaterialProperties::e_sup) { //TODO: remove this to enable endcap supports
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
    }
    
    tk2CMSSW::Composite tk2CMSSW::createComposite(std::string name, double density, MaterialProperties& mp) {
        Composite comp;
        comp.name = name;
        comp.density = density;
        comp.method = wt;
        for (unsigned int i = 0; i < mp.localMassCount(); i++) {
            std::pair<std::string, double> p;
            p.first = mp.getLocalTag(i);
            p.second = mp.getLocalMass(i);
            comp.elements.push_back(p);
        }
        for (unsigned int i = 0; i < mp.exitingMassCount(); i++) {
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
        }
        for (unsigned int i = 0; i < comp.elements.size(); i++)
            comp.elements.at(i).second = comp.elements.at(i).second / mp.getTotalMass();
        return comp;
    }
    
    std::vector<ModuleCap>::iterator tk2CMSSW::findPartnerModule(std::vector<ModuleCap>::iterator i,
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
    
    double tk2CMSSW::findDeltaR(std::vector<Module*>::iterator start,
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
        dr = (*mod1)->getMinRho() - (*mod2)->getMinRho() + (*mod1)->getThickness();
        return dr;
    }
    
    double tk2CMSSW::findDeltaZ(std::vector<Module*>::iterator start,
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
    
    int tk2CMSSW::findSpecParIndex(std::vector<SpecPar>& specs, std::string label) {
        int idx = 0, size = (int)(specs.size());
        while (idx < size) {
            if (specs.at(idx).name.compare(label) == 0) return idx;
            idx++;
        }
        return -1;
    }
    
    double tk2CMSSW::fromRim(double r, double w) {
        double s = asin(w / r);
        s = 1 - cos(s);
        s = s * r;
        return s;
    }
    
    double tk2CMSSW::compositeDensity(InactiveElement& ie) {
        double d = ie.getRWidth() + ie.getInnerRadius();
        d = d * d - ie.getInnerRadius() * ie.getInnerRadius();
        d = 1000 * ie.getTotalMass() / (PI * ie.getZLength() * d);
        return d;
    }
    
    double tk2CMSSW::compositeDensity(ModuleCap& mc) {
        double d = mc.getSurface() * mc.getModule().getThickness();
        d = 1000 * mc.getTotalMass() / d;
        return d;
    }
    
    int tk2CMSSW::Z(double x0, double A) {
        double d = 4 - 4 * (1.0 - 181.0 * A / x0);
        if (d > 0) return floor((sqrt(d) - 2.0) / 2.0 + 0.5);
        else return -1;
    }
    
    void tk2CMSSW::print() {
        std::cout << "tm2CMSSW internal status:" << std::endl;
        std::cout << "elements: " << elements.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < elements.size(); i++) {
            std::cout << "entry " << i << ": tag = " << elements.at(i).tag << ", density = " << elements.at(i).density << ", atomic number = ";
            std::cout << elements.at(i).atomic_number << ", atomic weight = " << elements.at(i).atomic_weight << std::endl;
        }
        std::cout << "composites: " << composites.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < composites.size(); i++) {
            std::cout << "entry " << i << ": name = " << composites.at(i).name << ", density = " << composites.at(i).density << ", method = ";
            switch (composites.at(i).method) {
                case wt: std::cout << "fraction by weight";
                break;
                case vl: std::cout << "fraction by volume";
                break;
                case ap: std::cout << "fraction by atomic proportion";
                break;
                default: std::cout << "unknown method";
            }
            std::cout << std::endl << "elements: ";
            std::vector<std::pair<std::string, double> >& elems = composites.at(i).elements;
            for (unsigned int j = 0; j < elems.size(); j++) std::cout << "(" << elems.at(j).first << ", " << elems.at(j).second << ") ";
            std::cout << std::endl;
        }
        std::cout << "logic: " << logic.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < logic.size(); i++) {
            std::cout << "name_tag = " << logic.at(i).name_tag << ", shape_tag = " << logic.at(i).shape_tag;
            std::cout << ", material_tag = " << logic.at(i).material_tag << std::endl;
        }
        std::cout << "shapes: " << shapes.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < shapes.size(); i++) {
            std::cout << "name_tag = " << shapes.at(i).name_tag << ", type = ";
            switch (shapes.at(i).type) {
                case bx: std::cout << "box, dx = " << shapes.at(i).dx << ", dy = " << shapes.at(i).dy << ", dz = ";
                std::cout << shapes.at(i).dz;
                break;
                case tb: std::cout << "tube, rmin = " << shapes.at(i).rmin << ", rmax = " << shapes.at(i).rmax;
                std::cout << ", dz = " << shapes.at(i).dz;
                break;
                case tp: std::cout << "trapezoid, dx = " << shapes.at(i).dx << ", dxx = " << shapes.at(i).dxx;
                std::cout << ", dy = " << shapes.at(i).dy << ", dz = " << shapes.at(i).dz;
                break;
                default: std::cout << "unknown shape";
            }
            std::cout << std::endl;
        }
        std::cout << "positions: " << positions.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < positions.size(); i++) {
            std::cout << "parent_tag = " << positions.at(i).parent_tag << ", child_tag = " << positions.at(i).child_tag;
            std::cout << ", rotation = (" << (positions.at(i).rot.name.empty() ? "[no name]": positions.at(i).rot.name) << ", ";
            std::cout  << positions.at(i).rot.phix << ", " << positions.at(i).rot.phiy << ", " << positions.at(i).rot.phiz << ", ";
            std::cout << positions.at(i).rot.thetax << ", " << positions.at(i).rot.thetay << ", " << positions.at(i).rot.thetaz;
            std::cout << "), translation = (" << positions.at(i).trans.dx << ", " << positions.at(i).trans.dy << ", ";
            std::cout << positions.at(i).trans.dz << ")" << std::endl;
        }
    }
}
