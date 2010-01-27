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
        // analyse tracker system and build up collection of elements, composites, hierarchy, shapes, positions, algorithms and topology
        analyse(mt, mb, elements, composites, logic, shapes, positions, algos, specs);
        // translate collected information to XML and write to buffers
        std::ostringstream gbuffer, tbuffer, pbuffer, sbuffer, mbuffer;
        gbuffer << xml_preamble << xml_const_section;
        tbuffer << xml_preamble;
        materialSection(xml_trackerfile, elements, composites, gbuffer);
        logicalPartSection(logic, xml_trackerfile, gbuffer);
        solidSection(shapes, xml_trackerfile, gbuffer);
        posPartSection(positions, algos, xml_trackerfile, gbuffer);
        specParSection(specs, xml_specpars_label, tbuffer);
        prodcuts(specs, pbuffer);
        trackersens(specs, sbuffer);
        recomaterial(specs, mbuffer);
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
        std::ofstream poutstream((outpath + xml_prodcutsfile).c_str());
        poutstream << pbuffer.str() << std::endl;
        poutstream.close();
        std::cout << "CMSSW prodcuts output has been written to " << outpath << xml_prodcutsfile << std::endl;
        std::ofstream soutstream((outpath + xml_trackersensfile).c_str());
        soutstream << sbuffer.str() << std::endl;
        soutstream.close();
        std::cout << "CMSSW sensor surface output has been written to " << outpath << xml_trackersensfile << std::endl;
        std::ofstream moutstream((outpath + xml_recomatfile).c_str());
        moutstream << mbuffer.str() << std::endl;
        moutstream.close();
        std::cout << "CMSSW reco material output has been written to " << outpath << xml_recomatfile << std::endl;
    }
    
    // protected
    void tk2CMSSW::prodcuts(std::vector<SpecParInfo>& t, std::ostringstream& stream) {
        stream << xml_preamble << xml_prodcuts_open;
        for (unsigned int i = 0; i < t.size(); i++) {
            if ((t.at(i).name.substr(0, xml_subdet_tobdet.size()).compare(xml_subdet_tobdet) == 0)
                    || (t.at(i).name.substr(0, xml_subdet_tiddet.size()).compare(xml_subdet_tiddet) == 0)) {
                for (unsigned int j = 0; j < t.at(i).partselectors.size(); j++) {
                    stream << xml_spec_par_selector << t.at(i).partselectors.at(j) << xml_general_endline;
                }
            }
        }
        stream << xml_prodcuts_close << xml_spec_par_close << xml_spec_par_section_close << xml_defclose;
    }
    
    void tk2CMSSW::trackersens(std::vector<SpecParInfo>& t, std::ostringstream& stream) {
        stream << xml_preamble << xml_trackersens_open;
        for (unsigned int i = 0; i < t.size(); i++) {
            if (t.at(i).name.compare(xml_subdet_tobdet + xml_par_tail) == 0) {
                for (unsigned int j = 0; j < t.at(i).partselectors.size(); j++) {
                    stream << xml_spec_par_selector << t.at(i).partselectors.at(j) << xml_general_endline;
                }
            }
        }
        stream << xml_trackersens_endtob << xml_spec_par_close;
        if (endcapsInTopology(t)) {
            stream << xml_trackersens_inter;
            for (unsigned int i = 0; i < t.size(); i++) {
                if (t.at(i).name.compare(xml_subdet_tiddet + xml_par_tail) == 0) {
                    for (unsigned int j = 0; j < t.at(i).partselectors.size(); j++) {
                        stream << xml_spec_par_selector << t.at(i).partselectors.at(j) << xml_general_endline;
                    }
                }
            }
            stream << xml_trackersens_endtid << xml_spec_par_close;
        }
        stream << xml_spec_par_section_close << xml_defclose;
    }
    
    void tk2CMSSW::recomaterial(std::vector<SpecParInfo>& t, std::ostringstream& stream) {
        std::vector<std::pair<std::string, std::vector<std::string> > > b;
        b = buildPaths(t, b);
        if (!b.empty()) {
            std::vector<std::pair<std::string, std::vector<std::string> > >::iterator iter, guard = b.end();
            stream << xml_preamble << xml_spec_par_section_open << xml_specpars_label << xml_general_inter;
            for (iter = b.begin(); iter != guard; iter++) {
                if ((!iter->first.empty()) && (!iter->second.empty())) {
                    std::vector<std::string>::iterator iiter, iguard = iter->second.end();
                    stream << xml_spec_par_open << iter->first << xml_eval_true;
                    for (iiter = iter->second.begin(); iiter != iguard; iiter++) {
                        stream << xml_spec_par_selector << *iiter << xml_general_endline;
                    }
                    stream << xml_recomat_parameters << xml_spec_par_close;
                }
            }
            stream << xml_spec_par_section_close << xml_defclose;
        }
    }
    
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
    
    void tk2CMSSW::specParSection(std::vector<SpecParInfo>& t, std::string label, std::ostringstream& stream) {
        std::vector<SpecParInfo>::iterator titer, tguard = t.end();
        stream << xml_spec_par_section_open << label << xml_general_inter;
        for (titer = t.begin(); titer != tguard; titer++) specPar(titer->name, titer->parameter, titer->partselectors, stream);
        stream << xml_spec_par_section_close;
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
            std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::vector<SpecParInfo>& t) {
        //TODO: finish endcaps
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
        shape.rzup.clear();
        shape.rzdown.clear();
        spec.partselectors.clear();
        // define top-level endcap volume containers (polycones)
        /*shape.name_tag = xml_tidf;
        analyseForwardEndcapContainer(tr, shape.rzup, shape.rzdown);//TODO
        s.push_back(shape);
        logic.name_tag = shape.name_tag;
        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
        l.push_back(logic);
        pos.parent_tag = xml_fileident + ":" + xml_tracker;
        pos.child_tag = logic.shape_tag;
        p.push_back(pos);
        shape.name_tag = xml_tidb;
        analyseBackwardEndcapContainer(tr, shape.rzup, shape.rzdown);//TODO
        s.push_back(shape);
        logic.name_tag = shape.name_tag;
        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
        l.push_back(logic);
        pos.parent_tag = xml_fileident + ":" + xml_tracker;
        pos.child_tag = logic.shape_tag;
        p.push_back(pos);
        spec.name = xml_tid_subdet + xml_par_tail;
        spec.parameter.first = xml_tkddd_structure;
        spec.parameter.second = xml_tid;
        spec.partselectors.push_back(xml_tidb);
        spec.partselectors.push_back(xml_tidf);
        t.push_back(spec);
        shape.rzup.clear();
        shape.rzdown.clear();
        spec.partselectors.clear();*/
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
    
    void tk2CMSSW::analyseLayers(MaterialTable& mt, std::vector<std::vector<ModuleCap> >& bc, Tracker& tr,
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
    
    void tk2CMSSW::analyseDiscs(MaterialTable& mt, std::vector<std::vector<ModuleCap> >& ec, Tracker& tr,
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
        SpecParInfo dspec, rspec, mspec;
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
         * ri.mthk = iiter->getModule().getModuleThickness();
         * ri.phi = iiter->getModule().getMeanPoint().Phi();
         * rinfo.insert(std::pair<int, RingInfo>(iiter->getModule().getRing(), ri));
         * // module trapezoid
         * shape.name_tag = mname.str();
         * shape.dx = static_cast<EndcapModule&>(iiter->getModule()).getWidthLo() / 2.0;
         * shape.dxx = static_cast<EndcapModule&>(iiter->getModule()).getWidthHi() / 2.0;
         * shape.dy = iiter->getModule().getHeight() / 2.0;
         * shape.dz = iiter->getModule().getModuleThickness() / 2.0;
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
         * specname << xml_apv_head << (iiter->getModule().getNStripsAcross() / 128) << xml_par_tail;
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
        if (!dspec.partselectors.empty()) t.push_back(dspec);
        if (!rspec.partselectors.empty()) t.push_back(rspec);
        if (!mspec.partselectors.empty()) t.push_back(mspec);
    }
    
    void tk2CMSSW::analyseBarrelServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
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
    
    void tk2CMSSW::analyseEndcapServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
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
    
    void tk2CMSSW::analyseSupports(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
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
    
    tk2CMSSW::Composite tk2CMSSW::createComposite(std::string name, double density, MaterialProperties& mp, bool nosensors) {
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
        dr = (*mod1)->getMinRho() - (*mod2)->getMinRho() + (*mod1)->getModuleThickness();
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
    
    int tk2CMSSW::findSpecParIndex(std::vector<SpecParInfo>& specs, std::string label) {
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
    
    double tk2CMSSW::compositeDensity(ModuleCap& mc, bool nosensors) {
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
    
    double tk2CMSSW::calculateSensorThickness(ModuleCap& mc, MaterialTable& mt) {
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
    
    int tk2CMSSW::Z(double x0, double A) {
        double d = 4 - 4 * (1.0 - 181.0 * A / x0);
        if (d > 0) return floor((sqrt(d) - 2.0) / 2.0 + 0.5);
        else return -1;
    }
    
    bool tk2CMSSW::endcapsInTopology(std::vector<SpecParInfo>& specs) {
        for (unsigned int i = 0; i < specs.size(); i++) {
            if (specs.at(i).name.compare(xml_subdet_tiddet + xml_par_tail) == 0) return true;
        }
        return false;
    }
    
    std::vector<std::pair<std::string, std::vector<std::string> > >& tk2CMSSW::buildPaths(std::vector<SpecParInfo>& specs,
            std::vector<std::pair<std::string, std::vector<std::string> > >& blocks) {
        std::vector<std::pair<std::string, std::vector<std::string> > >::iterator existing;
        int rrange, mrange, rindex, mindex, ioffset = 0;
        std::string prefix, postfix, spname;
        std::vector<std::string> paths;
        std::vector<Node> tree;
        Node n, nn;
        blocks.clear();
        // TOB
        // build tree: root
        n.name = xml_tob;
        rindex = findEntry(specs, xml_subdet_rod + xml_par_tail);
        if (rindex >= 0) rrange = specs.at(rindex).partselectors.size();
        else rrange = 0;
        // TOBSubDetRodPar loop
        for (int i = 0; i < rrange; i++) n.indices.push_back(i + 1);
        tree.push_back(n);
        for (unsigned int i = 0; i < n.indices.size(); i++) tree.push_back(nn);
        n.indices.clear();
        // build tree: layers and modules
        mindex = findEntry(specs, xml_subdet_tobdet + xml_par_tail);
        if (mindex >= 0) mrange = specs.at(mindex).partselectors.size();
        else mrange = 0;
        if (mindex >= 0) {
            ioffset = rrange + 1;
            for (int i = 1; i <= rrange; i++) {
                tree.at(i).name = specs.at(rindex).partselectors.at(i - 1);
                // TOBSubDetDetPar loop
                for (int j = 0; j < mrange; j++) {
                    int stop;
                    // extract layer number from both n.name and specs.at(mindex).partselectors.at(j)
                    std::string nrod = tree.at(i).name;
                    if ((nrod.size() > xml_plus.size()) &&
                            (nrod.substr(nrod.size() - xml_plus.size()).compare(xml_plus) == 0))
                        nrod = nrod.substr(0, nrod.size() - xml_plus.size());
                    if ((nrod.size() > xml_minus.size()) &&
                            (nrod.substr(nrod.size() - xml_minus.size()).compare(xml_minus) == 0))
                        nrod = nrod.substr(0, nrod.size() - xml_minus.size());
                    nrod = nrod.substr(xml_rod.size());
                    n.name = specs.at(mindex).partselectors.at(j);
                    std::string nmodule = n.name.substr(0, n.name.size() - xml_base_act.size());
                    nmodule = nmodule.substr(xml_barrel_module.size());
                    stop = findAlphaNumericPrefixSize(nmodule);
                    nmodule = nmodule.substr(stop + xml_layer.size());
                    // process a match
                    if (nrod.compare(nmodule) == 0) {
                        int nodenr = findNodeNamed(tree, n.name);
                        if (nodenr == -1) {
                            tree.at(i).indices.push_back(ioffset);
                            tree.push_back(n);
                            ioffset++;
                        }
                        else tree.at(i).indices.push_back(nodenr);
                    }
                }
            }
        }
        // traverse tree: tree.at(0).indices (children of root element TOB) loop
        for (unsigned int i = 0; i < tree.at(0).indices.size(); i++) {
            std::string number, plusminus;
            int stop, child, leaf;
            child = tree.at(0).indices.at(i);
            if ((tree.at(child).name.size() > xml_plus.size())
                    && (tree.at(child).name.substr(tree.at(child).name.size() - xml_plus.size()).compare(xml_plus) == 0))
                plusminus = tree.at(child).name.substr(tree.at(child).name.size() - xml_plus.size());
            if ((tree.at(child).name.size() > xml_minus.size())
                    && (tree.at(child).name.substr(tree.at(child).name.size() - xml_minus.size()).compare(xml_minus) == 0))
                plusminus = tree.at(child).name.substr(tree.at(child).name.size() - xml_minus.size());
            number = tree.at(child).name.substr(xml_rod.size());
            number = number.substr(0, number.size() - plusminus.size());
            spname = xml_tob_prefix + xml_layer + number;
            prefix = xml_tob + "/" + xml_layer + number + "/" + tree.at(child).name;
            // traverse tree: tree.at(i).indices (children of rods) loop
            for (unsigned int j = 0; j < tree.at(child).indices.size(); j++) {
                leaf = tree.at(child).indices.at(j);
                stop = tree.at(leaf).name.size() - xml_base_act.size();
                postfix = "/" + tree.at(leaf).name.substr(0, stop) + "/" + tree.at(leaf).name.substr(0, stop) + xml_base_waf + "/" + tree.at(leaf).name;
                paths.push_back(prefix + postfix);
            }
            existing = findEntry(spname, blocks);
            if (existing != blocks.end()) existing->second.insert(existing->second.end(), paths.begin(), paths.end());
            else blocks.push_back(std::pair<std::string, std::vector<std::string> >(spname, paths));
            paths.clear();
        }
        //TODO: endcaps
        return blocks;
    }
    
    int tk2CMSSW::findAlphaNumericPrefixSize(std::string s) {
        std::ostringstream st;
        int number;
        if (!s.empty()) {
            number = atoi(s.c_str());
            if (number == 0) return number;
            else {
                st << number;
                return st.str().size();
            }
        }
        return 0;
    }
    
    int tk2CMSSW::findNodeNamed(std::vector<Node>& tree, std::string name) {
        for (unsigned int i = 0; i < tree.size(); i++) {
            if (tree.at(i).name.compare(name) == 0) return i;
        }
        return -1;
    }
    
    int tk2CMSSW::findEntry(std::vector<SpecParInfo>& specs, std::string name) {
        int index = 0;
        while (index < (int)(specs.size())) {
            if (specs.at(index).name.compare(name) == 0) return index;
            index++;
        }
        return -1;
    }
    
    std::vector<std::pair<std::string, std::vector<std::string> > >::iterator tk2CMSSW::findEntry(std::string name,
            std::vector<std::pair<std::string, std::vector<std::string> > >& data) {
        std::vector<std::pair<std::string, std::vector<std::string> > >::iterator result = data.begin();
        std::vector<std::pair<std::string, std::vector<std::string> > >::iterator guard = data.end();
        while (result != guard) {
            if ((result->first).compare(name) == 0) return result;
            result++;
        }
        return result;
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
