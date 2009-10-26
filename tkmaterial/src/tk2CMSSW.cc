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
        // analyse tracker system and build up collection of elements, composites, hierarchy, shapes and position
        analyse(mt, mb, elements, composites, logic, shapes, positions, algos);
        // translate collected information to XML and write to buffer
        std::ostringstream buffer;
        buffer << xml_preamble << xml_const_section;
        materialSection(xml_trackerfile, elements, composites, buffer);
        logicalPartSection(logic, xml_trackerfile, buffer);
        solidSection(shapes, xml_trackerfile, buffer);
        posPartSection(positions, algos, xml_trackerfile, buffer);
        buffer << xml_defclose;
        // write contents of buffer to top-level file
        bfs::remove_all(outpath.c_str());
        bfs::create_directory(outpath);
        std::ofstream outstream((outpath + xml_trackerfile).c_str());
        outstream << buffer.str() << std::endl;
        outstream.close();
        std::cout << "CMSSW XML output has been written to " << outpath << xml_trackerfile << std::endl;
    }
    
    // protected
    void tk2CMSSW::materialSection(std::string name , std::vector<Element>& e, std::vector<Composite>& c, std::ostringstream& stream) {
        stream << xml_material_section_open << name << xml_material_section_inter;
        for (unsigned int i = 0; i < e.size(); i++) elementaryMaterial(e.at(i).tag, e.at(i).density, e.at(i).atomic_number, e.at(i).atomic_weight, stream);
        for (unsigned int i = 0; i < c.size(); i++) compositeMaterial(c.at(i).name, c.at(i).density, c.at(i).method, c.at(i).elements, stream);
        stream << xml_material_section_close;
    }
    
    void tk2CMSSW::logicalPartSection(std::vector<LogicalInfo>& l, std::string label, std::ostringstream& stream) {
        std::vector<LogicalInfo>::const_iterator iter, guard = l.end();
        stream << xml_logical_part_section_open << label << xml_logical_part_section_inter;
        for (iter = l.begin(); iter != guard; iter++) logicalPart(iter->name_tag, iter->shape_tag, iter->material_tag, stream);
        stream << xml_logical_part_section_close;
    }
    
    void tk2CMSSW::solidSection(std::vector<ShapeInfo>& s, std::string label, std::ostringstream& stream) {
        std::vector<ShapeInfo>::const_iterator iter, guard = s.end();
        stream << xml_solid_section_open << label << xml_solid_section_inter;
        for (iter = s.begin(); iter != guard; iter++) {
            switch (iter->type) {
                case bx : box(iter->name_tag, iter->dx, iter->dy, iter->dz, stream);
                break;
                case tp : trapezoid(iter->name_tag, iter->dx, iter->dxx, iter->dy, iter->dz, stream);
                break;
                case tb : tubs(iter->name_tag, iter->rmin, iter->rmax, iter->dz, stream);
                break;
                default: std::cerr << "solidSection(): unknown shape type found. Using box." << std::endl;
                box(iter->name_tag, iter->dx, iter->dy, iter->dz, stream);
            }
        }
        stream << xml_solid_section_close;
    }
    
    void tk2CMSSW::posPartSection(std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::string label, std::ostringstream& stream) {
        std::vector<PosInfo>::iterator piter, pguard = p.end();
        std::vector<AlgoInfo>::iterator aiter, aguard = a.end();
        stream << xml_pos_part_section_open << label << xml_pos_part_section_inter;
        for (piter = p.begin(); piter != pguard; piter++) posPart(piter->parent_tag, piter->child_tag, piter->rot, piter->trans, piter->copy, stream);
        for (aiter = a.begin(); aiter != aguard; aiter++) algorithm(aiter->name, aiter->parent, aiter->parameters, stream);
        stream << xml_pos_part_section_close;
    }
    
    void tk2CMSSW::algorithm(std::string name, std::string parent,
            std::vector<std::string>& params, std::ostringstream& stream) {
        stream << xml_algorithm_open << name << xml_algorithm_parent << parent << xml_algorithm_endline;
        for (unsigned int i = 0; i < params.size(); i++) stream << params.at(i);
        stream << xml_algorithm_close;
    }
    
    std::string tk2CMSSW::stringParam(std::string name, std::string value) {
        std::string res;
        res = xml_algorithm_string + name + xml_algorithm_value + value + xml_algorithm_endline;
        return res;
    }
    
    std::string tk2CMSSW::numericParam(std::string name, std::string value) {
        std::string res;
        res = xml_algorithm_numeric + name + xml_algorithm_value + value + xml_algorithm_endline;
        return res;
    }
    
    //                                    Vector: os << xml_algorithm_vector_open << name << xml_algorithm_vector_first_inter << type;
    //                                                os << xml_algorithm_vector_second_inter << entries << xml_algorithm_vector_third_inter;
    //                                                if (!vector.empty()) os << vector.at(0);
    //                                                for (unsigned int i = 1; i < vector.size(); i++) os << "," << vector.at(i);
    //                                                os << xml_algorithm_vector_close; (see xmltest.cpp for code)
    
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
        stream << xml_composite_material_third_inter;
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
    
    void tk2CMSSW::posPart(std::string parent, std::string child, Rotation& rot, Translation& trans, int copy, std::ostringstream& stream) {
        stream << xml_pos_part_open << copy << xml_pos_part_first_inter << parent;
        stream << xml_pos_part_second_inter << child << xml_pos_part_third_inter;
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
    
    // private
    void tk2CMSSW::analyse(MaterialTable& mt, MaterialBudget& mb, std::vector<Element>& e, std::vector<Composite>& c,
            std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<AlgoInfo>& a) {
        Tracker& tr = mb.getTracker();
        InactiveSurfaces& is = mb.getInactiveSurfaces();
        std::vector<std::vector<ModuleCap> >& bc = mb.getBarrelModuleCaps();
        std::vector<std::vector<ModuleCap> >& ec = mb.getEndcapModuleCaps();
        std::vector<std::vector<ModuleCap> >::iterator oiter, oguard;
        std::vector<ModuleCap>::iterator iiter, iguard;
        std::vector<InactiveElement>::iterator iter, guard;
        int counter, layer;
        e.clear();
        c.clear();
        l.clear();
        s.clear();
        p.clear();
        // inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        AlgoInfo alg;
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
        shape.rmin = 0.0;
        shape.rmax = outer_radius + volume_width;
        s.push_back(shape);
        logic.name_tag = xml_tracker;
        logic.shape_tag = xml_fileident + ":" + xml_tracker;
        logic.material_tag = xml_material_air;
        l.push_back(logic);
        // translate entries in mt to elementary materials
        for (unsigned int i = 0; i < mt.rowCount(); i++) {
            Element elem;
            MaterialRow& r = mt.getMaterial(i);
            elem.tag = r.tag;
            elem.density = r.density;
            elem.atomic_weight = pow((r.ilength / 35.), 3); // magic!
            elem.atomic_number = Z(r.rlength, elem.atomic_weight);
            e.push_back(elem);
        }
        // b_mod: one composite for every module position on rod
        // s and l: one entry for every module position on rod (box), one for every layer (tube), rods TBD
        // p: one entry for every layer, first rod+ and rod- of that layer, wafer and active of one +/-rod, module of rod+ and module of rod-
        //TODO: handle mezzanine layers, this one only manages regular layers for now!!!
        layer = 1;
        alg.name = xml_tobalgo;
        oguard = bc.end();
        for (oiter = bc.begin(); oiter != oguard; oiter++) {
            double rmin = tr.getBarrelLayers()->at(layer - 1)->getMinRho();
            double rmax = tr.getBarrelLayers()->at(layer - 1)->getMaxRho();
            double zmin = tr.getBarrelLayers()->at(layer - 1)->getMinZ();
            double zmax = tr.getBarrelLayers()->at(layer - 1)->getMaxZ();
            double deltar = findDeltaR(tr.getBarrelLayers()->at(layer - 1)->getModuleVector()->begin(),
                    tr.getBarrelLayers()->at(layer - 1)->getModuleVector()->end(), (rmin + rmax) / 2.0);
            double umin = findMinRhoUpperRod(tr.getBarrelLayers()->at(layer - 1)->getModuleVector()->begin(),
                    tr.getBarrelLayers()->at(layer - 1)->getModuleVector()->end(), (rmin + rmax) / 2.0, rmax);
            if ((deltar == 0.0) || (umin == 0.0)) continue;
            bool is_short = (zmax < 0.0) || (zmin > 0.0);
            bool is_relevant = !is_short || (zmin > 0.0);
            if (is_relevant) {
                std::ostringstream lname, rname, pconverter;
                lname << xml_layer << layer;
                rname << xml_rod << layer;
                counter = 0;
                iguard = oiter->end();
                for (iiter = oiter->begin(); iiter != iguard; iiter++) {
                    if (iiter->getModule().getRing() > counter) {
                        std::vector<ModuleCap>::iterator partner;
                        std::ostringstream matname, shapename;
                        matname << xml_base_actcomp << "L" << layer << "P" << iiter->getModule().getRing();
                        c.push_back(createComposite(matname.str(), compositeDensity(*iiter), *iiter));
                        shapename << iiter->getModule().getRing() << lname.str();
                        shape.type = bx;
                        shape.name_tag = xml_barrel_module + shapename.str();
                        shape.dx = iiter->getModule().getThickness() / 2.0;
                        shape.dy = iiter->getModule().getArea() / iiter->getModule().getHeight() / 2.0;
                        shape.dz = iiter->getModule().getHeight() / 2.0;
                        s.push_back(shape);
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                        logic.material_tag = xml_material_air;
                        l.push_back(logic);
                        partner = findPartnerModule(iiter, iguard, iiter->getModule().getRing());
                        pos.parent_tag = xml_fileident + ":" + rname.str();
                        pos.child_tag = logic.shape_tag;
                        if (((iiter->getModule().getMinRho() > rmin)
                                && (iiter->getModule().getMinRho() < rmin + (rmin + rmax) / 2.0))
                                || (iiter->getModule().getMaxRho() == rmax)) pos.trans.dx = deltar / 2.0 - shape.dx;
                        else pos.trans.dx = shape.dx - deltar / 2.0;
                        if (iiter->getModule().getMaxZ() > 0) {
                            pos.trans.dz = shape.dz + iiter->getModule().getMinZ();
                            p.push_back(pos);
                            if (partner != iguard) {
                                if (((partner->getModule().getMinRho() > rmin)
                                        && (partner->getModule().getMinRho() < rmin + (rmin + rmax) / 2.0))
                                        || (partner->getModule().getMaxRho() == rmax)) pos.trans.dx = deltar / 2.0 - shape.dx;
                                else pos.trans.dx = shape.dx - deltar / 2.0;
                                pos.trans.dz = shape.dz - iiter->getModule().getMinZ();
                                pos.copy = 2;
                                p.push_back(pos);
                                pos.copy = 1;
                            }
                        }
                        else {
                            pos.trans.dz = shape.dz - iiter->getModule().getMinZ();
                            pos.copy = 2;
                            p.push_back(pos);
                            pos.copy = 1;
                            if (partner != iguard) {
                                if (((partner->getModule().getMinRho() > rmin)
                                        && (partner->getModule().getMinRho() < rmin + (rmin + rmax) / 2.0))
                                        || (partner->getModule().getMaxRho() == rmax)) pos.trans.dx = deltar / 2.0 - shape.dx;
                                else pos.trans.dx = shape.dx - deltar / 2.0;
                                pos.trans.dz = shape.dz + iiter->getModule().getMinZ();
                                p.push_back(pos);
                            }
                        }
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
                        shape.name_tag = xml_barrel_module + shapename.str() + xml_base_act;
                        s.push_back(shape);
                        pos.parent_tag = logic.shape_tag;
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                        logic.material_tag = xml_fileident + ":" + matname.str();
                        l.push_back(logic);
                        pos.child_tag = logic.shape_tag;
                        p.push_back(pos);
                        counter++;
                    }
                }
                shape.name_tag = rname.str();
                shape.dx = deltar / 2.0;
                shape.dz = zmax;
                s.push_back(shape);
                logic.name_tag = shape.name_tag;
                logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                logic.material_tag = xml_material_air;
                l.push_back(logic);
                alg.parameters.push_back(stringParam(xml_childparam, logic.shape_tag));
                shape.type = tb;
                shape.dx = 0.0;
                shape.dxx = 0.0;
                shape.dy = 0.0;
                pos.trans.dx = 0.0;
                pos.trans.dz = 0.0;
                shape.name_tag = lname.str();
                shape.rmin = rmin;
                shape.rmax = rmax;
                shape.dz = zmax;
                s.push_back(shape);
                logic.name_tag = shape.name_tag;
                logic.shape_tag = xml_fileident + ":" + logic.name_tag;
                logic.material_tag = xml_material_air;
                l.push_back(logic);
                pos.child_tag = logic.shape_tag;
                p.push_back(pos);
                alg.parent = logic.shape_tag;
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
                pconverter << (umin + deltar / 2.0) << "*mm";
                alg.parameters.push_back(numericParam(xml_radiusout, pconverter.str()));
                pconverter.str("");
                alg.parameters.push_back(numericParam(xml_zposition, "0.0*mm"));
                pconverter << static_cast<BarrelLayer*>(tr.getBarrelLayers()->at(layer - 1))->getRods();
                alg.parameters.push_back(numericParam(xml_number, pconverter.str()));
                alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
                alg.parameters.push_back(numericParam(xml_incrcopyno, "1"));
                a.push_back(alg);
                alg.parameters.clear();
            }
            layer++;
        }
        // e_mod: one composite for every ring
        // s and l: one entry for every ring module, one for every ring in z+, one for every disc in z+
        // p: TBD, TODO
        /*layer = 1;
         * alg.name = xml_tobalgo; //TODO: placeholder, replace with endcap algo name
         * oguard = ec.end();
         * for (oiter = ec.begin(); oiter != oguard; oiter++) {
         * double rmin = tr.getEndcapLayers()->at(layer - 1)->getMinRho();
         * double rmax = tr.getEndcapLayers()->at(layer - 1)->getMaxRho();
         * double zmax = tr.getEndcapLayers()->at(layer - 1)->getMaxZ();
         * std::ostringstream dname;
         * dname << xml_disc << layer;
         * counter = 0;
         * iguard = oiter->end();
         * for (iiter = oiter->begin(); iiter != iguard; iiter++) {
         * if (iiter->getModule().getRing() > counter) {
         * std::ostringstream matname;
         * matname << xml_base_actcomp << "D" << layer << "R" << iiter->getModule().getRing();
         * c.push_back(createComposite(matname.str(), compositeDensity(*iiter), *iiter));
         * counter++;
         * }
         * }
         * //disc here
         * layer++;
         * }*/
        // all tubes from now on, only translations in z
        shape.type = tb;
        shape.dx = 0.0;
        shape.dxx = 0.0;
        shape.dy = 0.0;
        pos.trans.dx = 0.0;
        pos.trans.dz = 0.0;
        // b_ser, e_ser: one composite for every service volume on the z+ side
        // s, l and p: one entry per service volume
        std::vector<InactiveElement>& bs = is.getBarrelServices();
        std::vector<InactiveElement>& es = is.getEndcapServices();
        guard = bs.end();
        for (iter = bs.begin(); iter != guard; iter++) {
            std::ostringstream matname, shapename;
            matname << xml_base_serfcomp << iter->getCategory() << "R" << (int)(iter->getInnerRadius());
            shapename << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset()));
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
        /*guard = es.end();
         * for (iter = es.begin(); iter != guard; iter++) {
         * std::ostringstream matname, shapename;
         * matname << xml_base_serfcomp << iter->getCategory() << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0));
         * shapename << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0));
         * if ((iter->getZOffset() + iter->getZLength()) > 0) c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
         * shape.name_tag = shapename.str();
         * shape.dz = iter->getZLength() / 2.0;
         * shape.rmin = iter->getInnerRadius();
         * shape.rmax = shape.rmin + iter->getRWidth();
         * s.push_back(shape);
         * logic.name_tag = shapename.str();
         * logic.shape_tag = xml_fileident + ":" + shapename.str();
         * logic.material_tag = xml_fileident + ":" + matname.str();
         * l.push_back(logic);
         * pos.parent_tag = xml_fileident + ":" + xml_tracker;
         * pos.child_tag = logic.shape_tag;
         * pos.trans.dz = iter->getZOffset() + shape.dz;
         * p.push_back(pos);
         * }*/
        // b_sup, e_sup, o_sup, t_sup, u_sup: one composite per category
        // l, s and p: one entry per support part
        std::set<MaterialProperties::Category> found;
        std::set<MaterialProperties::Category>::iterator fres;
        std::vector<InactiveElement>& sp = is.getSupports();
        guard = sp.end();
        for (iter = sp.begin(); iter != guard; iter++) {
            if (iter->getCategory() != MaterialProperties::e_sup) { //TODO: remove later
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
                if (iter->getCategory() == MaterialProperties::o_sup) {
                    std::string tmp_l = xml_fileident + ":" + xml_layer;
                    std::string tmp_d = xml_fileident + ":" + xml_disc;
                    std::vector<PosInfo>::iterator iter, guard = p.end();
                    for (iter = p.begin(); iter != guard; iter++) {
                        if (((iter->child_tag).substr(0, tmp_l.size()).compare(tmp_l) == 0)
                                || ((iter->child_tag).substr(0, tmp_d.size()).compare(tmp_d) == 0)) iter->parent_tag = logic.shape_tag;
                    }
                }
                p.push_back(pos);
            }
        }
        //TODO: traverse collections in mb and fill up l, s and p
        //            find out how to handle <algorithm> blocks
        //FOR NOW: barrel only
        //            find out about DOWN configuration (inner barrel longer than outer one)
        //            => not necessarily a problem if going from tracker straight to layers
        //            find out about short layers and stacked layers
        //            => should be straightforward: treat stacked layers as two individual ones,
        //                  short layers have two copies instead of one and translations in +z and -z
        //loop through layers
        //        see notes for rest of hierarchy
        //loop through discs
        //        see XML files for hierarchy
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
        for (unsigned int i = 0; i < comp.elements.size(); i++) comp.elements.at(i).second = comp.elements.at(i).second / mp.getTotalMass();
        return comp;
    }
    
    std::vector<ModuleCap>::iterator tk2CMSSW::findPartnerModule(std::vector<ModuleCap>::iterator i,
            std::vector<ModuleCap>::iterator g, int ponrod) {
        std::vector<ModuleCap>::iterator res = i;
        if (i != g) {
            bool plus = i->getModule().getMaxZ() > 0;
            while (res != g) {
                if (res->getModule().getRing() == ponrod) {
                    if((plus && (res->getModule().getMaxZ() < 0)) || (!plus && (res->getModule().getMaxZ() > 0))) break;
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
    
    double tk2CMSSW::findMinRhoUpperRod(std::vector<Module*>::iterator start,
            std::vector<Module*>::iterator stop, double middle, double max) {
        std::vector<Module*>::iterator iter;
        double u = 0.0;
        for (iter = start; iter != stop; iter++) {
            if (((*iter)->getMeanPoint().Rho() > middle) && ((*iter)->getMaxRho() != max)) {
                u = (*iter)->getMinRho();
                break;
            }
        }
        return u;
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
