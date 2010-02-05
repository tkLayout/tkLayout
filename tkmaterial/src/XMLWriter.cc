/**
 * @file XMLWriter.cc
 * @brief
 */

#include <XMLWriter.h>
namespace insur {
    //public
    void XMLWriter::writeXML(std::vector<Element> e, std::vector<Composite> c, std::vector<LogicalInfo> l, std::vector<ShapeInfo> s,
                                 std::vector<PosInfo> p, std::vector<AlgoInfo> a, std::vector<SpecParInfo> t, std::ostringstream& gb,
                                 std::ostringstream& tb, std::ostringstream& pb, std::ostringstream& sb, std::ostringstream& mb) {
        materialSection(xml_trackerfile, e, c, gb);
        logicalPartSection(l, xml_trackerfile, gb);
        solidSection(s, xml_trackerfile, gb);
        posPartSection(p, a, xml_trackerfile, gb);
        specParSection(t, xml_specpars_label, tb);
        prodcuts(t, pb);
        trackersens(t, sb);
        recomaterial(t, mb);
    }
    //protected
    void XMLWriter::prodcuts(std::vector<SpecParInfo>& t, std::ostringstream& stream) {
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
    
    void XMLWriter::trackersens(std::vector<SpecParInfo>& t, std::ostringstream& stream) {
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
    
    void XMLWriter::recomaterial(std::vector<SpecParInfo>& t, std::ostringstream& stream) {
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
    
    void XMLWriter::materialSection(std::string name , std::vector<Element>& e, std::vector<Composite>& c, std::ostringstream& stream) {
        stream << xml_material_section_open << name << xml_general_inter;
        for (unsigned int i = 0; i < e.size(); i++) elementaryMaterial(e.at(i).tag, e.at(i).density, e.at(i).atomic_number, e.at(i).atomic_weight, stream);
        for (unsigned int i = 0; i < c.size(); i++) compositeMaterial(c.at(i).name, c.at(i).density, c.at(i).method, c.at(i).elements, stream);
        stream << xml_material_section_close;
    }
    
    void XMLWriter::logicalPartSection(std::vector<LogicalInfo>& l, std::string label, std::ostringstream& stream) {
        std::vector<LogicalInfo>::const_iterator iter, guard = l.end();
        stream << xml_logical_part_section_open << label << xml_general_inter;
        for (iter = l.begin(); iter != guard; iter++) logicalPart(iter->name_tag, iter->shape_tag, iter->material_tag, stream);
        stream << xml_logical_part_section_close;
    }
    
    void XMLWriter::solidSection(std::vector<ShapeInfo>& s, std::string label, std::ostringstream& stream) {
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
    
    void XMLWriter::posPartSection(std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::string label, std::ostringstream& stream) {
        std::vector<PosInfo>::iterator piter, pguard = p.end();
        std::vector<AlgoInfo>::iterator aiter, aguard = a.end();
        stream << xml_pos_part_section_open << label << xml_general_inter;
        for (piter = p.begin(); piter != pguard; piter++) posPart(piter->parent_tag, piter->child_tag, piter->rot, piter->trans, piter->copy, stream);
        for (aiter = a.begin(); aiter != aguard; aiter++) algorithm(aiter->name, aiter->parent, aiter->parameters, stream);
        stream << xml_pos_part_section_close;
    }
    
    void XMLWriter::specParSection(std::vector<SpecParInfo>& t, std::string label, std::ostringstream& stream) {
        std::vector<SpecParInfo>::iterator titer, tguard = t.end();
        stream << xml_spec_par_section_open << label << xml_general_inter;
        for (titer = t.begin(); titer != tguard; titer++) specPar(titer->name, titer->parameter, titer->partselectors, stream);
        stream << xml_spec_par_section_close;
    }
    
    void XMLWriter::algorithm(std::string name, std::string parent,
            std::vector<std::string>& params, std::ostringstream& stream) {
        stream << xml_algorithm_open << name << xml_algorithm_parent << parent << xml_general_endline;
        for (unsigned int i = 0; i < params.size(); i++) stream << params.at(i);
        stream << xml_algorithm_close;
    }
    
    void XMLWriter::elementaryMaterial(std::string tag, double density, int a_number,
            double a_weight, std::ostringstream& stream) {
        stream << xml_elementary_material_open << tag << xml_elementary_material_first_inter << tag;
        stream << xml_elementary_material_second_inter << a_number << xml_elementary_material_third_inter;
        stream << a_weight << xml_elementary_material_fourth_inter << density;
        stream << xml_elementary_material_close;
    }
    
    void XMLWriter::compositeMaterial(std::string name, double density, CompType method,
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
    
    void XMLWriter::logicalPart(std::string name, std::string solid, std::string material, std::ostringstream& stream) {
        stream << xml_logical_part_open << name << xml_logical_part_first_inter << solid;
        stream << xml_logical_part_second_inter << material << xml_logical_part_close;
    }
    
    void XMLWriter::box(std::string name, double dx, double dy, double dz, std::ostringstream& stream) {
        stream << xml_box_open << name << xml_box_first_inter << dx << xml_box_second_inter << dy;
        stream << xml_box_third_inter << dz << xml_box_close;
    }
    
    void XMLWriter::trapezoid(std::string name, double dx, double dxx, double dy, double dz, std::ostringstream& stream) {
        stream << xml_trapezoid_open << name << xml_trapezoid_first_inter << dx;
        stream << xml_trapezoid_second_inter << dxx << xml_trapezoid_third_inter << dy;
        stream << xml_trapezoid_fourth_inter << dy << xml_trapezoid_fifth_inter << dz;
        stream << xml_trapezoid_close;
    }
    
    void XMLWriter::tubs(std::string name, double rmin, double rmax, double dz, std::ostringstream& stream) {
        stream << xml_tubs_open << name << xml_tubs_first_inter << rmin << xml_tubs_second_inter << rmax;
        stream << xml_tubs_third_inter << dz << xml_tubs_close;
    }
    
    void XMLWriter::polycone(std::string name, std::vector<std::pair<double, double> >& rzu,
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
    
    void XMLWriter::posPart(std::string parent, std::string child, Rotation& rot, Translation& trans, int copy, std::ostringstream& stream) {
        stream << xml_pos_part_open << copy << xml_pos_part_first_inter << parent;
        stream << xml_pos_part_second_inter << child << xml_general_endline;
        if (!rot.name.empty()) rotation(rot.name, rot.phix, rot.phiy, rot.phiz, rot.thetax, rot.thetay, rot.thetaz, stream);
        if (!(trans.dx == 0.0 && trans.dy == 0.0 && trans.dz == 0.0)) translation(trans.dx, trans.dy, trans.dz, stream);
        stream << xml_pos_part_close;
    }
    
    void XMLWriter::rotation(std::string name, double phix, double phiy, double phiz,
            double thetax, double thetay, double thetaz, std::ostringstream& stream) {
        stream << xml_rotation_open << name << xml_rotation_first_inter << phix << xml_rotation_second_inter << phiy;
        stream << xml_rotation_third_inter << phiz << xml_rotation_fourth_inter << thetax << xml_rotation_fifth_inter;
        stream << thetay << xml_rotation_sixth_inter << thetaz << xml_rotation_close;
    }
    
    void XMLWriter::translation(double x, double y, double z, std::ostringstream& stream) {
        stream << xml_translation_open << x << xml_translation_first_inter << y << xml_translation_second_inter << z;
        stream << xml_translation_close;
    }
        
    void XMLWriter::specPar(std::string name, std::pair<std::string, std::string> param,
            std::vector<std::string>& partsel, std::ostringstream& stream) {
        stream << xml_spec_par_open << name << xml_general_inter;
        for (unsigned i = 0; i < partsel.size(); i++) {
            stream << xml_spec_par_selector << partsel.at(i) << xml_general_endline;
        }
        stream << xml_spec_par_parameter_first << param.first << xml_spec_par_parameter_second;
        stream << param.second << xml_spec_par_close;
    }
    //private
    std::vector<std::pair<std::string, std::vector<std::string> > >& XMLWriter::buildPaths(std::vector<SpecParInfo>& specs,
            std::vector<std::pair<std::string, std::vector<std::string> > >& blocks) {
        std::vector<std::pair<std::string, std::vector<std::string> > >::iterator existing;
        std::string prefix, postfix, spname;
        std::vector<std::string> paths;
        int dindex, rindex, mindex;
        blocks.clear();
        //TOB
        rindex = findEntry(specs, xml_subdet_rod + xml_par_tail);
        mindex = findEntry(specs, xml_subdet_tobdet + xml_par_tail);
        if ((rindex >= 0) && (mindex >= 0)) {
            // rod loop
            for (unsigned int i = 0; i < specs.at(rindex).partselectors.size(); i++) {
                std::string rnumber, mnumber, plusminus;
                std::string& rcurrent = specs.at(rindex).partselectors.at(i);
                if ((rcurrent.size() > xml_plus.size())
                        && (rcurrent.substr(rcurrent.size() - xml_plus.size()).compare(xml_plus) == 0))
                    plusminus = rcurrent.substr(rcurrent.size() - xml_plus.size());
                if ((rcurrent.size() > xml_minus.size())
                        && (rcurrent.substr(rcurrent.size() - xml_minus.size()).compare(xml_minus) == 0))
                    plusminus = rcurrent.substr(rcurrent.size() - xml_minus.size());
                rnumber = rcurrent.substr(xml_rod.size());
                rnumber = rnumber.substr(0, rnumber.size() - plusminus.size());
                spname = xml_tob_prefix + xml_layer + rnumber;
                prefix = xml_tob + "/" + xml_layer + rnumber + "/" + rcurrent;
                // module loop
                for (unsigned int j = 0; j < specs.at(mindex).partselectors.size(); j++) {
                    mnumber = specs.at(mindex).partselectors.at(j).substr(xml_barrel_module.size());
                    mnumber = mnumber.substr(0, mnumber.size() - xml_base_act.size());
                    mnumber = mnumber.substr(findNumericPrefixSize(mnumber) + xml_layer.size());
                    // matching layers
                    if (mnumber.compare(rnumber) == 0) {
                        postfix = specs.at(mindex).partselectors.at(j);
                        postfix = postfix.substr(0, postfix.size() - xml_base_act.size());
                        postfix = postfix + "/" + postfix + xml_base_waf + "/" + specs.at(mindex).partselectors.at(j);
                        paths.push_back(prefix + "/" + postfix);
                    }
                }
                existing = findEntry(spname, blocks);
                if (existing != blocks.end()) existing->second.insert(existing->second.end(), paths.begin(), paths.end());
                else blocks.push_back(std::pair<std::string, std::vector<std::string> >(spname, paths));
                paths.clear();
            }
        }
        //TID
        dindex = findEntry(specs, xml_subdet_wheel + xml_par_tail);
        rindex = findEntry(specs, xml_subdet_ring + xml_par_tail);
        if ((dindex >= 0) && (rindex >= 0)) {
            // disc loop
            for (unsigned int i = 0; i < specs.at(dindex).partselectors.size(); i++) {
                std::string dnumber, rnumber;
                bool plus;
                dnumber = specs.at(dindex).partselectors.at(i).substr(xml_disc.size());
                if ((int)specs.at(dindex).partselectors.size() / 2 < atoi(dnumber.c_str())) plus = true;
                else plus = false;
                spname = xml_tid_prefix + dnumber;
                if (plus) spname = spname + xml_forward;
                else spname = spname + xml_backward;
                if (plus) prefix = xml_tidf;
                else prefix = xml_tidb;
                prefix = prefix + "/" + xml_disc + dnumber;
                // ring loop
                for (unsigned int j = 0; j < specs.at(rindex).partselectors.size(); j++) {
                    std::string compstr = specs.at(rindex).partselectors.at(j);
                    compstr = compstr.substr(compstr.size() - specs.at(dindex).partselectors.at(i).size());
                    // matching discs
                    if (specs.at(dindex).partselectors.at(i).compare(compstr) == 0) {
                        rnumber = specs.at(rindex).partselectors.at(j).substr(xml_ring.size());
                        rnumber = rnumber.substr(0, findNumericPrefixSize(rnumber));
                        postfix = xml_endcap_module + rnumber + xml_disc + dnumber;
                        postfix = postfix + "/" + postfix + xml_base_waf + "/" + postfix + xml_base_act;
                        postfix = specs.at(rindex).partselectors.at(j) + "/" + postfix;
                        paths.push_back(prefix + "/" + postfix);
                    }
                }
                existing = findEntry(spname, blocks);
                if (existing != blocks.end()) existing->second.insert(existing->second.end(), paths.begin(), paths.end());
                else blocks.push_back(std::pair<std::string, std::vector<std::string> >(spname, paths));
                paths.clear();
            }
        }
        return blocks;
    }
    
    bool XMLWriter::endcapsInTopology(std::vector<SpecParInfo>& specs) {
        for (unsigned int i = 0; i < specs.size(); i++) {
            if (specs.at(i).name.compare(xml_subdet_tiddet + xml_par_tail) == 0) return true;
        }
        return false;
    }
        
    int XMLWriter::findNumericPrefixSize(std::string s) {
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
    
    int XMLWriter::findEntry(std::vector<SpecParInfo>& specs, std::string name) {
        int index = 0;
        while (index < (int)(specs.size())) {
            if (specs.at(index).name.compare(name) == 0) return index;
            index++;
        }
        return -1;
    }
    
    std::vector<std::pair<std::string, std::vector<std::string> > >::iterator XMLWriter::findEntry(std::string name,
            std::vector<std::pair<std::string, std::vector<std::string> > >& data) {
        std::vector<std::pair<std::string, std::vector<std::string> > >::iterator result = data.begin();
        std::vector<std::pair<std::string, std::vector<std::string> > >::iterator guard = data.end();
        while (result != guard) {
            if ((result->first).compare(name) == 0) return result;
            result++;
        }
        return result;
    }
}
