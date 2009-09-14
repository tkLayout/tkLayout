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
        //TODO: analyse tracker system and build up collection of shapes as well as dependency graphs for hierarchy and position
        analyse(mb, elements, composites, logic, shapes, positions);
        std::ostringstream buffer;
        //TODO: actual translation code starts here
        //            write translated xml to buffer first by calling logicalPartSection(), solidSection(),
        //            posPartSection() and algorithm() with the appropriate parameters
        //            find out about DOWN configuration (inner barrel longer than outer one)
        //            find out about short layers and stacked layers
        //            => should be straightforward: treat stacked layers as two individual ones,
        //                  short layers have two copies instead of one and translations in +z and -z
        
        // write results to top-level file
        bfs::remove_all(outpath.c_str());
        bfs::create_directory(outpath);
        std::ofstream outstream((outpath + trackerfile).c_str());
        outstream << buffer.str() << std::endl;
        outstream.close();        
        //dummy code
        std::cout << "tk2CMSSW::translate() has been called for output directory " << outpath << std::endl;
    }
    
    // protected
    void tk2CMSSW::materialSection(std::vector<Element>& e, std::vector<Composite>& c, std::ostringstream& stream) {
        //TODO: compute A, Z for each entry and write info to material block.
        //            PROBLEM: distinction between elementary and composite materials in CMSSW
        //            SOLUTION: treat everything as an elementary material???
    }
    
    void tk2CMSSW::logicalPartSection(std::vector<LogicalInfo>& l, std::string label, std::ostringstream& stream) {
        std::vector<LogicalInfo>::const_iterator iter, guard = l.end();
        stream << xml_logical_part_section_open << label << xml_logical_part_section_inter;
        //TODO: traverse the volume hierarchy graph and write the nodes to stream using a series of logicalPart() calls
        for (iter = l.begin(); iter != guard; iter++) logicalPart();
        stream << xml_logical_part_section_close;
    }
    
    void tk2CMSSW::solidSection(std::vector<ShapeInfo>& s, std::string label, std::ostringstream& stream) {
        std::string name;
        stream << xml_solid_section_open << label << xml_solid_section_inter;
        //TODO: traverse the shapes collection and write the nodes to stream using a series of box(), trapezoid() and tubs() calls
        stream << xml_solid_section_close;
    }
    
    void tk2CMSSW::posPartSection(std::vector<PosInfo>& p, std::string label, std::ostringstream& stream) {
        stream << xml_pos_part_section_open << label << xml_pos_part_section_inter;
        //TODO: position every single volume relative to its parent using a series of posPart() calls
        //            IMPORTANT: find out the role of algorithm() in reducing the number of posPart() calls
        stream << xml_pos_part_section_close;
    }
    
    void tk2CMSSW::algorithm(/*TODO: TBD*/ std::ostringstream& stream) {
        //TODO
    }
    
    void tk2CMSSW::elementaryMaterial(std::string tag, double density, int a_number,
                                                                                         double a_weight, std::ostringstream& stream) {
        stream << xml_elementary_material_open << tag << xml_elementary_material_first_inter << tag;
        stream << xml_elementary_material_second_inter << a_number << xml_elementary_material_third_inter;
        stream << a_weight << xml_elementary_material_fourth_inter << density;
        stream << xml_elementary_material_close;
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
    
    void tk2CMSSW::posPart(int copy, std::string parent, std::string child, Rotation& rot, Translation& trans, std::ostringstream& stream) {
        stream << xml_pos_part_open << copy << xml_pos_part_first_inter << parent;
        stream << xml_pos_part_second_inter << child << xml_pos_part_third_inter;
        rotation(rot.name, rot.phix, rot.phiy, rot.phiz, rot.thetax, rot.thetay, rot.thetaz, stream);
        translation(trans.dx, trans.dy, trans.dz, stream);
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
    void tk2CMSSW::analyse(MaterialBudget& mb, std::vector<Element>& elements, std::vector<Composite>& composites,
                                                                                  std::vector<LogicalInfo>& logic, std::vector<ShapeInfo>& shapes, std::vector<PosInfo>& positions) {
        //TODO: traverse collections in mb and fill up logic, shape and pos
        Tracker& tr = mb.getTracker();
        InactiveSurfaces& is = mb.getInactiveSurfaces();
        std::vector<std::vector<ModuleCap> >& bc = mb.getBarrelModuleCaps();
        std::vector<std::vector<ModuleCap> >& ec = mb.getEndcapModuleCaps();
        
    }
    
    void tk2CMSSW::print() {
        //TODO: print internal status
    }
}
