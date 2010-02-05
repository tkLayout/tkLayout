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
        ex.analyse(mt, mb, elements, composites, logic, shapes, positions, algos, specs);
        // translate collected information to XML and write to buffers
        std::ostringstream gbuffer, tbuffer, pbuffer, sbuffer, mbuffer;
        gbuffer << xml_preamble << xml_const_section;
        tbuffer << xml_preamble;
        wr.writeXML(elements, composites, logic, shapes, positions, algos, specs, gbuffer, tbuffer, pbuffer, sbuffer, mbuffer);
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
    
    // private
    void tk2CMSSW::print() { //TODO: complete
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
