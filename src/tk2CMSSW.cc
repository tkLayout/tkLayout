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
        std::string tmppath = default_xmlpath + "/" + xml_tmppath + "/";
        // analyse tracker system and build up collection of elements, composites, hierarchy, shapes, positions, algorithms and topology
        ex.analyse(mt, mb, data);
        // translate collected information to XML and write to buffers
        std::ifstream instream;
        std::ofstream outstream;
        // replace by calls to pixbar(), pixfwd() (if endcapsInTopology()), tracker(),
        // topology(), prodcuts(), trackersens(), recomaterial() - probably in that order
        //TODO: improve error handling! Rename original outpath directory, if it exists.
        //            Then create new directory and start writing files. If all goes well, remove
        //            temp directory. If not, remove new directory and restore original.
        // use try-catch block for error handling: throw runtime error if one of the output
        // streams comes back in fail state, handle rollback in catch block
        try {
            if (bfs::exists(outpath)) bfs::rename(outpath, tmppath);
            bfs::create_directory(outpath);
            //TODO: use instream and outstream instead, in a series of inits and re-inits until all files have been written.
            instream.open((default_xmlpath + "/" + xml_pixbarfile).c_str());
            outstream.open((outpath + xml_pixbarfile).c_str());
            if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the pixbar files.");
            wr.pixbar(data.shapes, instream, outstream);
            if (outstream.fail()) throw std::runtime_error("Error writing to pixbar file.");
            instream.close();
            instream.clear();
            outstream.close();
            outstream.clear();
            std::cout << "CMSSW modified pixel barrel has been written to " << outpath << xml_pixbarfile << std::endl;
            instream.open((default_xmlpath + "/" + xml_pixfwdfile).c_str());
            outstream.open((outpath + xml_pixfwdfile).c_str());
            if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the pixfwdn files.");
            wr.pixfwd(data.shapes, instream, outstream);
            if (outstream.fail()) throw std::runtime_error("Error writing to pixfwd file.");
            instream.close();
            instream.clear();
            outstream.close();
            outstream.clear();
            std::cout << "CMSSW modified pixel endcap has been written to " << outpath << xml_pixfwdfile << std::endl;
            outstream.open((outpath + xml_trackerfile).c_str());
            if (outstream.fail()) throw std::runtime_error("Error opening tracker file for writing.");
            wr.tracker(data, outstream);
            if (outstream.fail()) throw std::runtime_error("Error writing to tracker file.");
            outstream.close();
            outstream.clear();
            std::cout << "CMSSW tracker geometry output has been written to " << outpath << xml_trackerfile << std::endl;
            instream.open((default_xmlpath + "/" + xml_topologyfile).c_str());
            outstream.open((outpath + xml_topologyfile).c_str());
            if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the topology files.");
            wr.topology(data.specs, instream, outstream);
            if (outstream.fail()) throw std::runtime_error("Error writing to topology file.");
            instream.close();
            instream.clear();
            outstream.close();
            outstream.clear();
            std::cout << "CMSSW topology output has been written to " << outpath << xml_topologyfile << std::endl;
            instream.open((default_xmlpath + "/" + xml_prodcutsfile).c_str());
            outstream.open((outpath + xml_prodcutsfile).c_str());
            if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the prodcuts files.");
            wr.prodcuts(data.specs, instream, outstream);
            if (outstream.fail()) throw std::runtime_error("Error writing to prodcuts file.");
            instream.close();
            instream.clear();
            outstream.close();
            outstream.clear();
            std::cout << "CMSSW prodcuts output has been written to " << outpath << xml_prodcutsfile << std::endl;
            instream.open((default_xmlpath + "/" + xml_trackersensfile).c_str());
            outstream.open((outpath + xml_trackersensfile).c_str());
            if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the trackersens files.");
            wr.trackersens(data.specs, instream, outstream);
            if (outstream.fail()) throw std::runtime_error("Error writing trackersens to file.");
            instream.close();
            instream.clear();
            outstream.close();
            outstream.clear();
            std::cout << "CMSSW sensor surface output has been written to " << outpath << xml_trackersensfile << std::endl;
            instream.open((default_xmlpath + "/" + xml_recomatfile).c_str());
            outstream.open((outpath + xml_recomatfile).c_str());
            if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the recomaterial files.");
            wr.recomaterial(data.specs, data.lrilength, instream, outstream);
            if (outstream.fail()) throw std::runtime_error("Error writing recomaterial to file.");
            instream.close();
            instream.clear();
            outstream.close();
            outstream.clear();
            std::cout << "CMSSW reco material output has been written to " << outpath << xml_recomatfile << std::endl;
            bfs::remove_all(tmppath);
        }
        catch (std::runtime_error& e) {
            std::cerr << "Error writing files: " << e.what() << std::endl;
            if (bfs::exists(outpath)) bfs::remove_all(outpath);
            if (bfs::exists(tmppath)) bfs::rename(tmppath, outpath);
            std::cerr << "No files were changed." <<std::endl;
        }
    }
    
    // private
    /**
     * This prints the contents of the internal CMSSWBundle collection; used for debugging.
     */
    void tk2CMSSW::print() {
        std::cout << "tm2CMSSW internal status:" << std::endl;
        std::cout << "elements: " << data.elements.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < data.elements.size(); i++) {
            std::cout << "entry " << i << ": tag = " << data.elements.at(i).tag << ", density = " << data.elements.at(i).density << ", atomic number = ";
            std::cout << data.elements.at(i).atomic_number << ", atomic weight = " << data.elements.at(i).atomic_weight << std::endl;
        }
        std::cout << "composites: " << data.composites.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < data.composites.size(); i++) {
            std::cout << "entry " << i << ": name = " << data.composites.at(i).name << ", density = " << data.composites.at(i).density << ", method = ";
            switch (data.composites.at(i).method) {
                case wt: std::cout << "fraction by weight";
                break;
                case vl: std::cout << "fraction by volume";
                break;
                case ap: std::cout << "fraction by atomic proportion";
                break;
                default: std::cout << "unknown method";
            }
            std::cout << std::endl << "elements: ";
            std::vector<std::pair<std::string, double> >& elems = data.composites.at(i).elements;
            for (unsigned int j = 0; j < elems.size(); j++) std::cout << "(" << elems.at(j).first << ", " << elems.at(j).second << ") ";
            std::cout << std::endl;
        }
        std::cout << "rotations: " << data.rots.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < data.rots.size(); i++) {
            std::cout << "name = " << data.rots.at(i).name << ", thetax = " << data.rots.at(i).thetax << ", phix = ";
            std::cout << data.rots.at(i).phix << ", thetay = " << data.rots.at(i).thetay << ", phiy = " << data.rots.at(i).phiy;
            std::cout << ", thetaz = " << data.rots.at(i).thetaz << ", phiz = " << data.rots.at(i).phiz << std::endl;
        }
        std::cout << "logic: " << data.logic.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < data.logic.size(); i++) {
            std::cout << "name_tag = " << data.logic.at(i).name_tag << ", shape_tag = " << data.logic.at(i).shape_tag;
            std::cout << ", material_tag = " << data.logic.at(i).material_tag << std::endl;
        }
        std::cout << "shapes: " << data.shapes.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < data.shapes.size(); i++) {
            std::cout << "name_tag = " << data.shapes.at(i).name_tag << ", type = ";
            switch (data.shapes.at(i).type) {
                case bx: std::cout << "box, dx = " << data.shapes.at(i).dx << ", dy = " << data.shapes.at(i).dy << ", dz = ";
                std::cout << data.shapes.at(i).dz;
                break;
                case tb: std::cout << "tube, rmin = " << data.shapes.at(i).rmin << ", rmax = " << data.shapes.at(i).rmax;
                std::cout << ", dz = " << data.shapes.at(i).dz;
                break;
                case tp: std::cout << "trapezoid, dx = " << data.shapes.at(i).dx << ", dy = " << data.shapes.at(i).dy;
                std::cout << ", dyy = " << data.shapes.at(i).dyy << ", dz = " << data.shapes.at(i).dz;
                break;
                default: std::cout << "unknown shape";
            }
            std::cout << std::endl;
        }
        std::cout << "positions: " << data.positions.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < data.positions.size(); i++) {
            std::cout << "parent_tag = " << data.positions.at(i).parent_tag << ", child_tag = " << data.positions.at(i).child_tag;
            std::cout << ", rotref = " << (data.positions.at(i).rotref.empty() ? "[no name]": data.positions.at(i).rotref) << ", ";
            std::cout << ", translation = (" << data.positions.at(i).trans.dx << ", " << data.positions.at(i).trans.dy << ", ";
            std::cout << data.positions.at(i).trans.dz << ")" << std::endl;
        }
        std::cout << "algorithms: " << data.algos.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < data.algos.size(); i++ ) {
            std::cout << "name = " << data.algos.at(i).name << ", parent = " << data.algos.at(i).parent << std::endl;
            std::cout << "parameters:" << std::endl;
            for (unsigned int j = 0; j < data.algos.at(i).parameters.size(); j++) std::cout << data.algos.at(i).parameters.at(j) << std::endl;
        }
        SpecParInfo spec;
        std::cout << "topology: " << data.specs.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < data.specs.size(); i++) {
            std::cout << "name = " << data.specs.at(i).name << std::endl << "partselectors:" << std::endl;
            for (unsigned int j = 0; j < data.specs.at(i).partselectors.size(); j++) std::cout << data.specs.at(i).partselectors.at(j) << std::endl;
            std::cout << "parameter = (" << data.specs.at(i).parameter.first << ", " << data.specs.at(i).parameter.second << ")" << std::endl;
        }
    }
}
