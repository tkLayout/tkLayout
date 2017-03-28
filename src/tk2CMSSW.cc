/**
 * @file tk2CMSSW.cc
 * @brief This class provides the interface to analyse a tracker and write the results to XML files for CMSSW
 */

#include <SvnRevision.hh>
#include <tk2CMSSW.hh>

namespace insur {

    // public
    /**
     * This is the main translation function if a tkgeometry tracker model needs to be represented in an XML form that CMSSW
     * can understand. It splits the work into three main parts: filesystem gymnastics related to input and output, analysis of the
     * tkgeometry model, and generation of XML output. Analysis and output generation are delegated to internal instances of
     * <i>Extractor</i> and <i>XMLWriter</i>, respectively. This function deals mainly with the file system and makes sure
     * that the generated output files go where they are supposed to go. The user may specify the name of a subdirectory for the
     * output files. If that directory exists, all of its contents will be overwritten by the new files. If no name is given, a default is
     * used instead.
     * @mt A refernce to the global material table
     * @mb A reference to an existing material budget that serves as input to the translation
     * @xmlOutputName A string with the name of a subfolder for the output; empty by default.
     */
  void tk2CMSSW::translate(MaterialTable& mt, MaterialBudget& mb, XmlTags& trackerXmlTags, std::string xmlDirectoryPath, std::string xmlOutputPath, std::string xmlOutputName, bool wt) {

    bool isPixelTracker = mb.getTracker().isPixelTracker();

        // analyse tracker system and build up collection of elements, composites, hierarchy, shapes, positions, algorithms and topology
        // ex is an instance of Extractor class
        ex.analyse(mt, mb, trackerXmlTags, data, wt);

        std::stringstream simpleHeaderStream;
	std::string metdataFileName = xmlOutputName + "_" + currentDateTime(false) + ".cfg";
        writeSimpleHeader(simpleHeaderStream, metdataFileName);

	// wr is an instance of XMLWriter class
        wr.setSimpleHeader(simpleHeaderStream.str());
        
        // translate collected information to XML and write to buffers
        std::ifstream instream;
        std::ofstream outstream;

	    if (!isPixelTracker) {
	      outstream.open((xmlOutputPath + metdataFileName).c_str());
	      if (outstream.fail()) throw std::runtime_error("Error opening XML generation metadata file for writing.");
	      writeMetadata(outstream);
	      if (outstream.fail()) throw std::runtime_error("Error writing to XML generation metadata file.");
	      outstream.close();
	      outstream.clear();
	      std::cout << "CMSSW XML generation metadata has been written to " << xmlOutputPath << metdataFileName << std::endl;

	      if (!wt) {
		instream.open((xmlDirectoryPath + "/" + xml_pixbarfile).c_str());
		outstream.open((xmlOutputPath + xml_pixbarfile).c_str());
		if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the pixbar files.");
		//writeSimpleHeader(outstream);
		wr.pixbar(data.shapes, instream, outstream);
		if (outstream.fail()) throw std::runtime_error("Error writing to pixbar file.");
		instream.close();
		instream.clear();
		outstream.close();
		outstream.clear();
		std::cout << "CMSSW modified pixel barrel has been written to " << xmlOutputPath << xml_pixbarfile << std::endl;

		instream.open((xmlDirectoryPath + "/" + xml_pixfwdfile).c_str());
		outstream.open((xmlOutputPath + xml_pixfwdfile).c_str());
		if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the pixfwdn files.");
		//writeSimpleHeader(outstream);
		wr.pixfwd(data.shapes, instream, outstream);
		if (outstream.fail()) throw std::runtime_error("Error writing to pixfwd file.");
		instream.close();
		instream.clear();
		outstream.close();
		outstream.clear();
		std::cout << "CMSSW modified pixel endcap has been written to " << xmlOutputPath << xml_pixfwdfile << std::endl;
	      }
	    }

	    if (wt) outstream.open((xmlOutputPath + xml_newtrackerfile).c_str());
	    else outstream.open((xmlOutputPath + trackerXmlTags.trackerfile).c_str());
            if (outstream.fail()) throw std::runtime_error("Error opening tracker file for writing.");
            std::ifstream trackerVolumeTemplate((xmlDirectoryPath + "/" + xml_trackervolumefile).c_str());
            wr.tracker(data, outstream, trackerVolumeTemplate, isPixelTracker, trackerXmlTags, wt);
            if (outstream.fail()) throw std::runtime_error("Error writing to tracker file.");
            outstream.close();
            outstream.clear();
            std::cout << "CMSSW tracker geometry output has been written to " << xmlOutputPath << (wt ? xml_newtrackerfile : trackerXmlTags.trackerfile) << std::endl;

	    if (wt) instream.open((xmlDirectoryPath + "/" + xml_newtopologyfile).c_str());
	    else instream.open((xmlDirectoryPath + "/" + xml_topologyfile).c_str());
	    outstream.open((xmlOutputPath + trackerXmlTags.topologyfile).c_str());
            if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the topology files.");
            //writeSimpleHeader(outstream);
            wr.topology(data.specs, instream, outstream, isPixelTracker, trackerXmlTags);
            if (outstream.fail()) throw std::runtime_error("Error writing to topology file.");
            instream.close();
            instream.clear();
            outstream.close();
            outstream.clear();
	    std::cout << "CMSSW topology output has been written to " << xmlOutputPath << trackerXmlTags.topologyfile << std::endl;

	    instream.open((xmlDirectoryPath + "/" + trackerXmlTags.prodcutsfile).c_str());
	    outstream.open((xmlOutputPath + trackerXmlTags.prodcutsfile).c_str());
            if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the prodcuts files.");
            //writeSimpleHeader(outstream);
            wr.prodcuts(data.specs, instream, outstream, isPixelTracker, trackerXmlTags);
            if (outstream.fail()) throw std::runtime_error("Error writing to prodcuts file.");
            instream.close();
            instream.clear();
            outstream.close();
            outstream.clear();
            std::cout << "CMSSW prodcuts output has been written to " << xmlOutputPath << trackerXmlTags.prodcutsfile << std::endl;

	    instream.open((xmlDirectoryPath + "/" + trackerXmlTags.trackersensfile).c_str());
	    outstream.open((xmlOutputPath + trackerXmlTags.trackersensfile).c_str());
            if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the trackersens files.");
            //writeSimpleHeader(outstream);
            wr.trackersens(data.specs, instream, outstream, isPixelTracker, trackerXmlTags);
            if (outstream.fail()) throw std::runtime_error("Error writing trackersens to file.");
            instream.close();
            instream.clear();
            outstream.close();
            outstream.clear();
            std::cout << "CMSSW sensor surface output has been written to " << xmlOutputPath << trackerXmlTags.trackersensfile << std::endl;

	    if (!isPixelTracker) {
	      instream.open((xmlDirectoryPath + "/" + xml_recomatfile).c_str()); // for OT, takes template file.
	      outstream.open((xmlDirectoryPath + "/" + xml_tmppath + xml_recomatfile).c_str());
	    }
	    else {
	      instream.open((xmlDirectoryPath + "/" + xml_tmppath + xml_recomatfile).c_str()); // for PX, takes output file from OT already created.
	      outstream.open((xmlOutputPath + trackerXmlTags.recomatfile).c_str());
	    }
            if (instream.fail() || outstream.fail()) throw std::runtime_error("Error opening one of the recomaterial files.");
            //writeSimpleHeader(outstream);
            wr.recomaterial(data.specs, data.lrilength, instream, outstream, isPixelTracker, trackerXmlTags, wt);
            if (outstream.fail()) throw std::runtime_error("Error writing recomaterial to file.");
            instream.close();
            instream.clear();
            outstream.close();
            outstream.clear();
            std::cout << "CMSSW reco material output has been written to " << xmlOutputPath << trackerXmlTags.recomatfile << std::endl;
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
            std::map<std::string, double>& elements = data.composites.at(i).elements;
            for (const auto& elem : elements) std::cout << "(" << elem.first << ", " << elem.second << ") ";
            std::cout << std::endl;
        }
        std::cout << "rotations: " << data.rots.size() << " entries." << std::endl;
        for (auto const &it : data.rots) {
            std::cout << "name = " << it.second.name << ", thetax = " << it.second.thetax << ", phix = ";
            std::cout << it.second.phix << ", thetay = " << it.second.thetay << ", phiy = " << it.second.phiy;
            std::cout << ", thetaz = " << it.second.thetaz << ", phiz = " << it.second.phiz << std::endl;
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
	        case co: std::cout << "cone, rmin1 = " << data.shapes.at(i).rmin1 << ", rmax1 = " << data.shapes.at(i).rmax1;
	        std::cout << ", rmin2 = " << data.shapes.at(i).rmin2 << ", rmax2 = " << data.shapes.at(i).rmax2;
	        std::cout << ", dz = " << data.shapes.at(i).dz;
	        break;
                case tp: std::cout << "trapezoid, dx = " << data.shapes.at(i).dx << ", dy = " << data.shapes.at(i).dy;
                std::cout << ", dyy = " << data.shapes.at(i).dyy << ", dz = " << data.shapes.at(i).dz;
                break;
                default: std::cout << "unknown shape";
            }
            std::cout << std::endl;
        }
        std::cout << "operations on shapes: " << data.shapeOps.size() << " entries." << std::endl;
        for (unsigned int i = 0; i < data.shapeOps.size(); i++) {
            std::cout << "name_tag = " << data.shapeOps.at(i).name_tag << ", type = ";
            switch (data.shapeOps.at(i).type) {
	        case uni: std::cout << "union, rSolid1 = " << data.shapeOps.at(i).rSolid1 << ", rSolid2 = " << data.shapeOps.at(i).rSolid2;
	        break;
	        case intersec: std::cout << "intersection, rSolid1 = " << data.shapeOps.at(i).rSolid1 << ", rSolid2 = " << data.shapeOps.at(i).rSolid2;
	        break;
	        case substract: std::cout << "substraction, rSolid1 = " << data.shapeOps.at(i).rSolid1 << ", rSolid2 = " << data.shapeOps.at(i).rSolid2 << ", translation = (" << data.shapeOps.at(i).trans.dx << ", " << data.shapeOps.at(i).trans.dy << ", " << data.shapeOps.at(i).trans.dz << ")" << std::endl;
	        break;
	        default: std::cout << "unknown operation";
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



  void tk2CMSSW::writeSimpleHeader(std::ostream& os, std::string& metdataFileName) {
      os << "<!--" << std::endl;
      os << "============= XML GENERATION METADATA HEADER =============" << std::endl;
      os << "tkLayout revision: " << SvnRevision::revisionNumber << std::endl;
      os << "generated by: " << fullUserName() << std::endl;
      os << "generation date: " << currentDateTime(true) << std::endl;
      os << "note: see " << metdataFileName << " for full config files" << std::endl;
      os << "=========== END XML GENERATION METADATA HEADER ===========" << std::endl;
      os << "-->" << std::endl;
    }

  void tk2CMSSW::writeMetadata(std::ofstream& out) {
    out << "============= XML GENERATION METADATA =============" << std::endl;
    out << "tkLayout revision: " << SvnRevision::revisionNumber << std::endl;
    out << "generated by: " << fullUserName() << std::endl;
    out << "generation date: " << currentDateTime(true) << std::endl;
    for (std::vector<ConfigFile>::const_iterator it = getConfigFiles().begin(); it != getConfigFiles().end(); ++it) {
      out << std::endl << std::endl;
      out << "CONFIG FILE: " << it->name << std::endl;
      out << it->content << std::endl; 
    }
    out << "=========== END XML GENERATION METADATA ===========" << std::endl;
  }

    std::string tk2CMSSW::currentDateTime(bool withTime) const {
      time_t     now = time(0);
      struct tm  tstruct;
      char       buf[80];
      tstruct = *localtime(&now);
      if (withTime) strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
      else strftime(buf, sizeof(buf), "%Y_%m_%d", &tstruct);
      return buf;
    }

    std::string tk2CMSSW::fullUserName() const {
      struct passwd* pwd = getpwuid(getuid());
      char hn[256];
      hn[255] = 0;
      gethostname(hn, 255); // if hostname is too long it gets truncated and the string might or might not contain a terminating null byte, therefore we force the last char to be null by construction
      return std::string(pwd->pw_gecos) + " (" + pwd->pw_name + " on " + hn + ")";
    }
    
}
