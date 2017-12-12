/**
 * @file MatParser.cc
 * @brief This parses the global material list and the tracker-specific materials from files, initialising an instance of <i>MatCalc</i> in the process
 */

#include <MatParser.hh>
namespace insur {
    /**
     * Nothing to do for the constructor...
     */
    MatParser::MatParser() {}
    
    /**
     * Nothing to do for the destructor...
     */
    MatParser::~MatParser() {}
    

    bool MatParser::fillTable(std::string materialfile, MaterialTable2& mattab) {
      std::ifstream filein(materialfile.c_str());
      std::string line;
      while (std::getline(filein, line)) {
        mattab.parseMaterial(line);
      }
      return true;
    }




    /**
     * This function initialises the internal material table from a config file.
     * @param materialfile The name and, if necessary, path of the global material config file
     * @param mattab A reference to an empty material table
     * @return True if the config file was successfully parsed, false otherwise
     */
    bool MatParser::fillTable(std::string materialfile, MaterialTable& mattab) {
        // Material file name gymnastics
        if (materialfile.empty()) materialfile = default_mattabfile;
        bfs::path mpath(materialfile);
        if (bfs::exists(mpath)) {
            try {
                std::string line, word;
                std::ifstream infilestream;
                infilestream.open(materialfile.c_str());
                // material file line loop
                while (std::getline(infilestream, line)) {
                    // cosmetics and word extraction preparations
                    balgo::trim(line);
                    std::istringstream wordstream(line);
                    std::vector<std::string> tmp;
                    MaterialRow row;
                    // word loop
                    while (wordstream >> word) {
                        // save everything that is not a comment word for word in a temporary vector
                        if ((word.compare(0, c_comment.size(), c_comment) == 0) || (word.compare(0, shell_comment.size(), shell_comment) == 0)) break;
                        else tmp.push_back(word);
                    }
                    if (tmp.empty()) continue;
                    else {
                        // fill up necessary data with dummy values if there is too little information or complain if there is too much
                        if (tmp.size() < 4) {
                            while (tmp.size() < 4) tmp.push_back(dummy_value);
                        }
                        if (tmp.size() > 4) std::cerr << warning_too_many_values << std::endl;
                        // convert the information in the temporary vector to fill the fields in the material row
                        row.tag = tmp.at(0);
                        row.density = atof(tmp.at(1).c_str());
                        row.rlength = atof(tmp.at(2).c_str());
                        row.ilength = atof(tmp.at(3).c_str());
                        // add the completed struct to the internal material table
                        mattab.addMaterial(row);
                    }
                }
                infilestream.close();
            }
            catch (bfs::filesystem_error& bfe) {
	        std::cerr << bfe.what() << std::endl;
                return false;
            }
            catch (std::bad_alloc& ba) {
                std::cerr << ba.what() << std::endl;
                return false;
            }
            return true;
        }
        else logERROR(msg_no_mat_file);
        return false;
    }
 
    /**
     * This function initialises a material calculator from a config file.
     * @param configfile The name and, if necessary, path of a material config file for a tracker layout
     * @param calc A reference to an uninitialised material calculator
     * @param mattabdir The location of the global material table if it is not in the default place
     * @return True if the config file was successfully parsed, false otherwise
     */
     bool MatParser::initMatCalc(MatCalc& calc, std::string mattabdir ) {
        // fill up the global material table if necessary
       if (calc.getMaterialTable().empty()) {
	 if (mattabdir.empty()) mattabdir = default_mattabdir;
	 std::string filename(mattabdir + "/" + default_mattabfile);
	 if(!fillTable(filename, calc.getMaterialTable())) return false;
	 calc.initDone(true);
	 return true;
       }
       return false;
     }
    
}
