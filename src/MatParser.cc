#include <MatParser.h>
namespace insur {
    /**
     * Nothing to do for the constructor...
     */
    MatParser::MatParser() {}
    
    /**
     * Nothing to do for the destructor...
     */
    MatParser::~MatParser() {}
    
    /**
     * This function initialises the internal material table from a config file.
     * @param materialfile The name and, if necessary, path of the global material config file
     * @param mattab A reference to an empty material table
     * @return True if the config file was successfully parsed, false otherwise
     */
    bool MatParser::fillTable(std::string materialfile, MaterialTable& mattab) {
        // Material file name gymnastics
        std::cerr << "Trying to read material file "<<materialfile << std::endl;
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
        else std::cerr << "MatParser::fillTable()" << msg_no_mat_file << std::endl;
        return false;
    }
    
    /**
     * This is the core function that does the actual parsing during initialisation of a material calculator.
     * @param configfile The name and, if necessary, path of a material config file for a tracker layout
     * @param calc A reference to an uninitialised material calculator
     * @return True if the config file was successfully parsed, false otherwise
     */
    bool MatParser::readParameters(std::string configfile, MatCalc& calc) {
        // material file gymnastics
        bfs::path mpath(configfile);
        if (bfs::exists(mpath)) {
            try {
                std::string line, word, type;
                std::ifstream infilestream;
                infilestream.open(configfile.c_str());
                // material file line loop
                while (std::getline(infilestream, line)) {
                    // cosmetics and word extraction preparations
                    balgo::trim(line);
                    std::istringstream wordstream(line);
                    wordstream >> word;
                    // skip comments
                    if ((word.compare(0, c_comment.size(), c_comment) == 0)
                            || (word.compare(0, shell_comment.size(), shell_comment) == 0)) continue;
                    // line defining a module type
                    else if (word.compare(type_marker) == 0) {
                        // extract module type
                        word = getValue(line, false);
                        // set module type is stereo (special case)
                        if (word.compare(type_stereo) == 0) {
                            if (calc.typeRegistered(type_rphi)) {
                                type = word;
                                // record type, strips and segments in an internal data structure
                                calc.addTypeInfo(type, calc.getStripsAcross(type_rphi), calc.getSegmentsAlong(type_rphi));
                                // copy the existing material information from single sided to double sided as it forms the basis
                                calc.copyContents(type_rphi, type);
                            }
                            else std::cerr << "MatParser::readParameters(): " << msg_unexpected_type << std::endl;
                        }
                        // independent module type
                        else {
                            if (!calc.typeRegistered(word)) {
                                type = word;
                                // extract the number of strips and segments from the type definition
                                std::string strips, segs;
                                if (parseStripsSegs(infilestream, strips, segs)) {
                                    int chips, segments;
                                    chips = atoi(strips.c_str());
                                    segments = atoi(segs.c_str());
                                    // record type, strips and segments in an internal data structure
                                    calc.addTypeInfo(type, chips, segments);
                                }
                                else std::cerr << "MatParser::readParameters(): " << msg_readparam_failure << std::endl;
                            }
                            else std::cerr << "MatParser::readParameters(): " << msg_type_exists << std::endl;
                        }
                    }
                    // line defining a module material
                    else if (word.compare(m_line_delim) == 0) {
                        if (!parseMLine(line, type, calc)) std::cout << msg_m_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    // line defining a module-to-service mapping
                    else if(word.compare(d_line_delim) == 0) {
                        if (!parseDLine(line, calc)) std::cout << msg_d_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    // line defining a service material
                    else if (word.compare(s_line_delim) == 0) {
                        if (!parseSimpleLine(line, calc, s_line_delim)) std::cout << msg_s_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    // line defining a material for a barrel support disc
                    else if (word.compare(x_line_delim) == 0) {
                        if (!parseSimpleLine(line, calc, x_line_delim)) std::cout << msg_x_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    // line defining a material for an endcap support tube
                    else if (word.compare(y_line_delim) == 0) {
                        if (!parseSimpleLine(line, calc, y_line_delim)) std::cout << msg_y_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    // line defining a material for the outer support tube
                    else if (word.compare(z_line_delim) == 0) {
                        if (!parseSimpleLine(line, calc, z_line_delim)) std::cout << msg_z_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    // line defining a material for a user-defined barrel support
                    else if (word.compare(w_line_delim) == 0) {
                        if (!parseSimpleLine(line, calc, w_line_delim)) std::cout << msg_w_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    // line defining a material for a barrel support tube
                    else if (word.compare(v_line_delim) == 0) {
                        if (!parseSimpleLine(line, calc, v_line_delim)) std::cout << msg_v_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    // confusing line
                    else if (!line.empty()) {
                        std::cerr << "MatParser::readParameters(): Confusion detected...skipping line '" << line << "'." << std::endl;
                        break;
                    }
                    word.clear();
                }
                infilestream.close();
                return true;
            }
            catch (bfs::filesystem_error& bfe) {
                std::cerr << "MatParser::readParameters(): " << bfe.what() << std::endl;
                calc.clearTypeVector();
                calc.clearModVectors();
                return false;
            }
            catch (std::bad_alloc& ba) {
                std::cerr << "MatParser::readParameters(): " << ba.what() << std::endl;
                calc.clearTypeVector();
                calc.clearModVectors();
                return false;
            }
        }
        else std::cerr << "MatParser::readParameters(): " << msg_no_mat_file << std::endl;
        calc.clearTypeVector();
        calc.clearModVectors();
        return false;
    }
    
    /**
     * This function initialises a material calculator from a config file.
     * @param configfile The name and, if necessary, path of a material config file for a tracker layout
     * @param calc A reference to an uninitialised material calculator
     * @return True if the config file was successfully parsed, false otherwise
     */
  bool MatParser::initMatCalc(std::string configfile, MatCalc& calc, std::string mattabdir ) {
        // fill up the global material table if necessary
        if (calc.getMaterialTable().empty()) {
  	    if (mattabdir.empty()) mattabdir = default_mattabdir;
            std::string filename(mattabdir + "/" + default_mattabfile);
            if(!fillTable(filename, calc.getMaterialTable())) return false;
        }
        // read the material config file provided by the user
        calc.initDone(readParameters(configfile, calc));
        if (calc.initDone()) return true;
        return false;
    }
    
    // protected
    /**
     * This function parses the default dimensions of a module type and stores them in the provided string objects.
     * @param instream A reference to the input file stream pointing to the config file
     * @param strips A reference to the output string for the parsed number of strips
     * @param segs A reference to the output string for the parsed number of segments
     * @return True if there were no errors during parsing, false otherwise
     */
    bool MatParser::parseStripsSegs(std::ifstream& instream, std::string& strips, std::string& segs) {
        bool strips_done = false, segs_done = false;
        std::string line;
        // type block line loop
        while (!(strips_done && segs_done) && std::getline(instream, line)) {
            // cosmetics
            balgo::trim(line);
            // skip comments and empty lines
            if (line.empty()) continue;
            else if((line.compare(0, c_comment.size(), c_comment) == 0) || (line.compare(0, shell_comment.size(), shell_comment) == 0)) continue;
            // extract number of strips
            else if(!strips_done && (line.compare(0, strip_marker.size(), strip_marker) == 0)) {
                strips = readFromLine(line, strip_marker);
                strips_done = true;
            }
            // extract number of segments
            else if (!segs_done && (line.compare(0, seg_marker.size(), seg_marker) == 0)) {
                segs = readFromLine(line, seg_marker);
                segs_done = true;
            }
            else return false;
        }
        if (strips.empty() || segs.empty()) return false;
        return true;
    }
    
    /**
     * This function parses a line containing information about module materials and stores it in the
     * internal data structures of a <i>MatCalc</i> object.
     * @param line The input line that needs to be parsed
     * @param type The module type that the given material information refers to
     * @param calc The material calculator where the extracted information is transferred
     * @return True if there were no errors during parsing, false otherwise
     */
    bool MatParser::parseMLine(std::string line, std::string type, MatCalc& calc) {
        if (type.empty()) return false;
        // set starting and end points of information on line
        unsigned int start = line.find(m_line_delim) + m_line_delim.size();
        unsigned int stop = line.find(line_end_delim);
        if ((start >= line.size()) || (stop == std::string::npos)) return false;
        // clip line to information section and remove whitespace at beginning and end
        line = line.substr(start, stop - start);
        balgo::trim(line);
        // preparations for word extraction
        std::stringstream wordstream(line);
        std::string tag, tmp;
        if (!(wordstream >> tag)) return false;
        // parameter and unit retrieval: if no unit, assume grams
        double A, B, C, D;
        MatCalc::Matunit uA, uB, uC, uD;
        bool local;
        if (!(wordstream >> A)) return false;
        if (!(wordstream >> tmp)) return false;
        try {
            uA = getUnit(tmp);
            if (!(wordstream >> B)) return false;
        }
        catch(std::range_error&) {
            uA = MatCalc::gr;
            B = atof(tmp.c_str());
        }
        if (!(wordstream >> tmp)) return false;
        try {
            uB = getUnit(tmp);
            if (!(wordstream >> C)) return false;
        }
        catch(std::range_error&) {
            uB = MatCalc::gr;
            C = atof(tmp.c_str());
        }
        if (!(wordstream >> tmp)) return false;
        try {
            uC = getUnit(tmp);
            if (!(wordstream >> D)) return false;
        }
        catch(std::range_error&) {
            uC = MatCalc::gr;
            D = atof(tmp.c_str());
        }
        if (!(wordstream >> tmp)) return false;
        try {
            uD =  getUnit(tmp);
            if (!(wordstream >> tmp)) return false;
        }
        catch(std::range_error&) {
            uD = MatCalc::gr;
        }
        // local/exiting marker retrieval
        if (tmp.compare(local_marker) == 0) local = true;
        else if (tmp.compare(exit_marker) == 0) local = false;
        else return false;
        // convert module type from input parameter
        /*MatCalc::Modtype tp;
         * if (type.compare(type_rphi) == 0) tp = MatCalc::rphi;
         * else if (type.compare(type_stereo) == 0) tp = MatCalc::stereo;
         * else if (type.compare(type_pt) == 0) tp = MatCalc::pt;
         * else return false;*/
        // record material parameters, units and marker in an internal data structure
        calc.addModuleParameters(tag, type, A, uA, B, uB, C, uC, D, uD, local);
        return true;
    }
    
    /**
     * This function parses a line containing information about a material mapping at the layer/service
     * boundary and stores it in the internal data structures of a <i>MatCalc</i> object.
     * @param line The input line that needs to be parsed
     * @param calc The material calculator where the extracted information is transferred
     * @return True if there were no errors during parsing, false otherwise
     */
    bool MatParser::parseDLine(std::string line, MatCalc& calc) {
        // set starting and end points of information on line
        unsigned int start = line.find(d_line_delim) +d_line_delim.size();
        unsigned int stop = line.find(line_end_delim);
        if ((start >= line.size()) || (stop == std::string::npos)) return false;
        // clip line to information section and remove whitespace at beginning and end
        line = line.substr(start, stop - start);
        balgo::trim(line);
        // preparations for word extraction
        std::stringstream wordstream(line);
        std::string inTag, outTag, tmp;
        double in, out;
        MatCalc::Matunit uIn, uOut;
        bool local;
        // module material and unit retrieval: if no unit, assume grams
        if (!(wordstream >> in)) return false;
        if (!(wordstream >> tmp)) return false;
        else {
            try { uIn = getUnit(tmp); }
            catch(std::range_error& re) {
                std::cerr << re.what() << std::endl;
                return false;
            }
        }
        if (!(wordstream >> inTag)) return false;
        // service material and unit retrieval: if no unit, assume grams
        if (!(wordstream >> out)) return false;
        if (!(wordstream >> tmp)) return false;
        else {
            try { uOut = getUnit(tmp); }
            catch(std::range_error& re) {
                std::cerr << re.what() << std::endl;
                return false;
            }
        }
        if (!(wordstream >> outTag)) return false;
        if (!(wordstream >> tmp)) return false;
        // local/exiting marker retrieval
        else {
            if (tmp.compare(local_marker) == 0) local = true;
            else if (tmp.compare(exit_marker) == 0) local = false;
            else return false;
        }
        // record material parameters, units and marker in an internal data structure
        calc.addServiceParameters(inTag, in, uIn, outTag, out, uOut, local);
        return true;
    }
    
    /**
     * This function parses a line containing information about local materials within the detector and
     * stores it in the internal data structures of a <i>MatCalc</i> object.
     * @param line The input line that needs to be parsed
     * @param calc The material calculator where the extracted information is transferred
     * @param marker A string identifying the material category
     * @return True if there were no errors during parsing, false otherwise
     */
    bool MatParser::parseSimpleLine(std::string line, MatCalc& calc, std::string marker) {
        // set starting and end points of information on line
        unsigned int start = line.find(marker) + marker.size();
        unsigned int stop = line.find(line_end_delim);
        if ((start >= line.size()) || (stop == std::string::npos)) return false;
        // clip line to information section and remove whitespace at beginning and end
        line = line.substr(start, stop - start);
        balgo::trim(line);
        // preparations for word extraction
        std::stringstream wordstream(line);
        std::string tag, tmp;
        // material and unit retrieval: if no unit, assume grams
        if (!(wordstream >> tag)) return false;
        double val;
        MatCalc::Matunit uni;
        if (!(wordstream >> val)) return false;
        if (val != 0.0) {
            if (!(wordstream >> tmp)) return false;
            try { uni = getUnit(tmp); }
            catch (std::range_error& re) {
                std::cerr << re.what() << std::endl;
                return false;
            }
        }
        else uni = MatCalc::gr;
        // record material parameter and unit in local services internal data structure
        if (marker.compare(s_line_delim) == 0) calc.addServiceParameters(tag, val, uni);
        // record material parameter and unit under barrel support discs in internal data structure
        else if (marker.compare(x_line_delim) == 0) calc.addSupportParameters(tag, val, uni, MaterialProperties::b_sup);
        // record material parameter and unit under endcap support tubes in internal data structure
        else if (marker.compare(y_line_delim) == 0) calc.addSupportParameters(tag, val, uni, MaterialProperties::e_sup);
        // record material parameter and unit under outer support tube in internal data structure
        else if (marker.compare(z_line_delim) == 0) calc.addSupportParameters(tag, val, uni, MaterialProperties::o_sup);
        // record material parameter and unit under user-defined barrel supports in internal data structure
        else if (marker.compare(w_line_delim) == 0) calc.addSupportParameters(tag, val, uni, MaterialProperties::u_sup);
        // record material parameter and unit under barrel support tubes in internal data structure
        else if (marker.compare(v_line_delim) == 0) calc.addSupportParameters(tag, val, uni, MaterialProperties::t_sup);
        else return false;
        return true;
    }
    
    // private
    /**
     * This convenience function reads a parameter value for the given parameter name from a line in the config file.
     * Depending on whether the parameter name occurs in it or not, the function may return an empty string.
     * @param source The line as taken from the config file
     * @param paramname The parameter identifier that the requested value belongs to
     * @return The requested value in string form
     */
    std::string MatParser::readFromLine(std::string source, std::string paramname) {
        std::string value, tmp;
        std::istringstream wordstream(source);
        if ((wordstream >> tmp) && (tmp.compare(paramname) == 0)) {
            value = getValue(source);
        }
        return value;
    }
    
    /**
     * This convenience function extracts a value from a line found in the config file. Depending on what there is
     * to extract and if the line meets the syntax requirements or not, the return value may be an empty string.
     * @param source The line as taken from the config file
     * @param final_delim The end-of-line marker used in this line
     * @param delimiter The marker denoting the beginning of the substring that is to be extracted
     * @return The substring containing the requested value
     */
    std::string MatParser::getValue(std::string source, bool final_delim, std::string delimiter) {
        std::string value;
        unsigned int start = source.find(delimiter) + delimiter.size();
        if (start < source.size()) {
            if (final_delim) {
                unsigned int stop = source.find(line_end_delim);
                if (stop != std::string::npos) value = source.substr(start, stop - start);
            }
            else value = source.substr(start);
        }
        balgo::trim(value);
        return value;
    }
    
    /**
     * This convenience function turns a given string into a material unit of type <i>Matunit</i>. If the
     * contents of the string cannot be mapped to a unit, an exception is thrown.
     * @param source The string that needs to be converted
     * @return The unit described by the input string
     */
    MatCalc::Matunit MatParser::getUnit(std::string source) { // throws exception
        if (source.compare(gr_unit) == 0) return MatCalc::gr;
        else if(source.compare(mm_unit) == 0) return MatCalc::mm;
        else if(source.compare(mm3_unit) == 0) return MatCalc::mm3;
        else if(source.compare(grpm_unit) == 0) return MatCalc::grpm;
        else throw std::range_error(msg_unknown_unit);
    }
}
