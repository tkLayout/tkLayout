#include <MatParser.h>
namespace insur {
    MatParser::MatParser() {}
    
    MatParser::~MatParser() {}
    
    bool MatParser::fillTable(std::string materialfile, MaterialTable& mattab) {
        if (materialfile.empty()) materialfile = default_mattabfile;
        bfs::path mpath(materialfile);
        if (bfs::exists(mpath)) {
            try {
                std::string line, word;
                std::ifstream infilestream;
                infilestream.open(materialfile.c_str());
                while (std::getline(infilestream, line)) {
                    std::istringstream wordstream(line);
                    std::vector<std::string> tmp;
                    MaterialRow row;
                    while (wordstream >> word) {
                        if ((word.compare(0, c_comment.size(), c_comment) == 0) || (word.compare(0, shell_comment.size(), shell_comment) == 0)) break;
                        else tmp.push_back(word);
                    }
                    if (tmp.empty()) continue;
                    else {
                        if (tmp.size() < 4) {
                            while (tmp.size() < 4) tmp.push_back(dummy_value);
                        }
                        if (tmp.size() > 4) std::cerr << warning_too_many_values << std::endl;
                        row.tag = tmp.at(0);
                        row.density = atof(tmp.at(1).c_str());
                        row.rlength = atof(tmp.at(2).c_str());
                        row.ilength = atof(tmp.at(3).c_str());
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
        else std::cerr << msg_no_mat_file << std::endl;
        return false;
    }
    
    bool MatParser::readParameters(std::string configfile, MatCalc& calc) {
        bfs::path mpath(configfile);
        if (bfs::exists(mpath)) {
            try {
                std::string line, word, type;
                std::ifstream infilestream;
                infilestream.open(configfile.c_str());
                while (std::getline(infilestream, line)) {
                    balgo::trim(line);
                    std::istringstream wordstream(line);
                    wordstream >> word;
                    if ((word.compare(0, c_comment.size(), c_comment) == 0)
                            || (word.compare(0, shell_comment.size(), shell_comment) == 0)) continue;
                    else if (word.compare(type_marker) == 0) {
                        word = getValue(line, false);
                        if (word.compare(type_rphi) == 0) {
                            if (type.empty() || ((calc.registeredTypes() == 1) && (type.compare(type_pt)))) {
                                type = word;
                                std::string strips, segs;
                                if (parseStripsSegs(infilestream, strips, segs)) {
                                    int chips, segments;
                                    chips = atoi(strips.c_str());
                                    segments = atoi(segs.c_str());
                                    calc.addTypeInfo(MatCalc::rphi, chips, segments);
                                }
                                else std::cerr << msg_readparam_failure << std::endl;
                            }
                            else std::cerr << msg_unexpected_type << std::endl;
                        }
                        else if (word.compare(type_stereo) == 0) {
                            if (((calc.registeredTypes() > 1) && (type.compare(type_pt) == 0)) || (type.compare(type_rphi) == 0)) {
                                type = word;
                                calc.copyContents(MatCalc::rphi, MatCalc::stereo);
                                calc.addTypeInfo(MatCalc::stereo, calc.getStripsAcross(MatCalc::rphi), calc.getSegmentsAlong(MatCalc::rphi));
                            }
                            else std::cerr << msg_unexpected_type << std::endl;
                        }
                        else if (word.compare(type_pt) == 0) {
                            if (!calc.typeRegistered(MatCalc::pt)) {
                                type = word;
                                std::string strips, segs;
                                if (parseStripsSegs(infilestream, strips, segs)) {
                                    int chips, segments;
                                    chips = atoi(strips.c_str());
                                    segments = atoi(segs.c_str());
                                    calc.addTypeInfo(MatCalc::pt, chips, segments);
                                }
                                else std::cerr << msg_readparam_failure << std::endl;
                            }
                            else std::cerr << msg_unexpected_type << std::endl;
                        }
                        else {
                            std::cerr << msg_unknown_type << std::endl;
                            calc.clearTypeVector();
                            calc.clearModVectors();
                            return false;
                        }
                    }
                    else if (word.compare(m_line_delim) == 0) {
                        if (!parseMLine(line, type, calc)) std::cout << msg_m_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    else if(word.compare(d_line_delim) == 0) {
                        if (!parseDLine(line, calc)) std::cout << msg_d_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    else if (word.compare(s_line_delim) == 0) {
                        if (!parseSimpleLine(line, calc, s_line_delim)) std::cout << msg_s_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    else if (word.compare(x_line_delim) == 0) {
                        if (!parseSimpleLine(line, calc, x_line_delim)) std::cout << msg_x_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    else if (word.compare(y_line_delim) == 0) {
                        if (!parseSimpleLine(line, calc, y_line_delim)) std::cout << msg_y_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    else if (word.compare(z_line_delim) == 0) {
                        if (!parseSimpleLine(line, calc, z_line_delim)) std::cout << msg_z_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    else if (word.compare(w_line_delim) == 0) {
                        if (!parseSimpleLine(line, calc, w_line_delim)) std::cout << msg_w_line_err << std::endl << "Line was: " << line << std::endl;
                    }
                    else if (!line.empty()) {
                        std::cerr << "Confusion detected...skipping line '" << line << "'." << std::endl;
                        break;
                    }
                    word.clear();
                }
                infilestream.close();
                return true;
            }
            catch (bfs::filesystem_error& bfe) {
                std::cerr << bfe.what() << std::endl;
                calc.clearTypeVector();
                calc.clearModVectors();
                return false;
            }
            catch (std::bad_alloc& ba) {
                std::cerr << ba.what() << std::endl;
                calc.clearTypeVector();
                calc.clearModVectors();
                return false;
            }
        }
        else std::cerr << msg_no_mat_file << std::endl;
        calc.clearTypeVector();
        calc.clearModVectors();
        return false;
    }
    
    bool MatParser::initMatCalc(std::string configfile, MatCalc& calc) {
        std::string filename(default_mattabdir + "/" + default_mattabfile);
        if(fillTable(filename, calc.getMaterialTable())) {
            calc.initDone(readParameters(configfile, calc));
            if (calc.initDone()) return true;
            return false;
        }
        return false;
    }
    
// protected
    bool MatParser::parseStripsSegs(std::ifstream& instream, std::string& strips, std::string& segs) {
        bool strips_done = false, segs_done = false;
        std::string line;
        while (!(strips_done && segs_done) && std::getline(instream, line)) {
            balgo::trim(line);
            if (line.empty()) continue;
            else if((line.compare(0, c_comment.size(), c_comment) == 0) || (line.compare(0, shell_comment.size(), shell_comment) == 0)) continue;
            else if(line.compare(0, strip_marker.size(), strip_marker) == 0) {
                strips = readFromLine(line, strip_marker);
                strips_done = true;
            }
            else if (line.compare(0, seg_marker.size(), seg_marker) == 0) {
                segs = readFromLine(line, seg_marker);
                segs_done = true;
            }
            else return false;
        }
        if (strips.empty() || segs.empty()) return false;
        return true;
    }
    
    bool MatParser::parseMLine(std::string line, std::string type, MatCalc& calc) {
        if (type.empty()) return false;
        uint start = line.find(m_line_delim) + m_line_delim.size();
        uint stop = line.find(line_end_delim);
        if ((start >= line.size()) || (stop == std::string::npos)) return false;
        line = line.substr(start, stop - start);
        balgo::trim(line);
        std::stringstream wordstream(line);
        std::string tag, tmp;
        if (!(wordstream >> tag)) return false;
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
        if (tmp.compare(local_marker) == 0) local = true;
        else if (tmp.compare(exit_marker) == 0) local = false;
        else return false;
        MatCalc::Modtype tp;
        if (type.compare(type_rphi) == 0) tp = MatCalc::rphi;
        else if (type.compare(type_stereo) == 0) tp = MatCalc::stereo;
        else if (type.compare(type_pt) == 0) tp = MatCalc::pt;
        else return false;
        calc.addModuleParameters(tag, tp, A, uA, B, uB, C, uC, D, uD, local);
        return true;
    }
    
    bool MatParser::parseDLine(std::string line, MatCalc& calc) {
        uint start = line.find(d_line_delim) +d_line_delim.size();
        uint stop = line.find(line_end_delim);
        if ((start >= line.size()) || (stop == std::string::npos)) return false;
        line = line.substr(start, stop - start);
        balgo::trim(line);
        std::stringstream wordstream(line);
        std::string inTag, outTag, tmp;
        double in, out;
        MatCalc::Matunit uIn, uOut;
        bool local;
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
        else {
            if (tmp.compare(local_marker) == 0) local = true;
            else if (tmp.compare(exit_marker) == 0) local = false;
            else return false;
        }
        calc.addServiceParameters(inTag, in, uIn, outTag, out, uOut, local);
        return true;
    }
    
    bool MatParser::parseSimpleLine(std::string line, MatCalc& calc, std::string marker) {
        uint start = line.find(marker) + marker.size();
        uint stop = line.find(line_end_delim);
        if ((start >= line.size()) || (stop == std::string::npos)) return false;
        line = line.substr(start, stop - start);
        balgo::trim(line);
        std::stringstream wordstream(line);
        std::string tag, tmp;
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
        if (marker.compare(s_line_delim) == 0) calc.addServiceParameters(tag, val, uni);
        else if (marker.compare(x_line_delim) == 0) calc.addSupportParameters(tag, val, uni, MaterialProperties::b_sup);
        else if (marker.compare(y_line_delim) == 0) calc.addSupportParameters(tag, val, uni, MaterialProperties::e_sup);
        else if (marker.compare(z_line_delim) == 0) calc.addSupportParameters(tag, val, uni, MaterialProperties::t_sup);
        else if(marker.compare(w_line_delim) == 0) calc.addSupportParameters(tag, val, uni, MaterialProperties::u_sup);
        else return false;
        return true;
    }
// private
    std::string MatParser::readFromLine(std::string source, std::string paramname) {
        std::string value, tmp;
        std::istringstream wordstream(source);
        if ((wordstream >> tmp) && (tmp.compare(paramname) == 0)) {
            value = getValue(source);
        }
        return value;
    }
    
    std::string MatParser::getValue(std::string source, bool final_delim, std::string delimiter) {
        std::string value;
        uint start = source.find(delimiter) + delimiter.size();
        if (start < source.size()) {
            if (final_delim) {
                uint stop = source.find(line_end_delim);
                if (stop != std::string::npos) value = source.substr(start, stop - start);
            }
            else value = source.substr(start);
        }
        balgo::trim(value);
        return value;
    }
    
    MatCalc::Matunit MatParser::getUnit(std::string source) { // throws exception
        if (source.compare(gr_unit) == 0) return MatCalc::gr;
        else if(source.compare(mm_unit) == 0) return MatCalc::mm;
        else if(source.compare(mm3_unit) == 0) return MatCalc::mm3;
        else if(source.compare(grpm_unit) == 0) return MatCalc::grpm;
        else throw std::range_error(msg_unknown_unit);
    }
}
