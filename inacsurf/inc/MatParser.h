// 
// File:   MatParser.h
// Author: ndemaio
//
// Created on January 13, 2009, 12:27 PM
//

/**
 * @file MatParser.h
 * @brief This is the header file for the material list parser
 */

#ifndef _MATPARSER_H
#define	_MATPARSER_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <MatCalc.h>
#include <MaterialTable.h>
#include <global_constants.h>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/trim.hpp>
/**
 * A shorter alias for the filesystem library namespace
 */
namespace bfs = boost::filesystem;
namespace balgo = boost::algorithm;
namespace insur {
    /*
     * String constants that will be used during parsing.
     */
    static const std::string c_comment = "//";
    static const std::string shell_comment = "#";
    static const std::string dummy_value = "-1.0";
    static const std::string gr_unit = "g";
    static const std::string mm_unit = "mm";
    static const std::string mm3_unit = "mm3";
    static const std::string grpm_unit = "g/m";
    static const std::string local_marker = "L";
    static const std::string exit_marker = "E";
    static const std::string name_value_delim = "=";
    static const std::string m_line_delim = "M";
    static const std::string s_line_delim = "S";
    static const std::string d_line_delim = "D";
    static const std::string w_line_delim = "W";
    static const std::string x_line_delim = "X";
    static const std::string y_line_delim = "Y";
    static const std::string z_line_delim = "Z";
    static const std::string line_end_delim = ";";
    static const std::string type_marker = "type";
    static const std::string strip_marker = "nStripsAcross";
    static const std::string seg_marker = "nSegments";
    
    
    /*
     * Error messages that may be reported.
     */
    static const std::string msg_no_mat_file = "Material file does not exist.";
    static const std::string msg_invalid_type = "Invalid type description found. Skipping line.";
    static const std::string msg_unknown_type = "Unknown type description found.";
    static const std::string msg_unexpected_type = "Unexpected type description found.";
    static const std::string msg_unknown_unit = "Unknown unit found.";
    static const std::string msg_readparam_failure = "Error reading parameter.";
    static const std::string msg_m_line_err = "Error parsing M entry. Skipping line.";
    static const std::string msg_s_line_err = "Error parsing S entry. Skipping line.";
    static const std::string msg_d_line_err = "Error parsing D entry. Skipping line.";
    static const std::string msg_x_line_err = "Error parsing X entry. Skipping line";
    static const std::string msg_y_line_err = "Error parsing Y entry. Skipping line";
    static const std::string msg_z_line_err = "Error parsing Z entry. Skipping line";
    static const std::string msg_w_line_err = "Error parsing W entry. Skipping line";
    static const std::string warning_too_many_values = "Warning: ignoring values beyond the fourth.";
    
    /**
     * @class MatParser
     * @brief This class provides a parser for the - very simple - material config file format
     *
     * It encapsulates the parser code that reads the material list from file and initialises the material table.
     */
    class MatParser {
    public:
        MatParser();
        ~MatParser();
        bool fillTable(std::string materialfile, MaterialTable& mattab);
        bool readParameters(std::string configfile, MatCalc& calc);
        bool initMatCalc(std::string configfile, MatCalc& calc);
    protected:
        bool parseStripsSegs(std::ifstream& instream, std::string& strips, std::string& segs);
        bool parseMLine(std::string line, std::string type, MatCalc& calc);
        bool parseDLine(std::string line, MatCalc& calc);
        bool parseSimpleLine(std::string line, MatCalc& calc, std::string marker);
    private:
        std::string readFromLine(std::string source, std::string paramname);
        std::string getValue(std::string source, bool final_delim = true, std::string delimiter = name_value_delim);
        MatCalc::Matunit getUnit(std::string source); // throws exception
    };
}
#endif	/* _MATPARSER_H */

