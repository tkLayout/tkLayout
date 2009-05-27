#ifndef _GLOBAL_CONSTANTS_H
#define	_GLOBAL_CONSTANTS_H

/**
 * @file global_constants.h
 * @brief This header file lists the global constants of the application.
 */

#include <string>

namespace insur {
    static const double PI = 3.14159265358979323846;
    static const double E = 2.71828182845904523536;
    /**
     * Geometry constants; all length measurements are in mm
     * @param epsilon The standard distance between one solid object and the next
     * @param volume_width The standard geometrical thickness of an inactive volume
     * @param inner_radius The inner radius of the tracker; the inner support tube starts immediately inside, everything below is part of the pixel detector
     * @param outer_radius The outer radius of the tracker; the outer support tube starts immediately outside, everything above is part of ECAL
     * @param max_length The maximum length, -z to +z, available to place the tracker components
     */
    static const double epsilon = 0.1;
    static const double volume_width = 10.0;
    static const double inner_radius = 214.0;
    static const double outer_radius = 1190.0;
    static const double max_length = 2910.0;
    
    /**
     * Visualisation constants: material parameters for active surfaces, services and supports, plus top volume padding.
     * The selected materials are completely arbitrary and only meant to distinguish one type of surface from another visually.
     * @param a_silicon Silicon atomic number
     * @param z_silicon Silicon charge
     * @param d_silicon Silicon density
     * @param a_copper Copper atomic number
     * @param z_copper Copper charge
     * @param d_copper Copper density
     * @param a_carbon Carbon atomic number
     * @param z_carbon Carbon charge
     * @param d_carbon Carbon density
     * @param top_volume_pad The extra space that is added to the dimensions of the top volume on each side of the cube
     */
    static const double a_silicon = 28.0855;
    static const double z_silicon = 14;
    static const double d_silicon = 2.329;
    static const double a_copper = 63.546;
    static const double z_copper = 29;
    static const double d_copper = 8.96;
    static const double a_carbon = 12.0107;
    static const double z_carbon = 6;
    static const double d_carbon = 1.9;
    static const double top_volume_pad = 200;
    
    /**
     * Internal string constants
     */
    static const std::string type_rphi = "rphi";
    static const std::string type_stereo = "stereo";
    static const std::string type_pt = "pt";
    
    
    /**
     * Filename and path constants
     * @param default_mattabdir Relative path to the list of materials
     * @param default_mattabfile List of materials and of their properties as required by <i>MaterialTable</i>
     * @param default_rootfiledir Output directory for <i>.root</i> files
     * @param default_rootfile Default filename for <i>.root</i> geometry output
     * @param default_graphdir Output directory for neighbour graphs
     * @param default_graphfile Default filename for neighbour graph output
     */
    static const std::string default_mattabdir = "config";
    static const std::string default_mattabfile = "mattab.list";
    static const std::string default_rootfiledir = "rootfiles";
    static const std::string default_rootfile = "trackergeometry.root";
    static const std::string default_graphdir = "graphs";
    static const std::string default_graphfile = "neighbours.graph";
    static const std::string default_summarypath = "matsum";
    static const std::string default_summary = "profiles.html";
}
#endif /* _GLOBAL_CONSTANTS_H */
