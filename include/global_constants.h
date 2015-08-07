#ifndef _GLOBAL_CONSTANTS_H
#define	_GLOBAL_CONSTANTS_H

/**
 * @file global_constants.h
 * @brief This header file lists the global constants of the application.
 */

#include <string>
#include <vector>
#include <Rtypes.h>

namespace insur {
  static const double PI = 3.14159265358979323846;
  static const double E = 2.71828182845904523536;
  static const double magnetic_field = 3.8; // Tesla; CMS magnet field strength

  /**
   * Geometry constants; all length measurements are in mm
   * @param epsilon The standard distance between one solid object and the next
   * @param volume_width The standard geometrical thickness of an inactive volume
   * @param inner_radius The inner radius of the tracker; the inner support tube starts immediately inside, everything below is part of the pixel detector
   * @param outer_radius The outer radius of the tracker; the outer support tube starts immediately outside, everything above is part of ECAL
   * @param max_length The maximum length, in +z, available to place the tracker components
   */
  static const std::vector<std::string> name_eta_regions  = {""   ,"TRK-BRL","TRK-ENDCAP","TRK-FWD","TRK-IFWD"}; // Name tracker eta regions
  static const std::vector<double>      range_eta_regions = {0.001,1.5      ,2.5         ,4.0      ,6.0       }; // Name tracker eta regions

  static const double epsilon                    = 0.1;
  static const double volume_width               = 10.0;   // mm
  static const double inner_radius               = 800.0;  // mm
  static const double outer_radius               = 2500.0; // mm
  static const double max_length                 = 8000.0; // mm
  static const double pixel_radius               = 25.0;   // mm
  static const double z_threshold_service_zigzag = 100.0;
  static const double max_eta_coverage           = range_eta_regions[range_eta_regions.size()-1];   // Tracking performed from step_eta_epsilon to max_eta_coverage in steps
  static const int    n_eta_regions              = range_eta_regions.size();                        // Tracking performed in the following Number of eta regions

  static const int    default_n_tracks           = 120;                       // Default number of tracks simulated (max_eta_coverage/default_n_tracks = etaStep)
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
   * @param temperature_levels The number of different colour levels in 2D histogram plots
   */
  static const double a_silicon       = 28.0855;
  static const double z_silicon       = 14;
  static const double d_silicon       = 2.329;
  static const double a_copper        = 63.546;
  static const double z_copper        = 29;
  static const double d_copper        = 8.96;
  static const double a_carbon        = 12.0107;
  static const double z_carbon        = 6;
  static const double d_carbon        = 1.9;
  static const double top_volume_pad  = 200;
  static const int temperature_levels = 512;

  /**
   * Display formatting parameters - eta ticks displayed with short step in range 0 - short_eta_coverage, with long step in range
   * short_eta_coverage - long_eta_coverage
   */
  static const double step_eta_short     = 0.2;
  static const double step_eta_long      = 0.5;
  static const double step_eta_epsilon   = 0.001;
  static const double short_eta_coverage = range_eta_regions[1];
  static const double trk_eta_coverage   = range_eta_regions[2];
  static const double long_eta_coverage  = range_eta_regions[3];

  static const double max_dPtOverPt      = 500;  // [%]
  static const double min_dPtOverPt      = 0.001; // [%]
  static const double max_dZ0            = 5000.;
  static const double min_dZ0            = 1.;
  static const double max_dD0            = 5000.;
  static const double min_dD0            = 1.;
  static const double max_dPhi           = 100.;
  static const double min_dPhi           = 1E-4;
  static const double max_dCtgTheta      = 1.0;
  static const double min_dCtgTheta      = 1E-6;


  static const int    min_canvas_sizeX   = 600;
  static const int    std_canvas_sizeX   = 900;
  static const int    max_canvas_sizeX   =1800;
  static const int    min_canvas_sizeY   = 600;
  static const int    std_canvas_sizeY   = 900;
  static const int    max_canvas_sizeY   =1800;

  static const double default_eta_step   = 0.1;
  static const int    default_n_bins     = max_eta_coverage/default_eta_step;  // Default number of bins in histogram from eta=0  to max_eta_coverage

  /**
   * Internal string constants for standard one-sided and specialised double-sided, rotated types
   */
  static const std::string type_rphi   = "rphi";
  static const std::string type_stereo = "stereo";


  /**
   * Filename and path constants
   * @param default_mattabdir Relative path to the list of materials
   * @param default_mattabfile List of materials and of their properties as required by <i>MaterialTable</i>
   * @param default_rootfiledir Output directory for <i>.root</i> files
   * @param default_rootfile Default filename for <i>.root</i> geometry output
   * @param default_graphdir Output directory for neighbour graphs
   * @param default_graphfile Default filename for neighbour graph output
   * @param default_summarypath Output directory for material summaries
   * @param default_summary Default filename root for material summary
   * @param default_xmlpath Output base directory for CMSSW XML output
   * @param default_xml Default subdirectory name for CMSSW XML output
   */
  // TODO: make sure the following constants are only used in
  // mainConfigHandler
  static const std::string default_mattabdir                     = "config";
  static const std::string default_mattabfile                    = "mattab.list";
  static const std::string default_irradiationdir                = "config";
  static const std::vector<std::string> default_irradiationfiles = {"irradiation.map", "irradiationPixel.map"};
  //static const std::string default_irradiationfile = "irradiation.map";
  static const std::string default_materialsdir                  = "config";
  static const std::string default_tracker_materials_file        = "Materials.cfg";
  static const std::string default_pixel_materials_file          = "PixelMaterials.cfg";
  static const std::string suffix_tracker_material_file          = "_Materials.cfg";
  static const std::string suffix_pixel_material_file            = "_Materials.cfg.pix";
  static const std::string suffix_geometry_file                  = ".cfg";
  static const std::string suffix_types_file                     = "_Types.cfg";
  static const std::string default_rootfiledir                   = "rootfiles";
  static const std::string default_rootfile                      = "trackergeometry.root";
  static const std::string default_graphdir                      = "graphs";
  static const std::string default_graphfile                     = "neighbours.graph";
  static const std::string default_summary                       = "profiles.html";
  static const std::string default_xmlpath                       = "xml";
  static const std::string default_xml                           = "tk2CMSSWxml";
  static const std::string default_styledir                      = "style";
  static const std::string default_configdir                     = "config";
  static const std::string default_stdincludedir                 = "stdinclude";
  static const std::string default_geometriesdir                 = "geometries";
  static const std::string default_htmldir                       = "results";

  static const std::string csv_separator = ",";
  static const std::string csv_eol = "\n";

}
#endif /* _GLOBAL_CONSTANTS_H */
