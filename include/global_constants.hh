#ifndef _GLOBAL_CONSTANTS_H
#define	_GLOBAL_CONSTANTS_H

/**
 * @file global_constants.h
 * @brief This header file lists the global constants of the application.
 */

#include <string>
#include <cmath>
#include <vector>

namespace insur {
  //static const double magnetic_field = 3.8;                // Tesla; CMS magnet field strength -> Use definition in SimParms file
  static const double boltzmann_constant = 8.6173303E-05;  // eV/K
  static const double celsius_to_kelvin = 273.15;          // T(K) = T(°C) + celsius_to_kelvin
  static const double siliconEffectiveBandGap = 1.21;      // eV. Used in the Hamburg model (effect of temperature on sensor leakage current).

  /**
   * Geometry constants; all length measurements are in mm
   * @param epsilon The standard distance between one solid object and the next
   * @param volume_width The standard geometrical thickness of an inactive volume
   * @param inner_radius The inner radius of the tracker; the inner support tube starts immediately inside, everything below is part of the pixel detector
   * @param outer_radius The outer radius of the tracker; the outer support tube starts immediately outside, everything above is part of ECAL
   * @param max_length The maximum length, in +z, available to place the tracker components
   */
  static const std::vector<std::string> geom_name_eta_regions  = {""   ,"C","I","F","VF","VVF"};  // Name tracker eta regions
  static const std::vector<double>      geom_range_eta_regions = {0.001,0.8,1.6,2.4 ,3.2 ,4.0  }; // Name tracker eta regions

  static const double geom_epsilon                    = 0.1;
  static const double geom_inactive_volume_width      = 2.0;     // mm
  static const double geom_conversion_station_width   = 2.0;     // mm
  static const double geom_inner_pixel_radius         = 30.0;    // mm
  static const double geom_inner_strip_radius         = 150.0;   // mm
  static const double geom_outer_strip_radius         = 1200.0;  // mm
  static const double geom_z_threshold_service_zigzag = 100.0;
  static const double geom_top_volume_pad             = 200;     // mm

  static const double geom_support_margin_bottom      = 1;       // mm
  static const double geom_support_margin_top         = 2;       // mm

  static const double geom_safety_factor              = 1.1;

  static const double geom_min_radius                 = geom_inner_pixel_radius;
  static const double geom_max_radius                 = geom_outer_strip_radius;
  static const double geom_max_length                 = 2800.0; // mm

  static const double geom_max_eta_coverage           = geom_range_eta_regions[geom_range_eta_regions.size()-1]; // Tracking performed from step_eta_epsilon to max_eta_coverage in steps
  static const int    geom_n_eta_regions              = geom_range_eta_regions.size();                           // Tracking performed in the following Number of eta regions

  static const int    default_n_tracks                = 100;                       // Default number of tracks simulated (max_eta_coverage/default_n_tracks = etaStep)

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
  static const double mat_a_silicon          = 28.0855;
  static const double mat_z_silicon          = 14;
  static const double mat_d_silicon          = 2.329;
  static const double mat_a_copper           = 63.546;
  static const double mat_z_copper           = 29;
  static const double mat_d_copper           = 8.96;
  static const double mat_a_carbon           = 12.0107;
  static const double mat_z_carbon           = 6;
  static const double mat_d_carbon           = 1.9;

  static const int    vis_temperature_levels = 512;

  /**
   * Display formatting parameters - eta ticks displayed with short step in range 0 - short_eta_coverage, with long step in range
   * short_eta_coverage - long_eta_coverage
   */
  static const double vis_step_eta_short     = 0.2;
  static const double vis_step_eta_long      = 0.5;
  static const double vis_step_eta_epsilon   = 0.001;
  static const double vis_short_eta_coverage = geom_range_eta_regions[1];
  static const double vis_trk_eta_coverage   = geom_range_eta_regions[2];
  static const double vis_long_eta_coverage  = geom_range_eta_regions[3];

  static const double vis_max_dPtOverPt      = 100;  // [%]
  static const double vis_min_dPtOverPt      = 0.1; // [%]
  static const double vis_max_dZ0            = 5000.;
  static const double vis_min_dZ0            = 1.;
  static const double vis_max_dD0            = 5000.;
  static const double vis_min_dD0            = 1.;
  static const double vis_max_dPhi           = 100.;
  static const double vis_min_dPhi           = 1E-4;
  static const double vis_max_dCtgTheta      = 1.0;
  static const double vis_min_dCtgTheta      = 1E-6;
  static const double vis_max_dL             = 1E3;
  static const double vis_min_dL             = 1E0;
  static const double vis_max_beta           = M_PI/2;
  static const double vis_min_beta           = 0;
  static const double vis_max_omega          = +M_PI;
  static const double vis_min_omega          = -M_PI;

  static const double vis_safety_factor      = geom_safety_factor;

  static const int    vis_min_canvas_sizeX   = 600;
  static const int    vis_std_canvas_sizeX   = 900;
  static const int    vis_max_canvas_sizeX   =1800;
  static const int    vis_min_canvas_sizeY   = 600;
  static const int    vis_std_canvas_sizeY   = 900;
  static const int    vis_max_canvas_sizeY   =1800;

  static const double vis_eta_step           = 0.1;
  static const double vis_material_eta_step  = 0.05;
  static const int    vis_n_bins             = geom_max_eta_coverage/vis_eta_step;  // Default number of bins in histogram from eta=0  to max_eta_coverage

  /**
   * Internal string constants for standard one-sided and specialised double-sided, rotated types
   */
  static const std::string type_rphi   = "rphi";
  static const std::string type_stereo = "stereo";

  /**
   * Internal string constants for modules types.
   */
  static const std::string type_pixel   = "pixel";
  static const std::string type_timing = "timing";

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
  static const std::string default_fluence_file                  = "fluencePer1ppCol.map";
  static const std::string default_materialsdir                  = "config";
  static const std::string default_tracker_materials_file        = "Materials.cfg";
  static const std::string default_pixel_materials_file          = "PixelMaterials.cfg";
  static const std::string suffix_tracker_material_file          = "_Materials.cfg";
  static const std::string suffix_pixel_material_file            = "_Materials.cfg.pix";
  static const std::string suffix_geometry_file                  = ".cfg";
  static const std::string suffix_types_file                     = "_Types.cfg";
  static const std::string default_detidschemesdir               = "config";
  static const std::string default_detidschemesfile              = "det_id_schemes.list";
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

  static const std::string csv_separator = ",";
  static const std::string csv_eol       = "\n";

}
#endif /* _GLOBAL_CONSTANTS_H */
