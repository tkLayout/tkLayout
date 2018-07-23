#ifndef ANALYZERTOOLS_HH
#define ANALYZERTOOLS_HH

#include "global_constants.hh"

#include <TH2D.h>
#include <string>

namespace insur {
  // const specific to the Analyzer

  /**
   * Geometry coverage plots constants.
   **/

  // If there are more hits (or stubs) than the specified max, hits (or stubs) are counted in the max category.
  // Example: if plotNumberOfInnerTrackerStubs == 3, and a track has 4 stubs, it will be counted in the (>=3 stubs) category.
  static const int plotMaxNumberOfOuterTrackerHitsPerLayer = 5;
  static const int plotMaxNumberOfInnerTrackerHitsPerLayer = 4;
  static const int plotMaxNumberOfOuterTrackerStubs = 11;
  static const int plotMaxNumberOfInnerTrackerStubs = 3;

  // Simply set the plots MaxY.
  static const double plotNumberOfOuterTrackerHitModulesMaxY = 10.;
  static const double plotNumberOfInnerTrackerHitModulesMaxY = 20.;
  static const double plotNumberOfHitsMaxY = 20.;
  static const double plotNumberOfStubsMaxY = 10.;
  static const double plotNumberOfOuterTrackerHitLayersMaxY = 8.;
  static const double plotNumberOfInnerTrackerHitLayersMaxY = 13.;
  static const double plotNumberOfTriggeredPointsMaxY = 10.;

  /**
   * Warnings.
   **/
  static const std::string msg_module_warning = "Warning: tracker module with undefined subdetector type found.";


  /**
   * Categories.
   **/
  static const std::string outer_tracker_id = "Outer";
  static const std::string beam_pipe = "Beam Pipe";
  static const std::string services_under_pixel_tracking_volume = "Services under Pixel Tracking Volume";
  static const std::string services_in_pixel_tracking_volume = "Services in Pixel Tracking Volume";
  static const std::string supports_in_pixel_tracking_volume = "Supports in Pixel Tracking Volume";
  static const std::string services_and_supports_in_interstice = "Services and supports in interstice";
  static const std::string services_in_outer_tracking_volume = "Services in Outer Tracking Volume";
  static const std::string supports_in_outer_tracking_volume = "Supports in Outer Tracking Volume";


class AnalyzerTools {
protected:
  static void prepareTrackerMap(TH2D& myMap, const std::string& name, const std::string& title);
  static void prepareRadialTrackerMap(TH2D& myMap, const std::string& name, const std::string& title);
};

}

#endif
