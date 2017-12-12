#include <AnalyzerTools.hh>

namespace insur {

void AnalyzerTools::prepareTrackerMap(TH2D& myMap, const std::string& name, const std::string& title) { 
  int mapBinsY = int( (geom_max_radius + geom_inactive_volume_width) * geom_safety_factor / 10.); // every cm
  int mapBinsX = int( (geom_max_length) * geom_safety_factor / 10.); // every cm
  myMap.SetName(name.c_str());
  myMap.SetTitle(title.c_str());
  myMap.SetXTitle("z [mm]");
  myMap.SetYTitle("r [mm]");
  myMap.SetBins(mapBinsX, 0.0, geom_max_length*geom_safety_factor, mapBinsY, 0.0, (geom_max_radius + geom_inactive_volume_width) * geom_safety_factor);
  myMap.Reset();
}

void AnalyzerTools::prepareRadialTrackerMap(TH2D& myMap, const std::string& name, const std::string& title) { 
  int mapBinsY = int( (2*geom_max_radius) * geom_safety_factor / 10.); // every cm
  int mapBinsX = int( (2*geom_max_radius) * geom_safety_factor / 10.); // every cm
  myMap.SetName(name.c_str());
  myMap.SetTitle(title.c_str());
  myMap.SetXTitle("x [mm]");
  myMap.SetYTitle("y [mm]");
  myMap.SetBins(mapBinsX, -geom_max_radius*geom_safety_factor, geom_max_radius*geom_safety_factor, mapBinsY, -geom_max_radius*geom_safety_factor, geom_max_radius*geom_safety_factor);
  myMap.Reset();
}


}
