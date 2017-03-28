#ifndef ANALYZERTOOLS_HH
#define ANALYZERTOOLS_HH

#include "global_constants.hh"

#include <TH2D.h>
#include <string>

namespace insur {

class AnalyzerTools {
protected:
  static void prepareTrackerMap(TH2D& myMap, const std::string& name, const std::string& title);
  static void prepareRadialTrackerMap(TH2D& myMap, const std::string& name, const std::string& title);
};

}

#endif
