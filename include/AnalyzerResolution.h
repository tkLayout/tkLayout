/*
 * AnalyzerResolution.h
 *
 *  Created on: 20. 4. 2016
 *      Author: Drasal (CERN)
 */

#ifndef INCLUDE_ANALYZERRESOLUTION_H_
#define INCLUDE_ANALYZERRESOLUTION_H_

#include <string>
#include <vector>

#include "AnalyzerModule.h"

class Tracker;

class AnalyzerResolution : AnalyzerModule {

 public:
  // Constructor
  AnalyzerResolution(std::vector<const Tracker*> trackers);
  // Destructor
  ~AnalyzerResolution() {};

  bool init(int nTracks) {return true;}

  // Calculate track resolution
  bool analyze();

  // Visualize - add html page with all calculations & results
  bool visualize(RootWSite& webSite);

}; // Class

#endif /* INCLUDE_ANALYZERRESOLUTION_H_ */
