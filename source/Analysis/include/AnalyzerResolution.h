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

/*
 * @class AnalyzerResolution
 */
class AnalyzerResolution : public AnalyzerModule {

 public:

  //! Constructor
  AnalyzerResolution(std::vector<const Tracker*> trackers);

  //! Destructor
  ~AnalyzerResolution() {};

  //! Initialize - mostly histograms & other containers
  //! @return True if OK
  bool init(int nTracks);

  //! Inspect geometry layout (if init OK) -> collect data to histograms & tables
  //! @return True if OK
  bool analyze();

  //! Visualize geometry layout (if init & analysis OK) -> add html page with collected tables & created histograms
  //! @return True if OK
  bool visualize(RootWSite& webSite);

  //! Get number of used tracks
  int getNSimTracks() const { return m_nTracks;}

 private:

  int m_nTracks; //!< Number of simulation tracks to be used in the analysis

}; // Class

#endif /* INCLUDE_ANALYZERRESOLUTION_H_ */
