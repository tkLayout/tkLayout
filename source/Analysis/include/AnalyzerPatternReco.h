/*
 * AnalyzerPatternReco.h
 *
 *  Created on: 28. 11. 2016
 *      Author: drasal
 */

#ifndef INCLUDE_ANALYZERPATTERNRECO_H_
#define INCLUDE_ANALYZERPATTERNRECO_H_

#include <AnalyzerUnit.h>
#include "global_constants.h"

// Fwd declaration
class Detector;
class IrradiationMap;
class RootWSite;

/*
 * @class AnalyzerPatternReco
 * Analyze tracker pattern recognition capabilities for various combination of tracker layers/discs, tagged with given name: pixel (pixel only),
 * strip (strip only), fwd (forward only), etc. The tag is specified for each tracker component by a user in the geometry
 * configuration file. Data are then vizualized and printed out in a html formatted output. In order to simulate backround one loads in the occupancy
 * map.
 * Unique module name defined as "AnalyzerPatternReco".
 */
class AnalyzerPatternReco : public AnalyzerUnit {

 public:

  //! Constructor
  AnalyzerPatternReco(const Detector& detector);

  //! Destructor
  virtual ~AnalyzerPatternReco();

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

  //! Check that a file can be opened
  bool checkFile(const std::string& fileName, const std::string& filePath);

  int    m_nTracks; //!< Number of simulated tracks to be used in the analysis
  double m_etaMin;  //!< Minimum eta value;
  double m_etaMax;  //!< Maximum eta value

  IrradiationMap* m_chargedMap;   //!< Charged particles fluxes

  const int    c_nBins;

}; // Class

#endif /* INCLUDE_ANALYZERPATTERNRECO_H_ */
