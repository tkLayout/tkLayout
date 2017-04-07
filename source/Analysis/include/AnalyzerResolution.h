/*
 * AnalyzerResolution.h
 *
 *  Created on: 20. 4. 2016
 *      Author: Drasal (CERN)
 */

#ifndef INCLUDE_ANALYZERRESOLUTION_H_
#define INCLUDE_ANALYZERRESOLUTION_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <AnalyzerUnit.h>
#include "global_constants.h"
#include "Track.h"

class BarrelModule;
class CsvTextBuilder;
class Detector;
class EndcapModule;
class RootWContent;
class RootWPage;
class TProfile;

/*
 * @class AnalyzerResolution
 * Analyze tracker resolution for various combination of tracker layers/discs, tagged with given name: pixel (pixel only),
 * strip (strip only), fwd (forward only), etc. The tag is specified for each tracker component by a user in the geometry
 * configuration file. Data are then vizualized and printed out in a html formatted output.
 * Unique module name defined as "AnalyzerResolution".
 */
class AnalyzerResolution : public AnalyzerUnit {

 public:

  //! Constructor
  AnalyzerResolution(const Detector& detector);

  //! Destructor
  virtual ~AnalyzerResolution() {};

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

  //! Get Csv text output for const pt -> exception thrown if doesn't exist
  const CsvTextBuilder& getCsvResPt() const;

  //! Get Csv text output for const p -> exception thrown if doesn't exist
  const CsvTextBuilder& getCsvResP() const;

  //! Get Csv text output for hit collection-> exception thrown if doesn't exist
  const CsvTextBuilder& getCsvHitCol() const;

 private:

  //! Prepare plot: fill with data & set properties; varType specifies variable type to be filled: pT, p, d0, z0, phi0, cotgTheta & scenario const pT,  const p
  void preparePlot(std::vector<unique_ptr<TProfile>>& profHisArray, std::string varType, std::string scenario, const std::map<int, TrackCollection>& mapCollection);

  //! Prepare summary content table
  void prepareSummaryTable(std::string tag, std::string scenario, RootWPage& webPage, RootWContent& summaryContent, CsvTextBuilder& csvContainer);

  //! Calculate average values in defined regions for given profile histogram
  std::vector<double> averageHisValues(const TProfile& his, std::vector<double> regions);

  int    m_nTracks; //!< Number of simulation tracks to be used in the analysis
  double m_etaMin;  //!< Minimum eta value;
  double m_etaMax;  //!< Maximum eta value

  std::unique_ptr<CsvTextBuilder> m_csvResPt; //!< Csv containers -> keep final resolution for given pt in csv format
  std::unique_ptr<CsvTextBuilder> m_csvResP;  //!< Csv containers -> keep final resolution for given p in csv format

  std::map<std::string, std::map<int, TrackCollection>> m_taggedTrackPtCollectionMap;      //!< For given track tag -> map of track collection for given pT (full material)
  std::map<std::string, std::map<int, TrackCollection>> m_taggedTrackPtCollectionMapIdeal; //!< For given track tag -> map of track collection for given pT (ideal - no material)
  std::map<std::string, std::map<int, TrackCollection>> m_taggedTrackPCollectionMap;       //!< For given track tag -> map of track collection for given p (full material)
  std::map<std::string, std::map<int, TrackCollection>> m_taggedTrackPCollectionMapIdeal;  //!< For given track tag -> map of track collection for given p (ideal - no material)

  const double c_max_dPtOverPt = 1000;  // [%]
  const double c_min_dPtOverPt = 0.001; // [%]
  const double c_max_dZ0       = 5000.;
  const double c_min_dZ0       = 1.;
  const double c_max_dD0       = 5000.;
  const double c_min_dD0       = 1.;
  const double c_max_dPhi0     = 100.;
  const double c_min_dPhi0     = 1E-4;
  const double c_max_dCtgTheta = 1.0;
  const double c_min_dCtgTheta = 1E-6;
  const double c_min_dCTau     = 1.;
  const double c_max_dCTau     = 1E3;

  std::unique_ptr<CsvTextBuilder> m_csvHitCol;

  const int    c_nBins;

}; // Class

#endif /* INCLUDE_ANALYZERRESOLUTION_H_ */
