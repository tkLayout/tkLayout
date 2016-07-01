/*
 * AnalyzerResolution.h
 *
 *  Created on: 20. 4. 2016
 *      Author: Drasal (CERN)
 */

#ifndef INCLUDE_ANALYZERRESOLUTION_H_
#define INCLUDE_ANALYZERRESOLUTION_H_

#include <map>
#include <string>
#include <vector>

#include <AnalyzerUnit.h>
#include <Track.h>
#include <Visitor.h>

class BarrelModule;
class BeamPipe;
class EndcapModule;
class Tracker;

/*
 * @class AnalyzerResolution
 */
class AnalyzerResolution : public AnalyzerUnit {

 public:

  //! Constructor
  AnalyzerResolution(std::vector<const Tracker*> trackers, const BeamPipe* beamPipe);

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

 private:

  int    m_nTracks; //!< Number of simulation tracks to be used in the analysis
  double m_etaMin;  //!< Minimum eta value;
  double m_etaMax;  //!< Maximum eta value

  std::map<std::string, std::map<int, TrackCollection>> taggedTrackPtCollectionMap;      //!< For given track tag -> map of track collection for given pT (full material)
  std::map<std::string, std::map<int, TrackCollection>> taggedTrackPtCollectionMapIdeal; //!< For given track tag -> map of track collection for given pT (ideal - no material)
  std::map<std::string, std::map<int, TrackCollection>> taggedTrackPCollectionMap;       //!< For given track tag -> map of track collection for given p (full material)
  std::map<std::string, std::map<int, TrackCollection>> taggedTrackPCollectionMapIdeal;  //!< For given track tag -> map of track collection for given p (ideal - no material)

}; // Class

//! Helper class - Material track visitor (visitor pattern) - checks that modules, beam-pipe etc. hit by a track -> returns radiation & interaction length
class MatTrackVisitor : public ConstGeometryVisitor {

public:

  //! Constructor
  //! As a parameter provide a track, under which the material is studied
  MatTrackVisitor(Track& matTrack);

  //! Destructor
  ~MatTrackVisitor();

  //! Visit BeamPipe
  void visit(const BeamPipe& bp) override;

  //! Visit BarrelModule (no limits on Rods, Layers or Barrels)
  void visit(const BarrelModule& m) override;

  //! Visit EndcapModule (no limits on Rings, Discs or Endcaps)
  void visit(const EndcapModule& m) override;

private:

  //! Analyze if module crossed by given track & how much material is on the way
  void analyzeModuleMB(const DetectorModule& m);

  Track& m_matTrack;   //!< Shooting direction + origin encapsulated in a track class -> update track with hits once found

}; // Helper class

#endif /* INCLUDE_ANALYZERRESOLUTION_H_ */
