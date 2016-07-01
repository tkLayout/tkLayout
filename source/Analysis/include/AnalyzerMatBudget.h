/*
 * AnalyzerMatBudget.h
 *
 *  Created on: 9. 6. 2016
 *      Author: drasal
 */

#ifndef INCLUDE_ANALYZERMATBUDGET_H_
#define INCLUDE_ANALYZERMATBUDGET_H_

// System libraries
#include <string>
#include <vector>
#include <map>

#include <Visitor.h>
#include <AnalyzerUnit.h>
#include <CsvTextBuilder.h>
#include <DetectorModule.h>

// Forward declaration
class BarrelModule;
class BeamPipe;
class ConstGeometryVisitor;
class Disk;
class EndcapModule;
class ModuleCap;
class RILength;
class TCanvas;
class TGraph;
class TH1D;
class TH2D;
class Track;
class Tracker;
class RootWTable;

/*
 * @class AnalyzerMatBudget
 * Analyze material budget of given layout, vizualize data and print them out in a html formatted output.
 * Unique name defined as "AnalyzerMatBudget".
 */
class AnalyzerMatBudget : public AnalyzerUnit {

 public:

  //! Constructor
  AnalyzerMatBudget(std::vector<const Tracker*> trackers, const BeamPipe* beamPipe);

  //! Destructor
  virtual ~AnalyzerMatBudget();

  //! Initialize - mostly histograms & other containers
  //! @return True if OK
  bool init(int nMatTracks);

  //! Inspect geometry layout (if init OK) -> collect data to histograms & tables
  //! @return True if OK
  bool analyze();

  //! Visualize geometry layout (if init & analysis OK) -> add html page with collected tables & created histograms
  //! @return True if OK
  bool visualize(RootWSite& webSite);

  //! Get number of used material tracks
  int getNMatTracks() const { return m_nTracks;}

 private:

  //! Helper method creating profile histogram from histogram
  //TProfile* createProfileHis(const TH1D& sourceHistogram) const;

  //! Helper method calculating average value in the given cut regions
  std::vector<double> average(TGraph& myGraph, std::vector<double> cuts);

  //! Helper method modifying a TGraph, so that it looks like a histogram (can be filled)
  void closeGraph(TGraph& myGraph);

  //! Helper method calculating average histogram values in histogram range from lowest bin to cutoff
  double averageHistogramValues(const TH1D& his, double cutoff) const;

  //! Helper method calculating average histogram values in range cutoffStart - cutoffEnd
  double averageHistogramValues(const TH1D& his, double cutoffStart, double cutoffEnd) const;

  int    m_nTracks;     //!< Number of geometry tracks to be used in material analysis
  double m_etaSpan;     //!< Eta interval to be analyzed
  double m_etaMin;      //!< Minimum eta value;
  double m_etaMax;      //!< Maximum eta value

  const float c_etaSafetyMargin = 0.01;

  // Histogram output -> trkName stands for beam-pipe or given tracker name, componentName stands for "Barrel", "Endcap", "Supports", "Services" or "Total"
  std::map<std::string, std::map<std::string, TH1D>> m_radMB; //!< Individual contributions to material budget [in rad. lengths] for given tracker: m_radMBComponents[trkName][componentName] = TH1D
  std::map<std::string, std::map<std::string, TH1D>> m_intMB; //!< Individual contributions to material budget [in int. lengths] for given tracker: m_intMBComponents[trkName][componentName] = TH1D

  // Component histogram output -> trkName stands for beam-pipe or given tracker name, componentName stands for any defined module component (given in geom. config), "Supports", "Services" and "Total"
  std::map<std::string, std::map<std::string, TH1D>> m_radMBComp; //!< Individual contributions to material budget [in rad. lengths] for given tracker: m_radMBComponents[trkName][componentName] = TH1D
  std::map<std::string, std::map<std::string, TH1D>> m_intMBComp; //!< Individual contributions to material budget [in int. lengths] for given tracker: m_intMBComponents[trkName][componentName] = TH1D

  // Material maps for given trackers: "Total" corresponds to sum
  std::map<std::string, TH2D> m_radMap;      //!< Material map in radiation lengths for given tracker: m_radMap[trkName] = TH2D
  std::map<std::string, TH2D> m_radMapCount; //!< Counter for material map in terms of radiation length
  std::map<std::string, TH2D> m_intMap;      //!< Material map in interaction lengths for given tracker: m_radMap[trkName] = TH2D
  std::map<std::string, TH2D> m_intMapCount; //!< Counter for material map in terms of interaction length

  // Hadrons -> each variable for given tracker: "Total" corresponds to sum
  std::map<std::string, TGraph> m_hadronTotalHitsGraph;
  std::map<std::string, TGraph> m_hadronAverageHitsGraph;
  std::map<std::string, std::vector<double>> m_hadronNeededHitsFraction;
  std::map<std::string, std::vector<TGraph>> m_hadronGoodTracksFraction;

  // Csv output
  CsvTextBuilder       m_materialCsv;

}; // Class

//! Helper class - Material budget visitor (visitor pattern) - checks that modules, beam-pipe etc. hit by a track -> returns radiation & interaction length & preanalyze material related data
class MatBudgetVisitor : public ConstGeometryVisitor {

public:

  //! Constructor
  //! As a parameter provide a track, under which the material is studied, container filling in
  //! material for given subdetectors or components defined by name (Barrel, Endcap, ...) and
  //! radiation maps as TH2D to be filled
  MatBudgetVisitor(Track& matTrack, std::map<std::string, RILength>& matBudget, TH2D& radMap, TH2D& radMapCount, TH2D& intMap, TH2D& intMapCount);

  //! Destructor
  ~MatBudgetVisitor();

  //! Visit BeamPipe
  void visit(const BeamPipe& bp) override;

  //! Visit BarrelModule (no limits on Rods, Layers or Barrels)
  void visit(const BarrelModule& m) override;

  //! Visit EndcapModule (no limits on Rings, Discs or Endcaps)
  void visit(const EndcapModule& m) override;

  //! Post visit method called after the whole visit-accept pattern done to recalibrate histograms, etc.
  void postvisit();

private:

  //! Analyze if module crossed by given track & how much material is on the way
  void analyzeModuleMB(const DetectorModule& m);

  std::map<std::string, RILength>& m_matBudget;  //!< Material container for given subdetectors or components defined by name (Barrel, Endcap, ...)
  TH2D&                            m_radMap;     //!< Material map in terms of radiation length
  TH2D&                            m_radMapCount;//!< Counter for material map in terms of radiation length
  TH2D&                            m_intMap;     //!< Material map in terms of interaction length
  TH2D&                            m_intMapCount;//!< Counter for material map in terms of interaction length
  Track&                           m_matTrack;   //!< Shooting direction + origin encapsulated in a track class -> update track with hits once found
  int                              m_nEntries;   //!< Number of entries

}; // Helper class

#endif /* INCLUDE_ANALYZERMATBUDGET_H_ */
