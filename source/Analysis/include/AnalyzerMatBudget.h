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

// Forward declaration
class BeamPipe;
class ConstGeometryVisitor;
class TCanvas;
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
  virtual ~AnalyzerMatBudget() {};

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

  int    m_nTracks;     //!< Number of geometry tracks to be used in material analysis
  double m_etaSpan;     //!< Eta interval to be analyzed
  double m_etaMin;      //!< Minimum eta value;
  double m_etaMax;      //!< Maximum eta value

  // Histogram output

}; // Class

#endif /* INCLUDE_ANALYZERMATBUDGET_H_ */
