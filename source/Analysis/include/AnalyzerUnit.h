/*
 * AnalyzerUnit.h
 *
 * Base analyzis and visualization class -> inherit from this class whenever you want
 * to analyze data and visualize them in a user defined way. One has to implement the
 * content of two purly virtual methods: analyze() and visualize().
 *
 *  Created on: 27. 11. 2015
 *      Author: drasal
 */

#ifndef INCLUDE_ANALYZERUNIT_H_
#define INCLUDE_ANALYZERUNIT_H_

#include <vector>
#include <string>

class BeamPipe;
class Detector;
class MaterialBudget;
class RootWSite;
class TCanvas;
class Tracker;

/*
 * @class AnalyzerUnit
 * @brief Pure virtual class to be used as a base class for concrete implementation of analyzer modules (units).
 * @details Pure virtual class to be used as a base class for concrete implementation of analyzer modules (units).
 * Call init method to initialize variables, histograms, ... -> update m_isInitOK if OK
 * Call analyze method to perform calculations. Apply analyze method only if init didn't fail -> update m_isAnalysisOK if OK
 * Call visualize method to perform visualization. Apply visualize method only if init & analyze didn't fail.
 */
class AnalyzerUnit
{
 public:

  //! Constructor - set active trackers only
  AnalyzerUnit(std::string name, const Detector& detector);

  //! Virtual destructor
  virtual ~AnalyzerUnit();

  //! Pure virtual initialization method -> use to initialize various variables, histograms, ...
  virtual bool init(int nTracks) = 0;

  //! Pure virtual analysis method -> analyzes data
  virtual bool analyze() = 0;

  //! Pure virtual visualization method -> visualizes output
  virtual bool visualize(RootWSite& webSite) = 0;

  //! Is unit correctly initialized
  bool isInitOK() {return m_isInitOK;}

  //! Is unit correctly analyzed
  bool isAnalysisOK() {return m_isAnalysisOK;}

  //! Is unit correctly visualized
  bool isVisOK() {return m_isVisOK;}

  //! Pure virtual get name -> returns unique name of a unit
  virtual std::string getName() final {return m_name;}

  //

 protected:

  //! Set canvas standard properties - background color, border mode, ...
  void setCanvasProperties(TCanvas& canvas);

  //! Is correctly initialized
  bool m_isInitOK;

  //! Is correctly analyzed
  bool m_isAnalysisOK;

  //! Is correctly visualized
  bool m_isVisOK;

  //! Unique name
  std::string m_name;

  //! Vector of sub-trackers -> const references, one can't and shouldn't change its content, nor delete the pointers
  std::vector<const Tracker*> m_trackers;

  //! Beam pipe -> const reference, one can't and shouldn't change its content, nor delete the pointers
  const BeamPipe* m_beamPipe;
};

#endif /* INCLUDE_ANALYZERUNIT_H_ */
