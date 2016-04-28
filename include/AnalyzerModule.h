/*
 * AnalyzerModule.h
 *
 * Base analyzis and visualization class -> inherit from this class whenever you want
 * to analyze data and visualize them in a user defined way. One has to implement the
 * content of two purly virtual methods: analyze() and visualize().
 *
 *  Created on: 27. 11. 2015
 *      Author: drasal
 */

#ifndef INCLUDE_ANALYZERMODULE_H_
#define INCLUDE_ANALYZERMODULE_H_

#include <vector>
#include <string>

class RootWSite;
class Tracker;
class MaterialBudget;

/*
 * Pure virtual class to be used for Analyzer modules.
 * Call init method to initialized variables, histograms, ... -> update m_isInitOK if OK
 * Call analyze method to perform calculations. Apply analyze method only if init didn't fail -> update m_isAnalysisOK if OK
 * Call visualize method to perform visualization. Apply visualize method only if init & analyze didn't fail.
 */
class AnalyzerModule
{
 public:

  //! Constructor - set active trackers to be analyzed
  AnalyzerModule(std::string name, std::vector<const Tracker*> trackers);

  //! Virtual destructor
  virtual ~AnalyzerModule();

  //! Pure virtual initialization method -> use to initialize various variables, histograms, ...
  virtual bool init(int nTracks) = 0;

  //! Pure virtual analysis method -> analyzes data
  virtual bool analyze() = 0;

  //! Pure virtual visualization method -> visualizes output
  virtual bool visualize(RootWSite& webSite) = 0;

  //! Pure virtual get name
  virtual std::string getName() final {return m_name;}

 protected:

  //! Is correctly initialized
  bool m_isInitOK;

  //! Is correctly analyzed
  bool m_isAnalysisOK;

  //! Unique name
  std::string m_name;

  //! Vector of active trackers -> const pointer, one can't and shouldn't change its content
  std::vector<const Tracker*> m_trackers;
};

#endif /* INCLUDE_ANALYZERMODULE_H_ */
