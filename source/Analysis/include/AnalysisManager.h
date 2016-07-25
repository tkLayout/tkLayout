/*
 * AnalysisManager.h
 *
 *  Created on: 27. 4. 2016
 *      Author: drasal
 */

#ifndef INCLUDE_ANALYSISMANAGER_H_
#define INCLUDE_ANALYSISMANAGER_H_

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

class AnalyzerUnit;
class Detector;
class RootWSite;

/*
 * @class AnalysisManager
 * @brief The core analysis class, processing all analyzers to be used.
 * @details The core analysis class, processing all analyzers to be used. Each analyzer is identified by its unique
 * name and is called through initModule, analyzeModule & visualizeModule methods. Init method initializes parameters,
 * histograms etc. of given module (unit), analyze method performs the analysis and visualization outputs the results
 * to html format. All html web pages created by different modules (units) are collected in html format in a common web
 * container (main html page) and printed-out using makeWebSite() method. There are two other html methods called
 * by makeWebSite(). One to print out the log info: makeWebLogPage(). The other to collect and print out information
 * obtained by different modules (units): makeWebInfoPage(). A user is responsible to implement/change implementation of the
 * latter method only if extra info web page is required. Finally, once a new AnalyzerUnit (derived from pure
 * virtual AnalyzerUnit class) is created, user is responsible for updating the Manager constructor and creating
 * an instance of the module (unit) in there.
 */
class AnalysisManager {

 public:

  //! Constructor - create instances of all available analyzer units & prepare web container
  //! @param[in] detector - List of all active sub-trackers + its passive parts and the beam pipe
  AnalysisManager(const Detector& detector);

  //! Destructor
  ~AnalysisManager();

  //! Initialize required analyzer unit
  //! @param[in] nTracks Number of tracks to be simulated (if needed by module)
  //! @param[in] analyzerName Name of the unit to be used
  //! @return True if there were no errors during processing, false otherwise
  bool initUnit(int nTracks, std::string analyzerName);

  //! Analyze required analyzer unit
  //! @param[in] analyzerName Name of the unit to be used
  //! @return True if there were no errors during processing, false otherwise
  bool analyzeUnit(std::string analyzerName);

  //! Visualize required analyzer unit
  //! @param[in] analyzerName Name of the unit to be used
  //! @return True if there were no errors during processing, false otherwise
  bool visualizeUnit(std::string analyzerName);

  //! Make web site - publish all results in a html format using the results of visualize method of all used units
  //! @param[in] addInfoPage Print out the info page?
  //! @param[in] addLogPage  Print out the log page?
  //! @return True if there were no errors during processing, false otherwise
  bool makeWebSite(bool addInfoPage, bool addLogPage);

 private:

  //! Prepare web site (html container for all results)
  //! @param[in] layoutName Layout name
  //! @param[in] webDir     Directory in which the web content will be saved
  //! @return True if there were no errors during processing, false otherwise
  bool prepareWebSite(std::string layoutName, std::string webDir);

  //! Prepare extra info web page
  //! @return True if here were no errors during processing, false otherwise
  bool makeWebInfoPage();

  //! Prepare web log page for information
  //! @return True if any logs found
  bool makeWebLogPage();

  std::unique_ptr<RootWSite> m_webSite;         //!< Web container, where all analysis results will be available
  bool                       m_webSitePrepared; //!< Web container correctly prepared

  std::map<std::string, std::unique_ptr<AnalyzerUnit>> m_units; //!< List of all available analyzer units

}; // Class

#endif /* INCLUDE_ANALYSISMANAGER_H_ */
