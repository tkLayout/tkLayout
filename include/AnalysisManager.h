/*
 * AnalysisManager.h
 *
 *  Created on: 27. 4. 2016
 *      Author: drasal
 */

#ifndef INCLUDE_ANALYSISMANAGER_H_
#define INCLUDE_ANALYSISMANAGER_H_

#include <map>
#include <vector>
#include <string>

class AnalyzerModule;
class RootWSite;
namespace insur {
  class InactiveSurfaces;
}
class Support;
class Tracker;

/*
 * The core analysis class, containing all analyzers to be used. Each analyzer is identified by its unique name
 * and is called through initModule, analyzeModule & visualizeModule methods. Init method initializes parameters,
 * histograms etc. of given module, analyze method performs the analysis and visualization outputs the results
 * to html format. All html web pages created by different modules are collected in html format in a common web
 * container (main html page) and printed-out using makeWebSite() method. There are two other html methods called
 * by makeWebSite(). One to print out the log info: makeWebLogPage(). The other to collect and print out information
 * obtained by different modules: makeWebInfoPage(). A user is responsible to implement/change implementation of the
 * latter method only if extra infor web page is required. Finally, once a new AnalyzerModule is created, user is
 * responsible for updating the Manager constructor and creating an instance of the module in there.
 */
class AnalysisManager {

 public:

  //! Constructor - create instances of all available analyzer modules & prepare web container
  //! @param[in] layoutName     Layout name
  //! @param[in] webDir         Directory in which the web content will be saved
  //! @param[in] activeTrackers List of all active sub-trackers
  //! @param[in] pasiveTrackers List of pasives related to active sub-trackers
  //! @param[in] supports       List of independent support structures not directly related to active sub-trackers
  AnalysisManager(std::string layoutName, std::string webDir,
                  std::vector<const Tracker*> activeTrackers,
                  std::vector<const insur::InactiveSurfaces*> pasiveTrackers,
                  std::vector<const Support*> supports);

  //! Destructor -
  ~AnalysisManager();

  //! Initialize required analyzer module
  //! @param[in] nTracks Number of tracks to be simulated (if needed by module)
  //! @param[in] analyzerName Name of the module to be used
  //! @return True if there were no errors during processing, false otherwise
  bool initModule(int nTracks, std::string analyzerName);

  //! Analyze required analyzer module
  //! @param[in] analyzerName Name of the module to be used
  //! @return True if there were no errors during processing, false otherwise
  bool analyzeModule(std::string analyzerName);

  //! Visualize required analyzer module
  //! @param[in] analyzerName Name of the module to be used
  //! @return True if there were no errors during processing, false otherwise
  bool visualizeModule(std::string analyzerName);

  //! Make web site - publish all results in a html format using the results of visualize method of all used modules
  //! @param[in] addInfoPage Print out the info page?
  //! @param[in] addLogPage  Print out the log page?
  //! @return True if there were no errors during processing, false otherwise
  bool makeWebSite(bool addInfoPage, bool addLogPage);

  //! Set command line options passed over to program to analyze data
  //! @param[in] commandLine    Content of command line
  void setCommandLine(int argc, char* argv[]);

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

  RootWSite*  m_webSite;         //!< Web container, where all analysis results will be available
  bool        m_webSitePrepared; //!< Web container correctly prepared
  std::string m_webSiteDir;      //!< Web container directory
  std::string m_webLayoutName;   //!< Web layout name

  std::string m_commandLine;     //!< Command line options passed over to program to analyze data

  std::map<std::string, AnalyzerModule*> m_modules; //!< List of all available analyzer modules

}; // Class

#endif /* INCLUDE_ANALYSISMANAGER_H_ */
