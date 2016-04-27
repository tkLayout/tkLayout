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
 * The core analysis class, containing all analyzers to be used. Each analyzer is identified by its unique name.
 * Call analyzer init, analyze or visuzualize methods ...
 * The result are collected in html format in a web container and printed-out using makeWebSite method.
 */
class AnalysisManager {

 public:

  //! Constructor - create instances of all available analyzer modules & prepare web container
  //! @param[in] layoutName     Layout name
  //! @param[in] webDir         Directory in which the web content will be saved
  //! @param[in] activeTrackers List of all active sub-trackers
  //! @param[in] pasiveTrackers List of pasives related to active sub-trackers
  //! @param[in] supports       List of independent support structures not directly related to active sub-trackers
  AnalysisManager(std::string layoutName, std::string webDir, std::vector<Tracker*> activeTrackers, std::vector<insur::InactiveSurfaces*> pasiveTrackers, std::vector<Support*> supports);

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
  //! @param[in] addLogPage Print out the log page?
  //! @return True if there were no errors during processing, false otherwise
  bool makeWebSite(bool addLogPage);

 private:

  //! Prepare web site (html container for all results)
  //! @param[in] layoutName Layout name
  //! @param[in] webDir     Directory in which the web content will be saved
  //! @return True if there were no errors during processing, false otherwise
  bool prepareWebSite(std::string layoutName, std::string webDir);

  //! Prepare web log page for information
  //! @return True if any logs found
  bool makeWebLogPage();

  RootWSite*  m_webSite;         //!< Web container, where all analysis results will be available
  bool        m_webSitePrepared; //!< Web container correctly prepared
  std::string m_webSiteDir;      //!< Web container directory
  std::string m_webLayoutName;   //!< Web layout name

  std::map<std::string, AnalyzerModule*> m_modules; //!< List of all available analyzer modules

}; // Class

#endif /* INCLUDE_ANALYSISMANAGER_H_ */
