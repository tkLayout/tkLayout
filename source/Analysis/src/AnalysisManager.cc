/*
 * AnalysiManager.cc
 *
 *  Created on: 27. 4. 2016
 *      Author: drasal
 */
#include <AnalysisManager.h>

// List of units
#include <AnalyzerUnit.h>
#include <AnalyzerGeometry.h>
#include <AnalyzerMatBudget.h>
#include <AnalyzerResolution.h>

// Other include files
#include <BeamPipe.h>
#include <InactiveSurfaces.h>
#include <Tracker.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TH2I.h>
#include <SimParms.h>
#include <StopWatch.h>
#include <GitRevision.h>
#include <rootweb.h>

//
// Constructor - create instances of all available analyzer units & prepare web container
//
AnalysisManager::AnalysisManager(std::vector<const Tracker*> activeTrackers,
                                 std::vector<const insur::InactiveSurfaces*> passiveTrackers,
                                 const BeamPipe* beamPipe) :
 m_webSite(nullptr),
 m_webSitePrepared(false)
{
  AnalyzerUnit* unit = nullptr;

  // Create AnalyzerGeometry
  unit = new AnalyzerGeometry(activeTrackers, beamPipe);
  m_units[unit->getName()] = unit;

  // Create AnalyzerMatBudget
  unit = new AnalyzerMatBudget(activeTrackers, beamPipe);
  m_units[unit->getName()] = unit;

  // Prepare Web site
  auto simParms = SimParms::getInstance();
  m_webSitePrepared = prepareWebSite(simParms->getLayoutName(), simParms->getWebDir());
}

//
// Destructor - clean out memory
//
AnalysisManager::~AnalysisManager()
{
  for (auto iModule : m_units) {

    if (iModule.second==nullptr) delete m_units[iModule.first];
  }
}

//
// Initialize required analyzer unit
//
bool AnalysisManager::initUnit(int nTracks, std::string analyzerName)
{
  if (m_units.find(analyzerName)!=m_units.end()) {

    logINFO("Initializing " + analyzerName + " unit");
    return m_units[analyzerName]->init(nTracks);
  }
  else {

    logERROR("AnalysisManager::initModule: Module ("+analyzerName+") failed, no such unit.");
    return false;
  }
}

//
// Analyze required analyzer unit
//
bool AnalysisManager::analyzeUnit(std::string analyzerName)
{
  if (m_units.find(analyzerName)!=m_units.end()) {

    logINFO("Analyzing tracker by " + analyzerName + " unit");
    return m_units[analyzerName]->analyze();
  }
  else {

    logERROR("AnalysisManager::analyzerModule: Module ("+analyzerName+") failed, no such unit.");
    return false;
  }
}

//
// Visualize required analyzer unit
//
bool AnalysisManager::visualizeUnit(std::string analyzerName)
{
  if (m_units.find(analyzerName)!=m_units.end()) {

    logINFO("Visualizing output from " + analyzerName + " unit");
    return m_units[analyzerName]->visualize(*m_webSite);
  }
  else {

    logERROR("AnalysisManager::visualizeModule: Module ("+analyzerName+") failed, no such unit.");
    return false;
  }
}

//
// Make web site - publish all results in a html format using the results of visualize method of all used units
//
bool AnalysisManager::makeWebSite(bool addInfoPage, bool addLogPage) {

  if (!m_webSitePrepared) return false;
  else {

    startTaskClock(any2str("Creating website in: "+ SimParms::getInstance()->getWebDir()) );

    // Add info webPage
    if (addInfoPage) makeWebInfoPage();

    // Add log webPage
    if (addLogPage) makeWebLogPage();

    bool webOK = m_webSite->makeSite(false);
    stopTaskClock();

    return webOK;
  }
}

//
// Prepare web site (html container for all results)
//
bool AnalysisManager::prepareWebSite(std::string layoutName, std::string webDir)
{
  m_webSite = new RootWSite("TkLayout results");

  // Set directory
  if (webDir!="") m_webSite->setTargetDirectory(webDir);
  else {

    logERROR("AnalysisManager::prepareWebSite -> Web site directory not set!");
    return false;
  }

  // Set layout title
  if (layoutName!="") m_webSite->setTitle(layoutName);
  else {

    logERROR("AnalysisManager: Layout name not set!");
    return false;
  }

  // Set authors
  m_webSite->setComment("Layouts");
  m_webSite->setCommentLink("../");
  m_webSite->addAuthor("Giovanni Bianchi");
  m_webSite->addAuthor("Nicoletta De Maio");
  m_webSite->addAuthor("Stefano Martina");
  m_webSite->addAuthor("Stefano Mersi");
  m_webSite->setRevision(GitRevision::revisionNumber);
  return true;
}

//
// Prepare extra info web page
//
bool AnalysisManager::makeWebInfoPage()
{
  // Create web page: Log info
  RootWPage& myPage       = m_webSite->addPage("Info");
  RootWContent* myContent = nullptr;

  // Summary of parameters
  myContent = new RootWContent("Summary of tkLayout parameters");
  myPage.addContent(myContent);

  RootWInfo* myInfo;
  myInfo = new RootWInfo("Command line arguments");
  myInfo->setValue(SimParms::getInstance()->getCommandLine());
  myContent->addItem(myInfo);

  // Summary of used geometry & simulation tracks
  if (m_units.find("AnalyzerGeometry")!=m_units.end()) {

    AnalyzerGeometry* unit = dynamic_cast<AnalyzerGeometry*>(m_units["AnalyzerGeometry"]);
    myInfo = new RootWInfo("Number of tracks - geometry studies: ");
    myInfo->setValue(unit->getNGeomTracks());
    myContent->addItem(myInfo);
  }
  if (m_units.find("AnalyzerResolution")!=m_units.end()) {

    AnalyzerResolution* unit = dynamic_cast<AnalyzerResolution*>(m_units["AnalyzerResolution"]);
    myInfo = new RootWInfo("Number of tracks - resolution studies: ");
    myInfo->setValue(unit->getNSimTracks());
    myContent->addItem(myInfo);
  }

  // Summary of geometry config files
  auto simParms   = SimParms::getInstance();
  auto includeSet = simParms->getListOfConfFiles();
  if (!includeSet.empty()) {
    std::vector<std::string> includeNameSet;
    std::transform(includeSet.begin(), includeSet.end(), std::back_inserter(includeNameSet), [](const std::string& s) {

      auto pos = s.find_last_of('/');
      return (pos != string::npos ? s.substr(pos+1) : s);
    });
    RootWBinaryFileList* myBinaryFileList = new RootWBinaryFileList(includeNameSet.begin(), includeNameSet.end(),
                                                                    "Geometry configuration file(s)",
                                                                    includeSet.begin(), includeSet.end());
    myContent->addItem(myBinaryFileList);
  }

  return true;
}

//
// Prepare web log page for information
//
bool AnalysisManager::makeWebLogPage()
{
  // Found any log
  bool anyLogFound=false;

  // Create web page: Log info
  RootWPage& myPage = m_webSite->addPage("Log");

  // Check logs
  if (!MessageLogger::hasEmptyLog(MessageLogger::ERROR))        myPage.setAlert(1);
  else if (!MessageLogger::hasEmptyLog(MessageLogger::WARNING)) myPage.setAlert(0.5);

  // Print summary logs for each verbosity level
  for (auto iLevel=0; iLevel < MessageLogger::NumberOfLevels; ++iLevel) {

    if (!MessageLogger::hasEmptyLog(iLevel)) {

      bool defaultOpen = false;
      if (iLevel<=MessageLogger::WARNING) defaultOpen=true;

      anyLogFound = true;
      RootWContent& newContent = myPage.addContent(MessageLogger::getLevelName(iLevel), defaultOpen);
      newContent.addText("<pre>"+MessageLogger::getSummaryLog(iLevel)+"</pre>");
    }
  }

  return anyLogFound;
}
