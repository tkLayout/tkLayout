/*
 * AnalysiManager.cc
 *
 *  Created on: 27. 4. 2016
 *      Author: drasal
 */
#include <AnalysisManager.h>

// List of modules
#include <AnalyzerModule.h>
#include <AnalyzerGeometry.h>
#include <AnalyzerResolution.h>
#include <InactiveSurfaces.h>
#include <Support.h>
#include <Tracker.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TH2I.h>
#include <SimParms.h>
#include <StopWatch.h>
#include <SvnRevision.h>
#include <rootweb.hh>

//
// Constructor - create instances of all available analyzer modules & prepare web container
//
AnalysisManager::AnalysisManager(std::vector<const Tracker*> activeTrackers,
                                 std::vector<const insur::InactiveSurfaces*> pasiveTrackers,
                                 std::vector<const Support*> supports) :
 m_webSite(nullptr),
 m_webSitePrepared(false)
{
  // Create AnalyzerGeometry
  AnalyzerGeometry* module = new AnalyzerGeometry(activeTrackers);
  m_modules[module->getName()] = module;

  // Prepare Web site
  auto simParms = SimParms::getInstance();
  m_webSitePrepared = prepareWebSite(simParms->getLayoutName(), simParms->getWebDir());
}

//
// Destructor - clean out memory
//
AnalysisManager::~AnalysisManager()
{
  for (auto iModule : m_modules) {

    if (iModule.second==nullptr) delete m_modules[iModule.first];
  }
}

//
// Initialize required analyzer module
//
bool AnalysisManager::initModule(int nTracks, std::string analyzerName)
{
  if (m_modules.find(analyzerName)!=m_modules.end()) {

    logINFO("Initializing " + analyzerName + " module");
    return m_modules[analyzerName]->init(nTracks);
  }
  else {

    std::cerr << any2str("\nERROR: AnalysisManager::initModule -> Module ("+analyzerName+") failed, no such module.") << std::endl;
    logERROR("AnalysisManager::initModule: Module ("+analyzerName+") failed, no such module.");
    return false;
  }
}

//
// Analyze required analyzer module
//
bool AnalysisManager::analyzeModule(std::string analyzerName)
{
  if (m_modules.find(analyzerName)!=m_modules.end()) {

    logINFO("Analyzing tracker by " + analyzerName + " module");
    return m_modules[analyzerName]->analyze();
  }
  else {

    std::cerr << any2str("\nERROR: AnalysisManager::analyzeModule -> Module ("+analyzerName+") failed, no such module.") << std::endl;
    logERROR("AnalysisManager::analyzerModule: Module ("+analyzerName+") failed, no such module.");
    return false;
  }
}

//
// Visualize required analyzer module
//
bool AnalysisManager::visualizeModule(std::string analyzerName)
{
  if (m_modules.find(analyzerName)!=m_modules.end()) {

    logINFO("Visualizing output from " + analyzerName + " module");
    return m_modules[analyzerName]->visualize(*m_webSite);
  }
  else {

    std::cerr << any2str("\nERROR: AnalysisManager::visualizeModule -> Module ("+analyzerName+") failed, no such module.") << std::endl;
    logERROR("AnalysisManager::visualizeModule: Module ("+analyzerName+") failed, no such module.");
    return false;
  }
}

//
// Make web site - publish all results in a html format using the results of visualize method of all used modules
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

    std::cerr << any2str("\nERROR: AnalysisManager::prepareWebSite -> Web site directory not set!") << std::endl;
    logERROR("AnalysisManager::prepareWebSite -> Web site directory not set!");
    return false;
  }

  // Set layout title
  if (layoutName!="") m_webSite->setTitle(layoutName);
  else {

    std::cerr << any2str("\nERROR: AnalysisManager::prepareWebSite -> Layout name not set!") << std::endl;
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
  m_webSite->setRevision(SvnRevision::revisionNumber);
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
  if (m_modules.find("AnalyzerGeometry")!=m_modules.end()) {

    AnalyzerGeometry* module = dynamic_cast<AnalyzerGeometry*>(m_modules["AnalyzerGeometry"]);
    myInfo = new RootWInfo("Number of tracks - geometry studies: ");
    myInfo->setValue(module->getNGeomTracks());
    myContent->addItem(myInfo);
  }
  if (m_modules.find("AnalyzerResolution")!=m_modules.end()) {

    AnalyzerResolution* module = dynamic_cast<AnalyzerResolution*>(m_modules["AnalyzerResolution"]);
    myInfo = new RootWInfo("Number of tracks - resolution studies: ");
    myInfo->setValue(module->getNSimTracks());
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

  for (auto iLevel=0; iLevel < MessageLogger::NumberOfLevels; ++iLevel) {

    if (!MessageLogger::hasEmptyLog(iLevel)) {

      bool defaultOpen = false;
      if (iLevel<=MessageLogger::WARNING) defaultOpen=true;

      anyLogFound = true;
      RootWContent& newContent = myPage.addContent(MessageLogger::getLevelName(iLevel), defaultOpen);
      newContent.addText("<pre>"+MessageLogger::getLatestLog(iLevel)+"</pre>");
    }
  }

  return anyLogFound;
}
