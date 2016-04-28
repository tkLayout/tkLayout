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
#include <InactiveSurfaces.h>
#include <Support.h>
#include <Tracker.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TH2I.h>
#include <StopWatch.h>
#include <SvnRevision.h>
#include <rootweb.hh>

//
// Constructor - create instances of all available analyzer modules & prepare web container
//
AnalysisManager::AnalysisManager(std::string layoutName, std::string webDir,
                                 std::vector<const Tracker*> activeTrackers,
                                 std::vector<const insur::InactiveSurfaces*> pasiveTrackers,
                                 std::vector<const Support*> supports) :
 m_webSite(nullptr),
 m_webSitePrepared(false),
 m_webSiteDir(""),
 m_webLayoutName(""),
 m_commandLine("")
{
  // Create AnalyzerGeometry
  AnalyzerModule* module = (AnalyzerModule*)(new AnalyzerGeometry(activeTrackers));
  m_modules[module->getName()] = module;

  // Prepare Web site
  m_webSitePrepared = prepareWebSite(layoutName, webDir);
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
// Set command line options passed over to program to analyze data
//
void AnalysisManager::setCommandLine(int argc, char* argv[])
{
  m_commandLine = argv[1];
  for (auto i=2; i<argc; i++) m_commandLine += std::string(" ") + argv[i];
}

//
// Make web site - publish all results in a html format using the results of visualize method of all used modules
//
bool AnalysisManager::makeWebSite(bool addInfoPage, bool addLogPage) {

  if (!m_webSitePrepared) return false;
  else {

    startTaskClock(any2str("Creating website in: "+ m_webSiteDir));

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
  if (webDir!="") {

    m_webSiteDir = webDir;
    m_webSite->setTargetDirectory(m_webSiteDir);
  }
  else {

    std::cerr << any2str("\nERROR: AnalysisManager::prepareWebSite -> Web site directory not set!") << std::endl;
    logERROR("AnalysisManager::prepareWebSite -> Web site directory not set!");
    return false;
  }

  // Set layout title
  if (layoutName!="") {

    m_webLayoutName = layoutName;
    m_webSite->setTitle(m_webLayoutName);
  }
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
  RootWPage& myPage = m_webSite->addPage("Info");

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
