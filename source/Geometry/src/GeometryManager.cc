/*
 * GeometryManager.cc
 *
 *  Created on: 27. 4. 2016
 *      Author: drasal
 */
#include <GeometryManager.h>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <fstream>
#include <InactiveSurfaces.h>
#include <Materialway.h>
#include <MainConfigHandler.h>
#include <MessageLogger.h>
#include <SimParms.h>
#include <StopWatch.h>
#include <Support.h>
#include <Tracker.h>
#include <WeightDistributionGrid.h> // TODO: Dummy for now

//
// Constructor - initializes manager -> reads-in all geometry configuration files (using the mainConfigHandler)
//
GeometryManager::GeometryManager(std::string baseGeomFile) :
 m_initOK(false),
 m_geomTree(nullptr)
{
  // Set geometry file, directory, layout name ...
  setGeometryFile(baseGeomFile);

  // Geometry/config file not defined
  std::ifstream inputFile(baseGeomFile);
  if (inputFile.fail()) {

    std::cerr << any2str("n\ERROR: GeometryManager::GeometryManager -> Cannot open geometry file: "+baseGeomFile) << std::endl;
    logERROR(any2str("GeometryManager::GeometryManager: Cannot open geometry file: "+baseGeomFile));
  }
  else {

    m_initOK = true;
    logINFO("Reading base geometry file: "+baseGeomFile);

    //
    // Parse config file & build property tree
    std::stringstream outputConfig;
    m_includeSet = MainConfigHandler::getInstance().preprocessConfiguration(inputFile, outputConfig, baseGeomFile);

    m_geomTree = new boost::property_tree::ptree();
    boost::property_tree::info_parser::read_info(outputConfig, *m_geomTree);

    //
    // Fill SimParms singleton class (container for generic info necessary for simulation), store data to property tree
    auto simParms = SimParms::getInstance();

    simParms->store(getChild(*m_geomTree, "SimParms"));
    simParms->crosscheck();
  }
}

//
// Destructor - clean out the memory
//
GeometryManager::~GeometryManager()
{
  // Clear memory
  if (m_geomTree==nullptr) delete m_geomTree;

  // Active sub-trackers
  for (auto iTracker : m_activeTrackers) {

    if (iTracker==nullptr) delete iTracker;
  }

  // Standalone support
  for (auto iSupport : m_supports) {

    if (iSupport==nullptr) delete iSupport;
  }

  // Pasive components related to active sub-trackers
  for (auto iPasive : m_pasiveTrackers) {

    if (iPasive==nullptr) delete iPasive;
  }
}

//
// Build active part of the tracker (a bare-bones geometry of active modules). The resulting tracker object consists of individual sub-trackers
// (pixel, strip, central, forward, ...). This procedure replaces the previously registered tracker (applying correct memory managment), if such an object
// existed.
// Return True if there were no errors during processing, false otherwise
//
bool GeometryManager::buildActiveTracker()
{
  // Check that GeometryManager properly initialized
  if (!m_initOK) {

    std::cerr << any2str("\nERROR: GeometryManager::buildActiveTracker() -> Can't built active tracker if GeometryManager not properly initialized in constructor");
    logERROR(any2str("GeometryManager::buildActiveTracker(): Can't built active tracker if GeometryManager not properly initialized in constructor"));
    return false;
  }

  // Clear out content of trackers container
  for (auto iTracker : m_activeTrackers) {

    if (iTracker!=nullptr) delete iTracker;
  }
  m_activeTrackers.clear();

  // Build active tracker -> look for tag "Tracker"
  startTaskClock("Building active trackers");
  try {

    auto childRange = getChildRange(*m_geomTree, "Tracker");
    std::for_each(childRange.first, childRange.second, [&](const ptree::value_type& kv) {

      Tracker* trk = new Tracker(kv.second);
      trk->build();
      m_activeTrackers.push_back(trk);

    }); // Trackers

    // Declare unmatched properties
    std::set<string> unmatchedProperties = PropertyObject::reportUnmatchedProperties();
    if (!unmatchedProperties.empty()) {

      std::ostringstream message;
      message << "GeometryManager::buildActiveTracker(): The following unknown properties were ignored:" << std::endl;

      for (const string& property  : unmatchedProperties) {
        message << "  " << property << std::endl;
      }
      std::cerr << any2str("\nERROR: "+message.str()) << std::endl;
      logERROR(message);
    }
  } // Try
  catch (PathfulException& e) {

    std::cerr << any2str("\nERROR: " +e.path()+ " : " +e.what()) << std::endl;
    logERROR(e.path()+ " : " +e.what());
    stopTaskClock();
    return false;
  }

  stopTaskClock();
  return true;
}

//
// Build tracker support structures, which are independent on individual sub-trackers (the supports directly related to sub-trackers are built as pasive
// components of active sub-tracker). This procedure replaces the previously registered support (applying correct memory managment), if such an object
// existed.
//
bool GeometryManager::buildTrackerSupport()
{
  // Check that GeometryManager properly initialized
  if (!m_initOK) {

    std::cerr <<
    logERROR(any2str("\nERROR: GeometryManager::buildTrackerSupport() -> Can't built tracker supports if GeometryManager not properly initialized in constructor"));
    return false;
  }

  // Clear out content of supports container
  for (auto iSupport : m_supports) {

    if (iSupport!=nullptr) delete iSupport;
  }
  m_supports.clear();

  // Build tracker support -> look for tag "Support"
  startTaskClock("Building tracker support");
  try {

    // Look for tag "Support" not associated with a concrete Tracker and build supports
    auto childRange = getChildRange(*m_geomTree, "Support");
    std::for_each(childRange.first, childRange.second, [&](const ptree::value_type& kv) {

      Support* support = new Support();
      support->myid(kv.second.get_value(0));
      support->store(kv.second);
      support->build();
      m_supports.push_back(support);

    }); // Supports

  }// Try
  catch (PathfulException& e) {

    std::cerr << any2str("\nERROR: " +e.path()+ " : " +e.what()) << std::endl;
    logERROR(e.path()+ " : " +e.what());
    stopTaskClock();
    return false;
  }
  stopTaskClock();
  return true;
}

//
// Build all pasive components related to individual active sub-trackers. This procedure replaces the previously registered support (applying correct
// memory managment), if such an object existed.
//
bool GeometryManager::buildPasiveTracker()
{
  if (m_activeTrackers.size()>0) {

    startTaskClock("Building pasive materials of active trackers");

    for (auto iTracker : m_activeTrackers) {

      // Build new inactive surface
      insur::InactiveSurfaces*         iPasive = new insur::InactiveSurfaces();
      material::WeightDistributionGrid wdGrid(0.1); // TODO: Dummy for now


      // Get all materials in the way for given tracker
      material::Materialway materialway;
      materialway.build(*iTracker, *iPasive, wdGrid);
    }

    stopTaskClock();
    return true;
  }
  else return false;
}

//
// Get active sub-trackers
//
std::vector<const Tracker*> GeometryManager::getActiveTrackers() const
{
  std::vector<const Tracker*> activeTrackers;

  for (auto iTracker : m_activeTrackers) {

    const Tracker* trk = iTracker;
    activeTrackers.push_back(trk);
  }

  return activeTrackers;
}

//
// Get tracker supports, which are independent on individual sub-trackers
//
std::vector<const Support*> GeometryManager::getTrackerSupports() const
{
  std::vector<const Support*> trackerSupports;

  for (auto iSupport : m_supports) {

    const Support* support = iSupport;
    trackerSupports.push_back(support);
  }

  return trackerSupports;
}

//
// Get pasive components related to active sub-trackers
//
std::vector<const insur::InactiveSurfaces*> GeometryManager::getPasiveTrackers() const
{
  std::vector<const insur::InactiveSurfaces*> pasiveTrackers;

  for (auto iPasive : m_pasiveTrackers) {

    const insur::InactiveSurfaces* pasive = iPasive;
    pasiveTrackers.push_back(pasive);
  }

  return pasiveTrackers;
}

//
// Set base geometry file, directory, layout name, ...
//
void GeometryManager::setGeometryFile(std::string baseGeomFile)
{

  // Set file name according to geometry config file
  std::vector<std::string> info;
  boost::algorithm::split(info, baseGeomFile, boost::algorithm::is_any_of("/"));
  if (info.size()!=0) m_geomFile = info[info.size()-1];

  // Config directory
  m_baseDir = baseGeomFile;
  size_t pos = m_baseDir.find_last_of('/');
  if (pos != string::npos) { m_baseDir.erase(pos); }

  // Html directory
  m_htmlDir = m_baseDir + "/" + c_defaultHtmlDir;

  // Layout name
  m_layoutName = m_geomFile;
  pos = m_layoutName.find_last_of('.');
  if (pos != string::npos) { m_layoutName.erase(pos); }

  logINFO("Analysing layout: "         + m_layoutName);
  logINFO("Results will be saved in: " + m_htmlDir);
}
