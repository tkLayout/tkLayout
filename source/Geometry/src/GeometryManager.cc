/*
 * GeometryManager.cc
 *
 *  Created on: 27. 4. 2016
 *      Author: drasal
 */
#include <GeometryManager.h>

#include <BeamPipe.h>
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
#include <Tracker.h>

//
// Constructor - initializes manager -> reads-in all geometry configuration files (using the mainConfigHandler)
//
GeometryManager::GeometryManager(std::string baseGeomFile) :
 m_initOK(false),
 m_geomTree(nullptr),
 m_beamPipe(nullptr)
{
  // Set geometry file, directory, layout name ...
  setGeometryFile(baseGeomFile);

  // Geometry/config file not defined
  std::ifstream inputFile(baseGeomFile);
  if (inputFile.fail()) {

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

  // Passive components related to active sub-trackers
  for (auto iPassive : m_passiveTrackers) {

    if (iPassive==nullptr) delete iPassive;
  }

  // Beam pipe
  if (m_beamPipe!=nullptr) delete m_beamPipe;
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
    for (auto it=childRange.first; it!=childRange.second; ++it) {

      // Message
      if (it==childRange.first) std::cout << std::endl;
      std::cout << " " << it->second.data() << " tracker" << std::endl;

      // Create tracker
      Tracker* trk = new Tracker(it->second);
      trk->build();
      trk->setup();

      m_activeTrackers.push_back(trk);

    } // Trackers

    // Declare unmatched properties
    std::set<string> unmatchedProperties = PropertyObject::reportUnmatchedProperties();
    if (!unmatchedProperties.empty()) {

      std::ostringstream message;
      message << "GeometryManager::buildActiveTracker(): The following unknown properties were ignored:" << std::endl;

      for (const string& property  : unmatchedProperties) {
        message << "  " << property << std::endl;
      }
      logERROR(message);
    }
  } // Try
  catch (PathfulException& e) {

    logERROR(e.path()+ " : " +e.what());
    stopTaskClock();
    return false;
  }

  stopTaskClock();
  return true;
}

//
// Build all passive components related to individual active sub-trackers. This procedure replaces the previously registered support (applying correct
// memory managment), if such an object existed.
//
bool GeometryManager::buildPassiveTracker()
{
  if (m_activeTrackers.size()>0) {

    startTaskClock("Building passive materials of active trackers");

    for (auto iTracker : m_activeTrackers) {

      // Build new inactive surface
      insur::InactiveSurfaces* iPassive = new insur::InactiveSurfaces();

      // Get all materials in the way for given tracker
      material::Materialway materialway;
      materialway.build(*iTracker, *iPassive);
    }

    stopTaskClock();
    return true;
  }
  else return false;
}

//
// Build beam pipe. This procedure replaces the previously registered beam pipe
//
bool GeometryManager::buildBeamPipe()
{
  // Check that GeometryManager properly initialized
  if (!m_initOK) {

    logERROR(any2str("GeometryManager::buildBeamPipe(): Can't built beam pipe if GeometryManager not properly initialized in constructor"));
    return false;
  }

  // Clear out memory if beam pipe built before
  if (m_beamPipe!=nullptr) {

    delete m_beamPipe;
    m_beamPipe = nullptr;
  }

  // Build beam pipe -> look for tag "BeamPipe"
  startTaskClock("Building beam pipe");
  try {

    auto childRange = getChildRange(*m_geomTree, "BeamPipe");
    for (auto it=childRange.first; it!=childRange.second; ++it) {

      // Declare if more than one beam pipes defined
      if (m_beamPipe!=nullptr) {

        std::ostringstream message;
        message << "GeometryManager::buildBeamPipe(): More than one beam pipe defined -> only first one found used!" << std::endl;
        logERROR(message);
      }
      else {

        std::cout << "Building new beam pipe" << std::endl;
        m_beamPipe = new BeamPipe(it->second);
        m_beamPipe->build();
        m_beamPipe->setup();
      }

    }; // Beam pipes

    // Declare unmatched properties
    std::set<string> unmatchedProperties = PropertyObject::reportUnmatchedProperties();
    if (!unmatchedProperties.empty()) {

      std::ostringstream message;
      message << "GeometryManager::buildBeamPipe(): The following unknown properties were ignored:" << std::endl;

      for (const string& property  : unmatchedProperties) {
        message << "  " << property << std::endl;
      }
      logERROR(message);
    }
  } // Try
  catch (PathfulException& e) {

    logERROR(e.path()+ " : " +e.what());
    stopTaskClock();
    return false;
  }
  stopTaskClock();
  return true;
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
// Get passive components related to active sub-trackers
//
std::vector<const insur::InactiveSurfaces*> GeometryManager::getPassiveTrackers() const
{
  std::vector<const insur::InactiveSurfaces*> passiveTrackers;

  for (auto iPassive : m_passiveTrackers) {

    const insur::InactiveSurfaces* passive = iPassive;
    passiveTrackers.push_back(passive);
  }

  return passiveTrackers;
}

//
// Get beam pipe
//
const BeamPipe* GeometryManager::getBeamPipe() const
{
  const BeamPipe* beamPipe = m_beamPipe;
  return beamPipe;
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
