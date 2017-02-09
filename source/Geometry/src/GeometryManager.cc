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
#include "Detector.h"
#include <fstream>
#include <MainConfigHandler.h>
#include <MessageLogger.h>
#include <SimParms.h>
#include <StopWatch.h>

//
// Constructor - initializes manager -> reads-in all geometry configuration files (using the mainConfigHandler)
//
GeometryManager::GeometryManager(std::string baseGeomFile) :
 m_initOK(false)
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

    m_geomTree = std::unique_ptr<boost::property_tree::ptree>(new boost::property_tree::ptree());
    boost::property_tree::info_parser::read_info(outputConfig, *m_geomTree);

    //
    // Fill SimParms singleton class (container for generic info necessary for simulation), store data to property tree
    auto& simParms = SimParms::getInstance();

    simParms.store(getChild(*m_geomTree, "SimParms"));
    simParms.crosscheck();

    //
    // Create detector
    std::string detName = MainConfigHandler::getInstance().getProjectName();
    m_detector          = std::unique_ptr<Detector>(new Detector(detName));
  }
}

//
// Destructor - clean out the memory
//
GeometryManager::~GeometryManager()
{
  m_detector.reset();
  m_geomTree.reset();
}

//
// Build active part of the detector: individual trackers: pixel, strip, central, forward, ... (a bare-bones geometry of active modules).
// Return True if there were no errors during processing, false otherwise
//
bool GeometryManager::buildActiveDetector()
{
  // Check that GeometryManager properly initialized
  if (!m_initOK) {

    logERROR(any2str("GeometryManager::buildActiveDetector: Can't built active parts of detector if GeometryManager not properly initialized in constructor"));
    return false;
  }

  return m_detector->buildActiveTracker(*m_geomTree);
}

//
// Build all passive components related to individual active sub-trackers & beam-pipe.
//
bool GeometryManager::buildPassiveDetector()
{
  // Check that GeometryManager properly initialized
  if (!m_initOK) {

    logERROR(any2str("GeometryManager::buildPassiveDetector: Can't built passive parts of detector if GeometryManager not properly initialized in constructor"));
    return false;
  }

  return (m_detector->buildPassiveTracker() && m_detector->buildBeamPipe(*m_geomTree));
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
