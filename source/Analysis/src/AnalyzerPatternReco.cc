/*
 * AnalyzerPatternReco.cc
 *
 *  Created on: 28. 11. 2016
 *      Author: drasal
 */
#include "AnalyzerPatternReco.h"

#include "IrradiationMap.h"
#include "MainConfigHandler.h"
#include "MessageLogger.h"
#include "Palette.h"
#include "RootWSite.h"
#include "SimParms.h"

//
// AnalyzerPatternReco constructor
//
AnalyzerPatternReco::AnalyzerPatternReco(const Detector& detector) : AnalyzerUnit("AnalyzerPatternReco", detector),
 m_nTracks(0),
 m_etaMin(-1*SimParms::getInstance().getMaxEtaCoverage()),
 m_etaMax(+1*SimParms::getInstance().getMaxEtaCoverage()),
 c_nBins(SimParms::getInstance().getMaxEtaCoverage()/vis_eta_step)  // Default number of bins in histogram from eta=0  to max_eta_coverage)
{
  m_chargedMap = nullptr;

};

//
// AnalyzerPatternReco destructor
//
AnalyzerPatternReco::~AnalyzerPatternReco()
{
  if (m_chargedMap!=nullptr) delete m_chargedMap;
}

//
// AnalyzerPatternReco init method
//
bool AnalyzerPatternReco::init(int nTracks)
{
  // Get occupancy map
  std::string directory = MainConfigHandler::getInstance().getIrradiationDirectory();
  bool chargedMapOK     = checkFile(SimParms::getInstance().chargedMapFile(), directory);

  if (chargedMapOK) m_chargedMap  = new IrradiationMap(directory + "/" + SimParms::getInstance().chargedMapFile());
  else {

       std::cout << "BUUUUUU" << std::endl;
  }

  // Set nTracks
  m_nTracks = nTracks;

  if (m_nTracks<=0 || m_trackers.size()==0) {

    std::ostringstream message;
    message << "AnalyzerPatternReco::init(): Number of simulation tracks zero or negative: " << m_nTracks << " or zero number of trackers to be initialized";
    logERROR(message.str());
    return false;
  }
  else {

    m_isInitOK = chargedMapOK;
    return m_isInitOK;
  }
}

//
// AnalyzerPatternReco analyze method
//
bool AnalyzerPatternReco::analyze()
{
  // Check that initialization OK
  if (!m_isInitOK) return false;

  m_isAnalysisOK = true;
  return m_isAnalysisOK;
}

//
// AnalyzerPatternReco visualization method
//
bool AnalyzerPatternReco::visualize(RootWSite& webSite)
{
  // Check that initialization & analysis OK
  if (!m_isInitOK && !m_isAnalysisOK) return false;

  // Set Rainbow palette for drawing
  Palette::setRootPalette();

  m_isVisOK = true;
  return m_isVisOK;
}

//
// Check that a file can be opened
//
bool AnalyzerPatternReco::checkFile(const std::string& fileName, const std::string& filePath)
{
  fstream     file;
  std::string fullFileName(filePath+"/"+fileName);
  file.open(fullFileName);
  if (file.is_open()) {

    file.close();
    return true;
  }
  else {

    logERROR("AnalyzerPatternReco - failed opening file: " + fullFileName);
    return false;
  }
}
