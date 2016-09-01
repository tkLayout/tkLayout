/*
 * ExtractorCMSSW.cc
 *
 *  Created on: 1. 9. 2016
 */
#include "ExtractorCMSSW.h"

//
// Constructor
//
ExtractorCMSSW::ExtractorCMSSW(const Detector& detector) :
 AnalyzerUnit("ExtractorCMSSW", detector)
{}

//
// Destructor
//
ExtractorCMSSW::~ExtractorCMSSW()
{}

//
// Initialize - mostly histograms & other containers
// @return True if OK
bool ExtractorCMSSW::init(int nGeomTracks)
{
  //TODO: Transfer tkLayout CMSSW extractor init parts here

  m_isInitOK = true;
  return m_isInitOK;
}

//
// Gather information about the geometry layout (if init OK) & output the data to XML file
// @return True if OK
bool ExtractorCMSSW::analyze()
{
  //TODO: Transfer tkLayout CMSSW extractor analysis parts here

  m_isAnalysisOK = true;
  return m_isAnalysisOK;
}
