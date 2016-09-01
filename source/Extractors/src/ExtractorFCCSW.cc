/*
 * ExtractorFCCSW.cc
 *
 *  Created on: 1. 9. 2016
 *      Author: Z.Drasal (CERN)
 */
#include "ExtractorFCCSW.h"

//
// Constructor
//
ExtractorFCCSW::ExtractorFCCSW(const Detector& detector) :
 AnalyzerUnit("ExtractorFCCSW", detector)
{}

//
// Destructor
//
ExtractorFCCSW::~ExtractorFCCSW()
{}

//
// Initialize - mostly histograms & other containers
// @return True if OK
bool ExtractorFCCSW::init(int nGeomTracks)
{
  m_isInitOK = true;
  return m_isInitOK;
}

//
// Gather information about the geometry layout (if init OK) & output the data to XML file
// @return True if OK
bool ExtractorFCCSW::analyze()
{
  m_isAnalysisOK = true;
  return m_isAnalysisOK;
}
