/*
 * AnalyzerUnit.cpp
 *
 *  Created on: 27. 11. 2015
 *      Author: drasal
 */

#include "AnalyzerUnit.h"

#include "BeamPipe.h"
#include "Detector.h"
#include "RootWSite.h"
#include <TCanvas.h>
#include "Tracker.h"


//
// Constructor - set active trackers only
//
AnalyzerUnit::AnalyzerUnit(std::string name, const Detector& detector) :
 m_trackers(detector.getTrackers()),
 m_beamPipe(detector.getBeamPipe())
{
  // Initialization by default false
  m_isInitOK = false;

  // Analysis by default not done
  m_isAnalysisOK = false;

  // Visualizaton by default not done
  m_isVisOK = false;

  // Set unique name
  m_name = name;
}

//
// Destructor
//
AnalyzerUnit::~AnalyzerUnit() {}

//
// Set canvas standard properties
//
void AnalyzerUnit::setCanvasProperties(TCanvas& canvas)
{
  canvas.SetFillColor(kWhite);
  canvas.SetBorderMode(0);
  canvas.SetBorderSize(0);
}
