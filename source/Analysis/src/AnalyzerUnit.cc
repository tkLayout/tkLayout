/*
 * AnalyzerUnit.cpp
 *
 *  Created on: 27. 11. 2015
 *      Author: drasal
 */

#include "AnalyzerUnit.h"
#include "TCanvas.h"

//
// Constructor - set active trackers only
//
AnalyzerUnit::AnalyzerUnit(std::string name, std::vector<const Tracker*> trackers)
{
  // Initialization by default false
  m_isInitOK = false;

  // Analysis by default not done
  m_isAnalysisOK = false;

  // Set unique name
  m_name = name;

  // Set geometry, i.e. individual trackers
  for (auto it : trackers) m_trackers.push_back(it);

  // Set geometry of beam pipe
  m_beamPipe = nullptr;
}

//
// Constructor - set active trackers & beam pipe to be analyzed
//
AnalyzerUnit::AnalyzerUnit(std::string name, std::vector<const Tracker*> trackers, const BeamPipe* beamPipe)
{
  // Initialization by default false
  m_isInitOK = false;

  // Analysis by default not done
  m_isAnalysisOK = false;

  // Set unique name
  m_name = name;

  // Set geometry, i.e. individual trackers
  for (auto it : trackers) m_trackers.push_back(it);

  // Set geometry of beam pipe
  m_beamPipe = beamPipe;
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
