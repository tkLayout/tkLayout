/*
 * AnalyzerModule.cpp
 *
 *  Created on: 27. 11. 2015
 *      Author: drasal
 */

#include "AnalyzerModule.h"

AnalyzerModule::AnalyzerModule(std::string name, std::vector<const Tracker*> trackers)
{
  // Initialization by default false
  m_isInitOK = false;

  // Analysis by default not done
  m_isAnalysisOK = false;

  // Set unique name
  m_name = name;

  // Set geometry, i.e. individual trackers
  for (auto it : trackers) m_trackers.push_back(it);
}

AnalyzerModule::~AnalyzerModule() {}

