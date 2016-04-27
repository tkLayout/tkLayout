/*
 * AnalyzerResolution.cc
 *
 *  Created on: 20. 4. 2016
 *      Author: drasal
 */
#include <AnalyzerResolution.h>

#include <rootweb.hh>
#include <Tracker.h>


AnalyzerResolution::AnalyzerResolution(std::vector<Tracker*> trackers) : AnalyzerModule("AnalyzerResolution", trackers) {};

bool AnalyzerResolution::analyze()
{
  return true;
}

bool AnalyzerResolution::visualize(RootWSite& webSite)
{
  return true;
}


