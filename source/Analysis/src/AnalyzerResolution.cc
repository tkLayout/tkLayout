/*
 * AnalyzerResolution.cc
 *
 *  Created on: 20. 4. 2016
 *      Author: drasal
 */
#include <AnalyzerResolution.h>

#include <rootweb.hh>
#include <Tracker.h>


AnalyzerResolution::AnalyzerResolution(std::vector<const Tracker*> trackers) : AnalyzerUnit("AnalyzerResolution", trackers),
 m_nTracks(0)
{};

bool AnalyzerResolution::init(int nTracks)
{
  // Set nTracks
  m_nTracks = nTracks;

  if (m_nTracks<=0 || m_trackers.size()==0) {

    std::ostringstream message;
    message << "AnalyzerResolution::init(): Number of simulation tracks zero or negative: " << m_nTracks << " or zero number of trackers to be initialized";
    logERROR(message.str());
    return false;
  }
  else {

    m_isInitOK = true;
    return m_isInitOK;
  }
}

bool AnalyzerResolution::analyze()
{
  // Check that initialization OK
  if (!m_isInitOK) return false;

  return true;
}

bool AnalyzerResolution::visualize(RootWSite& webSite)
{
  // Check that initialization & analysis OK
  if (!m_isInitOK && !m_isAnalysisOK) return false;

  return true;
}


