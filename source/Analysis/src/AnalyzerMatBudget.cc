/*
 * AnalyzerMatBudget.cc
 *
 *  Created on: 9. 6. 2016
 *      Author: Z. Drasal (CERN)
 */
#include "AnalyzerMatBudget.h"

// Include files
#include <BeamPipe.h>
#include <global_constants.h>
#include <rootweb.h>
#include <ostream>
#include <Tracker.h>
#include <SimParms.h>

//
// AnalyzerMatBudget constructor
//
AnalyzerMatBudget::AnalyzerMatBudget(std::vector<const Tracker*> trackers, const BeamPipe* beamPipe) : AnalyzerUnit("AnalyzerMatBudget", trackers, beamPipe),
 m_nTracks(0),
 m_etaSpan(geom_max_eta_coverage - geom_max_eta_coverage),
 m_etaMin(-1*geom_max_eta_coverage),
 m_etaMax(+1*geom_max_eta_coverage)
{};

//
// AnalyzerMatBudget init method
//
bool AnalyzerMatBudget::init(int nMatTracks)
{
  // Set nTracks
  m_nTracks = nMatTracks;

  if (m_nTracks<=0 || m_trackers.size()==0) {

    std::ostringstream message;
    message << "AnalyzerMatBudget::init(): Number of material tracks zero or negative: " << m_nTracks << " or zero number of trackers to be initialized";
    logERROR(message.str());
    return false;
  }
  else {

    for (auto iTracker : m_trackers) {

      std::string trkName = iTracker->myid();

    } // For trackers

    m_isInitOK = true;
    return m_isInitOK;
  }
}

//
// AnalyzerMatBudget analysis method
//
bool AnalyzerMatBudget::analyze()
{
  // Check that initialization OK
  if (!m_isInitOK) return false;

  // Print
  std::cout << std::endl;

  // Go through all trackers
  for (auto iTracker : m_trackers) {

    std::string trkName =  iTracker->myid();

    // Print
    std::cout << " " << iTracker->myid() << " tracker"<< std::endl;

    // Analyze geometry layout
    for (auto iTrack=0; iTrack<m_nTracks; iTrack++) {

    }

  } // For trackers

  m_isAnalysisOK = true;
  return m_isAnalysisOK;
}

//
// AnalyzerMatBudget visualization method
//
bool AnalyzerMatBudget::visualize(RootWSite& webSite)
{
  // Check that initialization & analysis OK
  if (!m_isInitOK || !m_isAnalysisOK) return false;

  // Go through all trackers & prepare web content
  int webPriority         = 89;
  RootWPage*    myPage    = nullptr;
  RootWContent* myContent = nullptr;

  for (auto iTracker : m_trackers) {

    // Tracker name
    std::string trkName =  iTracker->myid();

    // Create dedicated web-page & its content
    std::string        pageTitle  = "MatBudget";
    if (trkName != "") pageTitle +=" (" + trkName + ")";
    myPage = new RootWPage(pageTitle);

    // Page address
    std::string pageAddress = "index"+trkName+".html";
    myPage->setAddress(pageAddress);

    webSite.addPage(myPage, webPriority);
    webPriority--;

  } // Trackers

  return true;
}
