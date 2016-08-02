/*
 * AnalyzerMatBudget.cc
 *
 *  Created on: 9. 6. 2016
 *      Author: Z. Drasal (CERN)
 */
#include "AnalyzerMatBudget.h"

// Include files
#include "BeamPipe.h"
#include "DetectorModule.h"
#include "Disk.h"
#include <global_constants.h>
#include <Hit.h>
#include <MaterialProperties.h>
#include <Math/Vector3D.h>
#include <ModuleCap.h>
#include <ostream>
#include <Palette.h>
#include "RootWContent.h"
#include "RootWImage.h"
#include "RootWPage.h"
#include "RootWSite.h"
#include "RootWTable.h"
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TPad.h>
#include <TProfile.h>
#include <Tracker.h>
#include <Track.h>
#include <TRandom3.h>
#include <Units.h>
#include <SimParms.h>

//
// AnalyzerMatBudget constructor
//
AnalyzerMatBudget::AnalyzerMatBudget(const Detector& detector) : AnalyzerUnit("AnalyzerMatBudget", detector),
 m_nTracks(0),
 m_etaSpan(2*SimParms::getInstance().getMaxEtaCoverage()),
 m_etaMin(-1*SimParms::getInstance().getMaxEtaCoverage()),
 m_etaMax(+1*SimParms::getInstance().getMaxEtaCoverage())
{};

//
// Destructor
//
AnalyzerMatBudget::~AnalyzerMatBudget()
{}

//
// AnalyzerMatBudget init method
//
bool AnalyzerMatBudget::init(int nMatTracks)
{
  // Set nTracks
  m_nTracks = nMatTracks;

  // Compute shooting intervals to analyze geometry
  double etaMin = -1*SimParms::getInstance().getMaxEtaCoverage();
  double etaMax = +1*SimParms::getInstance().getMaxEtaCoverage();

  float  safeMargin = c_etaSafetyMargin;
  m_etaSpan         = (etaMax - etaMin)*(1. + safeMargin);
  m_etaMax          = etaMax * (1 + safeMargin);

  if (m_nTracks<=0 || m_trackers.size()==0) {

    std::ostringstream message;
    message << "AnalyzerMatBudget::init(): Number of material tracks zero or negative: " << m_nTracks << " or zero number of trackers to be initialized";
    logERROR(message.str());
    return false;
  }
  else if (m_beamPipe==nullptr) {

    std::ostringstream message;
    message << "AnalyzerMatBudget::init(): No beam pipe defined for material budget studies";
    logERROR(message.str());
    return false;
  }
  else {

    // Calculate final plot range, increase by safety factor
    double maxZ = 0;
    double maxR = 0;
    for (auto iTrk : m_trackers) {

      maxZ = MAX(maxZ, iTrk->maxZ()*vis_safety_factor);
      maxR = MAX(maxR, iTrk->maxR()*vis_safety_factor);
    }

    // Prepare Csv containers -> keep final pT/p resolution in csv file
    m_csvMatBudget = std::unique_ptr<CsvTextBuilder>(new CsvTextBuilder());

    // Prepare histogram containers for each subdetector (components will be initialized on the fly)
    for (int iTrk=0; iTrk<=m_trackers.size(); iTrk++) {

      // Tracker name
      std::string trkName;
      if (iTrk==m_trackers.size()) trkName = "Tracker";
      else                         trkName = m_trackers[iTrk]->myid();

      // Maps
      double step = 5.*Units::mm;
      int materialMapBinsZ     = int( maxZ/step);
      int materialMapBinsR     = int( maxR/step);
      double materialMapRangeZ = materialMapBinsZ*step;
      double materialMapRangeR = materialMapBinsR*step;

      m_radMap[trkName].SetNameTitle("radMap", std::string(trkName +" Radiation length map;z [mm];r [mm]").c_str());
      m_radMap[trkName].SetBins( materialMapBinsZ, 0, materialMapRangeZ, materialMapBinsR, 0, materialMapRangeR);
      m_radMapCount[trkName].SetNameTitle("radMapCount", std::string(trkName +" Radiation length map (counter);z [mm];r [mm]").c_str());
      m_radMapCount[trkName].SetBins( materialMapBinsZ, 0, materialMapRangeZ, materialMapBinsR, 0, materialMapRangeR);
      m_intMap[trkName].SetNameTitle("intMap", std::string(trkName +" Interaction length map;z [mm];r [mm]").c_str());
      m_intMap[trkName].SetBins( materialMapBinsZ, 0, materialMapRangeZ, materialMapBinsR, 0, materialMapRangeR);
      m_intMapCount[trkName].SetNameTitle("intMapCount", std::string(trkName +" Interaction length map (counter);z [mm];r [mm]").c_str());
      m_intMapCount[trkName].SetBins( materialMapBinsZ, 0, materialMapRangeZ, materialMapBinsR, 0, materialMapRangeR);

      // Material distribution in barrels
      m_radMB[trkName]["Barrel"].Reset();
      m_radMB[trkName]["Barrel"].SetNameTitle(std::string("radMBBarrel"+trkName).c_str(), std::string(trkName +" Barrel Modules Radiation Length").c_str());
      m_radMB[trkName]["Barrel"].SetBins(m_nTracks, 0, m_etaMax);
      m_intMB[trkName]["Barrel"].Reset();
      m_intMB[trkName]["Barrel"].SetNameTitle(std::string("intMBBarrel"+trkName).c_str(), std::string(trkName +" Barrel Modules Interaction Length").c_str());
      m_intMB[trkName]["Barrel"].SetBins(m_nTracks, 0, m_etaMax);

      // Material distribution in endcaps
      m_radMB[trkName]["Endcap"].Reset();
      m_radMB[trkName]["Endcap"].SetNameTitle(std::string("radMBEndcap"+trkName).c_str(), std::string(trkName +" Endcap Modules Radiation Length").c_str());
      m_radMB[trkName]["Endcap"].SetBins(m_nTracks, 0, m_etaMax);
      m_intMB[trkName]["Endcap"].Reset();
      m_intMB[trkName]["Endcap"].SetNameTitle(std::string("intMBEndcap"+trkName).c_str(), std::string(trkName +" Endcap Modules Interaction Length").c_str());
      m_intMB[trkName]["Endcap"].SetBins(m_nTracks, 0, m_etaMax);

      // Material distribution in services
      m_radMB[trkName]["Services"].Reset();
      m_radMB[trkName]["Services"].SetNameTitle(std::string("radMBServices"+trkName).c_str(), std::string(trkName +" Services Radiation Length").c_str());
      m_radMB[trkName]["Services"].SetBins(m_nTracks, 0, m_etaMax);
      m_intMB[trkName]["Services"].Reset();
      m_intMB[trkName]["Services"].SetNameTitle(std::string("intMBServices"+trkName).c_str(), std::string(trkName +" Services Interaction Length").c_str());
      m_intMB[trkName]["Services"].SetBins(m_nTracks, 0, m_etaMax);

      // Material distribution in supports
      m_radMB[trkName]["Supports"].Reset();
      m_radMB[trkName]["Supports"].SetNameTitle(std::string("radMBSupports"+trkName).c_str(), std::string(trkName +" Supports Radiation Length").c_str());
      m_radMB[trkName]["Supports"].SetBins(m_nTracks, 0, m_etaMax);
      m_intMB[trkName]["Supports"].Reset();
      m_intMB[trkName]["Supports"].SetNameTitle(std::string("intMBSupports"+trkName).c_str(), std::string(trkName +" Supports Interaction Length").c_str());
      m_intMB[trkName]["Supports"].SetBins(m_nTracks, 0, m_etaMax);

      // Material distribution in total
      m_radMB[trkName]["Total"].Reset();
      m_radMB[trkName]["Total"].SetNameTitle(std::string("radMBTotal"+trkName).c_str(), "Total Radiation Length");
      m_radMB[trkName]["Total"].SetBins(m_nTracks, 0, m_etaMax);
      m_intMB[trkName]["Total"].Reset();
      m_intMB[trkName]["Total"].SetNameTitle(std::string("intMBTotal"+trkName).c_str(), "Total Interaction Length");
      m_intMB[trkName]["Total"].SetBins(m_nTracks, 0, m_etaMax);

      // Nuclear interactions
      m_hadronTotalHitsGraph[trkName].SetName(std::string("HadronTotalHitsGraphIn"+trkName).c_str());
      m_hadronAverageHitsGraph[trkName].SetName(std::string("HadronAverageHitsGraphIn"+trkName).c_str());

      // Clear the list of requested good hadron hits
      m_hadronNeededHitsFraction[trkName].clear();
      m_hadronGoodTracksFraction[trkName].clear();

      m_hadronNeededHitsFraction[trkName].push_back(0);
      m_hadronNeededHitsFraction[trkName].push_back(0.0001);
      //m_hadronNeededHitsFraction[trkName].push_back(.33);
      m_hadronNeededHitsFraction[trkName].push_back(.66);
      m_hadronNeededHitsFraction[trkName].push_back(1);

      std::sort(m_hadronNeededHitsFraction[trkName].begin(),
                m_hadronNeededHitsFraction[trkName].end());

      // Prepare the plots for the track survival fraction
      ostringstream tempSS;
      for (auto i=0; i<m_hadronNeededHitsFraction[trkName].size(); ++i) {
        tempSS.str("");
        tempSS << "hadronGoodTracksFraction" << trkName << "_at" << m_hadronNeededHitsFraction[trkName].at(i);
        TGraph myGraph;
        myGraph.SetName(tempSS.str().c_str());
        m_hadronGoodTracksFraction[trkName].push_back(myGraph);
      }
    }  // For trackers

    // Material distribution in beam-pipe
    m_radMB["Beampipe"]["Total"].Reset();
    m_radMB["Beampipe"]["Total"].SetNameTitle("radMBPipe", "BeamPipe Radiation Length");
    m_radMB["Beampipe"]["Total"].SetBins(m_nTracks, 0, m_etaMax);
    m_intMB["Beampipe"]["Total"].Reset();
    m_intMB["Beampipe"]["Total"].SetNameTitle("intMBPipe", "BeamPipe Interaction Length");
    m_intMB["Beampipe"]["Total"].SetBins(m_nTracks, 0, m_etaMax);

//    // Nuclear interactions
//    m_hadronTotalHitsGraph["Tracker"].SetName(std::string("HadronTotalHitsGraphInTracker").c_str());
//    m_hadronAverageHitsGraph["Tracker"].SetName(std::string("HadronAverageHitsGraphInTracker").c_str());
//
//    // Clear the list of requested good hadron hits
//    m_hadronNeededHitsFraction["Tracker"].clear();
//    m_hadronGoodTracksFraction["Tracker"].clear();
//
//    m_hadronNeededHitsFraction["Tracker"].push_back(0);
//    m_hadronNeededHitsFraction["Tracker"].push_back(0.0001);
//    //m_hadronNeededHitsFraction["Tracker"].push_back(.33);
//    m_hadronNeededHitsFraction["Tracker"].push_back(.66);
//    m_hadronNeededHitsFraction["Tracker"].push_back(1);
//
//    std::sort(m_hadronNeededHitsFraction["Tracker"].begin(),
//              m_hadronNeededHitsFraction["Tracker"].end());
//
//    // Prepare the plots for the track survival fraction
//    ostringstream tempSS;
//    for (auto i=0; i<m_hadronNeededHitsFraction["Tracker"].size(); ++i) {
//      tempSS.str("");
//      tempSS << "hadronGoodTracksFraction_at" << m_hadronNeededHitsFraction["Tracker"].at(i);
//      TGraph myGraph;
//      myGraph.SetName(tempSS.str().c_str());
//      m_hadronGoodTracksFraction["Tracker"].push_back(myGraph);
//    }

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

  TRandom3 myDice;
  myDice.SetSeed(random_seed);

  // Print
  std::cout << std::endl;

  // Study all trackers
  for (auto iTracker : m_trackers) {

    std::string trkName = iTracker->myid();
    std::cout << " " << trkName << " tracker" << std::endl;

    // Analyze material budget for individual subdetector & beam-pipe
    for (auto iTrack=0; iTrack<m_nTracks; iTrack++) {

      // Generate a track
      double phi   = myDice.Rndm() * 2 * M_PI;
      double eta   = 0 + m_etaMax/m_nTracks*(iTrack+0.5); // Currently assumed +-Z symmetry, to fill bin centre use iTrack + 0.5 // m_etaMin + m_etaSpan/m_nTracks*iTrack;
      double theta = 2*atan(exp(-1*eta));
      double pT    = 100*Units::TeV; // Arbitrarily high number

      Track matTrack;
      matTrack.setThetaPhiPt(theta, phi, pT);
      matTrack.setOrigin(0, 0, 0); // TODO: Not assuming z-error when calculating Material budget

      // Study beam pipe -> fill only once
      if (iTracker==*m_trackers.begin()) {

        m_radMB["Beampipe"]["Total"].Fill(eta, m_beamPipe->radLength()/sin(theta));
        m_intMB["Beampipe"]["Total"].Fill(eta, m_beamPipe->intLength()/sin(theta));

        m_radMB["Tracker"]["Total"].Fill(eta, m_beamPipe->radLength()/sin(theta));
        m_intMB["Tracker"]["Total"].Fill(eta, m_beamPipe->intLength()/sin(theta));
      }

      //
      // Get material related to detector modules, so-called module caps, beam-pipe etc.
      std::map<std::string, Material> matBudget;
      MatBudgetVisitor matVisitor(matTrack, matBudget, m_radMap[trkName], m_radMapCount[trkName], m_intMap[trkName], m_intMapCount[trkName]);
      iTracker->accept(matVisitor);  // Assign to material track hits corresponding to modules
      m_beamPipe->accept(matVisitor);// Assign to material track hit corresponding to beam-pipe

      // Given sub-tracker
      m_radMB[trkName]["Barrel"].Fill(eta, matBudget["Barrel"].radiation);
      m_intMB[trkName]["Barrel"].Fill(eta, matBudget["Barrel"].interaction);

      m_radMB[trkName]["Endcap"].Fill(eta, matBudget["Endcap"].radiation);
      m_intMB[trkName]["Endcap"].Fill(eta, matBudget["Endcap"].interaction);

      m_radMB[trkName]["Supports"].Fill(eta, matBudget["Supports"].radiation);
      m_intMB[trkName]["Supports"].Fill(eta, matBudget["Supports"].interaction);

      m_radMB[trkName]["Services"].Fill(eta, matBudget["Services"].radiation);
      m_intMB[trkName]["Services"].Fill(eta, matBudget["Services"].interaction);

      m_radMB[trkName]["Total"].Fill(eta, matBudget["Barrel"].radiation +
                                          matBudget["Endcap"].radiation +
                                          matBudget["Supports"].radiation +
                                          matBudget["Services"].radiation +
                                          m_beamPipe->radLength()/sin(theta));
      m_intMB[trkName]["Total"].Fill(eta, matBudget["Barrel"].interaction +
                                          matBudget["Endcap"].interaction +
                                          matBudget["Supports"].interaction +
                                          matBudget["Services"].interaction +
                                          m_beamPipe->intLength()/sin(theta));

      // Total sum to tracker
      m_radMB["Tracker"]["Barrel"].Fill(eta, matBudget["Barrel"].radiation);
      m_intMB["Tracker"]["Barrel"].Fill(eta, matBudget["Barrel"].interaction);

      m_radMB["Tracker"]["Endcap"].Fill(eta, matBudget["Endcap"].radiation);
      m_intMB["Tracker"]["Endcap"].Fill(eta, matBudget["Endcap"].interaction);

      m_radMB["Tracker"]["Supports"].Fill(eta, matBudget["Supports"].radiation);
      m_intMB["Tracker"]["Supports"].Fill(eta, matBudget["Supports"].interaction);

      m_radMB["Tracker"]["Services"].Fill(eta, matBudget["Services"].radiation);
      m_intMB["Tracker"]["Services"].Fill(eta, matBudget["Services"].interaction);

      m_radMB["Tracker"]["Total"].Fill(eta, matBudget["Barrel"].radiation +
                                            matBudget["Endcap"].radiation +
                                            matBudget["Supports"].radiation +
                                            matBudget["Services"].radiation);
      m_intMB["Tracker"]["Total"].Fill(eta, matBudget["Barrel"].interaction +
                                            matBudget["Endcap"].interaction +
                                            matBudget["Supports"].interaction +
                                            matBudget["Services"].interaction);

      // Create component contribution
      for (auto iMap : matBudget) {

        std::string compName = iMap.first;

        // Avoid Barrel or Endcap -> interested in Modules components + Services + Supports
        if (compName!="Barrel" && compName!="Endcap") {

          // Set histogram binning etc. if not yet initialized
          if (m_radMBComp[trkName].find(compName)==m_radMBComp[trkName].end()) {

            m_radMBComp[trkName][compName].Reset();
            m_radMBComp[trkName][compName].SetNameTitle(std::string("radMB"+compName).c_str(), std::string(trkName +" "+compName+" Radiation Length").c_str());
            m_radMBComp[trkName][compName].SetBins(m_nTracks, 0, m_etaMax);

            m_intMBComp[trkName][compName].Reset();
            m_intMBComp[trkName][compName].SetNameTitle(std::string("intMB"+compName).c_str(), std::string(trkName +" "+compName+" Interaction Length").c_str());
            m_intMBComp[trkName][compName].SetBins(m_nTracks, 0, m_etaMax);
          }
          if (m_radMBComp["Tracker"].find(compName)==m_radMBComp["Tracker"].end()) {

            m_radMBComp["Tracker"][compName].Reset();
            m_radMBComp["Tracker"][compName].SetNameTitle(std::string("radMB"+compName).c_str(), std::string("Tracker "+compName+" Radiation Length").c_str());
            m_radMBComp["Tracker"][compName].SetBins(m_nTracks, 0, m_etaMax);

            m_intMBComp["Tracker"][compName].Reset();
            m_intMBComp["Tracker"][compName].SetNameTitle(std::string("intMB"+compName).c_str(), std::string("Tracker "+compName+" Interaction Length").c_str());
            m_intMBComp["Tracker"][compName].SetBins(m_nTracks, 0, m_etaMax);
          }

          // Fill
          m_radMBComp[trkName][compName].Fill(eta, iMap.second.radiation);
          m_intMBComp[trkName][compName].Fill(eta, iMap.second.interaction);

          m_radMBComp["Tracker"][compName].Fill(eta, iMap.second.radiation);
          m_intMBComp["Tracker"][compName].Fill(eta, iMap.second.interaction);
        }
      }

      // Calculate efficiency plots
      if (!matTrack.hasNoHits()) {

        bool bySmallerR = true;
        matTrack.sortHits(bySmallerR);

        // Get hits
        int nActive = matTrack.getNActiveHits("all",true);
        if (nActive>0) {

          m_hadronTotalHitsGraph[trkName].SetPoint(m_hadronTotalHitsGraph[trkName].GetN(), eta, nActive);

          double probability;
          std::vector<double> probabilities = matTrack.getHadronActiveHitsProbability("all");

          double averageHits  = 0;
          double exactProb    = 0;
          double moreThanProb = 0;
          for (int i=probabilities.size()-1; i>=0; --i) {

            //if (nActive==10) { // debug
            //  std::cerr << "probabilities.at("
            //  << i << ")=" << probabilities.at(i)
            //  << endl;
            //}

            exactProb     = probabilities.at(i)-moreThanProb;
            averageHits  += (i+1)*exactProb;
            moreThanProb += exactProb;
          }

          m_hadronAverageHitsGraph[trkName].SetPoint(m_hadronAverageHitsGraph[trkName].GetN(), eta, averageHits);
          //m_hadronAverageHitsGraph.SetPointError(_hadronAverageHitsGraph.GetN()-1, 0, sqrt( averageSquaredHits - averageHits*averageHits) );

          unsigned int requiredHits;
          for (unsigned int i = 0; i<m_hadronNeededHitsFraction[trkName].size(); ++i) {

            requiredHits = int(ceil(double(nActive) * m_hadronNeededHitsFraction[trkName].at(i)));
            if      (requiredHits==0)                   probability = 1;
            else if (requiredHits>probabilities.size()) probability = 0;
            else                                        probability = probabilities.at(requiredHits-1);

            //if (probabilities.size()==10) { // debug
            //  std::cerr << "required " << requiredHits
            //              << " out of " << probabilities.size()
            //              << " == " << nActive
            //              << endl;
            // std::cerr << "      PROBABILITY = " << probability << endl << endl;
            //}
            m_hadronGoodTracksFraction[trkName].at(i).SetPoint(m_hadronGoodTracksFraction[trkName].at(i).GetN(), eta, probability);
          }
        }
      } // Calculate efficiency plots

    } // For nTracks

    //
    // Normalization calculation
    int nBins = m_radMap[trkName].GetNbinsX()*m_radMap[trkName].GetNbinsY();

    for (int iBin=1; iBin<=nBins; ++iBin) {

      // Update radiation map
      double content = m_radMap[trkName].GetBinContent(iBin);
      double count   = m_radMapCount[trkName].GetBinContent(iBin);

      if      (count==1) m_radMap[trkName].SetBinContent(iBin,content);
      else if (count>1)  m_radMap[trkName].SetBinContent(iBin,content/double(count));

      m_radMap["Tracker"].SetBinContent(iBin, m_radMap["Tracker"].GetBinContent(iBin) + m_radMap[trkName].GetBinContent(iBin));

      // Update interaction map (the same number of bins)
      content = m_intMap[trkName].GetBinContent(iBin);
      count   = m_intMapCount[trkName].GetBinContent(iBin);

      if      (count==1) m_intMap[trkName].SetBinContent(iBin,content);
      else if (count>1)  m_intMap[trkName].SetBinContent(iBin,content/double(count));

      m_intMap["Tracker"].SetBinContent(iBin, m_intMap["Tracker"].GetBinContent(iBin) + m_intMap[trkName].GetBinContent(iBin));
    }

  } // For nTrackers

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
  int webPriority         = web_priority_MB;

  // Set Rainbow palette for drawing
  Palette::setRootPalette(55);

  for (int iTrk=0; iTrk<=m_trackers.size(); iTrk++) {

    // Tracker name
    const Tracker* tracker = nullptr;
    std::string trkName;
    if (iTrk==m_trackers.size()) {

      tracker = nullptr;
      trkName = "Tracker";
    }
    else {

      tracker = m_trackers[iTrk];
      trkName = m_trackers[iTrk]->myid();
    }

    // Create dedicated web-page & its content
    std::string        pageTitle  = "MatBudget";
    if (trkName != "") pageTitle +=" (" + trkName + ")";

    std::string pageAddress = "indexMB"+trkName+".html";

    RootWPage& myPage = webSite.addPage(pageTitle, webPriority);
    myPage.setAddress(pageAddress);
    webPriority--;

    //
    // Material overview
    RootWContent& myContentOverview = myPage.addContent("Material overview", true);

    // Prepare canvas
    TCanvas myCanvasMat(std::string("MaterialInTrackingVolume"+trkName).c_str());
    myCanvasMat.SetFillColor(Palette::color_plot_background);
    myCanvasMat.Divide(2, 1);
    TPad* myPad = dynamic_cast<TPad*>(myCanvasMat.GetPad(0));
    myPad->SetFillColor(Palette::color_pad_background);

    // Rebin material histograms to readable values
    m_radMB[trkName]["Total"].SetFillColor(kGray + 2);
    m_radMB[trkName]["Total"].SetLineColor(kBlue + 2);
    m_radMB[trkName]["Total"].SetXTitle("#eta");
    double max = m_radMB[trkName]["Total"].GetMaximum();
    m_radMB[trkName]["Total"].GetYaxis()->SetRangeUser(0, vis_safety_factor*max);
    myPad = dynamic_cast<TPad*>(myCanvasMat.GetPad(1));
    myPad->cd();
    m_radMB[trkName]["Total"].SetStats(kFALSE);
    m_radMB[trkName]["Total"].Draw();

    m_intMB[trkName]["Total"].SetFillColor(kGray + 2);
    m_intMB[trkName]["Total"].SetLineColor(kBlue + 2);
    m_intMB[trkName]["Total"].SetXTitle("#eta");
    max = m_intMB[trkName]["Total"].GetMaximum();
    m_intMB[trkName]["Total"].GetYaxis()->SetRangeUser(0, vis_safety_factor*max);
    myPad = dynamic_cast<TPad*>(myCanvasMat.GetPad(2));
    myPad->cd();
    m_intMB[trkName]["Total"].SetStats(kFALSE);
    m_intMB[trkName]["Total"].Draw();

    // Write tracking volume plots to the web site
    RootWImage& myImageOverview = myContentOverview.addImage(myCanvasMat, 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImageOverview.setComment(std::string("Material overview in tracking volume: "+trkName));
    myImageOverview.setName("matOverview");

    RootWTable& myTableOverview = myContentOverview.addTable();
    char titleString[256];
    sprintf(titleString, std::string("Average radiation length [%] in tracking volume ("+web_etaLetter+" = [0, %.1f])").c_str(), SimParms::getInstance().getMaxEtaCoverage());
    myTableOverview.setContent(1, 1, titleString);
    sprintf(titleString, std::string("Average interaction length [%] in tracking volume ("+web_etaLetter+" = [0, %.1f])").c_str(), SimParms::getInstance().getMaxEtaCoverage());
    myTableOverview.setContent(2, 1, titleString);
    myTableOverview.setContent(1, 2, averageHistogramValues(m_radMB[trkName]["Total"], SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableOverview.setContent(2, 2, averageHistogramValues(m_intMB[trkName]["Total"], SimParms::getInstance().getMaxEtaCoverage())*100, 2);

    // Calculate summary table
    std::unique_ptr<RootWTable> materialSummaryTable(new RootWTable());

    double averageValue;
    materialSummaryTable->setContent(0,0,"Material    ");
    materialSummaryTable->setContent(1,0,"Rad. length [%]");
    materialSummaryTable->setContent(2,0,"Int. length [%]");
    materialSummaryTable->setContent(3,0,"Photon conv. prob.");

    // Material csv file
    if (!m_csvMatBudget->existCsvText("Label")) {

      m_csvMatBudget->addCsvElement("Label", "Tracker/Component name");
      m_csvMatBudget->addCsvElement("Label", "Type");
      for (unsigned int j=1; j< SimParms::getInstance().getNEtaRegions(); ++j) {

        ostringstream label;
        label << "eta(" << std::fixed << std::setprecision(1) << SimParms::getInstance().etaRegionRanges[j-1] << "-" << SimParms::getInstance().etaRegionRanges[j] << ")";
        m_csvMatBudget->addCsvElement("Label", label.str());
      }
      m_csvMatBudget->addCsvEOL("Label");
    }

    m_csvMatBudget->addCsvElement(trkName, trkName);
    m_csvMatBudget->addCsvElement(trkName, "Rad. length [%]");
    for (unsigned int j=1; j< SimParms::getInstance().getNEtaRegions(); ++j) {
      // First row: the cut name
      materialSummaryTable->setContent(0,j, SimParms::getInstance().etaRegionNames[j]);

      // Second row: the radiation length
      averageValue = averageHistogramValues(m_radMB[trkName]["Total"], SimParms::getInstance().etaRegionRanges[j-1], SimParms::getInstance().etaRegionRanges[j]);
      materialSummaryTable->setContent(1,j, averageValue*100 ,2);
      m_csvMatBudget->addCsvElement(trkName, averageValue*100);
    }
    m_csvMatBudget->addCsvEOL(trkName);

    m_csvMatBudget->addCsvElement(trkName, "");
    m_csvMatBudget->addCsvElement(trkName, "Int. length [%]");
    for (unsigned int j=1; j< SimParms::getInstance().getNEtaRegions(); ++j) {
      // Third row: the interaction length
      averageValue = averageHistogramValues(m_intMB[trkName]["Total"], SimParms::getInstance().etaRegionRanges[j-1], SimParms::getInstance().etaRegionRanges[j]);
      materialSummaryTable->setContent(2,j, averageValue*100 ,2);
      m_csvMatBudget->addCsvElement(trkName, averageValue*100);
    }
    m_csvMatBudget->addCsvEOL(trkName);

    m_csvMatBudget->addCsvElement(trkName, "");
    m_csvMatBudget->addCsvElement(trkName, "Photon conv. prob.");
    for (unsigned int j=1; j< SimParms::getInstance().getNEtaRegions(); ++j) {
      // Fourth row: the photon conversion probability
      averageValue  = averageHistogramValues(m_radMB[trkName]["Total"], SimParms::getInstance().etaRegionRanges[j-1], SimParms::getInstance().etaRegionRanges[j]);
      averageValue *= -7./9.;
      averageValue  = 1 - exp(averageValue);
      materialSummaryTable->setContent(3,j, averageValue ,4);
      m_csvMatBudget->addCsvElement(trkName, averageValue);
    }
    m_csvMatBudget->addCsvEOL(trkName);

    //
    // Detailed tracker material overview
    RootWContent& myContentCateg = myPage.addContent("Material overview by category", true);

    // Set variables and book histograms
    std::unique_ptr<THStack> radContainer(new THStack("rstack", "Radiation Length by Category"));
    std::unique_ptr<THStack> intContainer(new THStack("istack", "Interaction Length by Category"));

    // Radiation length in tracking volume by active, serving or passive
    m_radMB["Beampipe"]["Total"].SetFillColor(kGreen);
    m_radMB["Beampipe"]["Total"].SetLineColor(kGreen);
    m_radMB["Beampipe"]["Total"].SetXTitle("#eta");
    m_radMB[trkName]["Barrel"].SetFillColor(TColor::GetColor("#FFD21F"));
    m_radMB[trkName]["Barrel"].SetLineColor(TColor::GetColor("#FFD21F"));
    m_radMB[trkName]["Barrel"].SetXTitle("#eta");
    m_radMB[trkName]["Endcap"].SetFillColor(kRed);
    m_radMB[trkName]["Endcap"].SetLineColor(kRed);
    m_radMB[trkName]["Endcap"].SetXTitle("#eta");
    m_radMB[trkName]["Supports"].SetFillColor(kOrange+4);
    m_radMB[trkName]["Supports"].SetLineColor(kOrange+4);
    m_radMB[trkName]["Supports"].SetXTitle("#eta");
    m_radMB[trkName]["Services"].SetFillColor(kBlue);
    m_radMB[trkName]["Services"].SetLineColor(kBlue);
    m_radMB[trkName]["Services"].SetXTitle("#eta");

    // Interaction length in tracking volume by active, serving or passive
    m_intMB["Beampipe"]["Total"].SetFillColor(kGreen-2);
    m_intMB["Beampipe"]["Total"].SetLineColor(kGreen-2);
    m_intMB["Beampipe"]["Total"].SetXTitle("#eta");
    m_intMB[trkName]["Barrel"].SetFillColor(TColor::GetColor("#FFDE5C"));
    m_intMB[trkName]["Barrel"].SetLineColor(TColor::GetColor("#FFDE5C"));
    m_intMB[trkName]["Barrel"].SetXTitle("#eta");
    m_intMB[trkName]["Endcap"].SetFillColor(kRed+1);
    m_intMB[trkName]["Endcap"].SetLineColor(kRed+1);
    m_intMB[trkName]["Endcap"].SetXTitle("#eta");
    m_intMB[trkName]["Supports"].SetFillColor(kOrange+2);
    m_intMB[trkName]["Supports"].SetLineColor(kOrange+2);
    m_intMB[trkName]["Supports"].SetXTitle("#eta");
    m_intMB[trkName]["Services"].SetFillColor(kAzure-2);
    m_intMB[trkName]["Services"].SetLineColor(kAzure-2);
    m_intMB[trkName]["Services"].SetXTitle("#eta");

    // Write all material information into web page table
    RootWTable& myTableAllMat = myContentCateg.addTable();

    // Average values by active, service and passive
    sprintf(titleString, std::string("Average ("+web_etaLetter+" = [0, %.1f])").c_str(), SimParms::getInstance().getMaxEtaCoverage());
    myTableAllMat.setContent(0, 0, titleString);
    myTableAllMat.setContent(1, 0, "Beam pipe (green)");
    myTableAllMat.setContent(2, 0, "Barrel modules (yellow)");
    myTableAllMat.setContent(3, 0, "Endcap modules (red)");
    myTableAllMat.setContent(4, 0, "Supports (brown)");
    myTableAllMat.setContent(5, 0, "Services (blue)");
    myTableAllMat.setContent(6, 0, "Total");
    myTableAllMat.setContent(0, 1, "Radiation length [%]");
    myTableAllMat.setContent(0, 2, "Interaction length [%]");
    myTableAllMat.setContent(1, 1, averageHistogramValues(m_radMB["Beampipe"]["Total"], SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableAllMat.setContent(1, 2, averageHistogramValues(m_intMB["Beampipe"]["Total"], SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableAllMat.setContent(2, 1, averageHistogramValues(m_radMB[trkName]["Barrel"]  , SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableAllMat.setContent(2, 2, averageHistogramValues(m_intMB[trkName]["Barrel"]  , SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableAllMat.setContent(3, 1, averageHistogramValues(m_radMB[trkName]["Endcap"]  , SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableAllMat.setContent(3, 2, averageHistogramValues(m_intMB[trkName]["Endcap"]  , SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableAllMat.setContent(4, 1, averageHistogramValues(m_radMB[trkName]["Supports"], SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableAllMat.setContent(4, 2, averageHistogramValues(m_intMB[trkName]["Supports"], SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableAllMat.setContent(5, 1, averageHistogramValues(m_radMB[trkName]["Services"], SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableAllMat.setContent(5, 2, averageHistogramValues(m_intMB[trkName]["Services"], SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableAllMat.setContent(6, 1, averageHistogramValues(m_radMB[trkName]["Total"]   , SimParms::getInstance().getMaxEtaCoverage())*100, 2);
    myTableAllMat.setContent(6, 2, averageHistogramValues(m_intMB[trkName]["Total"]   , SimParms::getInstance().getMaxEtaCoverage())*100, 2);

    // Rebin histograms, draw them to a canvas and write the canvas to the web page
    TCanvas myCanvasMatCat(std::string("MaterialByCategoryIn"+trkName).c_str());
    myCanvasMatCat.SetFillColor(Palette::color_plot_background);
    myCanvasMatCat.Divide(2, 1);
    myPad = dynamic_cast<TPad*>(myCanvasMatCat.GetPad(0));
    myPad->SetFillColor(Palette::color_pad_background);

    myPad = dynamic_cast<TPad*>(myCanvasMatCat.GetPad(1));
    myPad->cd();
    radContainer->Add(&m_radMB["Beampipe"]["Total"]);
    radContainer->Add(&m_radMB[trkName]["Barrel"]);
    radContainer->Add(&m_radMB[trkName]["Endcap"]);
    radContainer->Add(&m_radMB[trkName]["Supports"]);
    radContainer->Add(&m_radMB[trkName]["Services"]);
    radContainer->Draw();
    radContainer->GetXaxis()->SetTitle("#eta");

    myPad = dynamic_cast<TPad*>(myCanvasMatCat.GetPad(2));
    myPad->cd();
    intContainer->Add(&m_intMB["Beampipe"]["Total"]);
    intContainer->Add(&m_intMB[trkName]["Barrel"]);
    intContainer->Add(&m_intMB[trkName]["Endcap"]);
    intContainer->Add(&m_intMB[trkName]["Supports"]);
    intContainer->Add(&m_intMB[trkName]["Services"]);
    intContainer->Draw();
    intContainer->GetXaxis()->SetTitle("#eta");

    RootWImage& myImageMatCat = myContentCateg.addImage(myCanvasMatCat, 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImageMatCat.setComment("Material overview by category in tracking volume: "+trkName);
    myImageMatCat.setName("riDistrCateg");

    //
    // Detailed material overview by component
    RootWContent& myContentComp = myPage.addContent("Material overview by component", false);

    RootWTable myTableComp = myContentComp.addTable();
    sprintf(titleString, std::string("Average ("+web_etaLetter+" = [0, %.1f])").c_str(), SimParms::getInstance().getMaxEtaCoverage());
    myTableComp.setContent(0, 0, titleString);
    myTableComp.setContent(0, 1, "Radiation length [%]");
    myTableComp.setContent(0, 2, "Interaction length [%]");

    // Set variables & book histograms
    THStack* radCompStack = new THStack("rcompstack", "Radiation Length by Component");
    THStack* intCompStack = new THStack("icompstack", "Interaction Length by Component");

    TLegend* compLegend = new TLegend(0.1,0.6,0.35,0.9);

    TCanvas myCanvasMatComp(std::string("MaterialByComponentIn"+trkName).c_str());
    myCanvasMatComp.SetFillColor(Palette::color_plot_background);
    myCanvasMatComp.Divide(2, 1);
    myPad = dynamic_cast<TPad*>(myCanvasMatComp.GetPad(0));
    myPad->SetFillColor(Palette::color_pad_background);

    // Radiation length for components
    myPad = dynamic_cast<TPad*>(myCanvasMatComp.GetPad(1));
    myPad->cd();
    int    compIndex      = 1;
    double totalRadLength = 0;

    m_radMB["Beampipe"]["Total"].SetLineColor(Palette::color(compIndex));
    m_radMB["Beampipe"]["Total"].SetFillColor(Palette::color(compIndex));
    m_radMB["Beampipe"]["Total"].SetXTitle("#eta");
    radCompStack->Add(&m_radMB["Beampipe"]["Total"]);
    compLegend->AddEntry(&m_radMB["Beampipe"]["Total"], "Beampipe");
    compIndex++;
    double avgValue = averageHistogramValues(m_radMB["Beampipe"]["Total"], SimParms::getInstance().getMaxEtaCoverage());
    myTableComp.setContent(compIndex, 0, "Beampipe");
    myTableComp.setContent(compIndex++, 1, avgValue*100, 2);
    totalRadLength += avgValue;

    for (auto& iComp : m_radMBComp[trkName]) {

      iComp.second.SetLineColor(Palette::color(compIndex));
      iComp.second.SetFillColor(Palette::color(compIndex));
      iComp.second.SetXTitle("#eta");
      radCompStack->Add(&(iComp.second));
      compLegend->AddEntry(&(iComp.second), iComp.first.c_str());
      avgValue = averageHistogramValues(iComp.second, SimParms::getInstance().getMaxEtaCoverage());
      myTableComp.setContent(compIndex, 0, iComp.first);
      myTableComp.setContent(compIndex++, 1, avgValue*100, 2);
      totalRadLength += avgValue;
    }
    myTableComp.setContent(compIndex, 0, "Total");
    myTableComp.setContent(compIndex, 1, totalRadLength*100, 2);
    radCompStack->Draw();
    radCompStack->GetXaxis()->SetTitle("#eta");
    compLegend->Draw();

    // Interaction length for components
    myPad = dynamic_cast<TPad*>(myCanvasMatComp.GetPad(2));
    myPad->cd();
    compIndex             = 1;
    double totalIntLength = 0;

    m_intMB["Beampipe"]["Total"].SetLineColor(Palette::color(compIndex));
    m_intMB["Beampipe"]["Total"].SetFillColor(Palette::color(compIndex));
    m_intMB["Beampipe"]["Total"].SetXTitle("#eta");
    intCompStack->Add(&m_intMB["Beampipe"]["Total"]);
    compIndex++;
    avgValue = averageHistogramValues(m_intMB["Beampipe"]["Total"], SimParms::getInstance().getMaxEtaCoverage());
    myTableComp.setContent(compIndex++, 2, avgValue*100, 2);
    totalIntLength += avgValue;

    for (auto& iComp : m_intMBComp[trkName]) {

      iComp.second.SetLineColor(Palette::color(compIndex));
      iComp.second.SetFillColor(Palette::color(compIndex));
      iComp.second.SetXTitle("#eta");
      intCompStack->Add(&(iComp.second));
      avgValue = averageHistogramValues(iComp.second, SimParms::getInstance().getMaxEtaCoverage());
      myTableComp.setContent(compIndex++, 2, avgValue*100, 2);
      totalIntLength += avgValue;
    }
    myTableComp.setContent(compIndex, 2, totalIntLength*100, 2);
    intCompStack->Draw();
    intCompStack->GetXaxis()->SetTitle("#eta");
    compLegend->Draw();

    RootWImage& myImageMatComp = myContentComp.addImage(myCanvasMatComp, 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImageMatComp.setComment(std::string("Material overview by component in tracking volume: "+trkName));
    myImageMatComp.setName("riDistrComp");

    //
    // Material map
    RootWContent& myContentMap = myPage.addContent("Material map", true);

    // Calculate final plot range, increase by safety factor
    double maxZ = 0;
    double maxR = 0;

    if (tracker!=nullptr) {

      maxZ = tracker->maxZ()*vis_safety_factor;
      maxR = tracker->maxR()*vis_safety_factor;
    }
    else {

      for (auto jTrk : m_trackers) {

        MAX(maxZ, jTrk->maxZ()*vis_safety_factor);
        MAX(maxR, jTrk->maxR()*vis_safety_factor);
      }
    }

    // Radiation length plot
    TCanvas myCanvasMapRad(std::string("RadMaterialMapIn"+trkName).c_str());
    myCanvasMapRad.SetFillColor(Palette::color_plot_background);
    myCanvasMapRad.cd();

    m_radMap[trkName].GetXaxis()->SetRangeUser(0, maxZ);
    m_radMap[trkName].GetYaxis()->SetRangeUser(0, maxR);
    m_radMap[trkName].SetContour(vis_temperature_levels, 0);
    m_radMap[trkName].GetYaxis()->SetTitleOffset(1.1);
    m_radMap[trkName].SetStats(kFALSE);
    m_radMap[trkName].Draw("COLZ");

    RootWImage& myImageMapRad = myContentMap.addImage(myCanvasMapRad, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImageMapRad.setComment("Radiation length material map");
    myImageMapRad.setName("matMapRad");

    // Interaction length plot
    TCanvas myCanvasMapInt(std::string("IntMaterialMapIn"+trkName).c_str());
    myCanvasMapInt.SetFillColor(Palette::color_plot_background);
    myCanvasMapInt.cd();

    m_intMap[trkName].GetXaxis()->SetRangeUser(0, maxZ);
    m_intMap[trkName].GetYaxis()->SetRangeUser(0, maxR);
    m_intMap[trkName].SetContour(vis_temperature_levels, 0);
    m_intMap[trkName].GetYaxis()->SetTitleOffset(1.1);
    m_intMap[trkName].SetStats(kFALSE);
    m_intMap[trkName].Draw("COLZ");

    RootWImage& myImageMapInt = myContentMap.addImage(myCanvasMapInt, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImageMapInt.setComment("Interaction length material map");
    myImageMapInt.setName("matMapInt");

    //
    // Hits occupancy: TODO Fix for full tracker
    if (trkName!="Tracker") {
      RootWContent& myContentOccup = myPage.addContent("Hits occupancy & track efficiency", false);

      // Set variables
      std::map<int, std::vector<double> > averages;

      // Number of hits
      TCanvas myCanvasOccup(std::string("HadronsHitsNumberIn"+trkName).c_str());
      myCanvasOccup.SetFillColor(Palette::color_plot_background);
      myCanvasOccup.Divide(2, 1);
      myPad = dynamic_cast<TPad*>(myCanvasOccup.GetPad(0));
      myPad->SetFillColor(Palette::color_pad_background);
      myPad = dynamic_cast<TPad*>(myCanvasOccup.GetPad(1));
      myPad->cd();

      m_hadronTotalHitsGraph[trkName].SetTitle("Maximum (black) & average (red) number of hits");
      m_hadronAverageHitsGraph[trkName].SetTitle("Maximum (black) & average (red) number of hits");
      m_hadronTotalHitsGraph[trkName].SetMarkerStyle(8);
      m_hadronTotalHitsGraph[trkName].SetMarkerColor(kBlack);
      m_hadronTotalHitsGraph[trkName].SetMinimum(0);
      m_hadronTotalHitsGraph[trkName].GetXaxis()->SetTitle("#eta");
      m_hadronTotalHitsGraph[trkName].Draw("alp");
      m_hadronAverageHitsGraph[trkName].SetMarkerStyle(8);
      m_hadronAverageHitsGraph[trkName].SetMarkerColor(kRed);
      m_hadronAverageHitsGraph[trkName].Draw("same lp");

      // Track fraction
      myPad = dynamic_cast<TPad*>(myCanvasOccup.GetPad(2));
      myPad->cd();
      std::unique_ptr<TLegend> myLegend(new TLegend(0.65, 0.16, .85, .40));
      // Old-style palette by Stefano, with custom-generated colors
      // Palette::prepare(hadronGoodTracksFraction.size()); // there was a 120 degree phase here
      // Replaced by the libreOffice-like palette
      std::unique_ptr<TH1D> ranger(new TH1D(std::string("HadTrackRangerIn"+trkName).c_str(),"Track efficiency with given fraction of hits ", 100, 0, SimParms::getInstance().getMaxEtaCoverage()));
      ranger->SetMaximum(1.);
      ranger->SetStats(kFALSE);
      ranger->GetXaxis()->SetTitle("#eta");
      //myAxis = ranger->GetYaxis();
      //myAxis->SetTitle("Tracks fraction");
      ranger->Draw();
      ostringstream tempSS;
      std::map<int, std::string> fractionTitles;

      // Prepare eta regions to vector
      std::vector<double> etaCuts;
      for (auto& iCut : SimParms::getInstance().etaRegionRanges) {
        etaCuts.push_back(iCut);
      }

      for (auto i=0; i<m_hadronGoodTracksFraction[trkName].size(); ++i) {

        TGraph& myGraph = m_hadronGoodTracksFraction[trkName].at(i);
        //std::cerr << "Good Hadrons fractions at (" << i <<") has " << myGraph.GetN() << " points" << std::endl;
        //double xx, yy;
        //myGraph.GetPoint(myGraph.GetN()-1, xx, yy);
        //std::cerr << "Last point (x,y) = ("<< xx <<", " << yy <<")" << std::endl;
        averages[i] = average(myGraph, etaCuts);
        closeGraph(myGraph);
        myGraph.SetFillColor(Palette::color(i+1));
        myGraph.Draw("same F");
        tempSS.str("");
        if (m_hadronNeededHitsFraction[trkName].at(i)!=0) {
          if (m_hadronNeededHitsFraction[trkName].at(i)==0.0001)
            tempSS << "1 hit required";
          else
            tempSS << int(m_hadronNeededHitsFraction[trkName].at(i)*100)
              << "% hits required";
          fractionTitles[i]=tempSS.str();
          myLegend->AddEntry(&myGraph, fractionTitles[i].c_str(), "F");
        }
      }
      ranger->Draw("sameaxis");
      myLegend->Draw();

      RootWImage& myImageOccup = myContentOccup.addImage(myCanvasOccup, 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImageOccup.setComment("Hits occupancy & track efficiency for hadrons");
      myImageOccup.setName("hadHitsTracks");
    }

    //
    // Summary table
    RootWContent& summaryContent = myPage.addContent("Summary", true);

    // Define variables
    ostringstream label;

    // Cuts summary table
    RootWTable&   cutsSummaryTable = summaryContent.addTable();
    cutsSummaryTable.setContent(0,0,"Region: ");
    cutsSummaryTable.setContent(1,0,"Min "+web_etaLetter+":");
    cutsSummaryTable.setContent(2,0,"Max "+web_etaLetter+":");

    RootWTable* myTable = &cutsSummaryTable;
    for (unsigned int iBorder=0; iBorder<SimParms::getInstance().getNEtaRegions()-1; ++iBorder) {
      myTable->setContent(0,iBorder+1,SimParms::getInstance().etaRegionNames[iBorder+1]);
      label.str(""); label << std::fixed << std::setprecision(1) << SimParms::getInstance().etaRegionRanges[iBorder];
      myTable->setContent(1,iBorder+1,label.str());
      label.str(""); label << std::fixed << std::setprecision(1) << SimParms::getInstance().etaRegionRanges[iBorder+1];
      myTable->setContent(2,iBorder+1,label.str());
    }

    // Calculate material bill
    //for (auto materialBudget : materialBudgets) {
    //  if (!m_materialBillCsv.existCsvText(materialBudget->getTracker().myid())) m_materialBillCsv.inspectTracker(*materialBudget);
    //}

    // Material summary table
    summaryContent.addItem(std::move(materialSummaryTable));

  } // Trackers

  m_isVisOK = true;
  return m_isVisOK;
}

//
// Helper method creating profile histogram from histogram
//
//TProfile* AnalyzerMatBudget::createProfileHis(const TH1D& sourceHistogram) const {
//
//  TProfile* resultProfile = nullptr;
//  resultProfile = new TProfile(Form("%s_profile",sourceHistogram.GetName()),
//                               sourceHistogram.GetTitle(),
//                               sourceHistogram.GetNbinsX(),
//                               sourceHistogram.GetXaxis()->GetXmin(),
//                               sourceHistogram.GetXaxis()->GetXmax());
//
//  for (int i=1; i<=sourceHistogram.GetNbinsX(); ++i) {
//    resultProfile->Fill(sourceHistogram.GetBinCenter(i), sourceHistogram.GetBinContent(i));
//  }
//  resultProfile->SetLineColor(sourceHistogram.GetLineColor());
//  resultProfile->SetLineWidth(sourceHistogram.GetLineWidth());
//  resultProfile->SetLineStyle(sourceHistogram.GetLineStyle());
//  resultProfile->SetFillColor(sourceHistogram.GetFillColor());
//  resultProfile->SetFillStyle(sourceHistogram.GetFillStyle());
//
//  return resultProfile;
//}

//
// Helper method calculating average value in the given cut regions
//
std::vector<double> AnalyzerMatBudget::average(TGraph& myGraph, std::vector<double> cuts) {

  std::vector<double> averages;
  if (cuts.size()<2) return averages;

  std::sort(cuts.begin(), cuts.end());
  int iBorder;
  int nBoarders=cuts.size();
  double valuesCount, valuesSum;
  double vx, vy;

  for (iBorder=0; iBorder<nBoarders-1; ++iBorder) {

    // Here we have a cut between
    // cuts[iBorder] and cuts[iBorder+1]
    //std::cerr << myGraph.GetTitle() << std::endl;
    //std::cerr << cuts[iBorder] << "< x <= "<< cuts[iBorder+1] << std::endl;

    // Average on points within the cut
    valuesSum=0;
    valuesCount=0;
    for (int iPoint=0; iPoint<myGraph.GetN(); ++iPoint) {
      myGraph.GetPoint(iPoint, vx, vy);
      if ((vx>=cuts[iBorder]) && (vx<cuts[iBorder+1])) {
        valuesCount++;
        valuesSum+=vy;
      }
    }
    averages.push_back(valuesSum/valuesCount);
  }
  return averages;
}

//
// Modifies a TGraph, so that it looks like a histogram (can be filled)
//
void AnalyzerMatBudget::closeGraph(TGraph& myGraph) {
  double x, y, x0, y0;
  myGraph.GetPoint(myGraph.GetN()-1, x, y);
  myGraph.GetPoint(0, x0, y0);
  myGraph.SetPoint(myGraph.GetN(), x,0);
  myGraph.SetPoint(myGraph.GetN(), x0,0);
}

//
// Helper method calculating average histogram values in histogram range from lowest bin to cutoff
//
double AnalyzerMatBudget::averageHistogramValues(const TH1D& his, double cutoff) const {

  double avg   = 0.0;
  int    cobin = 1;

  // Find last relevant bin
  while ((cobin < his.GetNbinsX()) && (his.GetBinLowEdge(cobin) < cutoff)) cobin++;

  // Calculate average
  for (auto i = 1; i <= cobin; i++) avg = avg + his.GetBinContent(i) / (double)cobin;

  return avg;
}

//
// Helper method calculating average histogram values in range cutoffStart - cutoffEnd
//
double AnalyzerMatBudget::averageHistogramValues(const TH1D& his, double cutoffStart, double cutoffEnd) const {

  double avg        = 0.0;
  int    coBinStart = 1;
  int    coBinEnd   = 1;

  if (cutoffStart >= cutoffEnd) return 0;

  // Find first relevant bin
  while ((coBinStart < his.GetNbinsX()) && (his.GetBinLowEdge(coBinStart) < cutoffStart)) coBinStart++;
  coBinEnd=coBinStart;

  // Find last relevant bin
  while ((coBinEnd < his.GetNbinsX()) && (his.GetBinLowEdge(coBinEnd) < cutoffEnd)) coBinEnd++;
  double coBinN=coBinEnd-coBinStart+1; // TODO: IMPORTANT check this

  // Calculate average
  if (coBinStart> his.GetNbinsX() - 1) coBinStart= his.GetNbinsX() - 1;
  if (coBinEnd > his.GetNbinsX() - 1) coBinEnd= his.GetNbinsX() - 1;
  for (int i = coBinStart; i <= coBinEnd; i++) avg = avg + his.GetBinContent(i) / coBinN;

  return avg;
}

//
// Material budget visitor - constructor
//
MatBudgetVisitor::MatBudgetVisitor(Track& matTrack, std::map<std::string, Material>& matBudget, TH2D& radMap, TH2D& radMapCount, TH2D& intMap, TH2D& intMapCount) :
    m_matTrack(matTrack),
    m_nEntries(0),
    m_matBudget(matBudget),
    m_radMap(radMap),
    m_radMapCount(radMapCount),
    m_intMap(intMap),
    m_intMapCount(intMapCount)
{
}

//
// Destructor
//
MatBudgetVisitor::~MatBudgetVisitor()
{
}

//
// Visit BeamPipe -> update track with beam pipe hit
//
void MatBudgetVisitor::visit(const BeamPipe& bp)
{
  // Add hit corresponding with beam-pipe
  double theta    = m_matTrack.getTheta();
  double distance = (bp.radius()+bp.thickness())/2./sin(theta);
  HitPtr hit(new Hit(distance));
  hit->setOrientation(HitOrientation::Horizontal);
  hit->setObjectKind(HitKind::Inactive);

  Material material;
  material.radiation   = bp.radLength()/sin(theta);
  material.interaction = bp.intLength()/sin(theta);
  hit->setCorrectedMaterial(material);
  hit->setBeamPipe(true);
  m_matTrack.addHit(std::move(hit));
}

//
// Visit BarrelModule (no limits on Rods, Layers or Barrels)
//
void MatBudgetVisitor::visit(const BarrelModule& m)
{
  analyzeModuleMB(m);
}

//
// Visit EndcapModule (no limits on Rings or Endcaps)
//
void MatBudgetVisitor::visit(const EndcapModule& m)
{
  analyzeModuleMB(m);
}

//
// Post visit method called after the whole visit-accept pattern done to recalibrate histograms, etc.
//
void MatBudgetVisitor::postvisit()
{
}

//
// Analyze if module crossed by given track & how much material is in the way
//
void MatBudgetVisitor::analyzeModuleMB(const DetectorModule& m)
{
  // Collision detection: material tracks being shot in z+ only, so consider only modules that lie on +Z side
  if (m.maxZ() > 0) {

    XYZVector direction(m_matTrack.getDirection());

    auto pair    = m.checkTrackHits(m_matTrack.getOrigin(), direction);
    auto hitRho  = pair.first.rho();
    auto hitType = pair.second;

    if (hitType!=HitType::NONE) {

      // Number of entries
      m_nEntries++;

      Material material;
      material.radiation   = m.getModuleCap().getRadiationLength();
      material.interaction = m.getModuleCap().getInteractionLength();

      // Fill material map
      double theta = m_matTrack.getTheta();
      double rho   = hitRho;
      double z     = rho/tan(theta);

      if (material.radiation>0){

        m_radMap.Fill(z,rho,material.radiation);
        m_radMapCount.Fill(z,rho);
      }
      if (material.interaction>0) {

        m_intMap.Fill(z,rho,material.interaction);
        m_intMapCount.Fill(z,rho);
      }

      // Treat barrel & endcap modules separately
      double tiltAngle = m.tiltAngle();

      if (m.subdet() == BARREL) {

        material.radiation   /= sin(theta + tiltAngle);
        material.interaction /= sin(theta + tiltAngle);

        // Fill barrel container
        m_matBudget["Barrel"].radiation   += material.radiation;
        m_matBudget["Barrel"].interaction += material.interaction;
      }
      else if (m.subdet() == ENDCAP) {

        material.radiation   /= cos(theta + tiltAngle - M_PI/2); // Endcap has tiltAngle = pi/2
        material.interaction /= cos(theta + tiltAngle - M_PI/2); // Endcap has tiltAngle = pi/2

        // Fill endcap container
        m_matBudget["Endcap"].radiation   += material.radiation;
        m_matBudget["Endcap"].interaction += material.interaction;
      }
      else {
        logWARNING("AnalyzerMatBudget::analyzeModuleMB -> incorrectly scaled material, unknown module type. Neither barrel or endcap");
      }

      // Fill other components
      std::map<std::string, Material> components = m.getModuleCap().getComponentsRI();

      for (auto iComponent : components) {

        m_matBudget[iComponent.first].radiation   += iComponent.second.radiation / (m.subdet()==BARREL ? sin(theta + tiltAngle) : cos(theta + tiltAngle - M_PI/2));
        m_matBudget[iComponent.first].interaction += iComponent.second.interaction / (m.subdet()==BARREL ? sin(theta + tiltAngle) : cos(theta + tiltAngle - M_PI/2));
      }

      // Create Hit object with appropriate parameters, add to Track t
      HitPtr hit(new Hit(pair.first.R(), &m, hitType));
      hit->setCorrectedMaterial(material);
      m_matTrack.addHit(std::move(hit));
    }
  } // Z>0
}
