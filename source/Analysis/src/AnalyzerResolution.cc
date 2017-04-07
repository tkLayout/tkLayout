/*
 * AnalyzerResolution.cc
 *
 *  Created on: 20. 4. 2016
 *      Author: drasal
 */
#include <AnalyzerResolution.h>

#include "BeamPipe.h"
#include "CsvTextBuilder.h"
#include "DetectorModule.h"
#include "MainConfigHandler.h"
#include "MaterialProperties.h"
#include "MessageLogger.h"
#include "ModuleCap.h"
#include <Palette.h>
#include "RootWContent.h"
#include "RootWImage.h"
#include "RootWPage.h"
#include "RootWSite.h"
#include "RootWTable.h"
#include "SimParms.h"
#include <TCanvas.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TList.h>
#include "Tracker.h"
#include "Track.h"
#include <TStyle.h>
#include "VisitorMatTrack.h"
#include "Units.h"

//
// AnalyzerResolution constructor
//
AnalyzerResolution::AnalyzerResolution(const Detector& detector) : AnalyzerUnit("AnalyzerResolution", detector),
 m_nTracks(0),
 m_etaMin(-1*SimParms::getInstance().getMaxEtaCoverage()),
 m_etaMax(+1*SimParms::getInstance().getMaxEtaCoverage()),
 c_nBins(SimParms::getInstance().getMaxEtaCoverage()/vis_eta_step)  // Default number of bins in histogram from eta=0  to max_eta_coverage)
{};

//
// AnalyzerResolution init method
//
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

    // Prepare Csv containers -> keep final pT/p resolution in csv file
    m_csvResPt = std::unique_ptr<CsvTextBuilder>(new CsvTextBuilder());
    m_csvResP  = std::unique_ptr<CsvTextBuilder>(new CsvTextBuilder());
    m_csvHitCol= std::unique_ptr<CsvTextBuilder>(new CsvTextBuilder());

    m_csvResPt->addCsvElement("Label", "Tag name");
    m_csvResPt->addCsvElement("Label", "Resolution type");
    m_csvResPt->addCsvElement("Label", "Transverse momentum [GeV]");

    m_csvResP->addCsvElement("Label", "Tag name");
    m_csvResP->addCsvElement("Label", "Resolution type");
    m_csvResP->addCsvElement("Label", "Total momentum [GeV]");

    m_csvHitCol->addCsvElement("Label", "Eta");
    m_csvHitCol->addCsvElement("Label", "BPHit");
    m_csvHitCol->addCsvElement("Label", "NxDetHits (s_R:s_Z:X/X0@Eta)");

    for (unsigned int iBorder=0; iBorder<SimParms::getInstance().getNEtaRegions()-1; ++iBorder) {

      ostringstream label("");
      label << "Real geom. eta(" << std::setiosflags(ios::fixed) << std::setprecision(1) << SimParms::getInstance().etaRegionRanges[iBorder] << "-" << SimParms::getInstance().etaRegionRanges[iBorder+1] << ")";
      m_csvResPt->addCsvElement("Label", label.str());
      m_csvResP->addCsvElement("Label", label.str());
    }

    m_csvResPt->addCsvEOL("Label");
    m_csvResP->addCsvEOL("Label");
    m_csvHitCol->addCsvEOL("Label");

    m_isInitOK = true;
    return m_isInitOK;
  }
}

//
// AnalyzerResolution visualization method
//
bool AnalyzerResolution::analyze()
{
  // Check that initialization OK
  if (!m_isInitOK) return false;

  // Random generator
  TRandom3 myDice;
  myDice.SetSeed(random_seed);

  // Initialize
  double efficiency  = SimParms::getInstance().efficiency();

  for (int iTrack = 0; iTrack < m_nTracks; iTrack++) {

    // Define track
    Track matTrack;

    double eta   = 0 + m_etaMax/m_nTracks*(iTrack+0.5);
    double theta = 2 * atan(exp(-eta));
    double phi   = myDice.Rndm() * M_PI * 2.0;
    double pT    = 100*Units::TeV; // Arbitrarily high number

    matTrack.setThetaPhiPt(theta, phi, pT);

    // Use uniform distribution to simulate position of primary interaction (rphi negligible). May be also triangular or gaussian (depends on accelerator design)!
    double zPos = 0.;
    if (SimParms::getInstance().useLumiRegInAnalysis()) zPos = (myDice.Rndm()*2 - 1)*SimParms::getInstance().zErrorIP();
    matTrack.setOrigin(0, 0, zPos);

    // Assign material to the track
    VisitorMatTrack matVisitor(matTrack);
    m_beamPipe->accept(matVisitor);                                 // Assign to material track hit corresponding to beam-pipe
    for (auto iTracker : m_trackers) iTracker->accept(matVisitor);  // Assign to material track hits corresponding to modules

    // Output hits to a file for debuggin purposes
    m_csvHitCol->addCsvElement(std::string("LastEta="+any2str(eta, 2)), any2str(eta, 2));

    std::map<double, std::string> hitInfoRSorted;
    for (std::vector<std::unique_ptr<Hit>>::const_iterator iter=matTrack.getBeginHits(); iter!=matTrack.getEndHits(); iter++) {

      std::string hitInfo = any2str((*iter)->getRPos()/Units::mm,1)+":"+any2str((*iter)->getZPos()/Units::mm,1)+":"+any2str((*iter)->getCorrectedMaterial().radiation, 6);
      hitInfoRSorted[(*iter)->getRPos()] = hitInfo;
    }
    for (auto& iMap : hitInfoRSorted) m_csvHitCol->addCsvElement(std::string("LastEta="+any2str(eta, 2)), iMap.second);
    hitInfoRSorted.clear();
    m_csvHitCol->addCsvEOL(std::string("LastEta="+any2str(eta, 2)));

    // Go through tracks with non-zero number of hits
    if (!matTrack.hasNoHits()) {

      for (string tag : matTrack.getTags()) {

        // Analyze only hits within track coming from the defined bunch of detectors (the same tag)
        matTrack.keepTaggedHitsOnly(tag);

        // Print hits
        //matTrack.printHits();

        // Remove some hits randomly based on inefficiency parameter
        //if (efficiency!=1) matTrack.addEfficiency(pxlEfficiency, stripEfficiency);

        // For each momentum/transverse momentum compute the tracks error
        for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

          int    parameter = pIter/Units::MeV; // Store p or pT in MeV as int (key to the map)
          double momentum  = pIter;

          // Case I) Initial momentum is equal to pT
          double pT = momentum;

          // Active+passive material
          TrackPtr trackPt(new Track(matTrack));
          trackPt->resetPt(pT);

          // Remove tracks with less than 3 hits
          if (trackPt->getNActiveHits(tag, true)>2) {

            std::map<int, TrackCollection>& myMap        = m_taggedTrackPtCollectionMap[tag];
            TrackCollection&                myCollection = myMap[parameter];
            myCollection.push_back(std::move(trackPt));
          }

          // Ideal (no material)
          TrackPtr idealTrackPt(new Track(matTrack));
          idealTrackPt->resetPt(pT);

          // Remove tracks with less than 3 hits & remove material
          idealTrackPt->removeMaterial();

          if (idealTrackPt->getNActiveHits(tag, true)>2) {

            std::map<int, TrackCollection>& myMapIdeal        = m_taggedTrackPtCollectionMapIdeal[tag];
            TrackCollection&                myCollectionIdeal = myMapIdeal[parameter];
            myCollectionIdeal.push_back(std::move(idealTrackPt));
          }

          // Case II) Initial momentum is equal to p
          pT = momentum*sin(theta);

          // Active+passive material
          TrackPtr trackP(new Track(matTrack));
          trackP->resetPt(pT);

          // Remove tracks with less than 3 hits
          if (trackP->getNActiveHits(tag, true)>2) {

            std::map<int, TrackCollection>& myMapII        = m_taggedTrackPCollectionMap[tag];
            TrackCollection&                myCollectionII = myMapII[parameter];
            myCollectionII.push_back(std::move(trackP));
          }

          // Ideal (no material)
          TrackPtr idealTrackP(new Track(matTrack));
          idealTrackP->resetPt(pT);

          // Remove tracks with less than 3 hits & remove material
          idealTrackP->removeMaterial();

          if (idealTrackP->getNActiveHits(tag, true)>2) {

            std::map<int, TrackCollection>& myMapIdealII        = m_taggedTrackPCollectionMapIdeal[tag];
            TrackCollection&                myCollectionIdealII = myMapIdealII[parameter];
            myCollectionIdealII.push_back(std::move(idealTrackP));
          }
        } // For momenta
      } // For tags
    }
  } // For tracks

  m_isAnalysisOK = true;
  return m_isAnalysisOK;
}

//
// AnalyzerResolution visualization method
//
bool AnalyzerResolution::visualize(RootWSite& webSite)
{
  // Check that initialization & analysis OK
  if (!m_isInitOK && !m_isAnalysisOK) return false;

  // Set Rainbow palette for drawing
  Palette::setRootPalette();

  // Go through all trackers & prepare web content
  for (auto& key : m_taggedTrackPtCollectionMap) {

    std::string tag       = key.first;
    std::string pageTitle = "Resolution";
    std::string additionalSummaryTag;
    double verticalScale  = 1;

    // Correct naming...
    std::string wName = "";
    if      (tag=="beampipe")                wName = "BP";
    else if (tag=="pixel")                   wName = "Pixel";
    else if (tag=="trigger" || tag=="strip") wName = "Strip";
    else if (tag=="barrel")                  wName = "BRL";
    else if (tag=="endcap")                  wName = "ECAP";
    else if (tag=="forward")                 wName = "FWD";
    else if (tag=="tracker")                 wName = "Tracker";
    else                                     wName = tag;

    pageTitle              += " ("+wName+")";
    additionalSummaryTag    = "_"+wName+"_";
    verticalScale           = 10;
    std::string pageAddress = "indexResol" + wName + ".html";

    int webPriority         = 0;
    if      (wName=="PIXEL")  webPriority = web_priority_Resol;
    else if (wName=="STRIP")  webPriority = web_priority_Resol-1;
    else if (wName=="FWD")    webPriority = web_priority_Resol-2;
    else if (wName=="BARREL") webPriority = web_priority_Resol-3;
    else if (wName=="ENDCAP") webPriority = web_priority_Resol-4;
    else if (wName=="TRK")    webPriority = web_priority_Resol-5;

    RootWPage& myPage = webSite.addPage(pageTitle, webPriority);
    myPage.setAddress(pageAddress);

    // Canvases
    gStyle->SetGridStyle(style_grid);
    gStyle->SetGridColor(Palette::color_hard_grid);
    gStyle->SetOptStat(0);

    // 4 scenarios - const pt with/without material + const p with/without material
    for (int i=0; i<4; i++) {

      const std::map<int,TrackCollection>* taggedTrackCollectionMap = nullptr;
      std::string scenario, scenarioName, webScenarioName;

      if (i==0) {
        taggedTrackCollectionMap = &m_taggedTrackPtCollectionMap[tag];
        scenario                 = "withMS_Pt";
        scenarioName             = "const p_{T}";
        webScenarioName          = "const pt";
      }
      if (i==1) {
        taggedTrackCollectionMap = &m_taggedTrackPtCollectionMapIdeal[tag];
        scenario                 = "noMS_Pt";
        scenarioName             = "const p_{T}";
        webScenarioName          = "const pt";
      }
      if (i==2) {
        taggedTrackCollectionMap = &m_taggedTrackPCollectionMap[tag];
        scenario                 = "withMS_P";
        scenarioName             = "const p";
        webScenarioName          = scenarioName;
      }
      if (i==3) {
        taggedTrackCollectionMap = &m_taggedTrackPCollectionMapIdeal[tag];
        scenario                 = "noMS_P";
        scenarioName             = "const p";
        webScenarioName          = scenarioName;
      }

      // Histogram arrays[momenta] - array of profile histograms for different momenta
      std::vector<unique_ptr<TProfile>> profHisArray_Pt;
      std::vector<unique_ptr<TProfile>> profHisArray_P;
      std::vector<unique_ptr<TProfile>> profHisArray_D0;
      std::vector<unique_ptr<TProfile>> profHisArray_Z0;
      std::vector<unique_ptr<TProfile>> profHisArray_Phi0;
      std::vector<unique_ptr<TProfile>> profHisArray_CotgTheta;
      std::vector<unique_ptr<TProfile>> profHisArray_CTau;

      preparePlot(profHisArray_Pt       , "pT"       , scenarioName, *taggedTrackCollectionMap);
      preparePlot(profHisArray_P        , "p"        , scenarioName, *taggedTrackCollectionMap);
      preparePlot(profHisArray_D0       , "d0"       , scenarioName, *taggedTrackCollectionMap);
      preparePlot(profHisArray_Z0       , "z0"       , scenarioName, *taggedTrackCollectionMap);
      preparePlot(profHisArray_Phi0     , "phi0"     , scenarioName, *taggedTrackCollectionMap);
      preparePlot(profHisArray_CotgTheta, "cotgTheta", scenarioName, *taggedTrackCollectionMap);
      preparePlot(profHisArray_CTau     , "ctau"     , scenarioName, *taggedTrackCollectionMap);

      std::string contentName = "";
      bool        contentVis  = true;

      if      (i==0) {
        contentName = "Track resolution for const Pt across "+web_etaLetter+" (active+pasive material)";
        contentVis  = true;
      }
      else if (i==1) {
        contentName = "Track resolution for const Pt across "+web_etaLetter+" (ideal - no material)";
        contentVis  = true;
      }
      else if (i==2) {
        contentName = "Track resolution for const P across "+web_etaLetter+" (active+pasive material)";
        contentVis  = false;
      }
      else if (i==3) {
        contentName = "Track resolution for const P across "+web_etaLetter+" (ideal - no material)";
        contentVis  = false;
      }
      else {

        std::ostringstream message;
        message << "AnalyzerResolution::visualize(): Badly implemented algorithm, should never get here, check!!!";
        logERROR(message.str());
        return false;
      }
      RootWContent& myContentPlots = myPage.addContent(contentName, contentVis);

      // a) Resolution in Pt
      TCanvas canvasResPtLin(std::string("ResPtLin_"+scenario+"_"+tag).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
      canvasResPtLin.SetGrid(1,1);
      canvasResPtLin.SetLogy(0);
      canvasResPtLin.SetFillColor(Palette::color_plot_background);
      canvasResPtLin.SetObjectStat(false);
      for (auto itHis=profHisArray_Pt.begin(); itHis!=profHisArray_Pt.end(); itHis++) {

        if (itHis==profHisArray_Pt.begin()) (*itHis)->Draw("PE1");
        else                                (*itHis)->Draw("PE1 SAME");
      }

      RootWImage& myImageLinPt = myContentPlots.addImage(canvasResPtLin, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImageLinPt.setComment("Transverse momentum resolution vs. "+web_etaLetter+" (lin. scale) - " + webScenarioName.c_str() + " across "+web_etaLetter);
      myImageLinPt.setName(Form("linptres_%s_%s", tag.c_str(), scenario.c_str()));

      TCanvas canvasResPtLog(std::string("ResPtLog_"+scenario+"_"+tag).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
      canvasResPtLog.SetGrid(1,1);
      canvasResPtLog.SetLogy(1);
      canvasResPtLog.SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_Pt.begin(); itHis!=profHisArray_Pt.end(); itHis++) {

        if (itHis==profHisArray_Pt.begin()) (*itHis)->Draw("PE1");
        else                                (*itHis)->Draw("PE1 SAME");
      }
      RootWImage& myImageLogPt = myContentPlots.addImage(canvasResPtLog, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImageLogPt.setComment("Transverse momentum resolution vs. "+web_etaLetter+" (log. scale) - " + webScenarioName.c_str() + " across "+web_etaLetter);
      myImageLogPt.setName(Form("logptres_%s_%s", tag.c_str(), scenario.c_str()));

      profHisArray_Pt.clear();

      // b) Resolution in P
      TCanvas canvasResP(std::string("ResP_"+scenario+"_"+tag).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
      canvasResP.SetGrid(1,1);
      canvasResP.SetLogy(1);
      canvasResP.SetFillColor(Palette::color_plot_background);

      for (auto itHis=profHisArray_P.begin(); itHis!=profHisArray_P.end(); itHis++) {

        if (itHis==profHisArray_P.begin()) (*itHis)->Draw("PE1");
        else                               (*itHis)->Draw("PE1 SAME");
      }

      RootWImage& myImageP = myContentPlots.addImage(canvasResP, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImageP.setComment("Momentum resolution vs. "+web_etaLetter+" - " + webScenarioName.c_str() + " across "+web_etaLetter);
      myImageP.setName(Form("pres_%s_%s", tag.c_str(), scenario.c_str()));

      profHisArray_P.clear();

      // c) Resolution in D0
      TCanvas canvasResD0(std::string("ResD0_"+scenario+"_"+tag).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
      canvasResD0.SetGrid(1,1);
      canvasResD0.SetLogy(1);
      canvasResD0.SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_D0.begin(); itHis!=profHisArray_D0.end(); itHis++) {

        if (itHis==profHisArray_D0.begin()) (*itHis)->Draw("PE1");
        else                                (*itHis)->Draw("PE1 SAME");
      }

      RootWImage& myImageD0 = myContentPlots.addImage(canvasResD0, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImageD0.setComment("d0 resolution vs. "+web_etaLetter+" - " + webScenarioName.c_str() + " across "+web_etaLetter);
      myImageD0.setName(Form("d0res_%s_%s", tag.c_str(), scenario.c_str()));

      profHisArray_D0.clear();

      // d) Resolution in Z0
      TCanvas canvasResZ0(std::string("ResZ0_"+scenario+"_"+tag).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
      canvasResZ0.SetGrid(1,1);
      canvasResZ0.SetLogy(1);
      canvasResZ0.SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_Z0.begin(); itHis!=profHisArray_Z0.end(); itHis++) {

        if (itHis==profHisArray_Z0.begin()) (*itHis)->Draw("PE1");
        else                                (*itHis)->Draw("PE1 SAME");
      }

      RootWImage& myImageZ0 = myContentPlots.addImage(canvasResZ0, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImageZ0.setComment("z0 resolution vs. "+web_etaLetter+" - " + webScenarioName.c_str() + " across "+web_etaLetter);
      myImageZ0.setName(Form("z0res_%s_%s", tag.c_str(), scenario.c_str()));

      profHisArray_Z0.clear();

      // e) Resolution in Phi0
      TCanvas canvasResPhi0(std::string("ResPhi0_"+scenario+"_"+tag).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
      canvasResPhi0.SetGrid(1,1);
      canvasResPhi0.SetLogy(1);
      canvasResPhi0.SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_Phi0.begin(); itHis!=profHisArray_Phi0.end(); itHis++) {

        if (itHis==profHisArray_Phi0.begin()) (*itHis)->Draw("PE1");
        else                                  (*itHis)->Draw("PE1 SAME");
      }

      RootWImage& myImagePhi0 = myContentPlots.addImage(canvasResPhi0, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePhi0.setComment(web_phiLetter + "0 resolution vs. "+web_etaLetter+" - " + webScenarioName.c_str() + " across "+web_etaLetter);
      myImagePhi0.setName(Form("phi0res_%s_%s", tag.c_str(), scenario.c_str()));

      profHisArray_Phi0.clear();

      // f) Resolution in cotg(theta)
      TCanvas canvasResCotgTh(std::string("ResCotgTh_"+scenario+"_"+tag).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
      canvasResCotgTh.SetGrid(1,1);
      canvasResCotgTh.SetLogy(1);
      canvasResCotgTh.SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_CotgTheta.begin(); itHis!=profHisArray_CotgTheta.end(); itHis++) {

        if (itHis==profHisArray_CotgTheta.begin()) (*itHis)->Draw("PE1");
        else                                       (*itHis)->Draw("PE1 SAME");
      }

      RootWImage& myImageCtg = myContentPlots.addImage(canvasResCotgTh, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImageCtg.setComment("Ctg("+web_thetaLetter+") resolution vs. "+web_etaLetter+" - " + webScenarioName.c_str() + " across "+web_etaLetter);
      myImageCtg.setName(Form("cotgThres_%s_%s", tag.c_str(), scenario.c_str()));

      profHisArray_CotgTheta.clear();

      // g) Resolution in ctau
      TCanvas canvasResCTau(std::string("ResCTau_"+scenario+"_"+tag).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
      canvasResCTau.SetGrid(1,1);
      canvasResCTau.SetLogy(1);
      canvasResCTau.SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_CTau.begin(); itHis!=profHisArray_CTau.end(); itHis++) {

        if (itHis==profHisArray_CTau.begin()) (*itHis)->Draw("PE1");
        else                                  (*itHis)->Draw("PE1 SAME");
      }

      RootWImage& myImageCTau = myContentPlots.addImage(canvasResCTau, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImageCTau.setComment("c"+web_tauLetter+" resolution vs. "+web_etaLetter+" - " + webScenarioName.c_str() + " across "+web_etaLetter);
      myImageCTau.setName(Form("ctauRes_%s_%s", tag.c_str(), scenario.c_str()));

      profHisArray_CTau.clear();

    } // For scenarios

    // Set Summary content for const Pt
    RootWContent& summaryContent_Pt = myPage.addContent("Summary - const Pt across "+web_etaLetter);
    prepareSummaryTable(tag, "Pt", myPage, summaryContent_Pt, *m_csvResPt);

    // Set Summary content for const P
    RootWContent& summaryContent_P = myPage.addContent("Summary - const P across "+web_etaLetter, false);
    prepareSummaryTable(tag, "P", myPage, summaryContent_P, *m_csvResP);
  } // For tag

  m_isVisOK = true;
  return m_isVisOK;
}

//
// Get Csv text output for const pt -> exception thrown if doesn't exist
//
const CsvTextBuilder& AnalyzerResolution::getCsvResPt() const
{
  // Check if exist
  if (m_csvResPt) return *m_csvResPt;

  // Otherwise throw exception
  else throw std::invalid_argument( "CsvTextBuilder::getCsvResPt() - csvResPt not defined (null reference), check!!!" );
}

//
// Get Csv text output for const p -> exception thrown if doesn't exist
//
const CsvTextBuilder& AnalyzerResolution::getCsvResP() const
{
  // Check if exist
  if (m_csvResP) return *m_csvResP;

  // Otherwise throw exception
  else throw std::invalid_argument( "CsvTextBuilder::getCsvResP() - csvResP not defined (null reference), check!!!" );
}

//
// Get Csv text output for hit collection-> exception thrown if doesn't exist
//
const CsvTextBuilder& AnalyzerResolution::getCsvHitCol() const
{
  // Check if exist
  if (m_csvHitCol) return *m_csvHitCol;

  // Otherwise throw exception
  else throw std::invalid_argument( "CsvTextBuilder::getCsvHitCol() - csvHitCol not defined (null reference), check!!!" );
}

//
// Prepare plot: fill with data & set properties; varType specifies variable type to be filled: pT, p, d0, z0, phi0, cotgTheta
//
void AnalyzerResolution::preparePlot(std::vector<unique_ptr<TProfile>>& profHisArray, std::string varType, std::string scenario, const std::map<int, TrackCollection>& mapCollection)
{
  // Momentum counter for color setting
  unsigned int iMomentum =0;

  // Histogram the data - for different momenta
  for (auto& col : mapCollection) {

    // Prepare histogram for given momentum
    int momentum = col.first;

    std::string title("");
    std::string name("");
    if (varType=="pT") {
      name  = "pT_vs_eta"+any2str(momentum/Units::GeV);
      title = "p_{T} resolution versus #eta - "+scenario+" across #eta;#eta;#delta p_{T}/p_{T} [%]";
    }
    else if (varType=="p") {
      name  = "p_vs_eta"+any2str(momentum/Units::GeV);
      title = "p resolution versus #eta - "+scenario+" across #eta;#eta;#delta p/p [%]";
    }
    else if (varType=="d0") {
      name  = "d0_vs_eta"+any2str(momentum/Units::GeV);
      title = "Transverse impact parameter error - "+scenario+" across #eta;#eta;#delta d_{0} [#mum]";
    }
    else if (varType=="z0") {
      name  = "z0_vs_eta"+any2str(momentum/Units::GeV);
      title = "Longitudinal impact parameter error - "+scenario+" across #eta;#eta;#delta z_{0} [#mum]";
    }
    else if (varType=="phi0") {
      name  = "phi0_vs_eta"+any2str(momentum/Units::GeV);
      title = "Track azimuthal angle error - "+scenario+" across #eta;#eta;#delta #phi [deg]";
    }
    else if (varType=="cotgTheta") {
      name  = "cotgTh_vs_eta"+any2str(momentum/Units::GeV);
      title = "Track polar angle error - "+scenario+" across #eta;#eta;#delta ctg(#theta)";
    }
    else if (varType=="ctau") {
      name  = "ctau_vs_eta"+any2str(momentum/Units::GeV);
      title = "c#tau resolution - "+scenario+" across #eta;#eta;#delta c#tau [#mum]";
    }
    std::unique_ptr<TProfile> profHis(new TProfile(name.c_str(), title.c_str(), c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));

    // Set style
    profHis->SetLineWidth(2.);
    profHis->SetMarkerStyle(21);
    profHis->SetMarkerSize(1.);
    profHis->SetLineColor(Palette::colorMomenta(iMomentum));
    profHis->SetMarkerColor(Palette::colorMomenta(iMomentum));

    for (auto& track : col.second) {

      double xVal = track->getEta();
      double yVal = 0;

      // Evaluate at reference point [r,z] =[0,0]
      double rPos = 0.0;
      if (varType=="pT")        yVal = track->getDeltaPtOverPt(rPos)*100; // In percent
      if (varType=="p")         yVal = track->getDeltaPOverP(rPos)*100;   // In percent
      if (varType=="d0")        yVal = track->getDeltaD0()/Units::um;
      if (varType=="z0")        yVal = track->getDeltaZ0()/Units::um;
      if (varType=="phi0")      yVal = track->getDeltaPhi0()/M_PI*180.; // In degerees
      if (varType=="cotgTheta") yVal = track->getDeltaCtgTheta(rPos);
      if (varType=="ctau")      yVal = track->getDeltaCTau()/Units::um;

      if (yVal>0) profHis->Fill(xVal, yVal);

    } // For tracks

    // Set axis range
    if (varType=="pT") {
      profHis->SetMaximum(c_max_dPtOverPt);
      profHis->SetMinimum(c_min_dPtOverPt);
    }
    if (varType=="p") {
      profHis->SetMaximum(c_max_dPtOverPt);
      profHis->SetMinimum(c_min_dPtOverPt);
    }
    if (varType=="d0") {
      profHis->SetMaximum(c_max_dD0);
      profHis->SetMinimum(c_min_dD0);
    }
    if (varType=="z0") {
      profHis->SetMaximum(c_max_dZ0);
      profHis->SetMinimum(c_min_dZ0);
    }
    if (varType=="phi0") {
      profHis->SetMaximum(c_max_dPhi0);
      profHis->SetMinimum(c_min_dPhi0);
    }
    if (varType=="cotgTheta") {
      profHis->SetMaximum(c_max_dCtgTheta);
      profHis->SetMinimum(c_min_dCtgTheta);
    }
    if (varType=="ctau") {
      profHis->SetMaximum(c_max_dCTau);
      profHis->SetMinimum(c_min_dCTau);
    }

    // Add histogram for given momentum to array
    profHisArray.push_back(std::move(profHis));

    // Increase counter
    iMomentum++;

  } // For momenta

}

//
// Prepare summary content table
//
void AnalyzerResolution::prepareSummaryTable(std::string tag, std::string scenario, RootWPage& webPage, RootWContent& summaryContent, CsvTextBuilder& csvContainer)
{
  RootWTable&   cutsSummaryTable  = summaryContent.addTable();
  RootWTable&   momSummaryTable   = summaryContent.addTable();

  std::map<std::string, RootWTable*> tableMap; // individual plot tables

  std::vector<std::string> plotNames;
  plotNames.push_back(web_deltaLetter+"pt/pt [%]:           ");
  plotNames.push_back(web_deltaLetter+"p/p [%]:             ");
  plotNames.push_back(web_deltaLetter+"d0 ["+web_muLetter+"m]:  ");
  plotNames.push_back(web_deltaLetter+"z0 ["+web_muLetter+"m]:  ");
  plotNames.push_back(web_deltaLetter+web_phiLetter+"0:          ");
  plotNames.push_back(web_deltaLetter+"ctg("+web_thetaLetter+"):");

  for (auto& name : plotNames) {

    RootWTable* myTable = &(summaryContent.addTable());
    tableMap[name] = myTable;
    tableMap[name]->setContent(0,0,name);
  }

  // Prepare the cuts for the averages
  ostringstream label;

  // Table explaining the cuts
  cutsSummaryTable.setContent(0,0,"Region: ");
  cutsSummaryTable.setContent(1,0,"Min "+web_etaLetter+":");
  cutsSummaryTable.setContent(2,0,"Max "+web_etaLetter+":");
  for (unsigned int iBorder=0; iBorder<SimParms::getInstance().getNEtaRegions()-1; ++iBorder) {

    cutsSummaryTable.setContent(0,iBorder+1,SimParms::getInstance().etaRegionNames[iBorder+1]);
    label.str(""); label << SimParms::getInstance().etaRegionRanges[iBorder];
    cutsSummaryTable.setContent(1,iBorder+1,label.str());
    label.str(""); label << SimParms::getInstance().etaRegionRanges[iBorder+1];
    cutsSummaryTable.setContent(2,iBorder+1,label.str());
  }

  // Table explaining momenta
  std::vector<double> momenta;
  for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) momenta.push_back(pIter);
  std::sort(momenta.begin(), momenta.end());

  momSummaryTable.setContent(0,0,"Particle momenta in GeV:  " );
  for (unsigned int iMom=0; iMom<momenta.size(); ++iMom) {

    label.str("");
    std::string color = Palette::colorMomentaNames(iMom);

    if (iMom!=momenta.size()-1) label << momenta[iMom]/Units::GeV << " (" << color << "),";
    else                        label << momenta[iMom]/Units::GeV << " (" << color << ").";
    momSummaryTable.setContent(0,iMom+1,label.str());
    momSummaryTable.setColor(0,iMom+1,Palette::colorMomenta(iMom));
  }

  // Cycle over the different measurements
  for (auto plotName : plotNames) {

    RootWTable& myTable = *(tableMap[plotName]);

    // Fill the table with the values, first the heading of momentum
    int baseColumn;
    int myColor = kBlack;

    // Set csv content
    std::string csvPlotName = plotName;
    std::string delta       = "&delta;";
    std::string ampersand   = "&";
    std::string semicolon   = ";";
    if (csvPlotName.find(delta)!= std::string::npos)    csvPlotName.replace(csvPlotName.find(delta), delta.length()    , "d");
    while (csvPlotName.find(ampersand)!=std::string::npos) csvPlotName.erase(csvPlotName.find(ampersand), ampersand.length());
    while (csvPlotName.find(semicolon)!=std::string::npos) csvPlotName.erase(csvPlotName.find(semicolon), semicolon.length());

    // Get values for all analyzed momenta
    for (unsigned int iMom=0; iMom<momenta.size(); ++iMom) {

      // Plot average values for different eta regions
      std::vector<double> averagesReal;
      std::vector<double> averagesIdeal;

      // Set table - 1.row & 1.column
      baseColumn = (SimParms::getInstance().getNEtaRegions()-1)*iMom + 1;
      myTable.setContent(0, baseColumn, momenta[iMom]/Units::GeV,0);
      myTable.setColor(0, baseColumn, Palette::colorMomenta(iMom));
      myTable.setContent(2, 0, "Real:      ");
      myTable.setContent(3, 0, "Ideal:     ");
      myTable.setContent(4, 0, "Real/Ideal:");

      // Analyze the canvas corresponding to given plotName -> print out average values
      std::string scenarioNameReal = "Track resolution for const "+scenario+" across "+web_etaLetter+" (active+pasive material)";
      std::string scenarioNameIdeal= "Track resolution for const "+scenario+" across "+web_etaLetter+" (ideal - no material)";

      const TCanvas* canvasReal = nullptr;
      const TCanvas* canvasIdeal= nullptr;

      std::string profileName = "";

      if      (plotName==std::string(web_deltaLetter+"pt/pt [%]:           ")) {
        canvasReal = &(webPage.findContent(scenarioNameReal).findImage("linptres_"+tag+"_withMS_"+scenario));
        canvasIdeal= &(webPage.findContent(scenarioNameIdeal).findImage("linptres_"+tag+"_noMS_"+scenario));
        profileName= "pT_vs_eta"+any2str(momenta[iMom]/Units::GeV,0);
      }
      else if (plotName==std::string(web_deltaLetter+"p/p [%]:             ")) {
        canvasReal = &(webPage.findContent(scenarioNameReal).findImage("pres_"+tag+"_withMS_"+scenario));
        canvasIdeal= &(webPage.findContent(scenarioNameIdeal).findImage("pres_"+tag+"_noMS_"+scenario));
        profileName= "p_vs_eta"+any2str(momenta[iMom]/Units::GeV,0);
      }
      else if (plotName==std::string(web_deltaLetter+"d0 ["+web_muLetter+"m]:  ")) {
        canvasReal = &(webPage.findContent(scenarioNameReal).findImage("d0res_"+tag+"_withMS_"+scenario));
        canvasIdeal= &(webPage.findContent(scenarioNameIdeal).findImage("d0res_"+tag+"_noMS_"+scenario));
        profileName= "d0_vs_eta"+any2str(momenta[iMom]/Units::GeV,0);
      }
      else if (plotName==std::string(web_deltaLetter+"z0 ["+web_muLetter+"m]:  ")) {
        canvasReal = &(webPage.findContent(scenarioNameReal).findImage("z0res_"+tag+"_withMS_"+scenario));
        canvasIdeal= &(webPage.findContent(scenarioNameIdeal).findImage("z0res_"+tag+"_noMS_"+scenario));
        profileName= "z0_vs_eta"+any2str(momenta[iMom]/Units::GeV,0);
      }
      else if (plotName==std::string(web_deltaLetter+web_phiLetter+"0:          ")) {
        canvasReal = &(webPage.findContent(scenarioNameReal).findImage("phi0res_"+tag+"_withMS_"+scenario));
        canvasIdeal= &(webPage.findContent(scenarioNameIdeal).findImage("phi0res_"+tag+"_noMS_"+scenario));
        profileName= "phi0_vs_eta"+any2str(momenta[iMom]/Units::GeV,0);
      }
      else if (plotName==std::string(web_deltaLetter+"ctg("+web_thetaLetter+"):")) {
        canvasReal = &(webPage.findContent(scenarioNameReal).findImage("cotgThres_"+tag+"_withMS_"+scenario));
        canvasIdeal= &(webPage.findContent(scenarioNameIdeal).findImage("cotgThres_"+tag+"_noMS_"+scenario));
        profileName= "cotgTh_vs_eta"+any2str(momenta[iMom]/Units::GeV,0);
      }

      // Prepare eta cuts
      std::vector<double> etaCuts;
      for (auto& iCut : SimParms::getInstance().etaRegionRanges) {
        etaCuts.push_back(iCut);
      }

      // Real geometry (active+passive)
      if (canvasReal!=nullptr) for (int i=0; i<canvasReal->GetListOfPrimitives()->GetSize(); ++i) {
        if (std::string(canvasReal->GetListOfPrimitives()->At(i)->ClassName())=="TProfile") {

          // TProfile
          const TProfile* myProfile = (const TProfile*)canvasReal->GetListOfPrimitives()->At(i);

          // If profile histogram for given momentum analyze
          if (std::string(myProfile->GetName())==profileName) {

            averagesReal = averageHisValues(*myProfile,etaCuts);
          }
        }
      }

      // Ideal geometry (no material)
      if (canvasIdeal!=nullptr) for (int i=0; i<canvasIdeal->GetListOfPrimitives()->GetSize(); ++i) {
        if (std::string(canvasIdeal->GetListOfPrimitives()->At(i)->ClassName())=="TProfile") {

          // TProfile
          const TProfile* myProfile = (const TProfile*)canvasIdeal->GetListOfPrimitives()->At(i);

          // If profile histogram for given momentum analyze
          if (std::string(myProfile->GetName())==profileName.c_str()) {

            averagesIdeal = averageHisValues(*myProfile,etaCuts);
          }
        }
      }

      // Fill resolution for different eta regions to a table
      for (unsigned int j=0; j<(SimParms::getInstance().getNEtaRegions()-1); ++j) {

        myTable.setContent(1, baseColumn+j, SimParms::getInstance().etaRegionNames[j+1]);
        myTable.setColor(1, baseColumn+j, myColor);
        if (averagesReal.size() > j) {

          myTable.setContent(2, baseColumn+j,averagesReal[j],2);
          myTable.setColor(2, baseColumn+j, myColor);
        }
        if (averagesIdeal.size() > j) {

          myTable.setContent(3, baseColumn+j,averagesIdeal[j],2);
          myTable.setColor(3, baseColumn+j, myColor);
        }
        if ((averagesReal.size() > j)&&(averagesIdeal.size() > j)) {
          myTable.setContent(4, baseColumn+j,averagesReal[j]/averagesIdeal[j],2);
          myTable.setColor(4, baseColumn+j, myColor);
        }
      }

      // Csv content
      if (plotName==*plotNames.begin() && iMom==0) csvContainer.addCsvElement(tag, tag);
      else                                         csvContainer.addCsvElement(tag, "");

      if (iMom==0) csvContainer.addCsvElement(tag, csvPlotName);
      else         csvContainer.addCsvElement(tag, "");

      csvContainer.addCsvElement(tag, any2str(momenta[iMom]/Units::GeV,0));

      for (const auto& average : averagesReal) csvContainer.addCsvElement(tag, average);

      csvContainer.addCsvEOL(tag);

    } // For momenta
  } // For plot names: dpt/pt, dp/p, ...

}

//
// Calculate average values in defined regions for given profile histogram
//
std::vector<double> AnalyzerResolution::averageHisValues(const TProfile& his, std::vector<double> regions)
{
  std::vector<double> averages;

  std::sort(regions.begin(), regions.end());
  for (int iBorder=0; iBorder<regions.size()-1; ++iBorder) {

    // Average on points within the cut
    double valuesSum  =0;
    int    valuesCount=0;
    for (int iBin=1; iBin<his.GetNbinsX(); ++iBin) {

      if ((his.GetBinLowEdge(iBin)>=regions[iBorder]) && ( (his.GetBinLowEdge(iBin)+his.GetBinWidth(iBin))<regions[iBorder+1])) {
        valuesCount++;
        valuesSum += his.GetBinContent(iBin);
      }
      if (his.GetBinLowEdge(iBin)+his.GetBinWidth(iBin)>regions[iBorder+1]) break;
    }
    if (valuesCount>0) averages.push_back(valuesSum/valuesCount);
    else               averages.push_back(0);
  }
  return averages;
}
