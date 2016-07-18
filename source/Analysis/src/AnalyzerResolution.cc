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
#include "rootweb.h"
#include "SimParms.h"
#include <TProfile.h>
#include <TRandom3.h>
#include "Tracker.h"
#include "Track.h"
#include <TStyle.h>
#include "Units.h"


AnalyzerResolution::AnalyzerResolution(std::vector<const Tracker*> trackers, const BeamPipe* beamPipe) : AnalyzerUnit("AnalyzerResolution", trackers, beamPipe),
 m_nTracks(0),
 m_etaMin(-1*geom_max_eta_coverage),
 m_etaMax(+1*geom_max_eta_coverage)
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

  // Random generator
  TRandom3 myDice;
  myDice.SetSeed(random_seed);

  // Initialize
  double efficiency  = SimParms::getInstance()->efficiency();

  // Tracks pruned
  bool isPruned = false;

  for (int iTrack = 0; iTrack < m_nTracks; iTrack++) {

    // Define track
    Track matTrack;

    double eta   = 0 + m_etaMax/m_nTracks*(iTrack+0.5);
    double theta = 2 * atan(exp(-eta));
    double phi   = myDice.Rndm() * M_PI * 2.0;
    double pT    = 100*Units::TeV; // Arbitrarily high number

    matTrack.setThetaPhiPt(theta, phi, pT);
    matTrack.setOrigin(0, 0, 0); // TODO: Not assuming z-error when calculating Material budget

    // Assign material to the track
    MatTrackVisitor matVisitor(matTrack);
    for (auto iTracker : m_trackers) iTracker->accept(matVisitor);  // Assign to material track hits corresponding to modules
    m_beamPipe->accept(matVisitor);                                 // Assign to material track hit corresponding to beam-pipe

    // Go through tracks with non-zero number of hits
    if (!matTrack.hasNoHits()) {

      for (string tag : matTrack.getTags()) {

        // Analyze only hits within track coming from the defined bunch of detectors (the same tag)
        matTrack.keepTaggedOnly(tag);

        // Add IP constraint
        if (SimParms::getInstance()->useIPConstraint()) matTrack.addIPConstraint(SimParms::getInstance()->rError(), SimParms::getInstance()->zErrorCollider());

        // Sort hits
        bool smallerRadius = true;
        matTrack.sortHits(true);
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
          trackPt->setThetaPhiPt(theta, phi, pT);

          // Remove tracks with less than 3 hits
          bool pruned = trackPt->pruneHits();
          if (pruned) isPruned = true;

          if (trackPt->getNActiveHits(tag, true)>2) {

            trackPt->computeErrors();
            std::map<int, TrackCollection>& myMap        = taggedTrackPtCollectionMap[tag];
            TrackCollection&                myCollection = myMap[parameter];
            myCollection.push_back(std::move(trackPt));
          }

          // Ideal (no material)
          TrackPtr idealTrackPt(new Track(matTrack));
          idealTrackPt->setThetaPhiPt(theta, phi, pT);

          // Remove tracks with less than 3 hits
          pruned = idealTrackPt->pruneHits();
          if (pruned) isPruned = true;

          // Remove material
          idealTrackPt->removeMaterial();
          if (idealTrackPt->getNActiveHits(tag, true)>2) {

            idealTrackPt->computeErrors();
            std::map<int, TrackCollection>& myMapIdeal        = taggedTrackPtCollectionMapIdeal[tag];
            TrackCollection&                myCollectionIdeal = myMapIdeal[parameter];
            myCollectionIdeal.push_back(std::move(idealTrackPt));
          }

          // Case II) Initial momentum is equal to p
          pT = momentum*sin(theta);

          // Active+passive material
          TrackPtr trackP(new Track(matTrack));
          trackP->setThetaPhiPt(theta, phi, pT);

          // Remove tracks with less than 3 hits
          pruned = trackP->pruneHits();
          if (pruned) isPruned = true;

          if (trackP->getNActiveHits(tag, true)>2) {

            trackP->computeErrors();
            std::map<int, TrackCollection>& myMapII        = taggedTrackPCollectionMap[tag];
            TrackCollection&                myCollectionII = myMapII[parameter];
            myCollectionII.push_back(std::move(trackP));
          }

          // Ideal (no material)
          TrackPtr idealTrackP(new Track(matTrack));
          idealTrackP->setThetaPhiPt(theta, phi, pT);

          // Remove tracks with less than 3 hits
          pruned = idealTrackP->pruneHits();
          if (pruned) isPruned = true;

          // Remove material
          idealTrackP->removeMaterial();

          if (idealTrackP->getNActiveHits(tag, true)>2) {

            idealTrackP->computeErrors();
            std::map<int, TrackCollection>& myMapIdealII        = taggedTrackPCollectionMapIdeal[tag];
            TrackCollection&                myCollectionIdealII = myMapIdealII[parameter];
            myCollectionIdealII.push_back(std::move(idealTrackP));
          }
        } // For momenta
      } // For tags
    }
  } // For tracks

  // Log pruning procedure -> no result mode might occur if starting momenta wrongly set
  if (isPruned) {

    std::string message = std::string("Resolution - some tracks pruned! Hits that don't follow the parabolic approximation removed. Check momenta if no result appears!");
    logWARNING(message);
  }

  return true;
}

bool AnalyzerResolution::visualize(RootWSite& webSite)
{
  // Check that initialization & analysis OK
  if (!m_isInitOK && !m_isAnalysisOK) return false;

  // Csv containers -> keep final pT/p resolution in csv file
  CsvTextBuilder csvResPt;
  CsvTextBuilder csvResP;

  csvResPt.addCsvElement("Label", "Tag name");
  csvResPt.addCsvElement("Label", "Resolution type");
  csvResPt.addCsvElement("Label", "Real/Ideal geometry");
  csvResPt.addCsvElement("Label", "Transverse momentum [GeV]");

  csvResP.addCsvElement("Label", "Tag name");
  csvResP.addCsvElement("Label", "Resolution type");
  csvResP.addCsvElement("Label", "Real/Ideal geometry");
  csvResP.addCsvElement("Label", "Total momentum [GeV]");

  for (unsigned int iBorder=0; iBorder<geom_name_eta_regions.size()-1; ++iBorder) {

    ostringstream label("");
    label << "eta(" << std::setiosflags(ios::fixed) << std::setprecision(1) << geom_range_eta_regions[iBorder] << "-" << geom_range_eta_regions[iBorder+1] << ")";
    csvResPt.addCsvElement("Label", label.str());
    csvResPt.addCsvElement("Label", label.str());
  }

  csvResPt.addCsvEOL("Label");
  csvResP.addCsvEOL("Label");

  // Set Rainbow palette for drawing
  Palette::setRootPalette(55);

  // Go through all trackers & prepare web content
  int webPriority         = web_priority_Resol;
  RootWPage*    myPage    = nullptr;
  RootWContent* myContent = nullptr;
  RootWImage*   myImage   = nullptr;
  RootWTable*   myTable   = nullptr;

  for (auto& key : taggedTrackPtCollectionMap) {

    std::string tag       = key.first;
    std::string pageTitle = "Resolution";
    std::string additionalSummaryTag;
    double verticalScale  = 1;

    // Correct naming...
    std::string wName = "";
    if      (tag=="beampipe")                wName = "BP";
    else if (tag=="pixel")                   wName = "PIXEL";
    else if (tag=="trigger" || tag=="strip") wName = "STRIP";
    else if (tag=="barrel")                  wName = "BRL";
    else if (tag=="endcap")                  wName = "ECAP";
    else if (tag=="forward")                 wName = "FWD";
    else if (tag=="tracker")                 wName = "TRK";
    else                                     wName = tag;

    pageTitle           += " ("+wName+")";
    additionalSummaryTag = "_"+wName+"_";
    verticalScale        = 10;
    std::string pageAddress = "errors" + wName + ".html";

    myPage = new RootWPage(pageTitle);
    myPage->setAddress(pageAddress);
    if      (wName=="PIXEL")  webSite.addPage(myPage,webPriority);
    else if (wName=="STRIP")  webSite.addPage(myPage,webPriority-1);
    else if (wName=="FWD")    webSite.addPage(myPage,webPriority-2);
    else if (wName=="BARREL") webSite.addPage(myPage,webPriority-3);
    else if (wName=="ENDCAP") webSite.addPage(myPage,webPriority-4);
    else if (wName=="TRK")    webSite.addPage(myPage,webPriority-5);
    else                      webSite.addPage(myPage);

    // Canvases
    std::unique_ptr<TCanvas> canvasResPtLog, canvasResPtLin, canvasResP, canvasResD0, canvasResZ0, canvasResPhi0, canvasResCotgTh;

    gStyle->SetGridStyle(style_grid);
    gStyle->SetGridColor(Palette::color_hard_grid);
    gStyle->SetOptStat(0);

    // 4 scenarios - const pt with/without material + const p with/without material
    for (int i=0; i<4; i++) {

      const std::map<int,TrackCollection>* taggedTrackCollectionMap = nullptr;
      std::string scenario, scenarioName;

      if (i==0) {
        taggedTrackCollectionMap = &taggedTrackPtCollectionMap[tag];
        scenario                 = "withMS_Pt";
        scenarioName             = "const P_{T}";
      }
      if (i==1) {
        taggedTrackCollectionMap = &taggedTrackPtCollectionMapIdeal[tag];
        scenario                 = "noMS_Pt";
        scenarioName             = "const P_{T}";
      }
      if (i==2) {
        taggedTrackCollectionMap = &taggedTrackPCollectionMap[tag];
        scenario                 = "withMS_P";
        scenarioName             = "const P";
      }
      if (i==3) {
        taggedTrackCollectionMap = &taggedTrackPCollectionMapIdeal[tag];
        scenario                 = "noMS_P";
        scenarioName             = "const P";
      }

      // Histogram arrays[momenta] - array of profile histograms for different momenta
      std::vector<unique_ptr<TProfile>> profHisArray_Pt;
      std::vector<unique_ptr<TProfile>> profHisArray_P;
      std::vector<unique_ptr<TProfile>> profHisArray_D0;
      std::vector<unique_ptr<TProfile>> profHisArray_Z0;
      std::vector<unique_ptr<TProfile>> profHisArray_Phi0;
      std::vector<unique_ptr<TProfile>> profHisArray_CotgTheta;

      preparePlot(profHisArray_Pt       , "pT"       , scenarioName, *taggedTrackCollectionMap);
      preparePlot(profHisArray_P        , "p"        , scenarioName, *taggedTrackCollectionMap);
      preparePlot(profHisArray_D0       , "d0"       , scenarioName, *taggedTrackCollectionMap);
      preparePlot(profHisArray_Z0       , "z0"       , scenarioName, *taggedTrackCollectionMap);
      preparePlot(profHisArray_Phi0     , "phi0"     , scenarioName, *taggedTrackCollectionMap);
      preparePlot(profHisArray_CotgTheta, "cotgTheta", scenarioName, *taggedTrackCollectionMap);

      if      (i==0) myContent = new RootWContent("Track resolution for const Pt across "  +web_etaLetter+" (active+pasive material)", true);
      else if (i==1) myContent = new RootWContent("Track resolution for const Pt across "  +web_etaLetter+" (ideal - no material)", true);
      else if (i==2) myContent = new RootWContent("Track resolution for const P across "   +web_etaLetter+" (active+pasive material)", false);
      else if (i==3) myContent = new RootWContent("Track resolution for const P across "   +web_etaLetter+" (ideal - no material)", false);
      else {

        std::ostringstream message;
        message << "AnalyzerResolution::visualize(): Badly implemented algorithm, should never get here, check!!!";
        logERROR(message.str());
        return false;
      }
      myPage->addContent(myContent);

      // a) Resolution in Pt
      canvasResPtLin = std::unique_ptr<TCanvas>(new TCanvas);
      canvasResPtLin->SetGrid(1,1);
      canvasResPtLin->SetLogy(0);
      canvasResPtLin->SetFillColor(Palette::color_plot_background);
      canvasResPtLin->SetObjectStat(false);
      for (auto itHis=profHisArray_Pt.begin(); itHis!=profHisArray_Pt.end(); itHis++) {

        if (itHis==profHisArray_Pt.begin()) (*itHis)->Draw("PE1");
        else                                (*itHis)->Draw("PE1 SAME");
      }

      myImage = new RootWImage(canvasResPtLin.release(), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("Transverse momentum resolution vs. "+web_etaLetter+" (lin. scale) - " + scenarioName.c_str() + " across "+web_etaLetter);
      myImage->setName(Form("linptres_%s_%s", tag.c_str(), scenario.c_str()));
      myContent->addItem(myImage);

      canvasResPtLog = std::unique_ptr<TCanvas>(new TCanvas);;
      canvasResPtLog->SetGrid(1,1);
      canvasResPtLog->SetLogy(1);
      canvasResPtLog->SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_Pt.begin(); itHis!=profHisArray_Pt.end(); itHis++) {

        if (itHis==profHisArray_Pt.begin()) (*itHis)->Draw("PE1");
        else                                (*itHis)->Draw("PE1 SAME");
      }
      myImage = new RootWImage(canvasResPtLog.release(), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("Transverse momentum resolution vs. "+web_etaLetter+" (log. scale) - " + scenarioName.c_str() + " across "+web_etaLetter);
      myImage->setName(Form("logptres_%s_%s", tag.c_str(), scenario.c_str()));
      myContent->addItem(myImage);

      profHisArray_Pt.clear();

      // b) Resolution in P
      canvasResP = std::unique_ptr<TCanvas>(new TCanvas);
      canvasResP->SetGrid(1,1);
      canvasResP->SetLogy(1);
      canvasResP->SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_P.begin(); itHis!=profHisArray_P.end(); itHis++) {

        if (itHis==profHisArray_P.begin()) (*itHis)->Draw("PE1");
        else                               (*itHis)->Draw("PE1 SAME");
      }

      myImage = new RootWImage(canvasResP.release(), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("Momentum resolution vs. "+web_etaLetter+" - " + scenarioName.c_str() + " across "+web_etaLetter);
      myImage->setName(Form("pres_%s_%s", tag.c_str(), scenario.c_str()));
      myContent->addItem(myImage);

      profHisArray_P.clear();

      // c) Resolution in D0
      canvasResD0 = std::unique_ptr<TCanvas>(new TCanvas);
      canvasResD0->SetGrid(1,1);
      canvasResD0->SetLogy(1);
      canvasResD0->SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_D0.begin(); itHis!=profHisArray_D0.end(); itHis++) {

        if (itHis==profHisArray_D0.begin()) (*itHis)->Draw("PE1");
        else                                (*itHis)->Draw("PE1 SAME");
      }

      myImage = new RootWImage(canvasResD0.release(), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("d0 resolution vs. "+web_etaLetter+" - " + scenarioName.c_str() + " across "+web_etaLetter);
      myImage->setName(Form("d0res_%s_%s", tag.c_str(), scenario.c_str()));
      myContent->addItem(myImage);

      profHisArray_D0.clear();

      // d) Resolution in Z0
      canvasResZ0 = std::unique_ptr<TCanvas>(new TCanvas);
      canvasResZ0->SetGrid(1,1);
      canvasResZ0->SetLogy(1);
      canvasResZ0->SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_Z0.begin(); itHis!=profHisArray_Z0.end(); itHis++) {

        if (itHis==profHisArray_Z0.begin()) (*itHis)->Draw("PE1");
        else                                (*itHis)->Draw("PE1 SAME");
      }

      myImage = new RootWImage(canvasResZ0.release(), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("z0 resolution vs. "+web_etaLetter+" - " + scenarioName.c_str() + " across "+web_etaLetter);
      myImage->setName(Form("z0res_%s_%s", tag.c_str(), scenario.c_str()));
      myContent->addItem(myImage);

      profHisArray_Z0.clear();

      // e) Resolution in Phi0
      canvasResPhi0 = std::unique_ptr<TCanvas>(new TCanvas);
      canvasResPhi0->SetGrid(1,1);
      canvasResPhi0->SetLogy(1);
      canvasResPhi0->SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_Phi0.begin(); itHis!=profHisArray_Phi0.end(); itHis++) {

        if (itHis==profHisArray_Phi0.begin()) (*itHis)->Draw("PE1");
        else                                  (*itHis)->Draw("PE1 SAME");
      }

      myImage = new RootWImage(canvasResPhi0.release(), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment(web_phiLetter + "0 resolution vs. "+web_etaLetter+" - " + scenarioName.c_str() + " across "+web_etaLetter);
      myImage->setName(Form("phi0res_%s_%s", tag.c_str(), scenario.c_str()));
      myContent->addItem(myImage);

      profHisArray_Phi0.clear();

      // f) Resolution in cotg(theta)
      canvasResCotgTh = std::unique_ptr<TCanvas>(new TCanvas);
      canvasResCotgTh->SetGrid(1,1);
      canvasResCotgTh->SetLogy(1);
      canvasResCotgTh->SetFillColor(Palette::color_plot_background);
      for (auto itHis=profHisArray_CotgTheta.begin(); itHis!=profHisArray_CotgTheta.end(); itHis++) {

        if (itHis==profHisArray_CotgTheta.begin()) (*itHis)->Draw("PE1");
        else                                       (*itHis)->Draw("PE1 SAME");
      }

      myImage = new RootWImage(canvasResCotgTh.release(), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("Ctg("+web_thetaLetter+") resolution vs. "+web_etaLetter+" - " + scenarioName.c_str() + " across "+web_etaLetter);
      myImage->setName(Form("cotgThres_%s_%s", tag.c_str(), scenario.c_str()));
      myContent->addItem(myImage);

      profHisArray_CotgTheta.clear();

    } // For scenarios

    // Set Summary content for const Pt
    RootWContent& summaryContent_Pt = myPage->addContent("Summary - const Pt across "+web_etaLetter);
    prepareSummaryTable(tag, "Pt", *myPage, summaryContent_Pt, csvResPt);

    // Set Summary content for const P
    RootWContent& summaryContent_P = myPage->addContent("Summary - const P across "+web_etaLetter, false);
    prepareSummaryTable(tag, "P", *myPage, summaryContent_P, csvResP);
  } // For tag

  return true;
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
      name  = "pT_vs_eta"+any2str(momentum);
      title = "p_{T} resolution versus #eta - const "+scenario+" across #eta;#eta;#delta p_{T}/p_{T} [%]";
    }
    else if (varType=="p") {
      name  = "p_vs_eta"+any2str(momentum);
      title = "p resolution versus #eta - const "+scenario+" across #eta;#eta;#delta p/p [%]";
    }
    else if (varType=="d0") {
      name  = "d0_vs_eta"+any2str(momentum);
      title = "Transverse impact parameter error - const "+scenario+" across #eta;#eta;#delta d_{0} [#mum]";
    }
    else if (varType=="z0") {
      name  = "z0_vs_eta"+any2str(momentum);
      title = "Longitudinal impact parameter error - const "+scenario+" across #eta;#eta;#delta z_{0} [#mum]";
    }
    else if (varType=="phi0") {
      name  = "phi0_vs_eta"+any2str(momentum);
      title = "Track azimuthal angle error - const "+scenario+" across #eta;#eta;#delta #phi [deg]";
    }
    else if (varType=="cotgTheta") {
      name  = "cotgTh_vs_eta"+any2str(momentum);
      title = "Track polar angle error - const "+scenario+" across #eta;#eta;#delta ctg(#theta)";
    }
    std::unique_ptr<TProfile> profHis(new TProfile(name.c_str(), title.c_str(), c_nBins, 0, geom_max_eta_coverage));

    // Set style
    profHis->SetLineWidth(2.);
    profHis->SetMarkerStyle(21);
    profHis->SetMarkerSize(1.);
    profHis->SetLineColor(Palette::colorMomenta(iMomentum));
    profHis->SetMarkerColor(Palette::colorMomenta(iMomentum));

    for (auto& track : col.second) {

      double xVal = track->getEta();
      double yVal = 0;
      if (varType=="pT")        yVal = track->getDeltaPtOverPt()*100; // In percent
      if (varType=="p")         yVal = track->getDeltaPOverP()*100;   // In percent
      if (varType=="d0")        yVal = track->getDeltaD0()/Units::um;
      if (varType=="z0")        yVal = track->getDeltaZ0()/Units::um;
      if (varType=="phi0")      yVal = track->getDeltaPhi()/M_PI*180.; // In degerees
      if (varType=="cotgTheta") yVal = track->getDeltaCtgTheta();

      profHis->Fill(xVal, yVal);

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
    if (varType=="cotgThetat") {
      profHis->SetMaximum(c_max_dCtgTheta);
      profHis->SetMinimum(c_min_dCtgTheta);
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

    RootWTable* myTable = new RootWTable();
    summaryContent.addItem(myTable);
    tableMap[name] = myTable;
    tableMap[name]->setContent(0,0,name);
  }

  // Prepare the cuts for the averages
  ostringstream label;

  // Table explaining the cuts
  cutsSummaryTable.setContent(0,0,"Region: ");
  cutsSummaryTable.setContent(1,0,"Min "+web_etaLetter+":");
  cutsSummaryTable.setContent(2,0,"Max "+web_etaLetter+":");
  for (unsigned int iBorder=0; iBorder<geom_name_eta_regions.size()-1; ++iBorder) {

    cutsSummaryTable.setContent(0,iBorder+1,geom_name_eta_regions[iBorder+1]);
    label.str(""); label << geom_range_eta_regions[iBorder];
    cutsSummaryTable.setContent(1,iBorder+1,label.str());
    label.str(""); label << geom_range_eta_regions[iBorder+1];
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
//      std::ostringstream myLabel;
//
    // Get values for all analyzed momenta
    for (unsigned int iMom=0; iMom<momenta.size(); ++iMom) {

      // Plot average values for different eta regions
      std::vector<double> averagesReal;
      std::vector<double> averagesIdeal;

      // Set table - 1.row & 1.column
      baseColumn = (geom_name_eta_regions.size()-1)*iMom + 1;
      myTable.setContent(0, baseColumn, momenta[iMom]/Units::GeV,0);
      myTable.setColor(0, baseColumn, Palette::colorMomenta(iMom));
//        myIndex.p=momentum[i];
//        myIndex.ideal = false;
//        myGraph = myPlotMap_Pt[myIndex];
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
        canvasReal = webPage.findContent(scenarioNameReal)->findImage("linptres_"+tag+"_withMS_"+scenario);
        canvasIdeal= webPage.findContent(scenarioNameIdeal)->findImage("linptres_"+tag+"_noMS_"+scenario);
        profileName= "pT_vs_eta"+any2str(momenta[iMom],0);
      }
      else if (plotName==std::string(web_deltaLetter+"p/p [%]:             ")) {
        canvasReal = webPage.findContent(scenarioNameReal)->findImage("pres_"+tag+"_withMS_"+scenario);
        canvasIdeal= webPage.findContent(scenarioNameIdeal)->findImage("pres_"+tag+"_noMS_"+scenario);
        profileName= "p_vs_eta"+any2str(momenta[iMom],0);
      }
      else if (plotName==std::string(web_deltaLetter+"d0 ["+web_muLetter+"m]:  ")) {
        canvasReal = webPage.findContent(scenarioNameReal)->findImage("d0res_"+tag+"_withMS_"+scenario);
        canvasIdeal= webPage.findContent(scenarioNameIdeal)->findImage("d0res_"+tag+"_noMS_"+scenario);
        profileName= "d0_vs_eta"+any2str(momenta[iMom],0);
      }
      else if (plotName==std::string(web_deltaLetter+"z0 ["+web_muLetter+"m]:  ")) {
        canvasReal = webPage.findContent(scenarioNameReal)->findImage("z0res_"+tag+"_withMS_"+scenario);
        canvasIdeal= webPage.findContent(scenarioNameIdeal)->findImage("z0res_"+tag+"_noMS_"+scenario);
        profileName= "z0_vs_eta"+any2str(momenta[iMom],0);
      }
      else if (plotName==std::string(web_deltaLetter+web_phiLetter+"0:          ")) {
        canvasReal = webPage.findContent(scenarioNameReal)->findImage("phi0res_"+tag+"_withMS_"+scenario);
        canvasIdeal= webPage.findContent(scenarioNameIdeal)->findImage("phi0res_"+tag+"_noMS_"+scenario);
        profileName= "phi0_vs_eta"+any2str(momenta[iMom],0);
      }
      else if (plotName==std::string(web_deltaLetter+"ctg("+web_thetaLetter+"):")) {
        canvasReal = webPage.findContent(scenarioNameReal)->findImage("cotgThres_"+tag+"_withMS_"+scenario);
        canvasIdeal= webPage.findContent(scenarioNameIdeal)->findImage("cotgThres_"+tag+"_noMS_"+scenario);
        profileName= "cotgTh_vs_eta"+any2str(momenta[iMom],0);
      }

      // Real geometry (active+passive)
      if (canvasReal!=nullptr) for (int i=0; i<canvasReal->GetListOfPrimitives()->GetSize(); ++i) {
        if (std::string(canvasReal->GetListOfPrimitives()->At(i)->ClassName())=="TProfile") {

          // TProfile
          const TProfile* myProfile = (const TProfile*)canvasReal->GetListOfPrimitives()->At(i);

          // If profile histogram for given momentum analyze
          if (std::string(myProfile->GetName())==profileName) {

            averagesReal = averageHis(*myProfile,geom_range_eta_regions);
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

            averagesIdeal = averageHis(*myProfile,geom_range_eta_regions);
          }
        }
      }

      // Fill resolution for different eta regions to a table
      for (unsigned int j=0; j<(geom_name_eta_regions.size()-1); ++j) {

        myTable.setContent(1, baseColumn+j, geom_name_eta_regions[j+1]);
        myTable.setColor(1, baseColumn+j, myColor);
        if (averagesReal.size() > j) {

          myTable.setContent(2, baseColumn+j,averagesReal[j],1);
          myTable.setColor(2, baseColumn+j, myColor);
        }
        if (averagesIdeal.size() > j) {

          myTable.setContent(3, baseColumn+j,averagesIdeal[j],1);
          myTable.setColor(3, baseColumn+j, myColor);
        }
        if ((averagesReal.size() > j)&&(averagesIdeal.size() > j)) {
          myTable.setContent(4, baseColumn+j,averagesReal[j]/averagesIdeal[j],2);
          myTable.setColor(4, baseColumn+j, myColor);
        }
      }

    } // For momenta
  } // For plot names: dpt/pt, dp/p, ...

}

//
// Calculate average values in defined regions for given profile histogram
//
std::vector<double> AnalyzerResolution::averageHis(const TProfile& his, std::vector<double> regions)
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

//
// Material budget visitor - constructor
//
MatTrackVisitor::MatTrackVisitor(Track& matTrack) :
    m_matTrack(matTrack)
{
}

//
// Destructor
//
MatTrackVisitor::~MatTrackVisitor()
{
}

//
// Visit BeamPipe -> update track with beam pipe hit
//
void MatTrackVisitor::visit(const BeamPipe& bp)
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
void MatTrackVisitor::visit(const BarrelModule& m)
{
  analyzeModuleMB(m);
}

//
// Visit EndcapModule (no limits on Rings or Endcaps)
//
void MatTrackVisitor::visit(const EndcapModule& m)
{
  analyzeModuleMB(m);
}

//
// Analyze if module crossed by given track & how much material is in the way
//
void MatTrackVisitor::analyzeModuleMB(const DetectorModule& m)
{
  // Collision detection: material tracks being shot in z+ only, so consider only modules that lie on +Z side
  if (m.maxZ() > 0) {

    XYZVector direction(m_matTrack.getDirection());

    auto pair    = m.checkTrackHits(m_matTrack.getOrigin(), direction);
    auto hitRho  = pair.first.rho();
    auto hitType = pair.second;

    if (hitType!=HitType::NONE) {

      Material material;
      material.radiation   = m.getModuleCap().getRadiationLength();
      material.interaction = m.getModuleCap().getInteractionLength();

      // Fill material map
      double theta = m_matTrack.getTheta();
      double rho   = hitRho;
      double z     = rho/tan(theta);

      // Treat barrel & endcap modules separately
      double tiltAngle = m.tiltAngle();

      if (m.subdet() == BARREL) {

        material.radiation   /= sin(theta + tiltAngle);
        material.interaction /= sin(theta + tiltAngle);

      }
      else if (m.subdet() == ENDCAP) {

        material.radiation   /= cos(theta + tiltAngle - M_PI/2); // Endcap has tiltAngle = pi/2
        material.interaction /= cos(theta + tiltAngle - M_PI/2); // Endcap has tiltAngle = pi/2

      }
      else {
        logWARNING("MatTrackVisitor::analyzeModuleMB -> incorrectly scaled material, unknown module type. Neither barrel or endcap");
      }

      // Create Hit object with appropriate parameters, add to Track t
      HitPtr hit(new Hit(pair.first.R(), &m, hitType));
      hit->setCorrectedMaterial(material);
      m_matTrack.addHit(std::move(hit));
    }
  } // Z>0
}

