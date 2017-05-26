/*
 * AnalyzerPatternReco.cc
 *
 *  Created on: 28. 11. 2016
 *      Author: drasal
 */
#include "AnalyzerPatternReco.h"

#include "BeamPipe.h"
#include "IrradiationMap.h"
#include "MainConfigHandler.h"
#include "MessageLogger.h"
#include "Palette.h"
#include "RootWContent.h"
#include "RootWImage.h"
#include "RootWInfo.h"
#include "RootWPage.h"
#include "RootWSite.h"
#include "SimParms.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "Track.h"
#include "Tracker.h"
#include "TStyle.h"
#include "Units.h"
#include "VisitorMatTrack.h"

//
// AnalyzerPatternReco constructor
//
AnalyzerPatternReco::AnalyzerPatternReco(const Detector& detector) : AnalyzerUnit("AnalyzerPatternReco", detector),
 m_nTracks(0),
 m_etaMin(-1*SimParms::getInstance().getMaxEtaCoverage()),
 m_etaMax(+1*SimParms::getInstance().getMaxEtaCoverage()),
 c_nBins(SimParms::getInstance().getMaxEtaCoverage()/vis_eta_step * 2)  // Default number of bins in histogram from eta=0  to max_eta_coverage)
{
  m_chargedMap = nullptr;
};

//
// AnalyzerPatternReco destructor
//
AnalyzerPatternReco::~AnalyzerPatternReco()
{
  // Clear memory
  if (m_chargedMap!=nullptr) delete m_chargedMap;

  for (auto & iter : m_hisPtBkgContInOut) delete iter;
  m_hisPtBkgContInOut.clear();
  for (auto & iter : m_hisPtBkgContInnerInOut) delete iter;
  m_hisPtBkgContInnerInOut.clear();
  for (auto & iter : m_hisPtBkgContOutIn) delete iter;
  m_hisPtBkgContOutIn.clear();

  for (auto & iter : m_hisPBkgContInOut) delete iter;
  m_hisPBkgContInOut.clear();
  for (auto & iter : m_hisPBkgContInnerInOut) delete iter;
  m_hisPBkgContInnerInOut.clear();
  for (auto & iter : m_hisPBkgContOutIn) delete iter;
  m_hisPBkgContOutIn.clear();
}

//
// AnalyzerPatternReco init method
//
bool AnalyzerPatternReco::init(int nTracks)
{
  // Get occupancy map
  std::string directory = MainConfigHandler::getInstance().getIrradiationDirectory();
  bool chargedMapOK     = checkFile(SimParms::getInstance().chargedMapFile(), directory);

  std::cout << "Reading in: " << directory + "/" + SimParms::getInstance().chargedMapFile() << std::endl;
  if (chargedMapOK) m_chargedMap  = new IrradiationMap(directory + "/" + SimParms::getInstance().chargedMapFile());

  // Set nTracks
  m_nTracks = nTracks;

  if (m_nTracks<=0 || m_trackers.size()==0) {

    std::ostringstream message;
    message << "AnalyzerPatternReco::init(): Number of simulation tracks zero or negative: " << m_nTracks << " or zero number of trackers to be initialized";
    logERROR(message.str());
    return false;
  }
  else {

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      // In-Out
      std::string name = "hisBkgContInOut_pT"+any2str(pIter/Units::GeV);
      m_hisPtBkgContInOut.push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));
      name             = "hisBkgContInnerInOut_pT"+any2str(pIter/Units::GeV);
      m_hisPtBkgContInnerInOut.push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));
      name             = "hisBkgContInOut_p"+any2str(pIter/Units::GeV);
      m_hisPBkgContInOut.push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));
      name             = "hisBkgContInnerInOut_p"+any2str(pIter/Units::GeV);
      m_hisPBkgContInnerInOut.push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));

      // Out-in
      name             = "hisBkgContOutIn_pT"+any2str(pIter/Units::GeV);
      m_hisPtBkgContOutIn.push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));
      name             = "hisBkgContOutIn_p"+any2str(pIter/Units::GeV);
      m_hisPBkgContOutIn.push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));
    }

    m_isInitOK = chargedMapOK;
    return m_isInitOK;
  }
}

//
// AnalyzerPatternReco analyze method
//
bool AnalyzerPatternReco::analyze()
{
  // Check that initialization OK
  if (!m_isInitOK) return false;

  // Random generator
  TRandom3 myDice;
  myDice.SetSeed(random_seed);

  for (int iTrack = 0; iTrack <m_nTracks; iTrack++) {

    // Define track
    Track matTrack;

    double eta   = 0.0 + m_etaMax/m_nTracks*(iTrack+0.5);
    double theta = 2 * atan(exp(-eta));
    double phi   = myDice.Rndm() * M_PI * 2.0; // M_PI/2.;
    double pT    = 100*Units::TeV; // Arbitrarily high number

    //std::cout << "Eta: " << eta << std::endl;

    matTrack.setThetaPhiPt(theta, phi, pT);

    // Use uniform distribution to simulate position of primary interaction (rphi negligible). May be also triangular or gaussian (depends on accelerator design)!
    double zPos = 0.;
    if (SimParms::getInstance().useLumiRegInAnalysis()) zPos = (myDice.Rndm()*2 - 1)*SimParms::getInstance().zErrorIP();
    matTrack.setOrigin(0, 0, zPos);

    // Assign material to the track
    VisitorMatTrack matVisitor(matTrack);
    m_beamPipe->accept(matVisitor);                                 // Assign to material track hit corresponding to beam-pipe
    for (auto iTracker : m_trackers) {

      iTracker->accept(matVisitor);  // Assign to material track hits corresponding to modules
    }

    // Add IP constraint to the track to be used with in-out approach -> TODO: Switch on
    //if (SimParms::getInstance().useIPConstraint()) matTrack.addIPConstraint(SimParms::getInstance().rphiErrorIP(), SimParms::getInstance().zErrorIP());

    // For each momentum/transverse momentum compute
    int iMomentum = 0;
    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      // 2 options: pT & p
      for (int pOption=0;pOption<2;pOption++) {

        // InOut or OutIn approach
        for (int approachOption=0; approachOption<2; approachOption++) {

          bool propagOutIn = approachOption;

          //if (!propagOutIn) std::cout << "InOut approach" << std::endl;
          //else              std::cout << "OutIn approach" << std::endl;

          // Reset track total probability & calculate p/pT based on option
          double pNotContamTot      = 1;
          double pNotContamTotInner = 1;

          if (pOption==0) pT = pIter;             // pT option
          else            pT = pIter*sin(theta);  // p option

          // Set track & prune hits
          Track track(matTrack);
          track.resetPt(pT);

          //
          // Pattern recognition procedure
          bool useIP        = true;
          int nMeasuredHits = 0;
          if (!propagOutIn) nMeasuredHits = track.getNMeasuredHits("all", !useIP);
          else              nMeasuredHits = track.getNMeasuredHits("all", !useIP);

          //track.printHits();
          //std::cout << ">> " << nMeasuredHits << " dD0=" << track.getDeltaD(0.0)/Units::um << " dZ0=" << track.getDeltaZ(0.0)/Units::um << std::endl;

          // Start with first 3 hits and end with N-1 hits (C counting from zero)
          bool testTriplet = true;
          for (int iHit=2; iHit<nMeasuredHits-1; iHit++) {

            // Start with triplet coming from 3 different layers -> avoid using 2 close measurements from overlapping modules in 1 layer
            if (testTriplet && !isTripletFromDifLayers(track, iHit, propagOutIn)) continue;
            else testTriplet = false;

            int nHitsUsed = iHit+1; // (C counting from zero)

            // Keep first/last N hits (based on approach) and find paramaters of the next hit
            double nextRPos    = 0;
            double nextZPos    = 0;
            double nextHitTilt = 0;
            double A           = 0; // r_i/2R
            std::string iHitID = "";

            if (!propagOutIn) {

              nextRPos    = track.getMeasurableOrIPHit(iHit+1)->getRPos();
              nextZPos    = track.getMeasurableOrIPHit(iHit+1)->getZPos();
              nextHitTilt = track.getMeasurableOrIPHit(iHit+1)->getTilt();
              if (SimParms::getInstance().isMagFieldConst()) A = track.getMeasurableOrIPHit(iHit+1)->getRPos()/2./track.getRadius(track.getMeasurableOrIPHit(iHit+1)->getZPos());  // r_i/2R

              iHitID      = track.getMeasurableOrIPHit(iHit+1)->getDetName() +
                            std::string(track.getMeasurableOrIPHit(iHit+1)->isBarrel() ? "_L_" : "_D_") +
                            any2str(track.getMeasurableOrIPHit(iHit+1)->getLayerOrDiscID());
            }
            else {

              nextRPos    = track.getRMeasurableOrIPHit(iHit+1)->getRPos();
              nextZPos    = track.getRMeasurableOrIPHit(iHit+1)->getZPos();
              nextHitTilt = track.getRMeasurableOrIPHit(iHit+1)->getTilt();
              if (SimParms::getInstance().isMagFieldConst()) A = track.getRMeasurableOrIPHit(iHit+1)->getRPos()/2./track.getRadius(track.getRMeasurableOrIPHit(iHit+1)->getZPos());  // r_i/2R
              iHitID      = track.getRMeasurableOrIPHit(iHit+1)->getDetName() +
                            std::string(track.getRMeasurableOrIPHit(iHit+1)->isBarrel() ? "_L_" : "_D_") +
                            any2str(track.getRMeasurableOrIPHit(iHit+1)->getLayerOrDiscID());
            }

            // Calculate dRPhi & dZ
            int    nPU    = SimParms::getInstance().numMinBiasEvents();
            double flux   = m_chargedMap->calculateIrradiationZR(nextZPos,nextRPos)*nPU;

            double corrFactor = cos(nextHitTilt) + track.getCotgTheta()/sqrt(1-A*A)*sin(nextHitTilt);
            double dZProj     = track.getDeltaZ(nextRPos,propagOutIn)/corrFactor;
            double dDProj     = track.getDeltaD(nextRPos,propagOutIn);

            // Print info
            //std::cout << ">> " << " R=" << nextRPos/Units::mm << " dD=" << dDProj/Units::um << " " << " Z=" << nextZPos/Units::mm << " dZ=" << dZProj/Units::um << std::endl;

            // Calculate how many sigmas does one need to get in 2D Gauss. 95% coverge, prob=0.95, alpha (p-value)=0.05
            // F(mu + n*sigma) - F(mu - n*sigma) = erf(n/sqrt(2)) --> 1 dimension
            // (x/sigX)^2 + (y/sigY)^2 = N^2 --> Follows chi2 distr with r=2 degrees of freedom --> Calculate 0.95 confidence level -> nSigmaFactor is then sqrt(chi2(r=2,p=0.95))
            // Or we want to get contour for ellipse being (x/sigX)^2 + (y/sigY)^2 = N^2 --> Calculate 0.95 confidence level
            // Hence calculate integral: int(Gaus(0,sigX)|-a;+a * int(Gaus(0,sigY)|-b;+b --> Make transformation (x/sigX)^2 + (y/sigY)^2 = z^2
            // The integral is then: int(z.exp(-z^2)/2Pi)|0;N * int(dtheta)|0;2Pi = 1 - exp(-N^2/2) = 1 - alpha --> N^2 = -2*ln(alpha) --> nSigmaFactor is then sqrt(-2*ln(1-p))
            static double nSigmaFactor = sqrt(-2*log(0.05));

            // Calculate errorEllipse multiplied by nSigmaFactor(RPhi)*nSigmaFactor(Z)
            double errorEllipse = M_PI*dZProj*dDProj*nSigmaFactor*nSigmaFactor;

            // Accumulate probibilites accross all hits not to be contaminated anywhere -> final probability of being contaminated is 1 - pContam
            double pContam = flux*errorEllipse;
            if (pContam>1) pContam = 1.0;
            pNotContamTot *= 1-pContam;

            // Output warning that some results may be inconsistent
            if (dZProj<=0 || dDProj<=0 || (1-pContam)<0 || pContam>1) logWARNING("Track propagator provided dZ or dD or probability lower than zero, check technique!");

            //
            // Fill d0proj, z0proj & pContam for individual layers/discs to histograms

            // Probability of inner tracker only
            if (iHitID.find("Inner")!=std::string::npos) pNotContamTotInner *= 1-pContam;

            // Create artificial map identifier with extra characters to have innermost first, then outermost & then fwd
            std::string iHitIDMap = "";
            if      (iHitID.find("Inner")!=std::string::npos) iHitIDMap = "A_" + iHitID;
            else if (iHitID.find("Outer")!=std::string::npos) iHitIDMap = "B_" + iHitID;
            else                                              iHitIDMap = "C_" + iHitID;

            // Pt & InOut
            if (pOption==0 && approachOption==0) {

              // Create profile histograms if don't exist yet
              if (m_hisPtHitDProjInOut.find(iHitIDMap)==m_hisPtHitDProjInOut.end()) {
                for (int iMom=0; iMom<MainConfigHandler::getInstance().getMomenta().size(); iMom++) {

                  std::string name = "hisPt_" + any2str(iMom) + "_Hit_" + iHitID + "_DProjInOut";
                  m_hisPtHitDProjInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));

                  name = "hisPt" + any2str(iMom) + "_Hit_" + iHitID + "_ZProjInOut";
                  m_hisPtHitZProjInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));

                  name = "hisPt" + any2str(iMom) + "_Hit_" + iHitID + "_ProbContamInOut";
                  m_hisPtHitProbContamInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));
                }
              }
              m_hisPtHitDProjInOut[iHitIDMap][iMomentum]->Fill(eta, dDProj/Units::um);
              m_hisPtHitZProjInOut[iHitIDMap][iMomentum]->Fill(eta, dZProj/Units::um);
              m_hisPtHitProbContamInOut[iHitIDMap][iMomentum]->Fill(eta, pContam);
            }
            // P & InOut
            if (pOption==1 && approachOption==0) {

              // Create profile histograms if don't exist yet
              if (m_hisPHitDProjInOut.find(iHitIDMap)==m_hisPHitDProjInOut.end()) {
                for (int iMom=0; iMom<MainConfigHandler::getInstance().getMomenta().size(); iMom++) {

                  std::string name = "hisP_" + any2str(iMom) + "_Hit_" + iHitID + "_DProjInOut";
                  m_hisPHitDProjInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));

                  name = "hisP" + any2str(iMom) + "_Hit_" + iHitID + "_ZProjInOut";
                  m_hisPHitZProjInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));

                  name = "hisP" + any2str(iMom) + "_Hit_" + iHitID + "_ProbContamInOut";
                  m_hisPHitProbContamInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));
                }
              }
              // TODO: Check that above zero!!!
              m_hisPHitDProjInOut[iHitIDMap][iMomentum]->Fill(eta, dDProj/Units::um);
              m_hisPHitZProjInOut[iHitIDMap][iMomentum]->Fill(eta, dZProj/Units::um);
              m_hisPHitProbContamInOut[iHitIDMap][iMomentum]->Fill(eta, pContam);
            }
            // Pt & OutIn
            if (pOption==0 && approachOption==1) {

              // Create profile histograms if don't exist yet
              if (m_hisPtHitDProjOutIn.find(iHitIDMap)==m_hisPtHitDProjOutIn.end()) {
                for (int iMom=0; iMom<MainConfigHandler::getInstance().getMomenta().size(); iMom++) {

                  std::string name = "hisPt_" + any2str(iMom) + "_Hit_" + iHitID + "_DProjOutIn";
                  m_hisPtHitDProjOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));

                  name = "hisPt" + any2str(iMom) + "_Hit_" + iHitID + "_ZProjOutIn";
                  m_hisPtHitZProjOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));

                  name = "hisPt" + any2str(iMom) + "_Hit_" + iHitID + "_ProbContamOutIn";
                  m_hisPtHitProbContamOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));
                }
              }
              m_hisPtHitDProjOutIn[iHitIDMap][iMomentum]->Fill(eta, dDProj/Units::um);
              m_hisPtHitZProjOutIn[iHitIDMap][iMomentum]->Fill(eta, dZProj/Units::um);
              m_hisPtHitProbContamOutIn[iHitIDMap][iMomentum]->Fill(eta, pContam);
            }
            // P & OutIn
            if (pOption==1 && approachOption==1) {

              // Create profile histograms if don't exist yet
              if (m_hisPHitDProjOutIn.find(iHitIDMap)==m_hisPHitDProjOutIn.end()) {
                for (int iMom=0; iMom<MainConfigHandler::getInstance().getMomenta().size(); iMom++) {

                  std::string name = "hisP_" + any2str(iMom) + "_Hit_" + iHitID + "_DProjOutIn";
                  m_hisPHitDProjOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));

                  name = "hisP" + any2str(iMom) + "_Hit_" + iHitID + "_ZProjOutIn";
                  m_hisPHitZProjOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));

                  name = "hisP" + any2str(iMom) + "_Hit_" + iHitID + "_ProbContamOutIn";
                  m_hisPHitProbContamOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),c_nBins, 0, SimParms::getInstance().getMaxEtaCoverage()));
                }
              }
              // TODO: Check that above zero!!!
              m_hisPHitDProjOutIn[iHitIDMap][iMomentum]->Fill(eta, dDProj/Units::um);
              m_hisPHitZProjOutIn[iHitIDMap][iMomentum]->Fill(eta, dZProj/Units::um);
              m_hisPHitProbContamOutIn[iHitIDMap][iMomentum]->Fill(eta, pContam);
            }
          } // Pattern reco loop

          //
          // Calculate total contamination based on different options: p/pT, in-out/out-in
          double pContamTot      = 1-pNotContamTot;
          double pContamInnerTot = 1-pNotContamTotInner;
          if (pOption==0) {

            if (approachOption==0) {
              m_hisPtBkgContInOut[iMomentum]->Fill(eta,pContamTot);
              m_hisPtBkgContInnerInOut[iMomentum]->Fill(eta,pContamInnerTot);
            }
            else {
              m_hisPtBkgContOutIn[iMomentum]->Fill(eta,pContamTot);
            }
          }
          else  {

            if (approachOption==0) {
              m_hisPBkgContInOut[iMomentum]->Fill(eta,pContamTot);
              m_hisPBkgContInnerInOut[iMomentum]->Fill(eta,pContamInnerTot);
            }
            else {
              m_hisPBkgContOutIn[iMomentum]->Fill(eta,pContamTot);
            }
          }

        } //In-Out, Out-In approach - options loop
      } //pT/p option loop

      iMomentum++;
    } // Momenta loop

  } // Tracks loop

  m_isAnalysisOK = true;
  return m_isAnalysisOK;
}

//
// AnalyzerPatternReco visualization method
//
bool AnalyzerPatternReco::visualize(RootWSite& webSite)
{
  // Check that initialization & analysis OK
  if (!m_isInitOK && !m_isAnalysisOK) return false;

  // Set Rainbow palette for drawing
  Palette::setRootPalette();

  std::string pageTitle   = "PatternReco";
  std::string pageAddress = "indexPatternReco.html";

  int webPriority         = web_priority_PR;

  RootWPage& myPage = webSite.addPage(pageTitle, webPriority);
  myPage.setAddress(pageAddress);

  // Canvases
  gStyle->SetGridStyle(style_grid);
  gStyle->SetGridColor(Palette::color_hard_grid);
  gStyle->SetOptStat(0);

  //
  // Pt option
  RootWContent& myContentPlotsPt = myPage.addContent("Pt: Track purity & probability of contamination by bkg hits when propagating track through tracker", true);

  // Summary
  RootWInfo& myInfoPt = myContentPlotsPt.addInfo("Number of minimum bias events per bunch crossing");
  myInfoPt.setValue(SimParms::getInstance().numMinBiasEvents(), 0);

  // a) Bkg contamination probability -> In-Out approach
  TCanvas canvasPtBkgContInOut("canvasPtBkgContInOut","",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
  canvasPtBkgContInOut.SetGrid(1,1);
  canvasPtBkgContInOut.SetLogy(0);
  canvasPtBkgContInOut.SetFillColor(Palette::color_plot_background);
  canvasPtBkgContInOut.SetObjectStat(false);

  // For each momentum/transverse momentum compute
  int iMomentum = 0;
  gStyle->SetTitleW(0.9);
  TLegend* legendPtInOut = new TLegend(0.11,0.66,0.36,0.89,"p_{T} options (InOut+IP, full TRK):");
  legendPtInOut->SetTextSize(0.025);

  for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

    m_hisPtBkgContInOut[iMomentum]->SetNameTitle(std::string("hisPtBkgContInOut"+any2str(iMomentum)).c_str(),"In-Out: Bkg contamination prob. in 95% area of 2D error ellipse accumulated accross N layers;#eta;1 - #Pi_{i=1}^{N} (1-p^{i}_{bkg_95%})");
    m_hisPtBkgContInOut[iMomentum]->SetLineWidth(2.);
    m_hisPtBkgContInOut[iMomentum]->SetMarkerStyle(21);
    m_hisPtBkgContInOut[iMomentum]->SetMarkerSize(1.);
    m_hisPtBkgContInOut[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
    m_hisPtBkgContInOut[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

    if (iMomentum==0) m_hisPtBkgContInOut[iMomentum]->Draw("PE1");
    else              m_hisPtBkgContInOut[iMomentum]->Draw("SAME PE1");

    m_hisPtBkgContInOut[iMomentum]->GetYaxis()->SetRangeUser(0,1.3);
    legendPtInOut->AddEntry(m_hisPtBkgContInOut[iMomentum],std::string(any2str(pIter/Units::GeV)+"GeV").c_str(),"lp");

    iMomentum++;
  }
  legendPtInOut->Draw("SAME");

  RootWImage& myImagePtBkgContInOut = myContentPlotsPt.addImage(canvasPtBkgContInOut, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
  myImagePtBkgContInOut.setComment("In-Out approach for pT: Bkg contamination prob. in 95% area of 2D error ellipse accumulated across N layers");
  myImagePtBkgContInOut.setName("bkg_pt_pContam_inout");

  // b) Bkg contamination probability for inner tracker only -> In-Out approach
  TCanvas canvasPtBkgContInnerInOut("canvasPtBkgContInnerInOut","",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
  canvasPtBkgContInnerInOut.SetGrid(1,1);
  canvasPtBkgContInnerInOut.SetLogy(0);
  canvasPtBkgContInnerInOut.SetFillColor(Palette::color_plot_background);
  canvasPtBkgContInnerInOut.SetObjectStat(false);

  // For each momentum/transverse momentum compute
  iMomentum = 0;
  gStyle->SetTitleW(0.9);
  legendPtInOut->SetEntryLabel("p_{T} options (InOut+IP, inner TRK):");
  legendPtInOut->SetTextSize(0.025);

  for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

    m_hisPtBkgContInnerInOut[iMomentum]->SetNameTitle(std::string("hisPtBkgContInnerInOut"+any2str(iMomentum)).c_str(),"In-Out: Bkg contamination prob. in 95% area of 2D error ellipse accumulated accross N inner trk layers;#eta;1 - #Pi_{i=1}^{N} (1-p^{i}_{bkg_95%})");
    m_hisPtBkgContInnerInOut[iMomentum]->SetLineWidth(2.);
    m_hisPtBkgContInnerInOut[iMomentum]->SetMarkerStyle(21);
    m_hisPtBkgContInnerInOut[iMomentum]->SetMarkerSize(1.);
    m_hisPtBkgContInnerInOut[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
    m_hisPtBkgContInnerInOut[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

    if (iMomentum==0) m_hisPtBkgContInnerInOut[iMomentum]->Draw("PE1");
    else              m_hisPtBkgContInnerInOut[iMomentum]->Draw("SAME PE1");

    m_hisPtBkgContInnerInOut[iMomentum]->GetYaxis()->SetRangeUser(0,1.3);

    iMomentum++;
  }
  legendPtInOut->Draw("SAME");

  RootWImage& myImagePtBkgContInnerInOut = myContentPlotsPt.addImage(canvasPtBkgContInnerInOut, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
  myImagePtBkgContInnerInOut.setComment("In-Out approach for pT & inner tracker only: Bkg contamination prob. in 95% area of 2D error ellipse");
  myImagePtBkgContInnerInOut.setName("bkg_pt_pContam_inner_inout");

  // c) Bkg contamination probability -> Out-In approach
  TCanvas canvasPtBkgContOutIn("canvasPtBkgContOutIn","",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
  canvasPtBkgContOutIn.SetGrid(1,1);
  canvasPtBkgContOutIn.SetLogy(0);
  canvasPtBkgContOutIn.SetFillColor(Palette::color_plot_background);
  canvasPtBkgContOutIn.SetObjectStat(false);

  // For each momentum/transverse momentum compute
  iMomentum = 0;
  gStyle->SetTitleW(0.9);
  TLegend* legendPtOutIn = new TLegend(0.11,0.66,0.31,0.89,"p_{T} options (OutIn, full TRK):");
  legendPtOutIn->SetTextSize(0.025);

  for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

    m_hisPtBkgContOutIn[iMomentum]->SetNameTitle(std::string("hisPtBkgContOutIn"+any2str(iMomentum)).c_str(),"Out-In: Bkg contamination prob. in 95% area of 2D error ellipse accumulated accross N layers;#eta;1 - #Pi_{i=1}^{N} (1-p^{i}_{bkg_95%})");
    m_hisPtBkgContOutIn[iMomentum]->SetLineWidth(2.);
    m_hisPtBkgContOutIn[iMomentum]->SetMarkerStyle(21);
    m_hisPtBkgContOutIn[iMomentum]->SetMarkerSize(1.);
    m_hisPtBkgContOutIn[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
    m_hisPtBkgContOutIn[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

    if (iMomentum==0) m_hisPtBkgContOutIn[iMomentum]->Draw("PE1");
    else              m_hisPtBkgContOutIn[iMomentum]->Draw("SAME PE1");

    m_hisPtBkgContOutIn[iMomentum]->GetYaxis()->SetRangeUser(0,1.3);
    legendPtOutIn->AddEntry(m_hisPtBkgContOutIn[iMomentum],std::string(any2str(pIter/Units::GeV)+"GeV").c_str(),"lp");

    iMomentum++;
  }
  legendPtOutIn->Draw("SAME");

  RootWImage& myImagePtBkgContOutIn = myContentPlotsPt.addImage(canvasPtBkgContOutIn, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
  myImagePtBkgContOutIn.setComment("Out-In approach for pT: Bkg contamination prob. in 95% area of 2D error ellipse accumulated across N layers");
  myImagePtBkgContOutIn.setName("bkg_pt_pContam_outin");

  //
  // Detail on pt in-out studies

  // a) D0
  RootWContent& myContentPlotsPtD0InOut = myPage.addContent("Pt in-out - D: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPtHitDProjInOut) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPtD0InOut(std::string("canvasPtD0InOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPtD0InOut.SetGrid(1,1);
    canvasPtD0InOut.SetLogy(0);
    canvasPtD0InOut.SetFillColor(Palette::color_plot_background);
    canvasPtD0InOut.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPtD0InOut"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{R-#Phi} from previous layers/discs;#eta; #sigma_{R-#Phi} [#mum]").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

      iMomentum++;
    }
    legendPtInOut->SetX1(0.64);
    legendPtInOut->SetX2(0.89);
    legendPtInOut->SetEntryLabel("p_{T} options (InOut+IP, full TRK):");
    legendPtInOut->Draw("SAME");

    RootWImage& myImagePtD0InOut = myContentPlotsPtD0InOut.addImage(canvasPtD0InOut, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePtD0InOut.setComment(std::string("Detector: "+name+" - an extrapolated sigma in R-Phi from previous layers/discs in in-out approach.").c_str());
    myImagePtD0InOut.setName(std::string("bkg_pt_d0_inout_"+name).c_str());
  }

  // b) Z0
  RootWContent& myContentPlotsPtZ0InOut = myPage.addContent("Pt in-out - Z: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPtHitZProjInOut) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPtZ0InOut(std::string("canvasPtZ0InOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPtZ0InOut.SetGrid(1,1);
    canvasPtZ0InOut.SetLogy(0);
    canvasPtZ0InOut.SetFillColor(Palette::color_plot_background);
    canvasPtZ0InOut.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPtZ0InOut"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{Z} from previous layers/discs;#eta; #sigma_{Z} [#mum]").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

      iMomentum++;
    }
    legendPtInOut->Draw("SAME");

    RootWImage& myImagePtZ0InOut = myContentPlotsPtZ0InOut.addImage(canvasPtZ0InOut, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePtZ0InOut.setComment(std::string("Detector: "+name+" - an extrapolated sigma in Z from previous layers/discs in in-out approach.").c_str());
    myImagePtZ0InOut.setName(std::string("bkg_pt_z0_inout_"+name).c_str());
  }

  // c) pContamination
  RootWContent& myContentPlotsPtProbContamInOut = myPage.addContent("Pt in-out - Bkg contamination prob.: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPtHitProbContamInOut) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPtProbContamInOut(std::string("canvasPtProbContamInOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPtProbContamInOut.SetGrid(1,1);
    canvasPtProbContamInOut.SetLogy(0);
    canvasPtProbContamInOut.SetFillColor(Palette::color_plot_background);
    canvasPtProbContamInOut.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPtProbContamInOut"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: bkg contamination prob. as error ellipse extrap. from previous layers/discs;#eta; probability").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,0.1);

      iMomentum++;
    }
    legendPtInOut->SetX1(0.11);
    legendPtInOut->SetX2(0.36);
    legendPtInOut->Draw("SAME");

    RootWImage& myImagePtProbContamInOut = myContentPlotsPtProbContamInOut.addImage(canvasPtProbContamInOut, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePtProbContamInOut.setComment(std::string("Detector: "+name+" - bkg contamination prob. as error ellipse extrap. from previous layers/discs in in-out approach.").c_str());
    myImagePtProbContamInOut.setName(std::string("bkg_pt_pContam_inout_"+name).c_str());
  }

  //
  // Detail on pt out-in studies

  // a) D0
  RootWContent& myContentPlotsPtD0OutIn = myPage.addContent("Pt out-in - D: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPtHitDProjOutIn) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPtD0OutIn(std::string("canvasPtD0OutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPtD0OutIn.SetGrid(1,1);
    canvasPtD0OutIn.SetLogy(0);
    canvasPtD0OutIn.SetFillColor(Palette::color_plot_background);
    canvasPtD0OutIn.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPtD0OutIn"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{R-#Phi} from previous layers/discs;#eta; #sigma_{R-#Phi} [#mum]").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

      iMomentum++;
    }
    legendPtOutIn->SetX1(0.64);
    legendPtOutIn->SetX2(0.89);
    legendPtOutIn->SetEntryLabel("p_{T} options (OutIn, full TRK):");
    legendPtOutIn->Draw("SAME");

    RootWImage& myImagePtD0OutIn = myContentPlotsPtD0OutIn.addImage(canvasPtD0OutIn, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePtD0OutIn.setComment(std::string("Detector: "+name+" - an extrapolated sigma in R-Phi from previous layers/discs in out-in approach.").c_str());
    myImagePtD0OutIn.setName(std::string("bkg_pt_d0_outin_"+name).c_str());
  }

  // b) Z0
  RootWContent& myContentPlotsPtZ0OutIn = myPage.addContent("Pt out-in - Z: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPtHitZProjOutIn) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPtZ0OutIn(std::string("canvasPtZ0OutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPtZ0OutIn.SetGrid(1,1);
    canvasPtZ0OutIn.SetLogy(0);
    canvasPtZ0OutIn.SetFillColor(Palette::color_plot_background);
    canvasPtZ0OutIn.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPtZ0OutIn"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{Z} from previous layers/discs;#eta; #sigma_{Z} [#mum]").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

      iMomentum++;
    }
    legendPtOutIn->Draw("SAME");

    RootWImage& myImagePtZ0OutIn = myContentPlotsPtZ0OutIn.addImage(canvasPtZ0OutIn, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePtZ0OutIn.setComment(std::string("Detector: "+name+" - an extrapolated sigma in Z from previous layers/discs in out-in approach.").c_str());
    myImagePtZ0OutIn.setName(std::string("bkg_pt_z0_outin_"+name).c_str());
  }

  // c) pContamination
  RootWContent& myContentPlotsPtProbContamOutIn = myPage.addContent("Pt out-in - Bkg contamination prob.: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPtHitProbContamOutIn) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPtProbContamOutIn(std::string("canvasPtProbContamOutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPtProbContamOutIn.SetGrid(1,1);
    canvasPtProbContamOutIn.SetLogy(0);
    canvasPtProbContamOutIn.SetFillColor(Palette::color_plot_background);
    canvasPtProbContamOutIn.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPtProbContamOutIn"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: bkg contamination prob. as error ellipse extrap. from previous layers/discs;#eta; probability").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,0.1);

      iMomentum++;
    }
    legendPtOutIn->SetX1(0.11);
    legendPtOutIn->SetX2(0.36);
    legendPtOutIn->Draw("SAME");

    RootWImage& myImagePtProbContamOutIn = myContentPlotsPtProbContamOutIn.addImage(canvasPtProbContamOutIn, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePtProbContamOutIn.setComment(std::string("Detector: "+name+" - bkg contamination prob. as error ellipse extrap. from previous layers/discs in out-in approach.").c_str());
    myImagePtProbContamOutIn.setName(std::string("bkg_pt_pContam_outin_"+name).c_str());
  }

  //
  // P option
  RootWContent& myContentPlotsP = myPage.addContent("P: Track purity & probability of contamination by bkg hits when propagating track through tracker", true);

  // d) Summary
  RootWInfo& myInfoP = myContentPlotsP.addInfo("Number of minimum bias events per bunch crossing");
  myInfoP.setValue(SimParms::getInstance().numMinBiasEvents(), 0);

  // a) Bkg contamination probability - In-Out approach
  TCanvas canvasPBkgContInOut("canvasPBkgContInOut","",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
  canvasPBkgContInOut.SetGrid(1,1);
  canvasPBkgContInOut.SetLogy(0);
  canvasPBkgContInOut.SetFillColor(Palette::color_plot_background);
  canvasPBkgContInOut.SetObjectStat(false);

  // For each momentum/transverse momentum compute
  iMomentum = 0;
  gStyle->SetTitleW(0.9);
  TLegend* legendPInOut = new TLegend(0.11,0.66,0.36,0.89,"p options (InOut+IP, full TRK):");
  legendPInOut->SetTextSize(0.025);

  for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

    m_hisPBkgContInOut[iMomentum]->SetNameTitle(std::string("hisPBkgContInOut"+any2str(iMomentum)).c_str(),"In-Out: Bkg contamination prob. in 95% area of 2D error ellipse accumulated accross N layers;#eta;1 - #Pi_{i=1}^{N} (1-p^{i}_{bkg_95%})");
    m_hisPBkgContInOut[iMomentum]->SetLineWidth(2.);
    m_hisPBkgContInOut[iMomentum]->SetMarkerStyle(21);
    m_hisPBkgContInOut[iMomentum]->SetMarkerSize(1.);
    m_hisPBkgContInOut[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
    m_hisPBkgContInOut[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

    if (iMomentum==0) m_hisPBkgContInOut[iMomentum]->Draw("PE1");
    else              m_hisPBkgContInOut[iMomentum]->Draw("SAME PE1");

    m_hisPBkgContInOut[iMomentum]->GetYaxis()->SetRangeUser(0,1.3);
    legendPInOut->AddEntry(m_hisPBkgContInOut[iMomentum],std::string(any2str(pIter/Units::GeV)+"GeV").c_str(),"lp");

    iMomentum++;
  }
  legendPInOut->Draw("SAME");

  RootWImage& myImagePBkgContInOut = myContentPlotsP.addImage(canvasPBkgContInOut, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
  myImagePBkgContInOut.setComment("In-Out approach for p: Bkg contamination prob. in 95% area of 2D error ellipse accumulated across N layers");
  myImagePBkgContInOut.setName("bkg_p_pContam_inout");

  // b) Bkg contamination probability for inner tracker only - In-Out approach
  TCanvas canvasPBkgContInnerInOut("canvasPBkgContInnerInOut","",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
  canvasPBkgContInnerInOut.SetGrid(1,1);
  canvasPBkgContInnerInOut.SetLogy(0);
  canvasPBkgContInnerInOut.SetFillColor(Palette::color_plot_background);
  canvasPBkgContInnerInOut.SetObjectStat(false);

  // For each momentum/transverse momentum compute
  iMomentum = 0;
  gStyle->SetTitleW(0.9);
  legendPInOut->SetEntryLabel("p options (InOut+IP, inner TRK):");
  legendPInOut->SetTextSize(0.025);

  for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

    m_hisPBkgContInnerInOut[iMomentum]->SetNameTitle(std::string("hisPBkgContInnerInOut"+any2str(iMomentum)).c_str(),"In-Out: Bkg contamination prob. in 95% area of 2D error ellipse accumulated accross N inner trk layers;#eta;1 - #Pi_{i=1}^{N} (1-p^{i}_{bkg_95%})");
    m_hisPBkgContInnerInOut[iMomentum]->SetLineWidth(2.);
    m_hisPBkgContInnerInOut[iMomentum]->SetMarkerStyle(21);
    m_hisPBkgContInnerInOut[iMomentum]->SetMarkerSize(1.);
    m_hisPBkgContInnerInOut[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
    m_hisPBkgContInnerInOut[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

    if (iMomentum==0) m_hisPBkgContInnerInOut[iMomentum]->Draw("PE1");
    else              m_hisPBkgContInnerInOut[iMomentum]->Draw("SAME PE1");

    m_hisPBkgContInnerInOut[iMomentum]->GetYaxis()->SetRangeUser(0,1.3);

    iMomentum++;
  }
  legendPInOut->Draw("SAME");

  RootWImage& myImagePBkgContInnerInOut = myContentPlotsP.addImage(canvasPBkgContInnerInOut, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
  myImagePBkgContInnerInOut.setComment("In-Out approach for p & inner tracker only: Bkg contamination prob. in 95% area of 2D error ellipse");
  myImagePBkgContInnerInOut.setName("bkg_p_pContam_inner_inout");

  // c) Bkg contamination probability - Out-In approach
  TCanvas canvasPBkgContOutIn("canvasPBkgContOutIn","",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
  canvasPBkgContOutIn.SetGrid(1,1);
  canvasPBkgContOutIn.SetLogy(0);
  canvasPBkgContOutIn.SetFillColor(Palette::color_plot_background);
  canvasPBkgContOutIn.SetObjectStat(false);

  // For each momentum/transverse momentum compute
  iMomentum = 0;
  gStyle->SetTitleW(0.9);
  TLegend* legendPOutIn = new TLegend(0.11,0.66,0.31,0.89,"p options (OutIn, full TRK):");
  legendPInOut->SetTextSize(0.025);

  for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

    m_hisPBkgContOutIn[iMomentum]->SetNameTitle(std::string("hisPBkgContOutIn"+any2str(iMomentum)).c_str(),"In-Out: Bkg contamination prob. in 95% area of 2D error ellipse accumulated accross N layers;#eta;1 - #Pi_{i=1}^{N} (1-p^{i}_{bkg_95%})");
    m_hisPBkgContOutIn[iMomentum]->SetLineWidth(2.);
    m_hisPBkgContOutIn[iMomentum]->SetMarkerStyle(21);
    m_hisPBkgContOutIn[iMomentum]->SetMarkerSize(1.);
    m_hisPBkgContOutIn[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
    m_hisPBkgContOutIn[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

    if (iMomentum==0) m_hisPBkgContOutIn[iMomentum]->Draw("PE1");
    else              m_hisPBkgContOutIn[iMomentum]->Draw("SAME PE1");

    m_hisPBkgContOutIn[iMomentum]->GetYaxis()->SetRangeUser(0,1.3);
    legendPOutIn->AddEntry(m_hisPBkgContOutIn[iMomentum],std::string(any2str(pIter/Units::GeV)+"GeV").c_str(),"lp");

    iMomentum++;
  }
  legendPOutIn->Draw("SAME");

  RootWImage& myImagePBkgContOutIn = myContentPlotsP.addImage(canvasPBkgContOutIn, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
  myImagePBkgContOutIn.setComment("Out-In approach for p: Bkg contamination prob. in 95% area of 2D error ellipse accumulated across N layers");
  myImagePBkgContOutIn.setName("bkg_p_pContam_outin");

  //
  // Detail on p in-out studies

  // a) D0
  RootWContent& myContentPlotsPD0InOut = myPage.addContent("P in-out - D: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPHitDProjInOut) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPD0InOut(std::string("canvasPD0InOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPD0InOut.SetGrid(1,1);
    canvasPD0InOut.SetLogy(0);
    canvasPD0InOut.SetFillColor(Palette::color_plot_background);
    canvasPD0InOut.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPD0InOut"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{R-#Phi} from previous layers/discs;#eta; #sigma_{R-#Phi} [#mum]").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

      iMomentum++;
    }
    legendPInOut->SetX1(0.64);
    legendPInOut->SetX2(0.89);
    legendPInOut->SetEntryLabel("p options (InOut+IP, full TRK):");
    legendPInOut->Draw("SAME");

    RootWImage& myImagePD0InOut = myContentPlotsPD0InOut.addImage(canvasPD0InOut, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePD0InOut.setComment(std::string("Detector: "+name+" - an extrapolated sigma in R-Phi from previous layers/discs in in-out approach.").c_str());
    myImagePD0InOut.setName(std::string("bkg_p_d0_inout_"+name).c_str());
  }

  // b) Z0
  RootWContent& myContentPlotsPZ0InOut = myPage.addContent("P in-out - Z: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPHitZProjInOut) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPZ0InOut(std::string("canvasPZ0InOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPZ0InOut.SetGrid(1,1);
    canvasPZ0InOut.SetLogy(0);
    canvasPZ0InOut.SetFillColor(Palette::color_plot_background);
    canvasPZ0InOut.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPZ0InOut"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{Z} from previous layers/discs;#eta; #sigma_{Z} [#mum]").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

      iMomentum++;
    }
    legendPInOut->SetX1(0.64);
    legendPInOut->SetX2(0.89);
    legendPInOut->Draw("SAME");

    RootWImage& myImagePZ0InOut = myContentPlotsPZ0InOut.addImage(canvasPZ0InOut, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePZ0InOut.setComment(std::string("Detector: "+name+" - an extrapolated sigma in Z from previous layers/discs in in-out approach.").c_str());
    myImagePZ0InOut.setName(std::string("bkg_p_z0_inout_"+name).c_str());
  }

  // c) pContamination
  RootWContent& myContentPlotsPProbContamInOut = myPage.addContent("P in-out - Bkg contamination prob.: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPHitProbContamInOut) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPProbContamInOut(std::string("canvasPProbContamInOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPProbContamInOut.SetGrid(1,1);
    canvasPProbContamInOut.SetLogy(0);
    canvasPProbContamInOut.SetFillColor(Palette::color_plot_background);
    canvasPProbContamInOut.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPProbContamInOut"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: bkg contamination prob. as error ellipse extrap. from previous layers/discs;#eta; probability").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1.1);

      iMomentum++;
    }
    legendPInOut->SetX1(0.11);
    legendPInOut->SetX2(0.36);
    legendPInOut->Draw("SAME");

    RootWImage& myImagePProbContamInOut = myContentPlotsPProbContamInOut.addImage(canvasPProbContamInOut, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePProbContamInOut.setComment(std::string("Detector: "+name+" - bkg contamination prob. as error ellipse extrap. from previous layers/discs in in-out approach.").c_str());
    myImagePProbContamInOut.setName(std::string("bkg_p_pContam_inout_"+name).c_str());
  }
  //
  // Detail on p out-in studies

  // a) D0
  RootWContent& myContentPlotsPD0OutIn = myPage.addContent("P out-in - D: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPHitDProjOutIn) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPD0OutIn(std::string("canvasPD0OutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPD0OutIn.SetGrid(1,1);
    canvasPD0OutIn.SetLogy(0);
    canvasPD0OutIn.SetFillColor(Palette::color_plot_background);
    canvasPD0OutIn.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPD0OutIn"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{R-#Phi} from previous layers/discs;#eta; #sigma_{R-#Phi} [#mum]").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

      iMomentum++;
    }
    legendPOutIn->SetX1(0.64);
    legendPOutIn->SetX2(0.89);
    legendPOutIn->SetEntryLabel("p options (OutIn, full TRK):");
    legendPOutIn->Draw("SAME");

    RootWImage& myImagePD0OutIn = myContentPlotsPD0OutIn.addImage(canvasPD0OutIn, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePD0OutIn.setComment(std::string("Detector: "+name+" - an extrapolated sigma in R-Phi from previous layers/discs in out-in approach.").c_str());
    myImagePD0OutIn.setName(std::string("bkg_p_d0_outin_"+name).c_str());
  }

  // b) Z0
  RootWContent& myContentPlotsPZ0OutIn = myPage.addContent("P out-in - Z: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPHitZProjOutIn) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPZ0OutIn(std::string("canvasPZ0OutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPZ0OutIn.SetGrid(1,1);
    canvasPZ0OutIn.SetLogy(0);
    canvasPZ0OutIn.SetFillColor(Palette::color_plot_background);
    canvasPZ0OutIn.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPZ0OutIn"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{Z} from previous layers/discs;#eta; #sigma_{Z} [#mum]").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

      iMomentum++;
    }
    legendPOutIn->SetX1(0.64);
    legendPOutIn->SetX2(0.89);
    legendPOutIn->Draw("SAME");

    RootWImage& myImagePZ0OutIn = myContentPlotsPZ0OutIn.addImage(canvasPZ0OutIn, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePZ0OutIn.setComment(std::string("Detector: "+name+" - an extrapolated sigma in Z from previous layers/discs in out-in approach.").c_str());
    myImagePZ0OutIn.setName(std::string("bkg_p_z0_outin_"+name).c_str());
  }

  // c) pContamination
  RootWContent& myContentPlotsPProbContamOutIn = myPage.addContent("P out-in - Bkg contamination prob.: Details of track purity & probability of contamination by bkg", false);
  for (auto const & iterMap : m_hisPHitProbContamOutIn) {

    // Get name -> remove first 2 artificial characters used to sort correctly map
    std::string name = iterMap.first;
    name = name.erase(0,2);

    int         iMomentum = 0;

    TCanvas canvasPProbContamOutIn(std::string("canvasPProbContamOutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY);
    canvasPProbContamOutIn.SetGrid(1,1);
    canvasPProbContamOutIn.SetLogy(0);
    canvasPProbContamOutIn.SetFillColor(Palette::color_plot_background);
    canvasPProbContamOutIn.SetObjectStat(false);

    for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

      auto & profHis = iterMap.second;
      profHis[iMomentum]->SetNameTitle(std::string("canvasPProbContamOutIn"+name+any2str(iMomentum)).c_str(),std::string(name+" In-Out approach: bkg contamination prob. as error ellipse extrap. from previous layers/discs;#eta; probability").c_str());
      profHis[iMomentum]->SetLineWidth(2.);
      profHis[iMomentum]->SetMarkerStyle(21);
      profHis[iMomentum]->SetMarkerSize(1.);
      profHis[iMomentum]->SetLineColor(Palette::colorMomenta(iMomentum));
      profHis[iMomentum]->SetMarkerColor(Palette::colorMomenta(iMomentum));

      if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
      else              profHis[iMomentum]->Draw("SAME PE1");

      profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1.1);

      iMomentum++;
    }
    legendPOutIn->SetX1(0.11);
    legendPOutIn->SetX2(0.36);
    legendPOutIn->Draw("SAME");

    RootWImage& myImagePProbContamOutIn = myContentPlotsPProbContamOutIn.addImage(canvasPProbContamOutIn, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePProbContamOutIn.setComment(std::string("Detector: "+name+" - bkg contamination prob. as error ellipse extrap. from previous layers/discs in out-in approach.").c_str());
    myImagePProbContamOutIn.setName(std::string("bkg_p_pContam_outin_"+name).c_str());
  }

  m_isVisOK = true;
  return m_isVisOK;
}

//
// Check that a file can be opened
//
bool AnalyzerPatternReco::checkFile(const std::string& fileName, const std::string& filePath)
{
  fstream     file;
  std::string fullFileName(filePath+"/"+fileName);
  file.open(fullFileName);
  if (file.is_open()) {

    file.close();
    return true;
  }
  else {

    logERROR("AnalyzerPatternReco - failed opening file: " + fullFileName);
    return false;
  }
}

//
// Is starting triplet from different layers (avoid using overlapping modules in one layer)
//
bool AnalyzerPatternReco::isTripletFromDifLayers(Track& track, int iHit, bool propagOutIn) {

  std::map<std::string, bool> hitIDs;

  for (int i=0; i<=iHit; i++) {

    std::string hitID = "";

    if (!propagOutIn) hitID = track.getMeasurableOrIPHit(i)->getDetName() +
                              std::string(track.getMeasurableOrIPHit(i)->isBarrel() ? "_L_" : "_D_") +
                              any2str(track.getMeasurableOrIPHit(i)->getLayerOrDiscID());
    else              hitID = track.getRMeasurableOrIPHit(i)->getDetName() +
                              std::string(track.getRMeasurableOrIPHit(i)->isBarrel() ? "_L_" : "_D_") +
                              any2str(track.getRMeasurableOrIPHit(i)->getLayerOrDiscID());

    hitIDs[hitID] = true;
  }

  if (hitIDs.size()>=3) return true;
  else                  return false;
}
