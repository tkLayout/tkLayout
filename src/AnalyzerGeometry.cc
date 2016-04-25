/*
 * AnalyzerGeometry.cc
 *
 *  Created on: 20. 4. 2016
 *      Author: drasal
 */
#include <AnalyzerGeometry.h>

// Include files
#include <Barrel.h>
#include <DetectorModule.h>
#include <Disk.h>
#include <global_constants.h>
#include <Endcap.h>
#include <Layer.h>
#include <Math/Vector3D.h>
#include <PlotDrawer.h>
#include <rootweb.hh>
#include <ostream>
#include <TagMaker.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TH2D.h>
#include <TH2I.h>
#include <TPad.h>
#include <TProfile.h>
#include <Tracker.h>
#include <TRandom3.h>
#include <SimParms.h>

// Namespaces
using namespace insur;

//
// AnalyzerGeometry constructor
//
AnalyzerGeometry::AnalyzerGeometry(std::vector<Tracker*> trackers, int nTracks) : AnalyzerModule(trackers),
 m_nTracks(nTracks),
 m_layerNames(trackers),
 m_etaSpan(geom_max_eta_coverage - geom_max_eta_coverage),
 m_etaMin(-1*geom_max_eta_coverage),
 m_etaMax(+1*geom_max_eta_coverage)
{};

//
// AnalyzerGeometry init method
//
bool AnalyzerGeometry::init()
{
  if (m_nTracks<=0 || m_trackers.size()==0) {

    std::ostringstream message;
    message << "AnalyzerGeometry::init(): Number of tracks zero or negative: " << m_nTracks << " or zero number of trackers to be initialized";
    logERROR(message.str());
    return false;
  }
  else {

    // Number of tracks needs to be expressed as int * int for proper binning
    int sqrtNTracks  = int(sqrt(m_nTracks));
    int nBins        = int(sqrtNTracks/2.);
    int nBinsProfile = 100;
    m_nTracks        = sqrtNTracks*sqrtNTracks;

    // Compute shooting intervals to analyze geometry
    double etaMin = -1*geom_max_eta_coverage;
    double etaMax = +1*geom_max_eta_coverage;

    float  safeMargin = c_etaSafetyMargin;
    m_etaSpan         = (etaMax - etaMin)*(1. + safeMargin);
    m_etaMax          = etaMax * (1 + safeMargin);

    // Prepare histograms
    for (auto iTracker : m_trackers) {

      std::string trkName = iTracker->myid();

      if (m_hitMapPhiEta.find(trkName)==m_hitMapPhiEta.end()) m_hitMapPhiEta[trkName] = TH2D();
      m_hitMapPhiEta[trkName].Reset();
      m_hitMapPhiEta[trkName].SetNameTitle("HitMapPhiEta", "Number of hits;#phi;#eta");
      m_hitMapPhiEta[trkName].SetBins(nBins, -M_PI, M_PI, nBins, -m_etaMax, +m_etaMax);

      if (m_trackMapPhiEta.find(trkName)==m_trackMapPhiEta.end()) m_trackMapPhiEta[trkName] = TH2I();
      m_trackMapPhiEta[trkName].Reset();
      m_trackMapPhiEta[trkName].SetNameTitle("TrackMapPhiEta", "Number of tracks;#phi;#eta");
      m_trackMapPhiEta[trkName].SetBins(nBins, -M_PI, M_PI, nBins, -m_etaMax, +m_etaMax);

      if (m_moduleHitEtaProfile.find(trkName)==m_moduleHitEtaProfile.end()) m_moduleHitEtaProfile[trkName] = TProfile();
      m_moduleHitEtaProfile[trkName].Reset();
      m_moduleHitEtaProfile[trkName].SetNameTitle("ModuleHitEtaProfile","Number of modules with at least one hit;#eta;Number of hit modules");
      m_moduleHitEtaProfile[trkName].SetBins(nBinsProfile, 0, +m_etaMax);
      m_moduleHitEtaProfile[trkName].SetMarkerStyle(8);
      m_moduleHitEtaProfile[trkName].SetMarkerColor(1);
      m_moduleHitEtaProfile[trkName].SetMarkerSize(1.5);
      m_moduleHitEtaProfile[trkName].SetMinimum(0);

      if (m_sensorHitEtaProfile.find(trkName)==m_sensorHitEtaProfile.end()) m_sensorHitEtaProfile[trkName] = TProfile();
      m_sensorHitEtaProfile[trkName].Reset();
      m_sensorHitEtaProfile[trkName].SetNameTitle("ModuleHitEtaProfile","Number of hit sensors;#eta;Number of hit sensors");
      m_sensorHitEtaProfile[trkName].SetBins(nBinsProfile, 0, +m_etaMax);
      m_sensorHitEtaProfile[trkName].SetMarkerStyle(8);
      m_sensorHitEtaProfile[trkName].SetMarkerColor(1);
      m_sensorHitEtaProfile[trkName].SetMarkerSize(1.5);
      m_sensorHitEtaProfile[trkName].SetMinimum(0);

      if (m_stubHitEtaProfile.find(trkName)==m_stubHitEtaProfile.end()) m_stubHitEtaProfile[trkName] = TProfile();
      m_stubHitEtaProfile[trkName].Reset();
      m_stubHitEtaProfile[trkName].SetNameTitle("ModuleHitEtaProfile","Number of modules with a stub;#eta;Number of hit stubs");
      m_stubHitEtaProfile[trkName].SetBins(nBinsProfile, 0, +m_etaMax);
      m_stubHitEtaProfile[trkName].SetMarkerStyle(8);
      m_stubHitEtaProfile[trkName].SetMarkerColor(1);
      m_stubHitEtaProfile[trkName].SetMarkerSize(1.5);
      m_stubHitEtaProfile[trkName].SetMinimum(0);

      // Prepare eta coverage profiles for individual barrel layers & end-cap disks
      if (m_layerEtaCoverProfile.find(trkName)!=m_layerEtaCoverProfile.end()) m_layerEtaCoverProfile[trkName].clear();
      if (m_layerStubEtaCoverProfile.find(trkName)!=m_layerStubEtaCoverProfile.end()) m_layerStubEtaCoverProfile[trkName].clear();

      for (auto layerName : m_layerNames.m_data) {

        // Set names, binning, etc.
        TProfile layerProfile;
        layerProfile.SetNameTitle("LayerEtaCoverage", std::string("Eta coverage of layer: "+layerName+";#eta;Coverage").c_str());
        layerProfile.SetBins(2*nBinsProfile, -m_etaMax, +m_etaMax);
        if (m_layerEtaCoverProfile.find(trkName)!=m_layerEtaCoverProfile.end()) m_layerEtaCoverProfile[trkName].clear();
        m_layerEtaCoverProfile[trkName][layerName] = layerProfile;
        m_layerEtaCoverProfile[trkName][layerName].SetMinimum(0);

        TProfile stubProfile;
        stubProfile.SetNameTitle("LayerStubEtaCoverage", std::string("Eta coverage of layer stub: "+layerName+";#eta;Coverage").c_str());
        stubProfile.SetBins(2*nBinsProfile, -m_etaMax, +m_etaMax);
        m_layerStubEtaCoverProfile[trkName][layerName] = stubProfile;
        m_layerStubEtaCoverProfile[trkName][layerName].SetMinimum(0);
      }
    } // For trackers
    return true;
  }
}

//
// AnalyzerGeometry analysis method
//
bool AnalyzerGeometry::analyze()
{
  // Initialize random number generator, counters and histograms
  TRandom3 myDice;
  myDice.SetSeed(random_seed);

  // Print
  std::cout << std::endl;

  // Go through all trackers
  for (auto iTracker : m_trackers) {

    std::string trkName =  iTracker->myid();

    // Print
    std::cout << " " << iTracker->myid() << " tracker ..."<< std::endl;

    // Counters
    std::map <std::string, int> moduleTypeCount;  // counts hit per module type -> if any of the sensors (or both) are hit, it counts 1
    std::map <std::string, int> sensorTypeCount;  // counts hit per sensor on module type -> if one sensor is hit, it counts 1, if both sensors are hit, it counts 2
    std::map <std::string, int> stubTypeCount;    // counts stubs per module type -> if both sensors are hit and a stub is formed, it counts 1

    int color = -1;

    // Analyze geometry layout
    for (auto iTrack=0; iTrack<m_nTracks; iTrack++) {

      // Generate a straight track & track origin
      double phi   = myDice.Rndm() * 2 * M_PI;
      double eta   = myDice.Rndm() * m_etaSpan + m_etaMin;
      double theta = 2*atan(exp(-1*eta));

      ROOT::Math::XYZVector direction = ROOT::Math::XYZVector(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));

      double zError = SimParms::getInstance()->zErrorCollider();
      double zPos   = myDice.Rndm()*2 - 1;

      ROOT::Math::XYZVector origin = ROOT::Math::XYZVector(0, 0, zPos);

      // Check whether generated track hits a module -> collect the list of hit modules
      std::vector<std::pair<Module*, HitType>> hitModules;

      static const double zErrorSafetyMargin = 5. ; // track origin shift in units of zError to compute boundaries

      // Go through all modules and save only those that were hit
      for (auto& module : iTracker->modules()) {

        // CouldHit -> save CPU time by calculating a cross-section of hit and given module within track phi and eta
        // (assuming origin error within +-5 sigma). Beware that module phi extends from -pi to +3*pi to avoid troubles
        // with -pi & +pi cross-line
        if (module->couldHit(direction, zError*zErrorSafetyMargin)) {

          // Get hit
          auto hit = module->checkTrackHits(origin, direction);
          if (hit.second != HitType::NONE) hitModules.push_back(std::make_pair(module,hit.second));
        }
      }

      // Reset counters
      int numStubs = 0;
      int numHits  = 0;

      for (auto& module : hitModules) {

        // Inner
        if (module.second & HitType::INNER) {
          numHits++;
        }
        // Outer
        if (module.second & HitType::OUTER) {
          numHits++;
        }
        // Stubs
        if (module.second == HitType::STUB) {
          numStubs++;
        }
        // Set color only for the first run
        if (color==-1) color = module.first->plotColor();
      }

      // Fill plots
      m_hitMapPhiEta[trkName].Fill(direction.Phi(), direction.Eta(), hitModules.size());
      m_trackMapPhiEta[trkName].Fill(direction.Phi(), direction.Eta());

      m_moduleHitEtaProfile[trkName].Fill(fabs(direction.Eta()), hitModules.size());
      m_sensorHitEtaProfile[trkName].Fill(fabs(direction.Eta()), numHits);
      m_stubHitEtaProfile[trkName].Fill(fabs(direction.Eta()), numStubs);

      for (auto layerName : m_layerNames.m_data) {

        int layerHit  = 0;
        int layerStub = 0;
        for (auto module : hitModules) {

          // Unique reference of module
          UniRef uniref = module.first->uniRef();
          if (layerName == (uniref.cnt + " " + any2str(uniref.layer))) {

            layerHit = 1;
            if (module.second == HitType::STUB) layerStub=1;

            // In maximum one hit needed to calculate coverage (i.e efficiency)
            if (layerHit && layerStub) break;
          }
        }

        // Fill value
        m_layerEtaCoverProfile[trkName][layerName].Fill(direction.Eta(), layerHit);
        m_layerStubEtaCoverProfile[trkName][layerName].Fill(direction.Eta(), layerStub);
      }

    } // For nTracks

    // Set fill color
    m_moduleHitEtaProfile[trkName].SetMarkerColor(Palette::color(color));
    m_sensorHitEtaProfile[trkName].SetMarkerColor(Palette::color(color));
    m_stubHitEtaProfile[trkName].SetMarkerColor(Palette::color(color));

    // Normalize histograms
    for (auto xBin=0; xBin<=m_hitMapPhiEta[trkName].GetNbinsX()+1; xBin++) {
      for (auto yBin=0; yBin<=m_hitMapPhiEta[trkName].GetNbinsY()+1; yBin++) {

        int trackCount = m_trackMapPhiEta[trkName].GetBinContent(xBin, yBin);
        if (trackCount>0) {

          double hitCount=m_hitMapPhiEta[trkName].GetBinContent(xBin, yBin);
          m_hitMapPhiEta[trkName].SetBinContent(xBin, yBin, hitCount/trackCount);
        }
      }
    }

  } // For trackers
  return true;
}

//
// AnalyzerGeometry visualization method
//
bool AnalyzerGeometry::visualize(RootWSite& webSite)
{
  // Go through all trackers & prepare web content
  int webPriority         = 99;
  RootWPage*    myPage    = nullptr;
  RootWContent* myContent = nullptr;

  // List of tracker RZ, BRL XY, ECAP XY canvases
  std::map<std::string, RootWImage*> trackerRZImg;
  std::map<std::string, RootWImage*> trackerXYImgBRL;
  std::map<std::string, RootWImage*> trackerXYImgEC;
  PlotDrawer<YZ, Type>               fullRZDrawer;

  for (auto iTracker : m_trackers) {

    // Tracker name
    std::string trkName =  iTracker->myid();

    // Create dedicated web-page & its content
    std::string        pageTitle  = "Geometry";
    if (trkName != "") pageTitle +=" (" + trkName + ")";

    myPage = new RootWPage(pageTitle);

    // Page address
    std::string pageAddress = "index"+trkName+".html";
    myPage->setAddress(pageAddress);

    webSite.addPage(myPage, webPriority);
    webPriority--;

    //
    // Layer & disk layout
    myContent = new RootWContent("Layer & disk layout");
    myPage->addContent(myContent);

    // Fill in all information directly from Tracker geometry object
    LayerDiskSummaryVisitor geometryVisitor;
    geometryVisitor.preVisit();
    iTracker->accept(geometryVisitor);
    geometryVisitor.postVisit();

    // Print out layer & disk table
    myContent->addItem(geometryVisitor.m_layerTable);
    myContent->addItem(geometryVisitor.m_diskTable);
    myContent->addItem(geometryVisitor.m_ringTable);

    //
    // Modules
    myContent = new RootWContent("Modules info", false);
    myPage->addContent(myContent);

    myContent->addItem(geometryVisitor.m_moduleTable);

    //
    // Plots
    myContent = new RootWContent("Layout plots");
    myPage->addContent(myContent);

    RootWImage* myImage;

    // Beam pipe parameters
    double bpRadius = SimParms::getInstance()->bpRadius();
    double bpThick  = SimParms::getInstance()->bpThickness();

    // R-Z view - avoid drawing canvas if barrel and endcap missing
    PlotDrawer<YZ, Type> rzDrawer;
    bool foundModules = rzDrawer.addModulesType(iTracker->modules().begin(), iTracker->modules().end(), BARREL | ENDCAP);
    if (foundModules) {

      TCanvas RZCanvas("RZCanvas", "RZView Canvas", insur::vis_max_canvas_sizeX, insur::vis_min_canvas_sizeY );
      RZCanvas.cd();
      rzDrawer.drawFrame<SummaryFrameStyle>(RZCanvas, iTracker->isPixelType());
      rzDrawer.drawModules<ContourStyle>(RZCanvas);
      fullRZDrawer.addModulesType(iTracker->modules().begin(), iTracker->modules().end());

      // Beam-pipe
      TPolyLine beamPipeRZ;
      beamPipeRZ.SetPoint(0, 0                                 , bpRadius + bpThick/2.);
      beamPipeRZ.SetPoint(1, iTracker->maxZ()*vis_safety_factor, bpRadius + bpThick/2.);
      beamPipeRZ.SetLineColor(14);
      beamPipeRZ.SetLineWidth(2);
      if (iTracker->isPixelType()) {
        logINFO("Drawing beam pipe to RZ canvas");
        beamPipeRZ.Draw("same");
      }

      myImage = new RootWImage(RZCanvas, RZCanvas.GetWindowWidth(), RZCanvas.GetWindowHeight() );
      myImage->setComment("RZ plot of the tracker");
      trackerRZImg[trkName] = myImage;
      myContent->addItem(myImage);
    }

    // Projections - avoid drawing canvas if barrel or endcap missing
    PlotDrawer<XY, Type> xyBarrelDrawer;
    foundModules = xyBarrelDrawer.addModulesType(iTracker->modules().begin(), iTracker->modules().end(), BARREL);
    if (foundModules) {

      TCanvas XYCanvasBRL("XYCanvasBRL", "XYView Barrel", insur::vis_min_canvas_sizeX, insur::vis_min_canvas_sizeY );
      XYCanvasBRL.cd();
      xyBarrelDrawer.drawFrame<SummaryFrameStyle>(XYCanvasBRL);
      xyBarrelDrawer.drawModules<ContourStyle>(XYCanvasBRL);

      // Beam-pipe
      TEllipse beamPipeXY(0, 0, bpRadius + bpThick/2.);
      beamPipeXY.SetFillColor(18); // "grey18"
      beamPipeXY.SetFillStyle(1001);
      if (iTracker->isPixelType()) beamPipeXY.Draw("same");

      myImage = new RootWImage(XYCanvasBRL, XYCanvasBRL.GetWindowWidth(), XYCanvasBRL.GetWindowHeight() );
      myImage->setComment("XY cross section of barrel modules");
      trackerXYImgBRL[trkName] = myImage;
      myContent->addItem(myImage);
    }

    PlotDrawer<XY, Type> xyEndcapDrawer;
    foundModules = xyEndcapDrawer.addModulesType(iTracker->modules().begin(), iTracker->modules().end(), ENDCAP);
    if (foundModules) {

      TCanvas XYCanvasEC("XYCanvasEC", "XYView Endcap", insur::vis_min_canvas_sizeX, insur::vis_min_canvas_sizeY );
      XYCanvasEC.cd();
      xyEndcapDrawer.drawFrame<SummaryFrameStyle>(XYCanvasEC);
      xyEndcapDrawer.drawModules<ContourStyle>(XYCanvasEC);

      // Beam-pipe
      TEllipse beamPipeXY(0, 0, bpRadius + bpThick/2.);
      beamPipeXY.SetFillColor(18); // "grey18"
      beamPipeXY.SetFillStyle(1001);
      if (iTracker->isPixelType()) beamPipeXY.Draw("same");

      //trackerXYProjEC[trkName] = (TCanvas*)XYCanvasEC->Clone();

      myImage = new RootWImage(XYCanvasEC, XYCanvasEC.GetWindowWidth(), XYCanvasEC.GetWindowHeight() );
      myImage->setComment("XY projection of endcap(s) modules");
      trackerXYImgEC[trkName] = myImage;
      myContent->addItem(myImage);
    }

    //
    // Hit eta profile plots

    // Sensor hit map
    if (m_sensorHitEtaProfile.find(trkName)!=m_sensorHitEtaProfile.end()) {

      TCanvas etaCanvasHitSensors("NumberOfHitSensorEtaProfile", "Number of hit sensors versus #eta");
      etaCanvasHitSensors.SetFillColor(kWhite);
      etaCanvasHitSensors.SetBorderMode(0);
      etaCanvasHitSensors.SetBorderSize(0);
      m_sensorHitEtaProfile[trkName].Draw();

      myImage = new RootWImage(etaCanvasHitSensors, etaCanvasHitSensors.GetWindowWidth(), etaCanvasHitSensors.GetWindowHeight());
      myImage->setComment("Number of hit sensors versus "+insur::web_etaLetter);
      myContent->addItem(myImage);
    }

    // Module hit map
    if (m_moduleHitEtaProfile.find(trkName)!=m_moduleHitEtaProfile.end()) {

      TCanvas etaCanvasHitModules("NumberOfHitModulesEtaProfile", "Number of hit modules versus #eta");
      etaCanvasHitModules.SetFillColor(kWhite);
      etaCanvasHitModules.SetBorderMode(0);
      etaCanvasHitModules.SetBorderSize(0);
      m_moduleHitEtaProfile[trkName].Draw();

      myImage = new RootWImage(etaCanvasHitModules, etaCanvasHitModules.GetWindowWidth(), etaCanvasHitModules.GetWindowHeight());
      myImage->setComment("Number of hit modules versus "+insur::web_etaLetter);
      myContent->addItem(myImage);
    }

    // 2D hit map
    if (m_hitMapPhiEta.find(trkName)!=m_hitMapPhiEta.end()) {

      TCanvas etaCanvasHitMap("HitMap2D", "Hit map - #eta versus #phi");
      etaCanvasHitMap.SetFillColor(kWhite);
      etaCanvasHitMap.SetBorderMode(0);
      etaCanvasHitMap.SetBorderSize(0);
      m_hitMapPhiEta[trkName].Draw("colz");

      myImage = new RootWImage(etaCanvasHitMap, etaCanvasHitMap.GetWindowWidth(), etaCanvasHitMap.GetWindowHeight());
      myImage->setComment("Hit coverage in "+insur::web_etaLetter+", "+insur::web_phiLetter);
      myContent->addItem(myImage);
    }

    // Layer coverage plots
    if ((m_layerEtaCoverProfile.find(trkName)!=m_layerEtaCoverProfile.end()) && m_layerEtaCoverProfile[trkName].size()>0) {

      RootWContent* myContent = new RootWContent("Layer coverage (Hits)", false);
      myPage->addContent(myContent);

      int layerCount = 0;
      for (auto& iter : m_layerEtaCoverProfile[trkName]) {

        TProfile& plot = iter.second;
        layerCount++;

        TCanvas coverageCanvas(Form("LayerCoverage%s%s", iter.first.c_str(), "Hits"), "Layer eta coverage (Hits)", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
        coverageCanvas.cd();

        std::string upperName = std::string(coverageCanvas.GetName())+"_upper";
        TPad upperPad(upperName.c_str(), "upperPad", 0, 0.4, 1, 1);
        std::string lowerName = std::string(coverageCanvas.GetName())+"_lower";
        TPad lowerPad(lowerName.c_str(), "lowerPad", 0, 0, 1, 0.4);
        upperPad.Draw();
        lowerPad.Draw();

        plot.SetMinimum(0);
        plot.SetMaximum(1.05);
        plot.SetMarkerColor(Palette::color(1));
        plot.SetLineColor(Palette::color(1));
        plot.SetMarkerStyle(1);

        TProfile* zoomedPlot = (TProfile*) plot.Clone(std::string(std::string(plot.GetName())+"zoomed").c_str());
        zoomedPlot->SetMinimum(0.9);
        zoomedPlot->SetMaximum(1.01);
        zoomedPlot->SetTitle("");
        upperPad.cd();
        plot.Draw();
        lowerPad.cd();
        zoomedPlot->Draw();

        RootWImage* myImage = new RootWImage(coverageCanvas, coverageCanvas.GetWindowWidth(), coverageCanvas.GetWindowHeight());
        myImage->setComment("Layer coverage in "+insur::web_etaLetter+" for hits (multiple occurrence in a layer counted just once)");
        myContent->addItem(myImage);
      }
    }
  } // Trackers

  //
  // Create introductory geometry web page
  myPage = new RootWPage("Overview");

  // Page address
  std::string pageAddress = "index.html";
  myPage->setAddress(pageAddress);
  webSite.addPage(myPage, 100);

  // Full tracker layout
  myContent = new RootWContent("Tracker layout: ");
  myPage->addContent(myContent);

  TCanvas fullRZCanvas("RZCanvasFull", "RZView Canvas Full", insur::vis_max_canvas_sizeX, insur::vis_min_canvas_sizeY );
  fullRZDrawer.drawFrame<SummaryFrameStyle>(fullRZCanvas);
  fullRZDrawer.drawModules<ContourStyle>(fullRZCanvas);

  RootWImage* myImage = new RootWImage(fullRZCanvas, fullRZCanvas.GetWindowWidth(), fullRZCanvas.GetWindowHeight());
  myContent->addItem(myImage);

  // All subdetectors
  for (auto iter : trackerRZImg) {

    std::string trkName = iter.first;

    myContent = new RootWContent("Geometry(" + trkName + "): ");
    myPage->addContent(myContent);
    myContent->addItem(iter.second);
    if (trackerXYImgBRL.find(trkName)!=trackerXYImgBRL.end()) myContent->addItem(trackerXYImgBRL[trkName]);
    if (trackerXYImgEC.find(trkName) !=trackerXYImgEC.end())  myContent->addItem(trackerXYImgEC[trkName]);
  }

  return true;


//  //
//  // Inactive surfaces
//  double inactiveSurfacesTotalMass;
//  if (inactive) {
//    std::vector<InactiveElement>& inactiveBarrelServices = inactive->getBarrelServices();
//    std::vector<InactiveElement>& inactiveEndcapServices = inactive->getEndcapServices();
//    std::vector<InactiveElement>& inactiveSupports       = inactive->getSupports();
//    std::vector<InactiveElement> allInactives;
//
//    allInactives.reserve( inactiveBarrelServices.size() + inactiveEndcapServices.size() + inactiveSupports.size() );
//    allInactives.insert(allInactives.end(), inactiveBarrelServices.begin(), inactiveBarrelServices.end() );
//    allInactives.insert(allInactives.end(), inactiveEndcapServices.begin(), inactiveEndcapServices.end() );
//    allInactives.insert(allInactives.end(), inactiveSupports.begin(), inactiveSupports.end() );
//
//    inactiveSurfacesTotalMass = 0;
//    for (const auto elem : allInactives ) {
//      if (elem.getTotalMass()>0) inactiveSurfacesTotalMass += elem.getTotalMass();
//    }
//  }



}

//
// LayerNameVisitor constructor method
//
LayerNameVisitor::LayerNameVisitor(const std::vector<Tracker*>& trackers)
{
  for (auto tracker : trackers) tracker->accept(*this);
}

//
// LayerNameVisitor - visit barrel
//
void LayerNameVisitor::visit(const Barrel& b) { m_id = b.myid(); }

//
// LayerNameVisitor - visit endcap
//
void LayerNameVisitor::visit(const Endcap& e) { m_id = e.myid(); }

//
// LayerNameVisitor - visit layer
//
void LayerNameVisitor::visit(const Layer& l)  { m_data.insert(m_id + " " + any2str(l.myid())); }

//
// LayerNameVisitor - visit disk
//
void LayerNameVisitor::visit(const Disk& d)   { m_data.insert(m_id + " " + any2str(d.myid())); }

//
// LayerDiskSummaryVisitor preVisit method - initialize
//
void LayerDiskSummaryVisitor::preVisit() {

  // Initialize
  m_nBarrelLayers      = 0;
  m_nDisks             = 0;
  m_totalBarrelModules = 0;
  m_totalEndcapModules = 0;

  m_totalArea       = 0;
  m_totalModules    = 0;
  m_totalSensors    = 0;
  m_totalChannels   = 0;
  m_totalSensorPower= 0;

  m_layerTable = new RootWTable();
  m_diskTable  = new RootWTable();
  m_ringTable  = new RootWTable();
  m_moduleTable= new RootWTable();

  m_layerTable->setContent(0, 0, "Layer no                    : ");
  m_layerTable->setContent(1, 0, "Radius [mm]                 : ");
  m_layerTable->setContent(2, 0, "Z-min [mm]                  : ");
  m_layerTable->setContent(3, 0, "Z-max [mm]                  : ");
  m_layerTable->setContent(4, 0, "Number of rods              : ");
  m_layerTable->setContent(5, 0, "Number of modules per rod   : ");
  m_layerTable->setContent(6, 0, "Number of modules           : ");

  m_diskTable->setContent(0, 0, "Disk no                      : ");
  m_diskTable->setContent(1, 0, "Radius-min [mm]              : ");
  m_diskTable->setContent(2, 0, "Radius-max [mm]              : ");
  m_diskTable->setContent(3, 0, "Average Z pos. [mm]          : ");
  m_diskTable->setContent(4, 0, "Number of rings              : ");
  m_diskTable->setContent(5, 0, "Number of modules per disk   : ");

  m_ringTable->setContent(0, 0, "Ring no                      : ");
  m_ringTable->setContent(1, 0, "R-min [mm]                   : ");
  m_ringTable->setContent(2, 0, "R-max [mm]                   : ");
  m_ringTable->setContent(3, 0, "Number of modules per ring   : ");

  m_moduleTable->setContent( 0, 0, "Module in:                                   ");
  m_moduleTable->setContent( 1, 0, "Position:                                    ");
  m_moduleTable->setContent( 2, 0, "Type:                                        ");
  m_moduleTable->setContent( 3, 0, "Sensor area [mm"+insur::web_superStart+"2"+insur::web_superEnd+"]: ");
  m_moduleTable->setContent( 4, 0, "Total area [m"+insur::web_superStart+"2"+insur::web_superEnd+"]:   ");
//  m_moduleTable->setContent( 5, 0, "Service Weight [kg]:                         ");
//  m_moduleTable->setContent( 6, 0, "Total Weight [kg]:                           ");
  m_moduleTable->setContent( 5, 0, "Number of modules:                           ");
  m_moduleTable->setContent( 6, 0, "Number of sensors:                           ");
  m_moduleTable->setContent( 7, 0, "Number of channels (M):                      ");
//  m_moduleTable->setContent(10, 0, "Number of channels (R-Phi 1.side):           ");
//  m_moduleTable->setContent(11, 0, "Number of channels (R-Phi 2.side):           ");
//  m_moduleTable->setContent(12, 0, "Number of channels (Z 1.side):               ");
//  m_moduleTable->setContent(13, 0, "Number of channels (Z 2.side):               ");
  m_moduleTable->setContent( 8, 0, "Min-Max R-Phi resolution ("+insur::web_muLetter+"m):    ");
  m_moduleTable->setContent( 9, 0, "Min-Max Z resolution ("+insur::web_muLetter+"m):        ");
}

//
// LayerDiskSummaryVisitor visit layer
//
void LayerDiskSummaryVisitor::visit(const Layer& l) {

  if (l.maxZ() < 0.) return;

  // Update layer counter
  ++m_nBarrelLayers;

  // Update module counter
  int nModules        = l.totalModules();
  m_totalBarrelModules += nModules;

  // Set table
  m_layerTable->setContent(0, m_nBarrelLayers, m_nBarrelLayers);
  m_layerTable->setContent(1, m_nBarrelLayers, l.placeRadius(), c_coordPrecision);
  m_layerTable->setContent(2, m_nBarrelLayers, l.minZ(), c_coordPrecision);
  m_layerTable->setContent(3, m_nBarrelLayers, l.maxZ(), c_coordPrecision);
  m_layerTable->setContent(4, m_nBarrelLayers, l.numRods());
  m_layerTable->setContent(5, m_nBarrelLayers, l.numModulesPerRod());
  m_layerTable->setContent(6, m_nBarrelLayers, l.totalModules());
}

//
// LayerDiskSummaryVisitor visit disk
//
void LayerDiskSummaryVisitor::visit(const Disk& d) {

  if (d.averageZ() < 0.) return;

  // Update disk counter
  if (m_nDisks==0) m_ringNModules.resize(d.numRings());
  ++m_nDisks;

  // Update module counter
  int nModules = d.totalModules();
  m_totalEndcapModules += nModules;

  // Set table
  m_diskTable->setContent(0, m_nDisks, d.myid());
  m_diskTable->setContent(1, m_nDisks, d.minR(),     c_coordPrecision);
  m_diskTable->setContent(2, m_nDisks, d.maxR(),     c_coordPrecision);
  m_diskTable->setContent(3, m_nDisks, d.averageZ(), c_coordPrecision);
  m_diskTable->setContent(4, m_nDisks, d.numRings());
  m_diskTable->setContent(5, m_nDisks, nModules);
}

//
// LayerDiskSummaryVisitor visit ring
//
void LayerDiskSummaryVisitor::visit(const Ring& r) {

  if (r.averageZ() < 0.) return;

  m_ringNModules[r.myid()-1] = r.numModules();
}

//
// LayerDiskSummaryVisitor visit module (barrel + endcap)
//
void LayerDiskSummaryVisitor::visit(const Module& m) {

  // Get unique sensor tag
  TagMaker moduleTagMaker(m);
  std::string tag = moduleTagMaker.sensorGeoTag;

  // Relate tags to positions
  m_moduleTagToPositionsMap[tag].insert(moduleTagMaker.posTag);

  // Get number of minimum bias events
  int nMinBiasEvents = SimParms::getInstance()->numMinBiasEvents();

  //Update sensor counters
  m_moduleCount[tag]++;
  m_moduleChannels[tag]          += m.totalChannels();
  m_moduleMaxStripOccupancy[tag]  = MAX(m.stripOccupancyPerEvent()*nMinBiasEvents, m_moduleMaxStripOccupancy[tag]);
  m_moduleAvgStripOccupancy[tag] += m.stripOccupancyPerEvent()*nMinBiasEvents;
  m_moduleMaxHitOccupancy[tag]    = MAX(m.hitOccupancyPerEvent()*nMinBiasEvents  , m_moduleMaxHitOccupancy[tag]);
  m_moduleAvgHitOccupancy[tag]   += m.hitOccupancyPerEvent()*nMinBiasEvents;
  m_moduleMaxRphiResolution[tag]  = MAX(m.resolutionLocalX(), m_moduleMaxRphiResolution[tag]);
  m_moduleMinRphiResolution[tag]  = MIN(m.resolutionLocalX(), m_moduleMaxRphiResolution[tag]);
  m_moduleAvgRphiResolution[tag] += m.resolutionLocalX();
  m_moduleMaxZResolution[tag]     = MAX(m.resolutionLocalY(), m_moduleMaxZResolution[tag]);
  m_moduleMinZResolution[tag]     = MIN(m.resolutionLocalY(), m_moduleMaxZResolution[tag]);
  m_moduleAvgZResolution[tag]    += m.resolutionLocalY();
  m_moduleMaxPower[tag]           = MAX(m.irradiationPower(), m_moduleMaxPower[tag]);
  m_moduleAvgPower[tag]          += m.irradiationPower();

  // Update counters
  m_totalModules++;
  m_totalSensors     += m.numSensors();
  m_totalChannels    += m.totalChannels();
  m_totalArea        += m.area()*m.numSensors();
  m_totalSensorPower += m.irradiationPower();

  // Get all sensor types for given module - if not find -> assign
  if (m_modulePtrMap.find(tag)==m_modulePtrMap.end()) m_modulePtrMap[tag] = &m;
}

//
// LayerDiskSummaryVisitor visit end-cap module
//
void LayerDiskSummaryVisitor::visit(const EndcapModule& m) {

  // All disks symmetric - visit only the first one
  if (m.disk() != 1 && m.side() != 1) return;

  // Get all end-cap module types for given ring - if not find -> assign
  if (m_ringModuleMap.find(m.ring())==m_ringModuleMap.end()) m_ringModuleMap[m.ring()] = &m;
}

//
// LayerDiskSummaryVisitor postvisit method - finalize
//
void LayerDiskSummaryVisitor::postVisit() {

  // Normalize
  for (auto& iter : m_moduleCount) {

    auto tag   = iter.first;
    auto count = iter.second;

    if (m_moduleCount[tag]>0) {

      int normFactor = m_moduleCount[tag];

      m_moduleAvgHitOccupancy[tag]   /= normFactor;
      m_moduleAvgStripOccupancy[tag] /= normFactor;
      m_moduleAvgRphiResolution[tag] /= normFactor;
      m_moduleAvgZResolution[tag]    /= normFactor;
      m_moduleAvgPower[tag]          /= normFactor;
    }
  }

  // Fill module table
  int iType = 0;

  double trkTotArea        = 0;
  double trkTotNumModules  = 0;
  double trkTotNumChannels = 0;

  for (auto& iter : m_modulePtrMap) {

    iType++;

    // Module
    auto moduleType = iter.first;
    auto module     = iter.second;

    // Position
    std::ostringstream position("");
    if (dynamic_cast<const BarrelModule*>(module)) position << std::dec << "Barrel";
    if (dynamic_cast<const EndcapModule*>(module)) position << std::dec << "Endcap";

    // Tag
    std::ostringstream tag("");
    tag << insur::web_smallStart;
    for (auto iter : m_moduleTagToPositionsMap[moduleType]) tag << iter << "<br/> ";
    tag << insur::web_smallEnd;

    // Type
    std::ostringstream type(module->moduleType());

    // Thickness
    std::ostringstream thickness("");
    thickness << module->dsDistance();

    // Area
    std::ostringstream sensorArea("");
    std::ostringstream totalArea("");
    sensorArea << std::dec << std::fixed << std::setprecision(c_areaPrecision) << module->area();
    totalArea  << std::dec << std::fixed << std::setprecision(c_areaPrecision) << module->area()*module->numSensors()*m_moduleCount[moduleType]* 1e-6;
    trkTotArea += module->area()*module->numSensors()*m_moduleCount[moduleType]* 1e-6;

    // Occupancy (in percent)
    std::ostringstream stripOccupancy("");
    std::ostringstream hitOccupancy("");
    stripOccupancy << std::dec << std::fixed << std::setprecision(c_occupancyPrecision) <<  m_moduleMaxStripOccupancy[moduleType]*100 << "/" << m_moduleAvgStripOccupancy[moduleType]*100; // PercenmoduleTypee
    hitOccupancy   << std::dec << std::fixed << std::setprecision(c_occupancyPrecision) <<  m_moduleMaxHitOccupancy[moduleType]*100   << "/" << m_moduleAvgHitOccupancy[moduleType]*100; // PercenmoduleTypee

    // Resolution (in um)
    std::ostringstream rphiResolution("");
    std::ostringstream zResolution("");
    rphiResolution << std::dec << std::fixed << std::setprecision(c_resolutionPrecision) << m_moduleMinRphiResolution[moduleType]*1000<< "-" << m_moduleMaxRphiResolution[moduleType]*1000; // mm -> um
    zResolution    << std::dec << std::fixed << std::setprecision(c_resolutionPrecision) << m_moduleMinZResolution[moduleType]*1000   << "-" << m_moduleMaxZResolution[moduleType]*1000; // mm -> um

    // Number of modules
    std::ostringstream numberMod("");
    numberMod << std::dec << m_moduleCount[moduleType];
    trkTotNumModules += m_moduleCount[moduleType];

    // Number of sensors
    std::ostringstream numberSens("");
    numberSens << std::dec << m_moduleCount[moduleType] * module->numSensors();

    // Number of channels (in millions)
    std::ostringstream numberChan("");
    numberChan << std::fixed << std::setprecision(c_channelPrecision) << m_moduleChannels[moduleType]/ 1e6 ;
    trkTotNumChannels += m_moduleChannels[moduleType]/ 1e6 ;

    // Fill
    m_moduleTable->setContent(0, iType, position.str());
    m_moduleTable->setContent(1, iType, tag.str());
    m_moduleTable->setContent(2, iType, type.str());
    m_moduleTable->setContent(3, iType, sensorArea.str());
    m_moduleTable->setContent(4, iType, totalArea.str());
    m_moduleTable->setContent(5, iType, numberMod.str());
    m_moduleTable->setContent(6, iType, numberSens.str());
    m_moduleTable->setContent(7, iType, numberChan.str());
    m_moduleTable->setContent(8, iType, rphiResolution.str());
    m_moduleTable->setContent(9, iType, zResolution.str());
  }

  // Finalize tables
  std::ostringstream sTotBarrelModule;
  sTotBarrelModule   << insur::web_emphStart << m_totalBarrelModules << insur::web_emphEnd;
  std::ostringstream sTotEndcapModule;
  sTotEndcapModule   << insur::web_emphStart << m_totalEndcapModules << insur::web_emphEnd;
  std::ostringstream sTrkTotArea;
  sTrkTotArea        << insur::web_emphStart << std::dec << std::fixed << std::setprecision(c_areaPrecision) << trkTotArea << insur::web_emphEnd;
  std::ostringstream sTrkTotModules;
  sTrkTotModules     << insur::web_emphStart << trkTotNumModules << insur::web_emphEnd;
  std::ostringstream sTrkTotNumChannels;
  sTrkTotNumChannels << insur::web_emphStart << std::dec << std::fixed << std::setprecision(c_channelPrecision) << trkTotNumChannels << insur::web_emphEnd;

  m_layerTable->setContent( 0, m_nBarrelLayers+1, "Total");
  m_layerTable->setContent( 6, m_nBarrelLayers+1, sTotBarrelModule.str());
  m_diskTable->setContent(  0, m_nDisks+1       , "Total");
  m_diskTable->setContent(  5, m_nDisks+1       , sTotEndcapModule.str());
  m_moduleTable->setContent(0, iType+1          , "Total");
  m_moduleTable->setContent(4, iType+1          , sTrkTotArea.str());
  m_moduleTable->setContent(5, iType+1          , sTrkTotModules.str());
  m_moduleTable->setContent(7, iType+1          , sTrkTotNumChannels.str());

  // Ring table
  for (auto& iterMap : m_ringModuleMap) {

    int                 ring   = iterMap.first;
    const EndcapModule* module = iterMap.second;

    m_ringTable->setContent(0, ring, ring);
    m_ringTable->setContent(1, ring, module->minR()                 , c_coordPrecision);
    m_ringTable->setContent(2, ring, module->minR()+module->length(), c_coordPrecision);
    m_ringTable->setContent(3, ring, m_ringNModules[ring-1]);
  }
}
