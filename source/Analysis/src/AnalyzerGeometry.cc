/*
 * AnalyzerGeometry.cc
 *
 *  Created on: 20. 4. 2016
 *      Author: drasal
 */
#include <AnalyzerGeometry.h>

// Include files
#include "Barrel.h"
#include "BeamPipe.h"
#include "Detector.h"
#include "DetectorModule.h"
#include "Disk.h"
#include "global_constants.h"
#include "Endcap.h"
#include "Layer.h"
#include "MaterialProperties.h"
#include <Math/Vector3D.h>
#include "PlotDrawer.h"
#include "RootWContent.h"
#include "RootWImage.h"
#include "RootWPage.h"
#include "RootWSite.h"
#include "RootWTable.h"
#include <ostream>
#include "TagMaker.h"
#include <TColor.h>
#include <TEllipse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2I.h>
#include <TPad.h>
#include <TProfile.h>
#include "Tracker.h"
#include <TRandom3.h>
#include "Units.h"
#include "SimParms.h"

//
// AnalyzerGeometry constructor
//
AnalyzerGeometry::AnalyzerGeometry(const Detector& detector) : AnalyzerUnit("AnalyzerGeometry", detector),
 m_nTracks(0),
 m_layerNamesVisitor(new VisitorLayerName(m_trackers)),
 m_etaSpan(2*SimParms::getInstance().getMaxEtaCoverage()),
 m_etaMin(-1*SimParms::getInstance().getMaxEtaCoverage()),
 m_etaMax(+1*SimParms::getInstance().getMaxEtaCoverage())
{
};

//
// Destructor
//
AnalyzerGeometry::~AnalyzerGeometry()
{
  m_layerNamesVisitor.reset(nullptr);
}

//
// AnalyzerGeometry init method
//
bool AnalyzerGeometry::init(int nGeomTracks)
{
  // Set nTracks
  m_nTracks = nGeomTracks;

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
    int nBinsProfile = c_nBinsProfile;
    m_nTracks        = sqrtNTracks*sqrtNTracks;

    // Compute shooting intervals to analyze geometry
    double etaMin = -1*SimParms::getInstance().getMaxEtaCoverage();
    double etaMax = +1*SimParms::getInstance().getMaxEtaCoverage();

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
      m_hitMapPhiEta[trkName].SetStats(kFALSE);

      if (m_trackMapPhiEta.find(trkName)==m_trackMapPhiEta.end()) m_trackMapPhiEta[trkName] = TH2I();
      m_trackMapPhiEta[trkName].Reset();
      m_trackMapPhiEta[trkName].SetNameTitle("TrackMapPhiEta", "Number of tracks;#phi;#eta");
      m_trackMapPhiEta[trkName].SetBins(nBins, -M_PI, M_PI, nBins, -m_etaMax, +m_etaMax);
      m_trackMapPhiEta[trkName].SetStats(kFALSE);

      if (m_moduleHitEtaProfile.find(trkName)==m_moduleHitEtaProfile.end()) m_moduleHitEtaProfile[trkName] = TProfile();
      m_moduleHitEtaProfile[trkName].Reset();
      m_moduleHitEtaProfile[trkName].SetNameTitle("ModuleHitEtaProfile","Number of modules with at least one hit;#eta;Number of hit modules");
      m_moduleHitEtaProfile[trkName].SetBins(nBinsProfile, 0, +m_etaMax);
      m_moduleHitEtaProfile[trkName].SetMarkerStyle(8);
      m_moduleHitEtaProfile[trkName].SetMarkerColor(1);
      m_moduleHitEtaProfile[trkName].SetMarkerSize(1.5);
      m_moduleHitEtaProfile[trkName].SetMinimum(0);
      m_moduleHitEtaProfile[trkName].SetStats(kFALSE);

      if (m_sensorHitEtaProfile.find(trkName)==m_sensorHitEtaProfile.end()) m_sensorHitEtaProfile[trkName] = TProfile();
      m_sensorHitEtaProfile[trkName].Reset();
      m_sensorHitEtaProfile[trkName].SetNameTitle("ModuleHitEtaProfile","Number of hit sensors;#eta;Number of hit sensors");
      m_sensorHitEtaProfile[trkName].SetBins(nBinsProfile, 0, +m_etaMax);
      m_sensorHitEtaProfile[trkName].SetMarkerStyle(8);
      m_sensorHitEtaProfile[trkName].SetMarkerColor(1);
      m_sensorHitEtaProfile[trkName].SetMarkerSize(1.5);
      m_sensorHitEtaProfile[trkName].SetMinimum(0);
      m_sensorHitEtaProfile[trkName].SetStats(kFALSE);

      if (m_stubHitEtaProfile.find(trkName)==m_stubHitEtaProfile.end()) m_stubHitEtaProfile[trkName] = TProfile();
      m_stubHitEtaProfile[trkName].Reset();
      m_stubHitEtaProfile[trkName].SetNameTitle("ModuleHitEtaProfile","Number of modules with a stub;#eta;Number of hit stubs");
      m_stubHitEtaProfile[trkName].SetBins(nBinsProfile, 0, +m_etaMax);
      m_stubHitEtaProfile[trkName].SetMarkerStyle(8);
      m_stubHitEtaProfile[trkName].SetMarkerColor(1);
      m_stubHitEtaProfile[trkName].SetMarkerSize(1.5);
      m_stubHitEtaProfile[trkName].SetMinimum(0);
      m_stubHitEtaProfile[trkName].SetStats(kFALSE);

      // Set hit profile histograms for different module types
      for (auto iModule : iTracker->modules()) {

        m_moduleTypes[trkName].insert(iModule->moduleType());
        m_moduleTypeColor[trkName][iModule->moduleType()] = iModule->plotColor();
      }

      if (m_moduleTypes[trkName].size()>1) for (auto iType : m_moduleTypes[trkName] ) {

        if (m_moduleTypeHitEtaProfile.find(trkName)!=m_moduleTypeHitEtaProfile.end()) {
          if (m_moduleTypeHitEtaProfile[trkName].find(iType)==m_moduleTypeHitEtaProfile[trkName].end()) m_moduleTypeHitEtaProfile[trkName][iType] = TProfile();
        }
        else m_moduleTypeHitEtaProfile[trkName][iType] = TProfile();
        m_moduleTypeHitEtaProfile[trkName][iType].Reset();
        m_moduleTypeHitEtaProfile[trkName][iType].SetNameTitle(std::string(trkName+"_"+iType+"_ModuleHitEtaProfile").c_str(), "Number of modules with at least one hit;#eta;Number of hit modules");
        m_moduleTypeHitEtaProfile[trkName][iType].SetBins(nBinsProfile, 0, +m_etaMax);
        m_moduleTypeHitEtaProfile[trkName][iType].SetMarkerStyle(8);
        m_moduleTypeHitEtaProfile[trkName][iType].SetMarkerColor(1);
        m_moduleTypeHitEtaProfile[trkName][iType].SetMarkerSize(1.0);
        m_moduleTypeHitEtaProfile[trkName][iType].SetMinimum(0);
        m_moduleTypeHitEtaProfile[trkName][iType].SetStats(kFALSE);
      }

      if (m_moduleTypes[trkName].size()>1) for (auto iType : m_moduleTypes[trkName] ) {

        if (m_moduleTypeStubEtaProfile.find(trkName)!=m_moduleTypeStubEtaProfile.end()) {
          if (m_moduleTypeStubEtaProfile[trkName].find(iType)==m_moduleTypeStubEtaProfile[trkName].end()) m_moduleTypeStubEtaProfile[trkName][iType] = TProfile();
        }
        else m_moduleTypeStubEtaProfile[trkName][iType] = TProfile();
        m_moduleTypeStubEtaProfile[trkName][iType].Reset();
        m_moduleTypeStubEtaProfile[trkName][iType].SetNameTitle(std::string(trkName+"_"+iType+"_ModuleHitEtaProfile").c_str(), "Number of modules with at least one hit;#eta;Number of hit modules");
        m_moduleTypeStubEtaProfile[trkName][iType].SetBins(nBinsProfile, 0, +m_etaMax);
        m_moduleTypeStubEtaProfile[trkName][iType].SetMarkerStyle(8);
        m_moduleTypeStubEtaProfile[trkName][iType].SetMarkerColor(1);
        m_moduleTypeStubEtaProfile[trkName][iType].SetMarkerSize(1.0);
        m_moduleTypeStubEtaProfile[trkName][iType].SetMinimum(0);
        m_moduleTypeStubEtaProfile[trkName][iType].SetStats(kFALSE);
      }

      // Prepare eta coverage profiles for individual barrel layers & end-cap disks
      if (m_layerEtaCoverProfile.find(trkName)!=m_layerEtaCoverProfile.end()) m_layerEtaCoverProfile[trkName].clear();
      if (m_layerStubEtaCoverProfile.find(trkName)!=m_layerStubEtaCoverProfile.end()) m_layerStubEtaCoverProfile[trkName].clear();

      std::set<std::string> layerNames;

      if (m_layerNamesVisitor->getLayerNames(trkName, layerNames)) for (auto layerName : layerNames) {

        // Set names, binning, etc.
        m_layerEtaCoverProfile[trkName][layerName] = TProfile();
        m_layerEtaCoverProfile[trkName][layerName].SetNameTitle(std::string("LayerEtaCoverage_"+trkName+"_"+layerName).c_str(), std::string(trkName+": Eta coverage of "+layerName+";#eta;Coverage").c_str());
        m_layerEtaCoverProfile[trkName][layerName].SetBins(2*nBinsProfile, -m_etaMax, +m_etaMax);
        m_layerEtaCoverProfile[trkName][layerName].SetMinimum(0);
        m_layerEtaCoverProfile[trkName][layerName].SetStats(kFALSE);

        m_layerStubEtaCoverProfile[trkName][layerName] = TProfile();
        m_layerStubEtaCoverProfile[trkName][layerName].SetNameTitle(std::string("LayerStubEtaCoverage_"+trkName+"_"+layerName).c_str(), std::string(trkName+": Eta coverage of "+layerName+" stub;#eta;Coverage").c_str());
        m_layerStubEtaCoverProfile[trkName][layerName].SetBins(2*nBinsProfile, -m_etaMax, +m_etaMax);
        m_layerStubEtaCoverProfile[trkName][layerName].SetMinimum(0);
        m_layerStubEtaCoverProfile[trkName][layerName].SetStats(kFALSE);
      }

    } // For trackers

    m_isInitOK = true;
    return m_isInitOK;
  }
}

//
// AnalyzerGeometry analysis method
//
bool AnalyzerGeometry::analyze()
{
  // Check that initialization OK
  if (!m_isInitOK) return false;

  // Initialize random number generator, counters and histograms
  TRandom3 myDice;
  myDice.SetSeed(random_seed);

  // Print
  std::cout << std::endl;

  // Go through all trackers
  for (auto iTracker : m_trackers) {

    std::string trkName =  iTracker->myid();

    // Print
    std::cout << " " << iTracker->myid() << " tracker"<< std::endl;

    // Counters
    std::map <std::string, int> moduleTypeCount;  // counts hit per module type -> if any of the sensors (or both) are hit, it counts 1
    std::map <std::string, int> sensorTypeCount;  // counts hit per sensor on module type -> if one sensor is hit, it counts 1, if both sensors are hit, it counts 2
    std::map <std::string, int> stubTypeCount;    // counts stubs per module type -> if both sensors are hit and a stub is formed, it counts 1

    // Analyze geometry layout
    for (auto iTrack=0; iTrack<m_nTracks; iTrack++) {

      // Generate a straight track & track origin
      double phi   = myDice.Rndm() * 2 * M_PI;
      double eta   = myDice.Rndm() * m_etaSpan + m_etaMin;
      double theta = 2*atan(exp(-1*eta));

      ROOT::Math::XYZVector direction = ROOT::Math::XYZVector(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));

      // Use uniform distribution to simulate position of primary interaction (rphi negligible). May be also triangular or gaussian (depends on accelerator design)!
      double zErrorIP = 0.;
      if (SimParms::getInstance().useLumiRegInAnalysis()) zErrorIP = SimParms::getInstance().zErrorIP();
      double zPos = (myDice.Rndm()*2 - 1)*zErrorIP;
      auto origin = ROOT::Math::XYZVector(0, 0, zPos);

      // Check whether generated track hits a module -> collect the list of hit modules
      std::vector<std::pair<DetectorModule*, HitType>> hitModules;

      static const double zErrorSafetyMargin = 5. ; // track origin shift in units of zError to compute boundaries

      // Go through all modules and save only those that were hit
      for (auto& module : iTracker->modules()) {

        // CouldHit -> save CPU time by calculating a cross-section of hit and given module within track phi and eta
        // (assuming origin error within +-5 sigma). Beware that module phi extends from -pi to +3*pi to avoid troubles
        // with -pi & +pi cross-line
        if (module->couldHit(direction, zErrorIP*zErrorSafetyMargin)) {

          // Get hit
          XYZVector hitPos;
          Material  hitMaterial;
          HitType   hitType;
          if (module->checkTrackHits(origin, direction, hitMaterial, hitType, hitPos)) {

            hitModules.push_back(std::make_pair(module,hitType));;
          }
        }
      }

      // Reset counters
      int numStubs = 0;
      int numHits  = 0;

      std::map<std::string, int> typedNumHits;
      std::map<std::string, int> typedNumStubs;
      for (auto iType : m_moduleTypes[trkName]) {

        typedNumHits[iType] = 0;
        typedNumStubs[iType]= 0;
      }

      // Get number of hits/stubs
      for (auto& module : hitModules) {

        auto modType = module.first->moduleType();

        // Inner
        if (module.second & HitType::INNER) numHits++;

        // Outer
        if (module.second & HitType::OUTER) numHits++;

        // Stubs
        if (module.second == HitType::STUB) {

          typedNumStubs[modType]++;
          numStubs++;
        }

        // Count module hits for given type
        if ((module.second & HitType::INNER) || (module.second & HitType::OUTER)) typedNumHits[modType]++;
      }

      // Fill plots
      m_hitMapPhiEta[trkName].Fill(direction.Phi(), direction.Eta(), hitModules.size());
      m_trackMapPhiEta[trkName].Fill(direction.Phi(), direction.Eta());

      m_moduleHitEtaProfile[trkName].Fill(fabs(direction.Eta()), hitModules.size());
      m_sensorHitEtaProfile[trkName].Fill(fabs(direction.Eta()), numHits);
      m_stubHitEtaProfile[trkName].Fill(fabs(direction.Eta()), numStubs);

      if (m_moduleTypes[trkName].size()>1) for (auto iType : m_moduleTypes[trkName]) m_moduleTypeHitEtaProfile[trkName][iType].Fill(fabs(direction.Eta()), typedNumHits[iType]);
      if (m_moduleTypes[trkName].size()>1) for (auto iType : m_moduleTypes[trkName]) m_moduleTypeStubEtaProfile[trkName][iType].Fill(fabs(direction.Eta()), typedNumStubs[iType]);

      std::set<std::string> layerNames;

      if (m_layerNamesVisitor->getLayerNames(trkName, layerNames)) for (auto layerName : layerNames) {

        int layerHit  = 0;
        int layerStub = 0;
        for (auto module : hitModules) {

          // Unique reference of module
          UniRef uniref = module.first->uniRef();
          if (layerName == (uniref.cnt + "_" + any2str(uniref.layer))) {

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

    // Set fill color - get first if several module types available
    int color = 1;
    if (m_moduleTypes[trkName].size()==1) color = Palette::color(m_moduleTypeColor[trkName].begin()->second);
    m_moduleHitEtaProfile[trkName].SetMarkerColor(color);
    m_sensorHitEtaProfile[trkName].SetMarkerColor(color);
    m_stubHitEtaProfile[trkName].SetMarkerColor(color);

    if (m_moduleTypes[trkName].size()>1) for (auto iType : m_moduleTypes[trkName]) m_moduleTypeHitEtaProfile[trkName][iType].SetMarkerColor(Palette::color(m_moduleTypeColor[trkName][iType]));
    if (m_moduleTypes[trkName].size()>1) for (auto iType : m_moduleTypes[trkName]) m_moduleTypeStubEtaProfile[trkName][iType].SetMarkerColor(Palette::color(m_moduleTypeColor[trkName][iType]));

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

  m_isAnalysisOK = true;
  return m_isAnalysisOK;
}

//
// AnalyzerGeometry visualization method
//
bool AnalyzerGeometry::visualize(RootWSite& webSite)
{
  // Check that initialization & analysis OK
  if (!m_isInitOK || !m_isAnalysisOK) return false;

  // Go through all trackers & prepare web content
  int webPriority         = web_priority_Geom;

  // List of tracker RZ, BRL XY, ECAP XY canvases
  std::map<std::string, std::unique_ptr<RootWImage>> trackerRZImg;
  std::map<std::string, std::unique_ptr<RootWImage>> trackerXYImgBRL;
  std::map<std::string, std::unique_ptr<RootWImage>> trackerXYImgEC;
  PlotDrawer<RZFull, Type>           fullRZDrawer;

  for (auto iTracker : m_trackers) {

    // Tracker name
    std::string trkName =  iTracker->myid();

    // Create dedicated web-page & its content
    std::string        pageTitle  = "Geometry";
    if (trkName != "") pageTitle +=" (" + trkName + ")";

    // Page address
    std::string pageAddress = "index"+trkName+".html";

    RootWPage& myPage = webSite.addPage(pageTitle, webPriority);
    myPage.setAddress(pageAddress);
    webPriority--;

    //
    // Layer & disk layout
    RootWContent& myContentLayout = myPage.addContent("Layer & disk layout");

    // Fill in all information directly from Tracker geometry object
    VisitorLayerDiscSummary geometryVisitor;
    geometryVisitor.preVisit();
    iTracker->accept(geometryVisitor);
    geometryVisitor.postVisit();

    // Print out layer & disk table
    myContentLayout.addItem(std::move(geometryVisitor.m_layerTable));
    myContentLayout.addItem(std::move(geometryVisitor.m_diskTable));
    myContentLayout.addItem(std::move(geometryVisitor.m_ringTable));

    // Print out detailed info about tilted layers
    TiltedLayersVisitor tiltLayerVisitor;
    iTracker->accept(tiltLayerVisitor);

    // If the layer is tilted, print flat and tilted parts tables
    if (tiltLayerVisitor.m_nTiltedLayers > 0) {

      RootWContent& myContent = myPage.addContent("Modules layout details",false);
      for (int i = 0; i < tiltLayerVisitor.m_nTiltedLayers; i++) {

        myContent.addItem(std::move(tiltLayerVisitor.m_tiltedLayerNames[i]));
        myContent.addItem(std::move(tiltLayerVisitor.m_flatPartTables[i]));
        myContent.addItem(std::move(tiltLayerVisitor.m_tiltedPartTables[i]));
      }
    }

    //
    // Modules
    RootWContent& myContentMods = myPage.addContent("Modules info", false);
    myContentMods.addItem(std::move(geometryVisitor.m_moduleTable));

    //
    // Plots
    RootWContent& myContentPlots = myPage.addContent("Layout plots");

    RootWImage* myImage;

    // R-Z view - avoid drawing canvas if barrel and endcap missing
    PlotDrawer<RZ, Type> rzDrawer;
    bool foundModules = rzDrawer.addModulesType(iTracker->modules().begin(), iTracker->modules().end(), BARREL | ENDCAP);
    if (foundModules) {

      // View
      TCanvas RZCanvas(std::string("RZCanvas_"+trkName).c_str(), "RZView Canvas", vis_max_canvas_sizeX, vis_min_canvas_sizeY);
      setCanvasProperties(RZCanvas);
      rzDrawer.drawFrame<TicksFrameStyle>(RZCanvas, iTracker->isPixelType());
      rzDrawer.drawModules<ContourStyle>(RZCanvas);
      fullRZDrawer.addModulesType(iTracker->modules().begin(), iTracker->modules().end());
      drawBeamPipeRZ(RZCanvas, iTracker->maxZ());

      RootWImage& myImage = myContentPlots.addImage(RZCanvas, vis_max_canvas_sizeX, vis_min_canvas_sizeY);
      myImage.setComment("RZ (envelope Radius versus Z) plot of the tracker");

      // Clone image view to have an exact copy
      std::unique_ptr<RootWImage> myImageClone(new RootWImage(RZCanvas, vis_max_canvas_sizeX, vis_min_canvas_sizeY));
      myImageClone->setComment("RZ (envelope Radius versus Z) plot of the tracker");
      trackerRZImg[trkName] = std::move(myImageClone);
    }

    // Projections - avoid drawing canvas if barrel or endcap missing
    PlotDrawer<XY, Type> xyBarrelDrawer;
    foundModules = xyBarrelDrawer.addModulesType(iTracker->modules().begin(), iTracker->modules().end(), BARREL);
    if (foundModules) {

      TCanvas XYCanvasBRL(std::string("XYCanvasBRL_"+trkName).c_str(), "XYView Barrel", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      setCanvasProperties(XYCanvasBRL);
      xyBarrelDrawer.drawFrame<TicksFrameStyle>(XYCanvasBRL);
      xyBarrelDrawer.drawModules<ContourStyle>(XYCanvasBRL);
      drawBeamPipeXY(XYCanvasBRL);

      RootWImage& myImage = myContentPlots.addImage(XYCanvasBRL, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage.setComment("XY cross section of barrel modules");

      // Clone image view to have an exact copy
      std::unique_ptr<RootWImage> myImageClone(new RootWImage(XYCanvasBRL, vis_min_canvas_sizeX, vis_min_canvas_sizeY));
      myImageClone->setComment("XY cross section of barrel modules");
      trackerXYImgBRL[trkName] = std::move(myImageClone);
    }

    PlotDrawer<XY, Type> xyEndcapDrawer;
    foundModules = xyEndcapDrawer.addModulesType(iTracker->modules().begin(), iTracker->modules().end(), ENDCAP);
    if (foundModules) {

      TCanvas XYCanvasEC(std::string("XYCanvasEC_"+trkName).c_str(), "XYView Endcap", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      setCanvasProperties(XYCanvasEC);
      xyEndcapDrawer.drawFrame<TicksFrameStyle>(XYCanvasEC);
      xyEndcapDrawer.drawModules<ContourStyle>(XYCanvasEC);
      drawBeamPipeXY(XYCanvasEC);

      RootWImage& myImage = myContentPlots.addImage(XYCanvasEC, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage.setComment("XY projection of endcap(s) modules");

      // Clone image view to have an exact copy
      std::unique_ptr<RootWImage> myImageClone(new RootWImage(XYCanvasEC, vis_min_canvas_sizeX, vis_min_canvas_sizeY));
      myImageClone->setComment("XY projection of endcap(s) modules");
      trackerXYImgEC[trkName] = std::move(myImageClone);
    }

    //
    // Hit eta profile plots

    // Sensor hit map
    if (m_sensorHitEtaProfile.find(trkName)!=m_sensorHitEtaProfile.end()) {

      TCanvas etaCanvasHitSensors(std::string("NumberOfHitSensorEtaProfile_"+trkName).c_str(), "Number of hit sensors versus #eta", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      setCanvasProperties(etaCanvasHitSensors);
      m_sensorHitEtaProfile[trkName].Draw();

      RootWImage& myImage = myContentPlots.addImage(etaCanvasHitSensors, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage.setComment("Number of hit sensors versus "+web_etaLetter);
    }

    // Module hit map
    if (m_moduleHitEtaProfile.find(trkName)!=m_moduleHitEtaProfile.end()) {

      TCanvas etaCanvasHitModules(std::string("NumberOfHitModulesEtaProfile_"+trkName).c_str(), "Number of hit modules versus #eta", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      setCanvasProperties(etaCanvasHitModules);
      m_moduleHitEtaProfile[trkName].Draw();

      if (m_moduleTypes[trkName].size()>1) for (auto iType : m_moduleTypes[trkName]) m_moduleTypeHitEtaProfile[trkName][iType].Draw("same");

      RootWImage& myImage = myContentPlots.addImage(etaCanvasHitModules, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage.setComment("Number of hit modules versus "+web_etaLetter);
    }

    // 2D hit map
    if (m_hitMapPhiEta.find(trkName)!=m_hitMapPhiEta.end()) {

      TCanvas etaCanvasHitMap(std::string("HitMap2D_"+trkName).c_str(), "Hit map - #eta versus #phi", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      setCanvasProperties(etaCanvasHitMap);
      m_hitMapPhiEta[trkName].Draw("colz");

      RootWImage& myImage = myContentPlots.addImage(etaCanvasHitMap, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage.setComment("Hit coverage in "+web_etaLetter+", "+web_phiLetter);
    }

    // Layer coverage plots
    if ((m_layerEtaCoverProfile.find(trkName)!=m_layerEtaCoverProfile.end()) && m_layerEtaCoverProfile[trkName].size()>0) {

      RootWContent& myContentHits = myPage.addContent("Layer coverage (Hits)", false);

      int layerCount = 0;
      for (auto& iter : m_layerEtaCoverProfile[trkName]) {

        TProfile& plot = iter.second;
        layerCount++;

        TCanvas coverageCanvas(std::string(Form("LayerCoverage%s%s_", iter.first.c_str(), "Hits")+trkName).c_str(), "Layer eta coverage (Hits)", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
        setCanvasProperties(coverageCanvas);

        std::string upperName = std::string(coverageCanvas.GetName())+"_upper";
        TPad upperPad(upperName.c_str(), "upperPad", 0, 0.3, 1, 1);
        upperPad.Draw();
        std::string lowerName = std::string(coverageCanvas.GetName())+"_lower";
        TPad lowerPad(lowerName.c_str(), "lowerPad", 0, 0.0, 1, 0.3);
        lowerPad.Draw();

        plot.SetMinimum(0);
        plot.SetMaximum(1.01);
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

        RootWImage& myImage = myContentHits.addImage(coverageCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
        myImage.setComment("Layer coverage in "+web_etaLetter+" (efficiency)");
      }
    }
  } // Trackers

  // Create introductory geometry web page
  std::string pageAddress = "index.html";
  RootWPage& myPage = webSite.addPage("Overview", web_priority_Geom+1);
  myPage.setAddress(pageAddress);

  // Full tracker layout
  RootWContent& myContentFullLayout = myPage.addContent("Tracker layout: ");

  TCanvas fullRZCanvas("RZCanvasFull", "RZView Canvas Full", vis_max_canvas_sizeX, vis_min_canvas_sizeY);
  setCanvasProperties(fullRZCanvas);
  fullRZDrawer.drawFrame<TicksFrameStyle>(fullRZCanvas);
  fullRZDrawer.drawModules<ContourStyle>(fullRZCanvas);

  myContentFullLayout.addImage(fullRZCanvas, vis_max_canvas_sizeX, vis_min_canvas_sizeY);

  // Hits
  TCanvas hitCanvas("HitCanvas", "Number of modules with at least one hit;#eta;Number of hit modules", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
  setCanvasProperties(hitCanvas);

  TH1D moduleTotHitsEtaHis;
  moduleTotHitsEtaHis.SetNameTitle("ModuleHitEtaProfile","Number of modules with at least one hit;#eta;Number of hit modules");
  moduleTotHitsEtaHis.SetBins(c_nBinsProfile, 0, +m_etaMax);
  moduleTotHitsEtaHis.SetMarkerStyle(8);
  moduleTotHitsEtaHis.SetMarkerColor(1);
  moduleTotHitsEtaHis.SetMarkerSize(1.5);
  moduleTotHitsEtaHis.SetMinimum(0);
  moduleTotHitsEtaHis.SetStats(kFALSE);

  for (auto iTracker : m_trackers) {
    std::string trkName = iTracker->myid();
    moduleTotHitsEtaHis.Add(&m_moduleHitEtaProfile[trkName]);
  }
  moduleTotHitsEtaHis.SetMinimum(0);
  moduleTotHitsEtaHis.Draw();

  for (auto iTracker : m_trackers) {
    std::string trkName = iTracker->myid();
    if (m_moduleTypes[trkName].size()>1) for (auto iType : m_moduleTypes[trkName]) m_moduleTypeHitEtaProfile[trkName][iType].Draw("same");
  }

  myContentFullLayout.addImage(hitCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);

  // All subdetectors
  for (auto& iter : trackerRZImg) {

    std::string trkName = iter.first;

    RootWContent& myContentSubDet = myPage.addContent("Geometry(" + trkName + "): ");
    myContentSubDet.addItem(std::move(iter.second));
    if (trackerXYImgBRL.find(trkName)!=trackerXYImgBRL.end()) myContentSubDet.addItem(std::move(trackerXYImgBRL[trkName]));
    if (trackerXYImgEC.find(trkName) !=trackerXYImgEC.end())  myContentSubDet.addItem(std::move(trackerXYImgEC[trkName]));
  }

  m_isVisOK = true;
  return m_isVisOK;
}

//
// Draw beam-pipe in RZ to given canvas
//
void AnalyzerGeometry::drawBeamPipeRZ(TCanvas& canvas, double maxZ)
{
  if (m_beamPipe) {
    double bpRadius = m_beamPipe->radius();
    double bpThick  = m_beamPipe->thickness();

    // Beam-pipe in RZ
    TPolyLine beamPipeRZ;
    beamPipeRZ.SetPoint(0, 0                     , bpRadius + bpThick/2.);
    beamPipeRZ.SetPoint(1, maxZ*vis_safety_factor, bpRadius + bpThick/2.);
    beamPipeRZ.SetLineColor(14);
    beamPipeRZ.SetLineWidth(2);
    canvas.cd();
    beamPipeRZ.DrawClone("same");
  }
  else {
    logWARNING("DrawBeamPipeRZ failed: no beam pipe defined!");
  }
}

//
// Draw beam-pipe in XY to given canvas
//
void AnalyzerGeometry::drawBeamPipeXY(TCanvas& canvas)
{
  if (m_beamPipe) {
    double bpRadius = m_beamPipe->radius();
    double bpThick  = m_beamPipe->thickness();

    // Beam-pipe in XY
    TEllipse beamPipeXY(0, 0, bpRadius + bpThick/2.);
    beamPipeXY.SetFillColor(18); // "grey18"
    beamPipeXY.SetFillStyle(1001);
    canvas.cd();
    beamPipeXY.DrawClone("same");
  }
  else {
    logWARNING("DrawBeamPipeXY failed: no beam pipe defined!");
  }
}

//
// VisitorLayerName constructor method
//
VisitorLayerName::VisitorLayerName(std::vector<const Tracker*>& trackers)
{
  for (auto tracker : trackers) {

    m_idTRK = tracker->myid();
    tracker->accept(*this);
  }
}

//
// VisitorLayerName - fill container with layer names for defined tracker if tracker exists
//
bool VisitorLayerName::getLayerNames(string trkName, std::set<std::string>& layerNames)
{
  if (m_data.find(trkName) != m_data.end()) {

    for (auto iter : m_data[trkName]) layerNames.insert(iter);
    return true;
  }
  else return false;
}

//
// VisitorLayerName - visit barrel
//
void VisitorLayerName::visit(const Barrel& b) { m_idBRLorEC = b.myid(); }

//
// VisitorLayerName - visit endcap
//
void VisitorLayerName::visit(const Endcap& e) { m_idBRLorEC = e.myid(); }

//
// VisitorLayerName - visit layer
//
void VisitorLayerName::visit(const Layer& l)  { m_data[m_idTRK].insert(m_idBRLorEC + "_" + any2str(l.myid())); }

//
// VisitorLayerName - visit disk
//
void VisitorLayerName::visit(const Disk& d)   { m_data[m_idTRK].insert(m_idBRLorEC + "_" + any2str(d.myid())); }

//
// VisitorLayerName destructor
//
VisitorLayerDiscSummary::~VisitorLayerDiscSummary() {};

//
// VisitorLayerDiscSummary preVisit method - initialize
//
void VisitorLayerDiscSummary::preVisit() {

  // Initialize
  m_nBarrelLayers      = 0;
  m_nDisks             = 0;
  m_nRings             = 0;
  m_totalBarrelModules = 0;
  m_totalEndcapModules = 0;

  m_totalArea       = 0;
  m_totalModules    = 0;
  m_totalSensors    = 0;
  m_totalChannels   = 0;
  m_totalSensorPower= 0;

  m_layerTable = std::unique_ptr<RootWTable>(new RootWTable());
  m_diskTable  = std::unique_ptr<RootWTable>(new RootWTable());
  m_ringTable  = std::unique_ptr<RootWTable>(new RootWTable());
  m_moduleTable= std::unique_ptr<RootWTable>(new RootWTable());

  m_layerTable->setContent(0, 0, "Layer no                    : ");
  m_layerTable->setContent(1, 0, "Average radius [mm]         : ");
  m_layerTable->setContent(2, 0, "Radius-min [mm]             : ");
  m_layerTable->setContent(3, 0, "Radius-max [mm]             : ");
  m_layerTable->setContent(4, 0, "Z-min [mm]                  : ");
  m_layerTable->setContent(5, 0, "Z-max [mm]                  : ");
  m_layerTable->setContent(6, 0, "Number of rods              : ");
  m_layerTable->setContent(7, 0, "Number of modules per rod   : ");
  m_layerTable->setContent(8, 0, "Number of modules           : ");

  m_diskTable->setContent(0, 0, "Disk no                      : ");
  m_diskTable->setContent(1, 0, "Radius-min [mm]              : ");
  m_diskTable->setContent(2, 0, "Radius-max [mm]              : ");
  m_diskTable->setContent(3, 0, "Average Z pos. [mm]          : ");
  m_diskTable->setContent(4, 0, "Z-min [mm]                   : ");
  m_diskTable->setContent(5, 0, "Z-max [mm]                   : ");
  m_diskTable->setContent(6, 0, "Number of rings              : ");
  m_diskTable->setContent(7, 0, "Number of modules per disk   : ");

  m_ringTable->setContent(0, 0, "Ring no                      : ");
  m_ringTable->setContent(1, 0, "R-min [mm]                   : ");
  m_ringTable->setContent(2, 0, "R-max [mm]                   : ");
  m_ringTable->setContent(3, 0, "Number of modules per ring   : ");

  m_moduleTable->setContent( 0, 0, "Module in:                  ");
  m_moduleTable->setContent( 1, 0, "Position:                   ");
  m_moduleTable->setContent( 2, 0, "Type:                       ");
  m_moduleTable->setContent( 3, 0, "Sensor area [mm"+web_superStart+"2"+web_superEnd+"]: ");
  m_moduleTable->setContent( 4, 0, "Total area [m"+web_superStart+"2"+web_superEnd+"]:   ");
//  m_moduleTable->setContent( 5, 0, "Service Weight [kg]:                         ");
//  m_moduleTable->setContent( 6, 0, "Total Weight [kg]:                           ");
  m_moduleTable->setContent( 5, 0, "Number of modules:          ");
  m_moduleTable->setContent( 6, 0, "Number of sensors:          ");
  m_moduleTable->setContent( 7, 0, "Number of channels (M):     ");
  m_moduleTable->setContent( 8, 0, "Channels per module (R-Phi): ");
  m_moduleTable->setContent( 9, 0, "Channels per module (Z):     ");
  m_moduleTable->setContent(10, 0, "Read-out chips per module:  ");
//  m_moduleTable->setContent(12, 0, "Number of channels (Z 1.side):               ");
//  m_moduleTable->setContent(13, 0, "Number of channels (Z 2.side):               ");
  m_moduleTable->setContent(11, 0, "Min-Max R-Phi resolution ("+web_muLetter+"m): ");
  m_moduleTable->setContent(12, 0, "Min-Max Z(R) resolution ("+web_muLetter+"m):  ");
}

//
// VisitorLayerDiscSummary visit layer
//
void VisitorLayerDiscSummary::visit(const Layer& l) {

  if (l.maxZ() < 0.) return;

  // Update layer counter
  ++m_nBarrelLayers;

  // Update module counter
  int nModules        = l.totalModules();
  m_totalBarrelModules += nModules;

  // Set table
  m_layerTable->setContent(0, m_nBarrelLayers, m_nBarrelLayers);
  m_layerTable->setContent(1, m_nBarrelLayers, l.avgBuildRadius()/Units::mm, c_coordPrecision+1);
  m_layerTable->setContent(2, m_nBarrelLayers, l.minRAllMat()/Units::mm, c_coordPrecision+1);
  m_layerTable->setContent(3, m_nBarrelLayers, l.maxRAllMat()/Units::mm, c_coordPrecision+1);
  m_layerTable->setContent(4, m_nBarrelLayers, l.minZ()/Units::mm, c_coordPrecision);
  m_layerTable->setContent(5, m_nBarrelLayers, l.maxZ()/Units::mm, c_coordPrecision);
  m_layerTable->setContent(6, m_nBarrelLayers, l.numRods());
  m_layerTable->setContent(7, m_nBarrelLayers, l.numModulesPerRod());
  m_layerTable->setContent(8, m_nBarrelLayers, l.totalModules());
}

//
// VisitorLayerDiscSummary visit disk
//
void VisitorLayerDiscSummary::visit(const Disk& d) {

  if (d.averageZ() < 0.) return;

  // Update disk counter
  if (d.numRings()>m_nRings) m_ringNModules.resize(d.numRings());
  m_nRings = d.numRings();
  ++m_nDisks;

  // Update module counter -> counting only +Z discs
  int nModules = d.totalModules();
  m_totalEndcapModules += nModules;

  // Set table
  m_diskTable->setContent(0, m_nDisks, d.myid());
  m_diskTable->setContent(1, m_nDisks, d.minR()/Units::mm,       c_coordPrecision);
  m_diskTable->setContent(2, m_nDisks, d.maxR()/Units::mm,       c_coordPrecision);
  m_diskTable->setContent(3, m_nDisks, d.averageZ()/Units::mm,   c_coordPrecision);
  m_diskTable->setContent(4, m_nDisks, d.minZAllMat()/Units::mm, c_coordPrecision);
  m_diskTable->setContent(5, m_nDisks, d.maxZAllMat()/Units::mm, c_coordPrecision);
  m_diskTable->setContent(6, m_nDisks, d.numRings());
  m_diskTable->setContent(7, m_nDisks, nModules);
}

//
// VisitorLayerDiscSummary visit ring
//
void VisitorLayerDiscSummary::visit(const Ring& r) {

  if (r.averageZ() < 0.) return;

  m_ringNModules[r.myid()-1] = r.numModules();
}

//
// VisitorLayerDiscSummary visit module (barrel + endcap)
//
void VisitorLayerDiscSummary::visit(const DetectorModule& m) {

  // Get unique sensor tag
  TagMaker moduleTagMaker(m);
  std::string tag = moduleTagMaker.sensorResTag;

  // Relate tags to positions
  m_moduleTagToPositionsMap[tag].insert(moduleTagMaker.posTag);

  // Get number of minimum bias events
  int nMinBiasEvents = SimParms::getInstance().numMinBiasEvents();

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
  m_moduleAvgChannelsRPhi[tag]   += m.innerSensor().numStripsAcross();
  m_moduleAvgChannelsZ[tag]      += m.innerSensor().numSegments();
  m_moduleAvgROCs[tag]           += m.innerSensor().numROCX()*m.innerSensor().numROCY();
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
// VisitorLayerDiscSummary visit end-cap module
//
void VisitorLayerDiscSummary::visit(const EndcapModule& m) {

  // All disks symmetric - visit only the first one (+-z symmetry -> check only positive discs)
  if (m.disk() != 1 || m.side() != 1) return;

  // Get all end-cap module types for given ring - if not find -> assign
  if (m_ringModuleMap.find(m.ring())==m_ringModuleMap.end()) m_ringModuleMap[m.ring()] = &m;
}

//
// VisitorLayerDiscSummary postvisit method - finalize
//
void VisitorLayerDiscSummary::postVisit() {

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
      m_moduleAvgChannelsRPhi[tag]   /= normFactor;
      m_moduleAvgChannelsZ[tag]      /= normFactor;
      m_moduleAvgROCs[tag]           /= normFactor;
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
    tag << web_smallStart;
    for (auto iter : m_moduleTagToPositionsMap[moduleType]) tag << iter << "<br/> ";
    tag << web_smallEnd;

    // Type
    std::string type = module->moduleType();

    // Thickness
    std::string thickness = any2str(module->dsDistance());

    // Area
    std::string sensorArea = any2str(module->area()/Units::mm2, c_areaPrecision);
    std::string totalArea  = any2str(module->area()/Units::m2*module->numSensors()*m_moduleCount[moduleType], c_areaPrecision);
    trkTotArea += module->area()/Units::m2*module->numSensors()*m_moduleCount[moduleType];

    // Occupancy (in percent)
    std::string stripOccupancy = any2str(m_moduleMaxStripOccupancy[moduleType]/Units::percent, c_occupancyPrecision) + "/" + any2str(m_moduleAvgStripOccupancy[moduleType]/Units::percent, c_occupancyPrecision); // Percent
    std::string hitOccupancy   = any2str(m_moduleMaxHitOccupancy[moduleType]/Units::percent  , c_occupancyPrecision) + "/" + any2str(m_moduleAvgHitOccupancy[moduleType]/Units::percent  , c_occupancyPrecision); // Percent

    // Resolution (in um)
    std::string rphiResolution = any2str(m_moduleMinRphiResolution[moduleType]/Units::um, c_resolutionPrecision) + "-" + any2str(m_moduleMaxRphiResolution[moduleType]/Units::um, c_resolutionPrecision);
    std::string zResolution    = any2str(m_moduleMinZResolution[moduleType]/Units::um   , c_resolutionPrecision) + "-" + any2str(m_moduleMaxZResolution[moduleType]/Units::um   , c_resolutionPrecision);

    // Number of modules
    std::string numberMod = any2str(m_moduleCount[moduleType]);
    trkTotNumModules += m_moduleCount[moduleType];

    // Number of sensors
    std::string numberSens = any2str(m_moduleCount[moduleType] * module->numSensors());

    // Number of channels (in millions)
    std::string numberChan = any2str(m_moduleChannels[moduleType]/Units::M, c_channelPrecision);
    trkTotNumChannels     += m_moduleChannels[moduleType]/Units::M;

    // Number of channels per module
    std::string numberChanRPhi = any2str(m_moduleAvgChannelsRPhi[moduleType]);
    std::string numberChanZ    = any2str(m_moduleAvgChannelsZ[moduleType]);
    std::string numberROCs     = any2str(m_moduleAvgROCs[moduleType]);

    // Fill
    m_moduleTable->setContent(0, iType, position.str());
    m_moduleTable->setContent(1, iType, tag.str());
    m_moduleTable->setContent(2, iType, type);
    m_moduleTable->setContent(3, iType, sensorArea);
    m_moduleTable->setContent(4, iType, totalArea);
    m_moduleTable->setContent(5, iType, numberMod);
    m_moduleTable->setContent(6, iType, numberSens);
    m_moduleTable->setContent(7, iType, numberChan);
    m_moduleTable->setContent(8, iType, numberChanRPhi);
    m_moduleTable->setContent(9, iType, numberChanZ);
    m_moduleTable->setContent(10,iType, numberROCs);
    m_moduleTable->setContent(11,iType, rphiResolution);
    m_moduleTable->setContent(12,iType, zResolution);
  }

  // Finalize tables
  std::string sTotBarrelModules = web_emphStart + any2str(m_totalBarrelModules)        + web_emphEnd;
  std::string sTotEndcapModules = web_emphStart + any2str(2*m_totalEndcapModules)      + web_emphEnd; // Disk calculations done only for one side due to symmetry, hence factor of 2 needed
  std::string sTrkTotArea       = web_emphStart + any2str(trkTotArea, c_areaPrecision) + web_emphEnd;;
  std::string sTrkTotModules    = web_emphStart + any2str(trkTotNumModules)            + web_emphEnd;
  std::string sTrkTotNumChannels= web_emphStart + any2str(trkTotNumChannels, c_channelPrecision) + web_emphEnd;

  m_layerTable->setContent( 0, m_nBarrelLayers+1, "Total");
  m_layerTable->setContent( 8, m_nBarrelLayers+1, sTotBarrelModules);
  m_diskTable->setContent(  0, m_nDisks+1       , "Total (+Z "+web_ampersand+" -Z)");
  m_diskTable->setContent(  7, m_nDisks+1       , sTotEndcapModules);
  m_moduleTable->setContent(0, iType+1          , "Total");
  m_moduleTable->setContent(4, iType+1          , sTrkTotArea);
  m_moduleTable->setContent(5, iType+1          , sTrkTotModules);
  m_moduleTable->setContent(7, iType+1          , sTrkTotNumChannels);

  // Ring table
  for (auto& iterMap : m_ringModuleMap) {

    int                 ring   = iterMap.first;
    const EndcapModule* module = iterMap.second;

    m_ringTable->setContent(0, ring, ring);
    m_ringTable->setContent(1, ring, module->minR()/Units::mm, c_coordPrecision);
    m_ringTable->setContent(2, ring, module->maxR()/Units::mm, c_coordPrecision);
    m_ringTable->setContent(3, ring, m_ringNModules[ring-1]);
  }
}

//
// Visit tilted tracker and initialize variables
//
void TiltedLayersVisitor::visit(const Tracker& t) {

  m_nTiltedLayers = 0;
  m_nLayers       = 0;
}

//
// Visit tilted layers and gathers info on its flat and tilted parts.
//
void TiltedLayersVisitor::visit(const Layer& l) {

  m_nLayers++;

  // Only for tilted layers with automatic placement.
  if (l.isTilted() && l.isTiltedAuto()) {

    m_nTiltedLayers++;

    // Initializes layer name
    std::unique_ptr<RootWTable> tiltedLayerName = std::unique_ptr<RootWTable>(new RootWTable());
    tiltedLayerName->setContent(0, 0, "Layer " + any2str(m_nLayers) + " :");
    m_tiltedLayerNames.push_back(std::move(tiltedLayerName));

    //
    // Table summarizing tilted part
    std::unique_ptr<RootWTable> tiltedPartTable = std::unique_ptr<RootWTable>(new RootWTable());
    for (auto iRing=0; iRing < l.numTiltedRings(); iRing++) {

      int totRingNumber = l.buildNumModulesFlat() + iRing + 1;

      tiltedPartTable->setContent(0 , 0  , "Tilted part - Rings :");
      tiltedPartTable->setContent(0 , iRing+1, totRingNumber);
      tiltedPartTable->setContent(1 , 0  , web_thetaLetter + web_subStart + "tilt" + web_subEnd + " [deg] :");
      tiltedPartTable->setContent(1 , iRing+1, l.tiltedRing(iRing).tiltAngle()/Units::degrees, c_anglePrecision);
      tiltedPartTable->setContent(2 , 0  , "Ideal " + web_thetaLetter + web_subStart + "tilt" + web_subEnd + "-inner [deg] : ");
      tiltedPartTable->setContent(2 , iRing+1, l.tiltedRing(iRing).tiltAngleIdealInner()/Units::degrees, c_anglePrecision);
      tiltedPartTable->setContent(3 , 0  , "Ideal " + web_thetaLetter + web_subStart + "tilt" + web_subEnd + "-outer [deg] : ");
      tiltedPartTable->setContent(3 , iRing+1, l.tiltedRing(iRing).tiltAngleIdealOuter()/Units::degrees, c_anglePrecision);
      tiltedPartTable->setContent(4 , 0  , web_DeltaLetter + web_thetaLetter + web_subStart + "tilt" + web_subEnd +"-inner [deg] :");
      tiltedPartTable->setContent(4 , iRing+1, l.tiltedRing(iRing).deltaTiltIdealInner()/Units::degrees, c_anglePrecision);
      tiltedPartTable->setContent(5 , 0  , web_DeltaLetter + web_thetaLetter + web_subStart + "tilt" + web_subEnd +"-outer [deg] :");
      tiltedPartTable->setContent(5 , iRing+1, l.tiltedRing(iRing).deltaTiltIdealOuter()/Units::degrees, c_anglePrecision);
      tiltedPartTable->setContent(6 , 0  , web_thetaLetter + web_subStart + "g" + web_subEnd + " (centre"+web_subStart+"inner"+web_subEnd+"->center"+web_subStart+"outer"+web_subEnd+") [deg] :");
      tiltedPartTable->setContent(6 , iRing+1, l.tiltedRing(iRing).theta_g()/Units::degrees, c_anglePrecision);
      tiltedPartTable->setContent(7 , 0  , "R"+web_subStart + "centre" + web_subEnd + "-inner [mm] :");
      tiltedPartTable->setContent(7 , iRing+1, l.tiltedRing(iRing).innerRadius()/Units::mm, c_tiltedCoordPrecision);
      tiltedPartTable->setContent(8 , 0  , "R"+web_subStart + "centre" + web_subEnd + "-outer [mm] :");
      tiltedPartTable->setContent(8 , iRing+1, l.tiltedRing(iRing).outerRadius()/Units::mm, c_tiltedCoordPrecision);
      tiltedPartTable->setContent(9 , 0  , "R"+web_subStart + "centre" + web_subEnd + "-avg [mm] :");
      tiltedPartTable->setContent(9 , iRing+1, l.tiltedRing(iRing).averageR()/Units::mm, c_tiltedCoordPrecision);
      tiltedPartTable->setContent(10, 0  , web_DeltaLetter +"/2 (R+-"+web_DeltaLetter+" @"+web_thetaLetter + web_subStart + "tilt" + web_subEnd + "=0) [mm] :");
      tiltedPartTable->setContent(10, iRing+1, l.tiltedRing(iRing).gapR()/Units::mm, c_tiltedCoordPrecision);
      tiltedPartTable->setContent(11, 0  , "Z" + web_subStart + "centre" + web_subEnd +"-inner [mm] :");
      tiltedPartTable->setContent(11, iRing+1, l.tiltedRing(iRing).zInner()/Units::mm, c_tiltedCoordPrecision);
      tiltedPartTable->setContent(12, 0  , "Z" + web_subStart + "centre" + web_subEnd +"-outer [mm] :");
      tiltedPartTable->setContent(12, iRing+1, l.tiltedRing(iRing).zOuter()/Units::mm, c_tiltedCoordPrecision);
      tiltedPartTable->setContent(13, 0  , "Z"+web_subStart + "centre" + web_subEnd + "-avg [mm] :");
      tiltedPartTable->setContent(13, iRing+1, l.tiltedRing(iRing).averageZ()/Units::mm, c_tiltedCoordPrecision);
      tiltedPartTable->setContent(14, 0  , web_DeltaLetter + "Z-inner (Ring"+web_subStart+"i"+web_subEnd+"-Ring"+web_subStart+"i-1"+web_subEnd+") [mm] :");
      tiltedPartTable->setContent(14, iRing+1, l.tiltedRing(iRing).deltaZInner()/Units::mm, c_tiltedCoordPrecision);
      tiltedPartTable->setContent(15, 0  , web_DeltaLetter + "Z-outer (Ring"+web_subStart+"i"+web_subEnd+"-Ring"+web_subStart+"i-1"+web_subEnd+") [mm] :");
      tiltedPartTable->setContent(15, iRing+1, l.tiltedRing(iRing).deltaZOuter()/Units::mm, c_tiltedCoordPrecision);
      tiltedPartTable->setContent(16, 0  , "Overlap in R-"+ web_phiLetter + " [mm] :");
      tiltedPartTable->setContent(16, iRing+1, l.tiltedRing(iRing).phiOverlap()/Units::mm, c_tiltedCoordPrecision);
      tiltedPartTable->setContent(17, 0  , "Overlap in Z [mm] :");
      tiltedPartTable->setContent(17, iRing+1, l.tiltedRing(iRing).ringZOverlap()/Units::mm, c_zOverlapPrecision);

      // Inner modules in the tilted rings
      double zErrorInner = l.tiltedRing(iRing).zErrorInner();
      tiltedPartTable->setContent(18, 0, "zError-inner (Ring"+web_subStart+"i"+web_subEnd+"-Ring"+web_subStart+"i-1"+web_subEnd+") [mm] :");
      if (!std::isnan(zErrorInner)) tiltedPartTable->setContent(18, iRing+1, zErrorInner/Units::mm, c_tiltedCoordPrecision);
      else tiltedPartTable->setContent(18, iRing+1, "n/a");

      // Outer modules in the tilted rings
      double zErrorOuter = l.tiltedRing(iRing).zErrorOuter();
      tiltedPartTable->setContent(19, 0, "zError-outer (Ring"+web_subStart+"i"+web_subEnd+"-Ring"+web_subStart+"i-1"+web_subEnd+") [mm] :");
      if (!std::isnan(zErrorOuter)) tiltedPartTable->setContent(19, iRing+1, zErrorOuter/Units::mm, c_tiltedCoordPrecision);
      else tiltedPartTable->setContent(19, iRing+1, "n/a");
    }
    m_tiltedPartTables.push_back(std::move(tiltedPartTable));

    //
    // Table summarizing flat part
    std::unique_ptr<RootWTable> flatPartTable = std::unique_ptr<RootWTable>(new RootWTable());

    const auto& minusBigDeltaRod     = dynamic_cast<const RodPairStraight&>(l.m_bigParity() > 0 ? l.flatRods()[1] : l.flatRods()[0]);
    const auto& minusBigDeltaModules = minusBigDeltaRod.modules().first;
    const auto& plusBigDeltaRod      = dynamic_cast<const RodPairStraight&>(l.m_bigParity() > 0 ? l.flatRods()[0] : l.flatRods()[1]);
    const auto& plusBigDeltaModules  = plusBigDeltaRod.modules().first;

    int iRing = 0;
    for (const auto& m : minusBigDeltaModules) {

      int totRingNumber = iRing+1;

      flatPartTable->setContent(0, 0  , "Flat part - Rings:");
      flatPartTable->setContent(0, iRing+1, totRingNumber);
      flatPartTable->setContent(1, 0  , "R"+web_subStart + "centre" + web_subEnd + "-inner [mm] :");
      flatPartTable->setContent(1, iRing+1, m.center().Rho()/Units::mm, c_tiltedCoordPrecision);
      flatPartTable->setContent(3, 0  , "R"+web_subStart + "centre" + web_subEnd + "-avg [mm] :");
      flatPartTable->setContent(3, iRing+1, l.avgBuildRadiusFlat()/Units::mm, c_tiltedCoordPrecision);
      flatPartTable->setContent(4, 0  , web_DeltaLetter +" (R+-"+web_DeltaLetter+" in R-Phi) [mm] :");
      flatPartTable->setContent(4, iRing+1, l.bigDelta()/Units::mm, c_tiltedCoordPrecision);
      flatPartTable->setContent(5, 0  , web_deltaLetter +" (R+-"+web_deltaLetter+" in Z) [mm] :");
      flatPartTable->setContent(5, iRing+1, l.smallDelta()/Units::mm, c_tiltedCoordPrecision);
      flatPartTable->setContent(6, 0  , "Z"+web_subStart + "centre" + web_subEnd + "-avg [mm] :");
      flatPartTable->setContent(6, iRing+1, m.center().Z()/Units::mm, c_tiltedCoordPrecision);
      flatPartTable->setContent(7, 0  , "Overlap in R-"+ web_phiLetter + " [mm] :");
      flatPartTable->setContent(7, iRing+1, (((minusBigDeltaRod.smallParity() * pow(-1, (iRing%2))) > 0) ? l.flatPhiOverlapSmallDeltaPlus()/Units::mm : l.flatPhiOverlapSmallDeltaMinus()/Units::mm), c_tiltedCoordPrecision);

      // In case beamSpotCover == false, zOverlap is the only parameter used as a Z-coverage constraint in the geometry construction process.
      // (There is then no zError taken into account) -> As a result, it is interesting to display zOverlap !
      if (SimParms::getInstance().useLumiRegInGeomBuild()) {
        flatPartTable->setContent(8, 0  , "Overlap in Z [mm] :");
        flatPartTable->setContent(8, iRing+1, (dynamic_cast<const RodPairStraight&>(l.flatRods()[0])).zOverlap()/Units::mm, c_tiltedCoordPrecision);
	      // WARNING : zOverlap, in the geometry construction process, is one value common for all flat part (or straight rod)
      }
      // calculates zError for ring (i) with respect to ring (i-1)
      if (iRing > 0) {

        int extraLine = SimParms::getInstance().useLumiRegInGeomBuild() ? 1 : 0;

        // Inner modules in the tilted rings
        double zErrorInner = l.zErrorInnerFlat(iRing);
        flatPartTable->setContent(8 + extraLine, 0, "zError-inner (Ring"+web_subStart+"i"+web_subEnd+"-Ring"+web_subStart+"i-1"+web_subEnd+") [mm] :");
        if (!std::isnan(zErrorInner)) flatPartTable->setContent(8 + extraLine, iRing+1, zErrorInner/Units::mm, c_tiltedCoordPrecision);
        else                          flatPartTable->setContent(8 + extraLine, iRing+1, "n/a");

        // Outer modules in the tilted rings
        double zErrorOuter = l.zErrorOuterFlat(iRing);
        flatPartTable->setContent(9 + extraLine, 0, "zError-outer (Ring"+web_subStart+"i"+web_subEnd+"-Ring"+web_subStart+"i-1"+web_subEnd+") [mm] :");
        if (!std::isnan(zErrorOuter)) flatPartTable->setContent(9 + extraLine, iRing+1, zErrorOuter/Units::mm, c_tiltedCoordPrecision);
        else                          flatPartTable->setContent(9 + extraLine, iRing+1, "n/a");
      }
      iRing++;
    }

    iRing = 0;
    for (const auto& m : plusBigDeltaModules) {
      flatPartTable->setContent(2, 0  , "R"+web_subStart + "centre" + web_subEnd + "-outer [mm] :");
      flatPartTable->setContent(2, iRing+1, m.center().Rho()/Units::mm, c_tiltedCoordPrecision);
      iRing++;
    }
    m_flatPartTables.push_back(std::move(flatPartTable));
  } // end of 'fills flat part table'
}
