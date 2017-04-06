/*
 * AnalyzerOccupancy.cc
 *
 *  Created on: 9. 11. 2015
 *      Author: Drasal (CERN)
 */

#include <AnalyzerOccupancy.h>

#include <memory>
#include <global_constants.h>

#include <BeamPipe.h>
#include <BFieldMap.h>
#include <Disk.h>
#include <IrradiationMap.h>
#include <Layer.h>
#include "MainConfigHandler.h"
#include <Ring.h>
#include "RootWContent.h"
#include "RootWImage.h"
#include "RootWPage.h"
#include "RootWSite.h"
#include "RootWTable.h"
#include "SimParms.h"
#include <Tracker.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <Units.h>
#include <TLegend.h>

#include "TList.h"

AnalyzerOccupancy::AnalyzerOccupancy(const Detector& detector) :
 AnalyzerUnit("AnalyzerOccupancy", detector)
{
  // Read data from files to memory
  m_photonsMap    = nullptr;
  m_chargedMap    = nullptr;
  m_bFieldMap     = nullptr;
  m_hisChargedFlux= nullptr;
  m_hisPhotonsFlux= nullptr;
}

AnalyzerOccupancy::~AnalyzerOccupancy()
{
  if (m_bFieldMap !=nullptr) delete m_bFieldMap;
  if (m_photonsMap!=nullptr) delete m_photonsMap;
  if (m_chargedMap!=nullptr) delete m_chargedMap;
}

//! Init variables
bool AnalyzerOccupancy::init(int nTracks)
{
  std::string directory = MainConfigHandler::getInstance().getIrradiationDirectory();

  bool bFieldMapOK  = checkFile(SimParms::getInstance().bFieldMapFile(), directory);
  bool chargedMapOK = checkFile(SimParms::getInstance().chargedMapFile(), directory);
  bool photonsMapOK = checkFile(SimParms::getInstance().photonsMapFile(), directory);

  if (bFieldMapOK)  m_bFieldMap   = new BFieldMap(directory+"/"+SimParms::getInstance().bFieldMapFile());
  if (chargedMapOK) m_chargedMap  = new IrradiationMap(directory + "/" + SimParms::getInstance().chargedMapFile());
  if (photonsMapOK) m_photonsMap  = new IrradiationMap(directory + "/" + SimParms::getInstance().photonsMapFile());

  m_isInitOK = bFieldMapOK && chargedMapOK && photonsMapOK;

  return m_isInitOK;
}

bool AnalyzerOccupancy::analyze()
{
  // Check that initialization OK
  if (!m_isInitOK) return false;

  // Make & fill all flux histograms
  fillHistogram(m_chargedMap,   m_hisChargedFlux,   "ChargedFluxPerPP",   "Flux of charged particles [cm^{-2}] per pp collision");
  //fillHistogram(m_photonsMap,   m_hisPhotonsFlux,   "PhotonsFluxPerPP",   "Flux of photons [cm^{-2}] per pp collision");

  m_isAnalysisOK = true;
  return m_isAnalysisOK;
}

bool AnalyzerOccupancy::visualize(RootWSite& webSite)
{
  // Check that initialization & analysis OK
  if (!m_isInitOK || !m_isAnalysisOK) return false;

  RootWPage& myPage = webSite.addPage("Occupancy", web_priority_Occup);
  myPage.setAddress("indexOccupancy.html");


  // Draw magnetic fiel - map
  if (m_bFieldMap!=nullptr && m_bFieldMap->isOK()) {

    RootWContent& magFieldContent = myPage.addContent("Magnetic field map:", true);

    TCanvas canvasXZBField("canvasXZBField", "XZ view of B field [T] (Y=0)", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    TCanvas canvasYZBField("canvasYZBField", "YZ view of B field [T] (X=0)", vis_std_canvas_sizeX, vis_min_canvas_sizeY);

    double geomMaxRadius = 0.0;
    double geomMaxLength = 0.0;

    for (auto& iTrk : m_trackers) {

      auto maxR = iTrk->maxR();
      auto maxZ = iTrk->maxZ();

      if (maxR>geomMaxRadius) geomMaxRadius = maxR;
      if (maxZ>geomMaxLength) geomMaxLength = maxZ;
    }

    m_bFieldMap->drawXZBFieldProj(canvasXZBField, "XZ view of B field [T] (Y=0)", 0, geomMaxRadius, 0, geomMaxLength);
    m_bFieldMap->drawYZBFieldProj(canvasYZBField, "YZ view of B field [T] (X=0)", 0, geomMaxRadius, 0, geomMaxLength);

    RootWImage& anImageXZBField = magFieldContent.addImage(canvasXZBField, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    anImageXZBField.setComment("XZ view of B field [T] (Y=0)");
    RootWImage& anImageYZBField = magFieldContent.addImage(canvasYZBField, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    anImageYZBField.setComment("YZ view of B field [T] (X=0)");
  }

//  // Draw plots - photons
//  RootWContent& plotsPhotonsContent = myPage.addContent("Fluka simulation - photons fluxes per pp collision -> adding individual effects:", false);
//
//  TCanvas canvasPhotons;
//
//  if (drawHistogram(canvasPhotons, m_hisPhotonsFlux, m_photonsMap, "PhotonsCanvas", "RZ view of photons flux")) {
//    canvasPhotons.SetLogz();
//    m_hisPhotonsFlux->SetMinimum(c_fluxMin);
//    m_hisPhotonsFlux->SetMaximum(c_fluxMax);
//    RootWImage& anImagePhotons = plotsPhotonsContent.addImage(canvasPhotons);
//    anImagePhotons.setComment("RZ view of photons flux [cm^-2] in a tracker");
//  }

  // Draw plots - charged
  RootWContent& plotsChargedContent = myPage.addContent("Fluka simulation - charged particles fluxes per pp collision:", true);

  TCanvas canvasCharged;

  if (drawHistogram(canvasCharged, m_hisChargedFlux, m_chargedMap, "ChargedCanvas", "RZ view of charged particles flux")) {
    canvasCharged.SetLogz();
    m_hisChargedFlux->SetMinimum(c_fluxMin);
    m_hisChargedFlux->SetMaximum(c_fluxMax);
    RootWImage& anImageCharged = plotsChargedContent.addImage(canvasCharged);
    anImageCharged.setComment("RZ view of charged particles flux [cm^-2] in a tracker");
  }

  // Plot a table with calculated occupancies
  for (auto itTracker : m_trackers) {

    // Create table
    std::string trkName = itTracker->myid();
    RootWContent& occupancyBarrelContent = myPage.addContent("Occupancy - charged particles ("+trkName+"-barrel)", true);
    RootWContent& occupancyEndcapContent = myPage.addContent("Occupancy - charged particles ("+trkName+"-endcap)", true);

    IrradiationMap* usedChargedMap = m_chargedMap;
    IrradiationMap* usedPhotonsMap = m_photonsMap;

    OccupancyVisitor geometryVisitor(usedPhotonsMap, usedChargedMap);
    geometryVisitor.setMaxPileUp(trk_pile_up[trk_pile_up.size()-1]);
    geometryVisitor.setMaxColFreq(collision_freq);

    // set all values through visitor
    itTracker->accept(geometryVisitor);

    // Print out layer & disk table
    for (auto nPileUps : trk_pile_up) {

      std::ostringstream namePileUp;
      namePileUp << nPileUps;

      std::unique_ptr<RootWTable> layerHeaderTable(new RootWTable());
      layerHeaderTable->setContent(0, 0, "Number of pile-up events: ");
      layerHeaderTable->setContent(0, 1, namePileUp.str());

      occupancyBarrelContent.addItem(std::move(layerHeaderTable));
      occupancyBarrelContent.addItem(std::move(geometryVisitor.getLayerTable(nPileUps, itTracker->myid())));

      std::unique_ptr<RootWTable> ringHeaderTable(new RootWTable());
      ringHeaderTable->setContent(0, 0, "Number of pile-up events: ");
      ringHeaderTable->setContent(0, 1, namePileUp.str());

      occupancyEndcapContent.addItem(std::move(ringHeaderTable));
      occupancyEndcapContent.addItem(std::move(geometryVisitor.getRingTable(nPileUps, itTracker->myid())));
    }
  }

  m_isVisOK = true;
  return m_isVisOK;
}

//
// Check that a file can be opened
//
bool AnalyzerOccupancy::checkFile(const std::string& fileName, const std::string& filePath)
{
  fstream     file;
  std::string fullFileName(filePath+"/"+fileName);
  file.open(fullFileName);
  if (file.is_open()) {

    file.close();
    return true;
  }
  else {

    logERROR("AnalyzerOccupancy - failed opening file: " + fullFileName);
    return false;
  }
}

bool AnalyzerOccupancy::fillHistogram(const IrradiationMap* map, TH2D*& his, std::string name, std::string title)
{
  if (map!=nullptr) {

    // Make & fill flux histograms
    if (map->isOfTypeMesh()) {
      his = new TH2D(name.c_str(), title.c_str(),
                     map->getZNBins(), map->getZMin()-map->getZBinWidth()/2., map->getZMax()+map->getZBinWidth()/2.,
                     map->getRNBins(), map->getRMin()-map->getRBinWidth()/2., map->getRMax()+map->getRBinWidth()/2.);
    }
    else {
      his = new TH2D(name.c_str(), title.c_str(),
                     map->getZNBins(), map->getZMin(), map->getZMax(),
                     map->getRNBins(), map->getRMin(), map->getRMax());
    }

    for (int zBin=0; zBin<map->getZNBins(); zBin++) {
      for (int rBin=0; rBin<map->getRNBins(); rBin++) {

        double zPos = 0;
        double rPos = 0;
        if (map->isOfTypeMesh()) {
          zPos = map->getZMin() + zBin*map->getZBinWidth();
          rPos = map->getRMin() + rBin*map->getRBinWidth();
        }
        else {
          zPos = map->getZMin() + (zBin+1/2.)*map->getZBinWidth();
          rPos = map->getRMin() + (rBin+1/2.)*map->getRBinWidth();
        }

        // Map arranged in the format ZxR (THist binning starts from 1)
        his->SetBinContent(zBin+1, rBin+1, map->calculateIrradiationZR(zPos, rPos)/(1./Units::cm2));
//        his->SetBinContent(rBin+1, zBin+1, map->calculateIrradiationRZ(rPos, zPos)/(1./Units::cm2));
      }
    }
    return true;
  }
  else return false;
}

bool AnalyzerOccupancy::drawHistogram(TCanvas& canvas, TH2D* his, const IrradiationMap* map, std::string name, std::string title)
{
  if (his!=nullptr) {

    std::string canvasName  = name+"Canvas";
    std::string canvasTitle = title;
    canvas.SetName(canvasName.c_str());
    canvas.SetTitle(canvasTitle.c_str());
    canvas.SetWindowSize(vis_std_canvas_sizeX, vis_min_canvas_sizeY);

    canvas.cd();
    his->Draw("COLZ");
    his->SetStats(kFALSE);
    his->GetXaxis()->SetTitle(std::string("Z ["+map->getZUnit()+"]").c_str());
    his->GetXaxis()->SetTitleOffset(1.2);
    his->GetYaxis()->SetTitle(std::string("R ["+map->getRUnit()+"]").c_str());
    his->GetYaxis()->SetTitleOffset(1.2);
    his->GetYaxis()->SetRangeUser(m_beamPipe->radius(), his->GetYaxis()->GetXmax());
    return true;
  }
  else return false;
}

//
// Helper occupancy visitor constructor
//
OccupancyVisitor::OccupancyVisitor(IrradiationMap* photonsMap, IrradiationMap* chargedMap) {

  m_photonsMap = photonsMap;
  m_chargedMap = chargedMap;
  m_nLayers    = 0;
  m_nDisks     = 0;
  m_iRing      = 0;
  m_nRings     = 0;
  m_zPosStep   = std::min(m_chargedMap->getZBinWidth(),m_photonsMap->getZBinWidth());
  m_rPosStep   = std::min(m_chargedMap->getRBinWidth(),m_photonsMap->getRBinWidth());

  m_maxPileUp  = 1.;
  m_maxColFreq = 1.;
}

//
// Visit layer
//
void OccupancyVisitor::visit(const Layer& layer) {

  if (layer.maxZ() < 0.) return;

  double minFlux = std::numeric_limits<double>::max();
  double maxFlux = 0;
  double zPos    = 0.;
  double maxZPos = 0.;

  while (zPos<=layer.maxZ()) {

    double flux  = m_chargedMap->calculateIrradiationZR(zPos, layer.avgBuildRadius());
//    double flux  = m_chargedMap->calculateIrradiationRZ(layer.avgBuildRadius(), zPos);

     if (flux>maxFlux) {

       maxFlux = flux;
       maxZPos = zPos;
     }
     if (flux<minFlux) minFlux = flux;
     zPos += m_zPosStep;
   }

   std::vector<long> vecNHits, vecNChannels;
   m_layerNHits.push_back(vecNHits);
   m_layerNChannels.push_back(vecNChannels);
   m_layerRadii.push_back(layer.avgBuildRadius());
   m_layerMinFluxes.push_back(minFlux);
   m_layerMaxFluxes.push_back(maxFlux);
   m_layerMaxFluxZ.push_back(maxZPos);
   m_layerNRods.push_back(layer.numRods());
   m_layerNModules.push_back(0);
   std::vector<int> vecAddrSpar, vecAddrUnspar, vecNPixels;
   std::vector<double> vecPixelArea, vecSenArea;
   m_layerSenAddrSparSize.push_back(vecAddrSpar);
   m_layerSenAddrUnsparSize.push_back(vecAddrUnspar);
   m_layerSenNPixels.push_back(vecNPixels);
   m_layerSenPixelArea.push_back(vecPixelArea);
   m_layerSenArea.push_back(vecSenArea);
   m_layerNSensorsInMod.push_back(0);

   m_nLayers++;
 }

//
// Visit barrel module
//
 void OccupancyVisitor::visit(const BarrelModule& module) {

   int    iLayer   = m_nLayers-1;
   double zPos     = fabs((module.planarMaxZ()+module.planarMinZ())/2.);
   double rPos     = (module.planarMaxR()+module.planarMinR())/2.;
   long   nHits    = module.area() * m_chargedMap->calculateIrradiationZR(zPos, rPos)/Units::mm2 * m_maxPileUp;
   // long   nHits    = module.area() * m_chargedMap->calculateIrradiationRZ(rPos, zPos)/Units::mm2 * m_maxPileUp;

   if (m_layerSenArea[iLayer].size()==0) m_layerSenArea[iLayer].push_back(module.area());

   short iSensor = 0;
   for (const auto& s : module.sensors()) {

     int nSegments         = s.numSegments();
     int nStrips           = s.numStripsAcross();
     SensorType sensorType = s.type();

     if (iSensor>=2) {
       logWARNING("Occupancy studies - module contains more than 2 sensors -> check the code, data rate will be wrong!");
       continue;
     }
     else if (m_layerNSensorsInMod[iLayer]==0) m_layerNSensorsInMod[iLayer] = module.numSensors();

     // Allocate memory if needed - 2 sensors per module in maximum
     if (m_layerNChannels[iLayer].size()==0         && iSensor==0) m_layerNChannels[iLayer].push_back(0);
     if (m_layerNChannels[iLayer].size()==1         && iSensor==1) m_layerNChannels[iLayer].push_back(0);
     if (m_layerNHits[iLayer].size()==0             && iSensor==0) m_layerNHits[iLayer].push_back(0);
     if (m_layerNHits[iLayer].size()==1             && iSensor==1) m_layerNHits[iLayer].push_back(0);
     if (m_layerSenAddrSparSize[iLayer].size()==0   && iSensor==0) m_layerSenAddrSparSize[iLayer].push_back(0);
     if (m_layerSenAddrSparSize[iLayer].size()==1   && iSensor==1) m_layerSenAddrSparSize[iLayer].push_back(0);
     if (m_layerSenAddrUnsparSize[iLayer].size()==0 && iSensor==0) m_layerSenAddrUnsparSize[iLayer].push_back(0);
     if (m_layerSenAddrUnsparSize[iLayer].size()==1 && iSensor==1) m_layerSenAddrUnsparSize[iLayer].push_back(0);
     if (m_layerSenNPixels[iLayer].size()==0        && iSensor==0) m_layerSenNPixels[iLayer].push_back(0);
     if (m_layerSenNPixels[iLayer].size()==1        && iSensor==1) m_layerSenNPixels[iLayer].push_back(0);
     if (m_layerSenPixelArea[iLayer].size()==0      && iSensor==0) m_layerSenPixelArea[iLayer].push_back(0);
     if (m_layerSenPixelArea[iLayer].size()==1      && iSensor==1) m_layerSenPixelArea[iLayer].push_back(0);

     // This separation works only for true crosssing strips -> not for several times repeating strip-lets (segments)
     //if (sensorType==SensorType::Largepix || sensorType==SensorType::Pixel) {
       long nReadOutChannels = 0;
       if (nHits<=nSegments*nStrips) nReadOutChannels = nHits;
       else                          nReadOutChannels = nSegments*nStrips;

       m_layerNChannels[iLayer][iSensor] += nReadOutChannels;
       m_layerNHits[iLayer][iSensor]     += nHits;
     /*}
     else if (sensorType==SensorType::Strip) {
       long nReadOutChannelsX = 0;
       long nReadOutChannelsY = 0;

       if (nHits<=nSegments) nReadOutChannelsX = nHits;
       else                  nReadOutChannelsX = nSegments;
       if (nHits<=nStrips)   nReadOutChannelsY = nHits;
       else                  nReadOutChannelsY = nStrips;
       m_layerNChannels[iLayer][iSensor] += nReadOutChannelsX+nReadOutChannelsY;
       m_layerNHits[iLayer][iSensor]     += nHits;
     }
     else {
       logWARNING("Occupancy studies: Sensor type not define -> couldn't calculate the data rate");
     }*/

    int senAddrSpar   = std::ceil(log2(nSegments)) + std::ceil(log2(nStrips)) * Units::b;
    int senAddrUnspar = nSegments*nStrips * Units::b;
    if (m_layerSenAddrSparSize[iLayer][iSensor]!=0) {
      if (m_layerSenAddrSparSize[iLayer][iSensor]!=senAddrSpar) logWARNING("Occupancy studies - module types differ within a layer -> check the code, data rate will be wrong!");
    }
    else {
      m_layerSenAddrSparSize[iLayer][iSensor] = senAddrSpar;
      //std::cout << ">>Spar> " << module.numSensors() << " " << m_layerSenAddrSparSize[iLayer][iSensor] << std::endl;
    }
    if (m_layerSenAddrUnsparSize[iLayer][iSensor]!=0) {
      if (m_layerSenAddrUnsparSize[iLayer][iSensor]!=senAddrUnspar) logWARNING("Occupancy studies - module types differ within a layer -> check the code, data rate will be wrong!");
    }
    else {
      m_layerSenAddrUnsparSize[iLayer][iSensor] = senAddrUnspar;
      //std::cout << ">>UnS> " << module.numSensors() << " " << m_layerSenAddrUnsparSize[iLayer][iSensor] << std::endl;
    }
    if (m_layerSenNPixels[iLayer][iSensor]!=0) {
      if ((m_layerSenNPixels[iLayer][iSensor]!=nStrips+nSegments) && (m_layerSenNPixels[iLayer][iSensor]!=nStrips*nSegments)) logWARNING("Occupancy studies - module types differ within a layer -> check the code, data rate will be wrong!");
    }
    else {
      //if (sensorType==SensorType::Strip) m_layerSenNPixels[iLayer][iSensor] = nStrips+nSegments;
      //else                               m_layerSenNPixels[iLayer][iSensor] = nStrips*nSegments;
      m_layerSenNPixels[iLayer][iSensor] = nStrips*nSegments;
    }
    if (m_layerSenPixelArea[iLayer][iSensor]!=0) {
      if (m_layerSenPixelArea[iLayer][iSensor]!=module.resolutionLocalX()*module.resolutionLocalY()*12) logWARNING("Occupancy studies - module types differ within a layer -> check the code, data rate will be wrong!");
    }
    else {
      m_layerSenPixelArea[iLayer][iSensor] = module.resolutionLocalX()*module.resolutionLocalY()*12;
    }

    iSensor++;
  }

  m_layerNModules[iLayer]++;
}

//
// Visit disk
//
void OccupancyVisitor::visit(const Disk& disk) {

  if (disk.averageZ() < 0.) return;

  // Update disk counter
  if (disk.numRings()>m_nRings) {

    m_ringAvgRadii.resize(disk.numRings());
    m_ringMinFluxes.resize(disk.numRings());
    m_ringMaxFluxes.resize(disk.numRings());
    m_ringMaxFluxZ.resize(disk.numRings());
    m_ringNChannels.resize(disk.numRings());
    m_ringNHits.resize(disk.numRings());
    m_ringNModules.resize(disk.numRings());
    m_ringNSensorsInMod.resize(disk.numRings());
    m_ringSenAddrSparSize.resize(disk.numRings());
    m_ringSenAddrUnsparSize.resize(disk.numRings());
    m_ringSenNPixels.resize(disk.numRings());
    m_ringSenPixelArea.resize(disk.numRings());
    m_ringSenArea.resize(disk.numRings());

    // Initialize min & max flux values
    for (int i=m_nRings; i<disk.numRings(); i++) {
      m_ringMinFluxes[i] = std::numeric_limits<double>::max();
      m_ringMaxFluxes[i] = 0;
    }
  }
  m_nRings = disk.numRings();
  ++m_nDisks;
}

//
// Visit ring
//
void OccupancyVisitor::visit(const Ring& ring) {

  if (ring.averageZ()<0.) return;

  double minFlux = std::numeric_limits<double>::max();
  double maxFlux = 0;
  double rPos    = ring.minR();

  m_iRing = ring.myid()-1;

  while (rPos<=ring.maxR()) {

    double flux  = m_chargedMap->calculateIrradiationZR(ring.averageZ(), rPos);//*cosTheta;
    // double flux  = m_chargedMap->calculateIrradiationRZ(rPos, ring.averageZ());


    if (flux>maxFlux) maxFlux = flux;
    if (flux<minFlux) minFlux = flux;
    rPos += m_rPosStep;
  }

  m_ringAvgRadii[m_iRing] = (ring.minR()+ring.maxR())/2.;
  double newMinFlux = std::min(minFlux, m_ringMinFluxes[m_iRing]);
  double newMaxFlux = std::max(maxFlux, m_ringMaxFluxes[m_iRing]);

  if (newMinFlux==minFlux) m_ringMinFluxes[m_iRing] = newMinFlux;
  if (newMaxFlux==maxFlux) {

    m_ringMaxFluxes[m_iRing] = newMaxFlux;
    m_ringMaxFluxZ[m_iRing]  = ring.averageZ();
  }
}

void OccupancyVisitor::visit(const EndcapModule& module) {

  static int countModPls = 0;
  if ((module.planarMaxZ()+module.planarMinZ())/2.<0) return;

  int    iDisk    = m_nDisks-1;
  double zPos     = fabs((module.planarMaxZ()+module.planarMinZ())/2.);
  double rPos     = (module.planarMaxR()+module.planarMinR())/2.;
  long   nHits    = module.area() * m_chargedMap->calculateIrradiationZR(zPos, rPos)/Units::mm2 * m_maxPileUp;
  //  long   nHits    = module.area() * m_chargedMap->calculateIrradiationRZ(rPos, zPos)/Units::mm2 * m_maxPileUp;

  if (m_ringSenArea[m_iRing].size()==0) m_ringSenArea[m_iRing].push_back(module.area());

  short iSensor = 0;
  for (const auto& s : module.sensors()) {

    int nSegments         = s.numSegments();
    int nStrips           = s.numStripsAcross();
    SensorType sensorType = s.type();

    if (iSensor>=2) {
      logWARNING("Occupancy studies - module contains more than 2 sensors -> check the code, data rate will be wrong!");
      continue;
    }
    else m_ringNSensorsInMod[m_iRing] = module.numSensors();

    // Allocate memory if needed - 2 sensors per module in maximum
    if (iSensor==0 && m_ringNChannels[m_iRing].size()<1) {
      m_ringNChannels[m_iRing].resize(1,0);
      m_ringNHits[m_iRing].resize(1,0);
      m_ringSenAddrSparSize[m_iRing].resize(1,0);
      m_ringSenAddrUnsparSize[m_iRing].resize(1,0);
      m_ringSenNPixels[m_iRing].resize(1,0);
      m_ringSenPixelArea[m_iRing].resize(1,0);
    }
    else if (iSensor==1 && m_ringNChannels[m_iRing].size()<2) {
      m_ringNChannels[m_iRing].resize(2,0);
      m_ringNHits[m_iRing].resize(2,0);
      m_ringSenAddrSparSize[m_iRing].resize(2,0);
      m_ringSenAddrUnsparSize[m_iRing].resize(2,0);
      m_ringSenNPixels[m_iRing].resize(2,0);
      m_ringSenPixelArea[m_iRing].resize(2,0);
    }

    // This separation works only for true crosssing strips -> not for several times repeating strip-lets (segments)
    //if (sensorType==SensorType::Largepix || sensorType==SensorType::Pixel) {
    long nReadOutChannels = 0;
    if (nHits<=nSegments*nStrips) nReadOutChannels = nHits;
    else                          nReadOutChannels = nSegments*nStrips;

    m_ringNChannels[m_iRing][iSensor] += nReadOutChannels;
    m_ringNHits[m_iRing][iSensor]     += nHits;
    // }
    // else if (sensorType==SensorType::Strip) {
    // long nReadOutChannelsX = 0;
    // long nReadOutChannelsY = 0;
    // if (nHits<=nSegments) nReadOutChannelsX = nHits;
    // else                  nReadOutChannelsX = nSegments;
    // if (nHits<=nStrips)   nReadOutChannelsY = nHits;
    // else                  nReadOutChannelsY = nStrips;
    // m_ringNChannels[m_iRing][iSensor] += nReadOutChannelsX+nReadOutChannelsY;
    // m_ringNHits[m_iRing][iSensor]     += nHits;
    // }
    // else {
    // logWARNING("Occupancy studies: Sensor type not define -> couldn't calculate the data rate");
    // }

    int senAddrSpar   = std::ceil(log2(nSegments)) + std::ceil(log2(nStrips)) * Units::b;
    int senAddrUnspar = nSegments*nStrips * Units::b;

    if (m_ringSenAddrSparSize[m_iRing][iSensor]!=0 && m_ringSenAddrSparSize[m_iRing][iSensor]!=senAddrSpar) logWARNING("Occupancy studies - module types differ within a ring -> check the code, data rate will be wrong!");
    else m_ringSenAddrSparSize[m_iRing][iSensor] = senAddrSpar;

    if (m_ringSenAddrUnsparSize[m_iRing][iSensor]!=0 && m_ringSenAddrUnsparSize[m_iRing][iSensor]!=senAddrUnspar) logWARNING("Occupancy studies - module types differ within a ring -> check the code, data rate will be wrong!");
    else m_ringSenAddrUnsparSize[m_iRing][iSensor] = senAddrUnspar;

    if (m_ringSenNPixels[m_iRing][iSensor]!=0 && m_ringSenNPixels[m_iRing][iSensor]!=nStrips+nSegments && m_ringSenNPixels[m_iRing][iSensor]!=nStrips*nSegments) logWARNING("Occupancy studies - module types differ within a ring -> check the code, data rate will be wrong!");
    else m_ringSenNPixels[m_iRing][iSensor] = nStrips*nSegments;

    if (m_ringSenPixelArea[m_iRing][iSensor]!=0 && m_ringSenPixelArea[m_iRing][iSensor]!=module.resolutionLocalX()*module.resolutionLocalY()*12) logWARNING("Occupancy studies - module types differ within a layer -> check the code, data rate will be wrong!");
    else m_ringSenPixelArea[m_iRing][iSensor] = module.resolutionLocalX()*module.resolutionLocalY()*12;

    iSensor++;

  }

  m_ringNModules[m_iRing]++;
}

//
// Get created layer summary table
//
std::unique_ptr<RootWTable> OccupancyVisitor::getLayerTable(signed int nPileUps, std::string trkName) {


  std::unique_ptr<RootWTable> layerTable(new RootWTable());

  double precisionFlux      = 2*c_coordPrecision;
  double precisionArea      = 2*c_coordPrecision;
  double precisionOccupancy = 2*c_coordPrecision;
  if (trkName=="Inner") {

    precisionFlux = 1*c_coordPrecision;
    precisionArea = 4*c_coordPrecision;
  }
  else {

    precisionFlux = 2*c_coordPrecision;
    precisionArea = 3*c_coordPrecision;
  }

  double totDataRateTriggerSpar   = 0;
  double totDataRateUnTriggerSpar = 0;

  for (int iLayer=0; iLayer<m_nLayers; iLayer++) {

    double minFlux     = m_layerMinFluxes[iLayer]*nPileUps;
    double maxFlux     = m_layerMaxFluxes[iLayer]*nPileUps;
    double maxCellArea = trk_max_occupancy/maxFlux;
    double numSensors  = m_layerNSensorsInMod[iLayer];

    // Calculate data rates for layers
    std::vector<double> hitRate;
    std::vector<double> channelRate;
    std::vector<int>    senAddrSparSize;
    std::vector<int>    senAddrUnsparSize;
    hitRate.resize(numSensors);
    channelRate.resize(numSensors);
    senAddrSparSize.resize(numSensors);
    senAddrUnsparSize.resize(numSensors);

    double totHitRate        = 0;
    double totChannelRate    = 0;
    int    totAddrSparSize   = 0;
    int    totAddrUnsparSize = 0;
    int    addrSparClsWidth  = 2 * Units::b;

    double dataRateCollisionSpar = 0;
    double dataRateTriggerSpar   = 0;
    double dataRateUnTriggerSpar = 0;
    double moduleOccupancy       = 0;
    double pixelArea             = 0;
    double moduleArea            = m_layerSenArea[iLayer][0];

    for (int iSensor=0; iSensor<numSensors; iSensor++) {
      hitRate[iSensor]           = m_layerNHits[iLayer][iSensor];
      channelRate[iSensor]       = m_layerNChannels[iLayer][iSensor];
      senAddrSparSize[iSensor]   = m_layerSenAddrSparSize[iLayer][iSensor] + addrSparClsWidth;
      senAddrUnsparSize[iSensor] = m_layerSenAddrUnsparSize[iLayer][iSensor];

      totHitRate        += hitRate[iSensor];
      totChannelRate    += channelRate[iSensor];
      totAddrSparSize   += senAddrSparSize[iSensor];
      totAddrUnsparSize += senAddrUnsparSize[iSensor];

      moduleOccupancy = channelRate[iSensor]/m_layerSenNPixels[iLayer][iSensor]/m_layerNModules[iLayer];

      pixelArea = m_layerSenPixelArea[iLayer][iSensor];

      dataRateCollisionSpar += channelRate[iSensor]*senAddrSparSize[iSensor];
      dataRateTriggerSpar   += channelRate[iSensor]*senAddrSparSize[iSensor];
      dataRateUnTriggerSpar += channelRate[iSensor]*senAddrSparSize[iSensor];
    }

    dataRateTriggerSpar   *= trigger_freq;
    dataRateUnTriggerSpar *= collision_freq;

    totDataRateTriggerSpar   += dataRateTriggerSpar;
    totDataRateUnTriggerSpar += dataRateUnTriggerSpar;

    // Layer table
    layerTable->setContent(0, 0, "Layer no                                   : ");
    layerTable->setContent(1, 0, "Radius [mm]                                : ");
    layerTable->setContent(2, 0, "Min flux in Z [particles/cm^2]             : ");
    layerTable->setContent(3, 0, "Max flux in Z [particles/cm^2]             : ");
    layerTable->setContent(4, 0, "Z position [mm] related to max flux        : ");
    layerTable->setContent(5, 0, "Max cell area (1% occupancy) [mm^2]        : ");
    layerTable->setContent(6, 0, "Module max occupancy (max[sen1,sen2])[%]   : ");
    if (nPileUps==trk_pile_up[trk_pile_up.size()-1]) {
      layerTable->setContent(7 , 0, "#Hits per BX (bunch crossing)             : ");
      layerTable->setContent(8 , 0, "#Hit-channels per BX                      : ");
      layerTable->setContent(9 , 0, "#Hit-channels per module per BX           : ");
      layerTable->setContent(10, 0, "Module avg occupancy (max[sen1,sen2])[%]  : ");
      layerTable->setContent(11, 0, "Module bandwidth/(addr+clsWidth=2b[b]     : ");
      layerTable->setContent(12, 0, "Mod. bandwidth(#chnls*(addr+clsWidth)[kb] : ");
      layerTable->setContent(13, 0, "Mod. bandwidth (matrix*1b/channel) [kb]   : ");
      layerTable->setContent(14, 0, "Data rate per layer - 40MHz,spars [Tb/s]  : ");
      layerTable->setContent(15, 0, "Data rate per layer -  1MHz,spars [Tb/s]  : ");
      layerTable->setContent(16, 0, "Data rate per ladder - 40Mhz,spars [Gb/s] : ");
      layerTable->setContent(17, 0, "Data rate per ladder -  1Mhz,spars [Gb/s] : ");
      layerTable->setContent(18, 0, "Data rate per module - 40Mhz,spars [Gb/s] : ");
      layerTable->setContent(19, 0, "Data rate per module -  1Mhz,spars [Gb/s] : ");
      layerTable->setContent(20, 0, "<b>Data rate per cm^2 - 40Mhz,spars [Gb/s/cm^2]</b>: ");
      layerTable->setContent(21, 0, "<b>Data rate per cm^2 -  1Mhz,spars [Gb/s/cm^2]</b>: ");
    }

    layerTable->setContent(0, iLayer+1, iLayer+1);
    layerTable->setContent(1, iLayer+1, m_layerRadii[iLayer]/Units::mm   , c_coordPrecision);
    layerTable->setContent(2, iLayer+1, minFlux/(1./Units::cm2)          , precisionFlux);
    layerTable->setContent(3, iLayer+1, maxFlux/(1./Units::cm2)          , precisionFlux);
    layerTable->setContent(4, iLayer+1, m_layerMaxFluxZ[iLayer]/Units::mm, c_coordPrecision);
    layerTable->setContent(5, iLayer+1, maxCellArea/Units::mm2           , precisionArea);
    layerTable->setContent(6, iLayer+1, maxFlux*pixelArea*100            , precisionOccupancy);
    if (nPileUps==trk_pile_up[trk_pile_up.size()-1]) {
      layerTable->setContent(7 , iLayer+1, totHitRate                          );
      layerTable->setContent(8 , iLayer+1, totChannelRate                      );
      layerTable->setContent(9 , iLayer+1, totChannelRate/m_layerNModules[iLayer]);
      layerTable->setContent(10, iLayer+1, moduleOccupancy*100, precisionOccupancy); // maxFlux*pixelArea*100, precisionOccupancy);
      layerTable->setContent(11, iLayer+1, totAddrSparSize/Units::b);
      layerTable->setContent(12, iLayer+1, dataRateCollisionSpar/m_layerNModules[iLayer]/Units::kb, 2*c_coordPrecision);
      layerTable->setContent(13, iLayer+1, totAddrUnsparSize/Units::kb                            , 2*c_coordPrecision);
      layerTable->setContent(14, iLayer+1, dataRateUnTriggerSpar/(Units::Tb/Units::s), c_coordPrecision);
      layerTable->setContent(15, iLayer+1, dataRateTriggerSpar/(Units::Tb/Units::s)  , c_coordPrecision);
      layerTable->setContent(16, iLayer+1, dataRateUnTriggerSpar/m_layerNRods[iLayer]/(Units::Gb/Units::s), c_coordPrecision);
      layerTable->setContent(17, iLayer+1, dataRateTriggerSpar/m_layerNRods[iLayer]/(Units::Gb/Units::s)  , c_coordPrecision);
      layerTable->setContent(18, iLayer+1, dataRateUnTriggerSpar/m_layerNModules[iLayer]/(Units::Gb/Units::s), 2*c_coordPrecision);
      layerTable->setContent(19, iLayer+1, dataRateTriggerSpar/m_layerNModules[iLayer]/(Units::Gb/Units::s)  , 2*c_coordPrecision);
      layerTable->setContent(20, iLayer+1, dataRateUnTriggerSpar/m_layerNModules[iLayer]/moduleArea/(Units::Gb/(Units::s*Units::cm2)), 2*c_coordPrecision);
      layerTable->setContent(21, iLayer+1, dataRateTriggerSpar/m_layerNModules[iLayer]/moduleArea/(Units::Gb/(Units::s*Units::cm2))  , 2*c_coordPrecision);
    }
  }
  if (m_nLayers>0 && (nPileUps==trk_pile_up[trk_pile_up.size()-1])) {
    layerTable->setContent(0 , m_nLayers+1, "Total [TB/s]");
    layerTable->setContent(14, m_nLayers+1, totDataRateUnTriggerSpar/(Units::TB/Units::s), c_coordPrecision);
    layerTable->setContent(15, m_nLayers+1, totDataRateTriggerSpar/(Units::TB/Units::s)  , c_coordPrecision);
  }

  return std::move(layerTable);
}

//
// Get created ring summary table
//
std::unique_ptr<RootWTable> OccupancyVisitor::getRingTable(signed int nPileUps, std::string trkName) {

  std::unique_ptr<RootWTable> ringTable(new RootWTable());

  double precisionFlux      = 2*c_coordPrecision;
  double precisionArea      = 2*c_coordPrecision;
  double precisionOccupancy = 2*c_coordPrecision;

  if (trkName=="Inner") {

    precisionFlux = 1*c_coordPrecision;
    precisionArea = 4*c_coordPrecision;
  }
  else {

    precisionFlux = 2*c_coordPrecision;
    precisionArea = 3*c_coordPrecision;
  }

  double totDataRateTriggerSpar   = 0;
  double totDataRateUnTriggerSpar = 0;

  for (int iRing=0; iRing<m_nRings; iRing++) {

    double minFlux      = m_ringMinFluxes[iRing]*nPileUps;
    double maxFlux      = m_ringMaxFluxes[iRing]*nPileUps;
    double maxCellArea  = trk_max_occupancy/maxFlux;
    double numSensors   = m_ringNSensorsInMod[iRing];

    // Calculate data rates for rings
    std::vector<double> hitRate;
    std::vector<double> channelRate;
    std::vector<int>    senAddrSparSize;
    std::vector<int>    senAddrUnsparSize;
    hitRate.resize(numSensors);
    channelRate.resize(numSensors);
    senAddrSparSize.resize(numSensors);
    senAddrUnsparSize.resize(numSensors);

    double totHitRate        = 0;
    double totChannelRate    = 0;
    int    totAddrSparSize   = 0;
    int    totAddrUnsparSize = 0;
    int    addrSparClsWidth  = 2 * Units::b;

    double dataRateCollisionSpar = 0;
    double dataRateTriggerSpar   = 0;
    double dataRateUnTriggerSpar = 0;
    double moduleOccupancy       = 0;
    double pixelArea             = 0;
    double moduleArea            = m_ringSenArea[iRing][0];

    for (int iSensor=0; iSensor<numSensors; iSensor++) {
      hitRate[iSensor]           = m_ringNHits[iRing][iSensor];
      channelRate[iSensor]       = m_ringNChannels[iRing][iSensor];
      senAddrSparSize[iSensor]   = m_ringSenAddrSparSize[iRing][iSensor] + addrSparClsWidth;
      senAddrUnsparSize[iSensor] = m_ringSenAddrUnsparSize[iRing][iSensor];

      totHitRate        += hitRate[iSensor];
      totChannelRate    += channelRate[iSensor];
      totAddrSparSize   += senAddrSparSize[iSensor];
      totAddrUnsparSize += senAddrUnsparSize[iSensor];

      moduleOccupancy = channelRate[iSensor]/m_ringSenNPixels[iRing][iSensor]/m_ringNModules[iRing];

      pixelArea = m_ringSenPixelArea[iRing][iSensor];

      dataRateCollisionSpar += channelRate[iSensor]*senAddrSparSize[iSensor];
      dataRateTriggerSpar   += channelRate[iSensor]*senAddrSparSize[iSensor];
      dataRateUnTriggerSpar += channelRate[iSensor]*senAddrSparSize[iSensor];
    }

    dataRateTriggerSpar   *= trigger_freq;
    dataRateUnTriggerSpar *= collision_freq;

    totDataRateTriggerSpar   += dataRateTriggerSpar;
    totDataRateUnTriggerSpar += dataRateUnTriggerSpar;

    // Ring table
    ringTable->setContent(0, 0, "Ring no                                 : ");
    ringTable->setContent(1, 0, "Average radius [mm]                     : ");
    ringTable->setContent(2, 0, "Min flux in R [particles/cm^2]          : ");
    ringTable->setContent(3, 0, "Max flux in R [particles/cm^2]          : ");
    ringTable->setContent(4, 0, "Z position [mm] related to max flux     : ");
    ringTable->setContent(5, 0, "Max cell area (1% occupancy) [mm^2]     : ");
    ringTable->setContent(6, 0, "Module max occupancy (max[sen1,sen2])[%]: ");
    if (nPileUps==trk_pile_up[trk_pile_up.size()-1]) {
      ringTable->setContent(7 , 0, "#Hits per BX (bunch crossing)             : ");
      ringTable->setContent(8 , 0, "#Hit-channels per BX                      : ");
      ringTable->setContent(9 , 0, "#Hit-channels per module per BX           : ");
      ringTable->setContent(10, 0, "Module avg occupancy (max[sen1,sen2]) [%] : ");
      ringTable->setContent(11, 0, "Module bandwidth/(addr+clsWidth=2b[b]     : ");
      ringTable->setContent(12, 0, "Mod. bandwidth(#chnls*(addr+clsWidth)[kb] : ");
      ringTable->setContent(13, 0, "Mod. bandwidth (matrix*1b/channel) [kb]   : ");
      ringTable->setContent(14, 0, "Data rate per ringLayer-40MHz,spars [Tb/s]: ");
      ringTable->setContent(15, 0, "Data rate per ringLayer- 1MHz,spars [Tb/s]: ");
      ringTable->setContent(16, 0, "Data rate per ring - 40Mhz,spars [Gb/s]   : ");
      ringTable->setContent(17, 0, "Data rate per ring -  1Mhz,spars [Gb/s]   : ");
      ringTable->setContent(18, 0, "Data rate per module - 40Mhz,spars [Gb/s]: ");
      ringTable->setContent(19, 0, "Data rate per module -  1Mhz,spars [Gb/s]: ");
      ringTable->setContent(20, 0, "<b>Data rate per cm^2 - 40Mhz,spars [Gb/s/cm^2]</b>: ");
      ringTable->setContent(21, 0, "<b>Data rate per cm^2 -  1Mhz,spars [Gb/s/cm^2]</b>: ");
    }

    ringTable->setContent(0, iRing+1, iRing+1);
    ringTable->setContent(1, iRing+1, m_ringAvgRadii[iRing]/Units::mm , c_coordPrecision);
    ringTable->setContent(2, iRing+1, minFlux/(1./Units::cm2)         , precisionFlux);
    ringTable->setContent(3, iRing+1, maxFlux/(1./Units::cm2)         , precisionFlux);
    ringTable->setContent(4, iRing+1, m_ringMaxFluxZ[iRing]/Units::mm , c_coordPrecision);
    ringTable->setContent(5, iRing+1, maxCellArea/Units::mm2          , precisionArea);
    ringTable->setContent(6, iRing+1, maxFlux*pixelArea*100           , precisionOccupancy);

    if (nPileUps==trk_pile_up[trk_pile_up.size()-1]) {
      ringTable->setContent(7 , iRing+1, totHitRate*2                        ); // Factor 2 for positive + negative side (neg. side don't used in calculations)
      ringTable->setContent(8 , iRing+1, totChannelRate*2                    ); // Factor 2 for positive + negative side (neg. side don't used in calculations)
      ringTable->setContent(9 , iRing+1, totChannelRate/m_ringNModules[iRing]);
      ringTable->setContent(10, iRing+1, moduleOccupancy*100, precisionOccupancy); // maxFlux*pixelArea*100, precisionOccupancy);
      ringTable->setContent(11, iRing+1, totAddrSparSize/Units::b);
      ringTable->setContent(12, iRing+1, dataRateCollisionSpar/m_ringNModules[iRing]/Units::kb, 2*c_coordPrecision);
      ringTable->setContent(13, iRing+1, totAddrUnsparSize/Units::kb                          , 2*c_coordPrecision);
      ringTable->setContent(14, iRing+1, dataRateUnTriggerSpar/(Units::Tb/Units::s)*2, c_coordPrecision); // Factor 2 for positive + negative side (neg. side don't used in calculations)
      ringTable->setContent(15, iRing+1, dataRateTriggerSpar/(Units::Tb/Units::s)*2  , c_coordPrecision); // Factor 2 for positive + negative side (neg. side don't used in calculations)
      ringTable->setContent(16, iRing+1, dataRateUnTriggerSpar/m_nDisks/(Units::Gb/Units::s), c_coordPrecision);
      ringTable->setContent(17, iRing+1, dataRateTriggerSpar/m_nDisks/(Units::Gb/Units::s)  , c_coordPrecision);
      ringTable->setContent(18, iRing+1, dataRateUnTriggerSpar/m_ringNModules[iRing]/(Units::Gb/Units::s), 2*c_coordPrecision);
      ringTable->setContent(19, iRing+1, dataRateTriggerSpar/m_ringNModules[iRing]/(Units::Gb/Units::s)  , 2*c_coordPrecision);
      ringTable->setContent(20, iRing+1, dataRateUnTriggerSpar/m_ringNModules[iRing]/moduleArea/(Units::Gb/(Units::s*Units::cm2)), 2*c_coordPrecision);
      ringTable->setContent(21, iRing+1, dataRateTriggerSpar/m_ringNModules[iRing]/moduleArea/(Units::Gb/(Units::s*Units::cm2))  , 2*c_coordPrecision);
    }
    if (m_nRings>0 && (nPileUps==trk_pile_up[trk_pile_up.size()-1])) {
      ringTable->setContent(0 , m_nRings+1, "Total [TB/s]");
      ringTable->setContent(14, m_nRings+1, totDataRateUnTriggerSpar/(Units::TB/Units::s)*2, c_coordPrecision); // Factor 2 for positive + negative side (neg. side don't used in calculations)
      ringTable->setContent(15, m_nRings+1, totDataRateTriggerSpar/(Units::TB/Units::s)*2  , c_coordPrecision); // Factor 2 for positive + negative side (neg. side don't used in calculations)
    }
  }
  return ringTable;
}
