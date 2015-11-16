/*
 * AnalyzerOccupancy.cc
 *
 *  Created on: 9. 11. 2015
 *      Author: Drasal (CERN)
 */

#include <AnalyzerOccupancy.h>
#include <rootweb.hh>
#include <Tracker.h>
#include <IrradiationMap.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <global_constants.h>
#include <Units.h>
#include <SimParms.h>

using namespace insur;

AnalyzerOccupancy::AnalyzerOccupancy(std::string chargedFileName, std::string photonsFileName, std::vector<Tracker*> trackers)
{
  // Set geometry, i.e. individual trackers
  for (auto it : trackers) m_trackers.push_back(it);

  // Read data from files to memory
  m_photonsMap = new IrradiationMap(photonsFileName);
  m_chargedMap = new IrradiationMap(chargedFileName);

  m_hisChargedFlux = nullptr;
  m_hisPhotonsFlux = nullptr;
}

AnalyzerOccupancy::~AnalyzerOccupancy()
{
  delete m_photonsMap;
  delete m_chargedMap;
}

bool AnalyzerOccupancy::calculate(double etaStep)
{
  // Make & fill flux histograms
  if (m_chargedMap->isOfTypeMesh()) {
    m_hisChargedFlux = new TH2D("ChargedFluxPerPP", "Flux of charged particles [cm^{-2}] per pp collision",
                                m_chargedMap->getZNBins(), m_chargedMap->getZMin()-m_chargedMap->getZBinWidth()/2., m_chargedMap->getZMax()+m_chargedMap->getZBinWidth()/2.,
                                m_chargedMap->getRNBins(), m_chargedMap->getRMin()-m_chargedMap->getRBinWidth()/2., m_chargedMap->getRMax()+m_chargedMap->getRBinWidth()/2.);
  }
  else {
    m_hisChargedFlux = new TH2D("ChargedFluxPerPP", "Flux of charged particles [cm^{-2}] per pp collision",
                                m_chargedMap->getZNBins(), m_chargedMap->getZMin(), m_chargedMap->getZMax(),
                                m_chargedMap->getRNBins(), m_chargedMap->getRMin(), m_chargedMap->getRMax());
  }
  if (m_photonsMap->isOfTypeMesh()) {
    m_hisPhotonsFlux = new TH2D("PhotonsFluxPerPP", "Flux of photons [cm^{-2}] per pp collision",
                                m_photonsMap->getZNBins(), m_photonsMap->getZMin()-m_photonsMap->getZBinWidth()/2., m_photonsMap->getZMax()+m_photonsMap->getZBinWidth()/2.,
                                m_photonsMap->getRNBins(), m_photonsMap->getRMin()-m_photonsMap->getRBinWidth()/2., m_photonsMap->getRMax()+m_photonsMap->getRBinWidth()/2.);
  }
  else {
    m_hisPhotonsFlux = new TH2D("PhotonsFluxPerPP", "Flux of photons [cm^{-2}] per pp collision",
                                m_photonsMap->getZNBins(), m_photonsMap->getZMin(), m_photonsMap->getZMax(),
                                m_photonsMap->getRNBins(), m_photonsMap->getRMin(), m_photonsMap->getRMax());
  }

  for (int zBin=0; zBin<m_chargedMap->getZNBins(); zBin++) {
    for (int rBin=0; rBin<m_chargedMap->getRNBins(); rBin++) {

      double zPos = 0;
      double rPos = 0;
      if (m_chargedMap->isOfTypeMesh()) {
        zPos = m_chargedMap->getZMin() + zBin*m_chargedMap->getZBinWidth();
        rPos = m_chargedMap->getRMin() + rBin*m_chargedMap->getRBinWidth();
      }
      else {
        zPos = m_chargedMap->getZMin() + (zBin+1/2.)*m_chargedMap->getZBinWidth();
        rPos = m_chargedMap->getRMin() + (rBin+1/2.)*m_chargedMap->getRBinWidth();
      }

      // Map arranged in the format ZxR (THist binning starts from 1)
      m_hisChargedFlux->SetBinContent(zBin+1, rBin+1, m_chargedMap->calculateIrradiationZR(zPos, rPos)/(1./Units::cm2));
    }
  }
  for (int zBin=0; zBin<m_photonsMap->getZNBins(); zBin++) {
    for (int rBin=0; rBin<m_photonsMap->getRNBins(); rBin++) {

      double zPos = 0;
      double rPos = 0;
      if (m_photonsMap->isOfTypeMesh()) {
        zPos = m_photonsMap->getZMin() + zBin*m_photonsMap->getZBinWidth();
        rPos = m_photonsMap->getRMin() + rBin*m_photonsMap->getRBinWidth();
      }
      else {
        zPos = m_photonsMap->getZMin() + (zBin+1/2.)*m_photonsMap->getZBinWidth();
        rPos = m_photonsMap->getRMin() + (rBin+1/2.)*m_photonsMap->getRBinWidth();
      }

      // Map arranged in the format ZxR (THist binning starts from 1)
      m_hisPhotonsFlux->SetBinContent(zBin+1, rBin+1, m_photonsMap->calculateIrradiationZR(zPos, rPos)/(1./Units::cm2));
    }
  }
  return true;
}

bool AnalyzerOccupancy::visualize(RootWSite& webSite, const SimParms* simParms)
{
  RootWPage* myPage = new RootWPage("Occupancy");
  myPage->setAddress("occupancy.html");
  webSite.addPage(myPage);

  // Draw plots
  TCanvas* canvasPhotons = new TCanvas("PhotonsFluxCanvas", "RZ View of photons flux", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
  TCanvas* canvasCharged = new TCanvas("ChargedFluxCanvas", "RZ View of charged flux", vis_std_canvas_sizeX, vis_min_canvas_sizeY);

  canvasPhotons->cd();
  m_hisPhotonsFlux->Draw("COLZ");
  m_hisPhotonsFlux->GetXaxis()->SetTitle(std::string("Z ["+m_photonsMap->getZUnit()+"]").c_str());
  m_hisPhotonsFlux->GetXaxis()->SetTitleOffset(1.2);
  m_hisPhotonsFlux->GetYaxis()->SetTitle(std::string("R ["+m_photonsMap->getRUnit()+"]").c_str());
  m_hisPhotonsFlux->GetYaxis()->SetTitleOffset(1.2);
  m_hisPhotonsFlux->GetYaxis()->SetRangeUser(simParms->bpRadius(), m_hisPhotonsFlux->GetYaxis()->GetXmax());


  canvasCharged->cd();
  m_hisChargedFlux->Draw("COLZ");
  m_hisChargedFlux->GetXaxis()->SetTitle(std::string("Z ["+m_photonsMap->getZUnit()+"]").c_str());
  m_hisChargedFlux->GetXaxis()->SetTitleOffset(1.2);
  m_hisChargedFlux->GetYaxis()->SetTitle(std::string("R ["+m_photonsMap->getRUnit()+"]").c_str());
  m_hisChargedFlux->GetYaxis()->SetTitleOffset(1.2);
  m_hisChargedFlux->GetYaxis()->SetRangeUser(simParms->bpRadius(), m_hisChargedFlux->GetYaxis()->GetXmax());

  RootWContent* plotsContent   = new RootWContent("Fluka simulation - fluxes per pp collision:", true);
  myPage->addContent(plotsContent);

  RootWImage*   anImagePhotons = new RootWImage(canvasPhotons, canvasPhotons->GetWindowWidth(), canvasPhotons->GetWindowHeight());
  anImagePhotons->setComment("RZ view of photons flux [cm^-2] in a tracker)");
  plotsContent->addItem(anImagePhotons);

  RootWImage*   anImageCharged = new RootWImage(canvasCharged, canvasCharged->GetWindowWidth(), canvasCharged->GetWindowHeight());
  anImageCharged->setComment("RZ view of charged particles flux [cm^-2] in a tracker)");
  plotsContent->addItem(anImageCharged);

  // Plot a table with calculated occupancies
  for (auto itTracker : m_trackers) {

    // Create table
    std::string trkName = itTracker->myid();
    RootWContent* occupancyBarrelContent = new RootWContent("Occupancy - photons & charged particles ("+trkName+"-barrel)", false);
    RootWContent* occupancyEndcapContent = new RootWContent("Occupancy - photons & charged particles ("+trkName+"-endcap)", false);
    myPage->addContent(occupancyBarrelContent);
    myPage->addContent(occupancyEndcapContent);

    // Create visitor class & fill tables with data
    class OccupancyVisitor : public ConstGeometryVisitor {
     private:
      IrradiationMap* m_photonsMap;
      IrradiationMap* m_chargedMap;
      int             m_nLayers;
      int             m_nDisks;
      int             m_nRings;

      std::vector<double> m_layerRadii;
      std::vector<double> m_layerMinFluxes;
      std::vector<double> m_layerMaxFluxes;
      std::vector<double> m_layerMaxFluxZ;

      std::vector<double> m_ringAvgRadii;
      std::vector<double> m_ringMinFluxes;
      std::vector<double> m_ringMaxFluxes;
      std::vector<double> m_ringMaxFluxZ;

      double    m_zPosStep;
      double    m_rPosStep;
      const int c_coordPrecision= 1;

     public:
      OccupancyVisitor(IrradiationMap* photonsMap, IrradiationMap* chargedMap) {

        m_photonsMap = photonsMap;
        m_chargedMap = chargedMap;
        m_nLayers    = 0;
        m_nDisks     = 0;
        m_nRings     = 0;
        m_zPosStep   = std::min(m_chargedMap->getZBinWidth(),m_photonsMap->getZBinWidth());
        m_rPosStep   = std::min(m_chargedMap->getRBinWidth(),m_photonsMap->getRBinWidth());
      }

      virtual ~OccupancyVisitor() {};

      void visit(const Layer& layer) override {
        if (layer.maxZ() < 0.) return;

        double minFlux = std::numeric_limits<double>::max();
        double maxFlux = 0;
        double zPos    = 0.;
        double maxZPos = 0.;

        while (zPos<=layer.maxZ()) {

          double flux  = m_chargedMap->calculateIrradiationZR(zPos, layer.placeRadius());
                 flux += m_photonsMap->calculateIrradiationZR(zPos, layer.placeRadius());

          if (flux>maxFlux) {

            maxFlux = flux;
            maxZPos = zPos;
          }
          if (flux<minFlux) minFlux = flux;
          zPos += m_zPosStep;
        }

        m_layerRadii.push_back(layer.placeRadius());
        m_layerMinFluxes.push_back(minFlux);
        m_layerMaxFluxes.push_back(maxFlux);
        m_layerMaxFluxZ.push_back(maxZPos);
        m_nLayers++;
      }

      void visit(const Disk& disk) override {
        if (disk.averageZ() < 0.) return;
        m_nRings = 0;
        ++m_nDisks;
      }

      void visit(const Ring& ring) override {
        if (ring.averageZ()<0.) return;

        double minFlux = std::numeric_limits<double>::max();
        double maxFlux = 0;
        double rPos    = ring.minR();

        while (rPos<=ring.maxR()) {

          double flux  = m_chargedMap->calculateIrradiationZR(ring.averageZ(), rPos);
                 flux += m_photonsMap->calculateIrradiationZR(ring.averageZ(), rPos);
          if (flux>maxFlux) maxFlux = flux;
          if (flux<minFlux) minFlux = flux;
          rPos += m_rPosStep;
        }

        // Find minimum & maximum value across disks
        // First disk
        if (m_nDisks==1) {
          m_ringAvgRadii.push_back( (ring.minR()+ring.maxR())/2. );
          m_ringMinFluxes.push_back(minFlux);
          m_ringMaxFluxes.push_back(maxFlux);
          m_ringMaxFluxZ.push_back(ring.averageZ());
        }
        // Other disks (need to have the same number of rings
        else {
          if (m_nRings<m_ringMinFluxes.size()) {
            double newMinFlux = std::min(minFlux, m_ringMinFluxes[m_nRings]);
            double newMaxFlux = std::max(maxFlux, m_ringMaxFluxes[m_nRings]);
            if (newMinFlux==minFlux) m_ringMinFluxes[m_nRings] = newMinFlux;
            if (newMaxFlux==maxFlux) {

              m_ringMaxFluxes[m_nRings] = newMaxFlux;
              m_ringMaxFluxZ[m_nRings]  = ring.averageZ();
            }
          }
          else logERROR("Occupancy calculation algorithm failed for disks! Expects the same number of rings in all disks across the given tracker!");
        }

        ++m_nRings;
      }

      RootWTable* getLayerTable(signed int nPileUps, std::string trkName) {

        RootWTable* layerTable = new RootWTable();

        double precisionFlux = 2*c_coordPrecision;
        double precisionArea = 2*c_coordPrecision;
        if (trkName=="Inner") {

          precisionFlux = 1*c_coordPrecision;
          precisionArea = 4*c_coordPrecision;
        }
        if (trkName=="Outer") {

          precisionFlux = 2*c_coordPrecision;
          precisionFlux = 1*c_coordPrecision;
        }

        for (int iLayer=0; iLayer<m_nLayers; iLayer++) {

          double minFlux      = m_layerMinFluxes[iLayer]*nPileUps;
          double maxFlux      = m_layerMaxFluxes[iLayer]*nPileUps;
          double maxCellArea  = insur::trk_max_occupancy/maxFlux;

          layerTable->setContent(0, 0, "Layer no                                : ");
          layerTable->setContent(1, 0, "Radius [mm]                             : ");
          layerTable->setContent(2, 0, "Min flux in Z [particles/cm^-2]         : ");
          layerTable->setContent(3, 0, "Max flux in Z [particles/cm^-2]         : ");
          layerTable->setContent(4, 0, "Z position [mm] related to max flux     : ");
          layerTable->setContent(5, 0, "Max cell area in Z (1% occupancy) [mm^2]: ");

          layerTable->setContent(0, iLayer+1, iLayer+1);
          layerTable->setContent(1, iLayer+1, m_layerRadii[iLayer]/Units::mm   , c_coordPrecision);
          layerTable->setContent(2, iLayer+1, minFlux/(1./Units::cm2)          , precisionFlux);
          layerTable->setContent(3, iLayer+1, maxFlux/(1./Units::cm2)          , precisionFlux);
          layerTable->setContent(4, iLayer+1, m_layerMaxFluxZ[iLayer]/Units::mm, c_coordPrecision);
          layerTable->setContent(5, iLayer+1, maxCellArea/Units::mm2           , precisionArea);
        }
        return layerTable;
      }

      RootWTable* getRingTable(signed int nPileUps, std::string trkName) {

        RootWTable* ringTable = new RootWTable();

        double precisionFlux = 2*c_coordPrecision;
        double precisionArea = 2*c_coordPrecision;
        if (trkName=="Inner") {

          precisionFlux = 1*c_coordPrecision;
          precisionArea = 4*c_coordPrecision;
        }
        if (trkName=="Outer") {

          precisionFlux = 2*c_coordPrecision;
          precisionArea = 1*c_coordPrecision;
        }

        for (int iRing=0; iRing<m_nRings; iRing++) {

          double minFlux      = m_ringMinFluxes[iRing]*nPileUps;
          double maxFlux      = m_ringMaxFluxes[iRing]*nPileUps;
          double maxCellArea  = insur::trk_max_occupancy/maxFlux;

          ringTable->setContent(0, 0, "Ring no                                 : ");
          ringTable->setContent(1, 0, "Average radius [mm]                     : ");
          ringTable->setContent(2, 0, "Min flux in R [particles/cm^-2]         : ");
          ringTable->setContent(3, 0, "Max flux in R [particles/cm^-2]         : ");
          ringTable->setContent(4, 0, "Z position [mm] related to max flux     : ");
          ringTable->setContent(5, 0, "Max cell area in R (1% occupancy) [mm^2]: ");

          ringTable->setContent(0, iRing+1, iRing+1);
          ringTable->setContent(1, iRing+1, m_ringAvgRadii[iRing]/Units::mm , c_coordPrecision);
          ringTable->setContent(2, iRing+1, minFlux/(1./Units::cm2)         , precisionFlux);
          ringTable->setContent(3, iRing+1, maxFlux/(1./Units::cm2)         , precisionFlux);
          ringTable->setContent(4, iRing+1, m_ringMaxFluxZ[iRing]/Units::mm , c_coordPrecision);
          ringTable->setContent(5, iRing+1, maxCellArea/Units::mm2          , precisionArea);
        }
        return ringTable;
      }
    };

    OccupancyVisitor geometryVisitor(m_photonsMap, m_chargedMap);
    itTracker->accept(geometryVisitor);

    // Print out layer & disk table
    for (auto nPileUps : insur::trk_pile_up) {

      RootWTable*        pileUpTable = new RootWTable();
      std::ostringstream namePileUp;
      namePileUp << nPileUps;

      pileUpTable->setContent(0, 0, "Number of pile-up events: ");
      pileUpTable->setContent(0, 1, namePileUp.str());

      occupancyBarrelContent->addItem(pileUpTable);
      occupancyBarrelContent->addItem(geometryVisitor.getLayerTable(nPileUps, itTracker->myid()));
      occupancyEndcapContent->addItem(pileUpTable);
      occupancyEndcapContent->addItem(geometryVisitor.getRingTable(nPileUps, itTracker->myid()));
    }
  }

  return true;
}

