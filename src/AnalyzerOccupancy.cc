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

AnalyzerOccupancy::AnalyzerOccupancy(std::string chargedFileName, std::string photonsFileName, std::vector<Tracker*> trackers) : AnalyzerModule(trackers)
{
  // Set geometry, i.e. individual trackers
  //for (auto it : trackers) m_trackers.push_back(it);

  // Read data from files to memory
  m_photonsMap = new IrradiationMap(photonsFileName);
  m_chargedMap = new IrradiationMap(chargedFileName);

  m_photonsMapNoB      = nullptr;
  m_photonsMapNoBNoMat = nullptr;
  m_chargedMapNoB      = nullptr;
  m_chargedMapNoBNoMat = nullptr;

  m_hisChargedFlux         = nullptr;
  m_hisChargedNoBFlux      = nullptr;
  m_hisChargedNoBNoMatFlux = nullptr;
  m_hisChargedRatioMat     = nullptr;
  m_hisChargedRatioMatB    = nullptr;
  m_hisPhotonsFlux         = nullptr;
  m_hisPhotonsNoBFlux      = nullptr;
  m_hisPhotonsNoBNoMatFlux = nullptr;
  m_hisPhotonsRatioMat     = nullptr;
  m_hisPhotonsRatioMatB    = nullptr;
}

AnalyzerOccupancy::~AnalyzerOccupancy()
{
  delete m_photonsMap;
  delete m_chargedMap;

  if (m_photonsMapNoB     !=nullptr) delete m_photonsMapNoB;
  if (m_photonsMapNoBNoMat!=nullptr) delete m_photonsMapNoBNoMat;
  if (m_chargedMapNoB     !=nullptr) delete m_chargedMapNoB;
  if (m_chargedMapNoBNoMat!=nullptr) delete m_chargedMapNoBNoMat;
}

bool AnalyzerOccupancy::analyze()
{
  // Make & fill all flux histograms
  fillHistogram(m_chargedMap,         m_hisChargedFlux,         "ChargedFluxPerPP",         "Flux of charged particles [cm^{-2}] per pp collision");
  fillHistogram(m_photonsMap,         m_hisPhotonsFlux,         "PhotonsFluxPerPP",         "Flux of photons [cm^{-2}] per pp collision");
  fillHistogram(m_chargedMapNoB,      m_hisChargedNoBFlux,      "ChargedFluxPerPPNoB",      "Flux of charged particles [cm^{-2}] per pp collision - no mag. field");
  fillHistogram(m_photonsMapNoB,      m_hisPhotonsNoBFlux,      "PhotonsFluxPerPPNoB",      "Flux of photons [cm^{-2}] per pp collision - no mag. field");
  fillHistogram(m_chargedMapNoBNoMat, m_hisChargedNoBNoMatFlux, "ChargedFluxPerPPNoBNoMat", "Flux of charged particles [cm^{-2}] per pp collision - no mag. field & no material");
  fillHistogram(m_photonsMapNoBNoMat, m_hisPhotonsNoBNoMatFlux, "PhotonsFluxPerPPNoBNoMat", "Flux of photons [cm^{-2}] per pp collision - no mag. field & no material");
  return true;
}

bool AnalyzerOccupancy::visualize(RootWSite& webSite, const SimParms* simParms)
{
  RootWPage* myPage = new RootWPage("Occupancy");
  myPage->setAddress("occupancy.html");
  webSite.addPage(myPage);

  // Draw plots
  TCanvas* canvasPhotons = new TCanvas("PhotonsCanvas", "RZ View of photons flux", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
  TCanvas* canvasCharged = new TCanvas("ChargedCanvas", "RZ View of charged flux", vis_std_canvas_sizeX, vis_min_canvas_sizeY);

  TCanvas* canvasPhotonsNoB      = nullptr;
  TCanvas* canvasChargedNoB      = nullptr;
  TCanvas* canvasPhotonsNoBNoMat = nullptr;
  TCanvas* canvasChargedNoBNoMat = nullptr;

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

  if (m_hisPhotonsNoBFlux!=nullptr) {
    canvasPhotonsNoB = new TCanvas("PhotonsNoBCanvas", "RZ View of photons flux", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    canvasPhotonsNoB->cd();
    m_hisPhotonsNoBFlux->Draw("COLZ");
    m_hisPhotonsNoBFlux->GetXaxis()->SetTitle(std::string("Z ["+m_photonsMapNoB->getZUnit()+"]").c_str());
    m_hisPhotonsNoBFlux->GetXaxis()->SetTitleOffset(1.2);
    m_hisPhotonsNoBFlux->GetYaxis()->SetTitle(std::string("R ["+m_photonsMapNoB->getRUnit()+"]").c_str());
    m_hisPhotonsNoBFlux->GetYaxis()->SetTitleOffset(1.2);
    m_hisPhotonsNoBFlux->GetYaxis()->SetRangeUser(simParms->bpRadius(), m_hisPhotonsNoBFlux->GetYaxis()->GetXmax());
  }

  if (m_hisChargedNoBFlux!=nullptr) {
    canvasChargedNoB = new TCanvas("ChargedNoBCanvas", "RZ View of charged flux", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    canvasChargedNoB->cd();
    m_hisChargedNoBFlux->Draw("COLZ");
    m_hisChargedNoBFlux->GetXaxis()->SetTitle(std::string("Z ["+m_chargedMapNoB->getZUnit()+"]").c_str());
    m_hisChargedNoBFlux->GetXaxis()->SetTitleOffset(1.2);
    m_hisChargedNoBFlux->GetYaxis()->SetTitle(std::string("R ["+m_chargedMapNoB->getRUnit()+"]").c_str());
    m_hisChargedNoBFlux->GetYaxis()->SetTitleOffset(1.2);
    m_hisChargedNoBFlux->GetYaxis()->SetRangeUser(simParms->bpRadius(), m_hisChargedNoBFlux->GetYaxis()->GetXmax());
  }

  if (m_hisPhotonsNoBNoMatFlux!=nullptr) {
    canvasPhotonsNoBNoMat = new TCanvas("PhotonsNoBNoMatCanvas", "RZ View of photons flux", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    canvasPhotonsNoBNoMat->cd();
    m_hisPhotonsNoBNoMatFlux->Draw("COLZ");
    m_hisPhotonsNoBNoMatFlux->GetXaxis()->SetTitle(std::string("Z ["+m_photonsMapNoBNoMat->getZUnit()+"]").c_str());
    m_hisPhotonsNoBNoMatFlux->GetXaxis()->SetTitleOffset(1.2);
    m_hisPhotonsNoBNoMatFlux->GetYaxis()->SetTitle(std::string("R ["+m_photonsMapNoBNoMat->getRUnit()+"]").c_str());
    m_hisPhotonsNoBNoMatFlux->GetYaxis()->SetTitleOffset(1.2);
    m_hisPhotonsNoBNoMatFlux->GetYaxis()->SetRangeUser(simParms->bpRadius(), m_hisPhotonsNoBNoMatFlux->GetYaxis()->GetXmax());
  }

  if (m_hisChargedNoBNoMatFlux!=nullptr) {
    canvasChargedNoBNoMat = new TCanvas("ChargedNoBNoMatCanvas", "RZ View of charged flux", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    canvasChargedNoBNoMat->cd();
    m_hisChargedNoBNoMatFlux->Draw("COLZ");
    m_hisChargedNoBNoMatFlux->GetXaxis()->SetTitle(std::string("Z ["+m_chargedMapNoBNoMat->getZUnit()+"]").c_str());
    m_hisChargedNoBNoMatFlux->GetXaxis()->SetTitleOffset(1.2);
    m_hisChargedNoBNoMatFlux->GetYaxis()->SetTitle(std::string("R ["+m_chargedMapNoBNoMat->getRUnit()+"]").c_str());
    m_hisChargedNoBNoMatFlux->GetYaxis()->SetTitleOffset(1.2);
    m_hisChargedNoBNoMatFlux->GetYaxis()->SetRangeUser(simParms->bpRadius(), m_hisChargedNoBNoMatFlux->GetYaxis()->GetXmax());
  }

  RootWContent* plotsContent   = new RootWContent("Fluka simulation - fluxes per pp collision:", true);
  myPage->addContent(plotsContent);

  RootWImage* anImagePhotons = new RootWImage(canvasPhotons, canvasPhotons->GetWindowWidth(), canvasPhotons->GetWindowHeight());
  anImagePhotons->setComment("RZ view of photons flux [cm^-2] in a tracker");
  plotsContent->addItem(anImagePhotons);

  RootWImage* anImageCharged = new RootWImage(canvasCharged, canvasCharged->GetWindowWidth(), canvasCharged->GetWindowHeight());
  anImageCharged->setComment("RZ view of charged particles flux [cm^-2] in a tracker");
  plotsContent->addItem(anImageCharged);

  if (canvasPhotonsNoB!=nullptr) {
    RootWImage* anImagePhotonsNoB = new RootWImage(canvasPhotonsNoB, canvasPhotonsNoB->GetWindowWidth(), canvasPhotonsNoB->GetWindowHeight());
    anImagePhotonsNoB->setComment("RZ view of photons flux [cm^-2] in a tracker (B=0T)");
    plotsContent->addItem(anImagePhotonsNoB);
  }
  if (canvasChargedNoB!=nullptr) {
    RootWImage* anImageChargedNoB = new RootWImage(canvasChargedNoB, canvasChargedNoB->GetWindowWidth(), canvasChargedNoB->GetWindowHeight());
    anImageChargedNoB->setComment("RZ view of charged particles flux [cm^-2] in a tracker (B=0T)");
    plotsContent->addItem(anImageChargedNoB);
  }
  if (canvasPhotonsNoBNoMat!=nullptr) {
    RootWImage* anImagePhotonsNoBNoMat = new RootWImage(canvasPhotonsNoBNoMat, canvasPhotonsNoBNoMat->GetWindowWidth(), canvasPhotonsNoBNoMat->GetWindowHeight());
    anImagePhotonsNoBNoMat->setComment("RZ view of photons flux [cm^-2] from pp collision (B=0T))");
    plotsContent->addItem(anImagePhotonsNoBNoMat);
  }
  if (canvasChargedNoBNoMat!=nullptr) {
    RootWImage* anImageChargedNoBNoMat = new RootWImage(canvasChargedNoBNoMat, canvasChargedNoBNoMat->GetWindowWidth(), canvasChargedNoBNoMat->GetWindowHeight());
    anImageChargedNoBNoMat->setComment("RZ view of charged particles flux [cm^-2] from pp collision (B=0T))");
    plotsContent->addItem(anImageChargedNoBNoMat);
  }

  // Ratios
  TCanvas * canvasPhotonsRatioMat  = nullptr;
  TCanvas * canvasPhotonsRatioMatB = nullptr;
  TCanvas * canvasChargedRatioMat  = nullptr;
  TCanvas * canvasChargedRatioMatB = nullptr;
  RootWContent* plotsRatioContent  = nullptr;

  if ((m_hisPhotonsNoBFlux && m_hisPhotonsNoBNoMatFlux) || (m_hisChargedNoBFlux && m_hisChargedNoBNoMatFlux)) {
    plotsRatioContent   = new RootWContent("Fluka simulation - ratio of fluxes per pp collision:", true);
    myPage->addContent(plotsRatioContent);
  }

  if (m_hisPhotonsNoBFlux && m_hisPhotonsNoBNoMatFlux) {

    m_hisPhotonsRatioMat  = (TH2D*)(m_hisPhotonsNoBFlux->Clone("PhotonsFluxPerPPRatioMat"));
    m_hisPhotonsRatioMat->Divide(m_hisPhotonsNoBNoMatFlux);
    m_hisPhotonsRatioMat->SetTitle("Ratio of photon flux - Material/(No material+No mag.field)");

    m_hisPhotonsRatioMatB = (TH2D*)(m_hisPhotonsFlux->Clone("PhotonsFluxPerPPRatioMatB"));
    m_hisPhotonsRatioMatB->Divide(m_hisPhotonsNoBNoMatFlux);
    m_hisPhotonsRatioMatB->SetTitle("Ratio of photon flux - (Material+Mag.field)/(No material+No mag.field)");

    canvasPhotonsRatioMat = new TCanvas("PhotonsRatioMatCanvas", "RZ Ratio of photons fluxes", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    canvasPhotonsRatioMat->cd();
    m_hisPhotonsRatioMat->Draw("COLZ");
    m_hisPhotonsRatioMat->GetXaxis()->SetTitle(std::string("Z ["+m_photonsMapNoB->getZUnit()+"]").c_str());
    m_hisPhotonsRatioMat->GetXaxis()->SetTitleOffset(1.2);
    m_hisPhotonsRatioMat->GetYaxis()->SetTitle(std::string("R ["+m_photonsMapNoB->getRUnit()+"]").c_str());
    m_hisPhotonsRatioMat->GetYaxis()->SetTitleOffset(1.2);
    m_hisPhotonsRatioMat->GetYaxis()->SetRangeUser(simParms->bpRadius(), m_hisPhotonsNoBFlux->GetYaxis()->GetXmax());

    canvasPhotonsRatioMatB = new TCanvas("PhotonsRatioMatBCanvas", "RZ Ratio of photons fluxes", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    canvasPhotonsRatioMatB->cd();
    m_hisPhotonsRatioMatB->Draw("COLZ");
    m_hisPhotonsRatioMatB->GetXaxis()->SetTitle(std::string("Z ["+m_photonsMapNoBNoMat->getZUnit()+"]").c_str());
    m_hisPhotonsRatioMatB->GetXaxis()->SetTitleOffset(1.2);
    m_hisPhotonsRatioMatB->GetYaxis()->SetTitle(std::string("R ["+m_photonsMapNoBNoMat->getRUnit()+"]").c_str());
    m_hisPhotonsRatioMatB->GetYaxis()->SetTitleOffset(1.2);
    m_hisPhotonsRatioMatB->GetYaxis()->SetRangeUser(simParms->bpRadius(), m_hisPhotonsNoBNoMatFlux->GetYaxis()->GetXmax());

    RootWImage* anImagePhotonsRatioMat = new RootWImage(canvasPhotonsRatioMat, canvasPhotonsRatioMat->GetWindowWidth(), canvasPhotonsRatioMat->GetWindowHeight());
    anImagePhotonsRatioMat->setComment("RZ ratio of photons fluxes - Material/(No material+No mag.field)");
    plotsRatioContent->addItem(anImagePhotonsRatioMat);

    RootWImage* anImagePhotonsRatioMatB = new RootWImage(canvasPhotonsRatioMatB, canvasPhotonsRatioMatB->GetWindowWidth(), canvasPhotonsRatioMatB->GetWindowHeight());
    anImagePhotonsRatioMatB->setComment("RZ ratio of photons fluxes - (Material+Mag.field)/(No material+No mag.field)");
    plotsRatioContent->addItem(anImagePhotonsRatioMatB);
  }
  if (m_hisChargedNoBFlux && m_hisChargedNoBNoMatFlux) {
    m_hisChargedRatioMat  = (TH2D*)(m_hisChargedNoBFlux->Clone("ChargedFluxPerPPRatioMat"));
    m_hisChargedRatioMat->Divide(m_hisChargedNoBNoMatFlux);
    m_hisChargedRatioMat->SetTitle("Ratio of charged particles flux - Material/No material (B=0T)");

    m_hisChargedRatioMatB = (TH2D*)(m_hisChargedFlux->Clone("ChargedFluxPerPPRatioMatB"));
    m_hisChargedRatioMatB->Divide(m_hisChargedNoBNoMatFlux);
    m_hisChargedRatioMatB->SetTitle("Ratio of charged particles flux - (Material+Mag.field)/(No material+No mag.field)");

    canvasChargedRatioMat = new TCanvas("ChargedRatioMatCanvas", "RZ Ratio of charged particles fluxes", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    canvasChargedRatioMat->cd();
    m_hisChargedRatioMat->Draw("COLZ");
    m_hisChargedRatioMat->GetXaxis()->SetTitle(std::string("Z ["+m_chargedMapNoB->getZUnit()+"]").c_str());
    m_hisChargedRatioMat->GetXaxis()->SetTitleOffset(1.2);
    m_hisChargedRatioMat->GetYaxis()->SetTitle(std::string("R ["+m_chargedMapNoB->getRUnit()+"]").c_str());
    m_hisChargedRatioMat->GetYaxis()->SetTitleOffset(1.2);
    m_hisChargedRatioMat->GetYaxis()->SetRangeUser(simParms->bpRadius(), m_hisChargedNoBFlux->GetYaxis()->GetXmax());

    canvasChargedRatioMatB = new TCanvas("ChargedRatioMatBCanvas", "RZ Ratio of charged particles fluxes", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    canvasChargedRatioMatB->cd();
    m_hisChargedRatioMatB->Draw("COLZ");
    m_hisChargedRatioMatB->GetXaxis()->SetTitle(std::string("Z ["+m_chargedMapNoBNoMat->getZUnit()+"]").c_str());
    m_hisChargedRatioMatB->GetXaxis()->SetTitleOffset(1.2);
    m_hisChargedRatioMatB->GetYaxis()->SetTitle(std::string("R ["+m_chargedMapNoBNoMat->getRUnit()+"]").c_str());
    m_hisChargedRatioMatB->GetYaxis()->SetTitleOffset(1.2);
    m_hisChargedRatioMatB->GetYaxis()->SetRangeUser(simParms->bpRadius(), m_hisChargedNoBNoMatFlux->GetYaxis()->GetXmax());

    RootWImage* anImageChargedRatioMat = new RootWImage(canvasChargedRatioMat, canvasChargedRatioMat->GetWindowWidth(), canvasChargedRatioMat->GetWindowHeight());
    anImageChargedRatioMat->setComment("RZ ratio of charged particles fluxes - Material/No material (B=0T)");
    plotsRatioContent->addItem(anImageChargedRatioMat);

    RootWImage* anImageChargedRatioMatB = new RootWImage(canvasChargedRatioMatB, canvasChargedRatioMatB->GetWindowWidth(), canvasChargedRatioMatB->GetWindowHeight());
    anImageChargedRatioMatB->setComment("RZ ratio of charged particles fluxes - (Material+Mag.field)/(No material+No mag.field)");
    plotsRatioContent->addItem(anImageChargedRatioMatB);
  }

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

void AnalyzerOccupancy::readNoMagFieldMap(std::string chargedFileName, std::string photonsFileName)
{
  m_photonsMapNoB = new IrradiationMap(photonsFileName);
  m_chargedMapNoB = new IrradiationMap(chargedFileName);
}

void AnalyzerOccupancy::readNoMagFieldNoMaterialMap(std::string chargedFileName, std::string photonsFileName)
{
  m_photonsMapNoBNoMat = new IrradiationMap(photonsFileName);
  m_chargedMapNoBNoMat = new IrradiationMap(chargedFileName);
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
      }
    }
    return true;
  }
  else return false;
}
