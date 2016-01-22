/*
 * AnalyzerOccupancy.cc
 *
 *  Created on: 9. 11. 2015
 *      Author: Drasal (CERN)
 */

#include <AnalyzerOccupancy.h>
#include <BFieldMap.h>
#include <rootweb.hh>
#include <Tracker.h>
#include <IrradiationMap.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <global_constants.h>
#include <Units.h>
#include <SimParms.h>
#include <TLegend.h>

using namespace insur;

AnalyzerOccupancy::AnalyzerOccupancy(std::string chargedFileName, std::string photonsFileName, std::vector<Tracker*> trackers) : AnalyzerModule(trackers)
{
  // Set geometry, i.e. individual trackers
  //for (auto it : trackers) m_trackers.push_back(it);

  // Read data from files to memory
  m_photonsMapBOnMatOn   = new IrradiationMap(photonsFileName);
  m_chargedMapBOnMatOn   = new IrradiationMap(chargedFileName);

  m_photonsMapBOnMatOnLTh= nullptr;
  m_photonsMapBOffMatOn  = nullptr;
  m_photonsMapBOnMatOff  = nullptr;
  m_photonsMapBOffMatOff = nullptr;
  m_photonsMapBOffTrkOff = nullptr;

  m_chargedMapBOnMatOnLTh= nullptr;
  m_chargedMapBOffMatOn  = nullptr;
  m_chargedMapBOnMatOff  = nullptr;
  m_chargedMapBOffMatOff = nullptr;
  m_chargedMapBOffTrkOff = nullptr;

  m_bFieldMap            = nullptr;

  m_hisChargedFluxBOnMatOn   = nullptr;
  m_hisChargedFluxBOnMatOnLTh= nullptr;
  m_hisChargedFluxBOffMatOn  = nullptr;
  m_hisChargedFluxBOnMatOff  = nullptr;
  m_hisChargedFluxBOffMatOff = nullptr;
  m_hisChargedFluxBOffTrkOff = nullptr;
  m_hisChargedRatioLTh       = nullptr;
  m_hisChargedRatioECalMat   = nullptr;
  m_hisChargedRatioMat       = nullptr;
  m_hisChargedRatioB         = nullptr;
  m_hisChargedRatioTrkB      = nullptr;
  m_hisChargedRatioMatB      = nullptr;

  m_hisPhotonsFluxBOnMatOn   = nullptr;
  m_hisPhotonsFluxBOnMatOnLTh= nullptr;
  m_hisPhotonsFluxBOffMatOn  = nullptr;
  m_hisPhotonsFluxBOnMatOff  = nullptr;
  m_hisPhotonsFluxBOffMatOff = nullptr;
  m_hisPhotonsFluxBOffTrkOff = nullptr;
  m_hisPhotonsRatioLTh       = nullptr;
  m_hisPhotonsRatioECalMat   = nullptr;
  m_hisPhotonsRatioMat       = nullptr;
  m_hisPhotonsRatioB         = nullptr;
  m_hisPhotonsRatioTrkB      = nullptr;
  m_hisPhotonsRatioMatB      = nullptr;
}

AnalyzerOccupancy::~AnalyzerOccupancy()
{
  delete m_photonsMapBOnMatOn;
  delete m_chargedMapBOnMatOn;

  if (m_photonsMapBOffMatOn  !=nullptr) delete m_photonsMapBOffMatOn;
  if (m_photonsMapBOnMatOff  !=nullptr) delete m_photonsMapBOnMatOff;
  if (m_photonsMapBOffMatOff !=nullptr) delete m_photonsMapBOffMatOff;
  if (m_chargedMapBOffMatOn  !=nullptr) delete m_chargedMapBOffMatOn;
  if (m_chargedMapBOnMatOff  !=nullptr) delete m_chargedMapBOnMatOff;
  if (m_chargedMapBOffMatOff !=nullptr) delete m_chargedMapBOffMatOff;
}

bool AnalyzerOccupancy::analyze()
{
  // Make & fill all flux histograms
  fillHistogram(m_chargedMapBOnMatOn,   m_hisChargedFluxBOnMatOn,   "ChargedFluxPerPPBOnMatOn",   "Flux of charged particles [cm^{-2}] per pp collision - B on, all material on");
  fillHistogram(m_chargedMapBOnMatOnLTh,m_hisChargedFluxBOnMatOnLTh,"ChargedFluxPerPPBOnMatOnLTh","Flux of charged particles [cm^{-2}] per pp collision - B on, all material on (e-low thr.)");
  fillHistogram(m_chargedMapBOffMatOn,  m_hisChargedFluxBOffMatOn,  "ChargedFluxPerPPBOffMatOn",  "Flux of charged particles [cm^{-2}] per pp collision - B off, all material on");
  fillHistogram(m_chargedMapBOnMatOff,  m_hisChargedFluxBOnMatOff,  "ChargedFluxPerPPBOnMatOff",  "Flux of charged particles [cm^{-2}] per pp collision - B on, material off");
  fillHistogram(m_chargedMapBOffMatOff, m_hisChargedFluxBOffMatOff, "ChargedFluxPerPPBOffMatOff", "Flux of charged particles [cm^{-2}] per pp collision - B off, material off");
  fillHistogram(m_chargedMapBOffTrkOff, m_hisChargedFluxBOffTrkOff, "ChargedFluxPerPPBOffTrkOff", "Flux of charged particles [cm^{-2}] per pp collision - B off, tracker material off");
  fillHistogram(m_photonsMapBOnMatOn,   m_hisPhotonsFluxBOnMatOn,   "PhotonsFluxPerPPBOnMatOn",   "Flux of photons [cm^{-2}] per pp collision - B on, all material on");
  fillHistogram(m_photonsMapBOnMatOnLTh,m_hisPhotonsFluxBOnMatOnLTh,"PhotonsFluxPerPPBOnMatOnLTh","Flux of photons [cm^{-2}] per pp collision - B on, all material on (e-low thr.)");
  fillHistogram(m_photonsMapBOffMatOn,  m_hisPhotonsFluxBOffMatOn,  "PhotonsFluxPerPPBOffMatOn",  "Flux of photons [cm^{-2}] per pp collision - B off, all material on");
  fillHistogram(m_photonsMapBOnMatOff,  m_hisPhotonsFluxBOnMatOff,  "PhotonsFluxPerPPBOnMatOff",  "Flux of photons [cm^{-2}] per pp collision - B on, material off");
  fillHistogram(m_photonsMapBOffMatOff, m_hisPhotonsFluxBOffMatOff, "PhotonsFluxPerPPBOffMatOff", "Flux of photons [cm^{-2}] per pp collision - B off, material off");
  fillHistogram(m_photonsMapBOffTrkOff, m_hisPhotonsFluxBOffTrkOff, "PhotonsFluxPerPPBOffTrkOff", "Flux of photons [cm^{-2}] per pp collision - B off, tracker material off");
  return true;
}

bool AnalyzerOccupancy::visualize(RootWSite& webSite, const SimParms* simParms)
{
  RootWPage* myPage = new RootWPage("Occupancy");
  myPage->setAddress("occupancy.html");
  webSite.addPage(myPage);

  // Draw magnetic fiel - map
  if (m_bFieldMap!=nullptr && m_bFieldMap->isOK()) {

    RootWContent* magFieldContent   = new RootWContent("Magnetic field map:", true);
    myPage->addContent(magFieldContent);

    TCanvas* canvasXZBField = new TCanvas("canvasXZBField", "XZ view of B field [T] (Y=0)", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    TCanvas* canvasYZBField = new TCanvas("canvasYZBField", "YZ view of B field [T] (X=0)", vis_std_canvas_sizeX, vis_min_canvas_sizeY);

    m_bFieldMap->drawXZBFieldProj(canvasXZBField, "XZ view of B field [T] (Y=0)", 0, insur::geom_max_radius, 0, insur::geom_max_length);
    m_bFieldMap->drawYZBFieldProj(canvasYZBField, "YZ view of B field [T] (X=0)", 0, insur::geom_max_radius, 0, insur::geom_max_length);

    RootWImage* anImageXZBField = new RootWImage(canvasXZBField, canvasXZBField->GetWindowWidth(), canvasXZBField->GetWindowHeight());
    anImageXZBField->setComment("XZ view of B field [T] (Y=0)");
    magFieldContent->addItem(anImageXZBField);
    RootWImage* anImageYZBField = new RootWImage(canvasYZBField, canvasYZBField->GetWindowWidth(), canvasYZBField->GetWindowHeight());
    anImageYZBField->setComment("YZ view of B field [T] (X=0)");
    magFieldContent->addItem(anImageYZBField);
  }

  // Draw plots - photons
  RootWContent* plotsPhotonsContent   = new RootWContent("Fluka simulation - photons fluxes per pp collision -> adding individual effects:", false);
  myPage->addContent(plotsPhotonsContent);

  TCanvas* canvasPhotonsBOnMatOn   = nullptr;
  TCanvas* canvasPhotonsBOffMatOn  = nullptr;
  TCanvas* canvasPhotonsBOnMatOff  = nullptr;
  TCanvas* canvasPhotonsBOffMatOff = nullptr;
  TCanvas* canvasPhotonsBOffTrkOff = nullptr;

  if (drawHistogram(canvasPhotonsBOffMatOff, m_hisPhotonsFluxBOffMatOff, m_photonsMapBOffMatOff, simParms, "PhotonsCanvasBOffMatOff", "RZ view of photons flux")) {
    canvasPhotonsBOffMatOff->SetLogz();
    m_hisPhotonsFluxBOffMatOff->SetMinimum(c_fluxMin);
    m_hisPhotonsFluxBOffMatOff->SetMaximum(c_fluxMax);
    RootWImage* anImagePhotonsBOffMatOff = new RootWImage(canvasPhotonsBOffMatOff, canvasPhotonsBOffMatOff->GetWindowWidth(), canvasPhotonsBOffMatOff->GetWindowHeight());
    anImagePhotonsBOffMatOff->setComment("RZ view of photons flux [cm^-2] in a tracker - B off,  material off");
    plotsPhotonsContent->addItem(anImagePhotonsBOffMatOff);
  }
  if (drawHistogram(canvasPhotonsBOnMatOff, m_hisPhotonsFluxBOnMatOff, m_photonsMapBOnMatOff, simParms, "PhotonsCanvasBOnMatOff", "RZ view of photons flux")) {
    canvasPhotonsBOnMatOff->SetLogz();
    m_hisPhotonsFluxBOnMatOff->SetMinimum(c_fluxMin);
    m_hisPhotonsFluxBOnMatOff->SetMaximum(c_fluxMax);
    RootWImage* anImagePhotonsBOnMatOff = new RootWImage(canvasPhotonsBOnMatOff, canvasPhotonsBOnMatOff->GetWindowWidth(), canvasPhotonsBOnMatOff->GetWindowHeight());
    anImagePhotonsBOnMatOff->setComment("RZ view of photons flux [cm^-2] in a tracker - B on,  material off");
    plotsPhotonsContent->addItem(anImagePhotonsBOnMatOff);
  }
  if (drawHistogram(canvasPhotonsBOffTrkOff, m_hisPhotonsFluxBOffTrkOff, m_photonsMapBOffTrkOff, simParms, "PhotonsCanvasBOffTrkOff", "RZ view of photons flux")) {
    canvasPhotonsBOffTrkOff->SetLogz();
    m_hisPhotonsFluxBOffTrkOff->SetMinimum(c_fluxMin);
    m_hisPhotonsFluxBOffTrkOff->SetMaximum(c_fluxMax);
    RootWImage* anImagePhotonsBOffTrkOff = new RootWImage(canvasPhotonsBOffTrkOff, canvasPhotonsBOffTrkOff->GetWindowWidth(), canvasPhotonsBOffTrkOff->GetWindowHeight());
    anImagePhotonsBOffTrkOff->setComment("RZ view of photons flux [cm^-2] in a tracker - B off,  tracker material off (calorimeter on)");
    plotsPhotonsContent->addItem(anImagePhotonsBOffTrkOff);
  }
  if (drawHistogram(canvasPhotonsBOffMatOn, m_hisPhotonsFluxBOffMatOn, m_photonsMapBOffMatOn, simParms, "PhotonsCanvasBOffMatOn", "RZ view of photonsflux")) {
    canvasPhotonsBOffMatOn->SetLogz();
    m_hisPhotonsFluxBOffMatOn->SetMinimum(c_fluxMin);
    m_hisPhotonsFluxBOffMatOn->SetMaximum(c_fluxMax);
    RootWImage* anImagePhotonsBOffMatOn = new RootWImage(canvasPhotonsBOffMatOn, canvasPhotonsBOffMatOn->GetWindowWidth(), canvasPhotonsBOffMatOn->GetWindowHeight());
    anImagePhotonsBOffMatOn->setComment("RZ view of photons flux [cm^-2] in a tracker - B off,  all material on");
    plotsPhotonsContent->addItem(anImagePhotonsBOffMatOn);
  }
  if (drawHistogram(canvasPhotonsBOnMatOn, m_hisPhotonsFluxBOnMatOn, m_photonsMapBOnMatOn, simParms, "PhotonsCanvasBOnMatOn", "RZ view of photons flux")) {
    canvasPhotonsBOnMatOn->SetLogz();
    m_hisPhotonsFluxBOnMatOn->SetMinimum(c_fluxMin);
    m_hisPhotonsFluxBOnMatOn->SetMaximum(c_fluxMax);
    RootWImage* anImagePhotonsBOnMatOn = new RootWImage(canvasPhotonsBOnMatOn, canvasPhotonsBOnMatOn->GetWindowWidth(), canvasPhotonsBOnMatOn->GetWindowHeight());
    anImagePhotonsBOnMatOn->setComment("RZ view of photons flux [cm^-2] in a tracker - B on,  all material on");
    plotsPhotonsContent->addItem(anImagePhotonsBOnMatOn);
  }


  // Draw plots - charged
  RootWContent* plotsChargedContent   = new RootWContent("Fluka simulation - charged particles fluxes per pp collision -> adding individual effects:", false);
  myPage->addContent(plotsChargedContent);

  TCanvas* canvasChargedBOnMatOn   = nullptr;
  TCanvas* canvasChargedBOffMatOn  = nullptr;
  TCanvas* canvasChargedBOnMatOff  = nullptr;
  TCanvas* canvasChargedBOffMatOff = nullptr;
  TCanvas* canvasChargedBOffTrkOff = nullptr;

  if (drawHistogram(canvasChargedBOffMatOff, m_hisChargedFluxBOffMatOff, m_chargedMapBOffMatOff, simParms, "ChargedCanvasBOffMatOff", "RZ view of charged particles flux")) {
    canvasChargedBOffMatOff->SetLogz();
    m_hisChargedFluxBOffMatOff->SetMinimum(c_fluxMin);
    m_hisChargedFluxBOffMatOff->SetMaximum(c_fluxMax);
    RootWImage* anImageChargedBOffMatOff = new RootWImage(canvasChargedBOffMatOff, canvasChargedBOffMatOff->GetWindowWidth(), canvasChargedBOffMatOff->GetWindowHeight());
    anImageChargedBOffMatOff->setComment("RZ view of charged particles flux [cm^-2] in a tracker - B off,  material off");
    plotsChargedContent->addItem(anImageChargedBOffMatOff);
  }
  if (drawHistogram(canvasChargedBOnMatOff, m_hisChargedFluxBOnMatOff, m_chargedMapBOnMatOff, simParms, "ChargedCanvasBOnMatOff", "RZ view of charged particles flux")) {
    canvasChargedBOnMatOff->SetLogz();
    m_hisChargedFluxBOnMatOff->SetMinimum(c_fluxMin);
    m_hisChargedFluxBOnMatOff->SetMaximum(c_fluxMax);
    RootWImage* anImageChargedBOnMatOff = new RootWImage(canvasChargedBOnMatOff, canvasChargedBOnMatOff->GetWindowWidth(), canvasChargedBOnMatOff->GetWindowHeight());
    anImageChargedBOnMatOff->setComment("RZ view of charged particles flux [cm^-2] in a tracker - B on,  material off");
    plotsChargedContent->addItem(anImageChargedBOnMatOff);
  }
  if (drawHistogram(canvasChargedBOffTrkOff, m_hisChargedFluxBOffTrkOff, m_chargedMapBOffTrkOff, simParms, "ChargedCanvasBOffTrkOff", "RZ view of charged particles flux")) {
    canvasChargedBOffTrkOff->SetLogz();
    m_hisChargedFluxBOffTrkOff->SetMinimum(c_fluxMin);
    m_hisChargedFluxBOffTrkOff->SetMaximum(c_fluxMax);
    RootWImage* anImageChargedBOffTrkOff = new RootWImage(canvasChargedBOffTrkOff, canvasChargedBOffTrkOff->GetWindowWidth(), canvasChargedBOffTrkOff->GetWindowHeight());
    anImageChargedBOffTrkOff->setComment("RZ view of charged particles flux [cm^-2] in a tracker - B off,  tracker material off");
    plotsChargedContent->addItem(anImageChargedBOffTrkOff);
  }
  if (drawHistogram(canvasChargedBOffMatOn, m_hisChargedFluxBOffMatOn, m_chargedMapBOffMatOn, simParms, "ChargedCanvasBOffMatOn", "RZ view of charged particles flux")) {
    canvasChargedBOffMatOn->SetLogz();
    m_hisChargedFluxBOffMatOn->SetMinimum(c_fluxMin);
    m_hisChargedFluxBOffMatOn->SetMaximum(c_fluxMax);
    RootWImage* anImageChargedBOffMatOn = new RootWImage(canvasChargedBOffMatOn, canvasChargedBOffMatOn->GetWindowWidth(), canvasChargedBOffMatOn->GetWindowHeight());
    anImageChargedBOffMatOn->setComment("RZ view of charged particles flux [cm^-2] in a tracker - B off,  all material on");
    plotsChargedContent->addItem(anImageChargedBOffMatOn);
  }
  if (drawHistogram(canvasChargedBOnMatOn, m_hisChargedFluxBOnMatOn, m_chargedMapBOnMatOn, simParms, "ChargedCanvasBOnMatOn", "RZ view of charged particles flux")) {
    canvasChargedBOnMatOn->SetLogz();
    m_hisChargedFluxBOnMatOn->SetMinimum(c_fluxMin);
    m_hisChargedFluxBOnMatOn->SetMaximum(c_fluxMax);
    RootWImage* anImageChargedBOnMatOn = new RootWImage(canvasChargedBOnMatOn, canvasChargedBOnMatOn->GetWindowWidth(), canvasChargedBOnMatOn->GetWindowHeight());
    anImageChargedBOnMatOn->setComment("RZ view of charged particles flux [cm^-2] in a tracker - B on,  all material on");
    plotsChargedContent->addItem(anImageChargedBOnMatOn);
  }

  // Ratios
  TCanvas * canvasPhotonsRatioLTh     = nullptr;
  TCanvas * canvasPhotonsRatioECalMat = nullptr;
  TCanvas * canvasPhotonsRatioMat     = nullptr;
  TCanvas * canvasPhotonsRatioB       = nullptr;
  TCanvas * canvasPhotonsRatioTrkB    = nullptr;
  TCanvas * canvasPhotonsRatioMatB    = nullptr;
  TCanvas * canvasChargedRatioLTh     = nullptr;
  TCanvas * canvasChargedRatioECalMat = nullptr;
  TCanvas * canvasChargedRatioMat     = nullptr;
  TCanvas * canvasChargedRatioB       = nullptr;
  TCanvas * canvasChargedRatioTrkB    = nullptr;
  TCanvas * canvasChargedRatioMatB    = nullptr;

  RootWContent* plotsPhotonsRatioContent  = nullptr;
  RootWContent* plotsChargedRatioContent  = nullptr;
  // Page
  if ((m_hisPhotonsFluxBOffMatOn  && m_hisChargedFluxBOffMatOn) ||
      (m_hisPhotonsFluxBOnMatOff  && m_hisChargedFluxBOnMatOff) ||
      (m_hisPhotonsFluxBOffMatOff && m_hisChargedFluxBOffMatOff)) {
    plotsPhotonsRatioContent   = new RootWContent("Fluka simulation - ratio of photons fluxes per pp collision -> adding individual effects:", true);
    plotsChargedRatioContent   = new RootWContent("Fluka simulation - ratio of charged particles fluxes per pp collision -> adding individual effects:", true);
    myPage->addContent(plotsPhotonsRatioContent);
    myPage->addContent(plotsChargedRatioContent);
  }

  // Ratio plots - photons
  if (m_hisPhotonsFluxBOnMatOff && m_hisPhotonsFluxBOffMatOff) {
    m_hisPhotonsRatioB  = (TH2D*)(m_hisPhotonsFluxBOnMatOff->Clone("PhotonsFluxPerPPRatioB"));
    m_hisPhotonsRatioB->Divide(m_hisPhotonsFluxBOffMatOff);
    m_hisPhotonsRatioB->SetTitle("Ratio of photon flux - (Mag.field)/(No material+No mag.field)");
    m_hisPhotonsRatioB->SetMaximum(5);

    if (drawHistogram(canvasPhotonsRatioB, m_hisPhotonsRatioB, m_photonsMapBOnMatOff, simParms, "PhotonsRatioB", "RZ ratio of photons fluxes")) {
      RootWImage* anImagePhotonsRatioB = new RootWImage(canvasPhotonsRatioB, canvasPhotonsRatioB->GetWindowWidth(), canvasPhotonsRatioB->GetWindowHeight());
      anImagePhotonsRatioB->setComment("RZ ratio of photons fluxes - (Mag.field)/(No material+No mag.field)");
      plotsPhotonsRatioContent->addItem(anImagePhotonsRatioB);
    }
  }
  if (m_hisPhotonsFluxBOffTrkOff && m_hisPhotonsFluxBOffMatOff) {
    m_hisPhotonsRatioECalMat  = (TH2D*)(m_hisPhotonsFluxBOffTrkOff->Clone("PhotonsFluxPerPPRatioECalMat"));
    m_hisPhotonsRatioECalMat->Divide(m_hisPhotonsFluxBOffMatOff);
    m_hisPhotonsRatioECalMat->SetTitle("Ratio of photon flux - (ECal Material)/(No material+No mag.field)");
    m_hisPhotonsRatioECalMat->SetMaximum(100);

    if (drawHistogram(canvasPhotonsRatioECalMat, m_hisPhotonsRatioECalMat, m_photonsMapBOffTrkOff, simParms, "PhotonsRatioECalMat", "RZ ratio of photons fluxes")) {
      RootWImage* anImagePhotonsRatioECalMat = new RootWImage(canvasPhotonsRatioECalMat, canvasPhotonsRatioECalMat->GetWindowWidth(), canvasPhotonsRatioECalMat->GetWindowHeight());
      anImagePhotonsRatioECalMat->setComment("RZ ratio of photons fluxes - (ECal Material)/(No material+No mag.field)");
      plotsPhotonsRatioContent->addItem(anImagePhotonsRatioECalMat);
    }
  }
  if (m_hisPhotonsFluxBOffMatOn && m_hisPhotonsFluxBOffMatOff) {
    m_hisPhotonsRatioMat  = (TH2D*)(m_hisPhotonsFluxBOffMatOn->Clone("PhotonsFluxPerPPRatioMat"));
    m_hisPhotonsRatioMat->Divide(m_hisPhotonsFluxBOffMatOff);
    m_hisPhotonsRatioMat->SetTitle("Ratio of photon flux - (All Material)/(No material+No mag.field)");
    m_hisPhotonsRatioMat->SetMaximum(100);

    if (drawHistogram(canvasPhotonsRatioMat, m_hisPhotonsRatioMat, m_photonsMapBOffMatOn, simParms, "PhotonsRatioMat", "RZ ratio of photons fluxes")) {
      RootWImage* anImagePhotonsRatioMat = new RootWImage(canvasPhotonsRatioMat, canvasPhotonsRatioMat->GetWindowWidth(), canvasPhotonsRatioMat->GetWindowHeight());
      anImagePhotonsRatioMat->setComment("RZ ratio of photons fluxes - (Material)/(No material+No mag.field)");
      plotsPhotonsRatioContent->addItem(anImagePhotonsRatioMat);
    }
  }
  if (m_hisPhotonsFluxBOnMatOn && m_hisPhotonsFluxBOffTrkOff) {
    m_hisPhotonsRatioTrkB  = (TH2D*)(m_hisPhotonsFluxBOnMatOn->Clone("PhotonsFluxPerPPRatioTrkB"));
    m_hisPhotonsRatioTrkB->Divide(m_hisPhotonsFluxBOffTrkOff);
    m_hisPhotonsRatioTrkB->SetTitle("Ratio of photon flux - (All Material+Mag.field)/(No Tracker+No mag.field)");
    m_hisPhotonsRatioTrkB->SetMaximum(5);

    if (drawHistogram(canvasPhotonsRatioTrkB, m_hisPhotonsRatioTrkB, m_photonsMapBOnMatOn, simParms, "PhotonsRatioTrkB", "RZ ratio of photons fluxes")) {
      RootWImage* anImagePhotonsRatioTrkB = new RootWImage(canvasPhotonsRatioTrkB, canvasPhotonsRatioTrkB->GetWindowWidth(), canvasPhotonsRatioTrkB->GetWindowHeight());
      anImagePhotonsRatioTrkB->setComment("RZ ratio of photons fluxes - (Material+Mag.field)/(No Tracker+No mag.field)");
      plotsPhotonsRatioContent->addItem(anImagePhotonsRatioTrkB);
    }
  }
  if (m_hisPhotonsFluxBOnMatOn && m_hisPhotonsFluxBOffMatOff) {
    m_hisPhotonsRatioMatB  = (TH2D*)(m_hisPhotonsFluxBOnMatOn->Clone("PhotonsFluxPerPPRatioMatB"));
    m_hisPhotonsRatioMatB->Divide(m_hisPhotonsFluxBOffMatOff);
    m_hisPhotonsRatioMatB->SetTitle("Ratio of photon flux - (All Material+Mag.field)/(No material+No mag.field)");
    m_hisPhotonsRatioMatB->SetMaximum(100);

    if (drawHistogram(canvasPhotonsRatioMatB, m_hisPhotonsRatioMatB, m_photonsMapBOnMatOn, simParms, "PhotonsRatioMatB", "RZ ratio of photons fluxes")) {
      RootWImage* anImagePhotonsRatioMatB = new RootWImage(canvasPhotonsRatioMatB, canvasPhotonsRatioMatB->GetWindowWidth(), canvasPhotonsRatioMatB->GetWindowHeight());
      anImagePhotonsRatioMatB->setComment("RZ ratio of photons fluxes - (Material+Mag.field)/(No material+No mag.field)");
      plotsPhotonsRatioContent->addItem(anImagePhotonsRatioMatB);
    }

    TCanvas* canvasRatioMatBProj = new TCanvas("PhotonsRatioMatBProj", "RZ ratio of photons fluxes - projection Z=?mm", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    TLegend* legend = new TLegend(0.12,0.5,0.5,0.88);
    for (int i=0; i<8; i++) {

      std::ostringstream name;
      name << i;
      TH1D* his = m_hisPhotonsRatioMatB->ProjectionY(std::string("PhotonsRatioMatBProj"+name.str()).c_str(),int(m_hisPhotonsFluxBOnMatOn->GetNbinsX()/8*i)+1,int(m_hisPhotonsFluxBOnMatOn->GetNbinsX()/8.1*i)+2);
      his->SetLineColor(i+1);
      name.str("");
      name << m_hisPhotonsFluxBOnMatOn->GetXaxis()->GetBinCenter(int(m_hisPhotonsFluxBOnMatOn->GetNbinsX()/8.1*i)+1);
      his->SetTitle(std::string("Ratio of photons flux - (All Material+Mag.field)/(No material+No mag.field)").c_str());
      his->GetYaxis()->SetRangeUser(0, 100);
      legend->AddEntry(his,std::string("Z="+name.str()+"mm").c_str());
      if (i==0) his->Draw();
      else      his->Draw("SAME");
    }
    legend->Draw("SAME");
    RootWImage* anImagePhotonsRatioMatBProj = new RootWImage(canvasRatioMatBProj, canvasRatioMatBProj->GetWindowWidth(), canvasRatioMatBProj->GetWindowHeight());
    anImagePhotonsRatioMatBProj->setComment("RZ ratio of photons fluxes - (Material+Mag.field)/(No material+No mag.field) Z=?mm");
    plotsPhotonsRatioContent->addItem(anImagePhotonsRatioMatBProj);
  }
  if (m_hisPhotonsFluxBOnMatOn && m_hisPhotonsFluxBOnMatOnLTh) {
    m_hisPhotonsRatioLTh  = (TH2D*)(m_hisPhotonsFluxBOnMatOnLTh->Clone("PhotonsFluxPerPPRatioLTh"));
    m_hisPhotonsRatioLTh->Divide(m_hisPhotonsFluxBOnMatOn);
    m_hisPhotonsRatioLTh->SetTitle("Ratio of photon flux - (Cut_{e prod}=10keV, Cut_{#gamma prod}=1keV)/(Cut_{e prod}=1MeV, Cut_{#gamma prod}=100keV)");

    if (drawHistogram(canvasPhotonsRatioLTh, m_hisPhotonsRatioLTh, m_photonsMapBOnMatOn, simParms, "PhotonsRatioLTh", "RZ ratio of photons fluxes")) {
      RootWImage* anImagePhotonsRatioLTh = new RootWImage(canvasPhotonsRatioLTh, canvasPhotonsRatioLTh->GetWindowWidth(), canvasPhotonsRatioLTh->GetWindowHeight());
      anImagePhotonsRatioLTh->setComment("RZ ratio of photons fluxes - (Cut_{e prod}=10keV, Cut_{#gamma prod}=1keV)/(Cut_{e prod}=1MeV, Cut_{#gamma prod}=100keV)");
      plotsPhotonsRatioContent->addItem(anImagePhotonsRatioLTh);
    }
  }

  // Ratio plots - charged particles
  if (m_hisChargedFluxBOnMatOff && m_hisChargedFluxBOffMatOff) {
    m_hisChargedRatioB  = (TH2D*)(m_hisChargedFluxBOnMatOff->Clone("ChargedFluxPerPPRatioB"));
    m_hisChargedRatioB->Divide(m_hisChargedFluxBOffMatOff);
    m_hisChargedRatioB->SetTitle("Ratio of charged particles flux - (Mag.field)/(No material+No mag.field)");
    m_hisChargedRatioB->SetMaximum(5);

    if (drawHistogram(canvasChargedRatioB, m_hisChargedRatioB, m_chargedMapBOnMatOff, simParms, "ChargedRatioB", "RZ ratio of charged particle fluxes")) {
      RootWImage* anImageChargedRatioB = new RootWImage(canvasChargedRatioB, canvasChargedRatioB->GetWindowWidth(), canvasChargedRatioB->GetWindowHeight());
      anImageChargedRatioB->setComment("RZ ratio of charged particles fluxes - (Mag.field)/(No material+No mag.field)");
      plotsChargedRatioContent->addItem(anImageChargedRatioB);
    }
  }
  if (m_hisChargedFluxBOffTrkOff && m_hisChargedFluxBOffMatOff) {
    m_hisChargedRatioECalMat  = (TH2D*)(m_hisChargedFluxBOffTrkOff->Clone("ChargedFluxPerPPRatioECalMat"));
    m_hisChargedRatioECalMat->Divide(m_hisChargedFluxBOffMatOff);
    m_hisChargedRatioECalMat->SetTitle("Ratio of charged particles flux - (ECal Material)/(No material+No mag.field)");
    m_hisChargedRatioECalMat->SetMaximum(5);

    if (drawHistogram(canvasChargedRatioECalMat, m_hisChargedRatioECalMat, m_chargedMapBOffTrkOff, simParms, "ChargedRatioECalMat", "RZ ratio of charged particle fluxes")) {
      RootWImage* anImageChargedRatioECalMat = new RootWImage(canvasChargedRatioECalMat, canvasChargedRatioECalMat->GetWindowWidth(), canvasChargedRatioECalMat->GetWindowHeight());
      anImageChargedRatioECalMat->setComment("RZ ratio of charged particles fluxes - (ECal Material)/(No material+No mag.field)");
      plotsChargedRatioContent->addItem(anImageChargedRatioECalMat);
    }
  }
  if (m_hisChargedFluxBOffMatOn && m_hisChargedFluxBOffMatOff) {
    m_hisChargedRatioMat  = (TH2D*)(m_hisChargedFluxBOffMatOn->Clone("ChargedFluxPerPPRatioMat"));
    m_hisChargedRatioMat->Divide(m_hisChargedFluxBOffMatOff);
    m_hisChargedRatioMat->SetTitle("Ratio of charged particles flux - (All Material)/(No material+No mag.field)");
    m_hisChargedRatioMat->SetMaximum(10);

    if (drawHistogram(canvasChargedRatioMat, m_hisChargedRatioMat, m_chargedMapBOffMatOn, simParms, "ChargedRatioMat", "RZ ratio of charged particle fluxes")) {
      RootWImage* anImageChargedRatioMat = new RootWImage(canvasChargedRatioMat, canvasChargedRatioMat->GetWindowWidth(), canvasChargedRatioMat->GetWindowHeight());
      anImageChargedRatioMat->setComment("RZ ratio of charged particles fluxes - (Material)/(No material+No mag.field)");
      plotsChargedRatioContent->addItem(anImageChargedRatioMat);
    }
  }
  if (m_hisChargedFluxBOnMatOn && m_hisChargedFluxBOffTrkOff) {
    m_hisChargedRatioTrkB  = (TH2D*)(m_hisChargedFluxBOnMatOn->Clone("ChargedFluxPerPPRatioTrkB"));
    m_hisChargedRatioTrkB->Divide(m_hisChargedFluxBOffTrkOff);
    m_hisChargedRatioTrkB->SetTitle("Ratio of charged particles flux - (All Material+Mag.field)/(No Tracker+No mag.field)");
    m_hisChargedRatioTrkB->SetMaximum(10);

    if (drawHistogram(canvasChargedRatioTrkB, m_hisChargedRatioTrkB, m_chargedMapBOnMatOn, simParms, "ChargedRatioTrkB", "RZ ratio of charged particle fluxes")) {
      RootWImage* anImageChargedRatioTrkB = new RootWImage(canvasChargedRatioTrkB, canvasChargedRatioTrkB->GetWindowWidth(), canvasChargedRatioTrkB->GetWindowHeight());
      anImageChargedRatioTrkB->setComment("RZ ratio of charged particles fluxes - (Material+Mag.field)/(No Tracker+No mag.field)");
      plotsChargedRatioContent->addItem(anImageChargedRatioTrkB);
    }
  }
  if (m_hisChargedFluxBOnMatOn && m_hisChargedFluxBOffMatOff) {
    m_hisChargedRatioMatB  = (TH2D*)(m_hisChargedFluxBOnMatOn->Clone("ChargedFluxPerPPRatioMatB"));
    m_hisChargedRatioMatB->Divide(m_hisChargedFluxBOffMatOff);
    m_hisChargedRatioMatB->SetTitle("Ratio of charged particles flux - (All Material+Mag.field)/(No material+No mag.field)");
    m_hisChargedRatioMatB->SetMaximum(10);

    if (drawHistogram(canvasChargedRatioMatB, m_hisChargedRatioMatB, m_chargedMapBOnMatOn, simParms, "ChargedRatioMatB", "RZ ratio of charged particle fluxes")) {
      RootWImage* anImageChargedRatioMatB = new RootWImage(canvasChargedRatioMatB, canvasChargedRatioMatB->GetWindowWidth(), canvasChargedRatioMatB->GetWindowHeight());
      anImageChargedRatioMatB->setComment("RZ ratio of charged particles fluxes - (Material+Mag.field)/(No material+No mag.field)");
      plotsChargedRatioContent->addItem(anImageChargedRatioMatB);
    }

    TCanvas* canvasRatioMatBProj = new TCanvas("ChargedRatioMatBProj", "RZ ratio of charged particle fluxes - projection Z=?mm", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    TLegend* legend = new TLegend(0.5,0.5,0.88,0.88);
    for (int i=0; i<8; i++) {

      std::ostringstream name;
      name << i;
      TH1D* his = m_hisChargedRatioMatB->ProjectionY(std::string("ChargedRatioMatBProj"+name.str()).c_str(),int(m_hisChargedFluxBOnMatOn->GetNbinsX()/8*i)+1,int(m_hisChargedFluxBOnMatOn->GetNbinsX()/8.1*i)+2);
      his->SetLineColor(i+1);
      name.str("");
      name << m_hisChargedFluxBOnMatOn->GetXaxis()->GetBinCenter(int(m_hisChargedFluxBOnMatOn->GetNbinsX()/8.1*i)+1);
      his->SetTitle(std::string("Ratio of charged particles flux - (All Material+Mag.field)/(No material+No mag.field)").c_str());
      his->GetYaxis()->SetRangeUser(0, 20);
      legend->AddEntry(his,std::string("Z="+name.str()+"mm").c_str());
      if (i==0) his->Draw();
      else      his->Draw("SAME");
    }
    legend->Draw("SAME");
    RootWImage* anImageChargedRatioMatBProj = new RootWImage(canvasRatioMatBProj, canvasRatioMatBProj->GetWindowWidth(), canvasRatioMatBProj->GetWindowHeight());
    anImageChargedRatioMatBProj->setComment("RZ ratio of charged particles fluxes - (Material+Mag.field)/(No material+No mag.field) Z=?mm");
    plotsChargedRatioContent->addItem(anImageChargedRatioMatBProj);
  }
  if (m_hisChargedFluxBOnMatOn && m_hisChargedFluxBOnMatOnLTh) {
    m_hisChargedRatioLTh  = (TH2D*)(m_hisChargedFluxBOnMatOnLTh->Clone("ChargedFluxPerPPRatioLTh"));
    m_hisChargedRatioLTh->Divide(m_hisChargedFluxBOnMatOn);
    m_hisChargedRatioLTh->SetTitle("Ratio of charged particles flux - (Cut_{e prod}=10keV, Cut_{#gamma prod}=1keV)/(Cut_{e prod}=1MeV, Cut_{#gamma prod}=100keV)");

    if (drawHistogram(canvasChargedRatioLTh, m_hisChargedRatioLTh, m_chargedMapBOnMatOn, simParms, "ChargedRatioLTh", "RZ ratio of charged particle fluxes")) {
      RootWImage* anImageChargedRatioLTh = new RootWImage(canvasChargedRatioLTh, canvasChargedRatioLTh->GetWindowWidth(), canvasChargedRatioLTh->GetWindowHeight());
      anImageChargedRatioLTh->setComment("RZ ratio of charged particles fluxes - (Cut_{e prod}=10keV, Cut_{#gamma prod}=1keV)/(Cut_{e prod}=1MeV, Cut_{#gamma prod}=100keV)");
      plotsChargedRatioContent->addItem(anImageChargedRatioLTh);
    }
  }

  // Plot a table with calculated occupancies
  for (auto itTracker : m_trackers) {

    // Create table
    std::string trkName = itTracker->myid();
    RootWContent* occupancyBarrelContent = new RootWContent("Occupancy - charged particles ("+trkName+"-barrel)", true);
    RootWContent* occupancyEndcapContent = new RootWContent("Occupancy - charged particles ("+trkName+"-endcap)", true);
    myPage->addContent(occupancyBarrelContent);
    myPage->addContent(occupancyEndcapContent);

    // Create visitor class & fill tables with data
    class OccupancyVisitor : public ConstGeometryVisitor {
     private:
      const bool c_assumeFlowsFromIP = true;  // Assume that all particles come from the interaction point

      IrradiationMap* m_photonsMap;
      IrradiationMap* m_chargedMap;
      int             m_nLayers;
      int             m_nDisks;
      int             m_nRings;

      std::vector<double>            m_layerRadii;             // Radius of a given layer
      std::vector<double>            m_layerMinFluxes;         // Minimum flux in a layer
      std::vector<double>            m_layerMaxFluxes;         // Maximum flux in a layer
      std::vector<double>            m_layerMaxFluxZ;          // Z-pos of the module in a layer with a maximum flux
      std::vector<std::vector<long>> m_layerNChannels;         // Number of channels to be read-out in each layer - separately for each sensor type (either 1 or 2 types)
      std::vector<std::vector<long>> m_layerNHits;             // Number of hits to be read-out in each layer - separately for each sensor type (either 1 or 2 types)
      std::vector<int>               m_layerNRods;             // Number of rods in each layer
      std::vector<int>               m_layerNModules;          // Number of modules in each layer
      std::vector<short>             m_layerNSensorsInMod;     // Number of sensors in each module -> assuming all modules are of the same type in a layer (max 2 sensors)
      std::vector<std::vector<int>>  m_layerSenAddrSparSize;   // Channel address size in bits (assuming max 2 types of sensors in the layer) -> used for sparsified data
      std::vector<std::vector<int>>  m_layerSenAddrUnsparSize; // Sensor addressing size in bits (assuming max 2 types of sensors in the layer, each channel is either 0 or 1 - n channels x 1b) -> used for unsparsified data
      std::vector<std::vector<int>>  m_layerSenNPixels;        // Number of pixels (strips) in each module (assuming max 2 types of sensors in the layer)

      std::vector<double>            m_ringAvgRadii;
      std::vector<double>            m_ringMinFluxes;
      std::vector<double>            m_ringMaxFluxes;
      std::vector<double>            m_ringMaxFluxZ;
      std::vector<std::vector<long>> m_ringNChannels;         // Number of channels to be read-out in each ring - separately for each sensor type (either 1 or 2 types)
      std::vector<std::vector<long>> m_ringNHits;             // Number of hits to be read-out in each ring - separately for each sensor type (either 1 or 2 types)
      std::vector<int>               m_ringNModules;          // Number of modules in each ring
      std::vector<short>             m_ringNSensorsInMod;     // Number of sensors in each module -> assuming all modules are of the same type in a ring (max 2 sensors)
      std::vector<std::vector<int>>  m_ringSenAddrSparSize;   // Channel address size in bits (assuming max 2 types of sensors in the ring) -> used for sparsified data
      std::vector<std::vector<int>>  m_ringSenAddrUnsparSize; // Sensor addressing size in bits (assuming max 2 types of sensors in the ring, each channel is either 0 or 1 - n channels x 1b) -> used for unsparsified data
      std::vector<std::vector<int>>  m_ringSenNPixels;        // Number of pixels (strips) in each module (assuming max 2 types of sensors in the ring)

      double    m_maxPileUp;
      double    m_maxColFreq;

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

        m_maxPileUp  = 1.;
        m_maxColFreq = 1.;
      }

      virtual ~OccupancyVisitor() {};

      void setMaxPileUp(double maxPileUp)   { m_maxPileUp = maxPileUp;}
      void setMaxColFreq(double maxColFreq) { m_maxColFreq= maxColFreq;}

      void visit(const Layer& layer) override {
        if (layer.maxZ() < 0.) return;

        double minFlux = std::numeric_limits<double>::max();
        double maxFlux = 0;
        double zPos    = 0.;
        double maxZPos = 0.;

        while (zPos<=layer.maxZ()) {

          // Assuming that all particles come from the primary interaction point
          //double cosTheta = 1;
          // Correction not needed - calculated already on surface
          //if (c_assumeFlowsFromIP && layer.placeRadius()!=0) cosTheta = cos(atan(zPos/layer.placeRadius()));

          double flux  = m_chargedMap->calculateIrradiationZR(zPos, layer.placeRadius());//*cosTheta;
                 //flux += m_photonsMap->calculateIrradiationZR(zPos, layer.placeRadius());

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
        m_layerRadii.push_back(layer.placeRadius());
        m_layerMinFluxes.push_back(minFlux);
        m_layerMaxFluxes.push_back(maxFlux);
        m_layerMaxFluxZ.push_back(maxZPos);
        m_layerNRods.push_back(layer.numRods());
        m_layerNModules.push_back(0);
        std::vector<int> vecAddrSpar, vecAddrUnspar, vecNPixels;
        m_layerSenAddrSparSize.push_back(vecAddrSpar);
        m_layerSenAddrUnsparSize.push_back(vecAddrUnspar);
        m_layerSenNPixels.push_back(vecNPixels);
        m_layerNSensorsInMod.push_back(0);

        m_nLayers++;
      }

      void visit(const BarrelModule& module) override {

        int    iLayer   = m_nLayers-1;
        double zPos     = fabs((module.planarMaxZ()+module.planarMinZ())/2.);
        double rPos     = (module.planarMaxR()+module.planarMinR())/2.;
        long   nHits    = module.area() * m_chargedMap->calculateIrradiationZR(zPos, rPos)/Units::mm2 * m_maxPileUp;

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

          iSensor++;
        }

        m_layerNModules[iLayer]++;
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

          // Assuming that all particles come from the primary interaction point
          //double cosTheta = 1;
          // Correction not needed - calculated already on surface
          //if (c_assumeFlowsFromIP && rPos!=0) cosTheta = cos(atan(rPos/ring.averageZ()));

          double flux  = m_chargedMap->calculateIrradiationZR(ring.averageZ(), rPos);//*cosTheta;
                 //flux += m_photonsMap->calculateIrradiationZR(ring.averageZ(), rPos);


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

          std::vector<long> vecRNHits, vecRNChannels;
          m_ringNHits.push_back(vecRNHits);
          m_ringNChannels.push_back(vecRNChannels);
          m_ringNModules.push_back(0);
          std::vector<int> vecRAddrSpar, vecRAddrUnspar, vecRNPixels;
          m_ringSenAddrSparSize.push_back(vecRAddrSpar);
          m_ringSenAddrUnsparSize.push_back(vecRAddrUnspar);
          m_ringSenNPixels.push_back(vecRNPixels);
          m_ringNSensorsInMod.push_back(0);
        }
        // Other disks (need to have the same number of rings)
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

      void visit(const EndcapModule& module) override {
        static int countModPls = 0;
        if ((module.planarMaxZ()+module.planarMinZ())/2.<0) return;

        int    iRing    = m_nRings-1;
        int    iDisk    = m_nDisks-1;
        double zPos     = fabs((module.planarMaxZ()+module.planarMinZ())/2.);
        double rPos     = (module.planarMaxR()+module.planarMinR())/2.;
        long   nHits    = module.area() * m_chargedMap->calculateIrradiationZR(zPos, rPos)/Units::mm2 * m_maxPileUp;

        short iSensor = 0;
        for (const auto& s : module.sensors()) {
          int nSegments         = s.numSegments();
          int nStrips           = s.numStripsAcross();
          SensorType sensorType = s.type();

          if (iSensor>=2) {
            logWARNING("Occupancy studies - module contains more than 2 sensors -> check the code, data rate will be wrong!");
            continue;
          }
          else if (iDisk==0 && m_ringNSensorsInMod[iRing]==0) m_ringNSensorsInMod[iRing] = module.numSensors();

          // Allocate memory if needed - 2 sensors per module in maximum
          if (iDisk==0 && m_ringNChannels[iRing].size()==0         && iSensor==0) m_ringNChannels[iRing].push_back(0);
          if (iDisk==0 && m_ringNChannels[iRing].size()==1         && iSensor==1) m_ringNChannels[iRing].push_back(0);
          if (iDisk==0 && m_ringNHits[iRing].size()==0             && iSensor==0) m_ringNHits[iRing].push_back(0);
          if (iDisk==0 && m_ringNHits[iRing].size()==1             && iSensor==1) m_ringNHits[iRing].push_back(0);
          if (iDisk==0 && m_ringSenAddrSparSize[iRing].size()==0   && iSensor==0) m_ringSenAddrSparSize[iRing].push_back(0);
          if (iDisk==0 && m_ringSenAddrSparSize[iRing].size()==1   && iSensor==1) m_ringSenAddrSparSize[iRing].push_back(0);
          if (iDisk==0 && m_ringSenAddrUnsparSize[iRing].size()==0 && iSensor==0) m_ringSenAddrUnsparSize[iRing].push_back(0);
          if (iDisk==0 && m_ringSenAddrUnsparSize[iRing].size()==1 && iSensor==1) m_ringSenAddrUnsparSize[iRing].push_back(0);
          if (iDisk==0 && m_ringSenNPixels[iRing].size()==0        && iSensor==0) m_ringSenNPixels[iRing].push_back(0);
          if (iDisk==0 && m_ringSenNPixels[iRing].size()==1        && iSensor==1) m_ringSenNPixels[iRing].push_back(0);

          // This separation works only for true crosssing strips -> not for several times repeating strip-lets (segments)
          //if (sensorType==SensorType::Largepix || sensorType==SensorType::Pixel) {
            long nReadOutChannels = 0;
            if (nHits<=nSegments*nStrips) nReadOutChannels = nHits;
            else                          nReadOutChannels = nSegments*nStrips;

            m_ringNChannels[iRing][iSensor] += nReadOutChannels;
            m_ringNHits[iRing][iSensor]     += nHits;
          /*}
          else if (sensorType==SensorType::Strip) {
            long nReadOutChannelsX = 0;
            long nReadOutChannelsY = 0;
            if (nHits<=nSegments) nReadOutChannelsX = nHits;
            else                  nReadOutChannelsX = nSegments;
            if (nHits<=nStrips)   nReadOutChannelsY = nHits;
            else                  nReadOutChannelsY = nStrips;
            m_ringNChannels[iRing][iSensor] += nReadOutChannelsX+nReadOutChannelsY;
            m_ringNHits[iRing][iSensor]     += nHits;
          }
          else {
            logWARNING("Occupancy studies: Sensor type not define -> couldn't calculate the data rate");
          }*/

          int senAddrSpar   = std::ceil(log2(nSegments)) + std::ceil(log2(nStrips)) * Units::b;
          int senAddrUnspar = nSegments*nStrips * Units::b;
          if (m_ringSenAddrSparSize[iRing][iSensor]!=0) {
            if (m_ringSenAddrSparSize[iRing][iSensor]!=senAddrSpar) logWARNING("Occupancy studies - module types differ within a ring -> check the code, data rate will be wrong!");
          }
          else {
            m_ringSenAddrSparSize[iRing][iSensor] = senAddrSpar;
            //std::cout << ">>Spar> " << module.numSensors() << " " << m_ringSenAddrSparSize[iLayer][iSensor] << std::endl;
          }
          if (m_ringSenAddrUnsparSize[iRing][iSensor]!=0) {
            if (m_ringSenAddrUnsparSize[iRing][iSensor]!=senAddrUnspar) logWARNING("Occupancy studies - module types differ within a ring -> check the code, data rate will be wrong!");
          }
          else {
            m_ringSenAddrUnsparSize[iRing][iSensor] = senAddrUnspar;
            //std::cout << ">>UnS> " << module.numSensors() << " " << m_ringSenAddrUnsparSize[iLayer][iSensor] << std::endl;
          }
          if (m_ringSenNPixels[iRing][iSensor]!=0) {
            if ((m_ringSenNPixels[iRing][iSensor]!=nStrips+nSegments) && (m_ringSenNPixels[iRing][iSensor]!=nStrips*nSegments)) logWARNING("Occupancy studies - module types differ within a ring -> check the code, data rate will be wrong!");
          }
          else {
            //if (sensorType==SensorType::Strip) m_ringSenNPixels[iRing][iSensor] = nStrips+nSegments;
            //else                               m_ringSenNPixels[iRing][iSensor] = nStrips*nSegments;
            m_ringSenNPixels[iRing][iSensor] = nStrips*nSegments;
          }

          iSensor++;
        }

        m_ringNModules[iRing]++;
      }

      RootWTable* getLayerTable(signed int nPileUps, std::string trkName) {

        RootWTable* layerTable = new RootWTable();

        double precisionFlux      = 2*c_coordPrecision;
        double precisionArea      = 2*c_coordPrecision;
        double precisionOccupancy = 2*c_coordPrecision;
        if (trkName=="Inner") {

          precisionFlux = 1*c_coordPrecision;
          precisionArea = 4*c_coordPrecision;
        }
        if (trkName=="Outer") {

          precisionFlux = 2*c_coordPrecision;
          precisionFlux = 1*c_coordPrecision;
        }

        double totDataRateTriggerSpar   = 0;
        double totDataRateUnTriggerSpar = 0;

        for (int iLayer=0; iLayer<m_nLayers; iLayer++) {

          double minFlux     = m_layerMinFluxes[iLayer]*nPileUps;
          double maxFlux     = m_layerMaxFluxes[iLayer]*nPileUps;
          double maxCellArea = insur::trk_max_occupancy/maxFlux;
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

          for (int iSensor=0; iSensor<numSensors; iSensor++) {
            hitRate[iSensor]           = m_layerNHits[iLayer][iSensor];
            channelRate[iSensor]       = m_layerNChannels[iLayer][iSensor];
            senAddrSparSize[iSensor]   = m_layerSenAddrSparSize[iLayer][iSensor] + addrSparClsWidth;
            senAddrUnsparSize[iSensor] = m_layerSenAddrUnsparSize[iLayer][iSensor];

            std::cout << ">iLayer>>" << iLayer << " " << channelRate[iSensor] << " " << senAddrSparSize[iSensor] << "  " << m_layerNModules[iLayer] << std::endl;
            totHitRate        += hitRate[iSensor];
            totChannelRate    += channelRate[iSensor];
            totAddrSparSize   += senAddrSparSize[iSensor];
            totAddrUnsparSize += senAddrUnsparSize[iSensor];

            double occupancy = hitRate[iSensor]/m_layerSenNPixels[iLayer][iSensor]/m_layerNModules[iLayer];
            if (occupancy>moduleOccupancy) moduleOccupancy = occupancy;

            dataRateCollisionSpar += channelRate[iSensor]*senAddrSparSize[iSensor];
            dataRateTriggerSpar   += channelRate[iSensor]*senAddrSparSize[iSensor];
            dataRateUnTriggerSpar += channelRate[iSensor]*senAddrSparSize[iSensor];
          }

          std::cout << ">>> " << dataRateCollisionSpar << std::endl;
          dataRateTriggerSpar   *= insur::trigger_freq;
          dataRateUnTriggerSpar *= insur::collision_freq;

          totDataRateTriggerSpar   += dataRateTriggerSpar;
          totDataRateUnTriggerSpar += dataRateUnTriggerSpar;

          // Layer table
          layerTable->setContent(0, 0, "Layer no                                   : ");
          layerTable->setContent(1, 0, "Radius [mm]                                : ");
          layerTable->setContent(2, 0, "Min flux in Z [particles/cm^-2]            : ");
          layerTable->setContent(3, 0, "Max flux in Z [particles/cm^-2]            : ");
          layerTable->setContent(4, 0, "Z position [mm] related to max flux        : ");
          layerTable->setContent(5, 0, "Max cell area in Z (1% occupancy) [mm^2]   : ");
          if (nPileUps==insur::trk_pile_up[insur::trk_pile_up.size()-1]) {
            layerTable->setContent(6 , 0, "#Hits per BX (bunch crossing)             : ");
            layerTable->setContent(7 , 0, "#Hit-channels per BX                      : ");
            layerTable->setContent(8 , 0, "#Hit-channels per module per BX           : ");
            layerTable->setContent(9 , 0, "Module avg occupancy (max[sen1,sen2])[%]  : ");
            layerTable->setContent(10, 0, "Module bandwidth/(addr+clsWidth=2b[b]     : ");
            layerTable->setContent(11, 0, "Mod. bandwidth(#chnls*(addr+clsWidth)[kb] : ");
            layerTable->setContent(12, 0, "Mod. bandwidth (matrix*1b/channel) [kb]   : ");
            layerTable->setContent(13, 0, "Data rate per layer - 40MHz,spars [Tb/s]  : ");
            layerTable->setContent(14, 0, "Data rate per layer -  1MHz,spars [Tb/s]  : ");
            layerTable->setContent(15, 0, "Data rate per ladder - 40Mhz,spars [Gb/s] : ");
            layerTable->setContent(16, 0, "Data rate per ladder -  1Mhz,spars [Gb/s] : ");
            layerTable->setContent(17, 0, "<b>Data rate per module - 40Mhz,spars [Gb/s]</b>: ");
            layerTable->setContent(18, 0, "<b>Data rate per module -  1Mhz,spars [Gb/s]</b>: ");
          }

          layerTable->setContent(0, iLayer+1, iLayer+1);
          layerTable->setContent(1, iLayer+1, m_layerRadii[iLayer]/Units::mm   , c_coordPrecision);
          layerTable->setContent(2, iLayer+1, minFlux/(1./Units::cm2)          , precisionFlux);
          layerTable->setContent(3, iLayer+1, maxFlux/(1./Units::cm2)          , precisionFlux);
          layerTable->setContent(4, iLayer+1, m_layerMaxFluxZ[iLayer]/Units::mm, c_coordPrecision);
          layerTable->setContent(5, iLayer+1, maxCellArea/Units::mm2           , precisionArea);
          if (nPileUps==insur::trk_pile_up[insur::trk_pile_up.size()-1]) {
            layerTable->setContent(6 , iLayer+1, totHitRate                          );
            layerTable->setContent(7 , iLayer+1, totChannelRate                      );
            layerTable->setContent(8 , iLayer+1, totChannelRate/m_layerNModules[iLayer]);
            layerTable->setContent(9 , iLayer+1, moduleOccupancy*100           , precisionOccupancy);
            layerTable->setContent(10, iLayer+1, totAddrSparSize/Units::b);
            layerTable->setContent(11, iLayer+1, dataRateCollisionSpar/m_layerNModules[iLayer]/Units::kb, 2*c_coordPrecision);
            layerTable->setContent(12, iLayer+1, totAddrUnsparSize/Units::kb                            , 2*c_coordPrecision);
            layerTable->setContent(13, iLayer+1, dataRateUnTriggerSpar/(Units::Tb/Units::s));
            layerTable->setContent(14, iLayer+1, dataRateTriggerSpar/(Units::Tb/Units::s));
            layerTable->setContent(15, iLayer+1, dataRateUnTriggerSpar/m_layerNRods[iLayer]/(Units::Gb/Units::s));
            layerTable->setContent(16, iLayer+1, dataRateTriggerSpar/m_layerNRods[iLayer]/(Units::Gb/Units::s));
            layerTable->setContent(17, iLayer+1, dataRateUnTriggerSpar/m_layerNModules[iLayer]/(Units::Gb/Units::s), 2*c_coordPrecision);
            layerTable->setContent(18, iLayer+1, dataRateTriggerSpar/m_layerNModules[iLayer]/(Units::Gb/Units::s)  , 2*c_coordPrecision);
          }
        }
        if (m_nLayers>0 && (nPileUps==insur::trk_pile_up[insur::trk_pile_up.size()-1])) {
          layerTable->setContent(0 , m_nLayers+1, "Total [TB/s]");
          layerTable->setContent(13, m_nLayers+1, totDataRateUnTriggerSpar/(Units::TB/Units::s));
          layerTable->setContent(14, m_nLayers+1, totDataRateTriggerSpar/(Units::TB/Units::s));
        }

        return layerTable;
      }

      RootWTable* getRingTable(signed int nPileUps, std::string trkName) {

        RootWTable* ringTable = new RootWTable();

        double precisionFlux      = 2*c_coordPrecision;
        double precisionArea      = 2*c_coordPrecision;
        double precisionOccupancy = 2*c_coordPrecision;

        if (trkName=="Inner") {

          precisionFlux = 1*c_coordPrecision;
          precisionArea = 4*c_coordPrecision;
        }
        if (trkName=="Outer") {

          precisionFlux = 2*c_coordPrecision;
          precisionArea = 1*c_coordPrecision;
        }

        double totDataRateTriggerSpar   = 0;
        double totDataRateUnTriggerSpar = 0;

        for (int iRing=0; iRing<m_nRings; iRing++) {

          double minFlux      = m_ringMinFluxes[iRing]*nPileUps;
          double maxFlux      = m_ringMaxFluxes[iRing]*nPileUps;
          double maxCellArea  = insur::trk_max_occupancy/maxFlux;
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

          for (int iSensor=0; iSensor<numSensors; iSensor++) {
            hitRate[iSensor]           = m_ringNHits[iRing][iSensor];
            channelRate[iSensor]       = m_ringNChannels[iRing][iSensor];
            senAddrSparSize[iSensor]   = m_ringSenAddrSparSize[iRing][iSensor] + addrSparClsWidth;
            senAddrUnsparSize[iSensor] = m_ringSenAddrUnsparSize[iRing][iSensor];

            totHitRate        += hitRate[iSensor];
            totChannelRate    += channelRate[iSensor];
            totAddrSparSize   += senAddrSparSize[iSensor];
            totAddrUnsparSize += senAddrUnsparSize[iSensor];

            double occupancy = hitRate[iSensor]/m_ringNModules[iRing]/m_ringSenNPixels[iRing][iSensor];
            if (occupancy>moduleOccupancy) moduleOccupancy = occupancy;

            dataRateCollisionSpar += channelRate[iSensor]*senAddrSparSize[iSensor];
            dataRateTriggerSpar   += channelRate[iSensor]*senAddrSparSize[iSensor];
            dataRateUnTriggerSpar += channelRate[iSensor]*senAddrSparSize[iSensor];
          }

          dataRateTriggerSpar   *= insur::trigger_freq;
          dataRateUnTriggerSpar *= insur::collision_freq;

          totDataRateTriggerSpar   += dataRateTriggerSpar;
          totDataRateUnTriggerSpar += dataRateUnTriggerSpar;

          // Ring table
          ringTable->setContent(0, 0, "Ring no                                 : ");
          ringTable->setContent(1, 0, "Average radius [mm]                     : ");
          ringTable->setContent(2, 0, "Min flux in R [particles/cm^-2]         : ");
          ringTable->setContent(3, 0, "Max flux in R [particles/cm^-2]         : ");
          ringTable->setContent(4, 0, "Z position [mm] related to max flux     : ");
          ringTable->setContent(5, 0, "Max cell area in R (1% occupancy) [mm^2]: ");
          if (nPileUps==insur::trk_pile_up[insur::trk_pile_up.size()-1]) {
            ringTable->setContent(6 , 0, "#Hits per BX (bunch crossing)             : ");
            ringTable->setContent(7 , 0, "#Hit-channels per BX                      : ");
            ringTable->setContent(8 , 0, "#Hit-channels per module per BX           : ");
            ringTable->setContent(9 , 0, "Module avg occupancy (max[sen1,sen2]) [%] : ");
            ringTable->setContent(10, 0, "Module bandwidth/(addr+clsWidth=2b[b]     : ");
            ringTable->setContent(11, 0, "Mod. bandwidth(#chnls*(addr+clsWidth)[kb] : ");
            ringTable->setContent(12, 0, "Mod. bandwidth (matrix*1b/channel) [kb]   : ");
            ringTable->setContent(13, 0, "Data rate per ringLayer-40MHz,spars [Tb/s]: ");
            ringTable->setContent(14, 0, "Data rate per ringLayer- 1MHz,spars [Tb/s]: ");
            ringTable->setContent(15, 0, "Data rate per ring - 40Mhz,spars [Gb/s]   : ");
            ringTable->setContent(16, 0, "Data rate per ring -  1Mhz,spars [Gb/s]   : ");
            ringTable->setContent(17, 0, "<b>Data rate per module - 40Mhz,spars [Gb/s]</b>: ");
            ringTable->setContent(18, 0, "<b>Data rate per module -  1Mhz,spars [Gb/s]</b>: ");
          }

          ringTable->setContent(0, iRing+1, iRing+1);
          ringTable->setContent(1, iRing+1, m_ringAvgRadii[iRing]/Units::mm , c_coordPrecision);
          ringTable->setContent(2, iRing+1, minFlux/(1./Units::cm2)         , precisionFlux);
          ringTable->setContent(3, iRing+1, maxFlux/(1./Units::cm2)         , precisionFlux);
          ringTable->setContent(4, iRing+1, m_ringMaxFluxZ[iRing]/Units::mm , c_coordPrecision);
          ringTable->setContent(5, iRing+1, maxCellArea/Units::mm2          , precisionArea);
          if (nPileUps==insur::trk_pile_up[insur::trk_pile_up.size()-1]) {
            ringTable->setContent(6 , iRing+1, totHitRate*2                        ); // Factor 2 for positive + negative side (neg. side don't used in calculations)
            ringTable->setContent(7 , iRing+1, totChannelRate*2                    ); // Factor 2 for positive + negative side (neg. side don't used in calculations)
            ringTable->setContent(8 , iRing+1, totChannelRate/m_ringNModules[iRing]);
            ringTable->setContent(9 , iRing+1, moduleOccupancy*100           , precisionOccupancy);
            ringTable->setContent(10, iRing+1, totAddrSparSize/Units::b);
            ringTable->setContent(11, iRing+1, dataRateCollisionSpar/m_ringNModules[iRing]/Units::kb, 2*c_coordPrecision);
            ringTable->setContent(12, iRing+1, totAddrUnsparSize/Units::kb                          , 2*c_coordPrecision);
            ringTable->setContent(13, iRing+1, dataRateUnTriggerSpar/(Units::Tb/Units::s)*2); // Factor 2 for positive + negative side (neg. side don't used in calculations)
            ringTable->setContent(14, iRing+1, dataRateTriggerSpar/(Units::Tb/Units::s)*2);   // Factor 2 for positive + negative side (neg. side don't used in calculations)
            ringTable->setContent(15, iRing+1, dataRateUnTriggerSpar/m_nDisks/(Units::Gb/Units::s));
            ringTable->setContent(16, iRing+1, dataRateTriggerSpar/m_nDisks/(Units::Gb/Units::s));
            ringTable->setContent(17, iRing+1, dataRateUnTriggerSpar/m_ringNModules[iRing]/(Units::Gb/Units::s), 2*c_coordPrecision);
            ringTable->setContent(18, iRing+1, dataRateTriggerSpar/m_ringNModules[iRing]/(Units::Gb/Units::s)  , 2*c_coordPrecision);
          }
          if (m_nRings>0 && (nPileUps==insur::trk_pile_up[insur::trk_pile_up.size()-1])) {
            ringTable->setContent(0 , m_nRings+1, "Total [TB/s]");
            ringTable->setContent(13, m_nRings+1, totDataRateUnTriggerSpar/(Units::TB/Units::s)*2); // Factor 2 for positive + negative side (neg. side don't used in calculations)
            ringTable->setContent(14, m_nRings+1, totDataRateTriggerSpar/(Units::TB/Units::s)*2);   // Factor 2 for positive + negative side (neg. side don't used in calculations)
          }
        }
        return ringTable;
      }
    };

    IrradiationMap* usedChargedMap = nullptr;
    IrradiationMap* usedPhotonsMap = nullptr;

    if (m_photonsMapBOnMatOnLTh!=nullptr) usedPhotonsMap = m_photonsMapBOnMatOnLTh;
    else                                  usedPhotonsMap = m_photonsMapBOnMatOn;
    if (m_photonsMapBOnMatOnLTh!=nullptr) usedChargedMap = m_chargedMapBOnMatOnLTh;
    else                                  usedChargedMap = m_chargedMapBOnMatOn;

    OccupancyVisitor geometryVisitor(usedPhotonsMap, usedChargedMap);
    geometryVisitor.setMaxPileUp(insur::trk_pile_up[insur::trk_pile_up.size()-1]);
    geometryVisitor.setMaxColFreq(insur::collision_freq);

    // set all values through visitor
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

  std::cout << "End" << std::endl;
  return true;
}

bool AnalyzerOccupancy::readMagFieldMap(std::string directory, std::string bFieldFileName)
{
  if (bFieldFileName!="") {
    m_bFieldMap = new BFieldMap(directory+"/"+bFieldFileName);
    return true;
  }
  else return false;
}

bool AnalyzerOccupancy::readNoMagFieldIrradMap(std::string directory, std::string chargedFileName, std::string photonsFileName)
{
  if (photonsFileName!="" && chargedFileName!="") {
    m_photonsMapBOffMatOn = new IrradiationMap(directory+"/"+photonsFileName);
    m_chargedMapBOffMatOn = new IrradiationMap(directory+"/"+chargedFileName);
    return true;
  }
  else return false;
}

bool AnalyzerOccupancy::readNoMaterialIrradMap(std::string directory, std::string chargedFileName, std::string photonsFileName)
{
  if (photonsFileName!="" && chargedFileName!="") {
    m_photonsMapBOnMatOff = new IrradiationMap(directory+"/"+photonsFileName);
    m_chargedMapBOnMatOff = new IrradiationMap(directory+"/"+chargedFileName);
    return true;
  }
  else return false;
}

bool AnalyzerOccupancy::readNoMagFieldNoMaterialIrradMap(std::string directory, std::string chargedFileName, std::string photonsFileName)
{
  if (photonsFileName!="" && chargedFileName!="") {
    m_photonsMapBOffMatOff = new IrradiationMap(directory+"/"+photonsFileName);
    m_chargedMapBOffMatOff = new IrradiationMap(directory+"/"+chargedFileName);
    return true;
  }
  else return false;
}

bool AnalyzerOccupancy::readNoMagFieldNoTrackerIrradMap(std::string directory, std::string chargedFileName, std::string photonsFileName)
{
  if (photonsFileName!="" && chargedFileName!="") {
    m_photonsMapBOffTrkOff = new IrradiationMap(directory+"/"+photonsFileName);
    m_chargedMapBOffTrkOff = new IrradiationMap(directory+"/"+chargedFileName);
    return true;
  }
  else return false;
}

bool AnalyzerOccupancy::readLowThresholdIrradMap(std::string directory, std::string chargedFileName, std::string photonsFileName)
{
  if (photonsFileName!="" && chargedFileName!="") {
    m_photonsMapBOnMatOnLTh = new IrradiationMap(directory+"/"+photonsFileName);
    m_chargedMapBOnMatOnLTh = new IrradiationMap(directory+"/"+chargedFileName);
    return true;
  }
  else return false;
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

bool AnalyzerOccupancy::drawHistogram(TCanvas*& canvas, TH2D* his, const IrradiationMap* map, const SimParms* simParms, std::string name, std::string title)
{
  if (his!=nullptr) {

    std::string canvasName  = name+"Canvas";
    std::string canvasTitle = title;
    canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), vis_std_canvas_sizeX, vis_min_canvas_sizeY);

    canvas->cd();
    his->Draw("COLZ");
    his->GetXaxis()->SetTitle(std::string("Z ["+map->getZUnit()+"]").c_str());
    his->GetXaxis()->SetTitleOffset(1.2);
    his->GetYaxis()->SetTitle(std::string("R ["+map->getRUnit()+"]").c_str());
    his->GetYaxis()->SetTitleOffset(1.2);
    his->GetYaxis()->SetRangeUser(simParms->bpRadius(), his->GetYaxis()->GetXmax());
    return true;
  }
  else return false;
}
