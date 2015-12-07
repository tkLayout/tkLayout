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
  m_hisChargedRatioMat       = nullptr;
  m_hisChargedRatioB         = nullptr;
  m_hisChargedRatioMatB      = nullptr;

  m_hisPhotonsFluxBOnMatOn   = nullptr;
  m_hisPhotonsFluxBOnMatOnLTh= nullptr;
  m_hisPhotonsFluxBOffMatOn  = nullptr;
  m_hisPhotonsFluxBOnMatOff  = nullptr;
  m_hisPhotonsFluxBOffMatOff = nullptr;
  m_hisPhotonsFluxBOffTrkOff = nullptr;
  m_hisPhotonsRatioLTh       = nullptr;
  m_hisPhotonsRatioMat       = nullptr;
  m_hisPhotonsRatioB         = nullptr;
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
  fillHistogram(m_chargedMapBOnMatOn,   m_hisChargedFluxBOnMatOn,   "ChargedFluxPerPPBOnMatOn",   "Flux of charged particles [cm^{-2}] per pp collision - B on, material on");
  fillHistogram(m_chargedMapBOnMatOnLTh,m_hisChargedFluxBOnMatOnLTh,"ChargedFluxPerPPBOnMatOnLTh","Flux of charged particles [cm^{-2}] per pp collision - B on, material on (e-low thr.)");
  fillHistogram(m_chargedMapBOffMatOn,  m_hisChargedFluxBOffMatOn,  "ChargedFluxPerPPBOffMatOn",  "Flux of charged particles [cm^{-2}] per pp collision - B off, material on");
  fillHistogram(m_chargedMapBOnMatOff,  m_hisChargedFluxBOnMatOff,  "ChargedFluxPerPPBOnMatOff",  "Flux of charged particles [cm^{-2}] per pp collision - B on, material off");
  fillHistogram(m_chargedMapBOffMatOff, m_hisChargedFluxBOffMatOff, "ChargedFluxPerPPBOffMatOff", "Flux of charged particles [cm^{-2}] per pp collision - B off, material off");
  fillHistogram(m_chargedMapBOffTrkOff, m_hisChargedFluxBOffTrkOff, "ChargedFluxPerPPBOffTrkOff", "Flux of charged particles [cm^{-2}] per pp collision - B off, tracker off");
  fillHistogram(m_photonsMapBOnMatOn,   m_hisPhotonsFluxBOnMatOn,   "PhotonsFluxPerPPBOnMatOn",   "Flux of photons [cm^{-2}] per pp collision - B on, material on");
  fillHistogram(m_photonsMapBOnMatOnLTh,m_hisPhotonsFluxBOnMatOnLTh,"PhotonsFluxPerPPBOnMatOnLTh","Flux of photons [cm^{-2}] per pp collision - B on, material on (e-low thr.)");
  fillHistogram(m_photonsMapBOffMatOn,  m_hisPhotonsFluxBOffMatOn,  "PhotonsFluxPerPPBOffMatOn",  "Flux of photons [cm^{-2}] per pp collision - B off, material on");
  fillHistogram(m_photonsMapBOnMatOff,  m_hisPhotonsFluxBOnMatOff,  "PhotonsFluxPerPPBOnMatOff",  "Flux of photons [cm^{-2}] per pp collision - B on, material off");
  fillHistogram(m_photonsMapBOffMatOff, m_hisPhotonsFluxBOffMatOff, "PhotonsFluxPerPPBOffMatOff", "Flux of photons [cm^{-2}] per pp collision - B off, material off");
  fillHistogram(m_photonsMapBOffTrkOff, m_hisPhotonsFluxBOffTrkOff, "PhotonsFluxPerPPBOffTrkOff", "Flux of photons [cm^{-2}] per pp collision - B off, tracker off");
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
  RootWContent* plotsPhotonsContent   = new RootWContent("Fluka simulation - photons fluxes per pp collision:", false);
  myPage->addContent(plotsPhotonsContent);

  TCanvas* canvasPhotonsBOnMatOn   = nullptr;
  TCanvas* canvasPhotonsBOffMatOn  = nullptr;
  TCanvas* canvasPhotonsBOnMatOff  = nullptr;
  TCanvas* canvasPhotonsBOffMatOff = nullptr;
  TCanvas* canvasPhotonsBOffTrkOff = nullptr;

  if (drawHistogram(canvasPhotonsBOnMatOn, m_hisPhotonsFluxBOnMatOn, m_photonsMapBOnMatOn, simParms, "PhotonsCanvasBOnMatOn", "RZ view of photons flux")) {
    canvasPhotonsBOnMatOn->SetLogz();
    RootWImage* anImagePhotonsBOnMatOn = new RootWImage(canvasPhotonsBOnMatOn, canvasPhotonsBOnMatOn->GetWindowWidth(), canvasPhotonsBOnMatOn->GetWindowHeight());
    anImagePhotonsBOnMatOn->setComment("RZ view of photons flux [cm^-2] in a tracker - B on,  material on");
    plotsPhotonsContent->addItem(anImagePhotonsBOnMatOn);
  }
  if (drawHistogram(canvasPhotonsBOffMatOn, m_hisPhotonsFluxBOffMatOn, m_photonsMapBOffMatOn, simParms, "PhotonsCanvasBOffMatOn", "RZ view of photonsflux")) {
    canvasPhotonsBOffMatOn->SetLogz();
    RootWImage* anImagePhotonsBOffMatOn = new RootWImage(canvasPhotonsBOffMatOn, canvasPhotonsBOffMatOn->GetWindowWidth(), canvasPhotonsBOffMatOn->GetWindowHeight());
    anImagePhotonsBOffMatOn->setComment("RZ view of photons flux [cm^-2] in a tracker - B off,  material on");
    plotsPhotonsContent->addItem(anImagePhotonsBOffMatOn);
  }
  if (drawHistogram(canvasPhotonsBOnMatOff, m_hisPhotonsFluxBOnMatOff, m_photonsMapBOnMatOff, simParms, "PhotonsCanvasBOnMatOff", "RZ view of photons flux")) {
    canvasPhotonsBOnMatOff->SetLogz();
    RootWImage* anImagePhotonsBOnMatOff = new RootWImage(canvasPhotonsBOnMatOff, canvasPhotonsBOnMatOff->GetWindowWidth(), canvasPhotonsBOnMatOff->GetWindowHeight());
    anImagePhotonsBOnMatOff->setComment("RZ view of photons flux [cm^-2] in a tracker - B on,  material off");
    plotsPhotonsContent->addItem(anImagePhotonsBOnMatOff);
  }
  if (drawHistogram(canvasPhotonsBOffMatOff, m_hisPhotonsFluxBOffMatOff, m_photonsMapBOffMatOff, simParms, "PhotonsCanvasBOffMatOff", "RZ view of photons flux")) {
    canvasPhotonsBOffMatOff->SetLogz();
    RootWImage* anImagePhotonsBOffMatOff = new RootWImage(canvasPhotonsBOffMatOff, canvasPhotonsBOffMatOff->GetWindowWidth(), canvasPhotonsBOffMatOff->GetWindowHeight());
    anImagePhotonsBOffMatOff->setComment("RZ view of photons flux [cm^-2] in a tracker - B off,  material off");
    plotsPhotonsContent->addItem(anImagePhotonsBOffMatOff);
  }
  if (drawHistogram(canvasPhotonsBOffTrkOff, m_hisPhotonsFluxBOffTrkOff, m_photonsMapBOffTrkOff, simParms, "PhotonsCanvasBOffTrkOff", "RZ view of photons flux")) {
    canvasPhotonsBOffTrkOff->SetLogz();
    RootWImage* anImagePhotonsBOffTrkOff = new RootWImage(canvasPhotonsBOffTrkOff, canvasPhotonsBOffTrkOff->GetWindowWidth(), canvasPhotonsBOffTrkOff->GetWindowHeight());
    anImagePhotonsBOffTrkOff->setComment("RZ view of photons flux [cm^-2] in a tracker - B off,  tracker material off");
    plotsPhotonsContent->addItem(anImagePhotonsBOffTrkOff);
  }

  // Draw plots - charged
  RootWContent* plotsChargedContent   = new RootWContent("Fluka simulation - charged particles fluxes per pp collision:", false);
  myPage->addContent(plotsChargedContent);

  TCanvas* canvasChargedBOnMatOn   = nullptr;
  TCanvas* canvasChargedBOffMatOn  = nullptr;
  TCanvas* canvasChargedBOnMatOff  = nullptr;
  TCanvas* canvasChargedBOffMatOff = nullptr;
  TCanvas* canvasChargedBOffTrkOff = nullptr;

  if (drawHistogram(canvasChargedBOnMatOn, m_hisChargedFluxBOnMatOn, m_chargedMapBOnMatOn, simParms, "ChargedCanvasBOnMatOn", "RZ view of charged particles flux")) {
    canvasChargedBOnMatOn->SetLogz();
    RootWImage* anImageChargedBOnMatOn = new RootWImage(canvasChargedBOnMatOn, canvasChargedBOnMatOn->GetWindowWidth(), canvasChargedBOnMatOn->GetWindowHeight());
    anImageChargedBOnMatOn->setComment("RZ view of charged particles flux [cm^-2] in a tracker - B on,  material on");
    plotsChargedContent->addItem(anImageChargedBOnMatOn);
  }
  if (drawHistogram(canvasChargedBOffMatOn, m_hisChargedFluxBOffMatOn, m_chargedMapBOffMatOn, simParms, "ChargedCanvasBOffMatOn", "RZ view of charged particles flux")) {
    canvasChargedBOffMatOn->SetLogz();
    RootWImage* anImageChargedBOffMatOn = new RootWImage(canvasChargedBOffMatOn, canvasChargedBOffMatOn->GetWindowWidth(), canvasChargedBOffMatOn->GetWindowHeight());
    anImageChargedBOffMatOn->setComment("RZ view of charged particles flux [cm^-2] in a tracker - B off,  material on");
    plotsChargedContent->addItem(anImageChargedBOffMatOn);
  }
  if (drawHistogram(canvasChargedBOnMatOff, m_hisChargedFluxBOnMatOff, m_chargedMapBOnMatOff, simParms, "ChargedCanvasBOnMatOff", "RZ view of charged particles flux")) {
    canvasChargedBOnMatOff->SetLogz();
    RootWImage* anImageChargedBOnMatOff = new RootWImage(canvasChargedBOnMatOff, canvasChargedBOnMatOff->GetWindowWidth(), canvasChargedBOnMatOff->GetWindowHeight());
    anImageChargedBOnMatOff->setComment("RZ view of charged particles flux [cm^-2] in a tracker - B on,  material off");
    plotsChargedContent->addItem(anImageChargedBOnMatOff);
  }
  if (drawHistogram(canvasChargedBOffMatOff, m_hisChargedFluxBOffMatOff, m_chargedMapBOffMatOff, simParms, "ChargedCanvasBOffMatOff", "RZ view of charged particles flux")) {
    canvasChargedBOffMatOff->SetLogz();
    RootWImage* anImageChargedBOffMatOff = new RootWImage(canvasChargedBOffMatOff, canvasChargedBOffMatOff->GetWindowWidth(), canvasChargedBOffMatOff->GetWindowHeight());
    anImageChargedBOffMatOff->setComment("RZ view of charged particles flux [cm^-2] in a tracker - B off,  material off");
    plotsChargedContent->addItem(anImageChargedBOffMatOff);
  }
  if (drawHistogram(canvasChargedBOffTrkOff, m_hisChargedFluxBOffTrkOff, m_chargedMapBOffTrkOff, simParms, "ChargedCanvasBOffTrkOff", "RZ view of charged particles flux")) {
    canvasChargedBOffTrkOff->SetLogz();
    RootWImage* anImageChargedBOffTrkOff = new RootWImage(canvasChargedBOffTrkOff, canvasChargedBOffTrkOff->GetWindowWidth(), canvasChargedBOffTrkOff->GetWindowHeight());
    anImageChargedBOffTrkOff->setComment("RZ view of charged particles flux [cm^-2] in a tracker - B off,  tracker material off");
    plotsChargedContent->addItem(anImageChargedBOffTrkOff);
  }

  // Ratios
  TCanvas * canvasPhotonsRatioLTh  = nullptr;
  TCanvas * canvasPhotonsRatioMat  = nullptr;
  TCanvas * canvasPhotonsRatioB    = nullptr;
  TCanvas * canvasPhotonsRatioMatB = nullptr;
  TCanvas * canvasChargedRatioLTh  = nullptr;
  TCanvas * canvasChargedRatioMat  = nullptr;
  TCanvas * canvasChargedRatioB    = nullptr;
  TCanvas * canvasChargedRatioMatB = nullptr;

  RootWContent* plotsPhotonsRatioContent  = nullptr;
  RootWContent* plotsChargedRatioContent  = nullptr;
  // Page
  if ((m_hisPhotonsFluxBOffMatOn  && m_hisChargedFluxBOffMatOn) ||
      (m_hisPhotonsFluxBOnMatOff  && m_hisChargedFluxBOnMatOff) ||
      (m_hisPhotonsFluxBOffMatOff && m_hisChargedFluxBOffMatOff)) {
    plotsPhotonsRatioContent   = new RootWContent("Fluka simulation - ratio of photons fluxes per pp collision:", true);
    plotsChargedRatioContent   = new RootWContent("Fluka simulation - ratio of charged particles fluxes per pp collision:", true);
    myPage->addContent(plotsPhotonsRatioContent);
    myPage->addContent(plotsChargedRatioContent);
  }

  // Ratio plots - photons
  if (m_hisPhotonsFluxBOnMatOn && m_hisPhotonsFluxBOnMatOnLTh) {
    m_hisPhotonsRatioLTh  = (TH2D*)(m_hisPhotonsFluxBOnMatOnLTh->Clone("PhotonsFluxPerPPRatioLTh"));
    m_hisPhotonsRatioLTh->Divide(m_hisPhotonsFluxBOnMatOn);
    m_hisPhotonsRatioLTh->SetTitle("Ratio of photon flux - (Cut_{e prod}=100keV, Cut_{#gamma prod}=10keV)/(Cut_{e prod}=1MeV, Cut_{#gamma prod}=100keV)");

    if (drawHistogram(canvasPhotonsRatioLTh, m_hisPhotonsRatioLTh, m_photonsMapBOnMatOn, simParms, "PhotonsRatioLTh", "RZ ratio of photons fluxes")) {
      RootWImage* anImagePhotonsRatioLTh = new RootWImage(canvasPhotonsRatioLTh, canvasPhotonsRatioLTh->GetWindowWidth(), canvasPhotonsRatioLTh->GetWindowHeight());
      anImagePhotonsRatioLTh->setComment("RZ ratio of photons fluxes - (Cut_{e prod}=100keV, Cut_{#gamma prod}=10keV)/(Cut_{e prod}=1MeV, Cut_{#gamma prod}=100keV)");
      plotsPhotonsRatioContent->addItem(anImagePhotonsRatioLTh);
    }
  }
  if (m_hisPhotonsFluxBOffMatOn && m_hisPhotonsFluxBOffMatOff) {
    m_hisPhotonsRatioMat  = (TH2D*)(m_hisPhotonsFluxBOffMatOn->Clone("PhotonsFluxPerPPRatioMat"));
    m_hisPhotonsRatioMat->Divide(m_hisPhotonsFluxBOffMatOff);
    m_hisPhotonsRatioMat->SetTitle("Ratio of photon flux - (Material)/(No material+No mag.field)");

    if (drawHistogram(canvasPhotonsRatioMat, m_hisPhotonsRatioMat, m_photonsMapBOffMatOn, simParms, "PhotonsRatioMat", "RZ ratio of photons fluxes")) {
      RootWImage* anImagePhotonsRatioMat = new RootWImage(canvasPhotonsRatioMat, canvasPhotonsRatioMat->GetWindowWidth(), canvasPhotonsRatioMat->GetWindowHeight());
      anImagePhotonsRatioMat->setComment("RZ ratio of photons fluxes - (Material)/(No material+No mag.field)");
      plotsPhotonsRatioContent->addItem(anImagePhotonsRatioMat);
    }
  }
  if (m_hisPhotonsFluxBOnMatOff && m_hisPhotonsFluxBOnMatOff) {
    m_hisPhotonsRatioB  = (TH2D*)(m_hisPhotonsFluxBOnMatOff->Clone("PhotonsFluxPerPPRatioB"));
    m_hisPhotonsRatioB->Divide(m_hisPhotonsFluxBOffMatOff);
    m_hisPhotonsRatioB->SetTitle("Ratio of photon flux - (Mag.field)/(No material+No mag.field)");

    if (drawHistogram(canvasPhotonsRatioB, m_hisPhotonsRatioB, m_photonsMapBOnMatOff, simParms, "PhotonsRatioB", "RZ ratio of photons fluxes")) {
      RootWImage* anImagePhotonsRatioB = new RootWImage(canvasPhotonsRatioB, canvasPhotonsRatioB->GetWindowWidth(), canvasPhotonsRatioB->GetWindowHeight());
      anImagePhotonsRatioB->setComment("RZ ratio of photons fluxes - (Mag.field)/(No material+No mag.field)");
      plotsPhotonsRatioContent->addItem(anImagePhotonsRatioB);
    }
  }
  if (m_hisPhotonsFluxBOnMatOn && m_hisPhotonsFluxBOffMatOff) {
    m_hisPhotonsRatioMatB  = (TH2D*)(m_hisPhotonsFluxBOnMatOn->Clone("PhotonsFluxPerPPRatioMatB"));
    m_hisPhotonsRatioMatB->Divide(m_hisPhotonsFluxBOffMatOff);
    m_hisPhotonsRatioMatB->SetTitle("Ratio of photon flux - (Material+Mag.field)/(No material+No mag.field)");

    if (drawHistogram(canvasPhotonsRatioMatB, m_hisPhotonsRatioMatB, m_photonsMapBOnMatOn, simParms, "PhotonsRatioMatB", "RZ ratio of photons fluxes")) {
      RootWImage* anImagePhotonsRatioMatB = new RootWImage(canvasPhotonsRatioMatB, canvasPhotonsRatioMatB->GetWindowWidth(), canvasPhotonsRatioMatB->GetWindowHeight());
      anImagePhotonsRatioMatB->setComment("RZ ratio of photons fluxes - (Material+Mag.field)/(No material+No mag.field)");
      plotsPhotonsRatioContent->addItem(anImagePhotonsRatioMatB);
    }
  }

  // Ratio plots - charged particles
  if (m_hisChargedFluxBOnMatOn && m_hisChargedFluxBOnMatOnLTh) {
    m_hisChargedRatioLTh  = (TH2D*)(m_hisChargedFluxBOnMatOnLTh->Clone("ChargedFluxPerPPRatioLTh"));
    m_hisChargedRatioLTh->Divide(m_hisChargedFluxBOnMatOn);
    m_hisChargedRatioLTh->SetTitle("Ratio of charged particles flux - (Cut_{e prod}=100keV, Cut_{#gamma prod}=10keV)/(Cut_{e prod}=1MeV, Cut_{#gamma prod}=100keV)");

    if (drawHistogram(canvasChargedRatioLTh, m_hisChargedRatioLTh, m_chargedMapBOnMatOn, simParms, "ChargedRatioLTh", "RZ ratio of charged particle fluxes")) {
      RootWImage* anImageChargedRatioLTh = new RootWImage(canvasChargedRatioLTh, canvasChargedRatioLTh->GetWindowWidth(), canvasChargedRatioLTh->GetWindowHeight());
      anImageChargedRatioLTh->setComment("RZ ratio of charged particles fluxes - (Cut_{e prod}=100keV, Cut_{#gamma prod}=10keV)/(Cut_{e prod}=1MeV, Cut_{#gamma prod}=100keV)");
      plotsChargedRatioContent->addItem(anImageChargedRatioLTh);
    }
  }
  if (m_hisChargedFluxBOffMatOn && m_hisChargedFluxBOffMatOff) {
    m_hisChargedRatioMat  = (TH2D*)(m_hisChargedFluxBOffMatOn->Clone("ChargedFluxPerPPRatioMat"));
    m_hisChargedRatioMat->Divide(m_hisChargedFluxBOffMatOff);
    m_hisChargedRatioMat->SetTitle("Ratio of charged particles flux - (Material)/(No material+No mag.field)");

    if (drawHistogram(canvasChargedRatioMat, m_hisChargedRatioMat, m_chargedMapBOffMatOn, simParms, "ChargedRatioMat", "RZ ratio of charged particle fluxes")) {
      RootWImage* anImageChargedRatioMat = new RootWImage(canvasChargedRatioMat, canvasChargedRatioMat->GetWindowWidth(), canvasChargedRatioMat->GetWindowHeight());
      anImageChargedRatioMat->setComment("RZ ratio of charged particles fluxes - (Material)/(No material+No mag.field)");
      plotsChargedRatioContent->addItem(anImageChargedRatioMat);
    }
  }
  if (m_hisChargedFluxBOnMatOff && m_hisChargedFluxBOnMatOff) {
    m_hisChargedRatioB  = (TH2D*)(m_hisChargedFluxBOnMatOff->Clone("ChargedFluxPerPPRatioB"));
    m_hisChargedRatioB->Divide(m_hisChargedFluxBOffMatOff);
    m_hisChargedRatioB->SetTitle("Ratio of charged particles flux - (Mag.field)/(No material+No mag.field)");

    if (drawHistogram(canvasChargedRatioB, m_hisChargedRatioB, m_chargedMapBOnMatOff, simParms, "ChargedRatioB", "RZ ratio of charged particle fluxes")) {
      RootWImage* anImageChargedRatioB = new RootWImage(canvasChargedRatioB, canvasChargedRatioB->GetWindowWidth(), canvasChargedRatioB->GetWindowHeight());
      anImageChargedRatioB->setComment("RZ ratio of charged particles fluxes - (Mag.field)/(No material+No mag.field)");
      plotsChargedRatioContent->addItem(anImageChargedRatioB);
    }
  }
  if (m_hisChargedFluxBOnMatOn && m_hisChargedFluxBOffMatOff) {
    m_hisChargedRatioMatB  = (TH2D*)(m_hisChargedFluxBOnMatOn->Clone("ChargedFluxPerPPRatioMatB"));
    m_hisChargedRatioMatB->Divide(m_hisChargedFluxBOffMatOff);
    m_hisChargedRatioMatB->SetTitle("Ratio of charged particles flux - (Material+Mag.field)/(No material+No mag.field)");

    if (drawHistogram(canvasChargedRatioMatB, m_hisChargedRatioMatB, m_chargedMapBOnMatOn, simParms, "ChargedRatioMatB", "RZ ratio of charged particle fluxes")) {
      RootWImage* anImageChargedRatioMatB = new RootWImage(canvasChargedRatioMatB, canvasChargedRatioMatB->GetWindowWidth(), canvasChargedRatioMatB->GetWindowHeight());
      anImageChargedRatioMatB->setComment("RZ ratio of charged particles fluxes - (Material+Mag.field)/(No material+No mag.field)");
      plotsChargedRatioContent->addItem(anImageChargedRatioMatB);
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

    IrradiationMap* usedChargedMap = nullptr;
    IrradiationMap* usedPhotonsMap = nullptr;

    if (m_photonsMapBOnMatOnLTh!=nullptr) usedPhotonsMap = m_photonsMapBOnMatOnLTh;
    else                                  usedPhotonsMap = m_photonsMapBOnMatOn;
    if (m_photonsMapBOnMatOnLTh!=nullptr) usedChargedMap = m_chargedMapBOnMatOnLTh;
    else                                  usedChargedMap = m_chargedMapBOnMatOn;

    OccupancyVisitor geometryVisitor(usedPhotonsMap, usedChargedMap);
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
