/*
 * AnalyzerOccupancy.h
 *
 *  Created on: 9. 11. 2015
 *      Author: Drasal (CERN)
 */
#ifndef INCLUDE_ANALYZEROCCUPANCY_H_
#define INCLUDE_ANALYZEROCCUPANCY_H_

#include <memory>
#include <string>
#include <vector>

#include <AnalyzerUnit.h>

class BFieldMap;
class Detector;
class IrradiationMap;
class RootWSite;
class TH2D;
class TCanvas;

/*
 * @class AnalyzerOccupancy
 * Analyze Fluka simulated fluxes of charged hadrons & photons and calculate occupancies, optimal pitch etc.
 * Vizualize data and print them out in a html formatted output.
 */
class AnalyzerOccupancy : public AnalyzerUnit {

 public:

  //! Constructor
  AnalyzerOccupancy(std::string chargedFileName, std::string photonsFileName, const Detector& detector);

  //! Destructor
  virtual ~AnalyzerOccupancy();

  //! Init variables
  bool init(int nTracks) {return true;};

  //! Calculate occupancy & ideal pitch size
  bool analyze();

  //! Visualize - add html page with all calculations & results
  bool visualize(RootWSite& webSite);

  //! Read magnetic field map - just to display
  bool readMagFieldMap(std::string directory, std::string bFieldFileName);
  //! Read no magnetic field irradiation map
  bool readNoMagFieldIrradMap(std::string directory, std::string chargedFileName, std::string photonsFileName);
  //! Read no material irradiation map
  bool readNoMaterialIrradMap(std::string directory, std::string chargedFileName, std::string photonsFileName);
  //! Read no magnetic field & no material irradiation map
  bool readNoMagFieldNoMaterialIrradMap(std::string directory, std::string chargedFileName, std::string photonsFileName);
  //! Read no magnetic field & no tracker irradiation map
  bool readNoMagFieldNoTrackerIrradMap(std::string directory, std::string chargedFileName, std::string photonsFileName);
  //! Read low threshold irradiation map
  bool readLowThresholdIrradMap(std::string directory, std::string chargedFileName, std::string photonsFileName);

 private:

  // Constants
  double c_fluxMin = 0.5E-5;
  double c_fluxMax = 10;

  //! Fill histogram
  bool fillHistogram(const IrradiationMap* map, TH2D*& his, std::string name, std::string title);

  //! Draw histogram
  bool drawHistogram(TCanvas& canvas, TH2D* his, const IrradiationMap* map, std::string nameType, std::string nameParticles);

  // Analyzed geometry
  //std::vector<Tracker*> m_trackers;

  IrradiationMap* m_photonsMapBOnMatOn;   //!< Photons fluxes
  IrradiationMap* m_photonsMapBOnMatOnLTh;//!< Photons low threshold
  IrradiationMap* m_photonsMapBOffMatOn;  //!< Photons fluxes in a detector without magnetic field
  IrradiationMap* m_photonsMapBOnMatOff;  //!< Photons fluxes in a space without detector and with magnetic field
  IrradiationMap* m_photonsMapBOffMatOff; //!< Photons fluxes in a space without detector and without magnetic field
  IrradiationMap* m_photonsMapBOffTrkOff; //!< Photons fluxes in a space without tracker (calorimeter is there) and without magnetic field
  IrradiationMap* m_chargedMapBOnMatOn;   //!< Charged particles fluxes
  IrradiationMap* m_chargedMapBOnMatOnLTh;//!< Photons low threshold
  IrradiationMap* m_chargedMapBOffMatOn;  //!< Charged particles fluxes in a detector without magnetic field
  IrradiationMap* m_chargedMapBOnMatOff;  //!< Charged particles fluxes in a space without detector and with magnetic field
  IrradiationMap* m_chargedMapBOffMatOff; //!< Charged particles fluxes in a space without detector and without magnetic field
  IrradiationMap* m_chargedMapBOffTrkOff; //!< Charged particles fluxes in a space without tracker (calorimeter is there) and without magnetic field

  BFieldMap*      m_bFieldMap;            //!< 3D map of b field inside a detector

  // Visualisation
  TH2D* m_hisChargedFluxBOnMatOn;
  TH2D* m_hisChargedFluxBOnMatOnLTh;
  TH2D* m_hisChargedFluxBOffMatOn;
  TH2D* m_hisChargedFluxBOnMatOff;
  TH2D* m_hisChargedFluxBOffMatOff;
  TH2D* m_hisChargedFluxBOffTrkOff;

  TH2D* m_hisChargedRatioLTh;
  TH2D* m_hisChargedRatioECalMat;
  TH2D* m_hisChargedRatioMat;
  TH2D* m_hisChargedRatioB;
  TH2D* m_hisChargedRatioTrkB;
  TH2D* m_hisChargedRatioMatB;

  TH2D* m_hisPhotonsFluxBOnMatOn;
  TH2D* m_hisPhotonsFluxBOnMatOnLTh;
  TH2D* m_hisPhotonsFluxBOffMatOn;
  TH2D* m_hisPhotonsFluxBOnMatOff;
  TH2D* m_hisPhotonsFluxBOffMatOff;
  TH2D* m_hisPhotonsFluxBOffTrkOff;

  TH2D* m_hisPhotonsRatioLTh;
  TH2D* m_hisPhotonsRatioECalMat;
  TH2D* m_hisPhotonsRatioMat;
  TH2D* m_hisPhotonsRatioB;
  TH2D* m_hisPhotonsRatioTrkB;
  TH2D* m_hisPhotonsRatioMatB;
};

#endif /* INCLUDE_ANALYZEROCCUPANCY_H_ */
