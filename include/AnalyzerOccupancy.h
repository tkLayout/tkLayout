/*
 * AnalyzerOccupancy.h
 *
 *  Created on: 9. 11. 2015
 *      Author: Drasal (CERN)
 */
#ifndef INCLUDE_ANALYZEROCCUPANCY_H_
#define INCLUDE_ANALYZEROCCUPANCY_H_

#include <string>
#include <vector>

class RootWSite;
class Tracker;
class IrradiationMap;
class TH2D;
class SimParms;

/*
 * Analyze Fluka simulated fluxes of charged hadrons & photons and calculate occupancies, optimal pitch etc.
 * Vizualize data and print them out in a html formatted output.
 */
class AnalyzerOccupancy {

 public:
  // Constructor
  AnalyzerOccupancy(std::string chargedFileName, std::string photonsFileName, std::vector<Tracker*> trackers);
  // Destructor
  ~AnalyzerOccupancy();

  // Calculate occupancy & ideal pitch size
  bool calculate(double etaStep);
  // Visualize - add html page with all calculations & results
  bool visualize(RootWSite& webSite, const SimParms* simParms);
  // Read no magnetic field map
  void readNoMagFieldMap(std::string chargedFileName, std::string photonsFileName);
  // Read no magnetic field & no material map
  void readNoMagFieldNoMaterialMap(std::string chargedFileName, std::string photonsFileName);

 private:

  // Fill histogram
  bool fillHistogram(const IrradiationMap* map, TH2D*& his, std::string name, std::string title);

  // Analyzed geometry
  std::vector<Tracker*> m_trackers;

  IrradiationMap* m_photonsMap;         // Photons fluxes
  IrradiationMap* m_photonsMapNoB;      // Photons fluxes in a detector without magnetic field
  IrradiationMap* m_photonsMapNoBNoMat; // Photons fluxes in a space without detector and without magnetic field
  IrradiationMap* m_chargedMap;         // Charged particles fluxes
  IrradiationMap* m_chargedMapNoB;      // Charged particles fluxes in a detector without magnetic field
  IrradiationMap* m_chargedMapNoBNoMat; // Charged particles fluxes in a space without detector and without magnetic field

  // Visualisation
  TH2D* m_hisChargedFlux;
  TH2D* m_hisPhotonsFlux;

  TH2D* m_hisChargedNoBFlux;
  TH2D* m_hisChargedNoBNoMatFlux;
  TH2D* m_hisPhotonsNoBFlux;
  TH2D* m_hisPhotonsNoBNoMatFlux;

  TH2D* m_hisChargedRatioMat;
  TH2D* m_hisPhotonsRatioMat;
  TH2D* m_hisChargedRatioMatB;
  TH2D* m_hisPhotonsRatioMatB;
};

#endif /* INCLUDE_ANALYZEROCCUPANCY_H_ */
