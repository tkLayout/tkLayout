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

 private:

  // Analyzed geometry
  std::vector<Tracker*> m_trackers;
  // Photons fluxes
  IrradiationMap* m_photonsMap;
  // Charged particles fluxes
  IrradiationMap* m_chargedMap;

  // Visualisation
  TH2D* m_hisChargedFlux;
  TH2D* m_hisPhotonsFlux;
};

#endif /* INCLUDE_ANALYZEROCCUPANCY_H_ */
