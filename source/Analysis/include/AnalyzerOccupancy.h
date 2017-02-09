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
  AnalyzerOccupancy(const Detector& detector);

  //! Destructor
  virtual ~AnalyzerOccupancy();

  //! Init variables
  bool init(int nTracks);

  //! Calculate occupancy & ideal pitch size
  bool analyze();

  //! Visualize - add html page with all calculations & results
  bool visualize(RootWSite& webSite);

 private:

  // Constants
  double c_fluxMin = 0.5E-5;
  double c_fluxMax = 10;

  //! Check that a file can be opened
  bool checkFile(const std::string& fileName, const std::string& filePath);

  //! Fill histogram
  bool fillHistogram(const IrradiationMap* map, TH2D*& his, std::string name, std::string title);

  //! Draw histogram
  bool drawHistogram(TCanvas& canvas, TH2D* his, const IrradiationMap* map, std::string nameType, std::string nameParticles);

  IrradiationMap* m_photonsMap;   //!< Photons fluxes
  IrradiationMap* m_chargedMap;   //!< Charged particles fluxes

  BFieldMap*      m_bFieldMap;     //!< 3D map of b field inside a detector

  // Visualisation
  TH2D* m_hisChargedFlux;
  TH2D* m_hisPhotonsFlux;
};

#endif /* INCLUDE_ANALYZEROCCUPANCY_H_ */
