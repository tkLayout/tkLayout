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
#include <Visitor.h>

class BFieldMap;
class ConstGeometryVisitor;
class Detector;
class IrradiationMap;
class RootWSite;
class RootWTable;
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

//
//! Helper class -> Analyzer occupancy using visitor pattern class
//
class OccupancyVisitor : public ConstGeometryVisitor {

 public:

  OccupancyVisitor(IrradiationMap* photonsMap, IrradiationMap* chargedMap);

  virtual ~OccupancyVisitor() {};

  std::unique_ptr<RootWTable> getLayerTable(signed int nPileUps, std::string trkName);
  std::unique_ptr<RootWTable> getRingTable(signed int nPileUps, std::string trkName);

  void setMaxPileUp(double maxPileUp)   { m_maxPileUp = maxPileUp;}
  void setMaxColFreq(double maxColFreq) { m_maxColFreq= maxColFreq;}

  void visit(const Layer& layer) override;
  void visit(const BarrelModule& module) override;
  void visit(const Disk& disk) override;
  void visit(const Ring& ring) override;
  void visit(const EndcapModule& module) override;

 private:

  const bool c_assumeFlowsFromIP = true;  // Assume that all particles come from the interaction point

  IrradiationMap* m_photonsMap;
  IrradiationMap* m_chargedMap;
  int             m_nLayers;
  int             m_nDisks;
  int             m_nRings;
  int             m_iRing;

  std::vector<double>             m_layerRadii;             // Radius of a given layer
  std::vector<double>             m_layerMinFluxes;         // Minimum flux in a layer
  std::vector<double>             m_layerMaxFluxes;         // Maximum flux in a layer
  std::vector<double>             m_layerMaxFluxZ;          // Z-pos of the module in a layer with a maximum flux
  std::vector<std::vector<long>>  m_layerNChannels;         // Number of channels to be read-out in each layer - separately for each sensor type (either 1 or 2 types)
  std::vector<std::vector<long>>  m_layerNHits;             // Number of hits to be read-out in each layer - separately for each sensor type (either 1 or 2 types)
  std::vector<int>                m_layerNRods;             // Number of rods in each layer
  std::vector<int>                m_layerNModules;          // Number of modules in each layer
  std::vector<short>              m_layerNSensorsInMod;     // Number of sensors in each module -> assuming all modules are of the same type in a layer (max 2 sensors)
  std::vector<std::vector<int>>   m_layerSenAddrSparSize;   // Channel address size in bits (assuming max 2 types of sensors in the layer) -> used for sparsified data
  std::vector<std::vector<int>>   m_layerSenAddrUnsparSize; // Sensor addressing size in bits (assuming max 2 types of sensors in the layer, each channel is either 0 or 1 - n channels x 1b) -> used for unsparsified data
  std::vector<std::vector<int>>   m_layerSenNPixels;        // Number of pixels (strips) in each module (assuming max 2 types of sensors in the layer)
  std::vector<std::vector<double>>m_layerSenPixelArea;      // Pixel (strips) area in each module calculated from resolutions & assuming binary read-out
  std::vector<std::vector<double>>m_layerSenArea;           // Sensor area

  std::vector<double>             m_ringAvgRadii;
  std::vector<double>             m_ringMinFluxes;
  std::vector<double>             m_ringMaxFluxes;
  std::vector<double>             m_ringMaxFluxZ;
  std::vector<std::vector<long>>  m_ringNChannels;          // Number of channels to be read-out in each ring - separately for each sensor type (either 1 or 2 types)
  std::vector<std::vector<long>>  m_ringNHits;              // Number of hits to be read-out in each ring - separately for each sensor type (either 1 or 2 types)
  std::vector<int>                m_ringNModules;           // Number of modules in each ring
  std::vector<short>              m_ringNSensorsInMod;      // Number of sensors in each module -> assuming all modules are of the same type in a ring (max 2 sensors)
  std::vector<std::vector<int>>   m_ringSenAddrSparSize;    // Channel address size in bits (assuming max 2 types of sensors in the ring) -> used for sparsified data
  std::vector<std::vector<int>>   m_ringSenAddrUnsparSize;  // Sensor addressing size in bits (assuming max 2 types of sensors in the ring, each channel is either 0 or 1 - n channels x 1b) -> used for unsparsified data
  std::vector<std::vector<int>>   m_ringSenNPixels;         // Number of pixels (strips) in each module (assuming max 2 types of sensors in the ring)
  std::vector<std::vector<double>>m_ringSenPixelArea;       // Pixel (strips) area in each module calculated from resolutions & assuming binary read-out
  std::vector<std::vector<double>>m_ringSenArea;            // Sensor area

  double    m_maxPileUp;
  double    m_maxColFreq;

  double    m_zPosStep;
  double    m_rPosStep;
  const int c_coordPrecision= 1;

};

#endif /* INCLUDE_ANALYZEROCCUPANCY_H_ */
