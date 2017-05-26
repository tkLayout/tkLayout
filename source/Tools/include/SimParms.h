#ifndef SIMPARMS_H
#define SIMPARMS_H

#include <map>
#include <string>
#include <sstream>
#include <fstream>

#include "Property.h"
#include "Visitable.h"

// Forward declaration
class GeometryVisitor;
class ConstGeometryVisitor;
class IrradiationMapsManager;

/*
 * @class SimParms
 * @brief A singleton class containing generic information needed across the tkLayout toolkit.
 * @details A singleton class containing generic information needed across the tkLayout toolkit. Various generic
 * parameters (geometry layout name, main geometry config file, etc) as well as generic environment/geometry
 * parameters (not directly related to individual sub-trackers) are held here. Two mechanisms are used to
 * fill-in the SimParms container: the environment variables are read-in from SimParms config file, the rest
 * are saved through setter methods at different level of processing (e.g. after full geometry is read-in).
 */
class SimParms : public PropertyObject, public Buildable, public Visitable {

 public:

  //! SimParms access method -> get instance of singleton class SimParms
  static SimParms& getInstance();
  
  //! Destructor
  ~SimParms();

  //! SimParms visitable -> implemented accept method to call visitor pattern
  void accept(GeometryVisitor& v);

  //! SimParms visitable -> implememented const accept method to call visitor pattern
  void accept(ConstGeometryVisitor& v) const;

  //! Cross-check that Sim parameters correctly read-in from the file & set units
  void crosscheck();

  //! Set command line options passed over to program to analyze data
  //! @param[in] commandLine    Content of command line
  void setCommandLine(int argc, char* argv[]);

  //! Set geometry layout name
  //! @param[in] layout name
  void setLayoutName(std::string layoutName) {m_layoutName = layoutName;}

  //! Set default html directory, where all results are saved
  //! @param[in] directory address
  void setWebDir(std::string htmlDir) {m_htmlDir = htmlDir;}

  //! Set name of main geometry config file
  //! @param[in] file name
  void setBaseGeomFileName(std::string geomFileName) {m_geomFile = geomFileName;}

  //! Set path of tkLayout run directory
  //! @param[in] run directory path
  void setRunDirPath(std::string baseDir) {m_baseDir = baseDir;}

  //! Set list of all configuration files obtained from the base geometry file using @include command
  //! @param [in] set of all include files (strings)
  void setListOfConfFiles(std::set<std::string> includeSet) {m_includeSet = includeSet;}

  //! Get command line content
  //! @return line content as one string
  std::string getCommandLine() const {return m_commandLine;}

  //! Get geometry layout name
  //! @return layout name
  std::string getLayoutName() const {return m_layoutName;}

  //! Get default html directory, where all results are saved
  //! @return directory address
  std::string getWebDir() const {return m_htmlDir;}

  //! Get name of main geometry config file
  //! @return file name
  std::string getBaseGeomFileName() const {return m_geomFile;}

  //! Get path of tkLayout run directory
  //! @return run directory path
  std::string getRunDirPath() const {return m_baseDir;}

  //! Get list of all configuration files obtained from the base geometry file using @include command
  //! @return set of all include files (strings)
  std::set<std::string> getListOfConfFiles() const {return m_includeSet;}

  //! Get number of defined eta regions
  size_t getNEtaRegions() const { return etaRegionRanges.size(); }

  //! Get eta max value as defined by user in definition of various tracker eta regions
  double getMaxEtaCoverage() const { if (etaRegionRanges.size()>0) return etaRegionRanges[etaRegionRanges.size()-1]; else return 0.0; }

  //! Check whether magnetic field const. or defined as a function
  bool isMagFieldConst() const { if (magField.size()==1) return true; else return false;}

  //! Get number of defined mag. field regions
  size_t getNMagFieldRegions() const { return magFieldZRegions.size(); }

  //! Get reference to irradiation maps manager
  const IrradiationMapsManager& irradiationMapsManager() const { return *m_irradiationMapsManager;}

  // Variables to be read in by SimParms class from a configuration file
  ReadonlyProperty<int   , NoDefault> numMinBiasEvents;     //!< Number of minimum bias events in p-p collision (#pile-ups)
  ReadonlyProperty<double, NoDefault> zErrorIP;             //!< Defines typical size of luminous region in Z, typically sigmaZ_bunch/sqrt(2)
  ReadonlyProperty<double, NoDefault> rphiErrorIP;          //!< Defines typical size of luminous region in R-Phi
  ReadonlyProperty<bool  , Default>   useLumiRegInGeomBuild;//!< Apply luminous region constraints, when building geometry (build hermetic detector using z+-zErrorIP, rphiErrorIP negligible)
  ReadonlyProperty<bool  , Default>   useLumiRegInAnalysis; //!< Apply luminous region constraints, when analysing geometry (resolution, pattern reco, etc.)
  ReadonlyProperty<bool  , NoDefault> useIPConstraint;      //!< Use IP constraint in track fitting
  ReadonlyProperty<int   , NoDefault> ptCost;
  ReadonlyProperty<int   , NoDefault> stripCost;
  ReadonlyProperty<double, NoDefault> efficiency;
  ReadonlyProperty<double, NoDefault> pixelEfficiency;

  ReadonlyProperty<int   , NoDefault> bunchSpacingNs;

  ReadonlyProperty<double, Default>   triggerEtaCut;
  ReadonlyProperty<double, Default>   triggerPtCut;
  ReadonlyProperty<int   , Default>   numTriggerTowersEta, numTriggerTowersPhi;

  ReadonlyProperty<double, Default>   timeIntegratedLumi;
  ReadonlyProperty<double, Default>   operatingTemp;
  ReadonlyProperty<double, Default>   chargeDepletionVoltage;
  ReadonlyProperty<double, Default>   alphaParm;
  ReadonlyProperty<double, Default>   referenceTemp;

  ReadonlyPropertyVector<double, ','> magField;             //!< Magnetic field [T] as a vector of approximate values of B(z) in predefined regions 0-z0, z0-z1, ...
  ReadonlyPropertyVector<double, ','> magFieldZRegions;     //!< Regions [m], in which the magnetic field values (B(z)) are defined, regions defined as a sequence of z0 [m], z1 [m], ... -> intervals 0-z0; z0-z1, etc.

  PropertyVector<std::string, ','>    irradiationMapFiles;

  // Define eta regions & region names to be plotted when drawing geometry.
  // The last value of region ranges represents the maximum tracker eta coverage -> used by analysis modules
  ReadonlyPropertyVector<double     , ','> etaRegionRanges; //!< Set ordered eta regions, the last one being maximum tracker eta coverage (e.g. 0, 2.5, 4.0)
  ReadonlyPropertyVector<std::string, ','> etaRegionNames;  //!< Set names for ordered eta borders (e.g TRK-0, TRK-BRL, TRK-MAX)

  Property<std::string, Default> bFieldMapFile;           // Map of b field - currently forvisualization use only (not for tracking)
  Property<std::string, Default> chargedMapFile;          // Map of charged hadron fluxes, segmented in R x Z
  Property<std::string, Default> chargedMapLowThFile;     // Map of charged hadron fluxes - electrons with lower threshold, segmented in R x Z
  Property<std::string, Default> chargedMapBOffMatOnFile; // Map of charged hadron fluxes, segmented in R x Z, no mag. field applied
  Property<std::string, Default> chargedMapBOnMatOffFile; // Map of charged hadron fluxes, segmented in R x Z, no material
  Property<std::string, Default> chargedMapBOffMatOffFile;// Map of charged hadron fluxes, segmented in R x Z, no mag. field applied, no material
  Property<std::string, Default> chargedMapBOffTrkOffFile;// Map of charged hadron fluxes, segmented in R x Z, no mag. field applied, no tracker material
  Property<std::string, Default> photonsMapFile;          // Map of photon fluxes, segmented in R x Z
  Property<std::string, Default> photonsMapLowThFile;     // Map of photon fluxes - electrons with lower threshold, segmented in R x Z
  Property<std::string, Default> photonsMapBOffMatOnFile; // Map of photon fluxes, segmented in R x Z, no mag. field applied
  Property<std::string, Default> photonsMapBOnMatOffFile; // Map of photon fluxes, segmented in R x Z, no material
  Property<std::string, Default> photonsMapBOffMatOffFile;// Map of photon fluxes, segmented in R x Z, no mag. field applied, no material
  Property<std::string, Default> photonsMapBOffTrkOffFile;// Map of photon fluxes, segmented in R x Z, no mag. field applied, no tracker material

  Property<        double, NoDefault> minTracksEta, maxTracksEta;
  PropertyNode<std::string>           taggedTracking;

 private:

  //! Singleton private constructor -> initialize all variables to be read-out from configuration file & define if checker should be called after parsing
  SimParms();

  std::string            m_commandLine;  //!< Command line options passed over to program to analyze data
  std::string            m_geomFile;     //!< Name of main geometry config file
  std::string            m_baseDir;      //!< Base tkLayout run directory
  std::string            m_htmlDir;      //!< Default html directory, where all results are saved
  std::string            m_layoutName;   //!< Geometry layout name
  std::set<std::string>  m_includeSet;   //!< List of all configuration files obtained from the base geometry file using @include command

  std::unique_ptr<IrradiationMapsManager> m_irradiationMapsManager; //!< Irradiation maps manager



};

#endif
