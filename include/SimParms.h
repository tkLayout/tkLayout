#ifndef SIMPARMS_H
#define SIMPARMS_H

#include <map>
#include <string>
#include <sstream>
#include <fstream>

#include "global_constants.h"
#include "messageLogger.h"
#include "Property.h"
#include "capabilities.h"
#include "Visitor.h"
#include "IrradiationMapsManager.h"
#include "Visitable.h"

//typedef std::map<std::pair<int,int>, double> IrradiationMap;

class SimParms : public PropertyObject, public Buildable, public Visitable {

public:

  //! SimParms access method -> get instance of singleton class SimParms
  static SimParms* getInstance();
  
  //! Destructor
  ~SimParms() {};

  //! SimParms visitable -> implemented accept method to call visitor pattern
  void accept(GeometryVisitor& v) { v.visit(*this); }

  //! SimParms visitable -> implememented const accept method to call visitor pattern
  void accept(ConstGeometryVisitor& v) const { v.visit(*this); }

  //! Cross-check that Sim parameters correctly read-in from the file
  void crosscheck();

  //! Read-in irradiation maps
  void readIrradiationMaps();

  //! Get reference to irradiation maps manager
  const IrradiationMapsManager& irradiationMapsManager() const { return m_irradiationMapsManager;}

  // Variables to be read in by SimParms class from a configuration file
  ReadonlyProperty<int   , NoDefault> numMinBiasEvents;
  ReadonlyProperty<int   , NoDefault> zErrorCollider;
  ReadonlyProperty<int   , NoDefault> rError;
  ReadonlyProperty<bool  , NoDefault> useIPConstraint;
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

  ReadonlyProperty<double, Default>   magneticField;      // Magnetic field in Tesla

  // To include dipole region in a quick way
  ReadonlyProperty<double, Default>   dipoleMagneticField;// Integral magnetic field in [Tm]
  ReadonlyProperty<double, Default>   dipoleDPlResAt10TeV;// deltaPl/Pl resolution of dipole tracker at 10 TeV
  ReadonlyProperty<double, Default>   dipoleXToX0;        // [%]

  // Beam pipe radius, thickness, thickness in rad. length, in int. length
  ReadonlyProperty<double, Default>   bpRadius;
  ReadonlyProperty<double, Default>   bpThickness;
  ReadonlyProperty<double, Default>   bpRadLength;        // [%]
  ReadonlyProperty<double, Default>   bpIntLength;        // [%]

  PropertyVector<std::string, ','>    irradiationMapFiles;
  //std::vector<Property<std::string, NoDefault>> irradiationMapFiles;

  Property<std::string, Default> bFieldMapFile;           // Map of b field - not currently currently for tracking
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

  //! An instance of SimParms class - singleton pattern
  static SimParms * s_instance;

  //! Singleton private constructor -> initialize all variables to be read-out from configuration file & define if checker should be called after parsing
  SimParms() : 
      numMinBiasEvents(        "numMinBiasEvents"        , parsedAndChecked()),
      zErrorCollider(          "zErrorCollider"          , parsedAndChecked()),
      rError(                  "rError"                  , parsedAndChecked()),
      useIPConstraint(         "useIPConstraint"         , parsedAndChecked()),
      ptCost(                  "ptCost"                  , parsedAndChecked()),
      stripCost(               "stripCost"               , parsedAndChecked()),
      efficiency(              "efficiency"              , parsedAndChecked()),
      pixelEfficiency(         "pixelEfficiency"         , parsedAndChecked()),
      bunchSpacingNs(          "bunchSpacingNs"          , parsedAndChecked()),
      triggerEtaCut(           "triggerEtaCut"           , parsedOnly(), 2),
      triggerPtCut(            "triggerPtCut"            , parsedOnly(), 1),
      numTriggerTowersEta(     "numTriggerTowersEta"     , parsedOnly(), 1),
      numTriggerTowersPhi(     "numTriggerTowersPhi"     , parsedOnly(), 1),
      timeIntegratedLumi(      "timeIntegratedLumi"      , parsedOnly(), 3000),
      operatingTemp(           "operatingTemp"           , parsedOnly(), -20),
      chargeDepletionVoltage(  "chargeDepletionVoltage"  , parsedOnly(), 600),
      alphaParm(               "alphaParm"               , parsedOnly(), 4e-17),
      referenceTemp(           "referenceTemp"           , parsedOnly(), 20),
      magneticField(           "magneticField"           , parsedOnly(), insur::magnetic_field),
      dipoleMagneticField(     "dipoleMagneticField"     , parsedOnly(), 0.0),
      dipoleDPlResAt10TeV(     "dipoleDPlResAt10TeV"     , parsedOnly(), 0.1),
      dipoleXToX0(             "dipoleXToX0"             , parsedOnly(), 0.1),
      bpRadius(                "beamPipeRadius"          , parsedOnly(), 0.0),
      bpThickness(             "beamPipeThickness"       , parsedOnly(), 0.0),
      bpRadLength(             "beamPipeRadLength"       , parsedOnly(), 0.0),
      bpIntLength(             "beamPipeIntLength"       , parsedOnly(), 0.0),
      irradiationMapFiles(     "irradiationMapFiles"     , parsedAndChecked()),
      bFieldMapFile(           "bFieldMapFile"           , parsedOnly(), std::string("")),
      chargedMapFile(          "chargedMapFile"          , parsedOnly(), std::string("")),
      chargedMapLowThFile(     "chargedMapLowThFile"     , parsedOnly(), std::string("")),
      chargedMapBOffMatOnFile( "chargedMapBOffMatOnFile" , parsedOnly(), std::string("")),
      chargedMapBOnMatOffFile( "chargedMapBOnMatOffFile" , parsedOnly(), std::string("")),
      chargedMapBOffMatOffFile("chargedMapBOffMatOffFile", parsedOnly(), std::string("")),
      chargedMapBOffTrkOffFile("chargedMapBOffTrkOffFile", parsedOnly(), std::string("")),
      photonsMapFile(          "photonsMapFile"          , parsedOnly(), std::string("")),
      photonsMapLowThFile(     "photonsMapLowThFile"     , parsedOnly(), std::string("")),
      photonsMapBOffMatOnFile( "photonsMapBOffMatOnFile" , parsedOnly(), std::string("")),
      photonsMapBOnMatOffFile( "photonsMapBOnMatOffFile" , parsedOnly(), std::string("")),
      photonsMapBOffMatOffFile("photonsMapBOffMatOffFile", parsedOnly(), std::string("")),
      photonsMapBOffTrkOffFile("photonsMapBOffTrkOffFile", parsedOnly(), std::string("")),
      minTracksEta(            "minTracksEta"            , parsedOnly()),
      maxTracksEta(            "maxTracksEta"            , parsedOnly()),
      taggedTracking(          "TaggedTracking"          , parsedOnly())
  {}

  // Irradiation maps manager
  IrradiationMapsManager m_irradiationMapsManager;

};

#endif
