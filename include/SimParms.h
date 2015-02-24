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
  IrradiationMapsManager irradiationMapsManager_;
public:
  
  ReadonlyProperty<int, NoDefault> numMinBiasEvents;
  ReadonlyProperty<int, NoDefault> zErrorCollider;
  ReadonlyProperty<int, NoDefault> rError;
  ReadonlyProperty<bool, NoDefault> useIPConstraint;
//  Property<int> zError;
  ReadonlyProperty<int, NoDefault> ptCost;
  ReadonlyProperty<int, NoDefault> stripCost;
  ReadonlyProperty<double, NoDefault> efficiency;
  ReadonlyProperty<double, NoDefault> pixelEfficiency;

  ReadonlyProperty<int, NoDefault> bunchSpacingNs;

  ReadonlyProperty<double, Default> triggerEtaCut;
  ReadonlyProperty<double, Default> triggerPtCut;
  ReadonlyProperty<int, Default> numTriggerTowersEta, numTriggerTowersPhi;

  ReadonlyProperty<double, Default> timeIntegratedLumi;
  ReadonlyProperty<double, Default> operatingTemp;
  ReadonlyProperty<double, Default> chargeDepletionVoltage;
  ReadonlyProperty<double, Default> alphaParm;
  ReadonlyProperty<double, Default> referenceTemp;

  ReadonlyProperty<double, Default> magneticField;

  PropertyVector<std::string, ','> irradiationMapFiles;
  //std::vector<Property<std::string, NoDefault>> irradiationMapFiles;

  Property<double, NoDefault> minTracksEta, maxTracksEta;

  PropertyNode<std::string> taggedTracking;

  SimParms() : 
      numMinBiasEvents("numMinBiasEvents", parsedAndChecked()),
      zErrorCollider("zErrorCollider", parsedAndChecked()),
      rError("rError", parsedAndChecked()),
      useIPConstraint("useIPConstraint", parsedAndChecked()),
      ptCost("ptCost", parsedAndChecked()),
      stripCost("stripCost", parsedAndChecked()),
      efficiency("efficiency", parsedAndChecked()),
      pixelEfficiency("pixelEfficiency", parsedAndChecked()),
      bunchSpacingNs("bunchSpacingNs", parsedAndChecked()),
      triggerEtaCut("triggerEtaCut", parsedOnly(), 2),
      triggerPtCut("triggerPtCut", parsedOnly(), 1),
      numTriggerTowersEta("numTriggerTowersEta", parsedOnly(), 1),
      numTriggerTowersPhi("numTriggerTowersPhi", parsedOnly(), 1),
      timeIntegratedLumi("timeIntegratedLumi", parsedOnly(), 3000),
      operatingTemp("operatingTemp", parsedOnly(), -20),
      chargeDepletionVoltage("chargeDepletionVoltage", parsedOnly(), 600),
      alphaParm("alphaParm", parsedOnly(), 4e-17),
      referenceTemp("referenceTemp", parsedOnly(), 20),
      magneticField("magneticField", parsedOnly(), 3.8),
      irradiationMapFiles("irradiationMapFiles", parsedAndChecked()),
      //irradiationMapFile("irradiationMapFile", parsedAndChecked()),
      minTracksEta("minTracksEta", parsedOnly()),
      maxTracksEta("maxTracksEta", parsedOnly()),
      taggedTracking("TaggedTracking", parsedOnly())
  { }

  void build();

  void addIrradiationMapFile(std::string path);

  const IrradiationMapsManager& irradiationMapsManager() const { return irradiationMapsManager_; }

  void accept(GeometryVisitor& v) { v.visit(*this); }
  void accept(ConstGeometryVisitor& v) const { v.visit(*this); }

  double calcCost(int) { return 0.; } // CUIDADO FIX THIS does nothing

  double particleCurvatureR(double pt) const { return pt/(0.3*insur::magnetic_field) * 1e3; }
};

#endif
