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
  ReadonlyProperty<int, NoDefault> bunchSpacingNs;

  ReadonlyProperty<int, NoDefault> zErrorCollider;
  ReadonlyProperty<int, NoDefault> rError;
  ReadonlyProperty<bool, NoDefault> useIPConstraint;

  ReadonlyProperty<int, NoDefault> ptCost;
  ReadonlyProperty<int, NoDefault> stripCost;

  ReadonlyProperty<double, NoDefault> efficiency;
  ReadonlyProperty<double, NoDefault> pixelEfficiency;

  ReadonlyProperty<double, NoDefault> triggerEtaCut;
  ReadonlyProperty<double, NoDefault> triggerPtCut;
  ReadonlyProperty<int, NoDefault> numTriggerTowersEta, numTriggerTowersPhi;

  ReadonlyProperty<double, NoDefault> timeIntegratedLumi;
  ReadonlyProperty<double, NoDefault> operatingTemp;
  ReadonlyProperty<double, NoDefault> referenceTemp;
  ReadonlyProperty<double, NoDefault> chargeDepletionVoltage;
  ReadonlyProperty<double, NoDefault> alphaParam;
  ReadonlyProperty<double, NoDefault> magneticField;

  PropertyVector<std::string, ','> irradiationMapFiles;
  //std::vector<Property<std::string, NoDefault>> irradiationMapFiles;

  Property<double, NoDefault> minTracksEta, maxTracksEta;

  PropertyNode<std::string> taggedTracking;

  SimParms() : 
      numMinBiasEvents("numMinBiasEvents", parsedAndChecked()),
      bunchSpacingNs("bunchSpacingNs", parsedAndChecked()),
      zErrorCollider("zErrorCollider", parsedAndChecked()),
      rError("rError", parsedAndChecked()),
      useIPConstraint("useIPConstraint", parsedAndChecked()),
      ptCost("ptCost", parsedAndChecked()),
      stripCost("stripCost", parsedAndChecked()),
      efficiency("efficiency", parsedAndChecked()),
      pixelEfficiency("pixelEfficiency", parsedAndChecked()), 
      triggerEtaCut("triggerEtaCut", parsedAndChecked()),
      triggerPtCut("triggerPtCut", parsedAndChecked()),
      numTriggerTowersEta("numTriggerTowersEta", parsedAndChecked()),
      numTriggerTowersPhi("numTriggerTowersPhi", parsedAndChecked()),
      timeIntegratedLumi("timeIntegratedLumi", parsedAndChecked()),
      operatingTemp("operatingTemp", parsedAndChecked()),
      referenceTemp("referenceTemp", parsedAndChecked()),
      chargeDepletionVoltage("chargeDepletionVoltage", parsedAndChecked()),
      alphaParam("alphaParam", parsedAndChecked()),    // radiation-damage coefficient, A/cm
      magneticField("magneticField", parsedAndChecked()),
      irradiationMapFiles("irradiationMapFiles", parsedAndChecked()),
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
