#ifndef SIMPARMS_H
#define SIMPARMS_H

#include <map>
#include <string>
#include <sstream>
#include <fstream>

#include "global_constants.hh"
#include "MessageLogger.hh"
#include "Property.hh"
#include "capabilities.hh"
#include "Visitor.hh"
#include "IrradiationMapsManager.hh"
#include "Visitable.hh"

enum LumiRegShape { SPOT, FLAT, GAUSSIAN };

class SimParms : public PropertyObject, public Buildable, public Visitable {

public:
  
  //! SimParms access method -> get instance of singleton class SimParms
  static SimParms& getInstance();

  ReadonlyProperty<int, NoDefault> numMinBiasEvents;
  ReadonlyProperty<int, NoDefault> bunchSpacingNs;

  ReadonlyProperty<double, NoDefault> lumiRegZError;
  ReadonlyProperty<LumiRegShape, NoDefault> lumiRegShape;
  ReadonlyProperty<LumiRegShape, NoDefault> lumiRegShapeInMatBudgetAnalysis;

  ReadonlyProperty<bool, NoDefault> useIPConstraint;
  ReadonlyProperty<double, NoDefault> rphiErrorCollider;

  ReadonlyProperty<int, NoDefault> ptCost;
  ReadonlyProperty<int, NoDefault> stripCost;

  ReadonlyProperty<double, NoDefault> triggerEtaCut;
  ReadonlyProperty<double, NoDefault> triggerPtCut;
  ReadonlyProperty<int, NoDefault> numTriggerTowersEta, numTriggerTowersPhi;

  ReadonlyProperty<double, NoDefault> timeIntegratedLumi;
  ReadonlyProperty<double, NoDefault> referenceTemp;
  ReadonlyProperty<double, NoDefault> alphaParam;
  ReadonlyProperty<double, NoDefault> magField;

  PropertyVector<std::string, ','> irradiationMapFiles;

  Property<double, NoDefault> minTracksEta, maxTracksEta;

  PropertyNode<std::string> taggedTracking;

  void build();

  void addIrradiationMapFile(std::string path);

  //! Check whether magnetic field const. or defined as a function -> for now assumed const.
  bool isMagFieldConst() const { return true;}

  //! Get number of defined mag. field regions
  //size_t getNMagFieldRegions() const { return magFieldZRegions.size(); }

  //! Get reference to irradiation maps manager
  const IrradiationMapsManager& irradiationMapsManager() const { return irradiationMapsManager_; }

  void accept(GeometryVisitor& v) { v.visit(*this); }
  void accept(ConstGeometryVisitor& v) const { v.visit(*this); }

  double particleCurvatureR(double pt) const { return pt/(0.3*magField()) * 1e3; }

private:

  //! Singleton private constructor -> initialize all variables to be read-out from configuration file & define if checker should be called after parsing
  SimParms();

  IrradiationMapsManager irradiationMapsManager_;
};

#endif
