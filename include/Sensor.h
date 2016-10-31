#ifndef SENSOR_H
#define SENSOR_H

#include <string>
#include <exception>

#include "global_funcs.h"
#include "Polygon3d.h"
#include "Property.h"
#include "CoordinateOperations.h"
#include "Visitor.h"


enum ModuleSubdetector { BARREL = 1, ENDCAP = 2 };
enum SensorPosition { NO, LOWER, UPPER };
enum class SensorType { Pixel, Largepix, Strip, None };

class DetectorModule;

class Sensor : public PropertyObject, public Buildable, public Identifiable<int>, public DetIdentifiable {
  const DetectorModule* parent_;
  ModuleSubdetector subdet_;
  SensorPosition innerOuter_ = SensorPosition::NO;
  mutable const Polygon3d<4>* hitPoly_ = 0; 
  mutable const Polygon3d<4>* envPoly_ = 0; 
  Polygon3d<4>* buildOwnPoly(double polyOffset) const;
public:
  ReadonlyProperty<int, NoDefault> numStripsAcross;
  ReadonlyProperty<double, NoDefault> pitchEstimate;
  ReadonlyProperty<int, NoDefault> numSegments;
  ReadonlyProperty<double, NoDefault> stripLengthEstimate;
  ReadonlyProperty<int, NoDefault> numROCX, numROCY;
  ReadonlyProperty<double, Default> sensorThickness;
  ReadonlyProperty<SensorType, Default> type;
  ReadonlyProperty<double, Computable> minR, maxR; // CUIDADO min/maxR don't take into account the sensor thickness!
  ReadonlyProperty<double, Computable> minZ, maxZ; // ditto for min/maxZ

 Sensor() :
  numStripsAcross("numStripsAcross", parsedOnly()),
    pitchEstimate("pitchEstimate", parsedOnly()),
    numSegments("numSegments", parsedOnly()),
    stripLengthEstimate("stripLengthEstimate", parsedOnly()),
    numROCX("numROCX", parsedOnly()),
    numROCY("numROCY", parsedOnly()),
    sensorThickness("sensorThickness", parsedOnly(), 0.1),
    type("sensorType", parsedOnly(), SensorType::None)
      {}

  void parent(const DetectorModule* m) { parent_ = m; }

  ModuleSubdetector subdet(ModuleSubdetector s) { subdet_ = s; }
  ModuleSubdetector subdet() const { return subdet_; }
  SensorPosition innerOuter(SensorPosition pos) { innerOuter_ = pos; }
  SensorPosition innerOuter() const { return innerOuter_; }

  int numStripsAcrossEstimate() const;
  int numSegmentsEstimate() const;
  int numChannels() const { return numStripsAcrossEstimate() * numSegmentsEstimate(); }
  double minPitch() const;
  double maxPitch() const;
  double pitch() const;
  double stripLength() const;

  int numROCRows() const { return numStripsAcrossEstimate() / numROCX(); } 
  int numROCCols() const { return numSegmentsEstimate() / numROCY(); }

  int totalROCs() const { return numROCX() * numROCY(); }

  double normalOffset() const;

  std::pair<XYZVector, int> checkHitSegment(const XYZVector& trackOrig, const XYZVector& trackDir) const;

  void check() override;

  void build() { 
    try { check(); } 
    catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }
    cleanup(); 
  }

  void setup() {
    clearPolys();
    minR.setup([&]() { return CoordinateOperations::computeMinR(envelopePoly()); });
    maxR.setup([&]() { return CoordinateOperations::computeMaxR(envelopePoly()); });
    minZ.setup([&]() { return CoordinateOperations::computeMinZ(envelopePoly()); });
    maxZ.setup([&]() { return CoordinateOperations::computeMaxZ(envelopePoly()); });
    }


//  double minRVertex() const { double min = std::numeric_limits<double>::max(); for (auto v : *poly_) { min = MIN(min, v.Rho()); } return min; }
//  double minZVertex() const { double min = std::numeric_limits<double>::max(); for (auto v : *poly_) { min = MIN(min, v.Z()); } return min; }
  
  void clearPolys();
  const Polygon3d<4>& hitPoly() const;
  const Polygon3d<4>& envelopePoly() const;

  void accept(SensorGeometryVisitor& v) { 
    v.visit(*this);
  }
};

#endif
