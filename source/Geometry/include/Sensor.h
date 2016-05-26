#ifndef SENSOR_H
#define SENSOR_H

#include <string>
#include <exception>

#include "global_funcs.h"
#include "Polygon3d.h"
#include "Property.h"
#include "CoordinateOperations.h"

class DetectorModule;

enum class SensorType { Pixel, Largepix, Strip, None };

class Sensor : public PropertyObject, public Buildable, public Identifiable<int> {
  const DetectorModule* parent_;
  mutable const Polygon3d<4>* hitPoly_ = 0; 
  mutable const Polygon3d<4>* envPoly_ = 0; 
  Polygon3d<4>* buildOwnPoly(double polyOffset) const;
public:
  ReadonlyProperty<int, NoDefault> numSegments;
  ReadonlyProperty<int, NoDefault> numStripsAcross;
  ReadonlyProperty<int, NoDefault> numROCX, numROCY;
  ReadonlyProperty<double, Default> sensorThickness;
  ReadonlyProperty<SensorType, Default> type;
  ReadonlyProperty<double, UncachedComputable> minR, maxR; // CUIDADO min/maxR don't take into account the sensor thickness!
  ReadonlyProperty<double, UncachedComputable> minZ, maxZ; // ditto for min/maxZ

  Sensor() : 
      numSegments("numSegments", parsedOnly()),
      numStripsAcross("numStripsAcross", parsedOnly()),
      numROCX("numROCX", parsedOnly()),
      numROCY("numROCY", parsedOnly()),
      sensorThickness("sensorThickness", parsedOnly(), 0.1),
      type("sensorType", parsedOnly(), SensorType::None)
  {}

  void parent(const DetectorModule* m) { parent_ = m; }

  int numChannels() const { return numStripsAcross() * numSegments(); }
  double minPitch() const;
  double maxPitch() const;
  double pitch() const;
  double stripLength() const;

  int numROCRows() const { return numStripsAcross() / numROCX(); } 
  int numROCCols() const { return numSegments() / numROCY(); }

  int totalROCs() const { return numROCX() * numROCY(); }

  double normalOffset() const;

  std::pair<XYZVector, int> checkHitSegment(const XYZVector& trackOrig, const XYZVector& trackDir) const;

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
};

#endif
