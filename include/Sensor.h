#include <string>

#include "global_funcs.h"
#include "Polygon3d.h"
#include "Property.h"

enum class SensorType { Pixel, Largepix, Strip, None };

class Sensor : public PropertyObject, public Buildable, public Identifiable<int> {
  Polygon3d<4>* poly_ = 0;
public:
  ReadonlyProperty<int, NoDefault> numSegments;
  ReadonlyProperty<int, NoDefault> numStripsAcross;
  ReadonlyProperty<int, NoDefault> numROCX, numROCY;
  ReadonlyProperty<double, Default> sensorThickness;
  ReadonlyProperty<SensorType, Default> type;
  Property<double, NoDefault> length, minWidth, maxWidth;

  Sensor() : 
      numSegments("numSegments", parsedOnly()),
      numStripsAcross("numStripsAcross", parsedOnly()),
      numROCX("numROCX", parsedOnly()),
      numROCY("numROCY", parsedOnly()),
      sensorThickness("sensorThickness", parsedOnly(), 0.1),
      type("sensorType", parsedOnly(), SensorType::None),
      length("l", checkedOnly()),
      minWidth("w", checkedOnly()),
      maxWidth("W", checkedOnly())
  {}

  int numChannels() const { return numStripsAcross() * numSegments(); }
  double minPitch() const { return minWidth() / (double)numStripsAcross(); }
  double maxPitch() const { return maxWidth() / (double)numStripsAcross(); }
  double pitch() const { return (maxWidth() + minWidth()) / 2. / (double)numStripsAcross(); }
  double stripLength() const { return length() / numSegments(); }

  int numROCRows() const { return numStripsAcross() / numROCX(); } 
  int numROCCols() const { return numSegments() / numROCY(); }

  int totalROCs() const { return numROCX() * numROCY(); }

  std::pair<XYZVector, int> checkHitSegment(const XYZVector& trackOrig, const XYZVector& trackDir) const;

  void build() { 
    try { check(); } 
    catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }
    cleanup(); 
  }
  
  bool hasPoly() const { return poly_ != 0; }
  void clearPoly() { delete poly_; poly_ = 0; } 
  void assignPoly(Polygon3d<4>* const poly) { poly_ = poly; }
};

