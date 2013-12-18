#include <string>

#include "global_funcs.h"
#include "Polygon3d.h"
#include "Property.h"

enum class SensorType { Pixel, Strip, None };

class Sensor : public PropertyObject, public Buildable, public Identifiable<int> {
  Polygon3d<4> poly_;
public:
  ReadonlyProperty<int, AutoDefault> xElements;
  ReadonlyProperty<int, AutoDefault> yElements;
  ReadonlyProperty<int, Default> numSegments;
  ReadonlyProperty<int, Default> numStripsAcross;
  ReadonlyProperty<int, Default> numROCs;
  ReadonlyProperty<double, Default> sensorThickness;
  ReadonlyProperty<SensorType, Default> type;
  Property<double, NoDefault> length, minWidth, maxWidth;

  Sensor() : 
      numSegments("numSegments", parsedOnly(), 1),
      numStripsAcross("numStripsAcross", parsedOnly(), 1),
      numROCs("numROCs", parsedOnly(), 128),
      sensorThickness("sensorThickness", parsedOnly(), 0.1),
      type("sensorType", parsedOnly(), SensorType::None),
      length("l", checkedOnly()),
      minWidth("w", checkedOnly()),
      maxWidth("W", checkedOnly())
  {}

  Polygon3d<4>& poly() { return poly_; }

  int numChannels() const { return numStripsAcross() * numSegments(); }
  double minPitch() const { return minWidth() / (double)numStripsAcross(); }
  double maxPitch() const { return maxWidth() / (double)numStripsAcross(); }
  double pitch() const { return (maxWidth() + minWidth()) / 2. / (double)numStripsAcross(); }
  double stripLength() const { return length() / numSegments(); }


  void build() { 
    try { check(); } 
    catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }
    cleanup(); }
  
};

