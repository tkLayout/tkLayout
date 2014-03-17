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
  ReadonlyProperty<double, Computable> minR, maxR; // CUIDADO min/maxR don't take into account the sensor thickness!
  ReadonlyProperty<double, Computable> minZ, maxZ; // ditto for min/maxZ

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

  void setup() {
    minR.setup([&]() {
      XYZVector side[2];
      std::partial_sort_copy(poly_->begin(), poly_->end(), std::begin(side), std::end(side), [](const XYZVector& v1, const XYZVector& v2) { return v1.Rho() < v2.Rho(); });
      return std::min(minRVertex(), ((side[0]+side[1])/2).Rho());          
    });
    maxR.setup([&]() {
      double max = 0.;
      for (auto v : *poly_) max = MAX(max, v.Rho());
      return max;
    });
    minZ.setup([&]() {
      XYZVector side[2];
      std::partial_sort_copy(poly_->begin(), poly_->end(), std::begin(side), std::end(side), [](const XYZVector& v1, const XYZVector& v2) { return v1.Z() < v2.Z(); });
      return std::min(minZVertex(), ((side[0]+side[1])/2).Rho());
    });
    maxZ.setup([&]() {
      double max = 0.;
      for (auto v : *poly_) max = MAX(max, v.Z());
      return max;
    });
  }

  double minRVertex() const { double min = std::numeric_limits<double>::max(); for (auto v : *poly_) { min = MIN(min, v.Rho()); } return min; }
  double minZVertex() const { double min = std::numeric_limits<double>::max(); for (auto v : *poly_) { min = MIN(min, v.Z()); } return min; }
  
  bool hasPoly() const { return poly_ != 0; }
  void clearPoly() { delete poly_; poly_ = 0; } 
  void assignPoly(Polygon3d<4>* const poly) { poly_ = poly; }
};

