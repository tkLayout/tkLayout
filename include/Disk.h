#ifndef DISK_H
#define DISK_H

#include <vector>
#include <string>
#include <memory>

#include <boost/ptr_container/ptr_vector.hpp>

#include "global_funcs.h"
#include "Property.h"
#include "Ring.h"

class Disk : public PropertyObject, public Buildable, public Identifiable<int> {
public:
  typedef PtrVector<Ring> Container;
  //typedef boost::ptr_map<int, Ring> RingIndexMap;
  //typedef PtrMap<int, Ring> RingIndexMap;
  typedef std::map<int, Ring*> RingIndexMap;
private:
  Container rings_;
  RingIndexMap ringIndexMap_;
  

  Property<double, NoDefault> innerRadius;
  Property<double, NoDefault> outerRadius;
  Property<double, NoDefault> bigDelta;
  Property<double, Default> rOverlap;
  Property<int, Default> bigParity;

  PropertyNode<int> ringNode;

  inline double getDsDistance(const vector<double>& buildDsDistances, int rindex) const;
  void buildTopDown(const vector<double>& buildDsDistances);
  void buildBottomUp(const vector<double>& buildDsDistances);

  double averageZ_ = 0;
public:
  Property<int, NoDefault> numRings;
  Property<double, NoDefault> zError;
  Property<double, NoDefault> buildZ;
  Property<double, NoDefault> placeZ;

  ReadonlyProperty<double, Computable> minZ, maxZ, minR, maxR;
  ReadonlyProperty<int, Computable> totalModules;
  ReadonlyProperty<double, Computable> maxRingThickness;

  Disk() :
    numRings("numRings", parsedAndChecked()),
    innerRadius("innerRadius", parsedAndChecked()),
    outerRadius("outerRadius", parsedAndChecked()),
    bigDelta("bigDelta", parsedAndChecked()),
    zError("zError", parsedAndChecked()),
    rOverlap("rOverlap", parsedOnly(), 1.),
    bigParity("bigParity", parsedOnly(), 1),
    buildZ("buildZ", parsedOnly()),
    placeZ("placeZ", parsedOnly()),
    ringNode("Ring", parsedOnly())
  {}

  void setup() {
    minZ.setup([this]() { double min = 99999; for (const Ring& r : rings_) { min = MIN(min, r.minZ()); } return min; });
    maxZ.setup([this]() { double max = 0; for (const Ring& r : rings_) { max = MAX(max, r.maxZ()); } return max; });
    minR.setup([this]() { double min = 99999; for (const Ring& r : rings_) { min = MIN(min, r.minR()); } return min; });
    maxR.setup([this]() { double max = 0; for (const Ring& r : rings_) { max = MAX(max, r.maxR()); } return max; });
    maxRingThickness.setup([this]() { double max = 0; for (const Ring& r : rings_) { max = MAX(max, r.thickness()); } return max; });
    totalModules.setup([this]() { int cnt = 0; for (const Ring& r : rings_) { cnt += r.numModules(); } return cnt; });
  }

  void check() override;
  void build(const vector<double>& buildDsDistances);
  void translateZ(double z);
  void mirrorZ();
  void cutAtEta(double eta);

  double averageZ() const { return averageZ_; }
  double thickness() const { return bigDelta()*2 + maxRingThickness(); } 

  const Container& rings() const { return rings_; }
  const RingIndexMap& ringsMap() const { return ringIndexMap_; }

  void accept(GeometryVisitor& v) { 
    v.visit(*this); 
    for (auto& r : rings_) { r.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const { 
    v.visit(*this); 
    for (const auto& r : rings_) { r.accept(v); }
  }
};

#endif
