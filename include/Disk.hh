#ifndef DISK_H
#define DISK_H

#include <vector>
#include <string>
#include <memory>
#include <limits.h>

#include <boost/ptr_container/ptr_vector.hpp>

#include "global_funcs.hh"
#include "Property.hh"
#include "Ring.hh"
#include "Visitable.hh"
#include "MaterialObject.hh"

namespace material {
  class ConversionStation;
}

using material::MaterialObject;
using material::ConversionStation;

typedef std::pair<std::vector<double>, std::vector<double> > ScanDiskInfo;
typedef std::pair<ScanDiskInfo, ScanDiskInfo> ScanEndcapInfo;

class Disk : public PropertyObject, public Buildable, public Identifiable<int>, public Visitable {
public:
  typedef PtrVector<Ring> Container;
  //typedef boost::ptr_map<int, Ring> RingIndexMap;
  //typedef PtrMap<int, Ring> RingIndexMap;
  typedef std::map<int, Ring*> RingIndexMap;
private:
  Container rings_;
  RingIndexMap ringIndexMap_;
  MaterialObject materialObject_;
  ConversionStation* flangeConversionStation_;
  std::vector<ConversionStation*> secondConversionStations_;
  int diskNumber_;

  Property<double, NoDefault> innerRadius;
  Property<double, NoDefault> outerRadius;
  Property<double, NoDefault> bigDelta;
  Property<double, Default>   rOverlap;
  Property<int   , Default>   bigParity;

  PropertyNode<int> ringNode;
  PropertyNodeUnique<std::string> stationsNode;

  const std::vector<double> scanSmallDeltas() const;
  const std::vector<double> scanDsDistances() const;
  inline const double getRingInfo(const vector<double>& ringsInfo, int ringNumber) const;

  std::pair<double, double> computeStringentZ(int i, int parity, const ScanEndcapInfo& extremaDisksInfo);
  double computeNextRho(int parity, double lastZ, double newZ, double lastRho);
  void buildTopDown(const ScanEndcapInfo& extremaDisksInfo);

  double averageZ_ = 0;
public:
  Property<int, NoDefault>    numRings;
  Property<double, NoDefault> zError;
  Property<double, NoDefault> zHalfLength;
  Property<double, NoDefault> buildZ;
  Property<double, NoDefault> placeZ;

  ReadonlyProperty<double, Computable> minZ, maxZ, minR, maxR;
  Property<double, Computable> minZwithHybrids, maxZwithHybrids, minRwithHybrids, maxRwithHybrids;
  ReadonlyProperty<int, Computable> totalModules;
  ReadonlyProperty<double, Computable> maxRingThickness;

  Disk() :
    materialObject_(MaterialObject::LAYER),
    flangeConversionStation_(nullptr),
    numRings(    "numRings"   , parsedAndChecked()),
    innerRadius( "innerRadius", parsedAndChecked()),
    outerRadius( "outerRadius", parsedAndChecked()),
    bigDelta(    "bigDelta"   , parsedAndChecked()),
    zError(      "zError"     , parsedAndChecked()),
    zHalfLength( "zHalfLength", parsedAndChecked()),
    rOverlap(    "rOverlap"   , parsedOnly(), 0.),
    bigParity(   "bigParity"  , parsedOnly(), 1),
    buildZ(      "buildZ"     , parsedOnly()),
    placeZ(      "placeZ"     , parsedOnly()),
    ringNode(    "Ring"       , parsedOnly()),
    stationsNode("Station"    , parsedOnly())
  {}

  void setup() {
    minZ.setup([this]() { double min = std::numeric_limits<double>::max(); for (const Ring& r : rings_) { min = MIN(min, r.minZ()); } return min; });
    maxZ.setup([this]() { double max = std::numeric_limits<double>::lowest(); for (const Ring& r : rings_) { max = MAX(max, r.maxZ()); } return max; });
    minR.setup([this]() { double min = std::numeric_limits<double>::max(); for (const Ring& r : rings_) { min = MIN(min, r.minR()); } return min; });
    maxR.setup([this]() { double max = 0; for (const Ring& r : rings_) { max = MAX(max, r.maxR()); } return max; });

    minZwithHybrids.setup([this]() { double min = std::numeric_limits<double>::max(); for (const Ring& r : rings_) { min = MIN(min, r.minZwithHybrids()); } return min; });
    maxZwithHybrids.setup([this]() { double max = 0; for (const Ring& r : rings_) { max = MAX(max, r.maxZwithHybrids()); } return max; });
    minRwithHybrids.setup([this]() { double min = std::numeric_limits<double>::max(); for (const Ring& r : rings_) { min = MIN(min, r.minRwithHybrids()); } return min; });
    maxRwithHybrids.setup([this]() { double max = 0; for (const Ring& r : rings_) { max = MAX(max, r.maxRwithHybrids()); } return max; });

    maxRingThickness.setup([this]() { double max = 0; for (const Ring& r : rings_) { max = MAX(max, r.thickness()); } return max; });
    totalModules.setup([this]() { int cnt = 0; for (const Ring& r : rings_) { cnt += r.numModules(); } return cnt; });
  }
  const std::pair<std::vector<double>, std::vector<double> > scanPropertyTree() const;

  void check() override;
  void build(const ScanEndcapInfo& extremaDisksInfo);
  void translateZ(double z);
  void mirrorZ();
  void cutAtEta(double eta);

  double averageZ() const { return averageZ_; }
  bool side() const { return averageZ_ > 0.; }
  double thickness() const { return bigDelta()*2 + maxRingThickness(); } 

  const Container& rings() const { return rings_; }
  const RingIndexMap& ringsMap() const { return ringIndexMap_; }

  void diskNumber(int num) { diskNumber_ = num; }
  int diskNumber() const { return diskNumber_; }
  int numEmptyRings() const { return count_if(rings_.begin(), rings_.end(), [](const Ring& r) { return r.numModules() == 0; }); }

  void accept(GeometryVisitor& v) {
    v.visit(*this); 
    for (auto& r : rings_) { r.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const { 
    v.visit(*this); 
    for (const auto& r : rings_) { r.accept(v); }
  }
  void accept(SensorGeometryVisitor& v) {
    v.visit(*this); 
    for (auto& r : rings_) { r.accept(v); }
  }
  const MaterialObject& materialObject() const;
  ConversionStation* flangeConversionStation() const;
  const std::vector<ConversionStation*>& secondConversionStations() const;
  std::vector< std::set<const Module*> > getModuleSurfaces() const;
 };

#endif
