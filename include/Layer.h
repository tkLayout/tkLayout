#ifndef LAYER_H
#define LAYER_H

#include <vector>
#include <string>
#include <memory>

#include <boost/ptr_container/ptr_vector.hpp>

#include "global_funcs.h"
#include "Property.h"
#include "Module.h"
#include "RodPair.h"

using std::string;
using std::vector;
using std::pair;
using std::unique_ptr;


class Layer : public PropertyObject, public Buildable, public Identifiable<int> {
public:
  typedef boost::ptr_vector<RodPair> Container;
private:
  Container rods_;

  double calculatePlaceRadius(int numRods, double bigDelta, double smallDelta, double dsDistance, double moduleWidth, double overlap);
  pair<float, int> calculateOptimalLayerParms(const RodTemplate&);
  RodTemplate makeRodTemplate();

  Property<double, NoDefault> smallDelta, bigDelta;
  Property<int, Default> bigParity;
  Property<double, Default> rodOverlapPhi;
  Property<int, Default> phiSegments;

  PropertyNode<int> ringNode; // to grab properties for specific rod modules

  double placeRadius_;
  int numRods_;

  void buildStraight();
  void buildTilted();
public:
  Property<int, AutoDefault> buildNumModules;
  ReadonlyProperty<double, UncachedComputable> maxZ, minZ;
  ReadonlyProperty<double, Computable> maxR, minR;

  enum RadiusMode { SHRINK, ENLARGE, FIXED, AUTO };
  Property<RadiusMode, Default> radiusMode;
  Property<double, NoDefault> placeRadiusHint;

  Property<double, NoDefault> minBuildRadius;
  Property<double, NoDefault> maxBuildRadius;
  Property<bool, Default> sameParityRods;

  Property<string, AutoDefault> tiltedLayerSpecFile;

  Layer() :
            smallDelta     ("smallDelta"     , parsedAndChecked()),
            bigDelta       ("bigDelta"       , parsedAndChecked()),
            bigParity      ("bigParity"      , parsedOnly(), -1),
            rodOverlapPhi  ("rodOverlapPhi"  , parsedAndChecked(), 1.),
            phiSegments    ("phiSegments"    , parsedAndChecked(), 4),
            ringNode       ("Ring"           , parsedOnly()),
            buildNumModules("numModules"     , parsedOnly()),
            maxZ           ("maxZ"           , parsedOnly()),
            radiusMode     ("radiusMode"     , parsedAndChecked(), RadiusMode::AUTO),
            placeRadiusHint("placeRadiusHint", parsedOnly()),
            minBuildRadius ("minBuildRadius" , parsedOnly()),
            maxBuildRadius ("maxBuildRadius" , parsedOnly()),
            sameParityRods ("sameParityRods" , parsedAndChecked(), false),
            tiltedLayerSpecFile("tiltedLayerSpecFile", parsedOnly())
  {}

  void setup() {
    maxZ.setup([this]() { return rods_.front().maxZ(); });
    minZ.setup([this]() { return rods_.front().minZ(); });
    maxR.setup([this]() { double max = 0; for (const auto& r : rods_) { max = MAX(max, r.maxR()); } return max; });
    minR.setup([this]() { double min = 99999; for (const auto& r : rods_) { min = MIN(min, r.minR()); } return min; });
    for (auto& r : rods_) r.setup();
  }

  double placeRadius() const { return placeRadius_; }
  int numRods() const { return rods_.size(); }
  int numModulesPerRod() const { return rods_.front().numModules(); };
  int totalModules() const { return numModulesPerRod()*numRods(); };
  
  double tilt() const { return 0.0; }
  double startAngle() const { return 0.0; }

  void check() override;
  void build();

  const Container& rods() const { return rods_; }

  void cutAtEta(double eta);
  void rotateZ(double angle) { for (auto& r : rods_) r.rotateZ(angle); } 

  void accept(GeometryVisitor& v) { 
    v.visit(*this); 
    for (auto& r : rods_) { r.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const { 
    v.visit(*this); 
    for (const auto& r : rods_) { r.accept(v); }
  }
};



#endif
