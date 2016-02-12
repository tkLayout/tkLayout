#ifndef LAYER_H
#define LAYER_H

#include <vector>
#include <string>
#include <memory>
#include <limits.h>

#include "global_funcs.h"
#include "Property.h"
#include "Module.h"
#include "RodPair.h"
#include "Ring.h"
#include "Visitable.h"
#include "MaterialObject.h"

namespace material {
  class ConversionStation;
}

using std::string;
using std::vector;
using std::pair;
using std::unique_ptr;
using material::MaterialObject;
using material::ConversionStation;

typedef std::vector<TiltedRing*> TiltedRodTemplate;

class Layer : public PropertyObject, public Buildable, public Identifiable<int>, public Clonable<Layer>, public Visitable {
public:
  typedef PtrVector<RodPair> Container;
private:
  Container rods_;
  MaterialObject materialObject_;
  ConversionStation* flangeConversionStation_;
  std::vector<ConversionStation*> secondConversionStations_;
  double flatPartThetaEnd_;
 
  double calculatePlaceRadius(int numRods, double bigDelta, double smallDelta, double dsDistance, double moduleWidth, double overlap);
  pair<float, int> calculateOptimalLayerParms(const RodTemplate&);
  RodTemplate makeRodTemplate();
  TiltedRodTemplate makeTiltedRodTemplate(double flatPartThetaEnd);

  Property<double, NoDefault> smallDelta, bigDelta;
  Property<int, Default> bigParity;
  Property<double, Default> phiOverlap;
  Property<int, Default> phiSegments;
  Property<int, NoDefault> numModulesPhi;

  PropertyNode<int> ringNode; // to grab properties for specific rod modules
  PropertyNodeUnique<std::string> stationsNode;

  double placeRadius_;
  int numRods_;

  void buildStraight(bool isFlatPart);
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
  Property<double, Default> layerRotation;

  Property<int, NoDefault> buildNumModulesFlat;
  Property<int, NoDefault> buildNumModulesTilted;
  Property<bool, Default> isTilted;
  Property<bool, Default> isTiltedAuto;
  Property<string, AutoDefault> tiltedLayerSpecFile;

  Layer() :
            materialObject_(MaterialObject::LAYER),
            flangeConversionStation_(nullptr),
            smallDelta     ("smallDelta"     , parsedAndChecked()),
            bigDelta       ("bigDelta"       , parsedAndChecked()),
            bigParity      ("bigParity"      , parsedOnly(), -1),
            phiOverlap     ("phiOverlap"     , parsedAndChecked(), 1.),
            phiSegments    ("phiSegments"    , parsedAndChecked(), 4),
	    numModulesPhi  ("numModulesPhi"  , parsedOnly()),
            ringNode       ("Ring"           , parsedOnly()),
            stationsNode   ("Station"        , parsedOnly()),
            buildNumModules("numModules"     , parsedOnly()),
            maxZ           ("maxZ"           , parsedOnly()),
            radiusMode     ("radiusMode"     , parsedAndChecked(), RadiusMode::AUTO),
            placeRadiusHint("placeRadiusHint", parsedOnly()),
            minBuildRadius ("minBuildRadius" , parsedOnly()),
            maxBuildRadius ("maxBuildRadius" , parsedOnly()),
	    layerRotation  ("layerRotation",   parsedOnly(), 0.),
	    sameParityRods ("sameParityRods" , parsedAndChecked(), false),
	    buildNumModulesFlat("numModulesFlat"     , parsedOnly()),	    
	    buildNumModulesTilted("numModulesTilted"     , parsedOnly()),
	    isTilted       ("isTilted"       , parsedOnly(), false),
	    isTiltedAuto   ("isTiltedAuto"   , parsedOnly(), true),
            tiltedLayerSpecFile("tiltedLayerSpecFile", parsedOnly())
  { setup(); }

  void setup() {
    maxZ.setup([this]() { return rods_.front().maxZ(); });
    minZ.setup([this]() { return rods_.front().minZ(); });
    maxR.setup([this]() { double max = 0; for (const auto& r : rods_) { max = MAX(max, r.maxR()); } return max; });
    minR.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& r : rods_) { min = MIN(min, r.minR()); } return min; });
  }

  double placeRadius() const { return placeRadius_; }
  int numRods() const { return rods_.size(); }
  int numModulesPerRod() const { return rods_.front().numModules(); }
  int numModulesPerRodSide(int side) const { return rods_.front().numModulesSide(side); }
  int totalModules() const { return numModulesPerRod()*numRods(); }
  double rodThickness() const { return rods_.front().thickness(); }
  double flatPartRodThickness() const { return smallDelta()*2. + rods_.front().maxModuleThickness(); }
  //bool isTilted() const { return rods_.front().isTilted(); }
  
  //double tilt() const { return 0.0; }
  double startAngle() const { return 90.0; }

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
  const MaterialObject& materialObject() const;

  ConversionStation* flangeConversionStation() const;
  const std::vector<ConversionStation*>& secondConversionStations() const;
};



#endif
