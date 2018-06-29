#ifndef LAYER_H
#define LAYER_H

#include <vector>
#include <string>
#include <memory>
#include <limits.h>

#include "global_funcs.hh"
#include "Property.hh"
#include "Module.hh"
#include "RodPair.hh"
#include "Ring.hh"
#include "Visitable.hh"
#include "MaterialObject.hh"

namespace material {
  class ConversionStation;
}

using std::string;
using std::vector;
using std::pair;
using std::unique_ptr;
using material::MaterialObject;
using material::ConversionStation;

typedef std::map<int, TiltedRing*> TiltedRingsTemplate;

struct SkewedLayerPhiShifts {
  const double installationMinusBigDeltaRodCenterPhiShift;
  const double commonRodCenterPhiShift;
  const double skewedRodCenterPhiShift;
};

struct SkewedLayerInfo {  
  const double skewedModuleMinRho;
  const double skewedModuleCenterRho;
  const double skewedModuleMaxRho;
  const double skewAngle;

  const double unitPhiOverlapLength;
  const double installationHorizontalOverlapLength;

  const SkewedLayerPhiShifts phiShifts;
};


class FlatRingsGeometryInfo {
public:
  void calculateFlatRingsGeometryInfo(std::vector<StraightRodPair*> flatPartRods, double bigParity);
  std::map<int, double> zErrorInner() const { return zErrorInner_; }
  std::map<int, double> zErrorOuter() const { return zErrorOuter_; }
private:
  std::map<int, double> zErrorInner_;
  std::map<int, double> zErrorOuter_;
};

class TiltedRingsGeometryInfo {
public:
  TiltedRingsGeometryInfo() {};
  TiltedRingsGeometryInfo(int numModulesFlat, double, double, double, double, TiltedRingsTemplate tiltedRingsGeometry);
  std::map<int, double> deltaZInner() const { return deltaZInner_; }
  std::map<int, double> deltaZOuter() const { return deltaZOuter_; }
  std::map<int, double> zErrorInner() const { return zErrorInner_; }
  std::map<int, double> zErrorOuter() const { return zErrorOuter_; }
private:
  std::map<int, double> deltaZInner_;
  std::map<int, double> deltaZOuter_;
  std::map<int, double> zErrorInner_;
  std::map<int, double> zErrorOuter_;
};


class Layer : public PropertyObject, public Buildable, public Identifiable<int>, public Clonable<Layer>, public Visitable {
public:
  typedef PtrVector<RodPair> Container;

  Property<double, NoDefault> smallDelta, bigDelta;
  Property<int, Default> bigParity;

  Property<int, AutoDefault> buildNumModules;
  ReadonlyProperty<double, UncachedComputable> maxZ, minZ;
  ReadonlyProperty<double, Computable> maxR, minR;
  Property<double, Computable> minZwithHybrids, maxZwithHybrids, minRwithHybrids, maxRwithHybrids;

  enum RadiusMode { SHRINK, ENLARGE, FIXED, AUTO };
  Property<RadiusMode, Default> radiusMode;
  Property<double, NoDefault> placeRadiusHint;

  Property<double, NoDefault> minBuildRadius;
  Property<double, NoDefault> maxBuildRadius;
  Property<bool, Default> sameParityRods;
  Property<double, Default> layerRotation;

  Property<int, NoDefault> numRods;
  FixedSizeMultiProperty<std::vector<double>, 4,','> phiForbiddenRanges;

  Property<int, NoDefault> buildNumModulesFlat;
  Property<int, NoDefault> buildNumModulesTilted;
  Property<bool, Default> isTilted;
  Property<bool, NoDefault> isTiltedAuto;
  Property<string, AutoDefault> tiltedLayerSpecFile;

  Property<bool, Default> isSkewedForInstallation;
  Property<double, NoDefault> skewedModuleEdgeShift;
  Property<double, Default> installationOverlapRatio;
  
  Property<double, AutoDefault> skewAngle;
  Property<double, AutoDefault> skewedModuleMinRho; 
  Property<double, AutoDefault> skewedModuleCenterRho;
  Property<double, AutoDefault> skewedModuleMaxRho; 
 
  Property<double, AutoDefault> unitPhiOverlapLength;
  Property<double, AutoDefault> installationHorizontalOverlapLength; 

  Layer() :
            materialObject_(MaterialObject::LAYER),
            flangeConversionStation_(nullptr),
            smallDelta     ("smallDelta"     , parsedAndChecked()),
            bigDelta       ("bigDelta"       , parsedAndChecked()),
            bigParity      ("bigParity"      , parsedOnly(), -1),
	    phiOverlap     ("phiOverlap"     , parsedOnly()), // used to be parsedAndChecked()
	    phiSegments    ("phiSegments"    , parsedOnly()), // used to be parsedAndChecked(), and default value = 4
	    numRods        ("numRods"        , parsedOnly()),
	    phiForbiddenRanges("phiForbiddenRanges", parsedOnly()),
            ringNode       ("Ring"           , parsedOnly()),
            stationsNode   ("Station"        , parsedOnly()),
            buildNumModules("numModules"     , parsedOnly()),
            maxZ           ("maxZ"           , parsedOnly()),
            radiusMode     ("radiusMode"     , parsedAndChecked(), RadiusMode::AUTO),
            placeRadiusHint("placeRadiusHint", parsedOnly()),
            minBuildRadius ("minBuildRadius" , parsedOnly()),
            maxBuildRadius ("maxBuildRadius" , parsedOnly()),
	    layerRotation  ("layerRotation",   parsedOnly(), 0.),
	    sameParityRods ("sameParityRods" , parsedAndChecked(), true),
	    buildNumModulesFlat("numModulesFlat"     , parsedOnly()),
	    buildNumModulesTilted("numModulesTilted"     , parsedOnly()),
	    isTilted       ("isTilted"       , parsedOnly(), false),
	    isTiltedAuto   ("isTiltedAuto"   , parsedOnly()),
            tiltedLayerSpecFile("tiltedLayerSpecFile", parsedOnly()),
	    isSkewedForInstallation("isSkewedForInstallation", parsedOnly(), false),
	    skewedModuleEdgeShift("skewedModuleEdgeShift", parsedOnly()),
	    installationOverlapRatio("installationOverlapRatio", parsedOnly(), 2.) // remove default??
  { setup(); }

  void setup() {
    maxZ.setup([this]() { double max = 0; for (const auto& r : rods_) { max = MAX(max, r.maxZ()); } return max; });
    minZ.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& r : rods_) { min = MIN(min, r.minZ()); } return min; });
    maxR.setup([this]() { double max = 0; for (const auto& r : rods_) { max = MAX(max, r.maxR()); } return max; });
    minR.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& r : rods_) { min = MIN(min, r.minR()); } return min; });

    maxZwithHybrids.setup([this]() { double max = 0; for (const auto& r : rods_) { max = MAX(max, r.maxZwithHybrids()); } return max; });
    minZwithHybrids.setup([this]()  { double min = std::numeric_limits<double>::max(); for (const auto& r : rods_) { min = MIN(min, r.minZwithHybrids()); } return min; });
    maxRwithHybrids.setup([this]() { double max = 0; for (const auto& r : rods_) { max = MAX(max, r.maxRwithHybrids()); } return max; });
    minRwithHybrids.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& r : rods_) { min = MIN(min, r.minRwithHybrids()); } return min; });
  }

  void check() override;
  void build();

  void cutAtEta(double eta);
  void rotateZ(double angle) { for (auto& r : rods_) r.rotateZ(angle); }

  double placeRadius() const { return placeRadius_; }
  int numModulesPerRod() const { return rods_.front().numModules(); }
  int numModulesPerRodSide(int side) const { return rods_.front().numModulesSide(side); }
  int totalModules() const { return numModulesPerRod()*numRods(); }
  double rodThickness() const { return rods_.front().thickness(); }
  double flatPartRodThickness() const { return smallDelta()*2. + rods_.front().maxModuleThickness(); }

  const Container& rods() const { return rods_; }
  std::vector<StraightRodPair*> flatPartRods() const { return flatPartRods_; }
  double flatPartPhiOverlapSmallDeltaMinus() const { return flatPartPhiOverlapSmallDeltaMinus_; }
  double flatPartPhiOverlapSmallDeltaPlus() const { return flatPartPhiOverlapSmallDeltaPlus_; }
  double flatPartAverageR() const { return flatPartAverageR_; }
  FlatRingsGeometryInfo flatRingsGeometryInfo() const { return flatRingsGeometryInfo_; }
  TiltedRingsTemplate tiltedRingsGeometry() const { return tiltedRingsGeometry_; }
  TiltedRingsGeometryInfo tiltedRingsGeometryInfo() const { return tiltedRingsGeometryInfo_; }

  void layerNumber(int num) { layerNumber_ = num; }
  int layerNumber() const { return layerNumber_; }

  bool isTiming() const { return rods_.front().isTiming(); }

  void accept(GeometryVisitor& v) {
    v.visit(*this);
    for (auto& r : rods_) { r.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const { 
    v.visit(*this);
    for (const auto& r : rods_) { r.accept(v); }
  }
  void accept(SensorGeometryVisitor& v) {
    v.visit(*this);
    for (auto& r : rods_) { r.accept(v); }
  }
  const MaterialObject& materialObject() const { return materialObject_; }

  ConversionStation* flangeConversionStation() const { return flangeConversionStation_; }
  const std::vector<ConversionStation*>& secondConversionStations() const { return secondConversionStations_; }


private:
  void buildStraight();
  void buildTilted();

  RodTemplate makeRodTemplate(const double skewAngle = 0.);
  TiltedRingsTemplate makeTiltedRingsTemplate(double flatPartThetaEnd);

  double calculatePlaceRadius(int numRods, double bigDelta, double smallDelta, double dsDistance, double moduleWidth, double overlap);
  pair<float, int> calculateOptimalLayerParms(const RodTemplate&);
  
  const SkewedLayerPhiShifts buildSkewed();
  static const SkewedLayerInfo computeSkewedLayerInfo(const double layerCenterRho, const double bigDelta, const int numRods, const double moduleWidth, const double skewedModuleEdgeShift, const double installationOverlapRatio);

  Container rods_;
  MaterialObject materialObject_;
  ConversionStation* flangeConversionStation_;
  std::vector<ConversionStation*> secondConversionStations_;
  std::vector<StraightRodPair*> flatPartRods_;
  double flatPartPhiOverlapSmallDeltaMinus_;
  double flatPartPhiOverlapSmallDeltaPlus_;
  double flatPartAverageR_;
  FlatRingsGeometryInfo flatRingsGeometryInfo_;
  TiltedRingsTemplate tiltedRingsGeometry_;
  TiltedRingsGeometryInfo tiltedRingsGeometryInfo_;
  int layerNumber_;

  Property<double, NoDefault> phiOverlap;
  Property<int, NoDefault> phiSegments;

  PropertyNode<int> ringNode; // to grab properties for specific rod modules
  PropertyNodeUnique<std::string> stationsNode;

  double placeRadius_;
};



/*class TiltedLayer : public Layer, public PropertyObject, public Buildable, public Identifiable<int>, public Clonable<TiltedLayer>, public Visitable {

 private :
  TiltedRodTemplate tiltedRings_;

  }*/

#endif
