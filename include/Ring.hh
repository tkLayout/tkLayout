#ifndef RING_H
#define RING_H

#include <vector>
#include <string>
#include <memory>
#include <limits.h>

#include <boost/ptr_container/ptr_vector.hpp>

using std::vector;
using std::string;

#include "global_funcs.hh"
#include "Property.hh"
#include "Module.hh"
#include "Visitable.hh"

#define MAX_WEDGE_CALC_LOOPS 100



class TiltedRing : public PropertyObject, public Buildable, public Identifiable<int>, public Visitable {

 public:
  typedef PtrVector<BarrelModule> Container;
 private :
  Container modules_;
  //MaterialObject materialObject_;
  
  double thetaOuterUP_, thetaOuterDOWN_, thetaOuter_, tiltAngleIdealOuter_, deltaTiltIdealOuter_;
  double thetaInner_, tiltAngleIdealInner_, deltaTiltIdealInner_;

  double thetaStart_, thetaEnd_;

  double thetaStartInner_, thetaEndInner_;

  double numPhi_, phiOverlap_;

  double rStartOuter_REAL_, zStartOuter_REAL_, rEndOuter_REAL_, zEndOuter_REAL_;
  double rStartInner_REAL_, zStartInner_REAL_, rEndInner_REAL_, zEndInner_REAL_;


 public:
  Property<double, NoDefault> innerRadius;
  Property<double, NoDefault> outerRadius;
  Property<double, NoDefault> zInner;
  Property<double, NoDefault> zOuter;
  Property<double, NoDefault> tiltAngle;
  Property<double, NoDefault> theta_g;
  Property<double, NoDefault> ringZOverlap;

  const Container& modules() const { return modules_; }

 TiltedRing() :
  //materialObject_(MaterialObject::ROD),
    innerRadius           ("ringInnerRadius"       , parsedAndChecked()),
    outerRadius           ("ringOuterRadius"       , parsedAndChecked()),
    zInner                ("ringInnerZ"            , parsedOnly()),
    zOuter                ("ringOuterZ"            , parsedOnly()),
    tiltAngle             ("tiltAngle"             , parsedAndChecked()),
    theta_g               ("theta_g"               , parsedAndChecked()),
    ringZOverlap          ("ringZOverlap"          , parsedOnly())
      {}

  void build(double lastThetaEnd);
  void buildLeftRight(double lastThetaEnd);
  void check() override;

  void accept(GeometryVisitor& v) {
    v.visit(*this); 
    for (auto& m : modules_) { m.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const {
    v.visit(*this); 
    for (const auto& m : modules_) { m.accept(v); }
    }

  //const MaterialObject& materialObject() const { return materialObject_; };

  
  //double zOuter() const { return zOuter_; }
  //double zInner() const { return zInner_; }
  double thetaOuter() const { return thetaOuter_; }
  double thetaInner() const { return thetaInner_; }
  double thetaEnd() const { return thetaEnd_; }


  double tiltAngleIdealOuter() const { return tiltAngleIdealOuter_; }
  double deltaTiltIdealOuter() const { return deltaTiltIdealOuter_; }

  double tiltAngleIdealInner() const { return tiltAngleIdealInner_; }
  double deltaTiltIdealInner() const { return deltaTiltIdealInner_; }

  //double thetaStartInner() const { return thetaStartInner_; }
  //double thetaEndInner() const { return thetaEndInner_; }

  double averageR() const { return (innerRadius() + outerRadius()) / 2.; }
  double averageZ() const { return (zInner() + zOuter()) / 2.; }

  double deltaR() const { return outerRadius() - innerRadius(); }
  double gapR() const { return (outerRadius() - innerRadius()) / sin(theta_g() * M_PI / 180.); }
  double deltaZ() const { return zOuter() - zInner(); }

  void numPhi(double numPhi) { numPhi_ = numPhi; }
  double numPhi() const { return numPhi_; }
  double phiOverlap() const { return  phiOverlap_; }

  double rStartOuter_REAL() const { return rStartOuter_REAL_; }
  double zStartOuter_REAL() const { return zStartOuter_REAL_; } 
  double rEndOuter_REAL() const { return rEndOuter_REAL_; }
  double zEndOuter_REAL() const { return zEndOuter_REAL_; }
  double thetaEndOuter_REAL () const { return atan(rEndOuter_REAL_ / zEndOuter_REAL_); }

  double rStartInner_REAL() const { return rStartInner_REAL_; }
  double zStartInner_REAL() const { return zStartInner_REAL_; } 
  double rEndInner_REAL() const { return rEndInner_REAL_; }
  double zEndInner_REAL() const { return zEndInner_REAL_; }
  double thetaEndInner_REAL () const { return atan(rEndInner_REAL_ / zEndInner_REAL_); }

};









class Ring : public PropertyObject, public Buildable, public Identifiable<int>, public Visitable {
  typedef PtrVector<EndcapModule> Container;
  Container modules_;
  MaterialObject materialObject_;

  template<class T> int roundToOdd(T x) { return round((x-1)/2)*2+1; }
  double solvex(double y);
  double compute_l(double x, double y, double d);
  double compute_d(double x, double y, double l);
  double computeTentativePhiAperture(double moduleWaferDiameter, double minRadius);
  std::pair<double, int> computeOptimalRingParametersWedge(double moduleWaferDiameter, double minRadius);
  std::pair<double, int> computeOptimalRingParametersRectangle(double moduleWidth, double highRadius);

  void buildModules(EndcapModule* templ, int numMods, double smallDelta);
  void buildBottomUp();
  void buildTopDown();

  Property<ModuleShape, NoDefault> moduleShape;
  Property<double, Default> phiOverlap;
  Property<bool  , Default> requireOddModsPerSlice;
  Property<int   , Default> phiSegments;
  Property<int   , Default> additionalModules;
  Property<bool  , Default> alignEdges;
  Property<double, Default> ringGap;
  Property<int   , Default> smallParity;

  double minRadius_, maxRadius_;

public:
  enum BuildDirection { TOPDOWN, BOTTOMUP };

  ReadonlyProperty<double, NoDefault> smallDelta;
  ReadonlyProperty<double, Computable> maxModuleThickness;
  Property<BuildDirection, NoDefault> buildDirection;
  Property<int   , AutoDefault> disk;
  Property<double, NoDefault> buildStartRadius;
  Property<double, NoDefault> buildCropRadius;
  Property<double, Computable> minZ, maxZ;
  Property<double, Computable> minZwithHybrids, maxZwithHybrids, minRwithHybrids, maxRwithHybrids;
  Property<int   , NoDefault> numModules; // if set forces the number of modules (in phi) to be exactly numModules

  Property<double, Default> zRotation;
  Property<double, Default> ringOuterRadius;
  Property<double, Default> ringInnerRadius;

  double minR()      const { return minRadius_; }
  double maxR()      const { return maxRadius_; }
  double thickness() const { return smallDelta()*2 + maxModuleThickness(); }

  const Container& modules() const { return modules_; }

  int nModules() const { return modules_.size(); }
 
  Ring() :
      materialObject_(MaterialObject::ROD),
      moduleShape           ("moduleShape"           , parsedAndChecked()),
      phiOverlap            ("phiOverlap"            , parsedOnly(), 1.),
      requireOddModsPerSlice("requireOddModsPerSlice", parsedOnly(), false),
      phiSegments           ("phiSegments"           , parsedOnly(), 4),
      numModules            ("numModules"            , parsedOnly()),
      additionalModules     ("additionalModules"     , parsedOnly(), 0),
      alignEdges            ("alignEdges"            , parsedOnly(), true),
      ringGap               ("ringGap"               , parsedOnly(), 0.),
      smallParity           ("smallParity"           , parsedOnly(), 1),
      smallDelta            ("smallDelta"            , parsedAndChecked()),
      zRotation             ("zRotation"             , parsedOnly(), 0.),
      ringOuterRadius       ("ringOuterRadius"       , parsedOnly(), -1.),
      ringInnerRadius       ("ringInnerRadius"       , parsedOnly(), -1.)
  {}

  void setup() {
    minZ.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& m : modules_) min = MIN(min, m.minZ()); return min; });
    maxZ.setup([this]() { double max = 0; for (const auto& m : modules_) max = MAX(max, m.maxZ()); return max; });
    maxModuleThickness.setup([this]() { 
      double max = 0;
      for (const auto& m : modules_) { 
        max = MAX(max, m.thickness()); 
      } 
      return max;
    });


    minZwithHybrids.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& m : modules_) min = MIN(min, m.minZwithHybrids()); return min; });
    maxZwithHybrids.setup([this]() { double max = 0; for (const auto& m : modules_) max = MAX(max, m.maxZwithHybrids()); return max; });
    minRwithHybrids.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& m : modules_) min = MIN(min, m.minRwithHybrids()); return min; });
    maxRwithHybrids.setup([this]() { double max = 0; for (const auto& m : modules_) max = MAX(max, m.maxRwithHybrids()); return max; });

  }
  
  void build();
  void check() override;

  void translateZ(double z);
  void mirrorZ();
  double averageZ() const {
    double averageZ = 0;
    for (const auto& m : modules_) { averageZ = averageZ + m.center().Z(); } 
    averageZ /= numModules(); return averageZ;
  }

  void cutAtEta(double eta);

  void removeModules() { modules_.erase_if([](DetectorModule& m) { return (m.removeModule()); }); numModules(modules_.size()); }

  void accept(GeometryVisitor& v) { 
    v.visit(*this); 
    for (auto& m : modules_) { m.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const { 
    v.visit(*this); 
    for (const auto& m : modules_) { m.accept(v); }
  }
  void accept(SensorGeometryVisitor& v) { 
    v.visit(*this); 
    for (auto& m : modules_) { m.accept(v); }
  }
  const MaterialObject& materialObject() const;
};


#endif
