#ifndef RING_H
#define RING_H

#include <vector>
#include <string>
#include <memory>
#include <limits.h>

#include <boost/ptr_container/ptr_vector.hpp>

using std::vector;
using std::string;

#include "global_funcs.h"
#include "Property.h"
#include "Module.h"
#include "Visitable.h"

#define MAX_WEDGE_CALC_LOOPS 100


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
  std::pair<double, int> computeOptimalRingParametersRectangle(double moduleWidth, double maxRadius);

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
  double averageZ_ = 0;

public:
  enum BuildDirection { TOPDOWN, BOTTOMUP };

  ReadonlyProperty<double, NoDefault> smallDelta;
  ReadonlyProperty<double, Computable> maxModuleThickness;
  Property<BuildDirection, NoDefault> buildDirection;
  Property<int   , AutoDefault> disk;
  Property<double, NoDefault> buildStartRadius;
  Property<double, NoDefault> buildCropRadius;
  Property<double, Computable> minZ, maxZ;
  Property<int   , NoDefault> numModules; // if set forces the number of modules (in phi) to be exactly numModules

  Property<double, Default> zRotation;
  Property<double, Default> ringOuterRadius;
  Property<double, Default> ringInnerRadius;

  double minR()      const { return minRadius_; }
  double maxR()      const { return maxRadius_; }
  double thickness() const { return smallDelta()*2 + maxModuleThickness(); } 

  const Container& modules() const { return modules_; }


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
    minZ.setup([this]() { double min = INT_MAX; for (const auto& m : modules_) min = MIN(min, m.minZ()); return min; });
    maxZ.setup([this]() { double max = 0; for (const auto& m : modules_) max = MAX(max, m.maxZ()); return max; });
    maxModuleThickness.setup([this]() { 
      double max = 0; 
      for (const auto& m : modules_) { 
        max = MAX(max, m.thickness()); 
      } 
      return max; 
    });
  }
  
  void build();
  void check() override;

  void translateZ(double z);
  void mirrorZ();
  double averageZ() const { return averageZ_; }
  void cutAtEta(double eta);

  void accept(GeometryVisitor& v) { 
    v.visit(*this); 
    for (auto& m : modules_) { m.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const { 
    v.visit(*this); 
    for (const auto& m : modules_) { m.accept(v); }
  }
  const MaterialObject& materialObject() const;
};


#endif
