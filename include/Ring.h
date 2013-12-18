#ifndef RING_H
#define RING_H

#include <vector>
#include <string>
#include <memory>

#include <boost/ptr_container/ptr_vector.hpp>

using std::vector;
using std::string;

#include "global_funcs.h"
#include "Property.h"
#include "Module.h"

#define MAX_WEDGE_CALC_LOOPS 100

class Ring : public PropertyObject, public Buildable, public Identifiable<int> {

  boost::ptr_vector<EndcapModule> modules_;

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
  Property<double, Default> moduleOverlapPhi;
  Property<bool, Default> requireOddModsPerSlice;
  Property<int, Default> phiSegments;
  Property<int, Default> additionalModules;
  Property<bool, Default> alignEdges;
  Property<double, Default> ringGap;
  Property<int, Default> smallParity;

  double minRadius_, maxRadius_;

public:
  enum BuildDirection { TOPDOWN, BOTTOMUP };
  Property<int, AutoDefault> disk;
  ReadonlyProperty<double, NoDefault> smallDelta;
  Property<BuildDirection, NoDefault> buildDirection;
  Property<double, NoDefault> buildStartRadius;
  Property<double, NoDefault> buildCropRadius;
  Property<double, Computable> minZ, maxZ;

  double minR() const { return minRadius_; }
  double maxR() const { return maxRadius_; }
  int numModules() const { return modules_.size(); }

  Ring() :
      moduleShape           ("moduleShape"           , parsedAndChecked()),
      moduleOverlapPhi      ("moduleOverlapPhi"      , parsedOnly(), 1.),
      requireOddModsPerSlice("requireOddModsPerSlice", parsedOnly(), false),
      phiSegments           ("phiSegments"           , parsedOnly(), 4),
      additionalModules     ("additionalModules"     , parsedOnly(), 0),
      alignEdges            ("alignEdges"            , parsedOnly(), true),
      ringGap               ("ringGap"               , parsedOnly(), 0.),
      smallParity           ("smallParity"           , parsedOnly(), 1),
      smallDelta            ("smallDelta"            , parsedAndChecked())
  {}

  void setup() {
    minZ.setup([this]() { double min = 99999; for (const auto& m : modules_) min = MIN(min, m.minZ()); return min; });
    maxZ.setup([this]() { double max = 0; for (const auto& m : modules_) max = MAX(max, m.maxZ()); return max; });
    for (auto& m : modules_) m.setup();
  }
  
  void build();
  void check() override;

  void translateZ(double z);
  void mirrorZ();
  void cutAtEta(double eta);

  void accept(GeometryVisitor& v) { 
    v.visit(*this); 
    for (auto& m : modules_) { m.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const { 
    v.visit(*this); 
    for (const auto& m : modules_) { m.accept(v); }
  }
};


#endif
