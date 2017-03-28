#ifndef BARREL_H
#define BARREL_H

#include <vector>
#include <string>
#include <memory>
#include <limits.h>

#include <boost/ptr_container/ptr_vector.hpp>

#include "global_funcs.hh"
#include "Property.hh"
#include "Layer.hh"
#include "Visitable.hh"

namespace material {
  class SupportStructure;
}

class Barrel : public PropertyObject, public Buildable, public Identifiable<string>, Clonable<Barrel>, public Visitable {

 private:
  typedef boost::ptr_vector<Layer>                      Container;
  typedef boost::ptr_vector<material::SupportStructure> SupportStructures;

  Container         layers_;
  SupportStructures supportStructures_;

  Property<double, NoDefault> innerRadius;
  Property<double, NoDefault> outerRadius;
  Property<bool  , Default>   sameRods;
  Property<double, Default>   barrelRotation;
  Property<double, Default>   supportMarginOuter;
  Property<double, Default>   supportMarginInner;
  Property<bool  , Default>   innerRadiusFixed;
  Property<bool  , Default>   outerRadiusFixed;
  
  PropertyNode<int>               layerNode;
  PropertyNodeUnique<std::string> supportNode;

 public:
  Barrel() : 
      numLayers(         "numLayers"         , parsedAndChecked()),
      innerRadius(       "innerRadius"       , parsedAndChecked()),
      outerRadius(       "outerRadius"       , parsedAndChecked()),
      innerRadiusFixed(  "innerRadiusFixed"  , parsedAndChecked(), true),
      outerRadiusFixed(  "outerRadiusFixed"  , parsedAndChecked(), true),
      sameRods(          "sameRods"          , parsedAndChecked(), false),
      barrelRotation(    "barrelRotation"    , parsedOnly(), 0.),
      supportMarginOuter("supportMarginOuter", parsedOnly(), 2.),
      supportMarginInner("supportMarginInner", parsedOnly(), 2.),
      skipServices(      "skipServices"      , parsedOnly(), false), // broken, do not use
      layerNode(         "Layer"             , parsedOnly()),
      supportNode(       "Support"           , parsedOnly())
      {}
  void setup() {
    maxZ.setup([this]() { double max = 0;                                  for (const auto& l : layers_) { max = MAX(max, l.maxZ()); } return max; });
    minZ.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& l : layers_) { min = MIN(min, l.minZ()); } return min; });
    maxR.setup([this]() { double max = 0;                                  for (const auto& l : layers_) { max = MAX(max, l.maxR()); } return max; });
    minR.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& l : layers_) { min = MIN(min, l.minR()); } return min; });

    maxZwithHybrids.setup([this]() { double max = 0;                                  for (const auto& l : layers_) { max = MAX(max, l.maxZwithHybrids()); } return max; });
    minZwithHybrids.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& l : layers_) { min = MIN(min, l.minZwithHybrids()); } return min; });
    maxRwithHybrids.setup([this]() { double max = 0;                                  for (const auto& l : layers_) { max = MAX(max, l.maxRwithHybrids()); } return max; });
    minRwithHybrids.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& l : layers_) { min = MIN(min, l.minRwithHybrids()); } return min; });
  }
  void build(); 
  void cutAtEta(double eta);
  void accept(GeometryVisitor& v) {
    v.visit(*this); 
    for (auto& l : layers_) { l.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const {
    v.visit(*this); 
    for (const auto& l : layers_) { l.accept(v); }
  }
  void accept(SensorGeometryVisitor& v) {
    v.visit(*this); 
    for (auto& l : layers_) { l.accept(v); }
  }

  const Container& layers() const        { return layers_; }
  SupportStructures& supportStructures() { return supportStructures_; }

  Property<        int   , NoDefault>  numLayers;
  ReadonlyProperty<double, Computable> maxZ, minZ;
  ReadonlyProperty<double, Computable> maxR, minR;
  Property<double, Computable> minZwithHybrids, maxZwithHybrids, minRwithHybrids, maxRwithHybrids;
  ReadonlyProperty<bool  , Default>    skipServices;
};

#endif
