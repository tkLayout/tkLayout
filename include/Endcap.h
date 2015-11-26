#ifndef ENDCAP_H
#define ENDCAP_H

#include <vector>
#include <string>
#include <memory>
#include <limits.h>

#include <boost/ptr_container/ptr_vector.hpp>

#include "global_funcs.h"
#include "Property.h"
#include "Disk.h"
#include "Visitable.h"

namespace material {
  class SupportStructure;
}

class Endcap : public PropertyObject, public Buildable, public Identifiable<std::string>, public Visitable {

 private:
  typedef boost::ptr_vector<Disk>                       Container;
  typedef boost::ptr_vector<material::SupportStructure> SupportStructures;

  Container         disks_;
  SupportStructures supportStructures_;

  Property<double, NoDefault>     barrelGap;
  PropertyNode<int>               diskNode;
  PropertyNodeUnique<std::string> supportNode;

  vector<double> findMaxDsDistances();

 public:
  Endcap() :
      barrelGap(   "barrelGap"   , parsedOnly()),
      numDisks(    "numDisks"    , parsedAndChecked()),
      innerZ(      "minZ"        , parsedOnly()),
      outerZ(      "maxZ"        , parsedAndChecked()),
      skipServices("skipServices", parsedOnly(), false), // broken, do not use
      diskNode(    "Disk"        , parsedOnly()),
      supportNode( "Support"     , parsedOnly())
  {}
  void setup() {
    maxR.setup([&]() { double max = 0;                                  for (const auto& d : disks_) { max = MAX(max, d.maxR()); } return max; });
    minR.setup([&]() { double min = std::numeric_limits<double>::max(); for (const auto& d : disks_) { min = MIN(min, d.minR()); } return min; });
    maxZ.setup([&]() { double max = 0;                                  for (const auto& d : disks_) { if(d.maxZ() > 0 ) max = MAX(max, d.maxZ()); } return max; });
    minZ.setup([&]() { double min = std::numeric_limits<double>::max(); for (const auto& d : disks_) { if(d.minZ() > 0 ) min = MIN(min, d.minZ()); } return min; });
  }
  void build();
  void cutAtEta(double eta);
  void accept(GeometryVisitor& v) {
    v.visit(*this);
    for (auto& d : disks_) { d.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const {
    v.visit(*this);
    for (const auto& d : disks_) { d.accept(v); }
  }

  const Container& disks() const         { return disks_; }
  SupportStructures& supportStructures() { return supportStructures_; }

  Property<        int   , NoDefault>  numDisks;
  Property<        double, NoDefault>  barrelMaxZ;
  Property<        double, NoDefault>  innerZ;
  Property<        double, NoDefault>  outerZ;
  ReadonlyProperty<double, Computable> maxZ, minZ;
  ReadonlyProperty<double, Computable> maxR, minR;
  ReadonlyProperty<bool  , Default>    skipServices;
};



#endif
