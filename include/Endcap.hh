#ifndef ENDCAP_H
#define ENDCAP_H

#include <vector>
#include <string>
#include <memory>
#include <limits.h>

#include <boost/ptr_container/ptr_vector.hpp>

#include "global_funcs.hh"
#include "Property.hh"
#include "Disk.hh"
#include "Visitable.hh"

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

  const ScanDiskInfo scanDiskPropertyTree(int diskNumber) const;
  const ScanEndcapInfo scanPropertyTree() const;

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

    maxRwithHybrids.setup([&]() { double max = 0;                                  for (const auto& d : disks_) { max = MAX(max, d.maxRwithHybrids()); } return max; });
    minRwithHybrids.setup([&]() { double min = std::numeric_limits<double>::max(); for (const auto& d : disks_) { min = MIN(min, d.minRwithHybrids()); } return min; });
    maxZwithHybrids.setup([&]() { double max = 0;                                  for (const auto& d : disks_) { if(d.maxZ() > 0 ) max = MAX(max, d.maxZwithHybrids()); } return max; });
    minZwithHybrids.setup([&]() {
	double min = std::numeric_limits<double>::max();
	for (const auto& d : disks_) { if(d.minZ() > 0 ) min = MIN(min, d.minZwithHybrids()); }
	if (innerZ.state()) { min = MIN(min, innerZ()); } // if minZ (or barrelGap) is specified in cfg file, this should be taken into account as well !
	return min;
      });
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
  void accept(SensorGeometryVisitor& v) {
    v.visit(*this);
    for (auto& d : disks_) { d.accept(v); }
  }

  const Container& disks() const         { return disks_; }
  SupportStructures& supportStructures() { return supportStructures_; }

  Property<        int   , NoDefault>  numDisks;
  Property<        double, NoDefault>  barrelMaxZ;
  Property<        double, NoDefault>  innerZ;
  Property<        double, NoDefault>  outerZ;
  ReadonlyProperty<double, Computable> maxZ, minZ;
  ReadonlyProperty<double, Computable> maxR, minR;
  Property<double, Computable> minZwithHybrids, maxZwithHybrids, minRwithHybrids, maxRwithHybrids;
  ReadonlyProperty<bool  , Default>    skipServices;
};



#endif
