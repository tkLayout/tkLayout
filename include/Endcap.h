#ifndef ENDCAP_H
#define ENDCAP_H

#include <vector>
#include <string>
#include <memory>

#include <boost/ptr_container/ptr_vector.hpp>

#include "global_funcs.h"
#include "Property.h"
#include "Disk.h"


class Endcap : public PropertyObject, public Buildable, public Identifiable<std::string> {
  typedef boost::ptr_vector<Disk> Container;

  Container disks_;
  Property<double, NoDefault> barrelGap;
  PropertyNode<int> diskNode;
  
  vector<double> findMaxDsDistances();
public:
  Property<int, NoDefault> numDisks;
  Property<double, NoDefault> barrelMaxZ;
  Property<double, NoDefault> minZ;
  Property<double, NoDefault> maxZ;
  ReadonlyProperty<double, Computable> maxR, minR;
  ReadonlyProperty<bool, Default> skipServices;

  Endcap() :
      barrelGap("barrelGap", parsedOnly()),
      numDisks("numDisks", parsedAndChecked()),
      minZ("minZ", parsedOnly()),
      maxZ("maxZ", parsedAndChecked()),
      skipServices("skipServices", parsedOnly(), false), // broken, do not use
      diskNode("Disk", parsedOnly())
  {}

  void setup() {
    maxR.setup([&]() { double max = 0; for (const auto& d : disks_) { max = MAX(max, d.maxR()); } return max; });
    minR.setup([&]() { double min = 0; for (const auto& d : disks_) { min = MIN(min, d.minR()); } return min; });
    for (auto& d : disks_) d.setup();
  }

  void build();
  void cutAtEta(double eta);

  const Container& disks() const { return disks_; }

  void accept(GeometryVisitor& v) { 
    v.visit(*this); 
    for (auto& d : disks_) { d.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const { 
    v.visit(*this); 
    for (const auto& d : disks_) { d.accept(v); }
  }
};



#endif
