#ifndef MODULETYPE_H
#define MODULETYPE_H

#include <map>

using std::map;

#include "Property.hh"
/*
class ModuleType : public PropertyObject {
public:
  Property<int, Default> numFaces;
  Property<float, Default> sensorThickness;
  Property<float, Default> dsDistance;
  Property<string, Default> innerSensorType;
  Property<string, Default> outerSensorType;
  Property<float, Default> waferDiameter;

  ModuleType() :
    numFaces("numFaces", unchecked(), 1),
    sensorThickness("sensorThickness", unchecked(), 0.1),
    dsDistance("dsDistance", unchecked(), 0.),
    innerSensorType("innerSensorType", unchecked(), string("null")),
    outerSensorType("outerSensorType", unchecked(), string("null")),
    waferDiameter("waferDiameter", unchecked(), 1.5) // CUIDADO Set suitable default!!!!
  {}
};


class ModuleTypeRepo {

  map<string, ModuleType*> types_;
  ModuleType* defType_;

public:
  ModuleTypeRepo() : defType_(new ModuleType()) {}
  static ModuleTypeRepo& getInstance() {
    static ModuleTypeRepo typeRepo;
    return typeRepo;
  }

  void store(const PropertyTree& pt) {
    for (auto& c : pt.getChildren("ModuleType")) {
      ModuleType* mt = new ModuleType();
      mt->store(c);
      mt->check();
      types_[c.getValue()] = mt;
    }
  }

  ModuleType* get(const string& typestr) const { return types_.count(typestr) > 0 ? types_.at(typestr) : NULL; }
  ModuleType* getDefault() const { return defType_; }
};
*/



#endif
