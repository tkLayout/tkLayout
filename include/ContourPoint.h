#ifndef CONTOURPOINT_H
#define CONTOURPOINT_H

#include <string>
#include <exception>

#include "global_funcs.h"
#include "Polygon3d.h"
#include "Property.h"
#include "CoordinateOperations.h"

class DetectorModule;

class ContourPoint : public PropertyObject, public Buildable, public Identifiable<int>, public DetIdentifiable {
 public:
  ReadonlyProperty<double, NoDefault> pointX;
  ReadonlyProperty<double, NoDefault> pointY;

 ContourPoint() :
  pointX("pointX", parsedOnly()),
    pointY("pointY", parsedOnly())
      {}

};

#endif
