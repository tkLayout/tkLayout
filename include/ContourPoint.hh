#ifndef CONTOURPOINT_H
#define CONTOURPOINT_H

#include <exception>
#include <string>

#include "CoordinateOperations.hh"
#include "Polygon3d.hh"
#include "Property.hh"
#include "global_funcs.hh"

class DetectorModule;

class ContourPoint : public PropertyObject,
                     public Buildable,
                     public Identifiable<int>,
                     public DetIdentifiable {
public:
  ReadonlyProperty<double, NoDefault> pointX;
  ReadonlyProperty<double, NoDefault> pointY;

  ContourPoint()
      : pointX("pointX", parsedOnly()), pointY("pointY", parsedOnly()) {}
};

#endif
