/**
 * @file MaterialSection.h
 *
 * @date 17/Oct/2014
 * @author Stefano Martina
 */

#ifndef MATERIALSECTION_H_
#define MATERIALSECTION_H_

#include "MaterialObject.hh"

namespace material {
  class MaterialSection : public MaterialObject {
  public:
    MaterialSection(double newMinZ, double newMinR, double newMaxZ, double newMaxR, Direction newBearing, MaterialSection* nextSection);
    MaterialSection(double newMinZ, double newMinR, double newMaxZ, double newMaxR, Direction newBearing);
    nextSection(MaterialSection* nextSection);
    inactiveElement(InactiveElement* inactiveElement);
  protected:
    Property<double, noDefault> minZ;
    Property<double, noDefault> maxZ;
    Property<double, noDefault> minR;
    Property<double, noDefault> maxR;
    Property<Direction, noDefault> bearing;
    MaterialSection* nextSection_;
    InactiveElement* inactiveElement_;
  }
}

#endif /* MATERIALSECTION_H_ */

