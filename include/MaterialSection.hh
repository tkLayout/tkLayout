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
    virtual ~MaterialSection();

    double isHit(double z, double r, double end, Direction aDirection);
    nextSection(MaterialSection* nextSection);
    MaterialSection* nextSection();
    inactiveElement(InactiveElement* inactiveElement);
    InactiveElement* inactiveElement();

    virtual void getServicesAndPass(MaterialObject& source);
  protected:
    Property<double, noDefault> minZ;
    Property<double, noDefault> maxZ;
    Property<double, noDefault> minR;
    Property<double, noDefault> maxR;
    Property<Direction, noDefault> bearing;
    MaterialSection* nextSection_;
    InactiveElement* inactiveElement_;
  }

  class MaterialStation : public MaterialSection {
  public:
    MaterialStation(double minZ, double minR, double maxZ, double maxR, Direction bearing, MaterialSection* nextSection);
    MaterialStation(double minZ, double minR, double maxZ, double maxR, Direction bearing);
    virtual ~MaterialStation();

    void getServicesAndPass(MaterialObject& source);
  }
}

#endif /* MATERIALSECTION_H_ */

