/**
 * @file MaterialSection.cc
 *
 * @date 31/oct/2014
 * @author Stefano Martina
 */

#include "MaterialSection.hh"

namespace material {
  MaterialSection::MaterialSection(double newMinZ, double newMinR, double newMaxZ, double newMaxR, Direction newBearing, MaterialSection* nextSection) :
    nextSection_(nextSection) {
    minZ(newMinZ);
    maxZ(newMaxZ);
    minR(newMinR);
    maxR(newMaxR);
    bearing(newBearing);
  }
  MaterialSection::MaterialSection(double newMinZ, double newMinR, double newMaxZ, double newMaxR, Direction newBearing) :
    MaterialSection(newMinZ, newMinR, newMaxZ, newMaxR, newBearing, nullptr) {}
  MaterialSection::~MaterialSection {}

  MaterialSection::double isHit(double z, double r, double end, Direction aDirection) {
    if (aDirection==HORIZONTAL) {
      if ((minR() - Materialway::sectionWidth - Materialway::safetySpace < r)&&(maxR() + Materialway::sectionWidth + Materialway::safetySpace > r)) {
        if (minZ()>z){
          if (minZ() <= end + Materialway::safetySpace) {
            return minZ();
          }
        } else if (maxZ()>z) {
          return -1;
        }
      }
    } else {
      if ((minZ() - Materialway::sectionWidth - Materialway::safetySpace < z)&&(maxZ() + Materialway::sectionWidth + Materialway::safetySpace > z)) {
        if (minR()>r) {
          if (minR() <= end + Materialway::safetySpace) {
            return minR();
          }
        } else if (maxR()>r) {
          return -1;
        }
      }
    }
    return 0;
  }

}
