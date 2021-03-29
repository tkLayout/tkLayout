#include "capabilities.hh"

/**
 * This calculates a DetId !
 * A DetId is calculated thanks to all the Ids from the geometry hierarchy.
 * A geometry hierarchy level is, for example, an OT Barrel, or a Ring, as
 * definied in the scheme in the cfg files.
 */
void DetIdentifiable::buildDetId(std::map<int, uint32_t> geometryHierarchyIds,
                                 std::vector<int> geometryHierarchySizes) {
  geometryHierarchyIds_ = geometryHierarchyIds;

  int numLevels =
      geometryHierarchySizes.size(); // Number of geometry hierarchy levels, as
                                     // expected by the DetId scheme.
  int numComputedLevels =
      geometryHierarchyIds.size(); // Number of geometry hierarchy levels, as
                                   // computed by the geometry visitor.
  if (numComputedLevels == numLevels) {
    for (int level = 0; level < numLevels; level++) {
      uint32_t id = geometryHierarchyIds.at(level);
      int shift = geometryHierarchySizes.at(level);

      if (id <= pow(2, shift)) {
        myDetId_ <<= shift; // Shift by the hierarchy level's number of bits.
        myDetId_ |= id;     // Add the hierarchy level Id.
      }

      // Check that the hierarchy level's number of bits is not too small !
      else
        logWARNING(
            "buildDetId : At rank " + any2str(level) + ", id number " +
            any2str(id) +
            " is reached, while size allocated by the DetId scheme is 2^" +
            any2str(shift) + ".");
    }
  } else
    logWARNING("buildDetId : the geometry hierarchy does not match the DetId "
               "scheme specified in the cfg files. You did not specify an "
               "appropriate DetId scheme in the cfg files, or the DetId "
               "building code is not adapted to the present geometry.");
}
