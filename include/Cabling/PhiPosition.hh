#ifndef PHIPOSITION_HH
#define PHIPOSITION_HH

#include "Cabling/cabling_functions.hh"


/* This class contains several useful Phi identifiers.
   A phi Slice is used to name a phiSegment, a phiRegion, or a phiSector.
   * phiSegment : Phi slice delimited by 2 consecutive (centers of) modules in Phi.
   * phiRegion : Phi slice of size defined in constructor. It is used for the following:
   Barrel : Phi slice containing bundles which can be connected to the same DTC.
   Endcap : Phi Slice contaning modules which can be connected to the same bundle.
   * phiSector : Phi slice of size cabling_nonantWidth.

   There is the following classification, by the phi angle which is covered : phiSegment < phiRegion < phiSector.

   NB: These phi slices correspond to purely geometrical delimitations only.
   In exceptional cases, the modules can be artificially placed in a phi slice which does not correspond to its geometrical position : staggering.
   This is used to avoid that a bundle (or cable) connects to more modules than possible.
*/
class PhiPosition {
public:
  PhiPosition(const double phi, const int numPhiSegments, const bool isBarrel, const int layerDiskNumber, const std::string subDetectorName = "", const Category& bundleType = Category::UNDEFINED);

  // PHI SEGMENT
  const double phiSegmentWidth() const { return phiSegmentWidth_; }
  const double phiSegmentStart() const { return phiSegmentStart_; }
  const int phiSegmentRef() const { return phiSegmentRef_; }

  // COMPLEMENTARY PHI SEGMENT
  const double complementaryPhiSegmentStart() const { return complementaryPhiSegmentStart_; }
  const int complementaryPhiSegmentRef() const { return complementaryPhiSegmentRef_; }
  
  // PHI REGION
  const double phiRegionWidth() const { return phiRegionWidth_; }
  const double phiRegionStart() const { return phiRegionStart_; }
  const int phiRegionRef() const { return phiRegionRef_; }

  // PHI SECTOR
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const double phiSectorStart() const { return phiSectorStart_; }
  const int phiSectorRef() const { return phiSectorRef_; }

private:
  // PHI SEGMENT
  double phiSegmentWidth_;
  double phiSegmentStart_;
  int phiSegmentRef_;

  // COMPLEMETARY PHI SEGMENT
  double complementaryPhiSegmentStart_;
  int complementaryPhiSegmentRef_;

  // PHI REGION
  double phiRegionWidth_;
  double phiRegionStart_;
  int phiRegionRef_;

  // PHI SECTOR
  const double phiSectorWidth_ = cabling_nonantWidth;
  double phiSectorStart_;
  int phiSectorRef_;
};


#endif  // PHIPOSITION_HH
