#ifndef PHIPOSITION_HH
#define PHIPOSITION_HH

//#include <global_constants.hh>
//#include "global_funcs.hh"
//#include "Cabling/cabling_constants.hh"
#include "Cabling/cabling_functions.hh"


class PhiPosition {
public:
  PhiPosition(const double phi, const int numPhiSegments, const bool isPositiveCablingSide, const bool isBarrel, const int layerDiskNumber, const std::string subDetectorName = "", const std::string bundleType = "");

  const double phiSegmentWidth() const { return phiSegmentWidth_; }
  const double phiSegmentStart() const { return phiSegmentStart_; }
  const int phiSegmentRef() const { return phiSegmentRef_; }
  
  const double phiRegionWidth() const { return phiRegionWidth_; }
  const double phiRegionStart() const { return phiRegionStart_; }
  const int phiRegionRef() const { return phiRegionRef_; }

  const double phiSectorWidth() const { return phiSectorWidth_; }
  const double phiSectorStart() const { return phiSectorStart_; }
  const int phiSectorRef() const { return phiSectorRef_; }

private:
  double phiSegmentWidth_;
  double phiSegmentStart_;
  int phiSegmentRef_;

  double phiRegionWidth_;
  double phiRegionStart_;
  int phiRegionRef_;

  const double phiSectorWidth_ = cabling_nonantWidth;
  double phiSectorStart_;
  int phiSectorRef_;
};


#endif  // PHIPOSITION_HH
