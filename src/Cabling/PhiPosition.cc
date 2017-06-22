#include "Cabling/PhiPosition.hh"


PhiPosition::PhiPosition(const double phi, const int numPhiSegments, const bool isPositiveCablingSide, const bool isBarrel, const int layerDiskNumber, const std::string subDetectorName, const Category& bundleType) {

  if (isBarrel) {
    double rodPhi = phi;
    int numRods = numPhiSegments;

    phiSegmentWidth_ = (2.*M_PI) / numRods;
    phiSegmentStart_ = computePhiSegmentStart(rodPhi, phiSegmentWidth_, isPositiveCablingSide);
    phiSegmentRef_ = computePhiSegmentRef(rodPhi, phiSegmentStart_, phiSegmentWidth_, isPositiveCablingSide);
	
    phiRegionWidth_ = ((layerDiskNumber == 1 || layerDiskNumber == 2 || layerDiskNumber == 4) ? cabling_nonantWidth : cabling_semiNonantWidth);
    phiRegionStart_ = 0.;
    phiRegionRef_ = computePhiSliceRef(rodPhi, phiRegionStart_, phiRegionWidth_, isPositiveCablingSide);

    phiSectorStart_ = 0.;
    phiSectorRef_ = computePhiSliceRef(rodPhi, phiSectorStart_, phiSectorWidth_, isPositiveCablingSide);
  }

  else {
    double modPhi = phi;
    int numModulesInRing = numPhiSegments;

    phiSegmentWidth_ = (2.*M_PI) / numModulesInRing;
    phiSegmentStart_ = computePhiSegmentStart(modPhi, phiSegmentWidth_, isPositiveCablingSide);
    phiSegmentRef_ = computePhiSegmentRef(modPhi, phiSegmentStart_, phiSegmentWidth_, isPositiveCablingSide);
	
    phiRegionWidth_ = 0;	  
    phiRegionStart_ = 0.;
    if (bundleType == Category::PS10G || bundleType == Category::PS5GA ) {
      phiRegionWidth_ = cabling_nonantWidth;
    }
    else if (bundleType == Category::PS5GB ) {
      phiRegionWidth_ = cabling_semiNonantWidth;
    }
    else if (bundleType == Category::SS ) {
      phiRegionWidth_ = cabling_endcapStripStripPhiRegionWidth;
      if (subDetectorName == cabling_tedd1) phiRegionStart_ = cabling_tedd1StripStripPhiRegionStart;
      else phiRegionStart_ = cabling_tedd2StripStripPhiRegionStart;
    }
    phiRegionRef_ = computePhiSliceRef(modPhi, phiRegionStart_, phiRegionWidth_, isPositiveCablingSide);
    
    phiSectorStart_ = 0.;
    phiSectorRef_ = computePhiSliceRef(modPhi, phiSectorStart_, phiSectorWidth_, isPositiveCablingSide);
  }
}
