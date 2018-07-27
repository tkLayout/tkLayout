#include "OuterCabling/PhiPosition.hh"


PhiPosition::PhiPosition(const double phi, const int numPhiSegments, const bool isBarrel, const int layerDiskNumber, const std::string subDetectorName, const Category& bundleType) {

  // BARREL
  if (isBarrel) {
    double rodPhi = phi;
    int numRods = numPhiSegments;

    // PHI SEGMENT
    phiSegmentWidth_ = (2.*M_PI) / numRods;
    phiSegmentStart_ = computePhiSegmentStart(rodPhi, phiSegmentWidth_);
    phiSegmentRef_ = computePhiSegmentRef(rodPhi, phiSegmentStart_, phiSegmentWidth_);

    // STEREO PHI SEGMENT
    double stereoRodPhi = femod(M_PI - rodPhi, 2.*M_PI);
    stereoPhiSegmentStart_ = computePhiSegmentStart(stereoRodPhi, phiSegmentWidth_);
    stereoPhiSegmentRef_ = computePhiSegmentRef(stereoRodPhi, stereoPhiSegmentStart_, phiSegmentWidth_);

    // PHI REGION
    // Depending on the layer number, different phiRegionWidth are assigned.
    // This is because for several layers, there can be too many modules per DTC, hence the phi width is defined smaller.	
    phiRegionWidth_ = ((layerDiskNumber == 1 || layerDiskNumber == 2 || layerDiskNumber == 4) ? outer_cabling_nonantWidth : outer_cabling_semiNonantWidth);
    phiRegionStart_ = 0.;
    phiRegionRef_ = computePhiSliceRef(rodPhi, phiRegionStart_, phiRegionWidth_);

    // PHI SECTOR
    phiSectorStart_ = 0.;
    phiSectorRef_ = computePhiSliceRef(rodPhi, phiSectorStart_, phiSectorWidth_);
  }

  // ENDCAPS
  else {
    double modPhi = phi;
    int numModulesInRing = numPhiSegments;

    // PHI SEGMENT
    phiSegmentWidth_ = (2.*M_PI) / numModulesInRing;
    phiSegmentStart_ = computePhiSegmentStart(modPhi, phiSegmentWidth_);
    phiSegmentRef_ = computePhiSegmentRef(modPhi, phiSegmentStart_, phiSegmentWidth_);

    // STEREO PHI SEGMENT
    double stereoModPhi =  femod(M_PI - modPhi, 2.*M_PI);
    stereoPhiSegmentStart_ = computePhiSegmentStart(stereoModPhi, phiSegmentWidth_);
    stereoPhiSegmentRef_ = computePhiSegmentRef(stereoModPhi, stereoPhiSegmentStart_, phiSegmentWidth_);
	
    // PHI REGION
    // Depending on the disk number and cabling type, different phiRegionWidth are assigned.
    // This is because for several cases, there can be too many modules per bundle, hence the phi width is defined smaller.
    phiRegionWidth_ = 0;	  
    phiRegionStart_ = 0.;
    // PS10GA, PS10GB
    if (bundleType == Category::PS10GA || bundleType == Category::PS10GB ) {
      phiRegionWidth_ = outer_cabling_nonantWidth;
    }
    // PS5G
    else if (bundleType == Category::PS5G ) {
      phiRegionWidth_ = outer_cabling_semiNonantWidth;
    }
    // 2S
    else if (bundleType == Category::SS ) {
      phiRegionWidth_ = outer_cabling_endcapStripStripPhiRegionWidth;
      // Use an offset to define these phiRegions 
      // (so that number of modules per phiRegion end up consistent with connection to 1 bundle only).
      if (subDetectorName == outer_cabling_tedd1) phiRegionStart_ = outer_cabling_tedd1StripStripPhiRegionStart;
      else phiRegionStart_ = outer_cabling_tedd2StripStripPhiRegionStart;
    }
    phiRegionRef_ = computePhiSliceRef(modPhi, phiRegionStart_, phiRegionWidth_);
    
    // PHI SECTOR
    phiSectorStart_ = 0.;
    phiSectorRef_ = computePhiSliceRef(modPhi, phiSectorStart_, phiSectorWidth_);
  }
}
