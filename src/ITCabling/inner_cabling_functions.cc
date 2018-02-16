#include "ITCabling/inner_cabling_functions.hh"


namespace inner_cabling_functions {

  //const int computePhiUnitRef(const double phi, const double phiUnitStart, const double phiUnitWidth) {
  const int computePhiUnitRef(const double phi, const int numPhiSubUnits, const bool isPositiveZEnd) {
    const int phiSubUnitRef = computePhiSubUnitRef(phi, numPhiSubUnits, isPositiveZEnd);

    // TO DO: should be simplified
    const double phiUnitRefExact = phiSubUnitRef / 2.;
    int phiUnitRef = 0;
    if (fabs(phiUnitRefExact - round(phiUnitRefExact)) < inner_cabling_roundingTolerance) phiUnitRef = fabs(round(phiUnitRefExact));
    else phiUnitRef = std::floor(phiUnitRefExact);
    return phiUnitRef;
  }


  /* Compute the int n for which we have: phi ~= (n * phiSubUnitWidth + phiSubUnitStart).
   * n = 0 for the module with the lowest Phi > -Pi/2.
   * n = 0 for the module with the lowest Phi > Pi/2.
   */
  //const int computePhiSubUnitRef(const double phi, const double phiSubUnitStart, const double phiSubUnitWidth) {
  const int computePhiSubUnitRef(const double phi, const int numPhiSubUnits, const bool isPositiveZEnd) {

    const double stereoPhi = computeStereoPhi(phi, isPositiveZEnd);

    const double phiSubUnitWidth = computePhiSubUnitWidth(numPhiSubUnits);
    const double phiSubUnitStart = computePhiSubUnitStart(stereoPhi, phiSubUnitWidth);

    const double phiRelative = femod(stereoPhi - phiSubUnitStart + M_PI / 2., M_PI);
    const int phiSubUnitRef = round( phiRelative / phiSubUnitWidth );
    return phiSubUnitRef;
  }


  /* Compute the offset in Phi with respect to PhiSubUnitWidth.
   */
  const double computePhiSubUnitStart(const double phi, const double phiSubUnitWidth) {
    double phiSubUnitStart = femod(phi, phiSubUnitWidth);
    return phiSubUnitStart;
  }


  const double computePhiSubUnitWidth(const int numPhiSubUnits) {
    const double phiSubUnitWidth = (2.*M_PI) / numPhiSubUnits;
    return phiSubUnitWidth;
  }


  const double computeStereoPhi(const double phi, const bool isPositiveZEnd) {
    const double stereoPhi = (isPositiveZEnd ? phi : femod(M_PI - phi, 2.*M_PI) );
    return stereoPhi;
  }






  const bool isBarrel(const std::string subDetectorName) {
    if (subDetectorName == inner_cabling_tbpx) return true;
    else if (subDetectorName == inner_cabling_tfpx || subDetectorName == inner_cabling_tepx) return false;
    else { 
      logERROR(any2str("Unknown subDetector name : ")
	       + any2str(subDetectorName)
	       );
      return false;
    }
  }


  const int computeInnerTrackerQuarterIndex(const bool isPositiveZEnd, const bool isPositiveXSide) {
    int innerTrackerQuarterIndex = 0;
    if (isPositiveZEnd) {
      innerTrackerQuarterIndex  = (isPositiveXSide ? 1 : 2);
    }
    else {
      innerTrackerQuarterIndex  = (isPositiveXSide ? 3 : 4);
    }

    return innerTrackerQuarterIndex;
  }


  const int computeSubDetectorIndex(const std::string subDetectorName) {
    int subDetectorIndex = 0;
    if (subDetectorName == inner_cabling_tbpx) subDetectorIndex = 1;
    else if (subDetectorName == inner_cabling_tfpx) subDetectorIndex = 2;
    else if (subDetectorName == inner_cabling_tepx) subDetectorIndex = 3;
    else {  logINFO(any2str("Unknown subDetector name : ")
		    + any2str(subDetectorName)
		    );
    }

    return subDetectorIndex;
  }


  const int computeRingQuarterIndex(const int ringNumber, const bool isRingInnerEnd) {
    const int isRingInnerEndIndex = (!isRingInnerEnd);
    const int ringQuarterIndex = (ringNumber < 1 ? 0 : (ringNumber - 1) * 2 + isRingInnerEndIndex);
    return ringQuarterIndex;
  }


  const int computeRingNumber(const int ringQuarterIndex) {
    const int ringNumber = 1 + ringQuarterIndex / 2;
    return ringNumber;
  }


  const bool isRingInnerEnd(const int ringQuarterIndex) {
    const bool isRingInnerEnd = (ringQuarterIndex % 2 ? true : false);
    return isRingInnerEnd;
  }








  /* Compute the int n for which we have: (phiSliceStart + n*phiSliceWidth) <= phi < (phiSliceStart + (n+1)*phiSliceWidth).
   */
  /*
    const int computePhiSliceRef(const double phi, const double phiSliceStart, const double phiSliceWidth) {
    double phiSliceRefExact = femod(phi - phiSliceStart, 2.*M_PI) / phiSliceWidth;
    int phiSliceRef = 0;
    // In case phiSliceRefExact is an integer, round it to an int!
    if (fabs((phiSliceRefExact - round(phiSliceRefExact))) < inner_cabling_roundingTolerance) phiSliceRef = fabs(round(phiSliceRefExact));
    else phiSliceRef = std::floor(phiSliceRefExact);

    return phiSliceRef;
    }
  */


  /* Compute phiSliceRef + 1 modulo numPhiSlices.
   */
  /*
    const int computeNextPhiSliceRef(const int phiSliceRef, const int numPhiSlices) {
    int nextPhiSliceRef = femod( (phiSliceRef + 1), numPhiSlices);
    return nextPhiSliceRef;
    }
  */


  /* Compute phiSliceRef - 1 modulo numPhiSlices.
   */
  /*
    const int computePreviousPhiSliceRef(const int phiSliceRef, const int numPhiSlices) {
    int previousPhiSliceRef = femod( (phiSliceRef - 1), numPhiSlices);
    return previousPhiSliceRef;
    }
  */


}
