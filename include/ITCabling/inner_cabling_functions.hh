#ifndef INNER_CABLING_FUNCTIONS_HH
#define INNER_CABLING_FUNCTIONS_HH

#include <global_constants.hh>
#include "global_funcs.hh"
#include "ITCabling/inner_cabling_constants.hh"
#include "MessageLogger.hh"


namespace inner_cabling_functions {

  const int computePhiUnitRef(const double phi, const int numPhiSubUnits, const bool isPositiveZEnd);
  const int computePhiSubUnitRef(const double phi, const int numPhiSubUnits, const bool isPositiveZEnd);
  const double computePhiSubUnitStart(const double phi, const double phiSubUnitWidth);
  const double computePhiSubUnitWidth(const int numPhiSubUnits);
  const double computeStereoPhi(const double phi, const bool isPositiveZEnd);

  const bool isBarrel(const std::string subDetectorName);
  const int computeInnerTrackerQuarterIndex(const bool isPositiveZEnd, const bool isPositiveXSide);
  const int computeSubDetectorIndex(const std::string subDetectorName);
  const int computeRingQuarterIndex(const int ringNumber, const bool isRingInnerEnd);
  const int computeRingNumber(const int ringQuarterIndex);
  const bool isRingInnerEnd(const int ringQuarterIndex);

}




//const int computePhiSliceRef(const double phi, const double phiSliceStart, const double phiSliceWidth);
//const int computeNextPhiSliceRef(const int phiSliceRef, const int numPhiSlices);
//const int computePreviousPhiSliceRef(const int phiSliceRef, const int numPhiSlices);


#endif  // INNER_CABLING_FUNCTIONS_HH
