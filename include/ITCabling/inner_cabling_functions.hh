#ifndef INNER_CABLING_FUNCTIONS_HH
#define INNER_CABLING_FUNCTIONS_HH

#include <global_constants.hh>
#include "global_funcs.hh"
#include "ITCabling/inner_cabling_constants.hh"
#include "MessageLogger.hh"


namespace inner_cabling_functions {

  const int computePhiUnitRef(const double phi, const int numPhiUnits, const bool isPositiveZEnd);
  const double computePhiFromMinY(const double phi, const bool isPositiveZEnd);
  const double computeStereoPhi(const double phi, const bool isPositiveZEnd);
  const double computePhiUnitWidth(const int numPhiUnits);
  const double computePhiUnitStart(const double phi, const double phiUnitWidth);  

  const bool isBarrel(const std::string subDetectorName);
  const int computeInnerTrackerQuarterIndex(const bool isPositiveZEnd, const bool isPositiveXSide);
  const int computeSubDetectorIndex(const std::string subDetectorName);
  const int computeRingQuarterIndex(const int ringNumber, const bool isRingInnerEnd);
  const int computeRingNumber(const int ringQuarterIndex);
  const bool isRingInnerEnd(const int ringQuarterIndex);

  const int computeNumELinksPerModule(const std::string subDetectorName, const int layerOrRingNumber);
}


#endif  // INNER_CABLING_FUNCTIONS_HH
