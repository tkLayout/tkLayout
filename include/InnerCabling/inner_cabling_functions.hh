#ifndef INNER_CABLING_FUNCTIONS_HH
#define INNER_CABLING_FUNCTIONS_HH

#include <global_constants.hh>
#include "global_funcs.hh"
#include "InnerCabling/inner_cabling_constants.hh"
#include "MessageLogger.hh"

#include <cstddef>

namespace inner_cabling_functions {

  // Phi-related free functions (IT cabling map related).
  const int computePhiUnitRef(const double phi, const int numPhiUnits, const bool isPositiveZEnd);
  const double computePhiFromMinY(const double phi, const bool isPositiveZEnd);
  const double computeStereoPhi(const double phi, const bool isPositiveZEnd);
  const double computePhiUnitWidth(const int numPhiUnits);
  const double computePhiUnitStart(const double phi, const double phiUnitWidth);  

  // Free functions used to locate a module (IT cabling map related).
  const bool isBarrel(const std::string subDetectorName);
  const int computeInnerTrackerQuarterIndex(const bool isPositiveZEnd, const bool isPositiveXSide);
  const int computeSubDetectorIndex(const std::string subDetectorName);
  const int computeHalfRingIndex(const int ringNumber, const bool isSmallerAbsZHalfRing);
  const int computeRingNumber(const int halfRingIndex);
  const bool isSmallerAbsZHalfRing(const int halfRingIndex);

  // compute number of ELinks per module
  std::size_t computeNumELinksModule(const std::string& subDetectorName, std::size_t layerOrRingNumber, int phiRefInPowerChain);
}


#endif  // INNER_CABLING_FUNCTIONS_HH
