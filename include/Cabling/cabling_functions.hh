#ifndef CABLING_FUNCTIONS_HH
#define CABLING_FUNCTIONS_HH

#include <global_constants.hh>
#include "global_funcs.hh"
#include "Cabling/cabling_constants.hh"



const double computePhiSegmentStart(const double phi, const double phiSegmentWidth, const bool isPositiveCablingSide);
const int computePhiSegmentRef(const double phi, const double phiSegmentStart, const double phiSegmentWidth, const bool isPositiveCablingSide);
const int computePhiSliceRef(const double phi, const double phiSliceStart, const double phiSliceWidth, const bool isPositiveCablingSide);



#endif  // CABLING_FUNCTIONS_HH
