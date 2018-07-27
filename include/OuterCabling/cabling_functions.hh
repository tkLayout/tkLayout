#ifndef CABLING_FUNCTIONS_HH
#define CABLING_FUNCTIONS_HH

#include <global_constants.hh>
#include "global_funcs.hh"
#include "OuterCabling/cabling_constants.hh"


const double computePhiSegmentStart(const double phi, const double phiSegmentWidth);
const int computePhiSegmentRef(const double phi, const double phiSegmentStart, const double phiSegmentWidth);
const int computePhiSliceRef(const double phi, const double phiSliceStart, const double phiSliceWidth);

const int computeNextPhiSliceRef(const int phiSliceRef, const int numPhiSlices);
const int computePreviousPhiSliceRef(const int phiSliceRef, const int numPhiSlices);


#endif  // CABLING_FUNCTIONS_HH
