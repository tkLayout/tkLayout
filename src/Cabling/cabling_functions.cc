#include "Cabling/cabling_functions.hh"


const double computePhiSegmentStart(const double phi, const double phiSegmentWidth, const bool isPositiveCablingSide) {
  double stereoPhi = (isPositiveCablingSide ? phi : M_PI - phi);
  double phiSegmentStart = femod(stereoPhi, phiSegmentWidth);
  return phiSegmentStart;
}


const int computePhiSegmentRef(const double phi, const double phiSegmentStart, const double phiSegmentWidth, const bool isPositiveCablingSide) {
  double stereoPhi = (isPositiveCablingSide ? phi : M_PI - phi);
  int phiSegmentRef = round(femod(stereoPhi - phiSegmentStart, 2.*M_PI) / phiSegmentWidth);
  return phiSegmentRef;
}


const int computePhiSliceRef(const double phi, const double phiSliceStart, const double phiSliceWidth, const bool isPositiveCablingSide) {
  double stereoPhi = (isPositiveCablingSide ? phi : M_PI - phi);
  double phiSliceRefExact = femod(stereoPhi - phiSliceStart, 2.*M_PI) / phiSliceWidth;
  int phiSliceRef = 0;
  if (fabs((phiSliceRefExact - round(phiSliceRefExact))) < cabling_roundingTolerance) phiSliceRef = fabs(round(phiSliceRefExact));
  else phiSliceRef = std::floor(phiSliceRefExact);

  return phiSliceRef;
}


const int computeNextPhiSliceRef(const int phiSliceRef, const int numPhiSlices) {
  int nextPhiSliceRef = femod( (phiSliceRef + 1), numPhiSlices);
  return nextPhiSliceRef;
}


const int computePreviousPhiSliceRef(const int phiSliceRef, const int numPhiSlices) {
  int previousPhiSliceRef = femod( (phiSliceRef - 1), numPhiSlices);
  return previousPhiSliceRef;
}
