#include "OuterCabling/cabling_functions.hh"


/* Compute the offset in Phi with respect to Phi = 0.
 */
const double computePhiSegmentStart(const double phi, const double phiSegmentWidth) {
  double phiSegmentStart = femod(phi, phiSegmentWidth);
  return phiSegmentStart;
}


/* Compute the int n for which we have: phi ~= (n * phiSegmentWidth + phiSegmentStart).
 */
const int computePhiSegmentRef(const double phi, const double phiSegmentStart, const double phiSegmentWidth) {
  int phiSegmentRef = round(femod(phi - phiSegmentStart, 2.*M_PI) / phiSegmentWidth);
  return phiSegmentRef;
}


/* Compute the int n for which we have: (phiSliceStart + n*phiSliceWidth) <= phi < (phiSliceStart + (n+1)*phiSliceWidth).
 */
const int computePhiSliceRef(const double phi, const double phiSliceStart, const double phiSliceWidth) {
  double phiSliceRefExact = femod(phi - phiSliceStart, 2.*M_PI) / phiSliceWidth;
  int phiSliceRef = 0;
  // In case phiSliceRefExact is an integer, round it to an int!
  if (fabs((phiSliceRefExact - round(phiSliceRefExact))) < cabling_roundingTolerance) phiSliceRef = fabs(round(phiSliceRefExact));
  else phiSliceRef = std::floor(phiSliceRefExact);

  return phiSliceRef;
}


/* Compute phiSliceRef + 1 modulo numPhiSlices.
 */
const int computeNextPhiSliceRef(const int phiSliceRef, const int numPhiSlices) {
  int nextPhiSliceRef = femod( (phiSliceRef + 1), numPhiSlices);
  return nextPhiSliceRef;
}


/* Compute phiSliceRef - 1 modulo numPhiSlices.
 */
const int computePreviousPhiSliceRef(const int phiSliceRef, const int numPhiSlices) {
  int previousPhiSliceRef = femod( (phiSliceRef - 1), numPhiSlices);
  return previousPhiSliceRef;
}
