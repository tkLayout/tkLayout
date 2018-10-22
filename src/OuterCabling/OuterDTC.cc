#include "OuterCabling/OuterDTC.hh"


OuterDTC::OuterDTC(const std::string name, const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) :
  name_(name),
  phiSectorWidth_(phiSectorWidth),
  phiSectorRef_(phiSectorRef),
  type_(type),
  slot_(slot),
  isPositiveCablingSide_(isPositiveCablingSide)
{
  const int oneCablingSideId = computeOneCablingSideId(phiSectorRef, type, slot);
  myCMSSWId_ = computeCMSSWId(oneCablingSideId, isPositiveCablingSide);
  plotColor_ = computePlotColor(oneCablingSideId);
};


const int OuterDTC::computeCMSSWId(const int oneCablingSideId, const bool isPositiveCablingSide) const {
  int cmsswId = oneCablingSideId;
  if (!isPositiveCablingSide) { cmsswId += outer_cabling_numNonants * outer_cabling_maxNumDTCsPerNonantPerZEnd; }
  return cmsswId;
}


const int OuterDTC::computePlotColor(const int oneCablingSideId) const {
  const int plotColor = oneCablingSideId;
  return plotColor;
}


const int OuterDTC::computeOneCablingSideId(const int phiSectorRef, const Category& type, const int slot) const {
  int oneCablingSideId = slot;
  if (type == Category::SS) oneCablingSideId += outer_cabling_maxNumDTCsPerNonantPerZEnd / 2;
  oneCablingSideId += phiSectorRef * outer_cabling_maxNumDTCsPerNonantPerZEnd;
  return oneCablingSideId;
}
