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


/** Memo: DTC CMMSW Id numbering scheme.
    In each Outer Tracker (Z) end, there are 9 Phi-sectors. In each Phi-sector, there are 12 DTC slots.
    Hence, on each Outer Tracker (Z) end, there are a total of 9*12 = 108 DTCs.

    DTC Id = phiSectorIndex * 12 + DTC slot.

    phiSectorIndex: from 0 to 8.
    phiSectorIndex == 0 <-> Phi Sector (0 deg, 40 deg). (CMS frame of reference)
    phiSectorIndex == 1 <-> Phi Sector (40 deg, 80 deg).
    phiSectorIndex == 2 <-> Phi Sector (80 deg, 120 deg).
    ...
    phiSectorIndex == 8 <-> Phi Sector (320 deg, 360 deg).

    DTC slot index: from 1 to 12.
    DTC slot index from 1 to 3 <-> PS 10G.
    DTC slot index from 4 to 6 <-> PS 5G.
    DTC slot index from 7 to 12 <-> 2S.

    On Tracker (-Z) end, 108 is added to the DTC id.
    Hence the DTC ids range from 1 to 108 (Tracker (+Z) end) and 109 to 216 (Tracker (-Z) end).

    Example:
    DTC Id = 182.
    This DTC is connected to modules on (-Z) end.
    182 - 108 = 74.
    74 = 6 * 12 +2.
    Hence phiSectorIndex = 6 and DTC slot index = 2.
    Hence this DTC is connected to modules located in the (240 deg, 280 deg) phi sector, and in PS 10G category.
*/
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
