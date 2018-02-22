#include "ITCabling/InnerDTC.hh"


InnerDTC::InnerDTC(const int DTCId, const bool isPositiveZEnd, const bool isPositiveXSide) :
  isPositiveZEnd_(isPositiveZEnd),
  isPositiveXSide_(isPositiveXSide)
{
  myid(DTCId);

  //DTCType_ = computePowerChainType(isBarrel_, layerDiskNumber, ringNumber_);

  //plotColor_ = computePlotColor(isBarrel_, isPositiveZEnd, phiRef, ringQuarterIndex);
};



void InnerDTC::addBundle(InnerBundle* bundle) { 
  bundles_.push_back(bundle);
}



/*const int InnerDTC::computePlotColor(const bool isBarrel, const bool isPositiveZEnd, const int phiRef, const int ringQuarterIndex) const {
  int plotColor = 0;

  const int plotPhi = femod(phiRef, 2);

  if (isBarrel) {
    const int plotZEnd = (isPositiveZEnd ? 0 : 1);
    plotColor = plotZEnd * 2 + plotPhi + 6;
  }
  else {
    const int plotRingQuarter = femod(ringQuarterIndex, 6);
    plotColor = plotRingQuarter * 2 + plotPhi + 1;
  }

  return plotColor;
  }*/
