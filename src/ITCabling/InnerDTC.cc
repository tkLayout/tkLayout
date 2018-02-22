#include "ITCabling/InnerDTC.hh"


InnerDTC::InnerDTC(const int DTCId, const bool isPositiveZEnd, const bool isPositiveXSide) :
  isPositiveZEnd_(isPositiveZEnd),
  isPositiveXSide_(isPositiveXSide)
{
  myid(DTCId);

  //DTCType_ = computePowerChainType(isBarrel_, layerDiskNumber, ringNumber_);

  plotColor_ = computePlotColor(DTCId);
};



void InnerDTC::addBundle(InnerBundle* bundle) { 
  bundles_.push_back(bundle);
}



const int InnerDTC::computePlotColor(const int DTCId) const {
  const int plotColor = DTCId % 10 + 4;

  return plotColor;
}
