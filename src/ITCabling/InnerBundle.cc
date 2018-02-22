#include "ITCabling/InnerBundle.hh"
#include "ITCabling/InnerDTC.hh"


InnerBundle::InnerBundle(const int bundleId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int myBundleIndex) :
  isPositiveZEnd_(isPositiveZEnd),
  isPositiveXSide_(isPositiveXSide),
  subDetectorName_(subDetectorName),
  layerDiskNumber_(layerDiskNumber),
  myBundleIndex_(myBundleIndex)
{
  myid(bundleId);
  isBarrel_ = inner_cabling_functions::isBarrel(subDetectorName);

  //bundleType_ = computePowerChainType(isBarrel_, layerDiskNumber, ringNumber_);

  //plotColor_ = computePlotColor(isBarrel_, isPositiveZEnd, phiRef, ringQuarterIndex);
};


InnerBundle::~InnerBundle() {
  delete myDTC_;    // TO DO: switch to smart pointers and remove this!
  myDTC_ = nullptr;
}


void InnerBundle::addGBT(GBT* myGBT) { 
  GBTs_.push_back(myGBT);
}



/*const int InnerBundle::computePlotColor(const bool isBarrel, const bool isPositiveZEnd, const int phiRef, const int ringQuarterIndex) const {
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
