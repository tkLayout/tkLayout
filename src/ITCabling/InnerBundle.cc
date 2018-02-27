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

  plotColor_ = computePlotColor(myBundleIndex);
};


InnerBundle::~InnerBundle() {
  delete myDTC_;    // TO DO: switch to smart pointers and remove this!
  myDTC_ = nullptr;
}


void InnerBundle::addGBT(GBT* myGBT) { 
  GBTs_.push_back(myGBT);
}



const int InnerBundle::computePlotColor(const int bundleIndex) const {
  const int plotColor = femod(bundleIndex, 3) + 1;

  return plotColor;
}
