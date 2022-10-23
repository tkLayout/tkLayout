#include "InnerCabling/InnerBundle.hh"
#include "InnerCabling/InnerDTC.hh"


InnerBundle::InnerBundle(const int bundleId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int myBundleIndex) :
  isPositiveZEnd_(isPositiveZEnd),
  isPositiveXSide_(isPositiveXSide),
  subDetectorName_(subDetectorName),
  layerDiskNumber_(layerDiskNumber),
  myBundleIndex_(myBundleIndex)
{
  myid(bundleId);
  isBarrel_ = inner_cabling_functions::isBarrel(subDetectorName);

  plotColor_ = computePlotColor(isBarrel_, isPositiveZEnd, layerDiskNumber, myBundleIndex);
};


/*
 *  Connect a GBT to the Bundle.
 */
void InnerBundle::addGBT(GBT* myGBT) { 
  GBTs_.push_back(myGBT);
}



/*
 * Compute bundle color on website.
 * Bundles next to each other in space, must be of different colors.
 */
const int InnerBundle::computePlotColor(const bool isBarrel, const bool isPositiveZEnd, const int layerDiskNumber, const int bundleIndex) const {
  int plotColor = 0;

  if (isBarrel) {
    const int plotZEnd = (isPositiveZEnd ? 0 : 1);
    const int plotLayer = femod(layerDiskNumber, 4);
    const int plotIndex = femod(bundleIndex, 3);
    plotColor = plotZEnd * 6 + plotLayer * 3 + plotIndex + 1;  // ? plotZEnd * 7
  }
  else {
    const int plotLayer = femod(layerDiskNumber, 5);           //  ?
    const int plotIndex = femod(bundleIndex, 4);
    plotColor = plotLayer * 2 + plotIndex + 1;
  }

  return plotColor;
}
