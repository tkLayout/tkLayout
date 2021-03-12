#ifndef INNERBUNDLE_HH
#define INNERBUNDLE_HH

#include <string>
#include <vector>

#include "Property.hh"

#include "InnerCabling/GBT.hh"
#include "InnerCabling/inner_cabling_functions.hh"

namespace insur {
class InnerDTC;
}
using insur::InnerDTC;

/*
 * Inner Tracker Fiber Bundle class.
 * All GBTs connected to a given Bundle can be accessed, as well as the Fiber
 * Cable it is connected to. General info on the Fiber Bundle is also provided.
 */
class InnerBundle : public PropertyObject,
                    public Buildable,
                    public Identifiable<int> {
  typedef std::vector<GBT *> Container;

public:
  InnerBundle(const int bundleId, const bool isPositiveZEnd,
              const bool isPositiveXSide, const std::string subDetectorName,
              const int layerDiskNumber, const int myBundleIndex);

  // GBTS CONNECTED TO THE BUNDLE
  const Container &GBTs() const { return GBTs_; }
  const int numGBTs() const { return GBTs_.size(); }
  void addGBT(GBT *myGBT);

  // DTC TO WHICH THE BUNDLE IS CONNECTED
  void setDTC(InnerDTC *myDTC) { myDTC_ = myDTC; }
  const InnerDTC *getDTC() const {
    if (!myDTC_)
      throw PathfulException("myDTC_ is nullptr");
    return myDTC_;
  }

  // GENERAL INFO
  const bool isPositiveZEnd() const { return isPositiveZEnd_; }
  const bool isPositiveXSide() const { return isPositiveXSide_; }
  const std::string subDetectorName() const { return subDetectorName_; }
  const int layerDiskNumber() const { return layerDiskNumber_; }
  const int bundleIndex() const { return myBundleIndex_; }

  const bool isBarrel() const { return isBarrel_; }

  const int plotColor() const { return plotColor_; }

private:
  const int computePlotColor(const bool isBarrel, const bool isPositiveZEnd,
                             const int layerDiskNumber,
                             const int bundleIndex) const;

  Container GBTs_;

  InnerDTC *myDTC_ = nullptr;

  bool isPositiveZEnd_;
  bool isPositiveXSide_;
  std::string subDetectorName_;
  int layerDiskNumber_;
  int myBundleIndex_;

  bool isBarrel_;

  int plotColor_;
};

#endif
