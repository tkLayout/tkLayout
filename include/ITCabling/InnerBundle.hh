#ifndef INNERBUNDLE_HH
#define INNERBUNDLE_HH

#include <vector>
#include <string>

#include "Property.hh"
//#include "Module.hh"
#include "ITCabling/GBT.hh"
#include "ITCabling/inner_cabling_functions.hh"


//namespace insur { class InnerDTC; }
//using insur::InnerDTC;


class InnerBundle : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef PtrVector<GBT> Container; 

public:
  InnerBundle(const int bundleId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int myBundleIndex);
  ~InnerBundle();

  // GBTS CONNECTED TO THE InnerBundle
  const Container& GBTs() const { return GBTs_; }
  const int numGBTs() const { return GBTs_.size(); }
  void addGBT(GBT* myGBT);

  // DTC TO WHICH THE InnerBundle IS CONECTED
  /*void setBundle(InnerBundle* bundle) { myBundle_ = bundle; }
    const InnerBundle* getBundle() const {
    if (!myBundle_) throw PathfulException("myBundle_ is nullptr");
    return myBundle_;
    }*/


  // GENERAL INFO
  const bool isPositiveZEnd() const { return isPositiveZEnd_; }
  const bool isPositiveXSide() const { return isPositiveXSide_; }
  const std::string subDetectorName() const { return subDetectorName_; }
  const int layerDiskNumber() const { return layerDiskNumber_; }
  const int bundleIndex() const { return myBundleIndex_; }

  const bool isBarrel() const { return isBarrel_; }
  //const bool isRingInnerEnd() const { return isRingInnerEnd_; }
  //const InnerBundleType powerChainType() const { return powerChainType_; }  

  //const int plotColor() const { return plotColor_; }

private:
  //const int computePlotColor(const bool isBarrel, const bool isPositiveZEnd, const int phiRef, const int ringQuarterIndex) const;

  Container GBTs_;

  //DTC* myDTC_ = nullptr;

  bool isPositiveZEnd_;
  bool isPositiveXSide_;
  std::string subDetectorName_;
  int layerDiskNumber_;
  int myBundleIndex_;
  
  bool isBarrel_;

  //int plotColor_;
};



#endif
