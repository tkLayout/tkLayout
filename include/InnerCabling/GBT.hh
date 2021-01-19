#ifndef GBT_HH
#define GBT_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Module.hh"
#include "InnerCabling/PowerChain.hh"
#include "InnerCabling/inner_cabling_functions.hh"


namespace insur { class InnerBundle; }
using insur::InnerBundle;


/*
 * Low-Power Giga Bit Transceiver class.
 * All modules connected to a given GBT can be accessed, as well as the Fiber Bundle the GBT is connected to.
 * General info on the GBT is also provided.
 */

class GBT : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef std::vector<Module*> Container; 

public:
  GBT(PowerChain* myPowerChain, const std::string GBTId, const int myGBTIndexColor, const int numELinksPerModule);

  // MODULES CONNECTED TO THE GBT
  const Container& modules() const { return modules_; }
  const int numModules() const { return modules_.size(); }
  void addModule(Module* m);
  const int numELinks() const { return numELinksPerModule_ * numModules(); } 

  // BUNDLE TO WHICH THE GBT IS CONECTED
  void setBundle(InnerBundle* bundle) { myBundle_ = bundle; }
  const InnerBundle* getBundle() const {
    if (!myBundle_) throw PathfulException("myBundle_ is nullptr");
    return myBundle_;
  }

  // POWER CHAIN FOR ALL MODULES CONNECTED TO THE GBT
  const PowerChain* getPowerChain() const {
    if (!myPowerChain_) throw PathfulException("myPowerChain_ is nullptr");
    return myPowerChain_;
  }


  // GENERAL INFO ON THE GBT
  const std::string GBTId() const { return myGBTId_; }
  void setCMSSWId(const int cmsswId) { myGBTCMSSWId_ = cmsswId; }
  const int getCMSSWId() const { return myGBTCMSSWId_; }
  const int GBTIndexInPowerChain() const { return myGBTIndexInPowerChain_; }
  const int plotGBTIndexInPowerChain() const { return plotGBTIndexInPowerChain_; }
  const int numELinksPerModule() const { return numELinksPerModule_; }
  
  const bool isPositiveZEnd() const { return myPowerChain_->isPositiveZEnd(); }
  const bool isPositiveXSide() const { return myPowerChain_->isPositiveXSide(); }
  const bool isBarrel() const { return myPowerChain_->isBarrel(); }
  const bool isLongBarrel() const { return myPowerChain_->isLongBarrel(); }
  const std::string subDetectorName() const { return myPowerChain_->subDetectorName(); }
  const int layerDiskNumber() const { return myPowerChain_->layerDiskNumber(); }
  const int ringNumber() const { return myPowerChain_->ringNumber(); }
  const bool isSmallerAbsZHalfRing() const { return myPowerChain_->isSmallerAbsZHalfRing(); }
  const int halfRingIndex() const { return myPowerChain_->halfRingIndex(); }
  const bool isAtSmallerAbsZDeeInDoubleDisk() const { return myPowerChain_->isAtSmallerAbsZDeeInDoubleDisk(); }
  const bool isAtSmallerAbsZSideInDee() const { return myPowerChain_->isAtSmallerAbsZSideInDee(); }
  const int powerChainPhiRef() const { return myPowerChain_->phiRef(); }
  

  const int plotPowerChainColor() const { return plotPowerChainColor_; }

private:
  const int computePlotGBTIndexInPowerChain(const int myGBTIndexInPowerChain, PowerChain* myPowerChain) const;
  const int computePlotColor(const PowerChain* myPowerChain) const;

  Container modules_;

  PowerChain* myPowerChain_ = nullptr;
  InnerBundle* myBundle_ = nullptr;

  std::string myGBTId_;
  int myGBTCMSSWId_;
  int myGBTIndexInPowerChain_;
  int plotGBTIndexInPowerChain_;
  int numELinksPerModule_;

  int plotPowerChainColor_;
};



#endif
