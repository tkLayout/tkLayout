#ifndef GBT_HH
#define GBT_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Module.hh"
#include "ITCabling/inner_cabling_functions.hh"
//#include "ITCabling/ELink.hh"
#include "ITCabling/PowerChain.hh"


namespace insur { class HvLine; }
using insur::HvLine;


class GBT : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef PtrVector<Module> Container; 

public:
  GBT(PowerChain* myPowerChain, const std::string GBTId, const int myGBTPhiIndex, const int numELinksPerModule);
  ~GBT();

  // MODULES CONNECTED TO THE GBT
  const Container& modules() const { return modules_; }
  const int numModules() const { return modules_.size(); }
  void addModule(Module* m);
  const int numELinks() const { return numELinksPerModule_ * numModules(); } 

  // HIGH VOLTAGE LINE, TO WHICH THE MODULES OF THE POWER CHAIN ARE ALL CONNECTED
  /*const HvLine* getHvLine() const {
    if (!hvLine_) throw PathfulException("hvLine_ is nullptr");
    return hvLine_;
    }*/


  // GENERAL INFO ON THE POWERCHAIN
  const std::string GBTId() const { return myGBTId_; }
  const int GBTPhiIndex() const { return myGBTPhiIndex_; }
  const int numELinksPerModule() const { return numELinksPerModule_; }
  
  const bool isPositiveZEnd() const { return myPowerChain_->isPositiveZEnd(); }
  const bool isPositiveXSide() const { return myPowerChain_->isPositiveXSide(); }
  const bool isBarrel() const { return myPowerChain_->isBarrel(); }
  const bool isBarrelLong() const { return myPowerChain_->isBarrelLong(); }
  const std::string subDetectorName() const { return myPowerChain_->subDetectorName(); }
  const int layerDiskNumber() const { return myPowerChain_->layerDiskNumber(); }
  const int ringNumber() const { return myPowerChain_->ringNumber(); }
  const bool isRingInnerEnd() const { return myPowerChain_->isRingInnerEnd(); }
  const int ringQuarterIndex() const { return myPowerChain_->ringQuarterIndex(); }
  const int phiRef() const { return myPowerChain_->phiRef(); }
  

  //const int plotColor() const { return plotColor_; }

private:
  //const int computePlotColor(const bool isBarrel, const bool isPositiveZEnd, const int phiRef, const int ringQuarterIndex) const;

  Container modules_;

  PowerChain* myPowerChain_ = nullptr;

  std::string myGBTId_;
  int myGBTPhiIndex_;
  int numELinksPerModule_;

  //int plotColor_;
};



#endif
