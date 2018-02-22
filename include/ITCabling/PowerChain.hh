#ifndef POWERCHAIN_HH
#define POWERCHAIN_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Module.hh"
#include "ITCabling/inner_cabling_functions.hh"


namespace insur { class HvLine; }
using insur::HvLine;


class PowerChain : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef PtrVector<Module> Container; 

public:
  PowerChain(const int powerChainId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int ringQuarterIndex);
  ~PowerChain();

  // MODULES CONNECTED TO THE POWERCHAIN.
  const Container& modules() const { return modules_; }
  Container& modules() { return modules_; }
  const int numModules() const { return modules_.size(); }
  void addModule(Module* m);

  // HIGH VOLTAGE LINE, TO WHICH THE MODULES OF THE POWER CHAIN ARE ALL CONNECTED
  const HvLine* getHvLine() const {
    if (!hvLine_) throw PathfulException("hvLine_ is nullptr");
    return hvLine_;
  }

  // GENERAL INFO ON THE POWERCHAIN
  
  const bool isPositiveZEnd() const { return isPositiveZEnd_; }
  const bool isPositiveXSide() const { return isPositiveXSide_; }
  const std::string subDetectorName() const { return subDetectorName_; }
  const int layerDiskNumber() const { return layerDiskNumber_; }
  const int phiRef() const { return phiRef_; }
  const int ringQuarterIndex() const { return ringQuarterIndex_; }

  const bool isBarrel() const { return isBarrel_; }
  const int ringNumber() const { return ringNumber_; }
  const bool isRingInnerEnd() const { return isRingInnerEnd_; }

  const PowerChainType powerChainType() const { return powerChainType_; }

  const bool isBarrelLong() const {
    if (isBarrel()) return (numModules() == inner_cabling_maxNumModulesPerPowerChain);
    else return false;
  }

  const int plotColor() const { return plotColor_; }

private:
  const PowerChainType computePowerChainType(const bool isBarrel, const int layerDiskNumber, const int ringNumber) const;
  const int computePlotColor(const bool isBarrel, const bool isPositiveZEnd, const int phiRef, const int ringQuarterIndex) const;

  void buildHvLine(const int powerChainId);
  const std::string computeHvLineName(const int powerChainId) const;

  Container modules_;

  HvLine* hvLine_ = nullptr;

  bool isPositiveZEnd_;
  bool isPositiveXSide_;
  std::string subDetectorName_;
  int layerDiskNumber_;
  int phiRef_;
  int ringQuarterIndex_;
  
  bool isBarrel_;
  int ringNumber_;
  bool isRingInnerEnd_;

  PowerChainType powerChainType_;

  int plotColor_;
};



#endif
