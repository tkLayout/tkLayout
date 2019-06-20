#ifndef POWERCHAIN_HH
#define POWERCHAIN_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Module.hh"
#include "InnerCabling/inner_cabling_functions.hh"


namespace insur { class HvLine; }
using insur::HvLine;


/*
 * Power chain class.
 * All modules connected to a given Power Chain can be accessed, as well as the HV line they are connected to.
 * General info on the Power Chain is also provided.
 */
class PowerChain : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef std::vector<Module*> Container; 

public:
  PowerChain(const int powerChainId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const bool isLongBarrel, const int ringQuarterIndex);

  // MODULES CONNECTED TO THE POWER CHAIN.
  const Container& modules() const { return modules_; }
  Container& modules() { return modules_; }
  const int numModules() const { return modules_.size(); }
  void addModule(Module* m);

  // HIGH VOLTAGE LINE, TO WHICH THE MODULES OF THE POWER CHAIN ARE ALL CONNECTED
  const HvLine* getHvLine() const {
    if (!hvLine_) throw PathfulException("hvLine_ is nullptr");
    return hvLine_.get();
  }

  // GENERAL INFO ON THE POWER CHAIN
  const bool isPositiveZEnd() const { return isPositiveZEnd_; }
  const bool isPositiveXSide() const { return isPositiveXSide_; }
  const std::string subDetectorName() const { return subDetectorName_; }
  const int layerDiskNumber() const { return layerDiskNumber_; }
  const int phiRef() const { return phiRef_; }
  const int ringQuarterIndex() const { return ringQuarterIndex_; }

  const bool isBarrel() const { return isBarrel_; }
  const int ringNumber() const { return ringNumber_; }
  const bool isSmallerAbsZRingSide() const { return isSmallerAbsZRingSide_; }

  const PowerChainType powerChainType() const { return powerChainType_; }

  const bool isLongBarrel() const {
    if (isBarrel()) return isLongBarrel_;
    else return false;
  }

  const int plotColor() const { return plotColor_; }

private:
  const PowerChainType computePowerChainType(const bool isBarrel, const int layerDiskNumber, const int ringNumber) const;
  const int computePlotColor(const bool isBarrel, const bool isPositiveZEnd, const int phiRef, const int ringQuarterIndex) const;

  void buildHvLine(const int powerChainId);
  const std::string computeHvLineName(const int powerChainId) const;

  Container modules_;

  std::unique_ptr<HvLine> hvLine_; // PowerChain owns HvLine

  bool isPositiveZEnd_;
  bool isPositiveXSide_;
  std::string subDetectorName_;
  int layerDiskNumber_;
  int phiRef_;
  bool isLongBarrel_;
  int ringQuarterIndex_;
  
  bool isBarrel_;
  int ringNumber_;
  bool isSmallerAbsZRingSide_;

  PowerChainType powerChainType_;

  int plotColor_;
};



#endif
