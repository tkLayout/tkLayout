#ifndef OUTERBUNDLE_HH
#define OUTERBUNDLE_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Module.hh"
#include "OuterCabling/PhiPosition.hh"
#include "OuterCabling/ServicesChannel.hh"


namespace insur { class OuterCable; }
using insur::OuterCable;


class OuterBundle : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef std::vector<Module*> Container; 

public:
  OuterBundle(const int id, const int stereoBundleId, const Category& type, const std::string subDetectorName, const int layerDiskNumber, const PhiPosition& phiPosition, const bool isPositiveCablingSide, const bool isTiltedPart);
  //~OuterBundle();

  // MODULES CONNECTED TO THE BUNDLE.
  const Container& modules() const { return modules_; }
  Container& modules() { return modules_; }
  const int numModules() const { return modules_.size(); }
  void addModule(Module* m) { modules_.push_back(m); }

  // CABLE THE BUNDLE IS CONNECTED TO.
  const OuterCable* getCable() const {
    if (!cable_) throw PathfulException("cable_ is nullptr");
    return cable_;
  }
  void setCable(OuterCable* cable) { cable_ = cable; }

  // GENERAL INFO ON THE BUNDLE
  const Category& type() const { return type_; }
  const std::string subDetectorName() const { return subDetectorName_; }
  const int layerDiskNumber() const { return layerDiskNumber_; }
  const PhiPosition& phiPosition() const { return phiPosition_; }  
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }

  const bool isTiltedPart() const { return isTiltedPart_; }
  const bool isBarrel() const { return (subDetectorName_ == outer_cabling_tbps || subDetectorName_ == outer_cabling_tb2s); }
  const bool isBarrelPSFlatPart() const { return (isBarrel() && type_ != Category::SS && !isTiltedPart_); }

  const int plotColor() const { return plotColor_; }

  // PHI INFORMATION FROM MODULES CONNECTED TO THE BUNDLE
  void moveMaxPhiModuleFromOtherBundle(OuterBundle* otherBundle);
  void moveMinPhiModuleFromOtherBundle(OuterBundle* otherBundle);

  const double minPhi() const;
  const double maxPhi() const;
  const double meanPhi() const;

  Module* minPhiModule() const;
  Module* maxPhiModule() const; 

  // SERVICES CHANNELS INFORMATION
  // VERY IMPORTANT: connection scheme from modules to optical bundles = connection scheme from modules to power cables.
  // As a result, 1 single Bundle object is used for both schemes.
  // Regarding the connections to services channels, each Bundle is then assigned:
  // - 1 Optical Services Channel Section (considering the Bundle as an optical Bundle);
  // - 1 Power Services Channels section (making as if the Bundle is a power cable);

  // Optical
  const ChannelSection* opticalChannelSection() const;
  // Power
  const ChannelSection* powerChannelSection() const {
    if (!powerChannelSection_) throw PathfulException("powerChannelSection_ is nullptr");
    return powerChannelSection_.get(); 
  }
  void setPowerChannelSection(std::unique_ptr<const ChannelSection> powerChannelSection) {
    powerChannelSection_ = std::move(powerChannelSection);
  }
  // Used to compute the power channel section.
  const int tiltedBundleId() const;
  const int stereoBundleId() const; // Id of the bundle located on the other cabling side, by rotation of 180° around CMS_Y.
  void setIsPowerRoutedToBarrelLowerSemiNonant(const bool isLower);
  const bool isPowerRoutedToBarrelLowerSemiNonant() const;
  

private:
  const int computePlotColor(const int id, const bool isPositiveCablingSide) const;

  Container modules_;

  OuterCable* cable_ = nullptr;

  Category type_;
  std::string subDetectorName_;
  int layerDiskNumber_;
  PhiPosition phiPosition_;
  bool isPositiveCablingSide_;
  bool isTiltedPart_;

  int plotColor_;

  std::unique_ptr<const ChannelSection> powerChannelSection_;

  int stereoBundleId_;
  bool isPowerRoutedToBarrelLowerSemiNonant_;
};



#endif
