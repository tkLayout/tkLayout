#ifndef CABLE_HH
#define CABLE_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Bundle.hh"


namespace insur { class DTC; }
using insur::DTC;


class Cable : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef PtrVector<Bundle> Container;
public:
  Cable(const int id, const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);
  ~Cable();

  // BUNDLES CONNECTED TO THE CABLE.
  const Container& bundles() const { return bundles_; }

  void addBundle(Bundle* b) { bundles_.push_back(b); }
  int numBundles() const { return bundles_.size(); }

  void assignPowerServicesChannels();

  // DTC THE CABLE IS CONNECTED TO.
  const DTC* getDTC() const { return myDTC_; }

  const Category& type() const { return type_; }
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const int phiSectorRef() const { return phiSectorRef_; }
  const int slot() const { return slot_; }
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }
  const int servicesChannel() const { return opticalServicesChannel_; }
  const ChannelSection& servicesChannelSection() const { return opticalServicesChannelSection_; }
  const int servicesChannelPlotColor() const { return opticalServicesChannelPlotColor_; }


private:
  const std::tuple<int, ChannelSection, int> computeOpticalServicesChannel(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const;
  const int computeOpticalServicesChannelPlotColor(const int servicesChannel, const ChannelSection& servicesChannelSection) const;
  std::pair<int, ChannelSection> computePowerServicesChannel(const int semiPhiRegionRef, const bool isPositiveCablingSide);

  void buildDTC(const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);
  const std::string computeDTCName(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const;
  
  Container bundles_;

  DTC* myDTC_ = nullptr;

  double phiSectorWidth_;
  int phiSectorRef_;
  Category type_;
  int slot_;
  bool isPositiveCablingSide_;
  int opticalServicesChannel_;
  ChannelSection opticalServicesChannelSection_;
  int opticalServicesChannelPlotColor_;
};



#endif
