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

  // DTC THE CABLE IS CONNECTED TO.
  const DTC* getDTC() const { return myDTC_; }

  const Category& type() const { return type_; }
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const int phiSectorRef() const { return phiSectorRef_; }
  const int slot() const { return slot_; }
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }
  const int servicesChannel() const { return servicesChannel_; }
  const ChannelSection& servicesChannelSection() const { return servicesChannelSection_; }
  const int servicesChannelPlotColor() const { return servicesChannelPlotColor_; }


private:
  const std::tuple<int, ChannelSection, int> computeServicesChannel(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const;
  const int computeServicesChannelPlotColor(const int servicesChannel, const ChannelSection& servicesChannelSection) const;
  void buildDTC(const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);
  const std::string computeDTCName(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const;
  
  Container bundles_;

  DTC* myDTC_ = nullptr;

  double phiSectorWidth_;
  int phiSectorRef_;
  Category type_;
  int slot_;
  bool isPositiveCablingSide_;
  int servicesChannel_;
  ChannelSection servicesChannelSection_;
  int servicesChannelPlotColor_;
};



#endif
