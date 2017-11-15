#ifndef SERVICESCHANEL_HH
#define SERVICESCHANEL_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Cabling/cabling_constants.hh"


// !!! TO DO !!!
// FOR THE MOMENT, ChannelSection info is only accessible directly from the Bundle class.
// This is sufficient, but it would be good to also have a gathering of the info per ServicesChannel.
// A map of all Tracker ServicesChannels could then be assigned to the CablingMap.
// (as is the map of Bundles and the map of Cables).

// As a result, what needs to be done is: Add ServicesChannel class, containing 3 ChannelSections.
// Example: ServicesChannel with number 12, on (+Z) side, contains sections A, B, and C.
// Then, each ChannelSection should contain pointeres to the bundles which it contains.
// A map of all Tracker ServicesChannels could then be assigned to the CablingMap.



/*class ServicesChannel : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef PtrVector<ChannelSection> Container;
public:
  const int getNumber() const { return channelNumber_; }
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }
  const Container& getSections() const { return sections_; }

protected:
  int channelNumber_;
  bool isPositiveCablingSide_; 
  Container sections_;
  };*/


class ChannelSection : public PropertyObject, public Buildable {
  //typedef PtrVector<Bundle> Container;
public:
  const int channelNumber() const { return channelNumber_; }
  const ChannelSlot& channelSlot() const { return channelSlot_; }
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }
  const int plotColor() const { return plotColor_; }
  //const Container& getBundles() const { return bundles_; }

protected:
  void build(const int channelNumber, const ChannelSlot& channelSlot, const bool isPositiveCablingSide, const int plotColor);

  int channelNumber_;
  ChannelSlot channelSlot_ = ChannelSlot::UNKNOWN;
  bool isPositiveCablingSide_;  
  int plotColor_;
  //Container bundles_;
};



class OpticalSection : public ChannelSection {
public:
  OpticalSection(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);

private:
  const int computeChannelNumber(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const;
  const int computeChannelPlotColor(const int number) const;

};



class PowerSection : public ChannelSection {
public:
  PowerSection(const int semiPhiRegionRef, const bool isPositiveCablingSide);

private:
  std::pair<int, ChannelSlot> computeChannelNumberAndSlot(const int semiPhiRegionRef, const bool isPositiveCablingSide) const;
  const int computeChannelPlotColor(const int number, const ChannelSlot& channelSlot, const bool isPositiveCablingSide) const;

};







#endif
