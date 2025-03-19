#ifndef SERVICESCHANEL_HH
#define SERVICESCHANEL_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "OuterCabling/outer_cabling_constants.hh"


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
  std::string channelName_;
  std::string chanelNamePP1Convention_;
  bool isPositiveCablingSide_; 
  Container sections_;
  };*/


/* Section of a services channel.
   As all cables exit the Tracker, they are routed through a services channel.
 * A services channel is divided into 3 sections: A, B and C.
 * The optical bundles are always routed in section B.
 * The power cables are always routed in section A or C.
 * The cooling pipes are always routed in section A or C.
 * IMPORTANT NOTE: SECTION IN SERVICES CHANNEL = SECTION IN PP1 (Patch Panel 1).
 */
class ChannelSection : public PropertyObject, public Buildable {
  //typedef PtrVector<OuterBundle> Container;
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


/* Channel Section filled by optical bundles.
 */
class OpticalSection : public ChannelSection {
public:
  OpticalSection(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);

private:
  const int computeChannelNumber(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const;
  const int computeChannelPlotColor(const int number) const;
};


/* Channel Section filled by power cables. 
 * NB: In the CablingMap, 1 Bundle = 1 Power cable.
 */
class PowerSection : public ChannelSection {
public:
  PowerSection(const int semiPhiRegionRef, const bool isPositiveCablingSide);

private:
  std::pair<int, ChannelSlot> computeChannelNumberAndSlot(const int semiPhiRegionRef, const bool isPositiveCablingSide) const;
  const int computeChannelPlotColor(const int number, const ChannelSlot& channelSlot, const bool isPositiveCablingSide) const;
};



#endif
