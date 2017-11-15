#ifndef SERVICESCHANEL_HH
#define SERVICESCHANEL_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Cabling/cabling_constants.hh"


class ChannelSection : public PropertyObject, public Buildable {
public:
  const int channelNumber() const { return channelNumber_; }
  const ChannelSlot& channelSlot() const { return channelSlot_; }
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }
  const int plotColor() const { return plotColor_; }

protected:
  void build(const int channelNumber, const ChannelSlot& channelSlot, const bool isPositiveCablingSide, const int plotColor);

  int channelNumber_;
  ChannelSlot channelSlot_ = ChannelSlot::UNKNOWN;
  bool isPositiveCablingSide_;  
  int plotColor_;
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
