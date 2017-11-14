#ifndef SERVICESCHANEL_HH
#define SERVICESCHANEL_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Cabling/cabling_constants.hh"


class ServicesChannel : public PropertyObject, public Buildable, public Identifiable<int> {
public:
  ServicesChannel(const int id, const ChannelSection& section, const bool isPositiveCablingSide, const int plotColor);

  const ChannelSection& section() const { return section_; }
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }
 
  const int plotColor() const { return plotColor_; }

protected:
  ChannelSection section_ = ChannelSection::UNKNOWN;
  bool isPositiveCablingSide_;  
  int plotColor_;
};



class OpticalChannel : public ServicesChannel {
public:
  OpticalChannel(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);

private:
  const int computeOpticalChannelNumber(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const;
  int computeOpticalChannelPlotColor(const int number) const;

};



class PowerChannel : public ServicesChannel {
public:
  PowerChannel(const int semiPhiRegionRef, const bool isPositiveCablingSide);

private:
  std::pair<int, ChannelSection> computePowerChannelNumberAndSection(const int semiPhiRegionRef, const bool isPositiveCablingSide) const;
  int computePowerChannelPlotColor(const int number, const ChannelSection& section, const bool isPositiveCablingSide) const;

};







#endif
