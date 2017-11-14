#include "Cabling/ServicesChannel.hh"


void ServicesChannel::build(const int id, const ChannelSection& section, const bool isPositiveCablingSide, const int plotColor) {
  section_ = section;
  isPositiveCablingSide_ = isPositiveCablingSide;
  plotColor_ = plotColor;
  myid(id); 
};


OpticalChannel::OpticalChannel(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) {
  const ChannelSection& section = ChannelSection::B;

  const int number = computeChannelNumber(phiSectorRef, type, slot, isPositiveCablingSide);
  const int plotColor = computeChannelPlotColor(number);

  build(number, section, isPositiveCablingSide, plotColor);
};


/* Compute optical services channels.
 * They are the channels where the optical cables are routed when they exit the tracker.
 * They are closely related to the phiSector ref.
 */
const int OpticalChannel::computeChannelNumber(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const {
  int servicesChannel = 0;

  if (type == Category::PS10G) {
    if (phiSectorRef == 0) { servicesChannel = 2;  }
    else if (phiSectorRef == 1) { servicesChannel = 3; }
    else if (phiSectorRef == 2) { servicesChannel = 5; }
    else if (phiSectorRef == 3) { servicesChannel = 6; }
    else if (phiSectorRef == 4) { servicesChannel = 7; }
    else if (phiSectorRef == 5) { servicesChannel = 8; }
    else if (phiSectorRef == 6) { servicesChannel = 9; }
    else if (phiSectorRef == 7) { servicesChannel = 11; }
    else if (phiSectorRef == 8) { servicesChannel = 12; }   
  }

  else if (type == Category::PS5G) {
    if (phiSectorRef == 0) { servicesChannel = 1; }
    else if (phiSectorRef == 1) { servicesChannel = 2; }
    else if (phiSectorRef == 2) { servicesChannel = 4; }
    else if (phiSectorRef == 3) { servicesChannel = 5; }
    else if (phiSectorRef == 4) { servicesChannel = 7; }
    else if (phiSectorRef == 5) { servicesChannel = 8; }
    else if (phiSectorRef == 6) { servicesChannel = 10; }
    else if (phiSectorRef == 7) { servicesChannel = 11; }
    else if (phiSectorRef == 8) { servicesChannel = 12; }
  }

  else if (type == Category::SS) {
    if (slot == 1 || slot == 2) {
      if (phiSectorRef == 0) { servicesChannel = 1; }
      else if (phiSectorRef == 1) { servicesChannel = 3; }
      else if (phiSectorRef == 2) { servicesChannel = 4; }
      else if (phiSectorRef == 3) { servicesChannel = 5; }
      else if (phiSectorRef == 4) { servicesChannel = 6; }
      else if (phiSectorRef == 5) { servicesChannel = 8; }
      else if (phiSectorRef == 6) { servicesChannel = 9; }
      else if (phiSectorRef == 7) { servicesChannel = 10; }
      else if (phiSectorRef == 8) { servicesChannel = 11; }
    }
    else if (slot == 3) {
      if (phiSectorRef == 0) { servicesChannel = 1; }
      else if (phiSectorRef == 1) { servicesChannel = 2; }
      else if (phiSectorRef == 2) { servicesChannel = 3; }
      else if (phiSectorRef == 3) { servicesChannel = 4; }
      else if (phiSectorRef == 4) { servicesChannel = 6; }
      else if (phiSectorRef == 5) { servicesChannel = 7; }
      else if (phiSectorRef == 6) { servicesChannel = 9; }
      else if (phiSectorRef == 7) { servicesChannel = 10; }
      else if (phiSectorRef == 8) { servicesChannel = 12; }
    }
    else {
      if (phiSectorRef == 0) { servicesChannel = 1; }
      else if (phiSectorRef == 1) { servicesChannel = 2; }
      else if (phiSectorRef == 2) { servicesChannel = 3; }
      else if (phiSectorRef == 3) { servicesChannel = 4; }
      else if (phiSectorRef == 4) { servicesChannel = 6; }
      else if (phiSectorRef == 5) { servicesChannel = 7; }
      else if (phiSectorRef == 6) { servicesChannel = 9; }
      else if (phiSectorRef == 7) { servicesChannel = 10; }
      else if (phiSectorRef == 8) { servicesChannel = 12; }
    }
  }


  // NEGATIVE CABLING SIDE.
  // A given services channel is simplified as a straight line all along (Z).
  // THIS DEFINES THE CHANNEL NUMBERING ON THE (-Z) SIDE.

  // OPTION A: 
  // Channel 1A on (+Z) side becomes -1A on the (-Z) side, 1C on (+Z) side becomes -1C on (-Z) side, and so on.
  if (!isPositiveCablingSide) {
    servicesChannel *= -1;
  }

  // OPTION B (NOT PRESENTLY RETAINED)
  /*
  // This is the following transformation:
  // 1 -> 6
  // 2 -> 5
  // 3 -> 4
  // 7 -> 12
  // 8 -> 11
  // 9 -> 10
  // This is so that the numbering follows a rotation of 180 degrees around CMS_Y for the negative cabling side.
  // The services channel is then set to negative on negative cabling side.
  if (!isPositiveCablingSide) {
  double pivot = (servicesChannel <= 6 ? 3.5 : 9.5);
  servicesChannel = servicesChannel + round( 2. * (pivot - servicesChannel) );
  servicesChannel *= -1;
  servicesChannelSection = (servicesChannelSection == ChannelSection::A ? ChannelSection::C : ChannelSection::A);
  }
  */

  return servicesChannel;
}


/* Compute color associated to services channel.
 */
int OpticalChannel::computeChannelPlotColor(const int number) const {
  int plotColor = fabs(number);
  return plotColor;
}


PowerChannel::PowerChannel(const int semiPhiRegionRef, const bool isPositiveCablingSide) {

  const std::pair<int, ChannelSection> numberAndSection = computeChannelNumberAndSection(semiPhiRegionRef, isPositiveCablingSide);
  const int number = numberAndSection.first;
  const ChannelSection& section = numberAndSection.second;
  const int plotColor = computeChannelPlotColor(number, section, isPositiveCablingSide);

  build(number, section, isPositiveCablingSide, plotColor);
};


std::pair<int, ChannelSection> PowerChannel::computeChannelNumberAndSection(const int semiPhiRegionRef, const bool isPositiveCablingSide) const {

  int servicesChannel = 0;
  ChannelSection servicesChannelSection = ChannelSection::UNKNOWN;

  if (isPositiveCablingSide) {
    if (semiPhiRegionRef == 0) { servicesChannel = 1; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 1) { servicesChannel = 2; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 2) { servicesChannel = 2; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 3) { servicesChannel = 3; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 4) { servicesChannel = 4; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 5) { servicesChannel = 4; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 6) { servicesChannel = 5; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 7) { servicesChannel = 6; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 8) { servicesChannel = 6; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 9) { servicesChannel = 7; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 10) { servicesChannel = 8; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 11) { servicesChannel = 8; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 12) { servicesChannel = 9; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 13) { servicesChannel = 10; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 14) { servicesChannel = 10; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 15) { servicesChannel = 11; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 16) { servicesChannel = 12; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 17) { servicesChannel = 12; servicesChannelSection = ChannelSection::C; }
    else { std::cout << "ERROR: semiPhiRegionRef = " << semiPhiRegionRef << std::endl; }
  }

  else {
    if (semiPhiRegionRef == 0) { servicesChannel = -1; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 1) { servicesChannel = -1; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 2) { servicesChannel = -2; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 3) { servicesChannel = -3; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 4) { servicesChannel = -3; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 5) { servicesChannel = -4; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 6) { servicesChannel = -5; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 7) { servicesChannel = -5; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 8) { servicesChannel = -6; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 9) { servicesChannel = -7; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 10) { servicesChannel = -7; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 11) { servicesChannel = -8; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 12) { servicesChannel = -9; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 13) { servicesChannel = -9; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 14) { servicesChannel = -10; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 15) { servicesChannel = -11; servicesChannelSection = ChannelSection::A; }
    else if (semiPhiRegionRef == 16) { servicesChannel = -11; servicesChannelSection = ChannelSection::C; }
    else if (semiPhiRegionRef == 17) { servicesChannel = -12; servicesChannelSection = ChannelSection::A; }
    else { std::cout << "ERROR: semiPhiRegionRef = " << semiPhiRegionRef << std::endl; }
  }

  return std::make_pair(servicesChannel, servicesChannelSection);
}


int PowerChannel::computeChannelPlotColor(const int number, const ChannelSection& section, const bool isPositiveCablingSide) const {
  int plotColor = fabs(number);
  if ( (isPositiveCablingSide && section == ChannelSection::A)
       || (!isPositiveCablingSide && section == ChannelSection::C)
       ) {
    //if (section == ChannelSection::A) {
    plotColor += 12;
  }
  return plotColor;
}
