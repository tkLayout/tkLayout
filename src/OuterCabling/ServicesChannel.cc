#include "OuterCabling/ServicesChannel.hh"


void ChannelSection::build(const int channelNumber, const ChannelSlot& channelSlot, const bool isPositiveCablingSide, const int plotColor) {
  channelNumber_ = channelNumber;
  channelSlot_ = channelSlot;
  isPositiveCablingSide_ = isPositiveCablingSide;
  plotColor_ = plotColor;
};


/* Channel Section filled by optical bundles. 
 */
OpticalSection::OpticalSection(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) {
  const ChannelSlot& channelSlot = ChannelSlot::B;  // The optical bundles are always routed in section B.

  const int channelNumber = computeChannelNumber(phiSectorRef, type, slot, isPositiveCablingSide);
  const int plotColor = computeChannelPlotColor(channelNumber);

  build(channelNumber, channelSlot, isPositiveCablingSide, plotColor);
};


/* Compute the number of the Channel, through which the optical bundles are routed when they exit the Tracker.
 * This number is closely related to the phiSector ref.
 * This is done in a way to avoid crossings between fiber bundles, as they exit the tracker.
 */
const int OpticalSection::computeChannelNumber(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const {
  int channelNumber = 0;

  if (type == Category::PS10G) {
    if (phiSectorRef == 0) { channelNumber = 2;  }
    else if (phiSectorRef == 1) { channelNumber = 3; }
    else if (phiSectorRef == 2) { channelNumber = 5; }
    else if (phiSectorRef == 3) { channelNumber = 6; }
    else if (phiSectorRef == 4) { channelNumber = 7; }
    else if (phiSectorRef == 5) { channelNumber = 8; }
    else if (phiSectorRef == 6) { channelNumber = 9; }
    else if (phiSectorRef == 7) { channelNumber = 11; }
    else if (phiSectorRef == 8) { channelNumber = 12; }   
  }

  else if (type == Category::PS5G) {
    if (phiSectorRef == 0) { channelNumber = 1; }
    else if (phiSectorRef == 1) { channelNumber = 2; }
    else if (phiSectorRef == 2) { channelNumber = 4; }
    else if (phiSectorRef == 3) { channelNumber = 5; }
    else if (phiSectorRef == 4) { channelNumber = 7; }
    else if (phiSectorRef == 5) { channelNumber = 8; }
    else if (phiSectorRef == 6) { channelNumber = 10; }
    else if (phiSectorRef == 7) { channelNumber = 11; }
    else if (phiSectorRef == 8) { channelNumber = 12; }
  }

  else if (type == Category::SS) {
    if (slot == 1 || slot == 2) {
      if (phiSectorRef == 0) { channelNumber = 1; }
      else if (phiSectorRef == 1) { channelNumber = 3; }
      else if (phiSectorRef == 2) { channelNumber = 4; }
      else if (phiSectorRef == 3) { channelNumber = 5; }
      else if (phiSectorRef == 4) { channelNumber = 6; }
      else if (phiSectorRef == 5) { channelNumber = 8; }
      else if (phiSectorRef == 6) { channelNumber = 9; }
      else if (phiSectorRef == 7) { channelNumber = 10; }
      else if (phiSectorRef == 8) { channelNumber = 11; }
    }
    else if (slot == 3) {
      if (phiSectorRef == 0) { channelNumber = 1; }
      else if (phiSectorRef == 1) { channelNumber = 2; }
      else if (phiSectorRef == 2) { channelNumber = 3; }
      else if (phiSectorRef == 3) { channelNumber = 4; }
      else if (phiSectorRef == 4) { channelNumber = 6; }
      else if (phiSectorRef == 5) { channelNumber = 7; }
      else if (phiSectorRef == 6) { channelNumber = 9; }
      else if (phiSectorRef == 7) { channelNumber = 10; }
      else if (phiSectorRef == 8) { channelNumber = 12; }
    }
    else {
      if (phiSectorRef == 0) { channelNumber = 1; }
      else if (phiSectorRef == 1) { channelNumber = 2; }
      else if (phiSectorRef == 2) { channelNumber = 3; }
      else if (phiSectorRef == 3) { channelNumber = 4; }
      else if (phiSectorRef == 4) { channelNumber = 6; }
      else if (phiSectorRef == 5) { channelNumber = 7; }
      else if (phiSectorRef == 6) { channelNumber = 9; }
      else if (phiSectorRef == 7) { channelNumber = 10; }
      else if (phiSectorRef == 8) { channelNumber = 12; }
    }
  }

  // NEGATIVE CABLING SIDE.
  // A given services channel is simplified as a straight line all along (Z).
  // THIS DEFINES THE CHANNEL NUMBERING ON THE (-Z) SIDE.

  // OPTION A: 
  // Channel 1A on (+Z) side becomes -1A on the (-Z) side, 1C on (+Z) side becomes -1C on (-Z) side, and so on.
  if (!isPositiveCablingSide) {
    channelNumber *= -1;
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
  double pivot = (channelNumber <= 6 ? 3.5 : 9.5);
  channelNumber = channelNumber + round( 2. * (pivot - channelNumber) );
  channelNumber *= -1;
  servicesChannelSlot = (servicesChannelSlot == ChannelSlot::A ? ChannelSlot::C : ChannelSlot::A);
  }
  */

  return channelNumber;
}


/* Compute color associated to optical channel section.
 */
const int OpticalSection::computeChannelPlotColor(const int channelNumber) const {
  int plotColor = fabs(channelNumber);
  return plotColor;
}


/* Channel Section filled by power cables. 
 * NB: In the CablingMap, 1 Bundle = 1 Power cable.
 */
PowerSection::PowerSection(const int semiPhiRegionRef, const bool isPositiveCablingSide) {

  const std::pair<int, ChannelSlot> channelNumberAndSlot = computeChannelNumberAndSlot(semiPhiRegionRef, isPositiveCablingSide);
  const int channelNumber = channelNumberAndSlot.first;
  const ChannelSlot& channelSlot = channelNumberAndSlot.second;
  const int plotColor = computeChannelPlotColor(channelNumber, channelSlot, isPositiveCablingSide);

  build(channelNumber, channelSlot, isPositiveCablingSide, plotColor);
};


/* Compute the number of the Channel, through which the power cables are routed when they exit the Tracker.
 * Also compute the channel slot on which the power cables are routed.
 * This is done in a way that leaves space free for the cooling pipes routing.
 * Positive cabling side: 1C, 3C, 5C, 7C, 9C, 11C used for cooling pipes.
 * Negative cabling side: -2A, -4A, -6A, -8A, -10A, -12A used for cooling pipes.
 */
std::pair<int, ChannelSlot> PowerSection::computeChannelNumberAndSlot(const int semiPhiRegionRef, const bool isPositiveCablingSide) const {

  int channelNumber = 0;
  ChannelSlot channelSlot = ChannelSlot::UNKNOWN;

  if (isPositiveCablingSide) {
    if (semiPhiRegionRef == 0) { channelNumber = 1; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 1) { channelNumber = 2; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 2) { channelNumber = 2; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 3) { channelNumber = 3; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 4) { channelNumber = 4; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 5) { channelNumber = 4; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 6) { channelNumber = 5; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 7) { channelNumber = 6; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 8) { channelNumber = 6; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 9) { channelNumber = 7; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 10) { channelNumber = 8; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 11) { channelNumber = 8; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 12) { channelNumber = 9; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 13) { channelNumber = 10; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 14) { channelNumber = 10; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 15) { channelNumber = 11; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 16) { channelNumber = 12; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 17) { channelNumber = 12; channelSlot = ChannelSlot::C; }
    else { std::cout << "ERROR: semiPhiRegionRef = " << semiPhiRegionRef << std::endl; }
  }

  else {
    if (semiPhiRegionRef == 0) { channelNumber = -1; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 1) { channelNumber = -1; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 2) { channelNumber = -2; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 3) { channelNumber = -3; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 4) { channelNumber = -3; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 5) { channelNumber = -4; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 6) { channelNumber = -5; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 7) { channelNumber = -5; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 8) { channelNumber = -6; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 9) { channelNumber = -7; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 10) { channelNumber = -7; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 11) { channelNumber = -8; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 12) { channelNumber = -9; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 13) { channelNumber = -9; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 14) { channelNumber = -10; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 15) { channelNumber = -11; channelSlot = ChannelSlot::A; }
    else if (semiPhiRegionRef == 16) { channelNumber = -11; channelSlot = ChannelSlot::C; }
    else if (semiPhiRegionRef == 17) { channelNumber = -12; channelSlot = ChannelSlot::C; }
    else { std::cout << "ERROR: semiPhiRegionRef = " << semiPhiRegionRef << std::endl; }
  }

  return std::make_pair(channelNumber, channelSlot);
}


/* Compute color associated to power channel section.
 */
const int PowerSection::computeChannelPlotColor(const int channelNumber, const ChannelSlot& channelSlot, const bool isPositiveCablingSide) const {
  int plotColor = fabs(channelNumber);
  if ( (isPositiveCablingSide && channelSlot == ChannelSlot::C)
       || (!isPositiveCablingSide && channelSlot == ChannelSlot::A)
       ) {
    //if (channelSlot == ChannelSlot::A) {
    plotColor += 12;  // This is used to activate transparency.
  }
  return plotColor;
}
