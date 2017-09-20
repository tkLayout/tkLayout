#include "Cabling/Cable.hh"
#include "Cabling/DTC.hh"


Cable::Cable(const int id, const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) : 
  phiSectorWidth_(phiSectorWidth),
  phiSectorRef_(phiSectorRef),
  type_(type),
  slot_(slot),
  isPositiveCablingSide_(isPositiveCablingSide)
{ 
  myid(id);
  // ASSIGN A SERVICESCHANNEL TO THE CABLE
  const std::tuple<int, ChannelSection, int>& servicesChannelInfo = computeServicesChannel(phiSectorRef, type, slot, isPositiveCablingSide);
  servicesChannel_ = std::get<0>(servicesChannelInfo);
  servicesChannelSection_ = std::get<1>(servicesChannelInfo);
  servicesChannelPlotColor_ = std::get<2>(servicesChannelInfo);

  // BUILD DTC ASOCIATED TO THE CABLE
  buildDTC(phiSectorWidth, phiSectorRef, type, slot, isPositiveCablingSide);  
};


Cable::~Cable() {
  delete myDTC_;       // TO DO: switch to smart pointers and remove this!
  myDTC_ = nullptr;
}


/* Compute services channels.
 * They are the channels where the optical cables are routed when they exit the tracker.
 * They are closely related to the phiSector ref.
 */
const std::tuple<int, ChannelSection, int> Cable::computeServicesChannel(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const {
  int servicesChannel = 0;
  ChannelSection servicesChannelSection = ChannelSection::UNKNOWN;

  if (type == Category::PS10G) {
    if (phiSectorRef == 0) { servicesChannel = 1; servicesChannelSection = ChannelSection::C; }
    else if (phiSectorRef == 1) { servicesChannel = 2; servicesChannelSection = ChannelSection::C; }
    else if (phiSectorRef == 2) { servicesChannel = 3; servicesChannelSection = ChannelSection::C; }
    else if (phiSectorRef == 3) { servicesChannel = 4; servicesChannelSection = ChannelSection::C; }
    else if (phiSectorRef == 4) { servicesChannel = 6; servicesChannelSection = ChannelSection::C; }
    else if (phiSectorRef == 5) { servicesChannel = 8; servicesChannelSection = ChannelSection::C; }
    else if (phiSectorRef == 6) { servicesChannel = 9; servicesChannelSection = ChannelSection::C; }
    else if (phiSectorRef == 7) { servicesChannel = 10; servicesChannelSection = ChannelSection::C; }
    else if (phiSectorRef == 8) { servicesChannel = 12; servicesChannelSection = ChannelSection::C; }   
  }
  else if (type == Category::PS5G) {
    if (slot == 3) {
      if (phiSectorRef == 0) { servicesChannel = 1; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 1) { servicesChannel = 2; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 2) { servicesChannel = 3; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 3) { servicesChannel = 4; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 4) { servicesChannel = 6; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 5) { servicesChannel = 8; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 6) { servicesChannel = 9; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 7) { servicesChannel = 10; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 8) { servicesChannel = 12; servicesChannelSection = ChannelSection::C; }
    }
    else {
      if (phiSectorRef == 0) { servicesChannel = 1; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 1) { servicesChannel = 2; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 2) { servicesChannel = 4; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 3) { servicesChannel = 5; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 4) { servicesChannel = 6; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 5) { servicesChannel = 8; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 6) { servicesChannel = 10; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 7) { servicesChannel = 11; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 8) { servicesChannel = 12; servicesChannelSection = ChannelSection::C; }
    }
  }
  else if (type == Category::SS) {
    if (slot == 1 || slot == 2) {
      if (phiSectorRef == 0) { servicesChannel = 2; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 1) { servicesChannel = 3; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 2) { servicesChannel = 4; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 3) { servicesChannel = 6; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 4) { servicesChannel = 7; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 5) { servicesChannel = 8; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 6) { servicesChannel = 9; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 7) { servicesChannel = 10; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 8) { servicesChannel = 12; servicesChannelSection = ChannelSection::A; }
    }
    else if (slot == 3) {
      if (phiSectorRef == 0) { servicesChannel = 2; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 1) { servicesChannel = 3; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 2) { servicesChannel = 4; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 3) { servicesChannel = 6; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 4) { servicesChannel = 7; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 5) { servicesChannel = 8; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 6) { servicesChannel = 9; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 7) { servicesChannel = 10; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 8) { servicesChannel = 12; servicesChannelSection = ChannelSection::A; }
    }
    else {
      if (phiSectorRef == 0) { servicesChannel = 2; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 1) { servicesChannel = 4; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 2) { servicesChannel = 5; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 3) { servicesChannel = 6; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 4) { servicesChannel = 7; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 5) { servicesChannel = 8; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 6) { servicesChannel = 10; servicesChannelSection = ChannelSection::A; }
      else if (phiSectorRef == 7) { servicesChannel = 11; servicesChannelSection = ChannelSection::C; }
      else if (phiSectorRef == 8) { servicesChannel = 12; servicesChannelSection = ChannelSection::A; }
    }
  }


  // VERY IMPORTANT!
  // COMPUTE SERVICES CHANNEL PLOT COLOR, FOR BOTH SIDES, BEFORE THE CHANNEL NUMBERING IS MIRRORED FOR (-Z) SIDE.
  // THIS IS BECAUSE A GIVEN SERVICE CHANNEL IS IDENTICAL ALL ALONG Z (going from (+z) to (-Z) side).
  // HENCE, EVEN IF THE NUMBERING IS MIRRORED, WE DONT CARE, AND WE ACTUALLY WANT A GIVEN CHANNEL COLORED BY A UNIQUE COLOR!!
  int servicesChannelPlotColor = computeServicesChannelPlotColor(servicesChannel, servicesChannelSection);
  //std::cout << "servicesChannelPlotColor = " << servicesChannelPlotColor << std::endl;


  // NEGATIVE CABLING SIDE.
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

  return std::make_tuple(servicesChannel, servicesChannelSection, servicesChannelPlotColor);
}


const int Cable::computeServicesChannelPlotColor(const int servicesChannel, const ChannelSection& servicesChannelSection) const {
  int plotColor = 0;
  plotColor = servicesChannel;
  //if (servicesChannelSection == ChannelSection::C) plotColor += 12 * 1;
  return plotColor;
}


/* Build DTC asociated to the cable.
 */
void Cable::buildDTC(const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) {
  std::string dtcName = computeDTCName(phiSectorRef, type, slot, isPositiveCablingSide);
  DTC* dtc = GeometryFactory::make<DTC>(dtcName, phiSectorWidth, phiSectorRef, type, slot, isPositiveCablingSide);
  dtc->addCable(this);
  myDTC_ = dtc;
}


/* Compute DTC name.
 */
const std::string Cable::computeDTCName(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const {
  std::ostringstream dtcNameStream;
  if (!isPositiveCablingSide) dtcNameStream << cabling_negativePrefix << "_";
  dtcNameStream << phiSectorRef << "_" << any2str(type) << "_" << slot;
  const std::string dtcName = dtcNameStream.str();
  return dtcName;
}
