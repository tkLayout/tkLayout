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
  servicesChannel_ = computeServicesChannel(phiSectorRef, type, slot, isPositiveCablingSide);

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
const int Cable::computeServicesChannel(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const {
  int servicesChannel = 0;
  if (type == Category::PS10G) {
    if (phiSectorRef == 0) servicesChannel = 1;
    else if (phiSectorRef == 1) servicesChannel = 3;
    else if (phiSectorRef == 2) servicesChannel = 4;
    else if (phiSectorRef == 3) servicesChannel = 5;
    else if (phiSectorRef == 4) servicesChannel = 6;
    else if (phiSectorRef == 5) servicesChannel = 7;
    else if (phiSectorRef == 6) servicesChannel = 9;
    else if (phiSectorRef == 7) servicesChannel = 10;
    else if (phiSectorRef == 8) servicesChannel = 12;   
  }
  else if (type == Category::PS5G) {
    if (slot == 3) {
      if (phiSectorRef == 0) servicesChannel = 1;
      else if (phiSectorRef == 1) servicesChannel = 3;
      else if (phiSectorRef == 2) servicesChannel = 4;
      else if (phiSectorRef == 3) servicesChannel = 6;
      else if (phiSectorRef == 4) servicesChannel = 7;
      else if (phiSectorRef == 5) servicesChannel = 8;
      else if (phiSectorRef == 6) servicesChannel = 10;
      else if (phiSectorRef == 7) servicesChannel = 11;
      else if (phiSectorRef == 8) servicesChannel = 12;
    }
    else {
      if (phiSectorRef == 0) servicesChannel = 1;
      else if (phiSectorRef == 1) servicesChannel = 3;
      else if (phiSectorRef == 2) servicesChannel = 4;
      else if (phiSectorRef == 3) servicesChannel = 6;
      else if (phiSectorRef == 4) servicesChannel = 7;
      else if (phiSectorRef == 5) servicesChannel = 8;
      else if (phiSectorRef == 6) servicesChannel = 10;
      else if (phiSectorRef == 7) servicesChannel = 11;
      else if (phiSectorRef == 8) servicesChannel = 12;
    }
  }
  else if (type == Category::SS) {
    if (slot != 1 && slot != 2) {
      if (phiSectorRef == 0) servicesChannel = 1;
      else if (phiSectorRef == 1) servicesChannel = 2;
      else if (phiSectorRef == 2) servicesChannel = 4;
      else if (phiSectorRef == 3) servicesChannel = 5;
      else if (phiSectorRef == 4) servicesChannel = 6;
      else if (phiSectorRef == 5) servicesChannel = 7;
      else if (phiSectorRef == 6) servicesChannel = 9;
      else if (phiSectorRef == 7) servicesChannel = 10;
      else if (phiSectorRef == 8) servicesChannel = 12;
    }
    if (slot == 1 || slot == 2) {
      if (phiSectorRef == 0) servicesChannel = 1;
      else if (phiSectorRef == 1) servicesChannel = 2;
      else if (phiSectorRef == 2) servicesChannel = 4;
      else if (phiSectorRef == 3) servicesChannel = 6;
      else if (phiSectorRef == 4) servicesChannel = 7;
      else if (phiSectorRef == 5) servicesChannel = 8;
      else if (phiSectorRef == 6) servicesChannel = 10;
      else if (phiSectorRef == 7) servicesChannel = 11;
      else if (phiSectorRef == 8) servicesChannel = 12;
    }
  }

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
  }

  return servicesChannel;
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
