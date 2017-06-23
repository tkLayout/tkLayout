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
  // Assign a servicesChannel to the cable
  servicesChannel_ = computeServicesChannel(phiSectorRef, type, slot, isPositiveCablingSide);

  // Build DTC asociated to the cable
  buildDTC(phiSectorWidth, phiSectorRef, type, slot, isPositiveCablingSide);  
};


Cable::~Cable() {
  delete myDTC_;       // TO DO: switch to smart pointers and remove this!
  myDTC_ = nullptr;
}


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
    if (slot != 3) {
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
      else if (phiSectorRef == 3) servicesChannel = 5;
      else if (phiSectorRef == 4) servicesChannel = 6;
      else if (phiSectorRef == 5) servicesChannel = 7;
      else if (phiSectorRef == 6) servicesChannel = 9;
      else if (phiSectorRef == 7) servicesChannel = 10;
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

  if (!isPositiveCablingSide) {
    double pivot = (servicesChannel <= 6 ? 3.5 : 9.5);
    servicesChannel = servicesChannel + round( 2. * (pivot - servicesChannel) );
    servicesChannel *= -1;
  }

  return servicesChannel;
}


// Build DTC asociated to the cable
void Cable::buildDTC(const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) {
  std::string dtcName = computeDTCName(phiSectorRef, type, slot, isPositiveCablingSide);
  DTC* dtc = GeometryFactory::make<DTC>(dtcName, phiSectorWidth, phiSectorRef, type, slot, isPositiveCablingSide);
  dtc->addCable(this);
  myDTC_ = dtc;
}


const std::string Cable::computeDTCName(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const {
  std::ostringstream dtcNameStream;
  if (!isPositiveCablingSide) dtcNameStream << cabling_negativePrefix;
  dtcNameStream << "_" << phiSectorRef << "_" << any2str(type) << "_" << slot;
  const std::string dtcName = dtcNameStream.str();
  return dtcName;
}
