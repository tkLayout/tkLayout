#include "Cabling/Cable.hh"
#include "Cabling/DTC.hh"


//void Cable::build() {}



Cable::Cable(int id, const double phiSectorWidth, int phiSectorRef, std::string type, int slot) {
  // Generic cable characteristics
  myid(id);
  phiSectorWidth_ = phiSectorWidth;
  phiSectorRef_ = phiSectorRef;
  type_ = type;
  slot_ = slot;

  // Assign cable a servicesChannel
  servicesChannel_ = computeServicesChannel(phiSectorRef, type, slot);



  // Build DTC asociated to the cable
  std::ostringstream dtcNameStream;
  dtcNameStream << phiSectorRef << "_" << type << "_" << slot;
  std::string dtcName = dtcNameStream.str();

  DTC* dtc = GeometryFactory::make<DTC>(dtcName, phiSectorWidth, phiSectorRef, type, slot);
  dtc->addCable(this);
  
  myDTC_ = dtc;
  };



int Cable::computeServicesChannel(int phiSectorRef, std::string type, int slot) {
  int servicesChannel = 0;
  if (type == "PS10G") {
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
  else if (type == "PS5G") {
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
  else if (type == "2S") {
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

  return servicesChannel;
}
