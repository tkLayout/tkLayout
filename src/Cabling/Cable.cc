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
  // ASSIGN AN OPTICAL SERVICESCHANNEL TO THE CABLE
  ServicesChannel* opticalChannel_ = GeometryFactory::make<ServicesChannel>(phiSectorRef, type, slot, isPositiveCablingSide);

  // BUILD DTC ASOCIATED TO THE CABLE
  buildDTC(phiSectorWidth, phiSectorRef, type, slot, isPositiveCablingSide);  
};


Cable::~Cable() {
  delete myDTC_;       // TO DO: switch to smart pointers and remove this!
  myDTC_ = nullptr;
}


void Cable::assignPowerServicesChannels() {
  for (auto& myBundle : bundles_) {

    const double meanPhi = femod(myBundle.meanPhi(), 2.*M_PI);
    const int cablePhiSectorRef = phiSectorRef_;

    const bool isBarrel = myBundle.isBarrel();
    const Category& type = type_;

    bool isLower;
    if (isBarrel) { isLower = myBundle.isInLowerSemiPhiSectorStereo(); }
    else {
      const double phiShift = ((type == Category::SS) ? 
				cabling_powerChannelsTeddStripStripSemiNonantBoundaryShift 
				: cabling_powerChannelsTeddPixelStripSemiNonantBoundaryShift);
      const double semiNonantBoundary = cablePhiSectorRef * cabling_nonantWidth + cabling_semiNonantWidth + phiShift;
      isLower = moduloComp(meanPhi, semiNonantBoundary, 2.*M_PI);
    }

    const int semiPhiRegionIndex = (isLower ? 0 : 1);
    const int semiPhiRegionRef = 2 * cablePhiSectorRef + semiPhiRegionIndex;

    ServicesChannel* powerChannel = GeometryFactory::make<ServicesChannel>(semiPhiRegionRef, isPositiveCablingSide_);
    myBundle.setPowerServicesChannel(powerChannel);
  }
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
