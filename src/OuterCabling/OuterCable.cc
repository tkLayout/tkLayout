#include "OuterCabling/OuterCable.hh"
#include "OuterCabling/OuterDTC.hh"

using namespace insur;

OuterCable::OuterCable(const int id, const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) : 
  phiSectorWidth_(phiSectorWidth),
  phiSectorRef_(phiSectorRef),
  type_(type),
  slot_(slot),
  isPositiveCablingSide_(isPositiveCablingSide)
{ 
  myid(id);
  // ASSIGN AN OPTICAL SERVICESCHANNEL TO THE CABLE
  std::unique_ptr<const ChannelSection> myOpticalChannelSection (new OpticalSection(phiSectorRef, type, slot, isPositiveCablingSide));
  opticalChannelSection_ = std::move(myOpticalChannelSection);

  // BUILD DTC ASOCIATED TO THE CABLE
  buildDTC(phiSectorWidth, phiSectorRef, type, slot, isPositiveCablingSide);  
};


/* Build DTC asociated to the cable.
 */
void OuterCable::buildDTC(const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) {
  std::string dtcName = computeDTCName(phiSectorRef, type, slot, isPositiveCablingSide);
  OuterDTC* dtc = new OuterDTC(dtcName, phiSectorWidth, phiSectorRef, type, slot, isPositiveCablingSide);
  dtc->addCable(this);
  myDTC_ = dtc;
}


/* Compute DTC name.
 */
const std::string OuterCable::computeDTCName(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const {
  std::ostringstream dtcNameStream;
  if (!isPositiveCablingSide) dtcNameStream << outer_cabling_negativePrefix << "_";
  dtcNameStream << phiSectorRef << "_" << any2str(type) << "_" << slot;
  const std::string dtcName = dtcNameStream.str();
  return dtcName;
}


/* POWER SERVICES CHANNELS CONNECTIONS.
 * One use the equivalence 1 Bundle = 1 Power cable.
 * The power channel mapping is done so that all modules connected to the same DTC, 
 * have their power cables routed through 2 consecutive channels sections at most.
 */
void OuterCable::assignPowerChannelSections() {
  // Loop on all bundles connected to the same DTC.
  for (auto& myBundle : bundles_) {
    
    const bool isBarrel = myBundle->isBarrel();
    const int cablePhiSectorRef = phiSectorRef_;
    
    // COMPUTE ISLOWER
    // THIS DECIDES WETHER THE BUNDLE (= POWER CABLE) IS ASSIGNED TO LOWER OR UPPER SEMI-NONANT.
    // 'lower' and 'upper' are defined by 'smaller' or 'bigger' Phi, 
    // in the trigonometric sense in the (XY) plane in CMS global frame of reference.
    bool isLower;

    // BARREL: Use complex boundary relying on rotation of 180° around CMS_Y.
    if (isBarrel) { 
      isLower = myBundle->isPowerRoutedToBarrelLowerSemiNonant(); 
    }

    // ENDCAP: Define a Phi Semi-Nonant boundary and split bundles accordingly.
    else {
      // Average Phi of all modules the Bundle is connected to.
      const double meanPhi = femod(myBundle->meanPhi(), 2.*M_PI);

      const Category& type = type_;

      // Compute Phi Semi-Nonant boundary
      const double phiShift = ((type == Category::SS) ? 
				outer_cabling_powerChannelsTeddStripStripSemiNonantBoundaryShift 
				: outer_cabling_powerChannelsTeddPixelStripSemiNonantBoundaryShift);
      const double semiNonantBoundary = cablePhiSectorRef * outer_cabling_nonantWidth + outer_cabling_semiNonantWidth + phiShift;
      isLower = moduloComp(meanPhi, semiNonantBoundary, 2.*M_PI);
    }

    // COMPUTE AT WHICH SEMINONANT THE BUNDLE IS LOCATED.
    const int semiPhiRegionIndex = (isLower ? 0 : 1);
    const int semiPhiRegionRef = 2 * cablePhiSectorRef + semiPhiRegionIndex;

    // COMPUTE THE POWER CHANEL SECTION CORRESPONDING TO THE PHI SEMINONANT.
    std::unique_ptr<const ChannelSection> powerChannelSection (new PowerSection(semiPhiRegionRef, isPositiveCablingSide_));

    // ASSIGN THE POWER CHANNEL SECTION TO THE BUNDLE (bundle = power cable)
    myBundle->setPowerChannelSection(std::move(powerChannelSection));
  }
}
