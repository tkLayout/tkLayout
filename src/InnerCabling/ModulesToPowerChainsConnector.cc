#include "InnerCabling/ModulesToPowerChainsConnector.hh"
#include <Tracker.hh>


void ModulesToPowerChainsConnector::visit(Barrel& b) {
  barrelName_ = b.myid();	
}


void ModulesToPowerChainsConnector::visit(Layer& l) {
  layerNumber_ = l.layerNumber();
  numRods_ = l.numRods();

  if (numRods_ % 4 != 0) logINFO(any2str("Found ") + any2str(numRods_)
				    + any2str(" as total number of rods in a barrel layer.")
				    + any2str(" Total number of rods should be a multiple of 4.")
				    );
}


void ModulesToPowerChainsConnector::visit(RodPair& r) {
  rodPhi_ = r.Phi();
}
 
     
/*
 * Gather all module info, build Power Chain, and assign module to it.
 */     
void ModulesToPowerChainsConnector::visit(BarrelModule& m) {
  // TO DO: this could be placed at the rod level
  const double modCenterX = m.center().X();
  const bool isPositiveXSide = computeXSide(modCenterX);

  const int halfNumRods = numRods_ / 2;
  const std::pair<bool, bool>& barrelModuleZEnd = computeBarrelModuleZEnd(m.uniRef().side, m.uniRef().ring, layerNumber_);
  const bool isPositiveZEnd = barrelModuleZEnd.first;
  const bool isLongBarrel = barrelModuleZEnd.second;

  const int phiUnitRef = inner_cabling_functions::computePhiUnitRef(rodPhi_, halfNumRods, isPositiveZEnd);
  const int modulePhiRefInPowerChain = femod(inner_cabling_functions::computePhiUnitRef(rodPhi_, numRods_, isPositiveZEnd), 2);
  m.setPhiRefInPowerChain(modulePhiRefInPowerChain);

  // BUILD POWER CHAIN IF NECESSARY, AND CONNECT MODULE TO POWER CHAIN
  buildPowerChain(m, powerChains_, isPositiveZEnd, isPositiveXSide, barrelName_, layerNumber_, phiUnitRef, isLongBarrel);
}


/*
 * Gather all module info, build Power Chain, and assign module to it.
 */   
void ModulesToPowerChainsConnector::visit(Endcap& e) {
  endcapName_ = e.myid();
}


void ModulesToPowerChainsConnector::visit(Disk& d) {
  diskNumber_ = d.myid();
  endcapEnd_ = d.side();   // geometrical Z-side
}


void ModulesToPowerChainsConnector::visit(Ring& r)   { 
  ringNumber_ = r.myid();
  numModulesInRing_ = r.numModules();

  if (numModulesInRing_ % 4 != 0) logINFO(any2str("Found ") + any2str(numModulesInRing_)
				    + any2str(" as total number of modules in a forward ring.")
				    + any2str(" Total number of modules should be a multiple of 4.")
				    );
}


void ModulesToPowerChainsConnector::visit(EndcapModule& m) {
  const bool isPositiveZEnd = endcapEnd_;    // Alyways true in the Endcaps : cabling side and geometrical side are the same.

  const double modCenterX = m.center().X();
  const bool isPositiveXSide = computeXSide(modCenterX);

  const bool isSmallerAbsZHalfRing = m.isSmallerAbsZModuleInRing();
  
  const double modPhi = m.center().Phi();
  const std::pair<int, int> phiRefs = computeForwardModulePhiPowerChain(modPhi, numModulesInRing_, isPositiveZEnd);
  const int powerChainPhiRef = phiRefs.first;
  const int modulePhiRefInPowerChain = phiRefs.second;
  m.setPhiRefInPowerChain(modulePhiRefInPowerChain);

  const int halfRingIndex = inner_cabling_functions::computeHalfRingIndex(ringNumber_, isSmallerAbsZHalfRing);
  const bool isAtSmallerAbsZDeeInDoubleDisk = m.isAtSmallerAbsZDeeInDoubleDisk();
  const bool isAtSmallerAbsZSideInDee = m.isAtSmallerAbsZSideInDee();

  // BUILD POWER CHAIN IF NECESSARY, AND CONNECT MODULE TO POWER CHAIN
  const bool isLongBarrel = false;
  buildPowerChain(m, powerChains_, isPositiveZEnd, isPositiveXSide, endcapName_, diskNumber_, powerChainPhiRef, isLongBarrel, halfRingIndex, isAtSmallerAbsZDeeInDoubleDisk, isAtSmallerAbsZSideInDee);
}


void ModulesToPowerChainsConnector::postVisit() {
  // CHECK
  checkModulesToPowerChainsCabling(powerChains_);
}



/*
 * Computes the (X) side on which a module cabling must be located.
 */
const bool ModulesToPowerChainsConnector::computeXSide(const double modCenterX) const {   
  if (fabs(modCenterX) < inner_cabling_roundingTolerance) {
    logINFO(any2str("Found a module at X ~ 0. This is not handled.")
	    );
  }
  const bool isPositiveXSide = (modCenterX > 0.);    // Geometrical X-side (CMS frame of reference).

  return isPositiveXSide;
}


/*
 * Computes the (Z) end on which a module cabling must be located.
 */
const std::pair<bool, bool> ModulesToPowerChainsConnector::computeBarrelModuleZEnd(const int side, const int ring, const int layerNumber) const {
  const bool isBarrelCentralModuleAtPositiveZEnd = computeBarrelCentralModuleZEnd(layerNumber);

  // Compute the module's (Z) end
  bool isPositiveZEnd;
  // Non-central rings
  if (ring != 1) {
    isPositiveZEnd = (side > 0.);      // geometrical Z-side
  }
  // Central ring
  else {
    isPositiveZEnd = isBarrelCentralModuleAtPositiveZEnd;
  }

  // The module is on the long barrel (Z) end <-> the central module is connected to the same (Z) end.
  const bool isLongBarrel = (isPositiveZEnd == isBarrelCentralModuleAtPositiveZEnd);

  return std::make_pair(isPositiveZEnd, isLongBarrel);
}


/* Compute the (Z) end to which the central modules in BPIX are cabled.
   Modules are alternatively connected to the (+Z) end and the (-Z) end, depending on the layer of the module.
*/
const bool ModulesToPowerChainsConnector::computeBarrelCentralModuleZEnd(const int layerNumber) const {
  const bool isOddLayer = ((layerNumber % 2) == 1);

  const bool isPositiveZEnd = isOddLayer;

  return isPositiveZEnd;
}


/*
 * Assign a Forward Module to a power chain based on its phi position in the ring.
 */
const std::pair<int, int> ModulesToPowerChainsConnector::computeForwardModulePhiPowerChain(const double modPhi, const int numModulesInRing, const bool isPositiveZEnd) const {
  int powerChainPhiRef = 0;

  const int numModulesInRingEnd = numModulesInRing / 2;
  const int numModulesInRingQuarter = numModulesInRingEnd / 2;
  const int phiUnitRef = inner_cabling_functions::computePhiUnitRef(modPhi, numModulesInRingEnd, isPositiveZEnd);
  int modulePhiRefInPowerChain = phiUnitRef;

  if (numModulesInRingQuarter > inner_cabling_maxNumModulesPerPowerChain) {
    const int numModulesInPowerChain = numModulesInRingQuarter / 2;
    if (phiUnitRef <= (numModulesInPowerChain - 1) ) { // powerChainPhiRef starts numbering from 0
      powerChainPhiRef = 0;
      modulePhiRefInPowerChain = phiUnitRef; 
    } 
    else { 
      powerChainPhiRef = 1; 
      modulePhiRefInPowerChain = phiUnitRef - numModulesInPowerChain; 
    }
  }
  return std::make_pair(powerChainPhiRef, modulePhiRefInPowerChain);
}


/* Build bundle.
 * The index associated to the bundleType is computed.
 * phiSliceRef: phiSegmentRef in Barrel, phiRegionRef in Endcaps.
 * This allows to compute the bundle Id.
 * Then, the bundle is created, and stored in the bundles_ or negPowerChains_ containers.
 * Lastly, each module is connected to its bundle, and vice-versa.
 */
void ModulesToPowerChainsConnector::buildPowerChain(DetectorModule& m, std::map<int, PowerChain*>& powerChains, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const bool isLongBarrel, const int halfRingIndex, const bool isAtSmallerAbsZDeeInDoubleDisk, const bool isAtSmallerAbsZSideInDee) {
  // COMPUTE POWER CHAIN ID
  const int powerChainId = computePowerChainId(isPositiveZEnd, isPositiveXSide, subDetectorName, layerDiskNumber, phiRef, halfRingIndex);

  // CREATE POWER CHAIN IF NECESSARY
  PowerChain* powerChain = nullptr;
  auto found = powerChains.find(powerChainId);
  if (found == powerChains.end()) {
    powerChain = createAndStorePowerChain(powerChains, powerChainId, isPositiveZEnd, isPositiveXSide, subDetectorName, layerDiskNumber, phiRef, isLongBarrel, halfRingIndex, isAtSmallerAbsZDeeInDoubleDisk, isAtSmallerAbsZSideInDee);
  }
  else {
    powerChain = found->second;
  }

  // CONNECT MODULE TO POWERCHAIN
  connectModuleToPowerChain(m, powerChain);
}


/* Compute the Id associated to each bundle.
 */
const int ModulesToPowerChainsConnector::computePowerChainId(const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int halfRingIndex) const {

  const int innerTrackerQuarterIndex = inner_cabling_functions::computeInnerTrackerQuarterIndex(isPositiveZEnd, isPositiveXSide);
  const int subdetectorIndex = inner_cabling_functions::computeSubDetectorIndex(subDetectorName);

  const int powerChainId = innerTrackerQuarterIndex * 10000 + subdetectorIndex * 1000 + layerDiskNumber * 100 + halfRingIndex * 10 + phiRef;
  return powerChainId;
}


/*  Create a powerChain, if it does not exist yet.
 *  Store it in the powerChains_ container.
 */
PowerChain* ModulesToPowerChainsConnector::createAndStorePowerChain(std::map<int, PowerChain*>& powerChains, const int powerChainId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const bool isLongBarrel, const int halfRingIndex, const bool isAtSmallerAbsZSideInDee) {

  PowerChain* powerChain = new PowerChain(powerChainId, isPositiveZEnd, isPositiveXSide, subDetectorName, layerDiskNumber, phiRef, isLongBarrel, halfRingIndex, isAtSmallerAbsZDeeInDoubleDisk, isAtSmallerAbsZSideInDee);

  powerChains.insert(std::make_pair(powerChainId, powerChain));
 
  return powerChain;
}


/* Connect module to powerChain and vice-versa.
 */
void ModulesToPowerChainsConnector::connectModuleToPowerChain(DetectorModule& m, PowerChain* powerChain) const {
  powerChain->addModule(&m);
  m.setPowerChain(powerChain);
}


/* Check modules-powerChains connections.
 */
void ModulesToPowerChainsConnector::checkModulesToPowerChainsCabling(const std::map<int, PowerChain*>& powerChains) const {
  for (auto& b : powerChains) {

    // CHECK THE NUMBER OF MODULES PER POWER CHAIN.
    const int powerChainNumModules = b.second->numModules();
    if (powerChainNumModules > inner_cabling_maxNumModulesPerPowerChain) {
      logERROR(any2str("Building IT cabling map: ")
	       + "PowerChain "  + any2str(b.first) + " is connected to " + any2str(powerChainNumModules) + " modules."
	       + "Maximum number of modules per powerChain allowed is " + any2str(inner_cabling_maxNumModulesPerPowerChain)
	       );
    }

  }
}
