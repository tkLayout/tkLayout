#include "ITCabling/ModulesToPowerChainsConnector.hh"
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
  //bundleType_ = computeBundleType(isBarrel_, barrelName_, layerNumber_); // Bundle cabling type
  rodPhi_ = r.Phi();
}
      
     
void ModulesToPowerChainsConnector::visit(BarrelModule& m) {
  // TO DO: this could be placed at the rod level
  const double modCenterX = m.center().X();
  const bool isPositiveXSide = computeXSide(modCenterX);

  const int halfNumRods = numRods_ / 2;
  const bool isPositiveZEnd = computeBarrelModuleZEnd(m.uniRef().side, m.uniRef().ring, rodPhi_, halfNumRods, isPositiveXSide);

  const int phiUnitRef = inner_cabling_functions::computePhiUnitRef(rodPhi_, halfNumRods, isPositiveZEnd);
  // PHIPOSITION.
  //const InnerPhiPosition& modulePhiPosition = InnerPhiPosition(rodPhi_, numRods_, isPositiveZEndTemp);


  // BUILD POWER CHAIN IF NECESSARY, AND CONNECT MODULE TO POWER CHAIN
  buildPowerChain(m, powerChains_, isPositiveZEnd, isPositiveXSide, barrelName_, layerNumber_, phiUnitRef);
}


void ModulesToPowerChainsConnector::visit(Endcap& e) {
  endcapName_ = e.myid();
}


void ModulesToPowerChainsConnector::visit(Disk& d) {
  diskNumber_ = d.diskNumber();
  endcapEnd_ = d.side();   // geometrical Z-side
}


void ModulesToPowerChainsConnector::visit(Ring& r)   { 
  ringNumber_ = r.myid();
  numModulesInRing_ = r.numModules();

  if (numModulesInRing_ % 4 != 0) logINFO(any2str("Found ") + any2str(numModulesInRing_)
				    + any2str(" as total number of modules in a forward ring.")
				    + any2str(" Total number of modules should be a multiple of 4.")
				    );

  //bundleType_ = computeBundleType(isBarrel_, endcapName_, diskNumber_, ringNumber_);
}


void ModulesToPowerChainsConnector::visit(EndcapModule& m) {
  const bool isPositiveZEnd = endcapEnd_;    // Alyways true in the Endcaps : cabling side and geometrical side are the same.

  const double modCenterX = m.center().X();
  const bool isPositiveXSide = computeXSide(modCenterX);
 
  // NOW THAT ALL INFORMATION HAS BEEN GATHERED, COMPUTE PHIPOSITION.
  //const PhiPosition& modulePhiPosition = PhiPosition(modPhi, numModulesInRing_, isBarrel_, diskNumber_, endcapName_, bundleType_);

  // BUILD BUNDLE IF NECESSARY, AND CONNECT MODULE TO BUNDLE
  //buildBundle(m, bundles_, negPowerChains_, bundleType_, isBarrel_, endcapName_, diskNumber_, modulePhiPosition, isPositiveCablingSide)

  const bool isRingInnerEnd = ( (m.diskSurface() % 2 ) == 1);
  
  const double modPhi = m.center().Phi();
  const int phiRef = computeForwardModulePhiPowerChain(modPhi, numModulesInRing_, isPositiveZEnd);

  const int ringQuarterIndex = inner_cabling_functions::computeRingQuarterIndex(ringNumber_, isRingInnerEnd);

  // BUILD POWER CHAIN IF NECESSARY, AND CONNECT MODULE TO POWER CHAIN
  buildPowerChain(m, powerChains_, isPositiveZEnd, isPositiveXSide, endcapName_, diskNumber_, phiRef, ringQuarterIndex);
}


void ModulesToPowerChainsConnector::postVisit() {
  // STAGGER MODULES
  //staggerModules(bundles_);
  //staggerModules(negPowerChains_);

  // CHECK
  checkModulesToPowerChainsCabling(powerChains_);
  //checkModulesToPowerChainsCabling(negPowerChains_);
}



const bool ModulesToPowerChainsConnector::computeXSide(const double modCenterX) const {   
  if (fabs(modCenterX) < inner_cabling_roundingTolerance) {
    logINFO(any2str("Found a module at X ~ 0. This is not handled.")
	    );
  }
  const bool isPositiveXSide = (modCenterX > 0.);    // Geometrical X-side (CMS frame of reference).

  return isPositiveXSide;
}


const bool ModulesToPowerChainsConnector::computeBarrelModuleZEnd(const int side, const int ring, const double rodPhi, const int numRods, const bool isPositiveXSide) const {
  bool isPositiveZEnd;

  // Non-central rings
  if (ring != 1) {
    isPositiveZEnd = (side > 0.);      // geometrical Z-side
  }
  // Central ring
  else {
    isPositiveZEnd = computeBarrelCentralModuleZEnd(rodPhi, numRods, isPositiveXSide);
  }

  return isPositiveZEnd;
}


/* Compute the Z end to which the central modules in BPIX are cabled.
   Modules are alternatively connected to the (+Z) end and the (-Z) end, depending on the Phi of the module.
*/
const bool ModulesToPowerChainsConnector::computeBarrelCentralModuleZEnd(const double rodPhi, const int numRods, const bool isPositiveXSide) const {
  // Compute the phi position of the module
  const bool isPositiveZEndTemp = true;
  const int phiUnitRefTemp = inner_cabling_functions::computePhiUnitRef(rodPhi, numRods, isPositiveZEndTemp);
  const bool isPhiUnitRefEven = ((phiUnitRefTemp % 2) == 0);

  const bool isPositiveZEnd = (isPositiveXSide ? isPhiUnitRefEven : !isPhiUnitRefEven);

  return isPositiveZEnd;
}


const int ModulesToPowerChainsConnector::computeForwardModulePhiPowerChain(const double modPhi, const int numModulesInRing, const bool isPositiveZEnd) const {
  int phiRef = 0;

  const int numModulesInRingEnd = numModulesInRing / 2;
  const int numModulesInRingQuarter = numModulesInRingEnd / 2;
  if (numModulesInRingQuarter > inner_cabling_maxNumModulesPerPowerChain) {
    const int phiUnitRef = inner_cabling_functions::computePhiUnitRef(modPhi, numModulesInRingEnd, isPositiveZEnd);
    const int numModulesInPowerChain = numModulesInRingQuarter / 2;
    if (phiUnitRef <= (numModulesInPowerChain - 1) ) { phiRef = 0; } // phiUnitRef starts numbering from 0
    else { phiRef = 1; }
  }
  return phiRef;
}


/* Build bundle.
 * The index associated to the bundleType is computed.
 * phiSliceRef: phiSegmentRef in Barrel, phiRegionRef in Endcaps.
 * This allows to compute the bundle Id.
 * Then, the bundle is created, and stored in the bundles_ or negPowerChains_ containers.
 * Lastly, each module is connected to its bundle, and vice-versa.
 */
void ModulesToPowerChainsConnector::buildPowerChain(DetectorModule& m, std::map<int, PowerChain*>& powerChains, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int ringQuarterIndex) {
  // COMPUTE POWER CHAIN ID
  const int powerChainId = computePowerChainId(isPositiveZEnd, isPositiveXSide, subDetectorName, layerDiskNumber, phiRef, ringQuarterIndex);

  // CREATE POWER CHAIN IF NECESSARY
  PowerChain* powerChain = nullptr;
  auto found = powerChains.find(powerChainId);
  if (found == powerChains.end()) {
    powerChain = createAndStorePowerChain(powerChains, powerChainId, isPositiveZEnd, isPositiveXSide, subDetectorName, layerDiskNumber, phiRef, ringQuarterIndex);
  }
  else {
    powerChain = found->second;
  }

  // CONNECT MODULE TO POWERCHAIN
  connectModuleToPowerChain(m, powerChain);
}


/* Compute the index associated to each bundle type.
 */
/*const int ModulesToPowerChainsConnector::computeBundleTypeIndex(const bool isBarrel, const Category& bundleType, const int totalNumFlatRings, const bool isTilted, const bool isExtraFlatPart) const {
  int bundleTypeIndex;
  // BARREL
  if (isBarrel) {
    if (bundleType == Category::SS) bundleTypeIndex = 0;
    else {
      if (!isTilted) {
	if (totalNumFlatRings <= cabling_maxNumModulesPerBundle) bundleTypeIndex = 1;
	else {
	  if (!isExtraFlatPart) bundleTypeIndex = 1;
	  else bundleTypeIndex = 2;
	}
      } 
      else bundleTypeIndex = 0;
    }
  }
  // ENDCAPS
  else {
    if (bundleType == Category::PS10GA) bundleTypeIndex = 0;
    else if (bundleType == Category::PS10GB) bundleTypeIndex = 1;
    else if (bundleType == Category::PS5G) bundleTypeIndex = 2;
    else if (bundleType == Category::SS) bundleTypeIndex = 3;
  }
  return bundleTypeIndex;
  }
*/


/* Compute the Id associated to each bundle.
 */
const int ModulesToPowerChainsConnector::computePowerChainId(const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int ringQuarterIndex) const {

  const int innerTrackerQuarterIndex = inner_cabling_functions::computeInnerTrackerQuarterIndex(isPositiveZEnd, isPositiveXSide);
  const int subdetectorIndex = inner_cabling_functions::computeSubDetectorIndex(subDetectorName);

  const int powerChainId = innerTrackerQuarterIndex * 10000 + subdetectorIndex * 1000 + layerDiskNumber * 100 + ringQuarterIndex * 10 + phiRef;
  return powerChainId;
}


/*  Create a powerChain, if it does not exist yet.
 *  Store it in the powerChains_ container.
 */
PowerChain* ModulesToPowerChainsConnector::createAndStorePowerChain(std::map<int, PowerChain*>& powerChains, const int powerChainId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int ringQuarterIndex) {

  PowerChain* powerChain = GeometryFactory::make<PowerChain>(powerChainId, isPositiveZEnd, isPositiveXSide, subDetectorName, layerDiskNumber, phiRef, ringQuarterIndex);

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
