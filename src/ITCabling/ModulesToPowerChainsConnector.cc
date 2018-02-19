#include "ITCabling/ModulesToPowerChainsConnector.hh"
#include <Tracker.hh>


void ModulesToPowerChainsConnector::visit(Barrel& b) {
  //isBarrel_ = true;
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

  std::cout << "Barrel" << std::endl;
  std::cout << "layerNumber_ = " << layerNumber_ << std::endl;
  const int phiUnitRef = inner_cabling_functions::computePhiUnitRef(rodPhi_, halfNumRods, isPositiveZEnd);
  // PHIPOSITION.
  //const InnerPhiPosition& modulePhiPosition = InnerPhiPosition(rodPhi_, numRods_, isPositiveZEndTemp);


  // BUILD POWER CHAIN IF NECESSARY, AND CONNECT MODULE TO POWER CHAIN
  buildPowerChain(m, powerChains_, isPositiveZEnd, isPositiveXSide, barrelName_, layerNumber_, phiUnitRef);
}


void ModulesToPowerChainsConnector::visit(Endcap& e) {
  //isBarrel_ = false;
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
    std::cout << "Endcap" << std::endl;
    const int phiUnitRef = inner_cabling_functions::computePhiUnitRef(modPhi, numModulesInRingEnd, isPositiveZEnd);
    const int numModulesInPowerChain = numModulesInRingQuarter / 2;
    if (phiUnitRef <= (numModulesInPowerChain - 1) ) { phiRef = 0; } // phiUnitRef starts numbering from 0
    else { phiRef = 1; }
  }
  return phiRef;
}


/* Compute the bundle cabling type.
 */
/*const Category ModulesToPowerChainsConnector::computeBundleType(const bool isBarrel, const std::string subDetectorName, const int layerDiskNumber, const int ringNumber) const {
  Category bundleType = Category::UNDEFINED;

  //BARREL
  if (isBarrel) {
    // TB2S
    if (subDetectorName == cabling_tb2s) {
      bundleType = Category::SS;
    }
    // TBPS
    else if (subDetectorName == cabling_tbps) {
      bundleType = (layerDiskNumber == 1 ? Category::PS10G : Category::PS5G);
    }
  }

  // ENDCAPS
  else {
    // TEDD_1
    if (subDetectorName == cabling_tedd1) {
      if (ringNumber <= 4) bundleType = Category::PS10GA;
      else if (ringNumber >= 5 && ringNumber <= 7) bundleType = Category::PS10GB;
      else if (ringNumber >= 8 && ringNumber <= 10) bundleType = Category::PS5G;
      else if (ringNumber >= 11) bundleType = Category::SS;
    }

    // TEDD_2
    else if (subDetectorName == cabling_tedd2) {
      if (ringNumber <= 3) bundleType = Category::UNDEFINED;
      else if (ringNumber >= 4 && ringNumber <= 6) bundleType = Category::PS10GB;
      else if (ringNumber >= 7 && ringNumber <= 10) bundleType = Category::PS5G;
      else if (ringNumber >= 11) bundleType = Category::SS;
    }
  }

  return bundleType;
  }*/


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


/* Compute the Id associated to each stereo bundle.
 * The stereo bundle is the bundle located on the other cabling side, by rotation of 180° around CMS_Y.
 */
/*const int ModulesToPowerChainsConnector::computeStereoBundleId(const bool isBarrel, const bool isPositiveCablingSide, const int layerDiskNumber, const int stereoPhiRef, const int bundleTypeIndex) const {
  int stereoBundleId = 0;
  if (isBarrel) {
    const bool stereoSide = !isPositiveCablingSide;
    stereoBundleId = computeBundleId(isBarrel, stereoSide, layerDiskNumber, stereoPhiRef, bundleTypeIndex);
  }

  return stereoBundleId;
}
*/


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
 

/* Stagger modules.
   This is a very important step.
   In theory, each module is connected to the bundle corresponding to its position in Phi.
   Though, there can be more modules connected to a bundle than possible.
   If this is the case, one needs to stagger modules :
   the module closest to an adjacent phiRegion is removed and placed in the adjacent phiRegion.
*/
/*
void ModulesToPowerChainsConnector::staggerModules(std::map<int, Bundle*>& bundles) {

  for (auto& b : bundles) {
    // All this happens in Endcaps only
    if (b.second->subDetectorName() == cabling_tedd1 || b.second->subDetectorName() == cabling_tedd2) {
      const bool isBarrel = false;

      // Too many modules per bundle: staggering needed!
      while (b.second->numModules() > cabling_maxNumModulesPerBundle) {
	const bool isPositiveCablingSide = b.second->isPositiveCablingSide();
	const int diskNumber = b.second->layerDiskNumber();

	const Category& bundleType = b.second->type();
	const int bundleTypeIndex = computeBundleTypeIndex(isBarrel, bundleType);

	const PhiPosition& bundlePhiPosition = b.second->phiPosition();

	const double phiRegionStart = bundlePhiPosition.phiRegionStart();

	const double phiRegionWidth = bundlePhiPosition.phiRegionWidth();
	const int numPhiRegions = round(2 * M_PI / phiRegionWidth);

	// Compute the ref of the phiRegion, and of the adjacent phi regions.
	// 'next' or 'previous' refers to:
	// - positiveCablingSide: positive Phi order.
	// - negativeCablingSide: negative Phi order.
	const int phiRegionRef = bundlePhiPosition.phiRegionRef();
	const int nextPhiRegionRef = computeNextPhiSliceRef(phiRegionRef, numPhiRegions);
	const int previousPhiRegionRef = computePreviousPhiSliceRef(phiRegionRef, numPhiRegions);

	// Compute the associated bundles ids (so that the associated bundles can be accessed).
	const int bundleId = b.first;
	const int nextBundleId = computeBundleId(isBarrel, isPositiveCablingSide, diskNumber, nextPhiRegionRef, bundleTypeIndex);
	const int previousBundleId = computeBundleId(isBarrel, isPositiveCablingSide, diskNumber, previousPhiRegionRef, bundleTypeIndex);

	// Distance in Phi from the smallest Phi module to the smallest Phi boundary of the phiRegion.
	const double minPhiBorder = fabs( femod((b.second->minPhi() - phiRegionStart), phiRegionWidth) );
	// Distance in Phi from the greatest Phi module to the greatest Phi boundary of the phiRegion.
	const double maxPhiBorder = fabs( femod((b.second->maxPhi() - phiRegionStart), phiRegionWidth) - phiRegionWidth);
	      

	// Access bundles.
	auto previousBundleSearch = bundles.find(previousBundleId);
	auto nextBundleSearch = bundles.find(nextBundleId);
	if (previousBundleSearch != bundles.end() && nextBundleSearch != bundles.end()) {
	  Bundle* previousBundle = previousBundleSearch->second;
	  Bundle* nextBundle = nextBundleSearch->second;

	  const int previousBundleNumModules = previousBundle->numModules();
	  const int nextBundleNumModules = nextBundle->numModules();
	  // A MODULE FROM THE PHI REGION NEEDS TO BE MOVED. LOOK WHERE TO MOVE IT!

	  // Cannot assign the extra module : both neighbouring phi regions are full !
	  if (previousBundleNumModules >= cabling_maxNumModulesPerBundle && nextBundleNumModules >= cabling_maxNumModulesPerBundle) {
	    logERROR(any2str("Building cabling map : Staggering modules.")
		     + "I am a module in side " + any2str(isPositiveCablingSide)
		     + ", disk " + any2str(diskNumber)
		     + ", bundleType " + any2str(bundleType) 
		     + ", phiRegionRef " + any2str(phiRegionRef)
		     + ", phiRegionWidth " + any2str(phiRegionWidth)
		     + ". My phiRegion has more than " + any2str(cabling_maxNumModulesPerBundle) + " modules."
		     + " My 2 neighbouring phiRegions have also more than " + any2str(cabling_maxNumModulesPerBundle) + " modules."
		     + " Hence, I cannot be staggered to any other phiRegion :/ "
		     );
	    break;
	  }

	  // Assign a module to the next phi region.
	  else if (previousBundleNumModules >= cabling_maxNumModulesPerBundle || maxPhiBorder <= minPhiBorder) {
	    logINFO(any2str("Building cabling map : Staggering modules.")
		    + " I am a module in side " + any2str(isPositiveCablingSide)
		    + ", disk " + any2str(diskNumber)
		    + ", bundleType " + any2str(bundleType)
		    + ", phiRegionRef " + any2str(phiRegionRef)
		    + ". maxPhiBorder " + any2str((maxPhiBorder * 180. / M_PI))
		    + ". My region has " + any2str(b.second->numModules())
		    + " > maxNumModulesPerBundle = " + any2str(cabling_maxNumModulesPerBundle)
		    + ". I am moved to the next phiRegion, which presently has " + any2str(nextBundleNumModules) + " modules."
		    );
		    
	    Module* maxPhiMod = b.second->maxPhiModule();
	    maxPhiMod->setBundle(nextBundle);  
	    nextBundle->moveMaxPhiModuleFromOtherBundle(b.second);
	  }

	  // Assign a module to the previous phi region.
	  else if (nextBundleNumModules >= cabling_maxNumModulesPerBundle || minPhiBorder < maxPhiBorder) {
	    logINFO(any2str("Building cabling map : Staggering modules.")
		    + " I am a module in side " + any2str(isPositiveCablingSide)
		    + ", disk " + any2str(diskNumber)
		    + ", bundleType " + any2str(bundleType)
		    + ", phiRegionRef " + any2str(phiRegionRef)
		    + ". minPhiBorder " + any2str((minPhiBorder * 180. / M_PI))
		    + ". My region has " + any2str(b.second->numModules())
		    + " > maxNumModulesPerBundle = " + any2str(cabling_maxNumModulesPerBundle)
		    + ". I am moved to the previous phiRegion, which presently has " + any2str(previousBundleNumModules) + " modules."
		    );

	    Module* minPhiMod = b.second->minPhiModule();
	    minPhiMod->setBundle(previousBundle);	  
	    previousBundle->moveMinPhiModuleFromOtherBundle(b.second);
	  }
	}
	else { // PowerChains not found when trying to access them.
	  logERROR(any2str("Building cabling map : Staggering modules.")
		   + "Error building previousBundleId or nextBundleId"); 
	  break; 
	}

      }
    }
  }

}
*/


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
