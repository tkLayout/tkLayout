#include "Cabling/CablingMap.hh"
#include <Tracker.hh>


CablingMap::CablingMap(Tracker* tracker) {
  try {
    // CONNECT MODULES TO BUNDLES
    connectModulesToBundles(tracker);

    // CONNECT BUNDLES TO CABLES
    connectBundlesToCables(bundles_, cables_, DTCs_);
    connectBundlesToCables(negBundles_, negCables_, negDTCs_);

    assignBundlesStereoSemiBoundaries(bundles_, negBundles_);
    assignBundlesStereoSemiBoundaries(negBundles_, bundles_);

    // COMPUTE POWER SERVICES CHANNELS
    computePowerServicesChannels(bundles_, cables_);
    computePowerServicesChannels(negBundles_, negCables_);
  }

  catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }
}


/* MODULES TO BUNDLES CONNECTIONS.
 */
void CablingMap::connectModulesToBundles(Tracker* tracker) {
  ModulesToBundlesConnector bundlesBuilder;
  tracker->accept(bundlesBuilder);
  bundlesBuilder.postVisit();
  bundles_ = bundlesBuilder.getBundles();
  negBundles_ = bundlesBuilder.getNegBundles();
}


void CablingMap::assignBundlesStereoSemiBoundaries(std::map<int, Bundle*>& bundles, std::map<int, Bundle*>& complementaryBundles) {
  int phiSectorRefMarker = 0;
  int complementaryPhiSectorRefMarker = 0;

  for (auto& b : bundles) {
    Bundle* myBundle = b.second;
    const bool isBarrel = myBundle->isBarrel();
    const bool isPSFlatPart = myBundle->isPSFlatPart();

    if (isBarrel && !isPSFlatPart) {
      bool isLower;

      const int phiSectorRef = myBundle->getCable()->phiSectorRef();
      if (phiSectorRef != phiSectorRefMarker) {
	phiSectorRefMarker = phiSectorRef;
	complementaryPhiSectorRefMarker = -1;
	isLower = true;
      }

      const int complementaryBundleId = myBundle->complementaryBundleId();
      auto found = complementaryBundles.find(complementaryBundleId);
      if (found != complementaryBundles.end()) {
	Bundle* myComplementaryBundle = found->second;
	const int complementaryPhiSectorRef = myComplementaryBundle->getCable()->phiSectorRef();
	
	if (complementaryPhiSectorRefMarker != -1 && complementaryPhiSectorRefMarker != complementaryPhiSectorRef) isLower = false;
	myBundle->setIsInLowerSemiPhiSectorStereo(isLower);
	complementaryPhiSectorRefMarker = complementaryPhiSectorRef;
      }
      else { std::cout << "Coud not find complementary bundle id " << complementaryBundleId << std::endl; }
    }
  }


  for (auto& b : bundles) {
    Bundle* myBundle = b.second;
    const bool isBarrel = myBundle->isBarrel();
    const bool isPSFlatPart = myBundle->isPSFlatPart();

    if (isBarrel && isPSFlatPart) {
      const int bundleId = myBundle->myid();
      const int tiltedBundleId = bundleId - femod(bundleId, 10);
 
      auto found = bundles.find(tiltedBundleId);
      if (found != bundles.end()) {
	Bundle* myTiltedBundle = found->second;
	const bool isLower = myTiltedBundle->isInLowerSemiPhiSectorStereo();
	myBundle->setIsInLowerSemiPhiSectorStereo(isLower);
      }
      else { std::cout << "Coud not find tilted bundle id " << tiltedBundleId << std::endl; }
    }
  }


}



/* BUNDLES TO POWER SERVICE CHANNELS CONNECTIONS.
 */
void CablingMap::computePowerServicesChannels(std::map<int, Bundle*>& bundles, std::map<int, Cable*>& cables) {

  for (auto& c : cables) {
    c.second->assignPowerServicesChannels();

    /*
    const bool isPositiveCablingSide = b.second->isPositiveCablingSide();

    const double meanPhiOfficial = b.second->meanPhi();
    const double meanPhi = (isPositiveCablingSide ? meanPhiOfficial : (M_PI - meanPhiOfficial));
    const double semiPhiRegionStart = 0.;
    const int semiPhiRegionRef = computePhiSliceRef(meanPhi, semiPhiRegionStart, cabling_semiNonantWidth, true);

    std::pair<int, ChannelSection> powerServicesChannel = computePowerServicesChannel(semiPhiRegionRef, isPositiveCablingSide);
    b.second->setPowerServicesChannel(powerServicesChannel);*/
  }

  // CHECK POWER SERVICES CHANNELS
  checkBundlesToPowerServicesChannels(bundles);
}



/* BUNDLES TO CABLES CONNECTIONS.
 */
void CablingMap::connectBundlesToCables(std::map<int, Bundle*>& bundles, std::map<int, Cable*>& cables, std::map<const std::string, const DTC*>& DTCs) {

  for (auto& b : bundles) {
    // COLLECT ALL INFORMATION NEEDED TO BUILD CABLES
    const double phiSectorWidth = b.second->phiPosition().phiSectorWidth();

    const Category& bundleType = b.second->type();
    const Category& cableType = computeCableType(bundleType);

    const int bundleId = b.first;
    std::map<int, std::pair<int, int> > cablesPhiSectorRefAndSlot = computeCablesPhiSectorRefAndSlot(bundles);
    const int cablePhiSectorRef = cablesPhiSectorRefAndSlot.at(bundleId).first;
    const int slot = cablesPhiSectorRefAndSlot.at(bundleId).second;

    const int cableTypeIndex = computeCableTypeIndex(cableType);
    bool isPositiveCablingSide = b.second->isPositiveCablingSide();
    const int cableId = computeCableId(cablePhiSectorRef, cableTypeIndex, slot, isPositiveCablingSide);

    // BUILD CABLES and DTCS AND STORE THEM
    createAndStoreCablesAndDTCs(b.second, cables, DTCs, cableId, phiSectorWidth, cablePhiSectorRef, cableType, slot, isPositiveCablingSide);
  }

  // CHECK CABLES
  checkBundlesToCablesCabling(cables);
}


/* Compute cabling type associated to cable.
 */
const Category CablingMap::computeCableType(const Category& bundleType) const {
 Category cableType = bundleType;
 if (bundleType == Category::PS10GA || bundleType == Category::PS10GB) cableType = Category::PS10G;
 return cableType;
}


/* Compute phiSectorRef and slot.
 * Per phiSector, there are several DTCS.
 * 1 slot = 1 DTC.
 * A staggering of bundles is done on the fly, when the number of bundles per cable is too high.
 */
const std::map<int, std::pair<int, int> > CablingMap::computeCablesPhiSectorRefAndSlot(const std::map<int, Bundle*>& bundles) const {
  std::map<int, std::pair<int, int> > cablesPhiSectorRefAndSlot;

  // Used to stagger several bundles
  std::map<int, int> Layer3FlatPhiSectorsCounter;
  std::map<int, int> Layer3TiltedPhiSectorsCounter;
  std::map<int, int> Layer4PhiSectorsCounter;
  std::map<int, int> Layer5PhiSectorsCounter;
  std::map<int, int> Layer6PhiSectorsCounter;

  for (const auto& b : bundles) {
    const Bundle* myBundle = b.second;

    // COLLECT RELEVANT INFO
    const PhiPosition& bundlePhiPosition = myBundle->phiPosition();
    const int phiSegmentRef = bundlePhiPosition.phiSegmentRef(); 
    const double phiSectorWidth = bundlePhiPosition.phiSectorWidth();
    const int phiSectorRef = bundlePhiPosition.phiSectorRef();  
   
    const int numPhiSectors = round(2 * M_PI / phiSectorWidth);
    const int nextPhiSectorRef = computeNextPhiSliceRef(phiSectorRef, numPhiSectors);
    const int previousPhiSectorRef = computePreviousPhiSliceRef(phiSectorRef, numPhiSectors);

    const Category& bundleType = myBundle->type();
    const Category& cableType = computeCableType(bundleType);

    const std::string subDetectorName = myBundle->subDetectorName();
    const int layerDiskNumber = myBundle->layerDiskNumber();

    const bool isPositiveCablingSide = myBundle->isPositiveCablingSide();

    // by default
    int cablePhiSectorRef = phiSectorRef;
    int slot = 0;

    // PS10G
    if (cableType == Category::PS10G) {
      // BARREL FLAT PART + ENDCAPS DISKS 1, 3, 5
      if ((subDetectorName == cabling_tbps && !myBundle->isTiltedPart()) || (subDetectorName == cabling_tedd1 && layerDiskNumber == 1) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 3) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 5)) {
	slot = 1;
      }
      // BARREL TILTED PART + ENDCAPS DISKS 2, 4
      if ((subDetectorName == cabling_tbps && myBundle->isTiltedPart()) || (subDetectorName == cabling_tedd1 && layerDiskNumber == 2) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 4)) {
	slot = 2;
      }
    }

    // PS5G
    else if (cableType == Category::PS5G) {
      if (subDetectorName == cabling_tbps && layerDiskNumber == 2) {
	slot = 3;
      }

      else if ((subDetectorName == cabling_tedd1 && layerDiskNumber == 1) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 3) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 5)) {
	slot = 4;
      }

      // STAGGERING
      else if ( (subDetectorName == cabling_tbps && layerDiskNumber == 3) || (subDetectorName == cabling_tedd1 && layerDiskNumber == 2) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 4) ) {
	// TBPS
	if (subDetectorName == cabling_tbps) {
	  // TILTED PART
	  if (myBundle->isTiltedPart()) {
	    int& myPhiSectorCounter = Layer3TiltedPhiSectorsCounter[phiSectorRef];
	    myPhiSectorCounter += 1;
	    // In case already 4 bundles from tilted part, assign to next phi Sector
	    if (myPhiSectorCounter > 4) {
	      myPhiSectorCounter -= 1;
	      Layer3TiltedPhiSectorsCounter[nextPhiSectorRef] += 1;
	      cablePhiSectorRef = nextPhiSectorRef;
	    }
	    slot = 6;
	  }
	  // FLAT PART : assign TBPS bundles with TEDD bundles
	  else {
	    int& myPhiSectorCounter = Layer3FlatPhiSectorsCounter[phiSectorRef];
	    myPhiSectorCounter += 1;
	    // In case already 4 bundles from flat part, assign to next phi Sector
	    if (myPhiSectorCounter > 4) {
	      myPhiSectorCounter -= 1;
	      Layer3FlatPhiSectorsCounter[nextPhiSectorRef] += 1;
	      cablePhiSectorRef = nextPhiSectorRef;
	    }
	    slot = 5;
	  }
	}
	// TEDD_2
	else {
	  if (subDetectorName == cabling_tedd1 && layerDiskNumber == 2) slot = 5;
	  else if (subDetectorName == cabling_tedd2 && layerDiskNumber == 4) slot = 6;
	}
      }
    }

    // 2S
    else if (cableType == Category::SS) {
      int phiSectorRefThird = femod(phiSectorRef % 3, 3);

      // STAGGERING TB2S LAYER 4
      if (subDetectorName == cabling_tb2s && layerDiskNumber == 4) {
	int& myPhiSectorCounter = Layer4PhiSectorsCounter[phiSectorRef];
	myPhiSectorCounter += 1;
	// In a few cases, need to reduce to 5 bundles (additional bundle from Layer 5 will be added).
	// As a result, the first bundle in the Phi Sector is assigned to the previous phiSector.
	if (phiSectorRefThird == 0 && myPhiSectorCounter == 1) {
	  cablePhiSectorRef = previousPhiSectorRef;
	}
	slot = 1;
      }

      // STAGGERING TB2S LAYER 5
      else if (subDetectorName == cabling_tb2s && layerDiskNumber == 5) {
	int& myPhiSectorCounter = Layer5PhiSectorsCounter[phiSectorRef];
	myPhiSectorCounter += 1;
	// STAGGER BUNDLES : ASSIGN BUNDLES FROM LAYER 5 TO LAYER 4
	// In 2 cases out of 3, also one bundle from Layer 5 added.
	// if (phiSectorRefThird == 0) : should have 5 bundles in Layer 4 (+ 1 added from Layer 5)
	// if (phiSectorRefThird == 1) : should have 5 bundles in Layer 4 (+ 1 added from Layer 5)
	// if (phiSectorRefThird == 2) : should have 6 bundles in Layer 4 (no bundle added from Layer 5)
	if (phiSectorRefThird != 2 && myPhiSectorCounter == 4) {
	  slot = 1;
	}
	else {
	  slot = 2;
	}
      }

      else if ( (subDetectorName == cabling_tb2s && layerDiskNumber == 6) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 3) ) {
	// STAGGER BUNDLES : ASSIGN BUNDLES FROM LAYER 6 TO DISK 3
	if (subDetectorName == cabling_tb2s) {
	  int& myPhiSectorCounter = Layer6PhiSectorsCounter[phiSectorRef];
	  myPhiSectorCounter += 1;
	  if (myPhiSectorCounter == 1 || myPhiSectorCounter == 5 || myPhiSectorCounter == 8) slot = 4;
	  else slot = 3;
	}
	else slot = 4;
      }

      else if ( (subDetectorName == cabling_tedd1 && layerDiskNumber == 1) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 4) ) {
	slot = 5;
      }

      else if ( (subDetectorName == cabling_tedd1 && layerDiskNumber == 2) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 5) ) {
	slot = 6;
      }
    }
  
    if (slot == 0) {
      std::cout << "bundleType = "  << bundleType << " cableType = " << cableType <<  " subDetectorName  =" << subDetectorName << " layerDiskNumber = " << layerDiskNumber << " isPositiveCablingSide = " << isPositiveCablingSide << std::endl;
      logERROR("Connection from ribbon to cable : ribbon category is unknown. Slot was not defined properly.");
    }


    // COLLECT phiSectorRef and slots which have been computed.
    const std::pair<int, int> phiSectorRefAndSlot = std::make_pair(cablePhiSectorRef, slot);
    const int bundleId = b.first;
    cablesPhiSectorRefAndSlot.insert(std::make_pair(bundleId, phiSectorRefAndSlot));
  }

  return cablesPhiSectorRefAndSlot;
}


/* Create a cable and DTC, if do not exist yet.
 *  Store them in the cables or DTCs containers.
 */
void CablingMap::createAndStoreCablesAndDTCs(Bundle* myBundle, std::map<int, Cable*>& cables, std::map<const std::string, const DTC*>& DTCs, const int cableId, const double phiSectorWidth, const int cablePhiSectorRef, const Category& cableType, const int slot, const bool isPositiveCablingSide) {

  auto found = cables.find(cableId);
  if (found == cables.end()) {
    Cable* cable = GeometryFactory::make<Cable>(cableId, phiSectorWidth, cablePhiSectorRef, cableType, slot, isPositiveCablingSide);
    connectOneBundleToOneCable(myBundle, cable);
    cables.insert(std::make_pair(cableId, cable));

    const DTC* dtc = cable->getDTC();
    DTCs.insert(std::make_pair(dtc->name(), dtc));      
  }
  else {
    connectOneBundleToOneCable(myBundle, found->second);
  }
}


/* Compute cabling type associated to a cable.
 */
const int CablingMap::computeCableTypeIndex(const Category& cableType) const {
  int cableTypeIndex;
  if (cableType == Category::PS10G) cableTypeIndex = 0;
  else if (cableType == Category::PS5G) cableTypeIndex = 1;
  else if (cableType == Category::SS) cableTypeIndex = 2;
  return cableTypeIndex;
}


/* Compute Id associated to a cable.
 */
const int CablingMap::computeCableId(const int cablePhiSectorRef, const int cableTypeIndex, const int slot, const bool isPositiveCablingSide) const {
  int cablingSideIndex = (isPositiveCablingSide ? 0 : 1);

  const int cableId = cablingSideIndex * 1000 + cablePhiSectorRef * 100 + cableTypeIndex * 10 + slot;
  return cableId;
}


/* Connect bundle to cable and vice-versa.
 */
void CablingMap::connectOneBundleToOneCable(Bundle* bundle, Cable* cable) const {
  cable->addBundle(bundle);
  bundle->setCable(cable);
}


void CablingMap::checkBundlesToPowerServicesChannels(std::map<int, Bundle*>& bundles) {
  std::map<std::pair<const int, const ChannelSection >, int > channels;
  for (auto& b : bundles) {
    const int servicesChannel = b.second->powerServicesChannel();
    const ChannelSection servicesChannelSection = b.second->powerServicesChannelSection();

    if (fabs(servicesChannel) == 0 || fabs(servicesChannel) >= 13) std::cout << "ERROR: power servicesChannel = " << servicesChannel << std::endl;
    if (servicesChannelSection != ChannelSection::A && servicesChannelSection != ChannelSection::C) std::cout << "ERROR: power servicesChannelSection = " << servicesChannelSection << std::endl;
    std::pair<const int, const ChannelSection > myChannel = std::make_pair(servicesChannel, servicesChannelSection);
    channels[myChannel] += 1;
  }

  for (const auto& c : channels) { 
    if (c.second > 36) std::cout << "Power services channel " << c.first.first << " section " << c.first.second << " has " << c.second << " bundles." << std::endl;
  }
}


/* Check bundles-cables connections.
 */
void CablingMap::checkBundlesToCablesCabling(std::map<int, Cable*>& cables) {
  for (auto& c : cables) {

    // CHECK WHETHER THE PHI SLICES REF MAKE SENSE.
    const int phiSectorRef = c.second->phiSectorRef();
    if (phiSectorRef <= -1) {
      logERROR(any2str("Building cabling map : a cable was not correctly created. ")
	       + "Cable " + any2str(c.first) + ", with cableType = " + any2str(c.second->type())
	       + ", has phiSectorRef = " + any2str(phiSectorRef)
	       + ", slot = " + any2str(c.second->slot())
	       );
    }

    // CHECK THE NUMBER OF BUNDLES PER CABLE.
    const int numBundles = c.second->numBundles();
    if (numBundles > cabling_maxNumBundlesPerCable) {
      logERROR(any2str("Building cabling map : Staggering bundles. ")
	       + "Cable "  + any2str(c.first) + " is connected to " + any2str(numBundles) + " bundles."
	       + "Maximum number of bundles per cable allowed is " + any2str(cabling_maxNumBundlesPerCable)
	       );
    }

  }
}

