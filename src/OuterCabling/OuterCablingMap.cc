#include "OuterCabling/OuterCablingMap.hh"
#include <Tracker.hh>


OuterCablingMap::OuterCablingMap(Tracker* tracker) {
  try {
    // CONNECT MODULES TO BUNDLES
    connectModulesToBundles(tracker);

    // CONNECT BUNDLES TO CABLES
    connectBundlesToCables(bundles_, cables_, DTCs_);
    connectBundlesToCables(negBundles_, negCables_, negDTCs_);

    computeCMSSWIds(DTCs_);
    computeCMSSWIds(negDTCs_);

    // COMPUTE SERVICES CHANNELS ASSIGNMENTS OF POWER CABLES
    computePowerServicesChannels();
  }

  catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }
}



/* MODULES TO BUNDLES CONNECTIONS.
 */
void OuterCablingMap::connectModulesToBundles(Tracker* tracker) {
  // Build bundles
  ModulesToBundlesConnector bundlesBuilder;
  tracker->accept(bundlesBuilder);
  bundlesBuilder.postVisit();

  // Store bundles
  for (const auto& bundleIt : bundlesBuilder.getBundles()) {
    std::unique_ptr<OuterBundle> myBundle(bundleIt.second);
    bundles_.insert(std::make_pair(bundleIt.first, std::move(myBundle)));
  }
  for (const auto& bundleIt : bundlesBuilder.getNegBundles()) {
    std::unique_ptr<OuterBundle> myBundle(bundleIt.second);
    negBundles_.insert(std::make_pair(bundleIt.first, std::move(myBundle)));
  }

  //Store GBTs
  for (const auto& gbtIt : bundlesBuilder.getGBTs()){
   std::unique_ptr<OuterGBT> myGBT(gbtIt.second);
   gbts_.insert(std::make_pair(gbtIt.first, std::move(myGBT)));
  }
}


/* BUNDLES TO CABLES CONNECTIONS.
 */
void OuterCablingMap::connectBundlesToCables(std::map<const int, std::unique_ptr<OuterBundle> >& bundles, std::map<const int, std::unique_ptr<OuterCable> >& cables, std::map<const std::string, std::unique_ptr<const OuterDTC> >& DTCs) {

  for (auto& b : bundles) {
    // COLLECT ALL INFORMATION NEEDED TO BUILD CABLES
    OuterBundle* myBundle = b.second.get();
    const double phiSectorWidth = myBundle->phiPosition().phiSectorWidth();

    const Category& bundleType = myBundle->type();
    const Category& cableType = computeCableType(bundleType);

    const int bundleId = b.first;
    std::map<int, std::pair<int, int> > cablesPhiSectorRefAndSlot = computeCablesPhiSectorRefAndSlot(bundles);
    const int cablePhiSectorRef = cablesPhiSectorRefAndSlot.at(bundleId).first;
    const int slot = cablesPhiSectorRefAndSlot.at(bundleId).second;

    const int cableTypeIndex = computeCableTypeIndex(cableType);
    bool isPositiveCablingSide = myBundle->isPositiveCablingSide();
    const int cableId = computeCableId(cablePhiSectorRef, cableTypeIndex, slot, isPositiveCablingSide);

    // BUILD CABLES and DTCS AND STORE THEM
    createAndStoreCablesAndDTCs(myBundle, cables, DTCs, cableId, phiSectorWidth, cablePhiSectorRef, cableType, slot, isPositiveCablingSide);
  }

  // CHECK CABLES
  checkBundlesToCablesCabling(cables);
}


/* Compute cabling type associated to cable.
 */
const Category OuterCablingMap::computeCableType(const Category& bundleType) const {
 Category cableType = bundleType;
 if (bundleType == Category::PS10GA || bundleType == Category::PS10GB) cableType = Category::PS10G;
 return cableType;
}


/* Compute phiSectorRef and slot.
 * Per phiSector, there are several DTCS.
 * 1 slot = 1 DTC.
 * A staggering of bundles is done on the fly, when the number of bundles per cable is too high.
 */
const std::map<int, std::pair<int, int> > OuterCablingMap::computeCablesPhiSectorRefAndSlot(const std::map<const int, std::unique_ptr<OuterBundle> >& bundles) const {
  std::map<int, std::pair<int, int> > cablesPhiSectorRefAndSlot;

  // Used to stagger several bundles
  std::map<int, int> Layer3FlatPhiSectorsCounter;
  std::map<int, int> Layer3TiltedPhiSectorsCounter;
  std::map<int, int> Layer4PhiSectorsCounter;
  std::map<int, int> Layer5PhiSectorsCounter;
  std::map<int, int> Layer6PhiSectorsCounter;

  for (const auto& b : bundles) {
    const OuterBundle* myBundle = b.second.get();

    // COLLECT RELEVANT INFO
    const PhiPosition& bundlePhiPosition = myBundle->phiPosition();
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
      // BARREL LAYER 1 FLAT PART + ENDCAPS DISKS 1, 3, 5
      if ( (subDetectorName == outer_cabling_tbps && layerDiskNumber == 1 && !myBundle->isTiltedPart()) 
	   || (subDetectorName == outer_cabling_tedd1 && layerDiskNumber == 1)
	   || (subDetectorName == outer_cabling_tedd2 && layerDiskNumber == 3)	   
	   || (subDetectorName == outer_cabling_tedd2 && layerDiskNumber == 5)
	   ) {
	slot = 1;
      }
      // BARREL LAYER 1 TILTED PART + ENDCAPS DISKS 2B, 4
      else if ( (subDetectorName == outer_cabling_tbps && layerDiskNumber == 1 && myBundle->isTiltedPart()) 
		|| (subDetectorName == outer_cabling_tedd1 && layerDiskNumber == 2 && bundleType == Category::PS10GB)
		|| (subDetectorName == outer_cabling_tedd2 && layerDiskNumber == 4)
		) {
	slot = 2;
      }
      // BARREL FULL LAYER 2 + ENDCAPS DISK 2A
      else if ( (subDetectorName == outer_cabling_tbps && layerDiskNumber == 2) 
		|| (subDetectorName == outer_cabling_tedd1 && layerDiskNumber == 2 && bundleType == Category::PS10GA)
		) {
	slot = 3;
      }
    }

    // PS5G
    else if (cableType == Category::PS5G) {
      if ( (subDetectorName == outer_cabling_tedd1 && layerDiskNumber == 1) 
	   || (subDetectorName == outer_cabling_tedd2 && layerDiskNumber == 3) 
	   || (subDetectorName == outer_cabling_tedd2 && layerDiskNumber == 5)
	   ) {
	slot = 4;
      }

      // STAGGERING
      else if ( (subDetectorName == outer_cabling_tbps && layerDiskNumber == 3) 
		|| (subDetectorName == outer_cabling_tedd1 && layerDiskNumber == 2) 
		|| (subDetectorName == outer_cabling_tedd2 && layerDiskNumber == 4)
		) {
	// TBPS
	if (subDetectorName == outer_cabling_tbps) {
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
	  if (subDetectorName == outer_cabling_tedd1 && layerDiskNumber == 2) slot = 5;
	  else if (subDetectorName == outer_cabling_tedd2 && layerDiskNumber == 4) slot = 6;
	}
      }
    }

    // 2S
    else if (cableType == Category::SS) {
      int phiSectorRefThird = femod(phiSectorRef % 3, 3);

      // STAGGERING TB2S LAYER 4
      if (subDetectorName == outer_cabling_tb2s && layerDiskNumber == 4) {
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
      else if (subDetectorName == outer_cabling_tb2s && layerDiskNumber == 5) {
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

      else if ( (subDetectorName == outer_cabling_tb2s && layerDiskNumber == 6) 
		|| (subDetectorName == outer_cabling_tedd2 && layerDiskNumber == 3) 
		) {
	// STAGGER BUNDLES : ASSIGN BUNDLES FROM LAYER 6 TO DISK 3
	if (subDetectorName == outer_cabling_tb2s) {
	  int& myPhiSectorCounter = Layer6PhiSectorsCounter[phiSectorRef];
	  myPhiSectorCounter += 1;
	  if (myPhiSectorCounter == 1 || myPhiSectorCounter == 5 || myPhiSectorCounter == 8) slot = 4;
	  else slot = 3;
	}
	else slot = 4;
      }

      else if ( (subDetectorName == outer_cabling_tedd1 && layerDiskNumber == 1) 
		|| (subDetectorName == outer_cabling_tedd2 && layerDiskNumber == 4)
		) {
	slot = 5;
      }

      else if ( (subDetectorName == outer_cabling_tedd1 && layerDiskNumber == 2) 
		|| (subDetectorName == outer_cabling_tedd2 && layerDiskNumber == 5)
		) {
	slot = 6;
      }
    }
  
    if (slot == 0) {
      logERROR(any2str("Connection from ribbon to cable : ribbon category is unknown. Slot was not defined properly.")
	       + " bundleType = " + any2str(bundleType) 
	       + ", cableType = " + any2str(cableType) 
	       + ", subDetectorName  =" + any2str(subDetectorName) 
	       + ", layerDiskNumber = " + any2str(layerDiskNumber) 
	       + ", isPositiveCablingSide = " + any2str(isPositiveCablingSide)
	       );
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
void OuterCablingMap::createAndStoreCablesAndDTCs(OuterBundle* myBundle, std::map<const int, std::unique_ptr<OuterCable> >& cables, std::map<const std::string, std::unique_ptr<const OuterDTC> >& DTCs, const int cableId, const double phiSectorWidth, const int cablePhiSectorRef, const Category& cableType, const int slot, const bool isPositiveCablingSide) {

  auto found = cables.find(cableId);
  if (found == cables.end()) {
    std::unique_ptr<OuterCable> cable (new OuterCable(cableId, phiSectorWidth, cablePhiSectorRef, cableType, slot, isPositiveCablingSide));
    connectOneBundleToOneCable(myBundle, cable.get());
    std::unique_ptr<const OuterDTC> dtc (cable.get()->getDTC());

    cables.insert(std::make_pair(cableId, std::move(cable))); 
    DTCs.insert(std::make_pair(dtc->name(), std::move(dtc)));      
  }
  else {
    connectOneBundleToOneCable(myBundle, found->second.get());
  }
}


/* Compute cabling type associated to a cable.
 */
const int OuterCablingMap::computeCableTypeIndex(const Category& cableType) const {
  int cableTypeIndex = -1;
  if (cableType == Category::PS10G) cableTypeIndex = 0;
  else if (cableType == Category::PS5G) cableTypeIndex = 1;
  else if (cableType == Category::SS) cableTypeIndex = 2;
  if (cableTypeIndex == -1) logWARNING("Unexpected cableType " + any2str(cableType) + ": assigning an invalid cableTypeIndex = -1");
  return cableTypeIndex;
}


/* Compute Id associated to a cable.
 */
const int OuterCablingMap::computeCableId(const int cablePhiSectorRef, const int cableTypeIndex, const int slot, const bool isPositiveCablingSide) const {
  int cablingSideIndex = (isPositiveCablingSide ? 0 : 1);

  const int cableId = cablingSideIndex * 1000 + cablePhiSectorRef * 100 + cableTypeIndex * 10 + slot;
  return cableId;
}


/* Connect bundle to cable and vice-versa.
 */
void OuterCablingMap::connectOneBundleToOneCable(OuterBundle* bundle, OuterCable* cable) const {
  cable->addBundle(bundle);
  bundle->setCable(cable);
}


/* Check bundles-cables connections.
 */
void OuterCablingMap::checkBundlesToCablesCabling(const std::map<const int, std::unique_ptr<OuterCable> >& cables) {
  for (const auto& c : cables) {
    const OuterCable* myCable = c.second.get();

    // CHECK THE CABLES TYPES, PHI SLICES REF AND SLOTS.
    const Category& type = myCable->type();
    const int phiSectorRef = myCable->phiSectorRef();
    const int slot = myCable->slot();
    if (type == Category::UNDEFINED
	|| phiSectorRef <= -1 || phiSectorRef >= outer_cabling_numNonants
	|| (slot <= 0) || (slot >= (outer_cabling_maxNumDTCsPerNonantPerZEnd / 2 + 1))
	) {
      logERROR(any2str("Building cabling map : a cable was not correctly created. ")
	       + "OuterCable " + any2str(c.first) + ", with cableType = " + any2str(type)
	       + ", has phiSectorRef = " + any2str(phiSectorRef)
	       + ", slot = " + any2str(slot)
	       );
    }

    // CHECK THE NUMBER OF BUNDLES PER CABLE.
    const int numBundles = myCable->numBundles();
    if (numBundles > outer_cabling_maxNumBundlesPerCable) {
      logERROR(any2str("Building cabling map : Staggering bundles. ")
	       + "OuterCable "  + any2str(c.first) + " is connected to " + any2str(numBundles) + " bundles."
	       + "Maximum number of bundles per cable allowed is " + any2str(outer_cabling_maxNumBundlesPerCable)
	       );
    }

  }
}


/* BUNDLES TO POWER SERVICE CHANNELS CONNECTIONS.
 * VERY IMPORTANT: connection scheme from modules to optical bundles = connection scheme from modules to power cables.
 * As a result, 1 single Bundle object is used for both schemes.
 * Regarding the connections to services channels, each Bundle is then assigned:
 * - 1 Optical Services Channel Section (considering the Bundle as an optical Bundle);
 * - 1 Power Services Channels section (making as if the Bundle is a power cable);
 * 
 * The optical channel mapping is done so that all bundles connected to the same DTC,
 * are routed through the same channel.
 *
 * The assignments of power cables to services channels must be done after the bundles to DTCs connections are established.
 * Indeed, the power channel mapping is done so that all modules connected to the same DTC, 
 * have their power cables routed through 2 consecutive channels sections at most.
 */
void OuterCablingMap::computePowerServicesChannels() {
  for (bool isPositiveCablingSide : { true, false }) {

    // BARREL ONLY: IN VIEW OF THE POWER CHANNEL ASSIGNMENT, SPLIT EACH NONANT INTO 2 SEMI-NONANTS
    routeBarrelBundlesPoweringToSemiNonants(isPositiveCablingSide);
 
    // ASSIGN POWER SERVICES CHANNELS
    std::map<const int, std::unique_ptr<OuterCable> >& cables = (isPositiveCablingSide ? cables_ : negCables_);
    for (auto& c : cables) {
      c.second.get()->assignPowerChannelSections();
    }

    // CHECK POWER SERVICES CHANNELS
    const std::map<const int, std::unique_ptr<OuterBundle> >& bundles = (isPositiveCablingSide ? bundles_ : negBundles_);
    checkBundlesToPowerServicesChannels(bundles);
  }
}


/* Barrel only: in view of the power channel assignment, split each nonant into 2 semi-nonants.
   The cooling pipes design should be invariant by rotation of 180° around CMS_Y, to avoid different cooling designs on both (Z) side.
   The cooling pipes and power cables are assigned to similar channels slots (A or C).
   As a result, the channel assignmennt of power cables need to follow the same symmetry as the cooling pipes.
   Hence, the CHANNEL ASSIGNMENNT OF POWER CABLES NEED TO BE INVARIANT BY ROTATION OF 180DEG AROUND CMS_Y.
   This is a priori not trivial, since the power cables scheme follow the bundles scheme, hence the optical map mirror symmetry.
   This is only possible if, for each phi nonant, one define a semi-nonant Phi boundary in a certain way.
   This is what is done here: the semi-nonant Phi boundaries are defined
   so that ALL NONANTS AND SEMI-NONANTS PHI BOUNDARIES ARE INVARIANT BY ROTATION OF 180° AROUND CMS_Y.
 */
void OuterCablingMap::routeBarrelBundlesPoweringToSemiNonants(const bool isPositiveCablingSide) {
  std::map<const int, std::unique_ptr<OuterBundle> >& bundles = (isPositiveCablingSide ? bundles_ : negBundles_);
  const std::map<const int, std::unique_ptr<OuterBundle> >& stereoBundles = (isPositiveCablingSide ? negBundles_ : bundles_);

  // phiSectorRefMarker keeps track of the Phi nonant we are in.
  int phiSectorRefMarker = -1;
  // phiSectorRefMarker keeps track of the stereo Phi nonant we are in.
  // 'stereo' means on the other cabling side, by a rotation of 180° around CMS_Y.
  int stereoPhiSectorRefMarker = -1;

  // TILTED PART OF TBPS + TB2S
  // Loop on all bundles of a given cabling side (sorted by their BundleIds).
  for (auto& b : bundles) {
    OuterBundle* myBundle = b.second.get();
    const bool isBarrel = myBundle->isBarrel();
    const bool isBarrelPSFlatPart = myBundle->isBarrelPSFlatPart();

    // Only for Barrel: tilted TBPS, or TB2S.
    if (isBarrel && !isBarrelPSFlatPart) {
      // Should the bundle be assigned to the lower or upper semi-nonant ?
      // 'lower' and 'upper' are defined by 'smaller' or 'bigger' Phi, 
      // in the trigonometric sense in the (XY) plane in CMS global frame of reference.
      bool isLower = false; // what we want to compute!
      bool isLowerAssigned = false;

      // Identifier of the Phi nonant we are in.
      const int phiSectorRef = myBundle->getCable()->phiSectorRef();
      // In case of a switch to a different Phi nonant, initialize variables.
      if (phiSectorRef != phiSectorRefMarker) {	
	phiSectorRefMarker = phiSectorRef;
	// Starts by assigning to bundle to the lower semi-nonant.
	isLower = true;
	isLowerAssigned = true;
	stereoPhiSectorRefMarker = -1;	
      }

      // Get the bundle located on the other cabling side, by a rotation of 180° around CMS_Y.
      const int stereoBundleId = myBundle->stereoBundleId();
      auto found = stereoBundles.find(stereoBundleId);
      if (found != stereoBundles.end()) {
	const OuterBundle* myStereoBundle = found->second.get();
	// Get the Phi nonant in which the stereoBundle is located.
	const int stereoPhiSectorRef = myStereoBundle->getCable()->phiSectorRef();
	
	// Decisive point!! 
	// As soon as a change in the identifier of the stereoBundle Phi nonant is detected,
	// one assigns the bundle to the upper semi-nonant.
	if (stereoPhiSectorRefMarker != -1 && stereoPhiSectorRefMarker != stereoPhiSectorRef) {
		isLower = false;
		isLowerAssigned = true;
	}

	// Lastly, assign the semi-nonant attribution decision to the bundle.
	
	if (!isLowerAssigned) logERROR("I did not manage to assign the 'isLower' variable for this cable. I pick a random one just to carry on (false)");
	myBundle->setIsPowerRoutedToBarrelLowerSemiNonant(isLower);

	// Keeps track of the Phi nonant in which the stereoBundle is located.
	stereoPhiSectorRefMarker = stereoPhiSectorRef;
      }
      else {
	logERROR(any2str("Could not find stereo bundle, id ") 
		    + any2str(stereoBundleId)
		    );
      }
    }

    // NB: All this is possible because the bundles and stereoBundles are sorted by their Ids.
    // One relies on 2 characteristics of the BundleId scheme:
    // - all Bundles connected to the same DTC will have consecutive Ids.
    // - the Id increment is in Phi, in the (XY) plane in CMS global frame of reference.
  }


  // FLAT PART OF TBPS
  // For a given bundle, connected to untilted modules:
  // Take the same semi-nonant assignment as the bundle located at the same Phi and connected to the tilted modules.

  // Loop on all bundles of a given cabling side (sorted by their BundleIds).
  for (auto& b : bundles) {
    OuterBundle* myBundle = b.second.get();
    const bool isBarrelPSFlatPart = myBundle->isBarrelPSFlatPart();

    // Only for Barrel: flat part of TBPS
    if (isBarrelPSFlatPart) {
      const int tiltedBundleId = myBundle->tiltedBundleId();
 
      // Get the bundle located at the same Phi, but connected to the tilted modules.
      auto found = bundles.find(tiltedBundleId);
      if (found != bundles.end()) {
	const OuterBundle* myTiltedBundle = found->second.get();

	// Decisive point!!
	// Get the semi-nonant attribution of the tiltedBundle.
	const bool isLower = myTiltedBundle->isPowerRoutedToBarrelLowerSemiNonant();

	// Lastly, assign the semi-nonant attribution decision to the bundle.
	myBundle->setIsPowerRoutedToBarrelLowerSemiNonant(isLower);
      }
      else {
	logERROR(any2str("Could not find bundle connected to tilted modules, id ") 
		 + any2str(tiltedBundleId)
		 );
      }
    }
  }

}


/* Check services channels sections containing power cables.
 */
void OuterCablingMap::checkBundlesToPowerServicesChannels(const std::map<const int, std::unique_ptr<OuterBundle> >& bundles) {
  std::map<std::pair<const int, const ChannelSlot>, int > channels;

  for (const auto& b : bundles) {
    const OuterBundle* myBundle = b.second.get();

    const ChannelSection* mySection = myBundle->powerChannelSection();
    const int myChannelNumber = mySection->channelNumber();
    const ChannelSlot& myChannelSlot = mySection->channelSlot();

    // Check channel number.
    if (fabs(myChannelNumber) == 0 || fabs(myChannelNumber) > outer_cabling_numServicesChannels) {
      logERROR(any2str("Invalid channel number = ")
	       + any2str(myChannelNumber)
	       + any2str(". Should not be 0 or have abs value > ")
	       + any2str(outer_cabling_numServicesChannels) + any2str(".")
	       );
    }
    // Check channel slot.
    if (myChannelSlot != ChannelSlot::A && myChannelSlot != ChannelSlot::C) {
      logERROR(any2str("Power cable: Invalid channel slot = ") 
	       + any2str(myChannelSlot)
	       + any2str(". Should be ") + any2str(ChannelSlot::A)
	       + any2str(" or ") + any2str(ChannelSlot::C) + any2str(".")
	       );
    }

    // Compute number of power cables per channel section.
    std::pair<const int, const ChannelSlot> myChannel = std::make_pair(myChannelNumber, myChannelSlot);
    channels[myChannel] += 1;
  }

  // Check number of power cables per channel section.
  for (const auto& c : channels) { 
    if (c.second > outer_cabling_maxNumPowerCablesPerChannel) {
      logERROR(any2str("Services channel ") + any2str(c.first.first) 
	       + any2str(" section ") + any2str(c.first.second) 
	       + any2str(" has " ) + any2str(c.second) + any2str(" power cables.")
	       + any2str(" Max number of power cables per channel section is ") 
	       + any2str(outer_cabling_maxNumPowerCablesPerChannel)
	       );
    }
  }
}

/* Compute GBTs CMSSW Ids.
 * Want to have these as consectuive integers, ordered by detID, so done after the full map is created.
 */
void OuterCablingMap::computeCMSSWIds(std::map<const std::string, std::unique_ptr<const OuterDTC> >& DTCs) {
  std::string lastDTCId = "";

  // DTCs loop
  for (auto& it : DTCs) {
    const std::string myDTCId = it.first;
    const OuterDTC* myDTC = it.second.get();

    if (myDTCId != lastDTCId) {
      lastDTCId = myDTCId;
    }

    std::vector<int> detIds;
    std::map<int, int> detIdToGBT;
    int gbtCMSSWId=0;
    const std::vector<OuterCable*>& myCables = myDTC->cable();
    for (const auto& myCable : myCables){
       const std::vector<OuterBundle*>& myBundles = myCable->bundles();
       for (const auto& myBundle : myBundles){
         const std::vector<Module*> & myModules = myBundle->modules();
         for (const auto& myModule : myModules){
             detIds.push_back(myModule->myDetId());
         }
      }
    }
    //Sort detIDs and determine GBT IDs
    std::sort(detIds.begin(), detIds.end());
    for(auto id : detIds){
      detIdToGBT[id] = gbtCMSSWId;
      gbtCMSSWId++;
    }
   
    //Now loop again, extracting the GBT from the module and setting the GBT ID just stored in the detIdToGBT map.
    for (const auto& myCable : myCables){
       const std::vector<OuterBundle*>& myBundles = myCable->bundles();
       for (const auto& myBundle : myBundles){
         const std::vector<Module*> & myModules = myBundle->modules();
         for (const auto& myModule : myModules){
            OuterGBT* myGBT = myModule->getOuterGBT();
            myGBT->setCMSSWId(detIdToGBT[myModule->myDetId()]);
         }
      }
    }
  }
}

