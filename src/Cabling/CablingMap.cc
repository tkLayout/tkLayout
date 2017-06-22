#include "Cabling/CablingMap.hh"
#include <Tracker.hh>


CablingMap::CablingMap(Tracker* tracker) {
  try {
    // MODULES TO BUNDLES
    connectModulesToBundles(tracker);

    // BUNDLES TO CABLES
    connectBundlesToCables(bundles_, cables_, DTCs_);
    connectBundlesToCables(negBundles_, negCables_, negDTCs_);
  }

  catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }
}


// MODULES TO BUNDLES
void CablingMap::connectModulesToBundles(Tracker* tracker) {
  ModulesToBundlesConnector bundlesBuilder;
  tracker->accept(bundlesBuilder);
  bundlesBuilder.postVisit();
  bundles_ = bundlesBuilder.getBundles();
  negBundles_ = bundlesBuilder.getNegBundles();
}


// BUNDLES TO CABLES
void CablingMap::connectBundlesToCables(std::map<int, Bundle*>& bundles, std::map<int, Cable*>& cables, std::map<const std::string, const DTC*>& DTCs) {

  for (auto& b : bundles) {
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


const Category CablingMap::computeCableType(const Category& bundleType) const {
 Category cableType = bundleType;
 if (bundleType == Category::PS5GA || bundleType == Category::PS5GB) cableType = Category::PS5G;
 return cableType;
}


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

    int cablePhiSectorRef = phiSectorRef;

    // Used to build cableId
    int slot = 0;

    if (cableType == Category::PS10G) {
      if (subDetectorName == cabling_tbps || (subDetectorName == cabling_tedd1 && layerDiskNumber == 1) || (subDetectorName == cabling_tedd1 && layerDiskNumber == 2)) {
	slot = 1;
      }
    }


    else if (cableType == Category::PS5G) {
      if ( (subDetectorName == cabling_tbps && layerDiskNumber == 2) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 3 && bundleType == Category::PS5GA) ) {
	slot = 2;
      }

      else if ( (subDetectorName == cabling_tbps && layerDiskNumber == 3) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 3 && bundleType == Category::PS5GB) ) {
	if (subDetectorName == cabling_tbps) {
	  // Tilted part
	  if (myBundle->isTiltedPart()) {
	    int& myPhiSectorCounter = Layer3TiltedPhiSectorsCounter[phiSectorRef];
	    myPhiSectorCounter += 1;
	    // In case already 4 bundles from tilted part, assign to next phi Sector
	    if (myPhiSectorCounter > 4) {
	      myPhiSectorCounter -= 1;
	      Layer3TiltedPhiSectorsCounter[nextPhiSectorRef] += 1;
	      cablePhiSectorRef = nextPhiSectorRef;
	    }
	    slot = 3;
	  }
	  // Flat part : assign TBPS bundles with TEDD bundles
	  else {
	    int& myPhiSectorCounter = Layer3FlatPhiSectorsCounter[phiSectorRef];
	    myPhiSectorCounter += 1;
	    // In case already 4 bundles from flat part, assign to next phi Sector
	    if (myPhiSectorCounter > 4) {
	      myPhiSectorCounter -= 1;
	      Layer3FlatPhiSectorsCounter[nextPhiSectorRef] += 1;
	      cablePhiSectorRef = nextPhiSectorRef;
	    }
	    slot = 4;
	  }
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


    else if (cableType == Category::SS) {
      int phiSectorRefThird = femod(phiSectorRef % 3, 3);

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
  
    if (slot == 0) logERROR("Connection from ribbon to cable : ribbon category is unknown. Slot was not defined properly.");

    const std::pair<int, int> phiSectorRefAndSlot = std::make_pair(cablePhiSectorRef, slot);
    const int bundleId = b.first;
    cablesPhiSectorRefAndSlot.insert(std::make_pair(bundleId, phiSectorRefAndSlot));
  }

  return cablesPhiSectorRefAndSlot;
}


void CablingMap::createAndStoreCablesAndDTCs(Bundle* myBundle, std::map<int, Cable*>& cables, std::map<const std::string, const DTC*>& DTCs, const int cableId, const double phiSectorWidth, const int cablePhiSectorRef, const Category& cableType, const int slot, const bool isPositiveCablingSide) {

  auto found = cables.find(cableId);
  if (found == cables.end()) {
    Cable* cable = GeometryFactory::make<Cable>(cableId, phiSectorWidth, cablePhiSectorRef, cableType, slot, isPositiveCablingSide);
    connectBundleToCable(myBundle, cable);
    cables.insert(std::make_pair(cableId, cable));

    const DTC* dtc = cable->getDTC();
    DTCs.insert(std::make_pair(dtc->name(), dtc));      
  }
  else {
    connectBundleToCable(myBundle, found->second);
  }
}


const int CablingMap::computeCableTypeIndex(const Category& cableType) const {
  int cableTypeIndex;
  if (cableType == Category::PS10G) cableTypeIndex = 0;
  else if (cableType == Category::PS5G) cableTypeIndex = 1;
  else if (cableType == Category::SS) cableTypeIndex = 2;
  return cableTypeIndex;
}


const int CablingMap::computeCableId(const int cablePhiSectorRef, const int cableTypeIndex, const int slot, const bool isPositiveCablingSide) const {
  int cablingSideIndex = (isPositiveCablingSide ? 0 : 1);

  const int cableId = cablingSideIndex * 1000 + cablePhiSectorRef * 100 + cableTypeIndex * 10 + slot;
  return cableId;
}


void CablingMap::connectBundleToCable(Bundle* bundle, Cable* cable) const {
  cable->addBundle(bundle);
  bundle->setCable(cable);
}


void CablingMap::checkBundlesToCablesCabling(std::map<int, Cable*>& cables) {
  for (auto& c : cables) {

    const int phiSectorRef = c.second->phiSectorRef();
    if (phiSectorRef <= -1) {
      logERROR(any2str("Building cabling map : a cable was not correctly created. ")
	       + "Cable " + any2str(c.first) + ", with cableType = " + any2str(c.second->type())
	       + ", has phiSectorRef = " + any2str(phiSectorRef)
	       + ", slot = " + any2str(c.second->slot())
	       );
    }

    const int numBundles = c.second->numBundles();
    if (numBundles > cabling_maxNumBundlesPerCable) {
      logERROR(any2str("Building cabling map : Staggering bundles. ")
	       + "Cable "  + any2str(c.first) + " is connected to " + any2str(numBundles) + " bundles."
	       + "Maximum number of bundles per cable allowed is " + any2str(cabling_maxNumBundlesPerCable)
	       );
    }

  }
}
