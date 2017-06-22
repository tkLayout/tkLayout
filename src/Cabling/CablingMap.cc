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

  // Used to stagger several bundles
  std::map<int, int> Layer3FlatPhiSectorsCounter;
  std::map<int, int> Layer3TiltedPhiSectorsCounter;
  std::map<int, int> Layer4PhiSectorsCounter;
  std::map<int, int> Layer5PhiSectorsCounter;
  std::map<int, int> Layer6PhiSectorsCounter;

  for (auto& b : bundles) {
    const PhiPosition& bundlePhiPosition = b.second->phiPosition();

    const int phiSectorRef = bundlePhiPosition.phiSectorRef();  
    const double phiSectorWidth = bundlePhiPosition.phiSectorWidth();
    const int numPhiSectors = round(2 * M_PI / phiSectorWidth);

    const int nextPhiSectorRef = computeNextPhiSliceRef(phiSectorRef, numPhiSectors);
    const int previousPhiSectorRef = computePreviousPhiSliceRef(phiSectorRef, numPhiSectors);

    const Category& bundleType = b.second->type();
    Category cableType = bundleType;
    if (bundleType == Category::PS5GA || bundleType == Category::PS5GB) cableType = Category::PS5G;

    int cableTypeIndex;
    if (cableType == Category::PS10G) cableTypeIndex = 0;
    else if (cableType == Category::PS5G) cableTypeIndex = 1;
    else if (cableType == Category::SS) cableTypeIndex = 2;


    const std::string subDetectorName = b.second->subDetectorName();
    const int layerDiskNumber = b.second->layerDiskNumber();

    int phiSectorRefCable = phiSectorRef;

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
	  if (b.second->isTiltedPart()) {
	    Layer3TiltedPhiSectorsCounter[phiSectorRef] += 1;
	    // In case already 4 bundles from tilted part, assign to next phi Sector
	    auto found = Layer3TiltedPhiSectorsCounter.find(phiSectorRef);
	    if (found.second > 4) {
	      found.second -= 1;
	      Layer3TiltedPhiSectorsCounter[nextPhiSectorRef] += 1;
	      phiSectorRefCable = nextPhiSectorRef;
	    }
	    slot = 3;
	  }
	  // Flat part : assign TBPS bundles with TEDD bundles
	  else {
	    Layer3FlatPhiSectorsCounter[phiSectorRef] += 1;
	    auto found = Layer3FlatPhiSectorsCounter.find(phiSectorRef);
	    // In case already 4 bundles from flat part, assign to next phi Sector
	    if (found.second) > 4) {
	      found.second -= 1;
	      Layer3FlatPhiSectorsCounter[nextPhiSectorRef] += 1;
	      phiSectorRefCable = nextPhiSectorRef;
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
	Layer4PhiSectorsCounter[phiSectorRef] += 1;
	// In a few cases, need to reduce to 5 bundles (additional bundle from Layer 5 will be added).
	// As a result, the first bundle in the Phi Sector is assigned to the previous phiSector.
	if (phiSectorRefThird == 0 && Layer4PhiSectorsCounter[phiSectorRef] == 1) {
	  phiSectorRefCable = previousPhiSectorRef;
	}
	slot = 1;
      }

      else if (subDetectorName == cabling_tb2s && layerDiskNumber == 5) {
	Layer5PhiSectorsCounter[phiSectorRef] += 1;
	// STAGGER BUNDLES : ASSIGN BUNDLES FROM LAYER 5 TO LAYER 4
	// In 2 cases out of 3, also one bundle from Layer 5 added.
	// if (phiSectorRefThird == 0) : should have 5 bundles in Layer 4 (+ 1 added from Layer 5)
	// if (phiSectorRefThird == 1) : should have 5 bundles in Layer 4 (+ 1 added from Layer 5)
	// if (phiSectorRefThird == 2) : should have 6 bundles in Layer 4 (no bundle added from Layer 5)
	if (phiSectorRefThird != 2 && Layer5PhiSectorsCounter[phiSectorRef] == 4) {
	  slot = 1;
	}
	else {
	  slot = 2;
	}
      }

      else if ( (subDetectorName == cabling_tb2s && layerDiskNumber == 6) || (subDetectorName == cabling_tedd2 && layerDiskNumber == 3) ) {
	// STAGGER BUNDLES : ASSIGN BUNDLES FROM LAYER 6 TO DISK 3
	if (subDetectorName == cabling_tb2s) {
	  Layer6PhiSectorsCounter[phiSectorRef] += 1;
	  if (Layer6PhiSectorsCounter[phiSectorRef] == 1 || Layer6PhiSectorsCounter[phiSectorRef] == 5 || Layer6PhiSectorsCounter[phiSectorRef] == 8) slot = 4;
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

    // BUILD CABLE AND STORE IT
    int cableId = phiSectorRefCable * 100 + cableTypeIndex * 10 + slot;

    bool isPositiveCablingSide = b.second->isPositiveCablingSide();
    if (!isPositiveCablingSide) cableId += 1000;

    if (cables.count(cableId) == 0) {
      Cable* cable = GeometryFactory::make<Cable>(cableId, phiSectorWidth, phiSectorRefCable, cableType, slot, isPositiveCablingSide);
      cable->addBundle(b.second);
      b.second->setCable(cable);
      cables.insert(std::make_pair(cableId, cable));

      const DTC* dtc = cable->getDTC();
      DTCs.insert(std::make_pair(dtc->name(), dtc));      
    }
    else {
      cables[cableId]->addBundle(b.second);
      b.second->setCable(cables[cableId]);
    }
  }


  // CHECK CABLES
  checkBundlesToCablesCabling(cables);
}


void CablingMap::checkBundlesToCablesCabling(std::map<int, Cable*>& cables) {
  for (auto& c : cables) {

    const int phiSectorRef = c.second->phiSectorRef();
    if (phiSectorRef <= -1) {
      logERROR(any2str("Building cabling map : a cable was not correctly created. ")
	       + "Cable " + any2str(c.first) + ", with cableType = " + any2str(c.second->type())
	       + ", has phiSectorRef = " + any2str(phiSectorRef)
	       + ", slot = " << any2str(c.second->slot())
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
