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

    const int nextPhiSectorRef = femod( (phiSectorRef + 1), numPhiSectors);
    const int previousPhiSectorRef = femod( (phiSectorRef - 1), numPhiSectors);

    const std::string bundleType = b.second->type();
    std::string cableType = bundleType;
    if (cableType == "PS5GA" || cableType == "PS5GB") cableType = "PS5G";

    int cableTypeIndex;
    if (cableType == "PS10G") cableTypeIndex = 0;
    else if (cableType == "PS5G") cableTypeIndex = 1;
    else if (cableType == "2S") cableTypeIndex = 2;


    const std::string subDetectorName = b.second->subDetectorName();
    const int layerDiskNumber = b.second->layerDiskNumber();

    int phiSectorRefCable = phiSectorRef;

    // Used to build cableId
    int slot = 0;

    if (cableType == "PS10G") {
      if (subDetectorName == "TBPS" || (subDetectorName == "TEDD_1" && layerDiskNumber == 1) || (subDetectorName == "TEDD_1" && layerDiskNumber == 2)) {
	slot = 1;
      }
    }


    else if (cableType == "PS5G") {
      if ( (subDetectorName == "TBPS" && layerDiskNumber == 2) || (subDetectorName == "TEDD_2" && layerDiskNumber == 3 && bundleType == "PS5GA") ) {
	slot = 2;
      }

      else if ( (subDetectorName == "TBPS" && layerDiskNumber == 3) || (subDetectorName == "TEDD_2" && layerDiskNumber == 3 && bundleType == "PS5GB") ) {
	if (subDetectorName == "TBPS") {
	  // Tilted part
	  if (b.second->isTiltedPart()) {
	    Layer3TiltedPhiSectorsCounter[phiSectorRef] += 1;
	    // In case already 4 bundles from tilted part, assign to next phi Sector
	    if (Layer3TiltedPhiSectorsCounter.at(phiSectorRef) > 4) {
	      Layer3TiltedPhiSectorsCounter[phiSectorRef] -= 1;
	      Layer3TiltedPhiSectorsCounter[nextPhiSectorRef] += 1;
	      phiSectorRefCable = nextPhiSectorRef;
	    }
	    slot = 3;
	  }
	  // Flat part : assign TBPS bundles with TEDD bundles
	  else {
	    Layer3FlatPhiSectorsCounter[phiSectorRef] += 1;
	    // In case already 4 bundles from flat part, assign to next phi Sector
	    if (Layer3FlatPhiSectorsCounter.at(phiSectorRef) > 4) {
	      Layer3FlatPhiSectorsCounter[phiSectorRef] -= 1;
	      Layer3FlatPhiSectorsCounter[nextPhiSectorRef] += 1;
	      phiSectorRefCable = nextPhiSectorRef;
	    }
	    slot = 4;
	  }
	}
	else slot = 4;
      }

      else if ( (subDetectorName == "TEDD_1" && layerDiskNumber == 1) || (subDetectorName == "TEDD_2" && layerDiskNumber == 4) ) {
	slot = 5;
      }

      else if ( (subDetectorName == "TEDD_1" && layerDiskNumber == 2) || (subDetectorName == "TEDD_2" && layerDiskNumber == 5) ) {
	slot = 6;
      }
    }


    else if (cableType == "2S") {
      int phiSectorRefThird = femod(phiSectorRef % 3, 3);

      if (subDetectorName == "TB2S" && layerDiskNumber == 4) {
	Layer4PhiSectorsCounter[phiSectorRef] += 1;
	// In a few cases, need to reduce to 5 bundles (additional bundle from Layer 5 will be added).
	// As a result, the first bundle in the Phi Sector is assigned to the previous phiSector.
	if (phiSectorRefThird == 0 && Layer4PhiSectorsCounter[phiSectorRef] == 1) {
	  phiSectorRefCable = previousPhiSectorRef;
	}
	slot = 1;
      }

      else if (subDetectorName == "TB2S" && layerDiskNumber == 5) {
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

      else if ( (subDetectorName == "TB2S" && layerDiskNumber == 6) || (subDetectorName == "TEDD_2" && layerDiskNumber == 3) ) {
	// STAGGER BUNDLES : ASSIGN BUNDLES FROM LAYER 6 TO DISK 3
	if (subDetectorName == "TB2S") {
	  Layer6PhiSectorsCounter[phiSectorRef] += 1;
	  if (Layer6PhiSectorsCounter[phiSectorRef] == 1 || Layer6PhiSectorsCounter[phiSectorRef] == 5 || Layer6PhiSectorsCounter[phiSectorRef] == 8) slot = 4;
	  else slot = 3;
	}
	else slot = 4;
      }

      else if ( (subDetectorName == "TEDD_1" && layerDiskNumber == 1) || (subDetectorName == "TEDD_2" && layerDiskNumber == 4) ) {
	slot = 5;
      }

      else if ( (subDetectorName == "TEDD_1" && layerDiskNumber == 2) || (subDetectorName == "TEDD_2" && layerDiskNumber == 5) ) {
	slot = 6;
      }
    }
  
    if (slot == 0) std::cout << "Connection from ribbon to cable : ribbon category is unknown. Slot was not defined properly." << std::endl;

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
    if (c.second->numBundles() > maxNumBundlesPerCable_) {
      std::cout << "There was an error while staggering bundles. Cable " 
		<< c.first << " is connected to " << c.second->numBundles() << " bundles." 
		<< std::endl;
    }

    if (c.second->phiSectorRef() <= -1) {
      std::cout << "Error while creating cable. Cable " << c.first << " has phiSectorRef = " << c.second->phiSectorRef() << ". type = " << c.second->type() << ", slot = " <<  c.second->slot() << std::endl;
    }

  }
}
