#include "Cabling/ModulesToBundlesConnector.hh"
#include <Tracker.hh>


void ModulesToBundlesConnector::visit(Barrel& b) {
  barrelName = b.myid();	
}

void ModulesToBundlesConnector::visit(Layer& l) {
  layerNumber = l.layerNumber();
  numRods = l.numRods();
  totalNumFlatRings = l.buildNumModulesFlat() * 2 - 1; // Total number of flat rings on both +Z side and -Z side

  if (layerNumber == 1 || layerNumber == 2 || layerNumber == 4) phiRegionWidth = 40. * M_PI / 180.;
  else phiRegionWidth = 20. * M_PI / 180.;
}

void ModulesToBundlesConnector::visit(RodPair& r) {
  double rodPhi = r.Phi();

  double phiSegmentWidth = (2.*M_PI) / numRods;

  // Positive cabling side
  double phiSegmentStart = femod( rodPhi, phiSegmentWidth);
  phiSegmentRef = round(femod(rodPhi - phiSegmentStart, 2.*M_PI) / phiSegmentWidth);
	
  double phiRegionStart = 0.;
  int phiRegionRef = computePhiSliceRef(rodPhi, phiRegionStart, phiRegionWidth);

  double phiSectorStart = 0.;
  int phiSectorRef = computePhiSliceRef(rodPhi, phiSectorStart, phiSectorWidth);


  // Negative cabling side
  double negPhiSegmentStart = femod( M_PI - rodPhi, phiSegmentWidth);
  negPhiSegmentRef = round(femod(M_PI - rodPhi - negPhiSegmentStart, 2.*M_PI) / phiSegmentWidth);
	
  double negPhiRegionStart = 0.;
  int negPhiRegionRef = computePhiSliceRef(M_PI - rodPhi, negPhiRegionStart, phiRegionWidth);

  double negPhiSectorStart = 0.;
  int negPhiSectorRef = computePhiSliceRef(M_PI - rodPhi, negPhiSectorStart, phiSectorWidth);


  bool isPositiveCablingSide = true;

  // CREATE 2S BUNDLES
  if (barrelName == "TB2S") {
    type = "2S";
    // Positive cabling side
    isPositiveCablingSide = true;
    bundleId = 10000 + layerNumber * 1000 + phiSegmentRef * 10;	  
    bundle = GeometryFactory::make<Bundle>(bundleId, type, barrelName, layerNumber, phiSegmentWidth, phiSegmentRef, phiRegionStart, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef, isPositiveCablingSide);
    bundles_.insert(std::make_pair(bundleId, bundle));

    // Negative cabling side
    isPositiveCablingSide = false;
    negBundleId = 30000 + layerNumber * 1000 + negPhiSegmentRef * 10;
    negBundle = GeometryFactory::make<Bundle>(negBundleId, type, barrelName, layerNumber, phiSegmentWidth, negPhiSegmentRef, negPhiRegionStart, phiRegionWidth, negPhiRegionRef, phiSectorWidth, negPhiSectorRef, isPositiveCablingSide);
    negBundles_.insert(std::make_pair(negBundleId, negBundle));
  }

  // CREATE PS BUNDLES
  else if (barrelName == "TBPS") {
    type = (layerNumber == 1 ? "PS10G" : "PS5G");

    // FLAT PART
    // Positive cabling side	  
    if ( (phiSegmentRef % 2) == 1 ) {
      isPositiveCablingSide = true;
      // standard case
      bundleFlatId = 10000 + layerNumber * 1000 + phiSegmentRef * 10 + 1;
      bundleFlat = GeometryFactory::make<Bundle>(bundleFlatId, type, barrelName, layerNumber, phiSegmentWidth, phiSegmentRef, phiRegionStart, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef, isPositiveCablingSide);
      bundles_.insert(std::make_pair(bundleFlatId, bundleFlat));

      // For layer 3, need to add a second bundle for flat part
      if (totalNumFlatRings > 12) {
	bundleFlatIdB = 10000 + layerNumber * 1000 + phiSegmentRef * 10 + 2;
	bundleFlatB = GeometryFactory::make<Bundle>(bundleFlatIdB, type, barrelName, layerNumber, phiSegmentWidth, phiSegmentRef, phiRegionStart, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef, isPositiveCablingSide);
	bundles_.insert(std::make_pair(bundleFlatIdB, bundleFlatB));
      }
    }
    // Negative cabling side
    if ( (negPhiSegmentRef % 2) == 0 ) {
      isPositiveCablingSide = false;
      // standard case
      negBundleFlatId = 30000 + layerNumber * 1000 + negPhiSegmentRef * 10 + 1;
      negBundleFlat = GeometryFactory::make<Bundle>(negBundleFlatId, type, barrelName, layerNumber, phiSegmentWidth, negPhiSegmentRef, negPhiRegionStart, phiRegionWidth, negPhiRegionRef, phiSectorWidth, negPhiSectorRef, isPositiveCablingSide);
      negBundles_.insert(std::make_pair(negBundleFlatId, negBundleFlat));

      // For layer 3, need to add a second negBundle for flat part
      if (totalNumFlatRings > 12) {
	negBundleFlatIdB = 30000 + layerNumber * 1000 + negPhiSegmentRef * 10 + 2;
	negBundleFlatB = GeometryFactory::make<Bundle>(negBundleFlatIdB, type, barrelName, layerNumber, phiSegmentWidth, negPhiSegmentRef, negPhiRegionStart, phiRegionWidth, negPhiRegionRef, phiSectorWidth, negPhiSectorRef, isPositiveCablingSide);
	negBundles_.insert(std::make_pair(negBundleFlatIdB, negBundleFlatB));
      }
    }

    // TILTED PART
    // Positive cabling side
    isPositiveCablingSide = true;
    bundleTiltedId = 10000 + layerNumber * 1000 + phiSegmentRef * 10;	  
    bundleTilted = GeometryFactory::make<Bundle>(bundleTiltedId, type, barrelName, layerNumber, phiSegmentWidth, phiSegmentRef, phiRegionStart, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef, isPositiveCablingSide);
    bundles_.insert(std::make_pair(bundleTiltedId, bundleTilted));

    // Negative cabling side
    isPositiveCablingSide = false;
    negBundleTiltedId = 30000 + layerNumber * 1000 + negPhiSegmentRef * 10;	  
    negBundleTilted = GeometryFactory::make<Bundle>(negBundleTiltedId, type, barrelName, layerNumber, phiSegmentWidth, negPhiSegmentRef, negPhiRegionStart, phiRegionWidth, negPhiRegionRef, phiSectorWidth, negPhiSectorRef, isPositiveCablingSide);
    negBundles_.insert(std::make_pair(negBundleTiltedId, negBundleTilted));	 
  }
}
      
     
void ModulesToBundlesConnector::visit(BarrelModule& m) {
  side = (m.uniRef().side > 0.);

  // CONNECT MODULES TO 2S BUNDLES
  if (barrelName == "TB2S") {
    // Positive cabling side
    if (side) {
      m.setBundle(bundle);
      bundles_[bundleId]->addModule(&m);
    }
    // Negative cabling side
    else {
      m.setBundle(negBundle);
      negBundles_[negBundleId]->addModule(&m);
    }
  }

  // CONNECT MODULES TO PS BUNDLES
  else if (barrelName == "TBPS") {

    // FLAT MODULES
    if (!m.isTilted()) {
      // Positive cabling side
      if ( (phiSegmentRef % 2) == 1) {
	// standard case
	if (totalNumFlatRings <= 12) {
	  m.setBundle(bundleFlat);
	  bundles_[bundleFlatId]->addModule(&m);
	}
	// For layer 3, need to add a second bundle for flat part
	else {
	  if (side) {
	    m.setBundle(bundleFlat);
	    bundles_[bundleFlatId]->addModule(&m);
	  } else {
	    m.setBundle(bundleFlatB);
	    bundles_[bundleFlatIdB]->addModule(&m);
	  }
	}
      }

      // Negative cabling side
      if ( (phiSegmentRef % 2) == 0) {
	// standard case
	if (totalNumFlatRings <= 12) {
	  m.setBundle(negBundleFlat);
	  negBundles_[negBundleFlatId]->addModule(&m);
	}
	// For layer 3, need to add a second bundle for flat part
	else {
	  if (!side) {
	    m.setBundle(negBundleFlat);
	    negBundles_[negBundleFlatId]->addModule(&m);
	  } else {
	    m.setBundle(negBundleFlatB);
	    negBundles_[negBundleFlatIdB]->addModule(&m);
	  }
	}
      }

    }

    // TILTED MODULES
    else if (m.isTilted()) {
      // Positive cabling side
      if (side) {
	m.setBundle(bundleTilted);
	bundles_[bundleTiltedId]->addModule(&m);
	bundles_[bundleTiltedId]->setIsTiltedPart(true);
      }
      // Negative cabling side
      else {
	m.setBundle(negBundleTilted);
	negBundles_[negBundleTiltedId]->addModule(&m);
	negBundles_[negBundleTiltedId]->setIsTiltedPart(true);
      }
    }
	  

  }
      
}


void ModulesToBundlesConnector::visit(Endcap& e) { 
  endcapName = e.myid();
}

void ModulesToBundlesConnector::visit(Disk& d) {
  diskNumber = d.diskNumber();
  side = d.side();
}

void ModulesToBundlesConnector::visit(Ring& r)   { 
  ringNumber = r.myid();
  numModulesInRing = r.numModules();

  if (endcapName == "TEDD_1") {
    if (ringNumber <= 4) type = "PS10G";
    else if (ringNumber >= 5 && ringNumber <= 7) type = "PS5GA";
    else if (ringNumber >= 8 && ringNumber <= 10) type = "PS5GB";
    else if (ringNumber >= 11) type = "2S";
  }

  else if (endcapName == "TEDD_2") {
    if (ringNumber <= 3) type = "null";
    else if (ringNumber >= 4 && ringNumber <= 6) type = "PS5GA";
    else if (ringNumber >= 7 && ringNumber <= 10) type = "PS5GB";
    else if (ringNumber >= 11) type = "2S";
  }

  if (type == "PS10G") typeIndex = 0;
  else if (type == "PS5GA") typeIndex = 1;
  else if (type == "PS5GB") typeIndex = 2;
  else if (type == "2S") typeIndex = 3;

}


void ModulesToBundlesConnector::visit(EndcapModule& m) {
  double modPhi = m.center().Phi();

  double phiSegmentWidth = (2.*M_PI) / numModulesInRing;
		  
  double phiRegionStart = 0.;
  if (type == "PS10G" || type == "PS5GA") {
    phiRegionWidth = 40. * M_PI / 180.;
  }

  else if (type == "PS5GB") {
    phiRegionWidth = 20. * M_PI / 180.;
  }

  else if (type == "2S") {
    phiRegionWidth = 360. / 27. * M_PI / 180.;
    if (endcapName == "TEDD_1") phiRegionStart = -0.55 * M_PI / 180.;
    else phiRegionStart = -0.001 * M_PI / 180.;
  }

  double phiSectorStart = 0.;

  bool isPositiveCablingSide = (side ? true : false);

  // Positive cabling side
  if (isPositiveCablingSide) {
    double phiSegmentStart = femod( modPhi, phiSegmentWidth);
    phiSegmentRef = round(femod(modPhi - phiSegmentStart, 2.*M_PI) / phiSegmentWidth);

    int phiRegionRef = computePhiSliceRef(modPhi, phiRegionStart, phiRegionWidth);
    bundleId = 20000 + diskNumber * 1000 + phiRegionRef * 10 + typeIndex;

    int phiSectorRef = computePhiSliceRef(modPhi, phiSectorStart, phiSectorWidth);

    if (bundles_.count(bundleId) == 0) {
      Bundle* bundleEndcap = GeometryFactory::make<Bundle>(bundleId, type, endcapName, diskNumber, phiSegmentWidth, phiSegmentRef, phiRegionStart, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef, isPositiveCablingSide);
      bundleEndcap->addModule(&m);
      bundles_.insert(std::make_pair(bundleId, bundleEndcap));
      m.setBundle(bundleEndcap);
    }
    else { 
      bundles_[bundleId]->addModule(&m);
      m.setBundle(bundles_[bundleId]);
    }
  }

  // Negative cabling side
  else {
    double negPhiSegmentStart = femod(M_PI - modPhi, phiSegmentWidth);
    negPhiSegmentRef = round(femod(M_PI - modPhi - negPhiSegmentStart, 2.*M_PI) / phiSegmentWidth);

    int negPhiRegionRef = computePhiSliceRef(M_PI - modPhi, phiRegionStart, phiRegionWidth);
    negBundleId = 40000 + diskNumber * 1000 + negPhiRegionRef * 10 + typeIndex;

    int negPhiSectorRef = computePhiSliceRef(M_PI - modPhi, phiSectorStart, phiSectorWidth);

    if (negBundles_.count(negBundleId) == 0) {
      Bundle* negBundleEndcap = GeometryFactory::make<Bundle>(negBundleId, type, endcapName, diskNumber, phiSegmentWidth, negPhiSegmentRef, phiRegionStart, phiRegionWidth, negPhiRegionRef, phiSectorWidth, negPhiSectorRef, isPositiveCablingSide);
      negBundleEndcap->addModule(&m);
      negBundles_.insert(std::make_pair(negBundleId, negBundleEndcap));
      m.setBundle(negBundleEndcap);
    }
    else { 
      negBundles_[negBundleId]->addModule(&m);
      m.setBundle(negBundles_[negBundleId]);
    }
  }
	
}


int ModulesToBundlesConnector::computePhiSliceRef(const double phi, const double phiSliceStart, const double phiSliceWidth) const {
  double phiSliceRefExact = femod(phi - phiSliceStart, 2.*M_PI) / phiSliceWidth;
  int phiSliceRef = 0;
  if (fabs((phiSliceRefExact - round(phiSliceRefExact))) < 0.0001) phiSliceRef = fabs(round(phiSliceRefExact));
  else phiSliceRef = std::floor(phiSliceRefExact);

  return phiSliceRef;
}


void ModulesToBundlesConnector::postVisit() {
  // STAGGER MODULES
  staggerModules(bundles_);
  staggerModules(negBundles_);

  // CHECK
  checkModulesToBundlesCabling(bundles_);
  checkModulesToBundlesCabling(negBundles_);
}


void ModulesToBundlesConnector::staggerModules(std::map<int, Bundle*>& bundles) {

  for (auto& b : bundles) {
    if (b.second->subDetectorName() == "TEDD_1" || b.second->subDetectorName() == "TEDD_2") {

      while (b.second->numModules() > 12) {
	bool isPositiveCablingSide = b.second->isPositiveCablingSide();
	int diskNumber = b.second->layerDiskNumber();

	std::string type = b.second->type();
	if (type == "PS10G") typeIndex = 0;
	else if (type == "PS5GA") typeIndex = 1;
	else if (type == "PS5GB") typeIndex = 2;
	else if (type == "2S") typeIndex = 3;

	double phiRegionStart = b.second->phiRegionStart();

	double phiRegionWidth = b.second->phiRegionWidth();
	int numPhiRegions = round(2 * M_PI / phiRegionWidth);

	int phiRegionRef = b.second->phiRegionRef();
	int nextPhiRegionRef = femod( (phiRegionRef + 1), numPhiRegions);
	int previousPhiRegionRef = femod( (phiRegionRef - 1), numPhiRegions);

	int bundleId = b.first;
	int nextBundleId = 0;
	int previousBundleId = 0;
	if (isPositiveCablingSide) {
	  nextBundleId = 20000 + diskNumber * 1000 + nextPhiRegionRef * 10 + typeIndex;
	  previousBundleId = 20000 + diskNumber * 1000 + previousPhiRegionRef * 10 + typeIndex;
	}
	else {
	  nextBundleId = 40000 + diskNumber * 1000 + nextPhiRegionRef * 10 + typeIndex;
	  previousBundleId = 40000 + diskNumber * 1000 + previousPhiRegionRef * 10 + typeIndex;
	}

	double minPhiBorder = fabs( femod((b.second->minPhi() - phiRegionStart), phiRegionWidth) );
	double maxPhiBorder = fabs( femod((b.second->maxPhi() - phiRegionStart), phiRegionWidth) - phiRegionWidth);
	      

	if (bundles.count(previousBundleId) != 0 && bundles.count(nextBundleId) != 0) {
	  // Cannot assign the extra module : both neighbouring phi regions are full !
	  if (bundles[previousBundleId]->numModules() >= 12 && bundles[nextBundleId]->numModules() >= 12) {
	    std::cout << "I am a refugee module in side " << isPositiveCablingSide << ", disk " << diskNumber << ", type " << type 
		      << ", phiRegionRef " << phiRegionRef << ", phiRegionWidth " << phiRegionWidth
		      << ", which has already more than 12 modules, and none of my neighbouring regions wants to welcome me :/" 
		      << std::endl;
	    break;
	  }

	  // Assign a module to the next phi region
	  else if (bundles[previousBundleId]->numModules() >= 12 || maxPhiBorder <= minPhiBorder) {
	    std::cout << "Removing module in side " << isPositiveCablingSide << ", disk " << diskNumber << ", type " << type 
		      << " from phiRegionRef " << phiRegionRef << ", maxPhiBorder " << (maxPhiBorder * 180. / M_PI)
		      << ", adding it to the next region." 
		      << std::endl;
	    std::cout << "my region numModules = " << b.second->numModules() << std::endl;
	    std::cout << "bundles[nextBundleId]->numModules = " << bundles[nextBundleId]->numModules() << std::endl;
	    Module* maxPhiMod = b.second->maxPhiModule();
	    maxPhiMod->setBundle(bundles[nextBundleId]);  
	    bundles[nextBundleId]->moveMaxPhiModuleFromOtherBundle(b.second);
	    std::cout << "NOWWWWWWWW my region numModules = " << b.second->numModules() << std::endl; 		  
	  }

	  // Assign a module to the previous phi region
	  else if (bundles[nextBundleId]->numModules() >= 12 || minPhiBorder < maxPhiBorder) {
	    std::cout << "Removing module in side " << isPositiveCablingSide << ", disk " << diskNumber << ", type " << type 
		      << " from phiRegionRef " << phiRegionRef << ", minPhiBorder " << (minPhiBorder * 180. / M_PI)
		      << ", adding it to the previous region." 
		      << std::endl;
	    std::cout << "my region numModules = " << b.second->numModules() << std::endl;
	    std::cout << "bundles[previousBundleId]->numModules = " << bundles[previousBundleId]->numModules() << std::endl;
	    Module* minPhiMod = b.second->minPhiModule();
	    minPhiMod->setBundle(bundles[previousBundleId]);	  
	    bundles[previousBundleId]->moveMinPhiModuleFromOtherBundle(b.second);
	    std::cout << "NOWWWWWWWW my region numModules = " << b.second->numModules() << std::endl;		  
	  }
	}
	else { std::cout << "Error building previousBundleId or nextBundleId" << std::endl; break; }

      }
    }
  }

}



void ModulesToBundlesConnector::checkModulesToBundlesCabling(std::map<int, Bundle*>& bundles) {
  for (auto& b : bundles) {
    if (b.second->numModules() > 12) {
      std::cout << "There was an error while staggering modules. Bundle " 
		<< b.first << " is connected to " << b.second->numModules() << " modules." 
		<< std::endl;
    }

    if (b.second->phiSegmentRef() <= -1 || b.second->phiRegionRef() <= -1 || b.second->phiSectorRef() <= -1) {
      std::cout << "Error while creating bundle. Bundle " << b.first << " has phiSegmentRef = " << b.second->phiSegmentRef() << ", phiRegionRef = " << b.second->phiRegionRef() << ", phiSectorRef = " << b.second->phiSectorRef() << ", type = " << b.second->type() << std::endl;
    }
  }
}
