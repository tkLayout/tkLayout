#include "Cabling/ModulesToBundlesConnector.hh"
#include <Tracker.hh>


void ModulesToBundlesConnector::visit(Barrel& b) {
  barrelName_ = b.myid();	
}


void ModulesToBundlesConnector::visit(Layer& l) {
  layerNumber_ = l.layerNumber();
  numRods_ = l.numRods();
  totalNumFlatRings_ = l.buildNumModulesFlat() * 2 - 1; // Total number of flat rings on both (+Z) side and (-Z) side

  if (layerNumber_ == 1 || layerNumber_ == 2 || layerNumber_ == 4) phiRegionWidth_ = 40. * M_PI / 180.;
  else phiRegionWidth_ = 20. * M_PI / 180.;
}


void ModulesToBundlesConnector::visit(RodPair& r) {
  double rodPhi = r.Phi();

  double phiSegmentWidth = (2.*M_PI) / numRods_;

  // Positive cabling side
  bool isPositiveCablingSide = true;

  double phiSegmentStart = computePhiSegmentStart(rodPhi, phiSegmentWidth, isPositiveCablingSide);
  phiSegmentRef_ = computePhiSegmentRef(rodPhi, phiSegmentStart, phiSegmentWidth, isPositiveCablingSide);
	
  double phiRegionStart = 0.;
  int phiRegionRef = computePhiSliceRef(rodPhi, phiRegionStart, phiRegionWidth_, isPositiveCablingSide);

  double phiSectorStart = 0.;
  int phiSectorRef = computePhiSliceRef(rodPhi, phiSectorStart, phiSectorWidth_, isPositiveCablingSide);


  // Negative cabling side
  isPositiveCablingSide = false;

  double negPhiSegmentStart = computePhiSegmentStart(rodPhi, phiSegmentWidth, isPositiveCablingSide);
  negPhiSegmentRef_ = computePhiSegmentRef(rodPhi, negPhiSegmentStart, phiSegmentWidth, isPositiveCablingSide);
	
  double negPhiRegionStart = 0.;
  int negPhiRegionRef = computePhiSliceRef(rodPhi, negPhiRegionStart, phiRegionWidth_, isPositiveCablingSide);

  double negPhiSectorStart = 0.;
  int negPhiSectorRef = computePhiSliceRef(rodPhi, negPhiSectorStart, phiSectorWidth_, isPositiveCablingSide);


  bool isBarrel = true;
  isPositiveCablingSide = true;

  // CREATE 2S BUNDLES
  if (barrelName_ == "TB2S") {
    type_ = "2S";
    // Positive cabling side
    isPositiveCablingSide = true;
    typeIndex_ = 0;
    bundleId_ = computeBundleId(isBarrel, isPositiveCablingSide, layerNumber_, phiSegmentRef_, typeIndex_);
    createAndStoreBundle(bundles_, negBundles_, bundleId_, type_, barrelName_, layerNumber_, phiSegmentWidth, phiSegmentRef_, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide);

    // Negative cabling side
    isPositiveCablingSide = false;
    typeIndex_ = 0;
    negBundleId_ = computeBundleId(isBarrel, isPositiveCablingSide, layerNumber_, negPhiSegmentRef_, typeIndex_);
    createAndStoreBundle(bundles_, negBundles_, negBundleId_, type_, barrelName_, layerNumber_, phiSegmentWidth, negPhiSegmentRef_, negPhiRegionStart, phiRegionWidth_, negPhiRegionRef, phiSectorWidth_, negPhiSectorRef, isPositiveCablingSide);
  }

  // CREATE PS BUNDLES
  else if (barrelName_ == "TBPS") {
    type_ = (layerNumber_ == 1 ? "PS10G" : "PS5G");

    // FLAT PART
    // Positive cabling side	  
    if ( (phiSegmentRef_ % 2) == 1 ) {
      isPositiveCablingSide = true;
      // standard case
      typeIndex_ = 1;
      bundleFlatId_ = computeBundleId(isBarrel, isPositiveCablingSide, layerNumber_, phiSegmentRef_, typeIndex_);
      createAndStoreBundle(bundles_, negBundles_, bundleFlatId_, type_, barrelName_, layerNumber_, phiSegmentWidth, phiSegmentRef_, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide);

      // For layer 3, need to add a second bundle for flat part
      if (totalNumFlatRings_ > maxNumModulesPerBundle_) {
	typeIndex_ = 2;
	bundleFlatIdB_ = computeBundleId(isBarrel, isPositiveCablingSide, layerNumber_, phiSegmentRef_, typeIndex_);
	createAndStoreBundle(bundles_, negBundles_, bundleFlatIdB_, type_, barrelName_, layerNumber_, phiSegmentWidth, phiSegmentRef_, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide);
      }
    }
    // Negative cabling side
    if ( (negPhiSegmentRef_ % 2) == 0 ) {
      isPositiveCablingSide = false;
      // standard case
      typeIndex_ = 1;
      negBundleFlatId_ = computeBundleId(isBarrel, isPositiveCablingSide, layerNumber_, negPhiSegmentRef_, typeIndex_);
      createAndStoreBundle(bundles_, negBundles_, negBundleFlatId_, type_, barrelName_, layerNumber_, phiSegmentWidth, negPhiSegmentRef_, negPhiRegionStart, phiRegionWidth_, negPhiRegionRef, phiSectorWidth_, negPhiSectorRef, isPositiveCablingSide);

      // For layer 3, need to add a second negBundle for flat part
      if (totalNumFlatRings_ > maxNumModulesPerBundle_) {
	typeIndex_ = 2;
	negBundleFlatIdB_ = computeBundleId(isBarrel, isPositiveCablingSide, layerNumber_, negPhiSegmentRef_, typeIndex_);
	createAndStoreBundle(bundles_, negBundles_, negBundleFlatIdB_, type_, barrelName_, layerNumber_, phiSegmentWidth, negPhiSegmentRef_, negPhiRegionStart, phiRegionWidth_, negPhiRegionRef, phiSectorWidth_, negPhiSectorRef, isPositiveCablingSide);
      }
    }

    // TILTED PART
    // Positive cabling side
    isPositiveCablingSide = true;
    typeIndex_ = 0;
    bundleTiltedId_ = computeBundleId(isBarrel, isPositiveCablingSide, layerNumber_, phiSegmentRef_, typeIndex_);
    createAndStoreBundle(bundles_, negBundles_, bundleTiltedId_, type_, barrelName_, layerNumber_, phiSegmentWidth, phiSegmentRef_, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide);
    
    // Negative cabling side
    isPositiveCablingSide = false;
    typeIndex_ = 0;
    negBundleTiltedId_ = computeBundleId(isBarrel, isPositiveCablingSide, layerNumber_, negPhiSegmentRef_, typeIndex_);
    createAndStoreBundle(bundles_, negBundles_, negBundleTiltedId_, type_, barrelName_, layerNumber_, phiSegmentWidth, negPhiSegmentRef_, negPhiRegionStart, phiRegionWidth_, negPhiRegionRef, phiSectorWidth_, negPhiSectorRef, isPositiveCablingSide);
  }
}
      
     
void ModulesToBundlesConnector::visit(BarrelModule& m) {
  side_ = (m.uniRef().side > 0.);

  // CONNECT MODULES TO 2S BUNDLES
  if (barrelName_ == "TB2S") {
    // Positive cabling side
    if (side_) {
      Bundle* bundle = bundles_[bundleId_];
      m.setBundle(bundle);
      bundle->addModule(&m);
    }
    // Negative cabling side
    else {
      Bundle* negBundle = negBundles_[negBundleId_];
      m.setBundle(negBundle);
      negBundle->addModule(&m);
    }
  }

  // CONNECT MODULES TO PS BUNDLES
  else if (barrelName_ == "TBPS") {

    // FLAT MODULES
    if (!m.isTilted()) {
      // Positive cabling side
      if ( (phiSegmentRef_ % 2) == 1) {
	// standard case
	if (totalNumFlatRings_ <= maxNumModulesPerBundle_) {
	  Bundle* bundleFlat = bundles_[bundleFlatId_];
	  m.setBundle(bundleFlat);
	  bundleFlat->addModule(&m);
	}
	// For layer 3, need to add a second bundle for flat part
	else {
	  if (side_) {
	    Bundle* bundleFlat = bundles_[bundleFlatId_];
	    m.setBundle(bundleFlat);
	    bundleFlat->addModule(&m);
	  } else {
	    Bundle* bundleFlatB = bundles_[bundleFlatIdB_];
	    m.setBundle(bundleFlatB);
	    bundleFlatB->addModule(&m);
	  }
	}
      }

      // Negative cabling side
      if ( (phiSegmentRef_ % 2) == 0) {
	// standard case
	if (totalNumFlatRings_ <= maxNumModulesPerBundle_) {
	  Bundle* negBundleFlat = negBundles_[negBundleFlatId_];
	  m.setBundle(negBundleFlat);
	  negBundleFlat->addModule(&m);
	}
	// For layer 3, need to add a second bundle for flat part
	else {
	  if (!side_) {
	    Bundle* negBundleFlat = negBundles_[negBundleFlatId_];
	    m.setBundle(negBundleFlat);
	    negBundleFlat->addModule(&m);
	  } else {
	    Bundle* negBundleFlatB = negBundles_[negBundleFlatIdB_];
	    m.setBundle(negBundleFlatB);
	    negBundleFlatB->addModule(&m);
	  }
	}
      }

    }

    // TILTED MODULES
    else if (m.isTilted()) {
      // Positive cabling side
      if (side_) {
	Bundle* bundleTilted = bundles_[bundleTiltedId_];
	m.setBundle(bundleTilted);
	bundleTilted->addModule(&m);
	bundleTilted->setIsTiltedPart(true);
      }
      // Negative cabling side
      else {
	Bundle* negBundleTilted = negBundles_[negBundleTiltedId_];
	m.setBundle(negBundleTilted);
	negBundleTilted->addModule(&m);
	negBundleTilted->setIsTiltedPart(true);
      }
    }
	  

  }
      
}


void ModulesToBundlesConnector::visit(Endcap& e) { 
  endcapName_ = e.myid();
}


void ModulesToBundlesConnector::visit(Disk& d) {
  diskNumber_ = d.diskNumber();
  side_ = d.side();
}


void ModulesToBundlesConnector::visit(Ring& r)   { 
  ringNumber_ = r.myid();
  numModulesInRing_ = r.numModules();

  if (endcapName_ == "TEDD_1") {
    if (ringNumber_ <= 4) type_ = "PS10G";
    else if (ringNumber_ >= 5 && ringNumber_ <= 7) type_ = "PS5GA";
    else if (ringNumber_ >= 8 && ringNumber_ <= 10) type_ = "PS5GB";
    else if (ringNumber_ >= 11) type_ = "2S";
  }

  else if (endcapName_ == "TEDD_2") {
    if (ringNumber_ <= 3) type_ = "null";
    else if (ringNumber_ >= 4 && ringNumber_ <= 6) type_ = "PS5GA";
    else if (ringNumber_ >= 7 && ringNumber_ <= 10) type_ = "PS5GB";
    else if (ringNumber_ >= 11) type_ = "2S";
  }

  if (type_ == "PS10G") typeIndex_ = 0;
  else if (type_ == "PS5GA") typeIndex_ = 1;
  else if (type_ == "PS5GB") typeIndex_ = 2;
  else if (type_ == "2S") typeIndex_ = 3;

}


void ModulesToBundlesConnector::visit(EndcapModule& m) {
  double modPhi = m.center().Phi();

  double phiSegmentWidth = (2.*M_PI) / numModulesInRing_;
		  
  double phiRegionStart = 0.;
  if (type_ == "PS10G" || type_ == "PS5GA") {
    phiRegionWidth_ = 40. * M_PI / 180.;
  }

  else if (type_ == "PS5GB") {
    phiRegionWidth_ = 20. * M_PI / 180.;
  }

  else if (type_ == "2S") {
    phiRegionWidth_ = 360. / 27. * M_PI / 180.;
    if (endcapName_ == "TEDD_1") phiRegionStart = -0.55 * M_PI / 180.;
    else phiRegionStart = -0.001 * M_PI / 180.;
  }

  double phiSectorStart = 0.;

  bool isBarrel = false;
  bool isPositiveCablingSide = side_;

 
  double phiSegmentStart = computePhiSegmentStart(modPhi, phiSegmentWidth, isPositiveCablingSide);
  phiSegmentRef_ = computePhiSegmentRef(modPhi, phiSegmentStart, phiSegmentWidth, isPositiveCablingSide);

  int phiRegionRef = computePhiSliceRef(modPhi, phiRegionStart, phiRegionWidth_, isPositiveCablingSide);
  bundleId_ = computeBundleId(isBarrel, isPositiveCablingSide, diskNumber_, phiRegionRef, typeIndex_);

  int phiSectorRef = computePhiSliceRef(modPhi, phiSectorStart, phiSectorWidth_, isPositiveCablingSide);

  Bundle* bundleEndcap = nullptr;
  std::map<int, Bundle*>& bundles = (isPositiveCablingSide ? bundles_ : negBundles_);
  auto found = bundles.find(bundleId_);
  if (found == bundles.end()) {
    createAndStoreBundle(bundles_, negBundles_, bundleId_, type_, endcapName_, diskNumber_, phiSegmentWidth, phiSegmentRef_, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide);
    bundleEndcap = bundles[bundleId_];
  }
  else {
    bundleEndcap = found->second;
  }
  bundleEndcap->addModule(&m);
  m.setBundle(bundleEndcap);	
}



void ModulesToBundlesConnector::postVisit() {
  // STAGGER MODULES
  staggerModules(bundles_);
  staggerModules(negBundles_);

  // CHECK
  checkModulesToBundlesCabling(bundles_);
  checkModulesToBundlesCabling(negBundles_);
}


double ModulesToBundlesConnector::computePhiSegmentStart(const double phi, const double phiSegmentWidth, const bool isPositiveCablingSide) const {
  double stereoPhi = (isPositiveCablingSide ? phi : M_PI - phi);
  double phiSegmentStart = femod(stereoPhi, phiSegmentWidth);
  return phiSegmentStart;
}


int ModulesToBundlesConnector::computePhiSegmentRef(const double phi, const double phiSegmentStart, const double phiSegmentWidth, const bool isPositiveCablingSide) const {
  double stereoPhi = (isPositiveCablingSide ? phi : M_PI - phi);
  int phiSegmentRef = round(femod(stereoPhi - phiSegmentStart, 2.*M_PI) / phiSegmentWidth);
  return phiSegmentRef;
}


int ModulesToBundlesConnector::computePhiSliceRef(const double phi, const double phiSliceStart, const double phiSliceWidth, const bool isPositiveCablingSide) const {
  double stereoPhi = (isPositiveCablingSide ? phi : M_PI - phi);
  double phiSliceRefExact = femod(stereoPhi - phiSliceStart, 2.*M_PI) / phiSliceWidth;
  int phiSliceRef = 0;
  if (fabs((phiSliceRefExact - round(phiSliceRefExact))) < 0.0001) phiSliceRef = fabs(round(phiSliceRefExact));
  else phiSliceRef = std::floor(phiSliceRefExact);

  return phiSliceRef;
}


int ModulesToBundlesConnector::computeBundleId(const bool isBarrel, const bool isPositiveCablingSide, const int layerDiskNumber, const int phiRef, const int typeIndex) {
  int cablingSideIndex = 0;
  if (isBarrel) {
    cablingSideIndex = (isPositiveCablingSide ? 1 : 3);
  }
  else {
    cablingSideIndex = (isPositiveCablingSide ? 2 : 4);
  }

  int bundleId = cablingSideIndex * 10000 + layerDiskNumber * 1000 + phiRef * 10 + typeIndex;
  return bundleId;
}



void ModulesToBundlesConnector::createAndStoreBundle(std::map<int, Bundle*>& bundles, std::map<int, Bundle*>& negBundles, const int bundleId, const std::string type, const std::string subDetectorName, const int layerDiskNumber, const double phiSegmentWidth, const int phiSegmentRef, const double phiRegionStart, const double phiRegionWidth, const int phiRegionRef, const double phiSectorWidth, const int phiSectorRef, const bool isPositiveCablingSide) {

  Bundle* bundle = GeometryFactory::make<Bundle>(bundleId, type, subDetectorName, layerDiskNumber, phiSegmentWidth, phiSegmentRef, phiRegionStart, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef, isPositiveCablingSide);

  if (isPositiveCablingSide) {
    bundles.insert(std::make_pair(bundleId, bundle));
  }
  else {
    negBundles.insert(std::make_pair(bundleId, bundle));
  }
}


void ModulesToBundlesConnector::staggerModules(std::map<int, Bundle*>& bundles) {

  for (auto& b : bundles) {
    if (b.second->subDetectorName() == "TEDD_1" || b.second->subDetectorName() == "TEDD_2") {

      while (b.second->numModules() > maxNumModulesPerBundle_) {
	bool isPositiveCablingSide = b.second->isPositiveCablingSide();
	int diskNumber = b.second->layerDiskNumber();

	std::string type = b.second->type();
	int typeIndex;
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
	  if (bundles[previousBundleId]->numModules() >= maxNumModulesPerBundle_ && bundles[nextBundleId]->numModules() >= maxNumModulesPerBundle_) {
	    std::cout << "I am a refugee module in side " << isPositiveCablingSide << ", disk " << diskNumber << ", type " << type 
		      << ", phiRegionRef " << phiRegionRef << ", phiRegionWidth " << phiRegionWidth
		      << ", which has already more than " << maxNumModulesPerBundle_<< " modules, and none of my neighbouring regions wants to welcome me :/" 
		      << std::endl;
	    break;
	  }

	  // Assign a module to the next phi region
	  else if (bundles[previousBundleId]->numModules() >= maxNumModulesPerBundle_ || maxPhiBorder <= minPhiBorder) {
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
	  else if (bundles[nextBundleId]->numModules() >= maxNumModulesPerBundle_ || minPhiBorder < maxPhiBorder) {
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


void ModulesToBundlesConnector::checkModulesToBundlesCabling(const std::map<int, Bundle*>& bundles) const {
  for (auto& b : bundles) {
    if (b.second->numModules() > maxNumModulesPerBundle_) {
      std::cout << "There was an error while staggering modules. Bundle " 
		<< b.first << " is connected to " << b.second->numModules() << " modules." 
		<< std::endl;
    }

    if (b.second->phiSegmentRef() <= -1 || b.second->phiRegionRef() <= -1 || b.second->phiSectorRef() <= -1) {
      std::cout << "Error while creating bundle. Bundle " << b.first << " has phiSegmentRef = " << b.second->phiSegmentRef() << ", phiRegionRef = " << b.second->phiRegionRef() << ", phiSectorRef = " << b.second->phiSectorRef() << ", type = " << b.second->type() << std::endl;
    }
  }
}
