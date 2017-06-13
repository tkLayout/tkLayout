#include "Cabling/ModulesToBundlesConnector.hh"
#include <Tracker.hh>


void ModulesToBundlesConnector::visit(Barrel& b) {
  barrelName_ = b.myid();	
}

void ModulesToBundlesConnector::visit(Layer& l) {
  layerNumber_ = l.layerNumber();
  numRods_ = l.numRods();
  totalNumFlatRings_ = l.buildNumModulesFlat() * 2 - 1; // Total number of flat rings on both +Z side and -Z side

  if (layerNumber_ == 1 || layerNumber_ == 2 || layerNumber_ == 4) phiRegionWidth_ = 40. * M_PI / 180.;
  else phiRegionWidth_ = 20. * M_PI / 180.;
}

void ModulesToBundlesConnector::visit(RodPair& r) {
  double rodPhi = r.Phi();

  double phiSegmentWidth = (2.*M_PI) / numRods_;

  // Positive cabling side
  double phiSegmentStart = femod( rodPhi, phiSegmentWidth);
  phiSegmentRef_ = round(femod(rodPhi - phiSegmentStart, 2.*M_PI) / phiSegmentWidth);
	
  double phiRegionStart = 0.;
  int phiRegionRef = computePhiSliceRef(rodPhi, phiRegionStart, phiRegionWidth_);

  double phiSectorStart = 0.;
  int phiSectorRef = computePhiSliceRef(rodPhi, phiSectorStart, phiSectorWidth_);


  // Negative cabling side
  double negPhiSegmentStart = femod( M_PI - rodPhi, phiSegmentWidth);
  negPhiSegmentRef_ = round(femod(M_PI - rodPhi - negPhiSegmentStart, 2.*M_PI) / phiSegmentWidth);
	
  double negPhiRegionStart = 0.;
  int negPhiRegionRef = computePhiSliceRef(M_PI - rodPhi, negPhiRegionStart, phiRegionWidth_);

  double negPhiSectorStart = 0.;
  int negPhiSectorRef = computePhiSliceRef(M_PI - rodPhi, negPhiSectorStart, phiSectorWidth_);


  bool isPositiveCablingSide = true;

  // CREATE 2S BUNDLES
  if (barrelName_ == "TB2S") {
    type_ = "2S";
    // Positive cabling side
    isPositiveCablingSide = true;
    bundleId_ = 10000 + layerNumber_ * 1000 + phiSegmentRef_ * 10;	  
    bundle_ = GeometryFactory::make<Bundle>(bundleId_, type_, barrelName_, layerNumber_, phiSegmentWidth, phiSegmentRef_, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide);
    bundles_.insert(std::make_pair(bundleId_, bundle_));

    // Negative cabling side
    isPositiveCablingSide = false;
    negBundleId_ = 30000 + layerNumber_ * 1000 + negPhiSegmentRef_ * 10;
    negBundle_ = GeometryFactory::make<Bundle>(negBundleId_, type_, barrelName_, layerNumber_, phiSegmentWidth, negPhiSegmentRef_, negPhiRegionStart, phiRegionWidth_, negPhiRegionRef, phiSectorWidth_, negPhiSectorRef, isPositiveCablingSide);
    negBundles_.insert(std::make_pair(negBundleId_, negBundle_));
  }

  // CREATE PS BUNDLES
  else if (barrelName_ == "TBPS") {
    type_ = (layerNumber_ == 1 ? "PS10G" : "PS5G");

    // FLAT PART
    // Positive cabling side	  
    if ( (phiSegmentRef_ % 2) == 1 ) {
      isPositiveCablingSide = true;
      // standard case
      bundleFlatId_ = 10000 + layerNumber_ * 1000 + phiSegmentRef_ * 10 + 1;
      bundleFlat_ = GeometryFactory::make<Bundle>(bundleFlatId_, type_, barrelName_, layerNumber_, phiSegmentWidth, phiSegmentRef_, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide);
      bundles_.insert(std::make_pair(bundleFlatId_, bundleFlat_));

      // For layer 3, need to add a second bundle for flat part
      if (totalNumFlatRings_ > maxNumModulesPerBundle_) {
	bundleFlatIdB_ = 10000 + layerNumber_ * 1000 + phiSegmentRef_ * 10 + 2;
	bundleFlatB_ = GeometryFactory::make<Bundle>(bundleFlatIdB_, type_, barrelName_, layerNumber_, phiSegmentWidth, phiSegmentRef_, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide);
	bundles_.insert(std::make_pair(bundleFlatIdB_, bundleFlatB_));
      }
    }
    // Negative cabling side
    if ( (negPhiSegmentRef_ % 2) == 0 ) {
      isPositiveCablingSide = false;
      // standard case
      negBundleFlatId_ = 30000 + layerNumber_ * 1000 + negPhiSegmentRef_ * 10 + 1;
      negBundleFlat_ = GeometryFactory::make<Bundle>(negBundleFlatId_, type_, barrelName_, layerNumber_, phiSegmentWidth, negPhiSegmentRef_, negPhiRegionStart, phiRegionWidth_, negPhiRegionRef, phiSectorWidth_, negPhiSectorRef, isPositiveCablingSide);
      negBundles_.insert(std::make_pair(negBundleFlatId_, negBundleFlat_));

      // For layer 3, need to add a second negBundle for flat part
      if (totalNumFlatRings_ > maxNumModulesPerBundle_) {
	negBundleFlatIdB_ = 30000 + layerNumber_ * 1000 + negPhiSegmentRef_ * 10 + 2;
	negBundleFlatB_ = GeometryFactory::make<Bundle>(negBundleFlatIdB_, type_, barrelName_, layerNumber_, phiSegmentWidth, negPhiSegmentRef_, negPhiRegionStart, phiRegionWidth_, negPhiRegionRef, phiSectorWidth_, negPhiSectorRef, isPositiveCablingSide);
	negBundles_.insert(std::make_pair(negBundleFlatIdB_, negBundleFlatB_));
      }
    }

    // TILTED PART
    // Positive cabling side
    isPositiveCablingSide = true;
    bundleTiltedId_ = 10000 + layerNumber_ * 1000 + phiSegmentRef_ * 10;	  
    bundleTilted_ = GeometryFactory::make<Bundle>(bundleTiltedId_, type_, barrelName_, layerNumber_, phiSegmentWidth, phiSegmentRef_, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide);
    bundles_.insert(std::make_pair(bundleTiltedId_, bundleTilted_));

    // Negative cabling side
    isPositiveCablingSide = false;
    negBundleTiltedId_ = 30000 + layerNumber_ * 1000 + negPhiSegmentRef_ * 10;	  
    negBundleTilted_ = GeometryFactory::make<Bundle>(negBundleTiltedId_, type_, barrelName_, layerNumber_, phiSegmentWidth, negPhiSegmentRef_, negPhiRegionStart, phiRegionWidth_, negPhiRegionRef, phiSectorWidth_, negPhiSectorRef, isPositiveCablingSide);
    negBundles_.insert(std::make_pair(negBundleTiltedId_, negBundleTilted_));	 
  }
}
      
     
void ModulesToBundlesConnector::visit(BarrelModule& m) {
  side_ = (m.uniRef().side > 0.);

  // CONNECT MODULES TO 2S BUNDLES
  if (barrelName_ == "TB2S") {
    // Positive cabling side
    if (side_) {
      m.setBundle(bundle_);
      bundles_[bundleId_]->addModule(&m);
    }
    // Negative cabling side
    else {
      m.setBundle(negBundle_);
      negBundles_[negBundleId_]->addModule(&m);
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
	  m.setBundle(bundleFlat_);
	  bundles_[bundleFlatId_]->addModule(&m);
	}
	// For layer 3, need to add a second bundle for flat part
	else {
	  if (side_) {
	    m.setBundle(bundleFlat_);
	    bundles_[bundleFlatId_]->addModule(&m);
	  } else {
	    m.setBundle(bundleFlatB_);
	    bundles_[bundleFlatIdB_]->addModule(&m);
	  }
	}
      }

      // Negative cabling side
      if ( (phiSegmentRef_ % 2) == 0) {
	// standard case
	if (totalNumFlatRings_ <= maxNumModulesPerBundle_) {
	  m.setBundle(negBundleFlat_);
	  negBundles_[negBundleFlatId_]->addModule(&m);
	}
	// For layer 3, need to add a second bundle for flat part
	else {
	  if (!side_) {
	    m.setBundle(negBundleFlat_);
	    negBundles_[negBundleFlatId_]->addModule(&m);
	  } else {
	    m.setBundle(negBundleFlatB_);
	    negBundles_[negBundleFlatIdB_]->addModule(&m);
	  }
	}
      }

    }

    // TILTED MODULES
    else if (m.isTilted()) {
      // Positive cabling side
      if (side_) {
	m.setBundle(bundleTilted_);
	bundles_[bundleTiltedId_]->addModule(&m);
	bundles_[bundleTiltedId_]->setIsTiltedPart(true);
      }
      // Negative cabling side
      else {
	m.setBundle(negBundleTilted_);
	negBundles_[negBundleTiltedId_]->addModule(&m);
	negBundles_[negBundleTiltedId_]->setIsTiltedPart(true);
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

  bool isPositiveCablingSide = side_;

  // Positive cabling side
  if (isPositiveCablingSide) {
    double phiSegmentStart = femod( modPhi, phiSegmentWidth);
    phiSegmentRef_ = round(femod(modPhi - phiSegmentStart, 2.*M_PI) / phiSegmentWidth);

    int phiRegionRef = computePhiSliceRef(modPhi, phiRegionStart, phiRegionWidth_);
    bundleId_ = 20000 + diskNumber_ * 1000 + phiRegionRef * 10 + typeIndex_;

    int phiSectorRef = computePhiSliceRef(modPhi, phiSectorStart, phiSectorWidth_);

    if (bundles_.count(bundleId_) == 0) {
      Bundle* bundleEndcap = GeometryFactory::make<Bundle>(bundleId_, type_, endcapName_, diskNumber_, phiSegmentWidth, phiSegmentRef_, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide);
      bundleEndcap->addModule(&m);
      bundles_.insert(std::make_pair(bundleId_, bundleEndcap));
      m.setBundle(bundleEndcap);
    }
    else { 
      bundles_[bundleId_]->addModule(&m);
      m.setBundle(bundles_[bundleId_]);
    }
  }

  // Negative cabling side
  else {
    double negPhiSegmentStart = femod(M_PI - modPhi, phiSegmentWidth);
    negPhiSegmentRef_ = round(femod(M_PI - modPhi - negPhiSegmentStart, 2.*M_PI) / phiSegmentWidth);

    int negPhiRegionRef = computePhiSliceRef(M_PI - modPhi, phiRegionStart, phiRegionWidth_);
    negBundleId_ = 40000 + diskNumber_ * 1000 + negPhiRegionRef * 10 + typeIndex_;

    int negPhiSectorRef = computePhiSliceRef(M_PI - modPhi, phiSectorStart, phiSectorWidth_);

    if (negBundles_.count(negBundleId_) == 0) {
      Bundle* negBundleEndcap = GeometryFactory::make<Bundle>(negBundleId_, type_, endcapName_, diskNumber_, phiSegmentWidth, negPhiSegmentRef_, phiRegionStart, phiRegionWidth_, negPhiRegionRef, phiSectorWidth_, negPhiSectorRef, isPositiveCablingSide);
      negBundleEndcap->addModule(&m);
      negBundles_.insert(std::make_pair(negBundleId_, negBundleEndcap));
      m.setBundle(negBundleEndcap);
    }
    else { 
      negBundles_[negBundleId_]->addModule(&m);
      m.setBundle(negBundles_[negBundleId_]);
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



void ModulesToBundlesConnector::checkModulesToBundlesCabling(std::map<int, Bundle*>& bundles) {
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
