#include "Cabling/ModulesToBundlesConnector.hh"
#include <Tracker.hh>


void ModulesToBundlesConnector::visit(Barrel& b) {
  isBarrel_ = true;
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
  bundleType_ = computeBundleType(isBarrel_, barrelName_, layerNumber_);

  rodPhi_ = r.Phi();
}
      
     
void ModulesToBundlesConnector::visit(BarrelModule& m) {
  double phiSegmentWidth = (2.*M_PI) / numRods_;

  side_ = (m.uniRef().side > 0.);

  bool isPositiveCablingSide = side_;
  bool isTilted = false;
  bool isExtraFlatPart = false;

  // 2S BUNDLES
  if (barrelName_ == "TB2S") {

    isPositiveCablingSide = side_;
    isTilted = false;
    isExtraFlatPart = false;
  }

  // PS BUNDLES
  else if (barrelName_ == "TBPS") {

    // FLAT PART
    if (!m.isTilted()) {
      isPositiveCablingSide = computeBarrelFlatPartRodCablingSide(rodPhi_, phiSegmentWidth);
      isTilted = false;
      isExtraFlatPart = false;

      // For layer 3, need to add a second bundle for flat part
      if (isPositiveCablingSide && (totalNumFlatRings_ > maxNumModulesPerBundle_) && !side_) isExtraFlatPart = true;
      if (!isPositiveCablingSide && (totalNumFlatRings_ > maxNumModulesPerBundle_) && side_) isExtraFlatPart = true;
    }

    // TILTED PART
    else if (m.isTilted()) {

      isPositiveCablingSide = side_;
      isTilted = true;
      isExtraFlatPart = false;
    }       
  }

  double phiSegmentStart = computePhiSegmentStart(rodPhi_, phiSegmentWidth, isPositiveCablingSide);
  int phiSegmentRef = computePhiSegmentRef(rodPhi_, phiSegmentStart, phiSegmentWidth, isPositiveCablingSide);
	
  double phiRegionStart = 0.;
  int phiRegionRef = computePhiSliceRef(rodPhi_, phiRegionStart, phiRegionWidth_, isPositiveCablingSide);

  double phiSectorStart = 0.;
  int phiSectorRef = computePhiSliceRef(rodPhi_, phiSectorStart, phiSectorWidth_, isPositiveCablingSide);

  bundleTypeIndex_ = computeBundleTypeIndex(isBarrel_, bundleType_, totalNumFlatRings_, maxNumModulesPerBundle_, isTilted, isExtraFlatPart);
  int bundleId = computeBundleId(isBarrel_, isPositiveCablingSide, layerNumber_, phiSegmentRef, bundleTypeIndex_);

  // BUILD BUNDLE IF NECESSARY, AND CONNECT MODULE TO BUNDLE
  buildBundle(m, bundles_, negBundles_, bundleId, bundleType_, barrelName_, layerNumber_, phiSegmentWidth, phiSegmentRef, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide, isTilted);
}


void ModulesToBundlesConnector::visit(Endcap& e) {
  isBarrel_ = false;
  endcapName_ = e.myid();
}


void ModulesToBundlesConnector::visit(Disk& d) {
  diskNumber_ = d.diskNumber();
  side_ = d.side();
}


void ModulesToBundlesConnector::visit(Ring& r)   { 
  ringNumber_ = r.myid();
  numModulesInRing_ = r.numModules();

  bundleType_ = computeBundleType(isBarrel_, endcapName_, diskNumber_, ringNumber_);
  bundleTypeIndex_ = computeBundleTypeIndex(isBarrel_, bundleType_);
}


void ModulesToBundlesConnector::visit(EndcapModule& m) {
  double modPhi = m.center().Phi();

  double phiSegmentWidth = (2.*M_PI) / numModulesInRing_;
		  
  double phiRegionStart = 0.;
  if (bundleType_ == "PS10G" || bundleType_ == "PS5GA") {
    phiRegionWidth_ = 40. * M_PI / 180.;
  }

  else if (bundleType_ == "PS5GB") {
    phiRegionWidth_ = 20. * M_PI / 180.;
  }

  else if (bundleType_ == "2S") {
    phiRegionWidth_ = 360. / 27. * M_PI / 180.;
    if (endcapName_ == "TEDD_1") phiRegionStart = -0.55 * M_PI / 180.;
    else phiRegionStart = -0.001 * M_PI / 180.;
  }

  double phiSectorStart = 0.;

  bool isPositiveCablingSide = side_;

 
  double phiSegmentStart = computePhiSegmentStart(modPhi, phiSegmentWidth, isPositiveCablingSide);
  int phiSegmentRef = computePhiSegmentRef(modPhi, phiSegmentStart, phiSegmentWidth, isPositiveCablingSide);

  int phiRegionRef = computePhiSliceRef(modPhi, phiRegionStart, phiRegionWidth_, isPositiveCablingSide);
  int bundleId = computeBundleId(isBarrel_, isPositiveCablingSide, diskNumber_, phiRegionRef, bundleTypeIndex_);

  int phiSectorRef = computePhiSliceRef(modPhi, phiSectorStart, phiSectorWidth_, isPositiveCablingSide);


  // BUILD BUNDLE IF NECESSARY, AND CONNECT MODULE TO BUNDLE
  buildBundle(m, bundles_, negBundles_, bundleId, bundleType_, endcapName_, diskNumber_, phiSegmentWidth, phiSegmentRef, phiRegionStart, phiRegionWidth_, phiRegionRef, phiSectorWidth_, phiSectorRef, isPositiveCablingSide);
}



void ModulesToBundlesConnector::postVisit() {
  // STAGGER MODULES
  staggerModules(bundles_);
  staggerModules(negBundles_);

  // CHECK
  checkModulesToBundlesCabling(bundles_);
  checkModulesToBundlesCabling(negBundles_);
}


bool ModulesToBundlesConnector::computeBarrelFlatPartRodCablingSide(const double rodPhi, const double phiSegmentWidth) {
  double phiSegmentStartOneCablingSide = computePhiSegmentStart(rodPhi, phiSegmentWidth, true);
  int phiSegmentRefOneCablingSide = computePhiSegmentRef(rodPhi, phiSegmentStartOneCablingSide, phiSegmentWidth, true);  
  bool isPositiveCablingSide = ((phiSegmentRefOneCablingSide % 2) == 1);
  return isPositiveCablingSide;
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


std::string ModulesToBundlesConnector::computeBundleType(const bool isBarrel, const std::string subDetectorName, const int layerDiskNumber, const int ringNumber) const {
  std::string bundleType;

  if (isBarrel) {
    if (subDetectorName == "TB2S") {
      bundleType = "2S";
    }
    else if (subDetectorName == "TBPS") {
      bundleType = (layerDiskNumber == 1 ? "PS10G" : "PS5G");
    }
  }

  else {
    if (subDetectorName == "TEDD_1") {
      if (ringNumber <= 4) bundleType = "PS10G";
      else if (ringNumber >= 5 && ringNumber <= 7) bundleType = "PS5GA";
      else if (ringNumber >= 8 && ringNumber <= 10) bundleType = "PS5GB";
      else if (ringNumber >= 11) bundleType = "2S";
    }

    else if (subDetectorName == "TEDD_2") {
      if (ringNumber <= 3) bundleType = "null";
      else if (ringNumber >= 4 && ringNumber <= 6) bundleType = "PS5GA";
      else if (ringNumber >= 7 && ringNumber <= 10) bundleType = "PS5GB";
      else if (ringNumber >= 11) bundleType = "2S";
    }
  }

  return bundleType;
}


int ModulesToBundlesConnector::computeBundleTypeIndex(const bool isBarrel, const std::string bundleType, const int totalNumFlatRings, const int maxNumModulesPerBundle, const bool isTilted, const bool isExtraFlatPart) const {
  int bundleTypeIndex;
  if (isBarrel) {
    if (bundleType == "2S") bundleTypeIndex = 0;
    else {
      if (!isTilted) {
	if (totalNumFlatRings <= maxNumModulesPerBundle) bundleTypeIndex = 1;
	else {
	  if (!isExtraFlatPart) bundleTypeIndex = 1;
	  else bundleTypeIndex = 2;
	}
      } 
      else bundleTypeIndex = 0;
    }
  }
  else {
    if (bundleType == "PS10G") bundleTypeIndex = 0;
    else if (bundleType == "PS5GA") bundleTypeIndex = 1;
    else if (bundleType == "PS5GB") bundleTypeIndex = 2;
    else if (bundleType == "2S") bundleTypeIndex = 3;
  }
  return bundleTypeIndex;
}


int ModulesToBundlesConnector::computeBundleId(const bool isBarrel, const bool isPositiveCablingSide, const int layerDiskNumber, const int phiRef, const int bundleTypeIndex) const {
  int cablingSideIndex = 0;
  if (isBarrel) {
    cablingSideIndex = (isPositiveCablingSide ? 1 : 3);
  }
  else {
    cablingSideIndex = (isPositiveCablingSide ? 2 : 4);
  }

  int bundleId = cablingSideIndex * 10000 + layerDiskNumber * 1000 + phiRef * 10 + bundleTypeIndex;
  return bundleId;
}


void ModulesToBundlesConnector::buildBundle(DetectorModule& m, std::map<int, Bundle*>& bundles, std::map<int, Bundle*>& negBundles, const int bundleId, const std::string bundleType, const std::string subDetectorName, const int layerDiskNumber, const double phiSegmentWidth, const int phiSegmentRef, const double phiRegionStart, const double phiRegionWidth, int phiRegionRef, const double phiSectorWidth, const int phiSectorRef, const bool isPositiveCablingSide, const bool isTiltedPart) {
  
  std::map<int, Bundle*>& stereoBundles = (isPositiveCablingSide ? bundles : negBundles);

  Bundle* bundle = nullptr;
  auto found = stereoBundles.find(bundleId);
  if (found == stereoBundles.end()) {
    bundle = createAndStoreBundle(bundles, negBundles, bundleId, bundleType, subDetectorName, layerDiskNumber, phiSegmentWidth, phiSegmentRef, phiRegionStart, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef, isPositiveCablingSide, isTiltedPart);
  }
  else {
    bundle = found->second;
  }

  connectModuleToBundle(m, bundle);
}


Bundle* ModulesToBundlesConnector::createAndStoreBundle(std::map<int, Bundle*>& bundles, std::map<int, Bundle*>& negBundles, const int bundleId, const std::string bundleType, const std::string subDetectorName, const int layerDiskNumber, const double phiSegmentWidth, const int phiSegmentRef, const double phiRegionStart, const double phiRegionWidth, const int phiRegionRef, const double phiSectorWidth, const int phiSectorRef, const bool isPositiveCablingSide, const bool isTiltedPart) {

  Bundle* bundle = GeometryFactory::make<Bundle>(bundleId, bundleType, subDetectorName, layerDiskNumber, phiSegmentWidth, phiSegmentRef, phiRegionStart, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef, isPositiveCablingSide, isTiltedPart);

  if (isPositiveCablingSide) {
    bundles.insert(std::make_pair(bundleId, bundle));
  }
  else {
    negBundles.insert(std::make_pair(bundleId, bundle));
  }
  return bundle;
}


void ModulesToBundlesConnector::connectModuleToBundle(DetectorModule& m, Bundle* bundle) const {
  bundle->addModule(&m);
  m.setBundle(bundle);
}
 

void ModulesToBundlesConnector::staggerModules(std::map<int, Bundle*>& bundles) {

  for (auto& b : bundles) {
    if (b.second->subDetectorName() == "TEDD_1" || b.second->subDetectorName() == "TEDD_2") {

      while (b.second->numModules() > maxNumModulesPerBundle_) {
	bool isPositiveCablingSide = b.second->isPositiveCablingSide();
	int diskNumber = b.second->layerDiskNumber();

	std::string bundleType = b.second->type();
	int bundleTypeIndex;
	if (bundleType == "PS10G") bundleTypeIndex = 0;
	else if (bundleType == "PS5GA") bundleTypeIndex = 1;
	else if (bundleType == "PS5GB") bundleTypeIndex = 2;
	else if (bundleType == "2S") bundleTypeIndex = 3;

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
	  nextBundleId = 20000 + diskNumber * 1000 + nextPhiRegionRef * 10 + bundleTypeIndex;
	  previousBundleId = 20000 + diskNumber * 1000 + previousPhiRegionRef * 10 + bundleTypeIndex;
	}
	else {
	  nextBundleId = 40000 + diskNumber * 1000 + nextPhiRegionRef * 10 + bundleTypeIndex;
	  previousBundleId = 40000 + diskNumber * 1000 + previousPhiRegionRef * 10 + bundleTypeIndex;
	}

	double minPhiBorder = fabs( femod((b.second->minPhi() - phiRegionStart), phiRegionWidth) );
	double maxPhiBorder = fabs( femod((b.second->maxPhi() - phiRegionStart), phiRegionWidth) - phiRegionWidth);
	      

	if (bundles.count(previousBundleId) != 0 && bundles.count(nextBundleId) != 0) {
	  // Cannot assign the extra module : both neighbouring phi regions are full !
	  if (bundles[previousBundleId]->numModules() >= maxNumModulesPerBundle_ && bundles[nextBundleId]->numModules() >= maxNumModulesPerBundle_) {
	    std::cout << "I am a refugee module in side " << isPositiveCablingSide << ", disk " << diskNumber << ", bundleType " << bundleType 
		      << ", phiRegionRef " << phiRegionRef << ", phiRegionWidth " << phiRegionWidth
		      << ", which has already more than " << maxNumModulesPerBundle_<< " modules, and none of my neighbouring regions wants to welcome me :/" 
		      << std::endl;
	    break;
	  }

	  // Assign a module to the next phi region
	  else if (bundles[previousBundleId]->numModules() >= maxNumModulesPerBundle_ || maxPhiBorder <= minPhiBorder) {
	    std::cout << "Removing module in side " << isPositiveCablingSide << ", disk " << diskNumber << ", bundleType " << bundleType 
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
	    std::cout << "Removing module in side " << isPositiveCablingSide << ", disk " << diskNumber << ", bundleType " << bundleType 
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
      std::cout << "Error while creating bundle. Bundle " << b.first << " has phiSegmentRef = " << b.second->phiSegmentRef() << ", phiRegionRef = " << b.second->phiRegionRef() << ", phiSectorRef = " << b.second->phiSectorRef() << ", bundleType = " << b.second->type() << std::endl;
    }
  }
}
