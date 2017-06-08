#include "Tracker.hh"

std::pair<double, double> Tracker::computeMinMaxEta() const {
  double min = std::numeric_limits<double>::max(), max = 0;
  for (auto m : modules()) {
    min = MIN(min, m->minEta());
    max = MAX(max, m->maxEta());
  } 

  //return std::make_pair(-1*log(tan(min/2.)), -1*log(tan(max/2.)));
  //return std::make_pair(min, max);
  return std::make_pair(-4.0,4.0); // CUIDADO to make it equal to the extended pixel - make it better ASAP!!
}

void Tracker::build() {
  try {
    check();

    double barrelMaxZ = 0;

    // Build barrel(s)
    for (auto& mapel : barrelNode) {
      if (!containsOnly.empty() && containsOnly.count(mapel.first) == 0) continue;
      Barrel* b = GeometryFactory::make<Barrel>();
      b->myid(mapel.first);
      b->store(propertyTree());
      b->store(mapel.second);
      b->build();
      b->cutAtEta(etaCut());
      barrelMaxZ = MAX(b->maxZ(), barrelMaxZ);
      barrels_.push_back(b);
    }

    // Build endcap(s)
    for (auto& mapel : endcapNode) {
      if (!containsOnly.empty() && containsOnly.count(mapel.first) == 0) continue;
      Endcap* e = GeometryFactory::make<Endcap>();
      e->myid(mapel.first);
      e->barrelMaxZ(barrelMaxZ);
      e->store(propertyTree());
      e->store(mapel.second);
      e->build();
      e->cutAtEta(etaCut());
      endcaps_.push_back(e);
    }

    // Remove requested modules
    class ModuleRemover : public GeometryVisitor {
      void visit(RodPair& m_rodPair) { m_rodPair.removeModules(); }
      void visit(Ring& m_ring) { m_ring.removeModules(); }
    };
    ModuleRemover m_moduleRemover;
    accept(m_moduleRemover);

    // Build support structures within tracker
    for (auto& mapel : supportNode) {
      SupportStructure* s = new SupportStructure();
      s->store(propertyTree());
      s->store(mapel.second);
      s->buildInTracker();
      supportStructures_.push_back(s);
    }
  }
  catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  // Add all modules to a set, directly accessible from Tracker
  accept(moduleSetVisitor_);

  // Add geometry hierarchy information to modules
  addHierarchyInfoToModules();

  // Add Layers and Disks a global Tracker numbering.
  // All Tracker Layers are numbered by increasing radius.
  // All Tracker Disks are numbered by increasing fabs(Z). Numbering starts from 1 on (+Z) side, and from 1 on (-Z) side.
  addLayerDiskNumbers();
  
  // Build DetIds
  buildDetIds();

  cleanup();
  builtok(true);
}


/** 
 * Add geometry hierarchy information to modules.
 */
void Tracker::addHierarchyInfoToModules() {
 
  class HierarchicalNameVisitor : public GeometryVisitor {
    int cntId = 0;
    string cnt;
    int c1, c2;
  public:
    void visit(Barrel& b) { cnt = b.myid(); cntId++; }
    void visit(Endcap& e) { cnt = e.myid(); cntId++; }
    void visit(Layer& l)  { c1 = l.myid(); }
    void visit(Disk& d)   { c1 = d.myid(); }
    void visit(RodPair& r){ c2 = r.myid(); }
    void visit(Ring& r)   { c2 = r.myid(); }
    void visit(Module& m) { m.cntNameId(cnt, cntId); }
    void visit(BarrelModule& m) { m.layer(c1); m.rod(c2); }
    void visit(EndcapModule& m) { m.disk(c1); m.ring(c2); }
  };

  HierarchicalNameVisitor cntNameVisitor;
  accept(cntNameVisitor);
}


/**
 * Add Layers and Disks a global Tracker numbering.
 * All Tracker Layers are numbered by increasing radius.
 * All Tracker Disks are numbered by increasing fabs(Z). Numbering starts from 1 on (+Z) side, and from 1 on (-Z) side.
 */
void Tracker::addLayerDiskNumbers() {

  class LayerDiskNumberBuilder : public GeometryVisitor {
  public:
    void visit(Layer& l)  { trackerLayers_.push_back(&l); }
    void visit(Disk& d)   { trackerDisks_.push_back(&d); }

    void postVisit() {
      std::sort(trackerLayers_.begin(), trackerLayers_.end(), [] (const Layer* l1, const Layer* l2) { return l1->minR() < l2->minR(); });
      int layerNumber = 1;
      for (auto& l : trackerLayers_) {
	l->layerNumber(layerNumber); 
	layerNumber++; 
      }

      std::sort(trackerDisks_.begin(), trackerDisks_.end(), [] (const Disk* d1, const Disk* d2) { return fabs(d1->averageZ()) < fabs(d2->averageZ()); });
      int positiveZDiskNumber = 1;
      int negativeZDiskNumber = 1;
      for (auto& d : trackerDisks_) {
	if (d->side()) { 
	  d->diskNumber(positiveZDiskNumber); 
	  positiveZDiskNumber++; 
	} 
	else { 
	  d->diskNumber(negativeZDiskNumber); 
	  negativeZDiskNumber++; 
	} 
      }
    }
  private:
    std::vector<Layer*> trackerLayers_;
    std::vector<Disk*> trackerDisks_;
    //PtrVector<Layer> trackerLayers_;  // Would be good to use this
    //PtrVector<Disk> trackerDisks_;
  };

  LayerDiskNumberBuilder builder;
  accept(builder);
  builder.postVisit();
}


/** 
 * Build DetIds in the Tracker.
 */
void Tracker::buildDetIds() {
  // Barrel part
  std::vector<int> barrelScheme = mainConfigHandler::instance().getDetIdScheme(barrelDetIdScheme());
  if (barrelScheme.size() != 0) {
    BarrelDetIdBuilder v(isPixelTracker(), barrelScheme);
    accept(v);
  }
  // Endcap part
  std::vector<int> endcapScheme = mainConfigHandler::instance().getDetIdScheme(endcapDetIdScheme());
  if (endcapScheme.size() != 0) {
    EndcapDetIdBuilder v(isPixelTracker(), endcapScheme);
    accept(v);
  }
  checkDetIds();
}


/**
 * Check DetIds : DetIds should oviously be unique !
 * NB : If DetIds are replicated outside of the same Tracker, this will not be noticed. Though, this cannot happen by construction.
 */
void Tracker::checkDetIds() {
  std::set<int> moduleDetIds;
  std::set<int> sensorDetIds;

  const Modules& modules = moduleSetVisitor_.modules();

  for (const auto& m : modules) {
    // Check Modules DetIds are unique
    uint32_t myModuleDetId = m->myDetId();
    auto found = moduleDetIds.find(myModuleDetId);
    if (found != moduleDetIds.end() && myModuleDetId != 0) logWARNING(any2str(myid()) + ": Error while building DetIds !! DetId " + any2str(myModuleDetId) + " is duplicated.");
    else moduleDetIds.insert(myModuleDetId);

    // Check Sensors DetIds are unique
    for (const auto& s : m->sensors()) {
      uint32_t mySensorDetId = s.myDetId();
      auto found = sensorDetIds.find(mySensorDetId);
      if (found != sensorDetIds.end() && mySensorDetId != 0) logWARNING(any2str(myid()) + ": Error while building DetIds !! DetId " + any2str(mySensorDetId) + " is duplicated.");
      else sensorDetIds.insert(mySensorDetId);
    }
  }
}





void Tracker::buildCabling() {
  try {
    
    class ModulesToBundlesBuilder : public GeometryVisitor {  
      std::string barrelName;     
      int layerNumber;
      int numRods;
      int totalNumFlatRings;   

      std::string endcapName;
      int diskNumber;
      int ringNumber;
      int numModulesInRing;

      std::string type;
      int typeIndex;
      bool side;
   
      int phiSegmentRef;
      int negPhiSegmentRef;
      double phiRegionWidth;
      const double phiSectorWidth = 40. * M_PI / 180.;
  
      int bundleId;
      int bundleFlatId;   
      int bundleFlatIdB;      
      int bundleTiltedId;

      int negBundleId;
      int negBundleFlatId;
      int negBundleFlatIdB;
      int negBundleTiltedId;

      Bundle* bundle = NULL;
      Bundle* bundleFlat = NULL;
      Bundle* bundleFlatB = NULL;
      Bundle* bundleTilted = NULL;

      Bundle* negBundle = NULL;
      Bundle* negBundleFlat = NULL;
      Bundle* negBundleFlatB = NULL;
      Bundle* negBundleTilted = NULL;

      std::map<int, Bundle*> bundles_;
      std::map<int, Bundle*> negBundles_;

    public:
      void visit(Barrel& b) {
	barrelName = b.myid();	
      }

      void visit(Layer& l) {
	layerNumber = l.layerNumber();
	numRods = l.numRods();
	totalNumFlatRings = l.buildNumModulesFlat() * 2 - 1; // Total number of flat rings on both +Z side and -Z side

	if (layerNumber == 1 || layerNumber == 2 || layerNumber == 4) phiRegionWidth = 40. * M_PI / 180.;
	else phiRegionWidth = 20. * M_PI / 180.;
      }

      void visit(RodPair& r) {
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
      
     
      void visit(BarrelModule& m) {
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


      void visit(Endcap& e) { 
	endcapName = e.myid();
      }

      void visit(Disk& d) {
	diskNumber = d.diskNumber();
	side = d.side();
      }

      void visit(Ring& r)   { 
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


      void visit(EndcapModule& m) {
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


      int computePhiSliceRef(const double phi, const double phiSliceStart, const double phiSliceWidth) const {
	double phiSliceRefExact = femod(phi - phiSliceStart, 2.*M_PI) / phiSliceWidth;
	int phiSliceRef = 0;
	if (fabs((phiSliceRefExact - round(phiSliceRefExact))) < 0.0001) phiSliceRef = fabs(round(phiSliceRefExact));
	else phiSliceRef = std::floor(phiSliceRefExact);

	return phiSliceRef;
      }


      void postVisit() {
	// STAGGER MODULES
	staggerModules(bundles_);
	staggerModules(negBundles_);

	// CHECK
	checkModulesToBundlesConnections(bundles_);
	checkModulesToBundlesConnections(negBundles_);
      }


      void staggerModules(std::map<int, Bundle*>& bundles) {

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



      void checkModulesToBundlesConnections(std::map<int, Bundle*>& bundles) {
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


      std::map<int, Bundle*> getBundles() { return bundles_; }
      std::map<int, Bundle*> getNegBundles() { return negBundles_; }
    };


    // MODULES TO BUNDLES
    ModulesToBundlesBuilder bundlesBuilder;
    accept(bundlesBuilder);
    bundlesBuilder.postVisit();
    bundles_ = bundlesBuilder.getBundles();
    negBundles_ = bundlesBuilder.getNegBundles();


    // BUNDLES TO CABLES
    connectBundlesToCables(bundles_, cables_, DTCs_);
    connectBundlesToCables(negBundles_, negCables_, negDTCs_);
    checkBundlesToCablesConnections(cables_);
    checkBundlesToCablesConnections(negCables_);
  }

  catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }
}






void Tracker::connectBundlesToCables(std::map<int, Bundle*>& bundles, std::map<int, Cable*>& cables, std::map<const std::string, const DTC*>& DTCs) {

  // Used to stagger several bundles
  std::map<int, int> Layer3FlatPhiSectorsCounter;
  std::map<int, int> Layer3TiltedPhiSectorsCounter;
  std::map<int, int> Layer4PhiSectorsCounter;
  std::map<int, int> Layer5PhiSectorsCounter;
  std::map<int, int> Layer6PhiSectorsCounter;
 

  for (auto& b : bundles) {
    int phiSectorRef = b.second->phiSectorRef();
    int phiSectorRefCable = phiSectorRef;
    double phiSectorWidth = b.second->phiSectorWidth();

    int numPhiSectors = round(2 * M_PI / phiSectorWidth);
    int nextPhiSectorRef = femod( (phiSectorRef + 1), numPhiSectors);
    int previousPhiSectorRef = femod( (phiSectorRef - 1), numPhiSectors);

    std::string bundleType = b.second->type();
    std::string cableType = bundleType;
    if (cableType == "PS5GA" || cableType == "PS5GB") cableType = "PS5G";

    int cableTypeIndex;
    if (cableType == "PS10G") cableTypeIndex = 0;
    else if (cableType == "PS5G") cableTypeIndex = 1;
    else if (cableType == "2S") cableTypeIndex = 2;


    std::string subDetectorName = b.second->subDetectorName();
    int layerDiskNumber = b.second->layerDiskNumber();

    // Used to stagger several bundles
    int phiRegionRef = b.second->phiRegionRef();   

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
}


void Tracker::checkBundlesToCablesConnections(std::map<int, Cable*>& cables) {
  for (auto& c : cables) {
    if (c.second->numBundles() > 6) {
      std::cout << "There was an error while staggering bundles. Cable " 
		<< c.first << " is connected to " << c.second->numBundles() << " bundles." 
		<< std::endl;
    }

    if (c.second->phiSectorRef() <= -1) {
      std::cout << "Error while creating cable. Cable " << c.first << " has phiSectorRef = " << c.second->phiSectorRef() << ". type = " << c.second->type() << ", slot = " <<  c.second->slot() << std::endl;
    }

  }
}

