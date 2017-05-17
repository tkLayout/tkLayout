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

  accept(moduleSetVisitor_);

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
  } cntNameVisitor;

  accept(cntNameVisitor);


  // THis is uglyy, sorry i will completely rewrite this
  std::vector< std::pair<const Layer*, std::string> > sortedLayers;
  std::string barrelId;
  for (auto& b : barrels_) {
    barrelId = b.myid();
    for (const auto& l : b.layers()) { sortedLayers.push_back(std::make_pair(&l, barrelId)); }
  }
  std::sort(sortedLayers.begin(), sortedLayers.end(), [&] (std::pair<const Layer*, std::string> l1, std::pair<const Layer*, std::string> l2) { return l1.first->minR() < l2.first->minR(); });
  std::map< std::pair<std::string, int>, int > sortedLayersIds;
  int i = 1;
  for (const auto& l : sortedLayers) { sortedLayersIds.insert(std::make_pair(std::make_pair(l.second, l.first->myid()), i)); i++; }

  std::vector< std::pair<const Disk*, std::string> > sortedZPlusDisks;
  std::vector< std::pair<const Disk*, std::string> > sortedZMinusDisks;
  std::string endcapId;
  for (auto& e : endcaps_) {
    endcapId = e.myid();
    for (const auto& d : e.disks()) {
      if (d.side()) sortedZPlusDisks.push_back(std::make_pair(&d, endcapId));
      else sortedZMinusDisks.push_back(std::make_pair(&d, endcapId));
    }
  }
  std::sort(sortedZPlusDisks.begin(), sortedZPlusDisks.end(), [&] (std::pair<const Disk*, std::string> d1, std::pair<const Disk*, std::string> d2) { return d1.first->averageZ() < d2.first->averageZ(); });
  std::sort(sortedZMinusDisks.begin(), sortedZMinusDisks.end(), [&] (std::pair<const Disk*, std::string> d1, std::pair<const Disk*, std::string> d2) { return fabs(d1.first->averageZ()) < fabs(d2.first->averageZ()); });
  std::map< std::tuple<std::string, int, bool>, int > sortedDisksIds;
  i = 1;
  for (const auto& d : sortedZPlusDisks) { sortedDisksIds.insert(std::make_pair(std::make_tuple(d.second, d.first->myid(), d.first->side()), i)); i++; }
  i = 1;
  for (const auto& d : sortedZMinusDisks) { sortedDisksIds.insert(std::make_pair(std::make_tuple(d.second, d.first->myid(), d.first->side()), i)); i++; }
  // END OF THE UGLY PART





  class BarrelDetIdBuilder : public SensorGeometryVisitor {
  private:
    std::string schemeName;
    std::vector<int> schemeShifts;
    std::map< std::pair<std::string, int>, int > sortedLayersIds;

    std::map<int, uint32_t> detIdRefs;

    std::string barrelName;
    bool isTiltedLayer;
    int numRods;
    int numFlatRings;
    int numRings;
    bool isCentered;
    uint32_t phiRef;

  public:
    BarrelDetIdBuilder(std::string name, std::vector<int> shifts, std::map< std::pair<std::string, int>, int > layersIds) : schemeName(name), schemeShifts(shifts), sortedLayersIds(layersIds) {}

    void visit(Barrel& b) {
      barrelName = b.myid();

      detIdRefs[0] = 1;

      detIdRefs[1] = 205 % 100;
			   
      detIdRefs[2] = 0;
    }

    void visit(Layer& l) {
      std::pair<std::string, int> layerId;
      layerId = std::make_pair(barrelName, l.myid());
      l.layerNumber(sortedLayersIds.at(layerId));

      detIdRefs[3] = l.layerNumber();

      isTiltedLayer = l.isTilted();
      numRods = l.numRods();
      numFlatRings = l.buildNumModulesFlat();
      numRings = l.buildNumModules();
      if (!isTiltedLayer) numRings = numFlatRings;
    }

    void visit(RodPair& r) {
      double startAngle = femod( r.Phi(), (2 * M_PI / numRods));
      phiRef = 1 + round(femod(r.Phi() - startAngle, 2*M_PI) / (2*M_PI) * numRods);

      isCentered = (r.startZMode() == RodPair::StartZMode::MODULECENTER);
    }

    void visit(BarrelModule& m) {
      int side = m.uniRef().side;
      uint32_t ringRef;

      if (!m.isTilted()) {
	detIdRefs[4] = 3;
	detIdRefs[5] = phiRef;

	if (isCentered) ringRef = (side > 0 ? (m.uniRef().ring + numFlatRings - 1) : (1 + numFlatRings - m.uniRef().ring));
	else ringRef = (side > 0 ? (m.uniRef().ring + numFlatRings) : (1 + numFlatRings - m.uniRef().ring));
	detIdRefs[6] = ringRef;
      }

      else {
	uint32_t category = (side > 0 ? 2 : 1);
	detIdRefs[4] = category;

	ringRef = (side > 0 ? (m.uniRef().ring - numFlatRings) : (1 + numRings - m.uniRef().ring));
	detIdRefs[5] = ringRef;

	detIdRefs[6] = phiRef;
      }

      uint32_t sensorRef = 0;
      detIdRefs[7] = sensorRef;

      m.buildDetId(detIdRefs, schemeShifts);

      //modules_.insert(std::make_pair(m.myDetId(), m));
    }

    void visit(Sensor& s) {
      uint32_t sensorRef = (s.innerOuter() == SensorPosition::LOWER ? 1 : 2);
      if (s.subdet() == ModuleSubdetector::BARREL) {
	detIdRefs[7] = sensorRef;

	s.buildDetId(detIdRefs, schemeShifts);
      }  
    }

  };


  class EndcapDetIdBuilder : public SensorGeometryVisitor {
  private:
    std::string schemeName;
    std::vector<int> schemeShifts;
    std::map< std::tuple<std::string, int, bool>, int > sortedDisksIds;

    std::map<int, uint32_t> detIdRefs;

    std::string endcapName;
    int numEmptyRings;
    int numModules;

  public:
    EndcapDetIdBuilder(std::string name, std::vector<int> shifts, std::map< std::tuple<std::string, int, bool>, int > disksIds) : schemeName(name), schemeShifts(shifts), sortedDisksIds(disksIds) {}

    void visit(Endcap& e) {
      endcapName = e.myid();

      detIdRefs[0] = 1;

      detIdRefs[1] = 204 % 100;
    }

    void visit(Disk& d) {
      bool side = d.side();

      std::tuple<std::string, int, bool> diskId;
      diskId = std::make_tuple(endcapName, d.myid(), side);
      d.diskNumber(sortedDisksIds.at(diskId));  

      uint32_t sideRef = (side ? 2 : 1);
      detIdRefs[2] = sideRef;

      detIdRefs[3] = 0;

      uint32_t diskRef = d.diskNumber();
      detIdRefs[4] = diskRef;

      numEmptyRings = d.numEmptyRings();
    }
   
    void visit(Ring& r) {
      uint32_t ringRef = r.myid() - numEmptyRings;
      detIdRefs[5] = ringRef;

      detIdRefs[6] = 1;

      numModules = r.numModules();
    }
 
    void visit(EndcapModule& m) {
      double startAngle = femod( m.center().Phi(), (2 * M_PI / numModules));
      uint32_t phiRef = 1 + round(femod(m.center().Phi() - startAngle, 2*M_PI) / (2*M_PI) * numModules);
      detIdRefs[7] = phiRef;

      uint32_t sensorRef = 0;
      detIdRefs[8] = sensorRef;
      m.buildDetId(detIdRefs, schemeShifts);

      //modules_.insert(std::make_pair(m.myDetId(), m));
      //std::cout << "disk = " << m.uniRef().layer << "ring = " <<  m.uniRef().ring << "side = " << m.uniRef().side << std::endl;
    }

    void visit(Sensor& s) {
      uint32_t sensorRef = (s.innerOuter() == SensorPosition::LOWER ? 1 : 2);
      if (s.subdet() == ModuleSubdetector::ENDCAP) {
	detIdRefs[8] = sensorRef;


	/*for (int a = 0; a < detIdRefs.size(); a++) {
	  std::cout << "values = " << std::endl;
	  std::cout << detIdRefs.at(a) << std::endl;
	  std::cout << "scheme = " << std::endl;
	  std::cout << schemeShifts.at(a) << std::endl;
	  }*/

	s.buildDetId(detIdRefs, schemeShifts);

	//std::bitset<32> test(s.myDetId());
	//std::cout << s.myDetId() << " " << test << " " << "rho = " <<  s.hitPoly().getCenter().Rho() << " z = " <<  s.hitPoly().getCenter().Z() << " phi = " <<  (s.hitPoly().getCenter().Phi() * 180. / M_PI) << std::endl;
      }
    }

  };


  if (!isPixelTracker()) {
    if (barrelDetIdScheme() == "Phase2Subdetector5" && detIdSchemes_.count("Phase2Subdetector5") != 0) {
      BarrelDetIdBuilder v(barrelDetIdScheme(), detIdSchemes_["Phase2Subdetector5"], sortedLayersIds);
      accept(v);
    }
    else logWARNING("barrelDetIdScheme = " + barrelDetIdScheme() + ". This barrel detId scheme is empty or incorrect or not currently implemented within tkLayout. No detId for Barrel sensors will be calculated.");

    if (endcapDetIdScheme() == "Phase2Subdetector4" && detIdSchemes_.count("Phase2Subdetector4") != 0) {
      EndcapDetIdBuilder v(endcapDetIdScheme(), detIdSchemes_["Phase2Subdetector4"], sortedDisksIds);
      accept(v);
    }
    else logWARNING("endcapDetIdScheme = " + endcapDetIdScheme() + ". This endcap detId scheme is empty or incorrect or not currently implemented within tkLayout. No detId for Endcap sensors will be calculated.");
  }


  cleanup();
  builtok(true);
}


std::map<std::string, std::vector<int> > Tracker::detIdSchemes() {
  std::map<std::string, std::vector<int> > schemes;

  std::ifstream schemesStream(mainConfigHandler::instance().getDetIdSchemesDirectory() + "/" + insur::default_detidschemesfile);

  if (schemesStream.good()) {
    std::string line;
    while(getline(schemesStream, line).good()) {
      if (line.empty()) continue;
      auto schemeData = split(line, " ");
      std::string schemeName = schemeData.at(0);
      if (schemeData.size() < 2) logWARNING("DetId scheme " + schemeName + " : no data was entered." );
      else {
	std::vector<int> detIdShifts;
	int sum = 0;
	for (int i = 1; i < schemeData.size(); i++) {
	  int shift = str2any<int>(schemeData.at(i));
	  detIdShifts.push_back(shift);
	  sum += shift; 
	}
	if (sum != 32) logWARNING("DetId scheme " + schemeName + " : The sum of Det Id shifts is not equal to 32." );
	else schemes.insert(std::make_pair(schemeName, detIdShifts));
      }
    }
    schemesStream.close();
  } else logWARNING("No file defining a detId scheme has been found.");

  return schemes;
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
      int numPhiRegions;
      const double phiSectorWidth = 40. * M_PI / 180.;
      const int numPhiSectors = round(2 * M_PI / phiSectorWidth);
  
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
	numPhiRegions = round(2 * M_PI / phiRegionWidth);
      }

      void visit(RodPair& r) {
	double rodPhi = r.Phi();

	double phiSegmentWidth = (2.*M_PI) / numRods;

	// Positive cabling side
	double phiSegmentStart = femod( rodPhi, phiSegmentWidth);
	phiSegmentRef = round(femod(rodPhi - phiSegmentStart, 2.*M_PI) / phiSegmentWidth);
	
	double phiRegionStart = 0.;
	double phiRegionRefExact = femod(rodPhi - phiRegionStart, 2.*M_PI) / phiRegionWidth;
	int phiRegionRef = 0;
	if (fabs((phiRegionRefExact - round(phiRegionRefExact))) < 0.0001) phiRegionRef = fabs(round(phiRegionRefExact));
	else phiRegionRef = std::floor(phiRegionRefExact);

	double phiSectorRefExact = femod(rodPhi, 2.*M_PI) / phiSectorWidth;
	int phiSectorRef = 0;
	if (fabs((phiSectorRefExact - round(phiSectorRefExact))) < 0.0001) phiSectorRef = fabs(round(phiSectorRefExact));
	else phiSectorRef = std::floor(phiSectorRefExact);

	// Negative cabling side
	double negPhiSegmentStart = femod( M_PI - rodPhi, phiSegmentWidth);
	negPhiSegmentRef = round(femod(M_PI - rodPhi - negPhiSegmentStart, 2.*M_PI) / phiSegmentWidth);
	
	double negPhiRegionStart = 0.;
	double negPhiRegionRefExact = femod(M_PI - rodPhi - negPhiRegionStart, 2.*M_PI) / phiRegionWidth;
	int negPhiRegionRef = 0;
	if (fabs((negPhiRegionRefExact - round(negPhiRegionRefExact))) < 0.0001) negPhiRegionRef = fabs(round(negPhiRegionRefExact));
	else negPhiRegionRef = std::floor(negPhiRegionRefExact);

	double negPhiSectorRefExact = femod(M_PI - rodPhi, 2.*M_PI) / phiSectorWidth;
	int negPhiSectorRef = 0;
	if (fabs((negPhiSectorRefExact - round(negPhiSectorRefExact))) < 0.0001) negPhiSectorRef = fabs(round(negPhiSectorRefExact));
	else negPhiSectorRef = std::floor(negPhiSectorRefExact);

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
	if (side) {
	  bool isPositiveCablingSide = true;
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

	  double phiSegmentWidth = (2.*M_PI) / numModulesInRing;
	  double phiSegmentStart = femod( m.center().Phi(), phiSegmentWidth);
	  phiSegmentRef = round(femod(m.center().Phi() - phiSegmentStart, 2.*M_PI) / phiSegmentWidth);

	  int phiRegionRef = std::floor(femod(m.center().Phi() - phiRegionStart, 2.*M_PI) / phiRegionWidth);	  
	  bundleId = 20000 + diskNumber * 1000 + phiRegionRef * 10 + typeIndex;
	  //else { bundleId = 30000 + diskNumber * 1000 + phiRegionRef * 10 + typeIndex; }

	  int phiSectorRef = std::floor(femod(m.center().Phi(), 2.*M_PI) / phiSectorWidth);

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
      }



      void postVisit() {
	// STAGGER MODULES
	staggerModules(bundles_);

	// CHECK
	checkModulesToBundlesConnections(bundles_);
	checkModulesToBundlesConnections(negBundles_);
      }


      void staggerModules(std::map<int, Bundle*>& bundles) {

	for (auto& b : bundles) {
	  if (b.second->subDetectorName() == "TEDD_1" || b.second->subDetectorName() == "TEDD_2") {

	    while (b.second->numModules() > 12) {
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
	      int nextBundleId = 20000 + diskNumber * 1000 + nextPhiRegionRef * 10 + typeIndex;
	      int previousBundleId = 20000 + diskNumber * 1000 + previousPhiRegionRef * 10 + typeIndex;

	      double minPhiBorder = fabs( femod((b.second->minPhi() - phiRegionStart), phiRegionWidth) );
	      double maxPhiBorder = fabs( femod((b.second->maxPhi() - phiRegionStart), phiRegionWidth) - phiRegionWidth);


	      if (bundles.count(previousBundleId) != 0 && bundles.count(nextBundleId) != 0) {
		// Cannot assign the extra module : both neighbouring phi regions are full !
		if (bundles[previousBundleId]->numModules() >= 12 && bundles[nextBundleId]->numModules() >= 12) {
		  std::cout << "I am a refugee module from disk " << diskNumber << ", type " << type 
			    << ", phiRegionRef " << phiRegionRef << ", phiRegionWidth " << phiRegionWidth
			    << ", which has already more than 12 modules, and none of my neighbouring regions wants to welcome me :/" 
			    << std::endl;
		  break;
		}

		// Assign the module with the biggest phi to the next phi region
		else if (bundles[previousBundleId]->numModules() >= 12 || maxPhiBorder <= minPhiBorder) {
		  std::cout << "Removing module in disk " << diskNumber << ", type " << type 
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

		// Assign the module with the lowest phi to the previous phi region
		else if (bundles[nextBundleId]->numModules() >= 12 || minPhiBorder < maxPhiBorder) {
		  std::cout << "Removing module in disk " << diskNumber << ", type " << type 
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
    if (!isPositiveCablingSide) cableId *= -1;

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

