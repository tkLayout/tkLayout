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
    
    class ModulesToRibbonsBuilder : public GeometryVisitor {  
      std::string barrelName;     
      int layerNumber;
      int numRods;
      int numFlatRings;   

      std::string endcapName;
      int diskNumber;
      int ringNumber;

      std::string type;
      int typeIndex;
      int phiRegionRef;
      int phiSelect;
      int side;

      const double phiSectorWidth = 40. * M_PI / 180.;
  
      int ribbonId;
      int ribbonFlatId;
      int ribbonFlatIdB;
      int ribbonTiltedId;

      std::map<int, Ribbon*> ribbons_;

    public:
      void visit(Barrel& b) {
	barrelName = b.myid();	
      }

      void visit(Layer& l) {
	layerNumber = l.layerNumber();
	numRods = l.numRods();

	numFlatRings = l.buildNumModulesFlat();
      }

      void visit(RodPair& r) {
	double startPhi = femod( r.Phi(), (2 * M_PI / numRods));
	double phiRegionWidth = (2*M_PI) / numRods;
	phiRegionRef = round(femod(r.Phi() - startPhi, 2*M_PI) / phiRegionWidth);
	int phiSectorRef = round(femod(r.Phi() - startPhi, 2*M_PI) / phiSectorWidth);	

	// Create 2S ribbons
	if (barrelName == "TB2S") {
	  ribbonId = 10000 + layerNumber * 1000 + phiRegionRef * 10;
	  type = "2S";
	  Ribbon ribbon(ribbonId, type, barrelName, layerNumber, startPhi, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef);
	  ribbons_.insert(std::make_pair(ribbonId, &ribbon));
	}

	// Create PS ribbons
	else if (barrelName == "TBPS") {
	  type = (layerNumber == 1 ? "PS10G" : "PS5G");

	  // Flat part
	  phiSelect = layerNumber % 2;
	  // standard case
	  if ( (phiRegionRef + 1) % 2 == phiSelect) {
	    ribbonFlatId = 10000 + layerNumber * 1000 + phiRegionRef * 10 + 1;
	    Ribbon ribbonFlat(ribbonFlatId, type, barrelName, layerNumber, startPhi, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef);
	    ribbons_.insert(std::make_pair(ribbonFlatId, &ribbonFlat));

	    // For layer 3, need to add a second ribbon for flat part
	    if (numFlatRings > 12) {
	      ribbonFlatIdB = 10000 + layerNumber * 1000 + phiRegionRef * 10 + 2;
	      Ribbon ribbonFlatB(ribbonFlatIdB, type, barrelName, layerNumber, startPhi, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef);
	      ribbons_.insert(std::make_pair(ribbonFlatIdB, &ribbonFlatB));
	    }
	  }

	  // Tilted part
	  ribbonTiltedId = 10000 + layerNumber * 1000 + phiRegionRef * 10;	  
	  Ribbon ribbonTilted(ribbonTiltedId, type, barrelName, layerNumber, startPhi, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef);
	  ribbons_.insert(std::make_pair(ribbonTiltedId, &ribbonTilted));	 
	}
      }
      
      void visit(Module& m) {
	side = m.uniRef().side;
      }
     
      void visit(BarrelModule& m) {

	// Connect modules to 2S ribbons
	if (barrelName == "TB2S") {
	  if (side > 0) {
	    m.setRibbonId(ribbonId);
	    ribbons_[ribbonId]->addModule(&m);
	  }
	}

	// Connect modules to PS ribbons
	else if (barrelName == "TBPS") {

	  // flat modules
	  if (!m.isTilted() && ( (phiRegionRef + 1) % 2 == phiSelect)) {
	    // standard case
	    if (numFlatRings <= 12) {
	      m.setRibbonId(ribbonFlatId);
	      ribbons_[ribbonFlatId]->addModule(&m);
	    }
	    // For layer 3, need to add a second ribbon for flat part
	    else {
	      if (side > 0) {
		m.setRibbonId(ribbonFlatId);
		ribbons_[ribbonFlatId]->addModule(&m);
	      } else {
		m.setRibbonId(ribbonFlatIdB);
		ribbons_[ribbonFlatIdB]->addModule(&m);
	      }
	    }
	  }

	  // tilted modules
	  else if (m.isTilted() && side > 0) {
	    m.setRibbonId(ribbonTiltedId);
	    ribbons_[ribbonTiltedId]->addModule(&m);
	  }

	}
      }





      void visit(Endcap& e) { 
	endcapName = e.myid();
      }

      void visit(Disk& d) {
	diskNumber = d.diskNumber();
      }

      void visit(Ring& r)   { 
	ringNumber = r.myid();

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
	double startPhi;
	double phiRegionWidth;

	if (type == "PS10G" || type == "PS5GA") {
	  startPhi = 0.;
	  phiRegionWidth = 40. * M_PI / 180.;
	}

	else if (type == "PS5GB") {
	  startPhi = 0.;
	  phiRegionWidth = 20. * M_PI / 180.;
	}

	else if (type == "2S") {
	  startPhi = 0.;
	  phiRegionWidth = 360. / 27. * M_PI / 180.;
	}

	phiRegionRef = round(femod(m.center().Phi() - startPhi, 2*M_PI) / phiRegionWidth);	  
	ribbonId = 20000 + diskNumber * 1000 + phiRegionRef * 10 + typeIndex;

	int phiSectorRef = round(femod(m.center().Phi() - startPhi, 2*M_PI) / phiSectorWidth);

	if (ribbons_.count(ribbonId) == 0) {
	  Ribbon ribbon(ribbonId, type, endcapName, diskNumber, startPhi, phiRegionWidth, phiRegionRef, phiSectorWidth, phiSectorRef);
	  ribbon.addModule(&m);
	  ribbons_.insert(std::make_pair(ribbonId, &ribbon));
	}
	else ribbons_[ribbonId]->addModule(&m);	 
	m.setRibbonId(ribbonId);
      }





      void postVisit() {
	for (auto& r : ribbons_) {
	  if (r.second->subDetectorName() == "TEDD_1" || r.second->subDetectorName() == "TEDD_2") {

	    while (r.second->numModules() > 12) {
	      int diskNumber = r.second->layerDiskNumber();

	      std::string type = r.second->type();
	      if (type == "PS10G") typeIndex = 0;
	      else if (type == "PS5GA") typeIndex = 1;
	      else if (type == "PS5GB") typeIndex = 2;
	      else if (type == "2S") typeIndex = 3;

	      double startPhi = r.second->startPhi();
	      double phiRegionWidth = r.second->phiRegionWidth();
	      int numPhiRegions = round(2 * M_PI / phiRegionWidth);

	      int phiRegionRef = r.second->phiRegionRef();
	      int nextPhiRegionRef = femod( (phiRegionRef + 1), numPhiRegions);
	      int previousPhiRegionRef = femod( (phiRegionRef - 1), numPhiRegions);

	      int ribbonId = r.first;
	      int nextRibbonId = 20000 + diskNumber * 1000 + nextPhiRegionRef * 10 + typeIndex;
	      int previousRibbonId = 20000 + diskNumber * 1000 + previousPhiRegionRef * 10 + typeIndex;


	      

	      double minPhiBorder = fabs( femod(r.second->minPhi(), phiRegionWidth) - startPhi );
	      double maxPhiBorder = fabs( femod(r.second->maxPhi(), phiRegionWidth) - phiRegionWidth);


	      if (ribbons_[previousRibbonId]->numModules() > 12 && ribbons_[nextRibbonId]->numModules() > 12) {
		std::cout << "I am a refugee module from disk " << diskNumber << ", phiRegionRef " << phiRegionRef 
			  << ", which has already more than 12 modules, and none of my neighbouring regions wants to welcome me :/" 
			  << std::endl;
	      }

	      // Assign the module with the biggest phi to the next phi region
	      else if (ribbons_[previousRibbonId]->numModules() > 12 || maxPhiBorder <= minPhiBorder) {

		Module* maxPhiMod = r.second->maxPhiModule();
		ribbons_[ribbonId]->removeModule(maxPhiMod);
		ribbons_[nextRibbonId]->addModule(maxPhiMod);
		maxPhiMod->setRibbonId(nextRibbonId);  // !!!!!!! ERROR : obviously doesnt change the ribbonId in the tracker modules
	      }

	      // Assign the module with the lowest phi to the previous phi region
	      else if (ribbons_[nextRibbonId]->numModules() > 12 || minPhiBorder < maxPhiBorder) {

		Module* minPhiMod = r.second->minPhiModule();
		ribbons_[ribbonId]->removeModule(minPhiMod);
		ribbons_[previousRibbonId]->addModule(minPhiMod);
		minPhiMod->setRibbonId(previousRibbonId); // !!!!!!!!!! ERROR : obviously doesnt change the ribbonId in the tracker modules
	      }

	    }
	  }
	}



	for (auto& r : ribbons_) {
	  if (r.second->numModules() > 12) {
	    std::cout << "There was an error while staggering ribbons. Rbbon " 
		      << r.first << " is connected to " << r.second->numModules() << " modules." 
		      << std::endl;
	  }
	}



      }




      std::map<int, Ribbon*> getRibbons() { return ribbons_; }
    };



 ModulesToRibbonsBuilder ribbonsBuilder;
 accept(ribbonsBuilder);
 ribbonsBuilder.postVisit();
 ribbons_ = ribbonsBuilder.getRibbons();






    }

  catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }
}
