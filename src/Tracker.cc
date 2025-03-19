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
      s->buildInTracker(myid());
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
 * The Ids here have nothing to do with DetIds: they are the numbers which appear in the configuration files.
 */
void Tracker::addHierarchyInfoToModules() {
 
  class HierarchyInfoToModulesVisitor : public GeometryVisitor {  
  public:
    void visit(Barrel& b) { 
      subdetectorId_++; 
    }
    void visit(Endcap& e) { 
      subdetectorId_++; 
    }
    void visit(Layer& l)  {
      layerDiskId_ = l.myid();     
    }
    void visit(Disk& d)   {
      layerDiskId_ = d.myid();
    }
    void visit(RodPair& r){
      rodRingId_ = r.myid(); 
    }
    void visit(Ring& r) {
      rodRingId_ = r.myid();
      isRingOn4Dees_ = r.isRingOn4Dees();
      const int ringIncrement = (isRingOn4Dees_ ? 1 : 2);
      bigDeltaIndex_ = (r.isSmallerAbsZRingInDisk() ? 0 : ringIncrement);
    }
    void visit(Module& m) { 
      m.subdetectorId(subdetectorId_); 
    }
    void visit(BarrelModule& m) { 
      m.layer(layerDiskId_); 
      m.rod(rodRingId_); 
    }
    void visit(EndcapModule& m) { 
      m.disk(layerDiskId_);
      m.ring(rodRingId_);
      const int moduleIncrement = (isRingOn4Dees_? 2 : 1);
      int smallDeltaIndex = (m.isSmallerAbsZModuleInRing() ? 0 : moduleIncrement);
      // Get which disk surface the module belongs to.
      // Disk Surface 1 is the surface with the lowest |Z|.
      // Disk Surface 4 is the surface with the biggest |Z|.
      int diskSurfaceIndex = 1 + bigDeltaIndex_ + smallDeltaIndex;
      m.endcapDiskSurface(diskSurfaceIndex);
    }

  private:
    int subdetectorId_ = 0;
    int layerDiskId_;
    int rodRingId_;
    bool isRingOn4Dees_;
    int bigDeltaIndex_;
  };

  HierarchyInfoToModulesVisitor v;
  accept(v);
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

      std::sort(trackerDisks_.begin(), trackerDisks_.end(), [] (const Disk* d1, const Disk* d2) { return fabs(d1->centerZ()) < fabs(d2->centerZ()); });
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
    EndcapDetIdBuilder v(isPixelTracker(), hasSubDisks(), endcapScheme);
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
