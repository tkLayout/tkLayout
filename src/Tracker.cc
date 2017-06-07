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

  // Add modules to a set, directly accessible from Tracker
  accept(moduleSetVisitor_);

  // Add geometry hierarchy information to modules
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





  

  


  if (detIdSchemes_.count(barrelDetIdScheme()) != 0) {
    BarrelDetIdBuilder v(isPixelTracker(), detIdSchemes_[barrelDetIdScheme()], sortedLayersIds);
    accept(v);
  }
  else logWARNING("barrelDetIdScheme = " + barrelDetIdScheme() + ". This barrel detId scheme is empty or incorrect or not currently implemented within tkLayout. No detId for Barrel sensors will be calculated.");

  if (detIdSchemes_.count(endcapDetIdScheme()) != 0) {
    EndcapDetIdBuilder v(isPixelTracker(), detIdSchemes_[endcapDetIdScheme()], sortedDisksIds);
    accept(v);
  }
  else logWARNING("endcapDetIdScheme = " + endcapDetIdScheme() + ". This endcap detId scheme is empty or incorrect or not currently implemented within tkLayout. No detId for Endcap sensors will be calculated.");


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
	std::vector<int> geometryHierarchySizes;
	int sum = 0;
	for (int i = 1; i < schemeData.size(); i++) {
	  int size = str2any<int>(schemeData.at(i));
	  geometryHierarchySizes.push_back(size);
	  sum += size; 
	}
	if (sum != 32) logWARNING("DetId scheme " + schemeName + " : The sum of geometry hierarchy sizes is not equal to 32." );
	else schemes.insert(std::make_pair(schemeName, geometryHierarchySizes));
      }
    }
    schemesStream.close();
  } else logWARNING("No file defining a detId scheme has been found.");

  return schemes;
}
