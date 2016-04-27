#include "Tracker.h"

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

  cleanup();
  builtok(true);
}
