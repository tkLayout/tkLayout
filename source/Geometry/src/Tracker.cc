#include "Tracker.h"
#include "global_constants.h"

#include <Barrel.h>
#include <DetectorModule.h>
#include <Disk.h>
#include <Endcap.h>
#include <Layer.h>
#include <Ring.h>
#include <RodPair.h>
#include "SimParms.h"
#include <SupportStructure.h>

//
// Helper class: Name visitor - methods
//
void HierarchicalNameVisitor::visit(Barrel& b)        { cnt = b.myid(); cntId++; }
void HierarchicalNameVisitor::visit(Endcap& e)        { cnt = e.myid(); cntId++; }
void HierarchicalNameVisitor::visit(Layer& l)         { c1 = l.myid(); }
void HierarchicalNameVisitor::visit(Disk& d)          { c1 = d.myid(); }
void HierarchicalNameVisitor::visit(RodPair& r)       { c2 = r.myid(); }
void HierarchicalNameVisitor::visit(Ring& r)          { c2 = r.myid(); }
void HierarchicalNameVisitor::visit(DetectorModule& m){ m.cntNameId(cnt, cntId); }
void HierarchicalNameVisitor::visit(BarrelModule& m)  { m.layer(c1); m.rod(c2); }
void HierarchicalNameVisitor::visit(EndcapModule& m)  { m.disk(c1); m.ring(c2); }

//
// Constructor - parse geometry config file using boost property tree & read-in Barrel, Endcap & Support nodes
//
Tracker::Tracker(const PropertyTree& treeProperty) :
 minR(     string("minR")           ),
 maxR(     string("maxR")           ),
 maxZ(     string("maxZ")           ),
 minEta(   string("minEta")         ),
 maxEta(   string("maxEta")         ),
 etaCut(          "etaCut"          , parsedOnly(), SimParms::getInstance().getMaxEtaCoverage()),
 isPixelType(     "isPixelType"     , parsedOnly(), true),
 servicesForcedUp("servicesForcedUp", parsedOnly(), true),
 skipAllServices( "skipAllServices" , parsedOnly(), false),
 skipAllSupports( "skipAllSupports" , parsedOnly(), false),
 m_barrelNode(    "Barrel"          , parsedOnly()),
 m_endcapNode(    "Endcap"          , parsedOnly()),
 m_supportNode(   "Support"         , parsedOnly()),
 m_containsOnly(  "containsOnly"    , parsedOnly()),
 m_modulesSetVisitor(nullptr),
 m_cntNameVisitor(nullptr)
{
  // Set the geometry config parameters
  this->myid(treeProperty.data());
  this->store(treeProperty);

  m_modulesSetVisitor = new ModulesSetVisitor();
  m_cntNameVisitor    = new HierarchicalNameVisitor();
}

//
// Destructor
//
Tracker::~Tracker()
{
  if (m_modulesSetVisitor!=nullptr) delete m_modulesSetVisitor;
  if (m_cntNameVisitor   !=nullptr) delete m_cntNameVisitor;
}

//
// Build recursively individual subdetector systems: Barrels, Endcaps
//
void Tracker::build() {

  try {
    check();

    double barrelMaxZ = 0;

    // Build barrel tracker
    for (auto& iBarrel : m_barrelNode) {
      if (!m_containsOnly.empty() && m_containsOnly.count(iBarrel.first) == 0) continue;

      std::cout << "  " << iBarrel.first << std::endl;
      Barrel* barrel = GeometryFactory::make<Barrel>(iBarrel.first, iBarrel.second, propertyTree());
      barrel->build();
      barrel->cutAtEta(etaCut());

      m_barrels.push_back(barrel);

      // Calculate barrel maxZ position
      barrelMaxZ = MAX(barrel->maxZ(), barrelMaxZ);
    }

    // Build end-cap tracker
    for (auto& iEndcap : m_endcapNode) {
      if (!m_containsOnly.empty() && m_containsOnly.count(iEndcap.first) == 0) continue;

      std::cout << "  " << iEndcap.first << std::endl;
      Endcap* endcap = GeometryFactory::make<Endcap>(barrelMaxZ, iEndcap.first, iEndcap.second, propertyTree());
      endcap->build();
      endcap->cutAtEta(etaCut());
      m_endcaps.push_back(endcap);
    }

    // Build support structures within tracker
    for (auto& iSupport : m_supportNode) {

      std::cout << "  " << iSupport.first << std::endl;
      SupportStructure* support = new SupportStructure();
      support->store(propertyTree());
      support->store(iSupport.second);
      support->buildInTracker();
      m_supportStructures.push_back(support);
    }
  }
  catch (PathfulException& pe) {

    pe.pushPath(fullid(*this));
    throw;
  }

  // Search tracker and get its properties
  accept(*m_modulesSetVisitor);
  accept(*m_cntNameVisitor);

  cleanup();
  builtok(true);
}

//
// Setup: link functions to various tracker related properties (use setup functions for ReadOnly Computable properties)
//
void Tracker::setup() {

  maxR.setup([&]() {
    double max = 0;
    for (const auto& b : m_barrels) max = MAX(max, b.maxR());
    for (const auto& e : m_endcaps) max = MAX(max, e.maxR());
    return max;
  });
  minR.setup([&]() {
    double min = std::numeric_limits<double>::max();
    for (const auto& b : m_barrels) min = MIN(min, b.minR());
    for (const auto& e : m_endcaps) min = MIN(min, e.minR());
    return min;
  });
  maxZ.setup([&]() {
    double max = 0;
    for (const auto& b : m_barrels) max = MAX(max, b.maxZ());
    for (const auto& e : m_endcaps) max = MAX(max, e.maxZ());
    return max;
  });
  minEta.setup([&]() {
    double min = std::numeric_limits<double>::max();
    for (const auto& m : m_modulesSetVisitor->modules()) min = MIN(min, m->minEta());
    return min;
  });
  maxEta.setup([&]() {
    double max = -std::numeric_limits<double>::max();
    for (const auto& m : m_modulesSetVisitor->modules()) max = MAX(max, m->maxEta());
    return max;
  });
}

//
// GeometryVisitor pattern -> tracker visitable
//
void Tracker::accept(GeometryVisitor& v)
{
  v.visit(*this);
  for (auto& b : m_barrels) { b.accept(v); }
  for (auto& e : m_endcaps) { e.accept(v); }
  for (auto& s : m_supportStructures) { s.accept(v); }
}

//
// GeometryVisitor pattern -> tracker visitable (const. option)
//
void Tracker::accept(ConstGeometryVisitor& v) const {
  v.visit(*this);
  for (const auto& b : m_barrels) { b.accept(v); }
  for (const auto& e : m_endcaps) { e.accept(v); }
  for (const auto& s : m_supportStructures) { s.accept(v); }
}

//
// Return all tracker modules
//
const std::vector<DetectorModule*>& Tracker::modules() const { return m_modulesSetVisitor->modules(); }
