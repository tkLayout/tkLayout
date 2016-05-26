#include "Endcap.h"
#include "MessageLogger.h"
#include "SupportStructure.h"

using material::SupportStructure;

//
// Constructor - parse geometry config file using boost property tree & read-in Disk, Support nodes
//
Endcap::Endcap(double barrelOuterZ, const std::string& name, const PropertyTree& nodeProperty, const PropertyTree& treeProperty) :
 numDisks(     "numDisks"    , parsedAndChecked()),
 innerZ(       "innerZ"      , parsedOnly()),
 outerZ(       "outerZ"      , parsedAndChecked()),
 minZ(  string("minZ")       ),
 maxZ(  string("maxZ")       ),
 minR(  string("minR")       ),
 maxR(  string("maxR")       ),
 skipServices( "skipServices", parsedOnly(), false), // broken, do not use
 m_barrelGap(  "barrelGap"   , parsedOnly()),
 m_diskNode(   "Disk"        , parsedOnly()),
 m_supportNode("Support"     , parsedOnly())
{
  // Set barrel outerZ
  this->m_barrelOuterZ(barrelOuterZ);

  // Set the geometry config parameters
  this->myid(name);
  this->store(treeProperty);
  this->store(nodeProperty);
}

//
// Limit endcap geometry by eta cut
//
void Endcap::cutAtEta(double eta)
{
  for (auto& d : m_disks) d.cutAtEta(eta);
  m_disks.erase_if([](const Disk& d) { return d.numRings() == 0; });
  numDisks(m_disks.size());
}

//
// Build recursively individual subdetector systems: Disks -> rings -> modules -> private method called by constructor
//
void Endcap::build()
{
  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();

    // Calculate parameters
    if     (!innerZ.state())     innerZ(m_barrelOuterZ() + m_barrelGap());
    else if(m_barrelGap.state()) logWARNING("'innerZ' was set, ignoring 'barrelGap'");
    else                         logWARNING("'innerZ' was set, ignoring 'barrelGap' & 'barrelOuterZ'");

    if (innerZ() > outerZ()) throw PathfulException("Endcap innerZ higher than outerZ!");

    vector<double> maxDsDistances = findMaxDsDistances();
    double         alpha          = pow(outerZ()/innerZ(), 1/double(numDisks()-1)); // geometric progression factor

    // Build disks
    vector<Disk*> tdisks;

    for (int i = 1; i <= numDisks(); i++) {
      Disk* diskp = GeometryFactory::make<Disk>();
      diskp->myid(i);

      // Standard is to build & calculate parameters for the central disc (in the middle)
      diskp->buildZ((innerZ() + outerZ())/2);

      // Apply correct offset for each disc versus middle position
      double offset = pow(alpha, i-1) * innerZ();
      diskp->placeZ(offset);

      // Store parameters in a tree
      diskp->store(propertyTree());
      if (m_diskNode.count(i) > 0) diskp->store(m_diskNode.at(i));

      // To test the extreme cases -> one needs to test either first or last layer (based on parity)
      diskp->zHalfLength((outerZ()-innerZ())/2.);

      // Build
      diskp->build(maxDsDistances);

      // Mirror discs
      Disk* diskn = GeometryFactory::clone(*diskp);
      diskn->mirrorZ();

      tdisks.push_back(diskp);
      tdisks.push_back(diskn);
    }
    std::stable_sort(tdisks.begin(), tdisks.end(), [](Disk* d1, Disk* d2) { return d1->minZ() < d2->maxZ(); });
    for (Disk* d : tdisks) m_disks.push_back(d);
    
  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  // Supports defined within a Barrel
  for (auto& mapel : m_supportNode) {
    SupportStructure* s = new SupportStructure();
    s->store(propertyTree());
    s->store(mapel.second);
    s->buildInEndcap(*this);
    m_supportStructures.push_back(s);
  }

  cleanup();
  builtok(true);
}

//
// Calculate various endcap related properties -> private method called by constructor
//
void Endcap::setup()
{
  maxR.setup([this]() { double max = 0;                                  for (const auto& d : m_disks) { max = MAX(max, d.maxR()); } return max; });
  minR.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& d : m_disks) { min = MIN(min, d.minR()); } return min; });
  maxZ.setup([this]() { double max =-std::numeric_limits<double>::max(); for (const auto& d : m_disks) { if(d.maxZ() > 0 ) max = MAX(max, d.maxZ()); } return max; });
  minZ.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& d : m_disks) { if(d.minZ() > 0 ) min = MIN(min, d.minZ()); } return min; });
}

//
// GeometryVisitor pattern -> endcap visitable
//
void Endcap::accept(GeometryVisitor& v) {
  v.visit(*this);
  for (auto& d : m_disks) { d.accept(v); }
}

//
// GeometryVisitor pattern -> endcap visitable (const. option)
//
void Endcap::accept(ConstGeometryVisitor& v) const {
  v.visit(*this);
  for (const auto& d : m_disks) { d.accept(v); }
}

//
// Drill down into the property tree to find the maximum dsDistances
//
vector<double> Endcap::findMaxDsDistances()
{
  vector<double> maxDsDistances;
  double endcapDsDistance = propertyTree().get("dsDistance", 0.);

  PropertyNode<int> ringNode("");
  for (auto& tel : pair2range(propertyTree().equal_range("Ring"))) // scan Ring subtrees outside Disks
    ringNode.fromPtree(tel.second);
  for (auto& rnel : ringNode) {
    Property<double, NoDefault> ringDsDistance; //endcapDsDistance);
    for (auto& tel : pair2range(rnel.second.equal_range("dsDistance"))) ringDsDistance.fromPtree(tel.second);
    if (ringDsDistance.state()) {
      if (maxDsDistances.size() < rnel.first) maxDsDistances.resize(rnel.first);//, endcapDsDistance);
      maxDsDistances[rnel.first-1] = MAX(maxDsDistances[rnel.first-1], ringDsDistance());
    }
  }

  for (auto& dnel : m_diskNode) {
    double diskDsDistance = dnel.second.get("dsDistance", endcapDsDistance);
    ringNode.clear();
    for (auto& tel : pair2range(dnel.second.equal_range("Ring"))) ringNode.fromPtree(tel.second); // scan Ring subtrees inside Disks
    for (auto& rnel : ringNode) {
      Property<double, NoDefault> ringDsDistance; // (diskDsDistance);
      for (auto& tel : pair2range(rnel.second.equal_range("dsDistance"))) ringDsDistance.fromPtree(tel.second);
      if (maxDsDistances.size() < rnel.first) maxDsDistances.resize(rnel.first); //, diskDsDistance);
      if (ringDsDistance.state()) {
        maxDsDistances[rnel.first-1] = MAX(maxDsDistances[rnel.first-1], ringDsDistance());
      } else {
        maxDsDistances[rnel.first-1] = MAX(maxDsDistances[rnel.first-1], diskDsDistance);
      }
    }
  }
  maxDsDistances.push_back(endcapDsDistance); // adds a default element to be used in case disks need more rings than the vector specifies dsDistances for
  return maxDsDistances;
}


